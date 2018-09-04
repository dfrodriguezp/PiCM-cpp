#include "particle.h"

Index steps, seed, gp, N;
Real  L, dr, dt, vt, vd, Bx, By, Bz;
Str   samplefile;

int main(int argc, char const *argv[]) {

    if (argc < 2) {
        std::cout << "ERROR: please include a parameters json file." << std::endl;
    }

    Json::Value root;
    Json::Reader reader;
    Str jsonfile = argv[1];
    std::ifstream data(jsonfile, std::ifstream::binary);
    if (!reader.parse(data, root, false)) {
        std::cout << "json file not found." << std::endl;
    }

    steps = root.get("steps", 50).asInt();
    gp =    root.get("npoints", 8).asInt();
    dt =    root.get("dt", 0.1).asDouble();
    N =     root["N"].asInt();
    dr =    root["dr"].asDouble();
    vt =    root["vt"].asDouble();
    vd =    root["vd"].asDouble();
    Json::Value Bfield = root["B"];
    Bx =    Bfield.get("Bx", 0.0).asDouble();
    By =    Bfield.get("By", 0.0).asDouble();
    Bz =    Bfield.get("Bz", 0.0).asDouble();
    samplefile = root["output"].asString();

    std::vector<Particle> parts;

    std::ifstream sample(samplefile);
    Real x, y, z, vx, vy, vz, n_r, qm;
    Index move;
    
    for (Index i = 0; i < N; i++) {
        sample >> x >> y >> z >> vx >> vy >> vz >> n_r >> qm >> move;
        Array pos = {x, y, z};
        Array vel = {vx, vy, vz};
        parts.push_back(Particle(pos, vel, n_r, qm, move));
    }

    Array B = {Bx, By, Bz};
    L = dr * (gp - 1);
    std::vector<Particle> finalParts;
    
    std::vector<Str> folders = {"/positions", "/velocities", "/density", "/potential", "/Efield", "/energy"};

    for (Index f = 0; f < folders.size(); ++f) {
        std::system(("mkdir -p results" + folders[f]).c_str());
    }

    std::cout << "\nSimulation running...\n" << std::endl;
    std::clock_t t_0 = std::clock();
    Real simulationTime;
    Real diff;

    VecArr RHO, PHI, EFIELDp;
    VecVecArr EFIELDn;
    Real KE, FE;
    std::ofstream energy;
    std::ofstream positions;
    std::ofstream velocities;
    std::ofstream density;
    std::ofstream potential;
    std::ofstream Efield;
    energy.open("results/energy/energy.dat");

    for (Index step = 0; step < steps; ++step) {
        RHO = update_density(parts);
        PHI = update_potential(RHO);
        EFIELDn = field_n(PHI);
        EFIELDp = field_p(EFIELDn, parts);

        if (step == 0) {
            outphase(-1.0, EFIELDp, B, parts);
        }

        Boris(EFIELDp, B, parts);
        finalParts = parts;
        outphase(1.0, EFIELDp, B, finalParts);

        if (mod(step, 10) == 0) {
            positions.open("results/positions/step_" + std::to_string(step) + ".dat");
            velocities.open("results/velocities/step_" + std::to_string(step) + ".dat");
            density.open("results/density/step_" + std::to_string(step) + ".dat");
            potential.open("results/potential/step_" + std::to_string(step) + ".dat");
            Efield.open("results/Efield/step_" + std::to_string(step) + ".dat");
        }

        KE = 0.0;
        FE = 0.0;

        for (Index p = 0; p < N; ++p) {
            if (finalParts.at(p).move_) {
                if (mod(step, 10) == 0) {
                    positions << finalParts.at(p).position_[0] << " " << finalParts.at(p).position_[1] << "\n";
                    velocities << finalParts.at(p).velocity_[0] << " " << finalParts.at(p).velocity_[1] << "\n";
                }
                KE += finalParts.at(p).mass_ * norm(finalParts.at(p).velocity_) * norm(finalParts.at(p).velocity_);
            }
        }

        KE *= 0.5;

        for (Index i = 0; i < gp; ++i) {
            for (Index j = 0; j < gp; ++j) {
                if (mod(step, 10) == 0) {
                    density << i * dr << " " << j * dr << " " << RHO[i][j] << "\n";
                    potential << i * dr << " " << j * dr << " " << PHI[i][j] << "\n";
                    Efield << i * dr << " " << j * dr << " " << EFIELDn[i][j][0] << " " << EFIELDn[i][j][1] << "\n";
                }
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;
        energy << step << " " << KE << " " << FE << "\n";
        positions.close();
        velocities.close();
        density.close();
        potential.close();
        Efield.close();

        if (step == 0) {
            std::clock_t t_1 = std::clock();
            diff = Real(t_1 - t_0);
            simulationTime = (diff / CLOCKS_PER_SEC);
        }
        std::cout << "Aproximate time remaining: " << simulationTime * (steps - step) << " seconds..." << std::endl;
    }

    energy.close();
    std::cout << "\nSimulation finished!\n" << std::endl; 
    std::cout << "\n*****************************************************************************\n" << std::endl;
    // std::cout << "Writing data...\n" << std::endl;

    // H5Fcreate(output.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // writeData(output, "position_x", pos_x);
    // writeData(output, "position_y", pos_y);

    // writeData(output, "velocity_x", vel_x);
    // writeData(output, "velocity_y", vel_y);
    // writeData(output, "velocity_z", vel_z);

    // writeData(output, "mesh", mesh_write);
    // writeData(output, "rho", rho);
    // writeData(output, "phi", phi);
    // writeData(output, "field_x", field_x);
    // writeData(output, "field_y", field_y);
    // writeData(output, "energy", energy);

    // std::cout << "\nDone!" << std::endl;
    // std::cout << "\n*****************************************************************************\n" << std::endl;

    return 0;
}