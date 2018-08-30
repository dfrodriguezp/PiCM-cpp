#include "mesh.h"

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
    Mesh mesh = Mesh();

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
    
    // VecVec pos_x(N);
    // VecVec pos_y(N);
    // VecVec pos_z(N);
    // VecVec vel_x(N);
    // VecVec vel_y(N);
    // VecVec vel_z(N);
    // VecVec mesh_write;
    // VecVec rho(n * n);
    // VecVec phi(n * n);
    // VecVec field_x(n * n);
    // VecVec field_y(n * n);
    // VecVec energy;

    // // Write mesh 
    // for (Index i = 0; i < n; ++i) {
    //     for (Index j = 0; j < n; ++j)
    //         mesh_write.push_back({i * dr, j * dr});
    // }

    std::vector<Str> folders = {"/positions", "/velocities", "/density", "/potential", "/Efield", "/energy"};

    for (Index f = 0; f < folders.size(); ++f) {
        std::system(("mkdir -p results" + folders[f]).c_str());
    }

    std::cout << "\nSimulation running...\n" << std::endl;
    std::clock_t t_0 = std::clock();
    Real simulationTime;
    Real diff;

    Real KE, FE;
    std::ofstream energy;
    std::ofstream positions;
    std::ofstream velocities;
    std::ofstream density;
    std::ofstream potential;
    std::ofstream Efield;
    energy.open("results/energy/energy.dat");

    for (Index step = 0; step < steps; ++step) {
        mesh.updateDensity(parts);
        mesh.updatePotential();
        mesh.updateEfield();

        for (Index p = 0; p < parts.size(); ++p) {
            parts[p].fieldInfluence(mesh.getEfield());
            if (step == 0) {
                parts[p].outPhase(-1, B);
            }
            parts[p].update(B);
        }

        finalParts = parts;

        for (Index p = 0; p < finalParts.size(); ++p) {
            finalParts[p].outPhase(1, B);
        }

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
            if (finalParts.at(p).Move() == 1) {
                if (mod(step, 10) == 0) {
                    positions << finalParts.at(p).getPosition()[0] << " " << finalParts.at(p).getPosition()[1] << "\n";
                    velocities << finalParts.at(p).getVelocity()[0] << " " << finalParts.at(p).getVelocity()[1] << "\n";
                }
                KE += finalParts.at(p).getMass() * norm(finalParts.at(p).getVelocity()) * norm(finalParts.at(p).getVelocity());
            }
            // pos_x.at(p).push_back(finalParts.at(p).getPosition()[0]);
            // pos_y.at(p).push_back(finalParts.at(p).getPosition()[1]);
            // vel_x.at(p).push_back(finalParts.at(p).getVelocity()[0]);
            // vel_y.at(p).push_back(finalParts.at(p).getVelocity()[1]);
            // vel_z.at(p).push_back(finalParts.at(p).getVelocity()[2]);
        }

        KE *= 0.5;

        for (Index i = 0; i < gp; ++i) {
            for (Index j = 0; j < gp; ++j) {
                if (mod(step, 10) == 0) {
                    density << i * dr << " " << j * dr << " " << mesh.getRho()[i][j] << "\n";
                    potential << i * dr << " " << j * dr << " " << mesh.getPhi()[i][j] << "\n";
                    Efield << i * dr << " " << j * dr << " " << mesh.getEfield()[i][j][0] << " " << mesh.getEfield()[i][j][1] << "\n";
                }
                // Index index = i * n + j;
                // rho.at(index).push_back(mesh.getRho()[i][j]);
                // phi.at(index).push_back(mesh.getPhi()[i][j]);
                // field_x.at(index).push_back(mesh.getEfield()[i][j][0]);
                // field_y.at(index).push_back(mesh.getEfield()[i][j][1]);

                FE += mesh.getRho()[i][j] * mesh.getPhi()[i][j];
            }
        }

        FE *= 0.5;
        energy << step << " " << KE << " " << FE << "\n";
        // energy.push_back({KE, FE});
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
    std::cout << "\n*****************************************************************************\n" << std::endl;

    return 0;
}