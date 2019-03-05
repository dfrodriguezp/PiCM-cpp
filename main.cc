#include "init.h"

int main(int argc, char const *argv[])
{
    if (argc < 2) 
    {
        std::cout << "ERROR! Please include a parameters JSON file." << std::endl;
        return 1;
    }

    Json::Value root;
    Json::Reader reader;
    std::string jsonfile = argv[1];
    std::ifstream data(jsonfile, std::ifstream::binary);
    if (!reader.parse(data, root, false)) 
    {
        std::cout << "ERROR! JSON file not found or invalid." << std::endl;
        return 1;
    }

    Json::Value results = root["results"];
    bool    writePhaseSpace = isInArray("phase_space", results);
    bool    writeEfield = isInArray("electric_field", results);
    bool    writePhi = isInArray("electric_potential", results);
    bool    writeRho = isInArray("charge_density", results);
    bool    writeSpace = isInArray("space", results);

    std::string samplefile = root["sample"].asString();
    std::string outputName = root.get("output", "results").asString();

    Int N;
    if (root.isMember("N"))
        N = root["N"].asUInt();
    else
    {
        std::cout << "ERROR! Number of particles \"N\" not found in JSON file." << std::endl;
        return 1;
    }

    Int     steps = root.get("steps", 50).asUInt();
    Int     ss_freq = root.get("ss_frequency", 10).asUInt();
    Int     seed =  root.get("seed", 69696969).asUInt();

    Real    dt = root.get("dt", 0.1).asDouble();
    Real    dx = root.get("dx", 1.0).asDouble();
    Real    dy = root.get("dy", 1.0).asDouble();

    std::valarray<Real> B(3);

    Json::Value Bfield;
    if (root.isMember("Bfield"))
    {
        Bfield = root["Bfield"];
        if (Bfield.size() != 3)
        {
            std::cout << "ERROR! Magnetic field must have three components!" << std::endl;
            return 1;
        }
        else
        {
            B[0] = Bfield[0].asDouble();
            B[1] = Bfield[1].asDouble();
            B[2] = Bfield[2].asDouble();
        }
    }
    else
    {
        B = {0.0, 0.0, 0.0};
    }

    Int Nx, Ny;

    Json::Value gridSize;
    if (root.isMember("grid_size"))
    {
        gridSize = root["grid_size"];
        if (gridSize.size() != 2)
        {
            std::cout << "ERROR! Grid size must have two components!" << std::endl;
            return 1;
        }
        else
        {
            Nx = gridSize[0].asUInt();
            Ny = gridSize[1].asUInt();
        }
    }
    else
    {
        Nx = 16;
        Ny = 16;
    }

    Real Lx = dx * Real(Nx);
    Real Ly = dy * Real(Ny);
    VecArr positions;
    VecArr velocities;
    VecArr new_velocities;
    std::vector<Real> charges;
    std::vector<Real> QoverM;
    std::vector<Real> masses;
    std::vector<Int> moves;

    std::ifstream sample(samplefile);
    if (!sample.good())
    {
        std::cout << "ERROR! Sample file not found." << std::endl;
        return 1;
    }

    Real x, y, vx, vy, vz, q_m, charge;
    Int move;
    
    for (Int i = 0; i < N; ++i) 
    {
        sample >> x >> y >> vx >> vy >> vz >> q_m >> move;
        charge = (Lx * Ly * q_m) / N;
        std::valarray<Real> pos = {x, y}; // 2D
        std::valarray<Real> vel = {vx, vy, vz}; // 3V
        positions.push_back(pos);
        velocities.push_back(vel);
        QoverM.push_back(q_m);
        charges.push_back(charge);
        masses.push_back(charge / q_m);
        moves.push_back(move);
    }

    std::vector<std::string> folders = {"/energy"};

    if (writePhaseSpace) folders.push_back("/phase_space");
    if (writeEfield) folders.push_back("/Efield");
    if (writePhi) folders.push_back("/phi");
    if (writeRho) folders.push_back("/rho");
    if (writeSpace) folders.push_back("/space");

    for (Int f = 0; f < folders.size(); ++f)
        std::system(("mkdir -p " + outputName + folders[f]).c_str());

    std::cout << "Simulation running..." << std::endl;
    Real simulationTime;
    Real diff;
    std::ofstream energy;
    std::ofstream phaseSpace;
    std::ofstream electricField;
    std::ofstream electricPotential;
    std::ofstream chargeDensity;
    std::ofstream space;

    VecArr RHO, PHI, EFIELDp;
    VecVecArr EFIELDn;
    Real KE, FE;

    bool writeStep;
    energy.open(outputName + "/energy/energy_seed_" + std::to_string(seed) +"_.dat");

    for (Int step = 0; step < steps; ++step)
    {
        std::clock_t t_0 = std::clock();
        writeStep = mod(Real(step), Real(ss_freq)) == 0;
        RHO = density(positions, charges, dx, dy, Nx, Ny, N);
        PHI = potential(RHO, dx, dy, Nx, Ny);
        EFIELDn = fieldNodes(PHI, dx, dy, Nx, Ny);
        EFIELDp = fieldParticles(EFIELDn, positions, moves, dx, dy, Nx, Ny, N);

        if (step == 0)
            outphase(-1.0, velocities, QoverM, moves, EFIELDp, B, dt, N);

        update(positions, velocities, QoverM, moves, EFIELDp, B, Lx, Ly, dt, N);
        new_velocities = velocities;
        outphase(1.0, new_velocities, QoverM, moves, EFIELDp, B, dt, N);

        if (writePhaseSpace && writeStep)
            phaseSpace.open(outputName + "/phase_space/step_" + std::to_string(step) + "_seed_" + std::to_string(seed) + "_.dat");

        if (writeEfield && writeStep)
            electricField.open(outputName + "/Efield/step_" + std::to_string(step) + "_seed_" + std::to_string(seed) +  "_.dat");

        if (writePhi && writeStep)
            electricPotential.open(outputName + "/phi/step_" + std::to_string(step) + "_seed_" + std::to_string(seed) + "_.dat");

        if (writeRho && writeStep)
            chargeDensity.open(outputName + "/rho/step_" + std::to_string(step) + "_seed_" + std::to_string(seed) + "_.dat");

        if (writeSpace && writeStep)
            space.open(outputName + "/space/step_" + std::to_string(step) + "_seed_" + std::to_string(seed) + "_.dat");


        KE = 0.0;
        FE = 0.0;

        for (Int p = 0; p < N; ++p)
        {
            if (moves.at(p))
            {
                if (writePhaseSpace && writeStep)
                    phaseSpace << positions.at(p)[0] << " " << positions.at(p)[1] << " " <<
                                  new_velocities.at(p)[0] << " " << new_velocities.at(p)[1] << " " << new_velocities.at(p)[2] << "\n";
                KE += masses.at(p) * norm(new_velocities.at(p)) * norm(new_velocities.at(p));
            }          
            space << positions.at(p)[0] << " " << positions.at(p)[1] << "\n";
        }

        KE *= 0.5;

        for (Int i = 0; i < Nx; ++i)
        {
            for (Int j = 0; j < Ny; ++j)
            {
                if (writeEfield && writeStep)
                    electricField << i * dx << " " << j * dy << " " << EFIELDn[i][j][0] << " " << EFIELDn[i][j][1] << "\n";
                if (writePhi && writeStep)
                    electricPotential << i * dx << " " << j * dy << " " << PHI[i][j] << "\n";
                if (writeRho && writeStep)
                    chargeDensity << i * dx << " " << j * dy << " " << RHO[i][j] << "\n";
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;

        energy << step << " " << KE << " " << FE << "\n";
        
        if (writePhaseSpace && mod(step, ss_freq)) phaseSpace.close();
        if (writeEfield && mod(step, ss_freq)) electricField.close();
        if (writePhi && mod(step, ss_freq)) electricPotential.close();
        if (writeRho && mod(step, ss_freq)) chargeDensity.close();
        if (writeSpace && mod(step, ss_freq)) space.close();

        std::clock_t t_1 = std::clock();
        diff = Real(t_1 - t_0);
        simulationTime = (diff / CLOCKS_PER_SEC);

        if ((writePhaseSpace || writeEfield || writePhi || writeRho || writeSpace) && writeStep)
            std::cout << "Writing data of step " << step << "..." << std::endl;
        else
            std::cout << "ETR: " << simulationTime * (steps - step) << " seconds..." << std::endl;
    }
    
    energy.close();  

    std::cout << "Simulation finished." << std::endl; 
    return 0;
}