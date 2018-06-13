#include "initial.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>

int main(int argc, char const *argv[])
{
    srand(parameters::seed);
    std::vector<Particle> parts;
    const Index N = parameters::N;
    const Index gridPoints = parameters::gp;
    const Real dr = parameters::dr;
    const Index steps = parameters::steps;

    // External magnetic field
    const Real Bx = parameters::Bx;
    const Real By = parameters::By;
    const Real Bz = parameters::Bz;

    // External electric field
    const Real Ex = parameters::Ex;
    const Real Ey = parameters::Ey;
    const Real Ez = parameters::Ez;

    std::valarray<Real> Bext = {Bx, By, Bz};
    std::valarray<Real> E;

    if (parameters::system == "two stream")
        parts = two_stream();
    else if (parameters::system == "random")
        parts = random_particles();    

    std::vector<Index> mobileParticles;
    std::vector<Real> rho_c;
    for (Index p = 0; p < N; p++)
    {
        rho_c.push_back(parts[p].q_ / (dr * dr));

        if (parts.at(p).move_)
            mobileParticles.push_back(p);
    }

    const Index Nm = mobileParticles.size();

    // Data output
    std::vector<Real> spaceX;
    // VecVec spaceX(Nm);
    // VecVec spaceY(Nm);
    // VecVec spaceZ(Nm);
    // VecVec velocityX(Nm);
    // VecVec velocityY(Nm);
    // VecVec velocityZ(Nm);
    // VecVec mesh;
    // VecVec rho(gridPoints * gridPoints);
    // VecVec phi(gridPoints * gridPoints);
    // VecVec fieldX(gridPoints * gridPoints);
    // VecVec fieldY(gridPoints * gridPoints);
    // VecVec energy;

    // Write mesh (only once)
    // for (Index i = 0; i < gridPoints; ++i)
    // {
    //     for (Index j = 0; j < gridPoints; ++j)
    //     {
    //         std::vector<Real> point = {i * dr, j * dr};
    //         mesh.push_back(point);
    //     }
    // }

    std::vector<Particle> finalParts;

    VecVal RHO;
    VecVal PHI;
    VecVecVal EFIELDn;
    VecVal EFIELDp;
    Real KE;
    Real FE;
    std::cout << "Simulation running..." << std::endl;
    std::clock_t t_0 = std::clock();
    Real ssTime; // Single step time
    Real diff;

    H5Fcreate(parameters::dataOutput.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    createOutput(parameters::dataOutput, "space_x", Nm, steps);
    createOutput(parameters::dataOutput, "space_y", Nm, steps);
    createOutput(parameters::dataOutput, "space_z", Nm, steps);

    createOutput(parameters::dataOutput, "velocity_x", Nm, steps);
    createOutput(parameters::dataOutput, "velocity_y", Nm, steps);
    createOutput(parameters::dataOutput, "velocity_z", Nm, steps);

    createOutput(parameters::dataOutput, "mesh", gridPoints * gridPoints, 2);
    createOutput(parameters::dataOutput, "rho", gridPoints * gridPoints, steps);
    createOutput(parameters::dataOutput, "phi", gridPoints * gridPoints, steps);
    createOutput(parameters::dataOutput, "fieldX", gridPoints * gridPoints, steps);
    createOutput(parameters::dataOutput, "fieldY", gridPoints * gridPoints, steps);
    createOutput(parameters::dataOutput, "energy", steps, 2);

    for (Index step = 0; step < steps; ++step)
    {
        RHO = density(parts, rho_c, dr);
        PHI = potential(RHO, dr);
        EFIELDn = EField_GP(PHI, dr);
        EFIELDp = EField_P(EFIELDn, parts, dr);

        if (step < 0.2 * steps)
            E = {Ex, Ey, Ez};

        else
            E = {0.0, 0.0, 0.0};

        if (step == 0)
            outPhase(-1.0, EFIELDp, E, Bext, parts);

        Boris(EFIELDp, E, Bext, parts);
        finalParts = parts;
        outPhase(1.0, EFIELDp, E, Bext, finalParts);

        KE = 0.0;
        FE = 0.0;

        for (Index p = 0; p < Nm; ++p)
        {
            spaceX.push_back(finalParts.at(p).position_[0]);
            // spaceY.at(p).push_back(finalParts.at(p).position_[1]);
            // spaceZ.at(p).push_back(finalParts.at(p).position_[2]);
            // velocityX.at(p).push_back(finalParts.at(p).velocity_[0]);
            // velocityY.at(p).push_back(finalParts.at(p).velocity_[1]);
            // velocityZ.at(p).push_back(finalParts.at(p).velocity_[2]);
            KE += finalParts[p].m_ * norm(finalParts[p].velocity_) * norm(finalParts[p].velocity_);    
        }

        writeData(parameters::dataOutput, "space_x", spaceX, step);
        spaceX.clear();

        KE *= 0.5;

        for (Index i = 0; i < gridPoints; ++i)
        {
            for (Index j = 0; j < gridPoints; ++j)
            {
                Index index = i * gridPoints + j;
                // rho.at(index).push_back(RHO[i][j]);
                // phi.at(index).push_back(PHI[i][j]);
                // fieldX.at(index).push_back(EFIELDn[i][j][0]);
                // fieldY.at(index).push_back(EFIELDn[i][j][1]);
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;
        std::vector<Real> energies = {KE, FE};
        // energy.push_back(energies);
        
        if (step == 0)
        {
            std::clock_t t_1 = std::clock();
            diff = Real(t_1 - t_0);
            ssTime = (diff / CLOCKS_PER_SEC);
        }

        std::cout << "Aproximate time remaining: " << ssTime * (steps - step) << " seconds." << std::endl;
    }
    
    std::cout << "Simulation finished!" << std::endl; 
    std::cout << "\n************************************************\n" << std::endl;

    // writeData(parameters::dataOutput, "space_x", spaceX);
    // writeData(parameters::dataOutput, "space_y", spaceY);
    // writeData(parameters::dataOutput, "space_z", spaceZ);

    // writeData(parameters::dataOutput, "velocity_x", velocityX);
    // writeData(parameters::dataOutput, "velocity_y", velocityY);
    // writeData(parameters::dataOutput, "velocity_z", velocityZ);

    // writeData(parameters::dataOutput, "mesh", mesh);
    // writeData(parameters::dataOutput, "rho", rho);
    // writeData(parameters::dataOutput, "phi", phi);
    // writeData(parameters::dataOutput, "fieldX", fieldX);
    // writeData(parameters::dataOutput, "fieldY", fieldY);
    // writeData(parameters::dataOutput, "energy", energy);

    return 0;
}