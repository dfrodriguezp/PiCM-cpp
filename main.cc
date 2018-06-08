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

    std::valarray<Real> B = {Bx, By, Bz};
    std::valarray<Real> E;

    if (parameters::system == "two stream")
    {
        parts = two_stream();
    }
    else if (parameters::system == "random")
    {
        parts = random_particles();    
    }

    std::vector<Index> mobileParticles;
    for (Index p = 0; p < N; p++)
    {
        if (parts.at(p).move_)
        {
            mobileParticles.push_back(p);
        }
    }

    const Index Nm = mobileParticles.size();

    // Data output
    VecVal spaceX = valarraysVector(Nm, steps);
    VecVal spaceY = valarraysVector(Nm, steps);
    VecVal spaceZ = valarraysVector(Nm, steps);
    VecVal velocityX = valarraysVector(Nm, steps);
    VecVal velocityY = valarraysVector(Nm, steps);
    VecVal velocityZ = valarraysVector(Nm, steps);
    VecVal mesh;
    VecVal rho = valarraysVector(gridPoints*gridPoints, steps);
    VecVal phi = valarraysVector(gridPoints*gridPoints, steps);
    VecVal fieldX = valarraysVector(gridPoints*gridPoints, steps);
    VecVal fieldY = valarraysVector(gridPoints*gridPoints, steps);
    VecVal energy;

    std::vector<Real> rho_c;

    // Write mesh (only once)
    for (Index i = 0; i < gridPoints; ++i)
    {
        for (Index j = 0; j < gridPoints; ++j)
            mesh.push_back({i * dr, j * dr});
    }

    for (Index p = 0; p < N; ++p)
    {
        rho_c.push_back(parts[p].q_ / (dr * dr));
    }

    std::vector<Particle> finalParts;

    std::cout << "Simulation running..." << std::endl;
    std::clock_t t_0 = std::clock();
    Real simulationTime;
    Real diff;
    
    for (Index step = 0; step < steps; ++step)
    {
        VecVal RHO = density(parts, rho_c, dr);
        VecVal PHI = potential(RHO, dr);
        VecVecVal EFIELDn = EField_GP(PHI, dr);
        VecVal EFIELDp = EField_P(EFIELDn, parts, dr);

        if (step < Index(0.2 * steps)) 
        {
            E = {Ex, Ey, Ez};
        }

        else
        {
            E = {0.0, 0.0, 0.0};
        }

        if (step == 0)
        {
            outPhase(-1.0, EFIELDp, E, B, parts);
        }

        Boris(EFIELDp, E, B, parts);
        finalParts = parts;
        outPhase(1.0, EFIELDp, E, B, finalParts);

        Real KE = 0.0;
        Real FE = 0.0;

        for (size_t p = 0; p < Nm; ++p)
        {
            spaceX.at(p)[step] = finalParts.at(p).position_[0];
            spaceY.at(p)[step] = finalParts.at(p).position_[1];
            spaceZ.at(p)[step] = finalParts.at(p).position_[2];
            velocityX.at(p)[step] = finalParts.at(p).velocity_[0];
            velocityY.at(p)[step] = finalParts.at(p).velocity_[1];
            velocityZ.at(p)[step] = finalParts.at(p).velocity_[2];
            KE += finalParts[p].m_ * norm(finalParts[p].velocity_) * norm(finalParts[p].velocity_);    
        }

        KE *= 0.5;

        for (Index i = 0; i < gridPoints; ++i)
        {
            for (Index j = 0; j < gridPoints; ++j)
            {
                Index index = i * gridPoints + j;
                rho.at(index)[step] = RHO[i][j];
                phi.at(index)[step] = PHI[i][j];
                fieldX.at(index)[step] = EFIELDn[i][j][0];
                fieldY.at(index)[step] = EFIELDn[i][j][1];
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;

        energy.push_back({KE, FE});
        
        if (step == 0)
        {
            std::clock_t t_1 = std::clock();
            diff = Real(t_1 - t_0);
            simulationTime = (diff / CLOCKS_PER_SEC);
        }

        std::cout << "Aproximate time remaining: " << simulationTime * (steps - step) << " seconds." << std::endl;
    }
    
    std::cout << "Simulation finished." << std::endl; 
    std::cout << "\n************************************************\n" << std::endl;
    std::cout << "Writing data..." << std::endl;

    hid_t file_id;
    file_id = H5Fcreate(parameters::dataOutput.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    writeData(parameters::dataOutput, "space_x", spaceX);
    writeData(parameters::dataOutput, "space_y", spaceY);
    writeData(parameters::dataOutput, "space_z", spaceZ);

    writeData(parameters::dataOutput, "velocity_x", velocityX);
    writeData(parameters::dataOutput, "velocity_y", velocityY);
    writeData(parameters::dataOutput, "velocity_z", velocityZ);

    writeData(parameters::dataOutput, "mesh", mesh);
    writeData(parameters::dataOutput, "rho", rho);
    writeData(parameters::dataOutput, "phi", phi);
    writeData(parameters::dataOutput, "fieldX", fieldX);
    writeData(parameters::dataOutput, "fieldY", fieldY);
    writeData(parameters::dataOutput, "energy", energy);

    std::cout << "Done!" << std::endl;
    std::cout << "\n************************************************\n" << std::endl;

    return 0;
}