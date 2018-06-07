#include "cosa.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>

int main(int argc, char const *argv[])
{
    srand(parameters::seed);
    std::vector<Particle> parts;
    const int N = parameters::N;
    const int gridPoints = parameters::gp;
    const double dr = parameters::dr;
    const int steps = parameters::steps;

    // External magnetic field
    const double Bx = parameters::Bx;
    const double By = parameters::By;
    const double Bz = parameters::Bz;

    // External electric field
    const double Ex = parameters::Ex;
    const double Ey = parameters::Ey;
    const double Ez = parameters::Ez;

    std::valarray<double> B = {Bx, By, Bz};
    std::valarray<double> E;

    // if (int(std::sqrt(N)) * int(std::sqrt(N)) != N) 
    // {
    //     std::cout << "Something went wrong: N must be a power of 2." << std::endl;
    //     return 1;
    // }

    if (parameters::system == "two stream")
    {
        parts = two_stream();
    }
    else if (parameters::system == "random")
    {
        parts = random_particles();    
    }

    // Data output
    VecVal spaceX = valarraysVector(N, steps);
    VecVal spaceY = valarraysVector(N, steps);
    VecVal velocityX = valarraysVector(N, steps);
    VecVal velocityY = valarraysVector(N, steps);
    VecVal velocityZ = valarraysVector(N, steps);
    // std::vector<double> spaceX;
    // std::vector<double> spaceY;

    std::vector<double> rho_c;

    for (int p = 0; p < N; ++p)
    {
        rho_c.push_back(parts[p].q_ / (dr * dr));
    }

    // std::string mainFolder = parameters::folder;

    // This doesn't work on a Linux OS.
    // std::string directory = "~/Desktop/" + mainFolder;
    // std::string dataOutput = name;
    // std::vector<std::string> folders = {"/space", "/phaseSpace", "/velocities", "/rho", "/phi", "/Efield", "/energy"};
    
    // for (int f = 0; f < folders.size(); ++f)
    // {
    //     system(("mkdir -p " + directory + folders[f]).c_str());
    // }

    std::vector<Particle> finalParts;

    std::cout << "Simulation running..." << std::endl;
    std::clock_t t_0 = std::clock();
    double simulationTime;
    double diff;
    // std::ofstream energy;
    // energy.open(directory + "/energy/energy.dat");
    
    for (int step = 0; step < steps; ++step)
    {
        VecVal RHO = density(parts, rho_c, dr);
        VecVal PHI = potential(RHO, dr);
        VecVecVal EFIELDn = EField_GP(PHI, dr);
        VecVal EFIELDp = EField_P(EFIELDn, parts, dr);

        if (step < int(0.2 * steps)) 
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

        // std::ofstream phaseSpace;
        // std::ofstream space;
        // std::ofstream velocities;
        // std::ofstream electricField;
        // std::ofstream electricPotential;
        // std::ofstream chargeDensity;

        // phaseSpace.open(directory + "/phaseSpace/step" + std::to_string(step) + ".dat");
        // space.open(directory + "/space/step" + std::to_string(step) + ".dat");
        // velocities.open(directory + "/velocities/step" + std::to_string(step) + ".dat");
        // electricField.open(directory + "/Efield/step" + std::to_string(step) + ".dat");
        // electricPotential.open(directory + "/phi/step" + std::to_string(step) + ".dat");
        // chargeDensity.open(directory + "/rho/step" + std::to_string(step) + ".dat");

        double KE = 0.0;
        double FE = 0.0;

        for (size_t p = 0; p < N; ++p)
        {
            if (finalParts.at(p).move_)
            {
                spaceX.at(p)[step] = finalParts.at(p).position_[0];
                spaceY.at(p)[step] = finalParts.at(p).position_[1];
                velocityX.at(p)[step] = finalParts.at(p).velocity_[0];
                velocityY.at(p)[step] = finalParts.at(p).velocity_[1];
                velocityZ.at(p)[step] = finalParts.at(p).velocity_[2];
                // std::cout << "puta mierda" << std::endl;
                // phaseSpace << finalParts[p].position_[0] << " " << finalParts[p].velocity_[0] << "\n";
                // space << finalParts[p].position_[0] << " " << finalParts[p].position_[1] << "\n";
                // velocities << finalParts[p].velocity_[0] << " " << finalParts[p].velocity_[1] << " " << finalParts[p].velocity_[2] << "\n";
                // KE += finalParts[p].m_ * norm(finalParts[p].velocity_) * norm(finalParts[p].velocity_);
            }          
        }

        KE *= 0.5;

        // for (int i = 0; i < gridPoints; ++i)
        // {
        //     for (int j = 0; j < gridPoints; ++j)
        //     {
        //         electricField << i * dr << " " << j * dr << " " << EFIELDn[i][j][0] << " " << EFIELDn[i][j][1] << "\n";
        //         electricPotential << i * dr << " " << j * dr << " " << PHI[i][j] << "\n";
        //         chargeDensity << i * dr << " " << j * dr << " " << RHO[i][j] << "\n";
        //         FE += RHO[i][j] * PHI[i][j];
        //     }
        // }

        FE *= 0.5;

        // energy << step << " " << KE << " " << FE << "\n";
        
        // phaseSpace.close();
        // space.close();
        // velocities.close();
        // electricField.close();
        // electricPotential.close();
        // chargeDensity.close();
        
        if (step == 0)
        {
            std::clock_t t_1 = std::clock();
            diff = double(t_1 - t_0);
            simulationTime = (diff / CLOCKS_PER_SEC);
        }

        std::cout << "Aproximate time remaining: " << simulationTime * (steps - step) << " seconds." << std::endl;
    }
    
    // energy.close();  
    std::cout << "Simulation finished." << std::endl; 
    hid_t file_id;
    file_id = H5Fcreate(parameters::dataOutput.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    writeSpaceOrVelocities(parameters::dataOutput, "space", spaceX, spaceY, spaceX);
    writeSpaceOrVelocities(parameters::dataOutput, "velocities", velocityX, velocityY, velocityZ);
    return 0;
}