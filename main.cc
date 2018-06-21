#include "initial.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iostream>

void printVecVal(VecVal& A) {
    for (Index i = 0; i < A.size(); i++) {
        std::cout << "[";
        for (Index j = 0; j < A.at(i).size(); j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << "\b]" << std::endl;
    }
}

int main(int argc, char const *argv[]) {
    srand(params::seed);
    std::vector<Particle> parts;
    const Index N = params::N;
    const Index gp = params::gp;
    const Real dr = params::dr;
    const Index steps = params::steps;

    // External magnetic field
    const Real Bx = params::Bx;
    const Real By = params::By;
    const Real Bz = params::Bz;

    // External electric field
    const Real Ex = params::Ex;
    const Real Ey = params::Ey;
    const Real Ez = params::Ez;

    Array B = {Bx, By, Bz};
    Array E = {Ex, Ey, Ez};

    // External electric potential

    if (params::system == "two stream"){
        parts = two_stream();
    } else if (params::system == "random") {
        parts = random_particles();    
    }

    std::vector<Real> rho_c;
    std::vector<Index> mobileParticles;

    for (Index p = 0; p < N; p++) {
        rho_c.push_back(parts[p].q_ / (dr * dr));
        if (parts.at(p).move_) {
            mobileParticles.push_back(p);
        }
    }

    const Index Nm = mobileParticles.size();

    // Data output
    VecVec spaceX(Nm);
    VecVec spaceY(Nm);
    VecVec spaceZ(Nm);
    VecVec velocityX(Nm);
    VecVec velocityY(Nm);
    VecVec velocityZ(Nm);
    VecVec mesh;
    VecVec rho(gp * gp);
    VecVec phi(gp * gp);
    VecVec fieldX(gp * gp);
    VecVec fieldY(gp * gp);
    VecVec energy;

    // Write mesh (only once)
    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j)
            mesh.push_back({i * dr, j * dr});
    }

    std::vector<Particle> finalParts;

    std::cout << "\nSimulation running...\n" << std::endl;
    std::clock_t t_0 = std::clock();
    Real simulationTime;
    Real diff;

    VecVal      RHO, PHI, EFIELDp;
    VecVecVal   EFIELDn;
    Real        KE, FE;

    for (Index step = 0; step < steps; ++step) {
        RHO = density(parts, rho_c);
        PHI = potential(RHO);
        EFIELDn = EField_GP(PHI);
        EFIELDp = EField_P(EFIELDn, parts);

        if (step == 0) {
            outPhase(-1.0, EFIELDp, E, B, parts);
        }

        Boris(EFIELDp, E, B, parts);
        finalParts = parts;
        outPhase(1.0, EFIELDp, E, B, finalParts);
        KE = 0.0;
        FE = 0.0;

        for (size_t p = 0; p < Nm; ++p) {
            spaceX.at(p).push_back(finalParts.at(p).position_[0]);
            spaceY.at(p).push_back(finalParts.at(p).position_[1]);
            spaceZ.at(p).push_back(finalParts.at(p).position_[2]);
            velocityX.at(p).push_back(finalParts.at(p).velocity_[0]);
            velocityY.at(p).push_back(finalParts.at(p).velocity_[1]);
            velocityZ.at(p).push_back(finalParts.at(p).velocity_[2]);
            KE += finalParts[p].m_ * norm(finalParts[p].velocity_) * norm(finalParts[p].velocity_);    
        }

        KE *= 0.5;

        for (Index i = 0; i < gp; ++i) {
            for (Index j = 0; j < gp; ++j) {
                Index index = i * gp + j;
                rho.at(index).push_back(RHO[i][j]);
                phi.at(index).push_back(PHI[i][j]);
                fieldX.at(index).push_back(EFIELDn[i][j][0]);
                fieldY.at(index).push_back(EFIELDn[i][j][1]);
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;
        
        energy.push_back({KE, FE});
        
        if (step == 0) {
            std::clock_t t_1 = std::clock();
            diff = Real(t_1 - t_0);
            simulationTime = (diff / CLOCKS_PER_SEC);
        }

        std::cout << "Aproximate time remaining: " << simulationTime * (steps - step) << " seconds." << std::endl;
    }

    std::cout << "\nSimulation finished." << std::endl; 
    std::cout << "\n************************************************\n" << std::endl;
    std::cout << "Writing data...\n" << std::endl;

    H5Fcreate(params::dataOutput.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    writeData(params::dataOutput, "space_x", spaceX);
    writeData(params::dataOutput, "space_y", spaceY);
    writeData(params::dataOutput, "space_z", spaceZ);

    writeData(params::dataOutput, "velocity_x", velocityX);
    writeData(params::dataOutput, "velocity_y", velocityY);
    writeData(params::dataOutput, "velocity_z", velocityZ);

    writeData(params::dataOutput, "mesh", mesh);
    writeData(params::dataOutput, "rho", rho);
    writeData(params::dataOutput, "phi", phi);
    writeData(params::dataOutput, "fieldX", fieldX);
    writeData(params::dataOutput, "fieldY", fieldY);
    writeData(params::dataOutput, "energy", energy);

    std::cout << "\nDone!" << std::endl;
    std::cout << "\n************************************************\n" << std::endl;

    return 0;
}