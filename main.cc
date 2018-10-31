#include "functions.h"

int main(int argc, char const *argv[])
{
    if (argc < 2) 
    {
        std::cout << "ERROR: please include a parameters json file." << std::endl;
    }

    Json::Value root;
    Json::Reader reader;
    std::string jsonfile = argv[1];
    std::ifstream data(jsonfile, std::ifstream::binary);
    if (!reader.parse(data, root, false)) 
    {
        std::cout << "json file not found." << std::endl;
    }

    Json::Value Bfield = root["Bfield"];
    Json::Value Efield = root["Efield"];
    std::string samplefile = root["sample"].asString();
    std::string outputName = root["output"].asString();
    int     steps = root.get("steps", 50).asInt();
    int     gp =    root.get("npoints", 16).asInt();
    int     N =     root["N"].asInt();
    double  dt =    root.get("dt", 0.1).asDouble();
    double  dr =    root["dr"].asDouble();
    double  Bx =    Bfield.get("Bx", 0.0).asDouble();
    double  By =    Bfield.get("By", 0.0).asDouble();
    double  Bz =    Bfield.get("Bz", 0.0).asDouble();
    double  Ex =    Efield.get("Ex", 0.0).asDouble();
    double  Ey =    Efield.get("Ey", 0.0).asDouble();
    double  Ez =    Efield.get("Ez", 0.0).asDouble();

    std::vector<std::valarray<double>> positions;
    std::vector<std::valarray<double>> velocities;
    std::vector<std::valarray<double>> new_velocities;
    std::vector<double> charges;
    std::vector<double> QoverM;
    std::vector<double> masses;
    std::vector<int> moves;

    std::ifstream sample(samplefile);
    double x, y, z, vx, vy, vz, charge;
    int move;
    
    for (int i = 0; i < N; ++i) 
    {
        sample >> x >> y >> z >> vx >> vy >> vz >> charge >> move;
        std::valarray<double> pos = {x, y, z};
        std::valarray<double> vel = {vx, vy, vz};
        positions.push_back(pos);
        velocities.push_back(vel);
        charges.push_back(charge);
        QoverM.push_back(sign(charge));
        masses.push_back(charge / sign(charge));
        moves.push_back(move);
    }

    std::valarray<double> B = {Bx, By, Bz};
    std::valarray<double> E;
    double L = dr * double(gp);

    std::string directory = "/home/daniel/Desktop/" + outputName;
    std::vector<std::string> folders = {"/phaseSpace", "/phi", "/Efield", "/energy"};
    
    for (int f = 0; f < folders.size(); ++f)
    {
        std::system(("mkdir -p " + directory + folders[f]).c_str());
    }

    std::cout << "Simulation running..." << std::endl;
    double simulationTime;
    double diff;
    std::ofstream energy;
    energy.open(directory + "/energy/energy.dat");
    
    for (int step = 0; step < steps; ++step)
    {
        std::clock_t t_0 = std::clock();
        VecVal RHO = density(positions, charges, dr, gp, N);
        VecVal PHI = potential(RHO, dr, gp);
        VecVecVal EFIELDn = EField_GP(PHI, dr, gp);
        VecVal EFIELDp = EField_P(EFIELDn, positions, moves, dr, gp, N);

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
            outphase(velocities, QoverM, moves, -1.0, EFIELDp, E, B, dt, N);
        }

        Boris(positions, velocities, QoverM, moves, EFIELDp, E, B, L, dt, N);
        new_velocities = velocities;
        outphase(new_velocities, QoverM, moves, 1.0, EFIELDp, E, B, dt, N);

        std::ofstream phaseSpace;
        std::ofstream electricField;
        std::ofstream electricPotential;

        phaseSpace.open(directory + "/phaseSpace/step" + std::to_string(step) + ".dat");
        electricField.open(directory + "/Efield/step" + std::to_string(step) + ".dat");
        electricPotential.open(directory + "/phi/step" + std::to_string(step) + ".dat");

        double KE = 0.0;
        double FE = 0.0;

        for (int p = 0; p < N; ++p)
        {
            if (moves.at(p))
            {
                phaseSpace << positions.at(p)[0] << " " << new_velocities.at(p)[0] << "\n";
                KE += masses.at(p) * norm(new_velocities.at(p)) * norm(new_velocities.at(p));
            }          
        }

        KE *= 0.5;

        for (int i = 0; i < gp; ++i)
        {
            for (int j = 0; j < gp; ++j)
            {
                electricField << i * dr << " " << j * dr << " " << EFIELDn[i][j][0] << " " << EFIELDn[i][j][1] << "\n";
                electricPotential << i * dr << " " << j * dr << " " << PHI[i][j] << "\n";
                FE += RHO[i][j] * PHI[i][j];
            }
        }

        FE *= 0.5;

        energy << step << " " << KE << " " << FE << "\n";
        
        phaseSpace.close();
        electricField.close();
        electricPotential.close();
        
        if (step == 0)
        {
            std::clock_t t_1 = std::clock();
            diff = double(t_1 - t_0);
            simulationTime = (diff / CLOCKS_PER_SEC);
        }

        std::cout << "ETR: " << simulationTime * (steps - step) << " seconds..." << std::endl;
    }
    
    energy.close();  

    std::cout << "Simulation finished." << std::endl; 
    return 0;
}