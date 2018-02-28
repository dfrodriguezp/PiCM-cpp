#include <iostream>
#include "particle.cc"
#include <random>
#include <vector>
#include <cmath>

int main(int argc, char const *argv[])
{

    double vt = 0.1;
    double vd = 1.0;
    int N = 49;
    double L = 1.0;
    double dr = 0.1;

    if (int(std::sqrt(N)) * int(std::sqrt(N)) != N) {
        std::cout << "You are so sudaca !!!" << std::endl;
        return 1;
    }

    double n = N / (L * L);
    double Ne_Ni_ratio = 2.0;
    double Ne = N / Ne_Ni_ratio;
    double margin = dr / 10.0;

    std::mt19937_64 engine;
    std::normal_distribution<double> maxwell_right(vd, vt);
    std::normal_distribution<double> maxwell_left(-vd, vt);

    double deltaL = (L - 2.0 * margin) / (std::sqrt(N) - 1);

    for (double x = margin; x <= (L - margin); x += deltaL)
    {
        for (double y = margin; y <= (L - margin); y += deltaL)
        {
            std::cout << x << " " << y << std::endl;
        }
    }

    // Particle a({0.0, 0.0}, {1.0, 0.0}, 1.0, 1.0, true);
    std::vector<Particle> parts;




    return 0;
}