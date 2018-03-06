#include <iostream>
#include <cmath>
#include <array>
#include <valarray>
#include <vector>
#include <complex>
#include "particle.cc"
#include "parameters.h"


typedef std::vector<std::valarray<double>> grid2D;

grid2D density(std::vector<Particle> parts, std::vector<double> rho_c, const double dr)
// grid2D density(const double dr)
{
    grid2D rho;
    std::valarray<double> rho_j(0.0, parameters::gp);

    for (int i = 0; i < parameters::gp; ++i)
    {
        rho.push_back(rho_j);
    }

    for (int p = 0; p < parameters::N; ++p)
    {
        int i = std::floor(parts[p].position_[0] / dr);
        int j = std::floor(parts[p].position_[1] / dr);

        double hx = parts[p].position_[0] - (i * dr);
        double hy = parts[p].position_[1] - (j * dr);

        rho[i][j] += rho_c[p] * (dr - hx) * (dr - hy);
        rho[i][j+1] += rho_c[p] * (dr - hx) * hy;
        rho[i+1][j] += rho_c[p] * hx * (dr - hy);
        rho[i+1][j+1] += rho_c[p] * hx * hy;
    }

    for (int u = 0; u < parameters::gp; ++u)
    {
        rho[u][parameters::gp - 1] = (rho[u][parameters::gp - 1] + rho[u][0]) * 0.5;
        rho[u][0] = rho[u][parameters::gp - 1];
    }

    for (int u = 0; u < parameters::gp; ++u)
    {
        rho[parameters::gp - 1][u] = (rho[parameters::gp - 1][u] + rho[0][u]) * 0.5;
        rho[0][u] = rho[parameters::gp - 1][u];
    }

    for (int i = 0; i < parameters::gp; ++i)
    {
        rho[i] /= (dr * dr);
    }

    return rho;
}

grid2D potential(grid2D rho, const double dr)
{
    grid2D phi;
    grid2D phi_k;
    grid2D rho_k;

    

    return phi;
}
int main(int argc, char const *argv[])
{
    // grid2D RHO;
    // const double dr = 0.1;
    // RHO = density(dr);

    // std::cout << "[" << std::endl;
    // for (int i = 0; i < parameters::gp; ++i)
    // {   
    //     std::cout << "[";
    //     for (int j = 0; j < parameters::gp; ++j)
    //     {
    //         std::cout << RHO[i][j] << ", ";
    //     }
    //     std::cout << "]" << std::endl;
    // }
    // std::cout << "]" << std::endl;
	return 0;
}