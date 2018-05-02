#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

namespace parameters
{
    constexpr double Bx(0.0);
    constexpr double By(0.0);
    constexpr double Bz(0.1);
    constexpr int gp(64);
    constexpr int N(710 * 710);
    constexpr double dt(0.05);
    constexpr int steps(2000);
    constexpr int seed(561765);
    constexpr double vt(1.0);
    constexpr double vd(5.0);
    constexpr double dr(1.0);
    const std::string system = "random"; // Either "two stream", "random", etc.
    const std::string folder = "random_particles_Bz";
}

namespace constants
{
    constexpr double pi(3.141592653589793);
}
#endif