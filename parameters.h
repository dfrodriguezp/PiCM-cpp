#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

namespace parameters
{
    // External magnetic field
    constexpr double Bx(0.0);
    constexpr double By(0.0);
    constexpr double Bz(0.1);

    // External electric field
    constexpr double Ex(0.0);
    constexpr double Ey(0.05);
    constexpr double Ez(0.0);

    constexpr int gp(64);
    constexpr int N(1000 * 1000);
    constexpr double dt(0.05);
    constexpr int steps(2000);
    constexpr int seed(1886615);
    constexpr double vt(1.0);
    constexpr double vd(5.0);
    constexpr double dr(1.0);
    const std::string system = "random"; // Either "two stream", "random", etc.
    const std::string folder = "random2";
}

namespace constants
{
    constexpr double pi(3.141592653589793);
}
#endif