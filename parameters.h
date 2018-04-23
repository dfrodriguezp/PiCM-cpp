#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

namespace parameters
{
    constexpr double Bx(0.0);
    constexpr double By(0.0);
    constexpr double Bz(0.0);
    constexpr double Ex(0.0);
    constexpr double Ey(0.0);
    constexpr double Ez(0.0);
    constexpr int gp(64 * 64);
    constexpr int N(100 * 100);
    constexpr double dt(0.025);
    constexpr int steps(1500);
    constexpr int seed(615821);
    const std::string folder = "bistream_run4";
    constexpr double vt1(10.0);
    constexpr double vt2(1.0);
    constexpr double vd1(0.0);
    constexpr double vd2(20.0);
    constexpr double R(0.05);
    constexpr double dr(1.0);
}

namespace constants
{
    constexpr double pi(3.141592653589793);
}
#endif