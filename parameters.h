#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

namespace parameters {
    // External magnetic field
    constexpr double Bx(0.0);
    constexpr double By(0.0);
    constexpr double Bz(0.1);

    // External electric field
    constexpr double Ex(0.0);
    constexpr double Ey(0.0);
    constexpr double Ez(0.0);

    constexpr int gp(64);
    constexpr int N(100 * 100);
    constexpr double dt(0.05);
    constexpr int steps(1000);
    constexpr int seed(1886615);
    constexpr double vt(1.0);
    constexpr double vd(5.0);
    constexpr double dr(1.0);
    const std::string system = "two stream"; // Either "two stream", "random", etc.
    const std::string dataOutput = "prueba.h5";

}

namespace constants {
    constexpr double pi(3.141592653589793);
}

#endif