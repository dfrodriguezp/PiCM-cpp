#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

typedef double Real;
typedef size_t Index;
namespace params {
    // External magnetic field
    constexpr Real Bx(0.0);
    constexpr Real By(0.0);
    constexpr Real Bz(0.1);

    // External electric field
    constexpr Real Ex(0.0);
    constexpr Real Ey(0.0);
    constexpr Real Ez(0.0);

    constexpr Index gp(8);
    constexpr Index N(10 * 10);
    constexpr Real dt(0.05);
    constexpr Index steps(10);
    constexpr Index seed(1886615);
    constexpr Real vt(1.0);
    constexpr Real vd(5.0);
    constexpr Real dr(1.0);
    const std::string system = "two stream"; // Either "two stream", "random", etc.
    const std::string dataOutput = "pruebac.h5";

}

namespace constants {
    constexpr Real pi(3.141592653589793);
}

#endif