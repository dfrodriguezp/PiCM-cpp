#ifndef PARAMETERS_H
#define PARAMETERS_H

namespace parameters
{
	constexpr double bx(0.0);
	constexpr double by(0.0);
	constexpr int gp(8);
	constexpr int N(10 * 10);
	constexpr double dt(0.1);
	constexpr double p_vt(0.001526);
	constexpr double p_vd(0.03125);
	constexpr double L(1.0);
	constexpr int steps(1);
	constexpr int seed(69);
}

namespace constants
{
	constexpr double pi(3.141592653589793);
}
#endif