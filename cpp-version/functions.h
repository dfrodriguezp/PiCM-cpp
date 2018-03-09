#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <array>
#include <valarray>
#include <vector>
#include <complex>
#include <random>
#include "particle.cc"
#include "parameters.h"

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<std::valarray<double>> VecVal;
typedef std::vector<std::vector<std::valarray<double>>> VecVecVal;

void fft(CArray &x);
void ifft(CArray& x);
VecVal valarraysVector(const int rows, const int cols);
VecVal density(std::vector<Particle> parts, std::vector<double> rho_c, const double dr);
VecVal potential(VecVal rho, const double dr);
VecVecVal EField_GP(VecVal phi, const double dr);
VecVal EField_P(VecVecVal field, std::vector<Particle> parts, const double dr);
inline double norm(std::valarray<double> Array);
inline double cross(std::valarray<double> A, std::valarray<double> B);
inline double mod(double a, double b);
void Boris(VecVal E, std::valarray<double> B, std::vector<Particle>& parts);
void rewind(const double direction, VecVal E, std::valarray<double> B, std::vector<Particle>& parts);

#endif