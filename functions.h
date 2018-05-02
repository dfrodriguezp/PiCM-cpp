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

std::vector<Particle> two_stream();
std::vector<Particle> random_particles();
void fft(CArray& x);
void ifft(CArray& x);
VecVal valarraysVector(const int& rows, const int& cols);
VecVal density(std::vector<Particle>& parts, const std::vector<double>& rho_c, const double& dr);
VecVal potential(const VecVal& rho, const double& dr);
VecVecVal EField_GP(const VecVal& phi, const double& dr);
VecVal EField_P(const VecVecVal& field, const std::vector<Particle>& parts, const double& dr);
double norm(const std::valarray<double>& Array);
inline std::valarray<double> cross(const std::valarray<double>& A, const std::valarray<double>& B);
inline double mod(const double& a, const double& b);
void Boris(const VecVal& E, const std::valarray<double>& B, std::vector<Particle>& parts);
void rewind(const double& direction, const VecVal& E, const std::valarray<double>& B, std::vector<Particle>& parts);

#endif