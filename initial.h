#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <array>
#include <valarray>
#include <vector>
#include <complex>
#include <random>
#include <string>
#include "/usr/include/hdf5/serial/hdf5.h"
#include "particle.cc"
#include "parameters.h"

typedef size_t Index;
typedef double Real;
typedef std::vector<std::valarray<Real>> VecVal;
typedef std::vector<std::vector<std::valarray<Real>>> VecVecVal;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;

std::vector<Particle> two_stream();
std::vector<Particle> random_particles();
void fft(CArray& x);
void ifft(CArray& x);
VecVal valarraysVector(const Index& rows, const Index& cols);
VecVal density(std::vector<Particle>& parts, const std::vector<Real>& rho_c, const Real& dr);
VecVal potential(const VecVal& rho, const Real& dr);
VecVecVal EField_GP(const VecVal& phi, const Real& dr);
VecVal EField_P(const VecVecVal& field, const std::vector<Particle>& parts, const Real& dr);
Real norm(const std::valarray<Real>& Array);
inline std::valarray<Real> cross(const std::valarray<Real>& A, const std::valarray<Real>& B);
inline Real mod(const Real& a, const Real& b);
void Boris(const VecVal& E, const std::valarray<Real>& extE, const std::valarray<Real>& B, std::vector<Particle>& parts);
void outPhase(const Real& direction, const VecVal& E, const std::valarray<Real>& extE, const std::valarray<Real>& B, std::vector<Particle>& parts);

void writeData(std::string filename, std::string sname, const VecVal& data);

#endif
