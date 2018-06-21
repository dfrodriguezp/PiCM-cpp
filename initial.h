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
#include "/usr/include/jsoncpp/json/json.h"
#include "parameters.h"
#include "particle.cc"

typedef size_t Index;
typedef double Real;
typedef std::valarray<Real> Array;
typedef std::vector<Array> VecVal;
typedef std::vector<VecVal> VecVecVal;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<std::vector<Real>> VecVec;

std::vector<Particle> two_stream();
std::vector<Particle> random_particles();
void fft(CArray& x);
void ifft(CArray& x);
VecVal valarraysVector(const Index& rows, const Index& cols, const Real& value);
VecVal density(std::vector<Particle>& parts, const std::vector<Real>& rho_c);
VecVal potential(const VecVal& rho);
VecVecVal EField_GP(const VecVal& phi);
VecVal EField_P(const VecVecVal& field, const std::vector<Particle>& parts);
Real norm(const Array& A);
inline Array cross(const Array& A, const Array& B);
inline Real mod(const Real& a, const Real& b);
void Boris(const VecVal& E, const Array& extE, const Array& B, std::vector<Particle>& parts);
void outPhase(const Real& direction, const VecVal& E, const Array& extE, const Array& B, std::vector<Particle>& parts);

void writeData(std::string filename, std::string sname, const VecVec& data);

#endif
