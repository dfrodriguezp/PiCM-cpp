#ifndef INITIAL_H
#define INITIAL_H

#include "/usr/include/jsoncpp/json/json.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <valarray>
#include <complex>
#include <ctime>
#include <cstdlib>

typedef double Real;
typedef size_t Index;
typedef std::string Str;
typedef std::valarray<Real> Array;
typedef std::vector<Array> VecArr;
typedef std::vector<VecArr> VecVecArr;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<std::vector<Real>> VecVec;

extern Index steps, seed, gp, N;
extern Real  L, dr, dt, vt, vd, Bx, By, Bz;
extern Str   samplefile;

VecArr arraysVector(const Index& rows, const Index& cols, const Real& value);
VecVecArr arraysVectorVector(const Index& rows, const Index& cols, const Real& value);
Real norm(const Array& A);
inline Array cross(const Array& A, const Array& B);
Real mod(const Real& a, const Real& b);
void fft(CArray& x);
void ifft(CArray& x);
VecArr update_potential(const VecArr& rho);
VecVecArr field_n(const VecArr& phi);


#endif