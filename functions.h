#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <random>
#include <fstream>
#include <string>
#include "/usr/include/jsoncpp/json/json.h"


typedef size_t Int;
typedef double Real;
typedef std::valarray<Real> Array;
typedef std::vector<Array> VecArr;
typedef std::vector<VecArr> VecVecArr;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;

Real sign(Real& x);

Real norm(const Array& A);

inline Array cross(const Array& A, const Array& B);

Real mod(const Real& a, const Real& b);

VecArr valarraysVector(const Int& rows, const Int& cols);

VecArr density(const VecArr& positions, 
               const std::vector<Real>& charges, const Real& dx, 
               const Real& dy, const Int& Nx, const Int& Ny, const Int& N);

VecArr potential(const VecArr& rho, const Real& dx, const Real& dy, 
                 const Int& Nx, const Int& Ny);

VecVecArr fieldNodes(const VecArr& phi, const Real& dx, const Real& dy, 
                     const Int& Nx, const Int& Ny);

VecArr fieldParticles(const VecVecArr& field, 
                      const VecArr& positions,
                      const std::vector<Int>& moves, 
                      const Real& dx, const Real& dy, 
                      const Int& Nx, const Int& Ny, const Int& N);

void boris(VecArr& velocities,
           const std::vector<Real>& QoverM,
           const std::vector<Int>& moves,
           const VecArr& E, const Array& B, 
           const Real& dt, const Int& N);

void update(VecArr& positions,
           VecArr& velocities,
           const std::vector<Real>& QoverM,
           const std::vector<Int>& moves,
           const VecArr& E, const Array& B, 
           const Real& Lx, const Real& Ly, 
           const Real& dt, const Int& N);

void outphase(const Real& direction,
              VecArr& velocities,
              const std::vector<Real>& QoverM,
              const std::vector<Int>& moves,
              const VecArr& E, const Array& B, 
              const Real& dt, const Int& N);

void fft(CArray& x);

void ifft(CArray& x);

#endif