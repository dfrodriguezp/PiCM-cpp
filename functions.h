#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <array>
#include <valarray>
#include <vector>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <random>
#include <fstream>
#include <string>
#include "/usr/include/jsoncpp/json/json.h"


typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<std::valarray<double>> VecVal;
typedef std::vector<std::vector<std::valarray<double>>> VecVecVal;

void fft(CArray& x);
void ifft(CArray& x);
double sign(double& x);
VecVal valarraysVector(const int& rows, const int& cols);
VecVal density(const std::vector<std::valarray<double>>& positions, 
               const std::vector<double>& charges, const double& dr, 
               const int& gp, const int& N);

VecVal potential(const VecVal& rho, const double& dr, const int& gp);
VecVecVal EField_GP(const VecVal& phi, const double& dr, const int& gp);
VecVal EField_P(const VecVecVal& field, 
                const std::vector<std::valarray<double>>& positions,
                const std::vector<int>& moves, 
                const double& dr, const int& gp, const int& N);

double norm(const std::valarray<double>& Array);
inline std::valarray<double> cross(const std::valarray<double>& A, const std::valarray<double>& B);
inline double mod(const double& a, const double& b);
void Boris(std::vector<std::valarray<double>>& positions, 
           std::vector<std::valarray<double>>& velocities,
           const std::vector<double>& QoverM,
           const std::vector<int>& moves,
           const VecVal& E, const std::valarray<double>& extE, 
           const std::valarray<double>& B, 
           const double& L, const double& dt, const int& N);

void outphase(std::vector<std::valarray<double>>& velocities,
              const std::vector<double>& QoverM,
              const std::vector<int>& moves,
              const double& direction, const VecVal& E, 
              const std::valarray<double>& extE, 
              const std::valarray<double>& B, 
              const double& dt, const int& N);

#endif