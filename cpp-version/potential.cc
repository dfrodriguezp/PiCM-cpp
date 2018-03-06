#include <complex>
#include <cmath>
#include <valarray>
#include "FFT.cc"
#include "parameters.h"

int main()
{   
    // std::cout << "FFT solution of Poisson's equation\n"
    //           << "----------------------------------\n";
    // std::cout << "Enter number of points in x or y: ";
    // int N;
    // std::cin >> N;
    // int n = 1;
    // while(n < N)
    // {
    //     n *= 2;
    // }
    // if (n != N)
    // {
    //     std::cout << "must be a power of 2, using " << (N=n) << std::endl;
    // }
    
    const int N = 64;
    double h = 1 / double(N - 1);

    double q = 10;
    Complex rho[N][N];

    for (int j = 0; j < N; ++j)
    {
        for (int k = 0; k < N; ++k)
        {
            if (j == N/2 && k == N/2)
            {
                rho[j][k] = q / (h * h);
            }
            else
            {
                rho[j][k] = 0.0;
            }
        }
    }    

    CArray f(N);

    for (int j = 0; j < N; ++j)
    {
        for (int k = 0; k < N; ++k)
        {
            f[k] = rho[j][k];
        }
        fft(f);
        for (int k = 0; k < N; ++k)
        {
            rho[j][k] = f[k];
        }        
    }    

    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {    
            f[j] = rho[j][k];
        }
        fft(f);
        for (int j = 0; j < N; ++j)
        {
            rho[j][k] = f[j];
        }
    }

    Complex i(0.0, 1.0);
    Complex W = std::exp(2.0 * constants::pi * i / double(N));
    Complex Wm = 1.0, Wn = 1.0;

    for (int m = 0; m < N; ++m)
    {
        for (int n = 0; n < N; ++n)
        {
            Complex denom = 4.0;
            denom -= Wm + 1.0/Wm + Wn + 1.0/Wn;
            if (denom != 0.0)
            {
                rho[m][n] *= h * h / denom;
            }
            Wn *= W;
        }
        Wm *= W;
    }

    for (int j = 0; j < N; ++j)
    {
        for (int k = 0; k < N; ++k)
        {
            f[k] = rho[j][k];
        }
        ifft(f);
        for (int k = 0; k < N; ++k)
        {
            rho[j][k] = f[k];
        }
    }

    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            f[j] = rho[j][k];
        }
        ifft(f);
        for (int j = 0; j < N; ++j)
        {
            rho[j][k] = f[j];
        }
    }

    for (int j = 0; j < N; ++j)
    {
        double x = j * h;
        for (int k = 0; k < N; ++k)
        {
            double y = k * h;
            std::cout << x << " " << y << " " << std::real(rho[j][k]) << "\n";
        }
        std::cout << "\n";
    }
    // std::clock_t t_0 = std::clock();
 //    const Complex test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
 //    CArray data(test, 8);
 
 //    // forward fft
 //    fft(data);
 
 //    std::cout << "fft" << std::endl;
 //    for (int i = 0; i < 8; ++i)
 //    {
 //        std::cout << data[i] << std::endl;
 //    }
 
 //    // inverse fft
 //    ifft(data);
 
 //    std::cout << std::endl << "ifft" << std::endl;
 //    for (int i = 0; i < 8; ++i)
 //    {
 //        std::cout << data[i] << std::endl;
 //    }
 //    std::clock_t t_1 = std::clock();
 //    std::cout << "CPU time: " << double(t_1 - t_0) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}