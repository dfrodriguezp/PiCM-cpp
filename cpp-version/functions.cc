#include "functions.h"

typedef std::vector<std::valarray<double>> VecVal;
typedef std::vector<std::vector<std::valarray<double>>> VecVecVal;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

VecVal valarraysVector(const int rows, const int cols)
{
    VecVal f;
    std::valarray<double> zeros(0.0, cols);
    for (int i = 0; i < rows; ++i)
    {
        f.push_back(zeros);
    }

    return f;
}

VecVal density(std::vector<Particle> parts, std::vector<double> rho_c, const double dr)
{
    const int gp = parameters::gp;
    const int N = parameters::N;
    VecVal rho = valarraysVector(gp, gp);

    for (int p = 0; p < N; ++p)
    {
        int i = std::floor(parts[p].position_[0] / dr);
        int j = std::floor(parts[p].position_[1] / dr);
        double hx = parts[p].position_[0] - (i * dr);
        double hy = parts[p].position_[1] - (j * dr);

        rho[i][j] += rho_c[p] * (dr - hx) * (dr - hy);
        rho[i][j+1] += rho_c[p] * (dr - hx) * hy;
        rho[i+1][j] += rho_c[p] * hx * (dr - hy);
        rho[i+1][j+1] += rho_c[p] * hx * hy;
    }
    for (int u = 0; u < gp; ++u)
    {
        rho[gp - 1][u] = (rho[gp - 1][u] + rho[0][u]) * 0.5;
        rho[0][u] = rho[gp - 1][u];
    }

    for (int u = 0; u < gp; ++u)
    {
        rho[u][gp - 1] = (rho[u][gp - 1] + rho[u][0]) * 0.5;
        rho[u][0] = rho[u][gp - 1];
    }


    for (int i = 0; i < gp; ++i)
    {
        rho[i] /= (dr * dr);
    }

    return rho;
}

VecVal potential(VecVal rho, const double dr)
{   
    const int gp = parameters::gp;
    Complex rho_k[gp][gp];

    for (int i = 0; i < gp; ++i)
       {
           for (int j = 0; j < gp; ++j)
           {
               rho_k[i][j] = rho[i][j];
           }
       }   

    CArray f(gp);

    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            f[j] = rho_k[i][j];
        }
        fft(f);
        for (int j = 0; j < gp; ++j)
        {
            rho_k[i][j] = f[j];
        }
    }

    for (int j = 0; j < gp; ++j)
    {
        for (int i = 0; i < gp; ++i)
        {    
            f[i] = rho_k[i][j];
        }
        fft(f);
        for (int i = 0; i < gp; ++i)
        {
            rho_k[i][j] = f[i];
        }
    }

    Complex phi_k[gp][gp];

    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            phi_k[i][j] = rho_k[i][j];
        }
    }

    Complex i(0.0, 1.0);
    Complex W = std::exp(2.0 * constants::pi * i / double(gp));
    Complex Wm = 1.0;
    Complex Wn = 1.0;

    for (int m = 0; m < gp; ++m)
    {
        for (int n = 0; n < gp; ++n)
        {
            Complex denom = 4.0;
            denom -= Wm + 1.0/Wm + Wn + 1.0/Wn;
            if (denom != 0.0)
            {
                phi_k[m][n] *= (dr * dr) / denom;
            }
            Wn *= W;
        }
        Wm *= W;
    }

    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            f[j] = phi_k[i][j];
        }
        ifft(f);
        for (int j = 0; j < gp; ++j)
        {
            phi_k[i][j] = f[j];
        }
    }

    for (int j = 0; j < gp; ++j)
    {
        for (int i = 0; i < gp; ++i)
        {
            f[i] = phi_k[i][j];
        }
        ifft(f);
        for (int i = 0; i < gp; ++i)
        {
            phi_k[i][j] = f[i];
        }
    }

    VecVal phi = valarraysVector(gp, gp);
    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            phi[i][j] = std::real(phi_k[i][j]);
        }
    }

    return phi;
}

VecVecVal EField_GP(VecVal phi, const double dr)
{
    const int gp = parameters::gp;
    VecVecVal E;
    std::valarray<double> zeros(0.0, 2);
    for (int i = 0; i < gp; ++i)
    {   
        VecVal rows;
        for (int j = 0; j < gp; ++j)
        {
            rows.push_back(zeros);
        }
        E.push_back(rows);  
    }

    for (int j = 0; j < gp; ++j)
    {
        for (int i = 0; i < gp; ++i)
        {
            int nxt_i = (i < gp - 1) ? (i + 1) : 0;
            int prv_i = (i > 0) ? (i - 1) : (gp - 1);

            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (2 * dr);
        }
    }

    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            int nxt_j = (j < gp - 1) ? (j + 1) : 0;
            int prv_j = (j > 0) ? (j - 1) : (gp - 1);

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (2 * dr);
        }
    }

    return E;
}

VecVal EField_P(VecVecVal field, std::vector<Particle> parts, const double dr)
{
    const int N = parameters::N;
    VecVal E = valarraysVector(N, 2);

    int index = 0;

    for (int p = 0; p < N; ++p)
    {
        if (parts[p].move_)
        {
            int i = std::floor(parts[p].position_[0] / dr);
            int j = std::floor(parts[p].position_[1] / dr);
            double hx = parts[p].position_[0] - (i * dr);
            double hy = parts[p].position_[1] - (j * dr);

            double A = (dr - hx) * (dr - hy);
            double B = (dr - hx) * hy;
            double C = hx * (dr - hy);
            double D = hx * hy;

            E[index][0] = field[i][j][0] * A + field[i][j+1][0] * B + field[i+1][j][0] * C + field[i+1][j+1][0] * D;
            E[index][1] = field[i][j][1] * A + field[i][j+1][1] * B + field[i+1][j][1] * C + field[i+1][j+1][1] * D;
        }
        index++;
    }

    for (int i = 0; i < N; ++i)
    {
        E[i] /= (dr * dr);
    }

    return E;
}

inline double norm(std::valarray<double> Array)
{
    return std::sqrt((Array * Array).sum());
}

inline double cross(std::valarray<double> A, std::valarray<double> B)
{
    return (A[0] * B[1]) - (A[1] * B[0]);
}

inline double mod(double a, double b)
{
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}

void Boris(VecVal E, std::valarray<double> B, std::vector<Particle>& parts)
{
    int index = 0;
    const double L = parameters::L;
    const double dt = parameters::dt;
    const int N = parameters::N;
    std::valarray<double> t;
    double t_2;
    std::valarray<double> s;
    std::valarray<double> v_minus;
    std::valarray<double> v_prime;
    std::valarray<double> v_plus;

    for (int p = 0; p < N; ++p)
    {
        if (parts[p].move_)
        {
            t = 0.5 * (parts[p].qm_) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = parts[p].velocity_ + 0.5 * (parts[p].qm_) * E[index] * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            parts[p].velocity_ = v_plus + 0.5 * (parts[p].qm_) * E[index] * dt;

            parts[p].position_ += parts[p].velocity_ * dt;
            parts[p].position_[0] = mod(parts[p].position_[0], L);
            parts[p].position_[1] = mod(parts[p].position_[1], L);
        }
        index++;
    }
}

void rewind(const double direction, VecVal E, std::valarray<double> B, std::vector<Particle>& parts)
{
    int index = 0;
    const double dt = direction * 0.5 * parameters::dt;
    const int N = parameters::N;
    std::valarray<double> t;
    double t_2;
    std::valarray<double> s;
    std::valarray<double> v_minus;
    std::valarray<double> v_prime;
    std::valarray<double> v_plus;

    for (int p = 0; p < N; ++p)
    {
        if (parts[p].move_)
        {
            t = 0.5 * (parts[p].qm_) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = parts[p].velocity_ + 0.5 * (parts[p].qm_) * E[index] * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            parts[p].velocity_ = v_plus + 0.5 * (parts[p].qm_) * E[index] * dt;
        }
        index++;
    }
}
