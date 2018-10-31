#include "functions.h"

typedef std::vector<std::valarray<double>> VecVal;
typedef std::vector<std::vector<std::valarray<double>>> VecVecVal;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

double sign(double& x)
{
    if (x < 0) return -1.0;
    else if (x > 0) return 1.0;
    else return 0.0;
}

double norm(const std::valarray<double>& Array)
{
    return std::sqrt((Array * Array).sum());
}

inline std::valarray<double> cross(const std::valarray<double>& A, const std::valarray<double>& B)
{
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

inline double mod(const double& a, const double& b)
{
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}

VecVal valarraysVector(const int& rows, const int& cols)
{
    VecVal f;
    std::valarray<double> zeros(0.0, cols);
    for (int i = 0; i < rows; ++i)
    {
        f.push_back(zeros);
    }

    return f;
}


VecVal density(const std::vector<std::valarray<double>>& positions, 
               const std::vector<double>& charges, const double& dr, 
               const int& gp, const int& N)
{
    VecVal rho = valarraysVector(gp, gp);

    for (int p = 0; p < N; ++p)
    {
        int i = std::floor(positions.at(p)[0] / dr);
        int j = std::floor(positions.at(p)[1] / dr);
        double hx = positions.at(p)[0] - (i * dr);
        double hy = positions.at(p)[1] - (j * dr);
        int nxt_i = mod(i + 1, gp);
        int nxt_j = mod(j + 1, gp);

        rho[i][j] += charges.at(p) * (dr - hx) * (dr - hy);
        rho[i][nxt_j] += charges.at(p) * (dr - hx) * hy;
        rho[nxt_i][j] += charges.at(p) * hx * (dr - hy);
        rho[nxt_i][nxt_j] += charges.at(p) * hx * hy;
    }

    for (int i = 0; i < gp; ++i)
    {
        rho[i] /= (dr * dr * dr * dr);
    }

    return rho;
}


VecVal potential(const VecVal& rho, const double& dr, const int& gp)
{   
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
    Complex W = std::exp(2.0 * std::acos(-1.0) * i / double(gp));
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


VecVecVal EField_GP(const VecVal& phi, const double& dr, const int& gp)
{
    VecVecVal E;
    std::valarray<double> zeros(0.0, 3);
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
            int nxt_i = mod(i + 1, gp);
            int prv_i = mod(i - 1, gp);
            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (2 * dr);
        }
    }

    for (int i = 0; i < gp; ++i)
    {
        for (int j = 0; j < gp; ++j)
        {
            int nxt_j = mod(j + 1, gp);
            int prv_j = mod(j - 1, gp);

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (2 * dr);
        }
    }

    return E;
}


VecVal EField_P(const VecVecVal& field, 
                const std::vector<std::valarray<double>>& positions,
                const std::vector<int>& moves, 
                const double& dr, const int& gp, const int& N)
{
    VecVal E = valarraysVector(N, 3);

    int index = 0;
    for (int p = 0; p < N; ++p)
    {
        if (moves.at(p))
        {
            int i = std::floor(positions.at(p)[0] / dr);
            int j = std::floor(positions.at(p)[1] / dr);
            double hx = positions.at(p)[0] - (i * dr);
            double hy = positions.at(p)[1] - (j * dr);
            int nxt_i = mod(i + 1, gp);
            int nxt_j = mod(j + 1, gp);

            double A = (dr - hx) * (dr - hy);
            double B = (dr - hx) * hy;
            double C = hx * (dr - hy);
            double D = hx * hy;

            E[index][0] = field[i][j][0] * A + field[i][nxt_j][0] * B + field[nxt_i][j][0] * C + field[nxt_i][nxt_j][0] * D;
            E[index][1] = field[i][j][1] * A + field[i][nxt_j][1] * B + field[nxt_i][j][1] * C + field[nxt_i][nxt_j][1] * D;
        }
        index++;
    }

    for (int i = 0; i < N; ++i)
    {
        E[i] /= (dr * dr);
    }

    return E;
}


void Boris(std::vector<std::valarray<double>>& positions, 
           std::vector<std::valarray<double>>& velocities,
           const std::vector<double>& QoverM,
           const std::vector<int>& moves,
           const VecVal& E, const std::valarray<double>& extE, 
           const std::valarray<double>& B, 
           const double& L, const double& dt, const int& N)
{
    int index = 0;
    std::valarray<double> t;
    double t_2;
    std::valarray<double> s;
    std::valarray<double> v_minus;
    std::valarray<double> v_prime;
    std::valarray<double> v_plus;

    for (int p = 0; p < N; ++p)
    {
        if (moves.at(p))
        {
            t = 0.5 * QoverM.at(p) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = velocities.at(p) + 0.5 * QoverM.at(p) * (E[index] + extE) * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            velocities.at(p) = v_plus + 0.5 * QoverM.at(p) * (E[index] + extE) * dt;

            positions.at(p) += velocities.at(p) * dt;
            positions.at(p)[0] = mod(positions.at(p)[0], L);
            positions.at(p)[1] = mod(positions.at(p)[1], L);
        }
        index++;
    }
}


void outphase(std::vector<std::valarray<double>>& velocities,
              const std::vector<double>& QoverM,
              const std::vector<int>& moves,
              const double& direction, const VecVal& E, 
              const std::valarray<double>& extE, 
              const std::valarray<double>& B, 
              const double& dt, const int& N)
{
    int index = 0;
    const double dT = direction * 0.5 * dt;
    std::valarray<double> t;
    double t_2;
    std::valarray<double> s;
    std::valarray<double> v_minus;
    std::valarray<double> v_prime;
    std::valarray<double> v_plus;

    for (int p = 0; p < N; ++p)
    {
        if (moves.at(p))
        {
            t = 0.5 * QoverM.at(p) * B * dT;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = velocities.at(p) + 0.5 * QoverM.at(p) * (E[index] + extE) * dT;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            velocities.at(p) = v_plus + 0.5 * QoverM.at(p) * (E[index] + extE) * dT;
        }
        index++;
    }
}
