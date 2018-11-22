#include "functions.h"

Real sign(Real& x)
{
    if (x < 0) return -1.0;
    else if (x > 0) return 1.0;
    else return 0.0;
}


Real norm(const Array& A)
{
    return std::sqrt((A * A).sum());
}


inline Array cross(const Array& A, const Array& B)
{
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}


Real mod(const Real& a, const Real& b)
{
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}


VecArr valarraysVector(const Int& rows, const Int& cols)
{
    VecArr f;
    Array zeros(0.0, cols);
    for (Int i = 0; i < rows; ++i)
    {
        f.push_back(zeros);
    }

    return f;
}


VecArr density(const VecArr& positions, 
               const std::vector<Real>& charges, const Real& dx, 
               const Real& dy, const Int& Nx, const Int& Ny, const Int& N)
{
    VecArr rho = valarraysVector(Nx, Ny);

    for (Int p = 0; p < N; ++p)
    {
        Int i = std::floor(positions.at(p)[0] / dx);
        Int j = std::floor(positions.at(p)[1] / dy);
        Real hx = positions.at(p)[0] - (i * dx);
        Real hy = positions.at(p)[1] - (j * dy);
        Int nxt_i = mod(i + 1, Nx);
        Int nxt_j = mod(j + 1, Ny);

        rho[i][j] += charges.at(p) * (dx - hx) * (dy - hy);
        rho[i][nxt_j] += charges.at(p) * (dx - hx) * hy;
        rho[nxt_i][j] += charges.at(p) * hx * (dy - hy);
        rho[nxt_i][nxt_j] += charges.at(p) * hx * hy;
    }

    for (Int i = 0; i < Nx; ++i)
    {
        rho[i] /= (dx * dx * dy * dy);
    }

    return rho;
}


VecArr potential(const VecArr& rho, const Real& dx, const Real& dy, 
                 const Int& Nx, const Int& Ny)
{   
    Complex rho_k[Nx][Ny];

    for (Int i = 0; i < Nx; ++i)
       {
           for (Int j = 0; j < Ny; ++j)
           {
               rho_k[i][j] = rho[i][j];
           }
       }   

    CArray f(Nx);
    CArray g(Ny);

    for (Int i = 0; i < Nx; ++i)
    {
        for (Int j = 0; j < Ny; ++j)
        {
            g[j] = rho_k[i][j];
        }
        fft(g);
        for (Int j = 0; j < Ny; ++j)
        {
            rho_k[i][j] = g[j];
        }
    }

    for (Int j = 0; j < Ny; ++j)
    {
        for (Int i = 0; i < Nx; ++i)
        {    
            f[i] = rho_k[i][j];
        }
        fft(f);
        for (Int i = 0; i < Nx; ++i)
        {
            rho_k[i][j] = f[i];
        }
    }

    Complex phi_k[Nx][Ny];

    for (Int i = 0; i < Nx; ++i)
    {
        for (Int j = 0; j < Ny; ++j)
        {
            phi_k[i][j] = rho_k[i][j];
        }
    }

    Complex i(0.0, 1.0);
    Complex Wx = std::exp(2.0 * std::acos(-1.0) * i / Real(Nx));
    Complex Wy = std::exp(2.0 * std::acos(-1.0) * i / Real(Ny));
    Complex Wn = 1.0;
    Complex Wm = 1.0;
    Real dx_2 = dx * dx;
    Real dy_2 = dy * dy;

    for (Int n = 0; n < Nx; ++n)
    {
        for (Int m = 0; m < Ny; ++m)
        {
            Complex denom = dy_2 * (2.0 - Wn - 1.0/Wn) + dx_2 * (2.0 - Wm - 1.0/Wm);
            if (denom != 0.0)
            {
                phi_k[n][m] *= (dx_2 * dy_2) / denom;
            }
            Wm *= Wy;
        }
        Wn *= Wx;
    }

    for (Int i = 0; i < Nx; ++i)
    {
        for (Int j = 0; j < Ny; ++j)
        {
            g[j] = phi_k[i][j];
        }
        ifft(g);
        for (Int j = 0; j < Ny; ++j)
        {
            phi_k[i][j] = g[j];
        }
    }

    for (Int j = 0; j < Ny; ++j)
    {
        for (Int i = 0; i < Nx; ++i)
        {
            f[i] = phi_k[i][j];
        }
        ifft(f);
        for (Int i = 0; i < Nx; ++i)
        {
            phi_k[i][j] = f[i];
        }
    }

    VecArr phi = valarraysVector(Nx, Ny);
    for (Int i = 0; i < Nx; ++i)
    {
        for (Int j = 0; j < Ny; ++j)
        {
            phi[i][j] = std::real(phi_k[i][j]);
        }
    }

    return phi;
}


VecVecArr fieldNodes(const VecArr& phi, const Real& dx, const Real& dy, 
                     const Int& Nx, const Int& Ny)
{
    VecVecArr E;
    Array zeros(0.0, 3);
    for (Int i = 0; i < Nx; ++i)
    {   
        VecArr rows;
        for (Int j = 0; j < Ny; ++j)
        {
            rows.push_back(zeros);
        }
        E.push_back(rows);  
    }

    for (Int j = 0; j < Ny; ++j)
    {
        for (Int i = 0; i < Nx; ++i)
        {
            Int nxt_i = mod(i + 1, Nx);
            Int prv_i = mod(i - 1, Nx);
            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (dx * 2);
        }
    }

    for (Int i = 0; i < Nx; ++i)
    {
        for (Int j = 0; j < Ny; ++j)
        {
            Int nxt_j = mod(j + 1, Ny);
            Int prv_j = mod(j - 1, Ny);

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (dy * 2);
        }
    }

    return E;
}


VecArr fieldParticles(const VecVecArr& field, 
                      const VecArr& positions,
                      const std::vector<Int>& moves, 
                      const Real& dx, const Real& dy, 
                      const Int& Nx, const Int& Ny, const Int& N)
{
    VecArr E = valarraysVector(N, 3);

    for (Int p = 0; p < N; ++p)
    {
        if (moves.at(p))
        {
            Int i = std::floor(positions.at(p)[0] / dx);
            Int j = std::floor(positions.at(p)[1] / dy);
            Real hx = positions.at(p)[0] - (i * dx);
            Real hy = positions.at(p)[1] - (j * dy);
            Int nxt_i = mod(i + 1, Nx);
            Int nxt_j = mod(j + 1, Ny);

            Real A = (dx - hx) * (dy - hy);
            Real B = (dx - hx) * hy;
            Real C = hx * (dy - hy);
            Real D = hx * hy;

            E[p][0] = field[i][j][0] * A + field[i][nxt_j][0] * B + field[nxt_i][j][0] * C + field[nxt_i][nxt_j][0] * D;
            E[p][1] = field[i][j][1] * A + field[i][nxt_j][1] * B + field[nxt_i][j][1] * C + field[nxt_i][nxt_j][1] * D;
        }
    }

    for (Int i = 0; i < N; ++i)
    {
        E[i] /= (dx * dy);
    }

    return E;
}


void boris(VecArr& velocities,
           const std::vector<Real>& QoverM,
           const std::vector<Int>& moves,
           const VecArr& E, const Array& B, 
           const Real& dt, const Int& N)
{
    Array t;
    Real t_2;
    Array s;
    Array v_minus;
    Array v_prime;
    Array v_plus;

    for (Int p = 0; p < N; ++p)
    {
        if (moves.at(p))
        {
            t = 0.5 * QoverM.at(p) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = velocities.at(p) + 0.5 * QoverM.at(p) * E[p] * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            velocities.at(p) = v_plus + 0.5 * QoverM.at(p) * E[p] * dt;
        }
    }
}


void update(VecArr& positions,
           VecArr& velocities,
           const std::vector<Real>& QoverM,
           const std::vector<Int>& moves,
           const VecArr& E, const Array& B, 
           const Real& Lx, const Real& Ly, 
           const Real& dt, const Int& N)
{
    boris(velocities, QoverM, moves, E, B, dt, N);
    for (Int p = 0; p < N; ++p)
    {
        positions.at(p) += velocities.at(p) * dt;
        positions.at(p)[0] = mod(positions.at(p)[0], Lx);
        positions.at(p)[1] = mod(positions.at(p)[1], Ly);
    }

}


void outphase(const Real& direction,
              VecArr& velocities,
              const std::vector<Real>& QoverM,
              const std::vector<Int>& moves,
              const VecArr& E, const Array& B, 
              const Real& dt, const Int& N)
{
    Real dT = 0.5 * direction * dt;
    boris(velocities, QoverM, moves, E, B, dT, N);
}
