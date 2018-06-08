#include "initial.h"

typedef size_t Index;
typedef double Real;
typedef std::vector<std::valarray<Real>> VecVal;
typedef std::vector<std::vector<std::valarray<Real>>> VecVecVal;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;

VecVal valarraysVector(const Index& rows, const Index& cols)
{
    VecVal f;
    std::valarray<Real> zeros(0.0, cols);
    for (Index i = 0; i < rows; ++i)
    {
        f.push_back(zeros);
    }

    return f;
}

// INITIALIZATION OF THE SYSTEM

std::vector<Particle> two_stream()
{
    const int N = parameters::N;
    const Real vt = parameters::vt;
    const Real vd = parameters::vd;
    const int gridPoints = parameters::gp;
    const Real dr = parameters::dr;
    const Real L = dr * Real(gridPoints - 1);
    const Real n = N / (L * L);
    const Real margin = dr / 10.0;
    std::mt19937_64 engine;
    std::normal_distribution<Real> vel_left(vd, vt);
    std::normal_distribution<Real> vel_right(-vd, vt);
    std::vector<std::valarray<Real>> pos;
    std::vector<Real> indexes;
    std::vector<Real> right;
    std::vector<Real> left;
    std::vector<Particle> parts;

    Real deltaL = (L - 2.0 * margin) / (std::sqrt(N) - 1);

    for (Real x = margin; x <= L; x += deltaL)
    {
        for (Real y = margin; y <= L; y += deltaL)
        {
            pos.push_back({x, y, 0.0});
        }
    }

    for (Index i = 0; i < N; ++i)
    {
        indexes.push_back(i);
    }

    std::random_shuffle(indexes.begin(), indexes.end());
    Index electronBeam = Index(N / 4);

    for (Index i = N-1; i >= N-electronBeam; --i)
    {
        right.push_back(indexes[i]);
        indexes.pop_back();
    }

    for (Index i = N-electronBeam-1; i >= Index(N/2); --i)
    {
        left.push_back(indexes[i]);
        indexes.pop_back();
    }

    for (auto i = left.begin(); i != left.end(); ++i)
    {
        parts.push_back(Particle(pos[*i], {vel_left(engine), 0.0, 0.0}, n, -1.0, true));
    }

    for (auto i = right.begin(); i != right.end(); ++i)
    {
        parts.push_back(Particle(pos[*i], {vel_right(engine), 0.0, 0.0}, n, -1.0, true));
    }

    for (auto i = indexes.begin(); i != indexes.end(); ++i)
    {
        parts.push_back(Particle(pos[*i], {0.0, 0.0, 0.0}, n, 1.0, false));
    }

    return parts;
}

std::vector<Particle> random_particles()
{
    const Index N = parameters::N;
    const Real vd = parameters::vd;
    const Index gridPoints = parameters::gp;
    const Real dr = parameters::dr;
    const Real L = dr * Real(gridPoints - 1);
    const Real n = N / (L * L);
    std::mt19937_64 engine;
    std::uniform_real_distribution<Real> vel(-vd, vd);
    std::uniform_real_distribution<Real> posi(0, L);
    std::vector<std::valarray<Real>> pos;
    std::vector<Particle> parts;

    for (Index i = 0; i < N; ++i)
    {
        if (i < (N / 2))
        {
            parts.push_back(Particle({posi(engine), posi(engine), 0.0}, {0.0, 0.0, 0.0}, n, 1.0, true));
        }

        else
        {
            parts.push_back(Particle({posi(engine), posi(engine), 0.0}, {0.0, 0.0, 0.0}, n, -1.0, true));
        }
    }

    return parts;
}

VecVal density(std::vector<Particle>& particles, const std::vector<Real>& rho_c, const Real& dr)
{
    const Index gp = parameters::gp;
    const Index N = parameters::N;
    VecVal rho = valarraysVector(gp, gp);

    for (Index p = 0; p < N; ++p)
    {
        Index i = std::floor(particles[p].position_[0] / dr);
        Index j = std::floor(particles[p].position_[1] / dr);
        Real hx = particles[p].position_[0] - (i * dr);
        Real hy = particles[p].position_[1] - (j * dr);

        rho[i][j] += rho_c[p] * (dr - hx) * (dr - hy);
        rho[i][j+1] += rho_c[p] * (dr - hx) * hy;
        rho[i+1][j] += rho_c[p] * hx * (dr - hy);
        rho[i+1][j+1] += rho_c[p] * hx * hy;
    }

    for (Index u = 0; u < gp; ++u)
    {
        rho[gp - 1][u] = (rho[gp - 1][u] + rho[0][u]) * 0.5;
        rho[0][u] = rho[gp - 1][u];
    }

    for (Index u = 0; u < gp; ++u)
    {
        rho[u][gp - 1] = (rho[u][gp - 1] + rho[u][0]) * 0.5;
        rho[u][0] = rho[u][gp - 1];
    }

    for (Index i = 0; i < gp; ++i)
    {
        rho[i] /= (dr * dr);
    }

    return rho;
}

VecVal potential(const VecVal& rho, const Real& dr)
{   
    const Index gp = parameters::gp;
    Complex rho_k[gp][gp];

    for (Index i = 0; i < gp; ++i)
       {
           for (Index j = 0; j < gp; ++j)
           {
               rho_k[i][j] = rho[i][j];
           }
       }   

    CArray f(gp);

    for (Index i = 0; i < gp; ++i)
    {
        for (Index j = 0; j < gp; ++j)
        {
            f[j] = rho_k[i][j];
        }
        fft(f);
        for (Index j = 0; j < gp; ++j)
        {
            rho_k[i][j] = f[j];
        }
    }

    for (Index j = 0; j < gp; ++j)
    {
        for (Index i = 0; i < gp; ++i)
        {    
            f[i] = rho_k[i][j];
        }
        fft(f);
        for (Index i = 0; i < gp; ++i)
        {
            rho_k[i][j] = f[i];
        }
    }

    Complex phi_k[gp][gp];

    for (Index i = 0; i < gp; ++i)
    {
        for (Index j = 0; j < gp; ++j)
        {
            phi_k[i][j] = rho_k[i][j];
        }
    }

    Complex i(0.0, 1.0);
    Complex W = std::exp(2.0 * constants::pi * i / Real(gp));
    Complex Wm = 1.0;
    Complex Wn = 1.0;

    for (Index m = 0; m < gp; ++m)
    {
        for (Index n = 0; n < gp; ++n)
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

    for (Index i = 0; i < gp; ++i)
    {
        for (Index j = 0; j < gp; ++j)
        {
            f[j] = phi_k[i][j];
        }
        ifft(f);
        for (Index j = 0; j < gp; ++j)
        {
            phi_k[i][j] = f[j];
        }
    }

    for (Index j = 0; j < gp; ++j)
    {
        for (Index i = 0; i < gp; ++i)
        {
            f[i] = phi_k[i][j];
        }
        ifft(f);
        for (Index i = 0; i < gp; ++i)
        {
            phi_k[i][j] = f[i];
        }
    }

    VecVal phi = valarraysVector(gp, gp);
    for (Index i = 0; i < gp; ++i)
    {
        for (Index j = 0; j < gp; ++j)
        {
            phi[i][j] = std::real(phi_k[i][j]);
        }
    }

    return phi;
}

VecVecVal EField_GP(const VecVal& phi, const Real& dr)
{
    const Index gp = parameters::gp;
    VecVecVal E;
    std::valarray<Real> zeros(0.0, 3);
    for (Index i = 0; i < gp; ++i)
    {   
        VecVal rows;
        for (Index j = 0; j < gp; ++j)
        {
            rows.push_back(zeros);
        }
        E.push_back(rows);  
    }

    for (Index j = 0; j < gp; ++j)
    {
        for (Index i = 0; i < gp; ++i)
        {
            Index nxt_i = (i < gp - 1) ? (i + 1) : 0;
            Index prv_i = (i > 0) ? (i - 1) : (gp - 1);

            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (2 * dr);
        }
    }

    for (Index i = 0; i < gp; ++i)
    {
        for (Index j = 0; j < gp; ++j)
        {
            Index nxt_j = (j < gp - 1) ? (j + 1) : 0;
            Index prv_j = (j > 0) ? (j - 1) : (gp - 1);

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (2 * dr);
        }
    }

    return E;
}

VecVal EField_P(const VecVecVal& field, const std::vector<Particle>& particles, const Real& dr)
{
    const Index N = parameters::N;
    VecVal E = valarraysVector(N, 3);

    Index index = 0;

    for (Index p = 0; p < N; ++p)
    {
        if (particles[p].move_)
        {
            Index i = std::floor(particles[p].position_[0] / dr);
            Index j = std::floor(particles[p].position_[1] / dr);
            Real hx = particles[p].position_[0] - (i * dr);
            Real hy = particles[p].position_[1] - (j * dr);

            Real A = (dr - hx) * (dr - hy);
            Real B = (dr - hx) * hy;
            Real C = hx * (dr - hy);
            Real D = hx * hy;

            E[index][0] = field[i][j][0] * A + field[i][j+1][0] * B + field[i+1][j][0] * C + field[i+1][j+1][0] * D;
            E[index][1] = field[i][j][1] * A + field[i][j+1][1] * B + field[i+1][j][1] * C + field[i+1][j+1][1] * D;
        }
        index++;
    }

    for (Index i = 0; i < N; ++i)
    {
        E[i] /= (dr * dr);
    }

    return E;
}

Real norm(const std::valarray<Real>& Array)
{
    return std::sqrt((Array * Array).sum());
}

inline std::valarray<Real> cross(const std::valarray<Real>& A, const std::valarray<Real>& B)
{
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

inline Real mod(const Real& a, const Real& b)
{
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}

void Boris(const VecVal& E, const std::valarray<Real>& extE, const std::valarray<Real>& B, std::vector<Particle>& particles)
{
    const Real L = parameters::dr * Real(parameters::gp - 1);
    const Real dt = parameters::dt;
    const Index N = parameters::N;
    std::valarray<Real> t;
    Real t_2;
    std::valarray<Real> s;
    std::valarray<Real> v_minus;
    std::valarray<Real> v_prime;
    std::valarray<Real> v_plus;

    for (Index p = 0; p < N; ++p)
    {
        if (particles[p].move_)
        {
            t = 0.5 * (particles[p].qm_) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = particles[p].velocity_ + 0.5 * (particles[p].qm_) * (E[p] + extE) * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            particles[p].velocity_ = v_plus + 0.5 * (particles[p].qm_) * (E[p] + extE) * dt;

            particles[p].position_ += particles[p].velocity_ * dt;
            particles[p].position_[0] = mod(particles[p].position_[0], L);
            particles[p].position_[1] = mod(particles[p].position_[1], L);
            particles[p].position_[2] = mod(particles[p].position_[2], L);
        }
    }
}

void outPhase(const Real& direction, const VecVal& E, const std::valarray<Real>& extE, const std::valarray<Real>& B, std::vector<Particle>& particles)
{
    const Real dt = direction * 0.5 * parameters::dt;
    const Index N = parameters::N;
    std::valarray<Real> t;
    Real t_2;
    std::valarray<Real> s;
    std::valarray<Real> v_minus;
    std::valarray<Real> v_prime;
    std::valarray<Real> v_plus;

    for (Index p = 0; p < N; ++p)
    {
        if (particles[p].move_)
        {
            t = 0.5 * (particles[p].qm_) * B * dt;
            t_2 = norm(t) * norm(t);
            s = (2.0 * t) / (1.0 + t_2);
            v_minus = particles[p].velocity_ + 0.5 * (particles[p].qm_) * (E[p] + extE) * dt;
            v_prime = v_minus + cross(v_minus, t);
            v_plus = v_minus + cross(v_prime, s);
            particles[p].velocity_ = v_plus + 0.5 * (particles[p].qm_) * (E[p] + extE) * dt;
        }
    }
}


/**********************************************************************************************************/

/*                                       DATA WRITING FUNTION                                            */

/**********************************************************************************************************/


void writeData(std::string filename, std::string sname, const VecVal& data)
{
    hid_t   file_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    Index   RANK = 2;
    Real    sdata[data.size()][data.at(0).size()];

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    dims[0] = data.size();
    dims[1] = data.at(0).size();

    dataspace_id = H5Screate_simple(RANK, dims, NULL);
    dataset_id = H5Dcreate2(file_id, sname.c_str(), H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (Index i = 0; i < data.size(); i++)
    {
        for (Index j = 0; j < data.at(0).size(); j++)
            sdata[i][j] = data[i][j];
    }

    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
}
