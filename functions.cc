#include "initial.h"

typedef size_t Index;
typedef double Real;
typedef std::valarray<Real> Array;
typedef std::vector<Array> VecVal;
typedef std::vector<VecVal> VecVecVal;
typedef std::complex<Real> Complex;
typedef std::valarray<Complex> CArray;
typedef std::vector<std::vector<Real>> VecVec;

VecVal valarraysVector(const Index& rows, const Index& cols, const Real& value) {
    VecVal f;
    Array zeros(value, cols);
    for (Index i = 0; i < rows; ++i) {
        f.push_back(zeros);
    }
    return f;
}

// INITIALIZATION OF THE SYSTEM

std::vector<Particle> two_stream() {
    const int N = params::N;
    const Real vt = params::vt;
    const Real vd = params::vd;
    const int gridPoints = params::gp;
    const Real dr = params::dr;
    const Real L = dr * Real(gridPoints - 1);
    const Real n = N / (L * L);
    const Real margin = dr / 10.0;
    std::mt19937_64 engine;
    std::normal_distribution<Real> vel_left(vd, vt);
    std::normal_distribution<Real> vel_right(-vd, vt);
    VecVal pos;
    std::vector<Real> neutro;
    std::vector<Real> right;
    std::vector<Real> left;
    std::vector<Particle> parts;

    Real deltaL = (L - 2.0 * margin) / (std::sqrt(N) - 1);

    for (Real x = margin; x <= L; x += deltaL) {
        for (Real y = margin; y <= L; y += deltaL) {
            pos.push_back({x, y, 0.0});
        }
    }

    // for (auto p : pos) {
    //     for (auto i : p) {
    //         std::cout << i << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // for (Index i = 0; i < N; ++i) {
    //     indexes.push_back(i);
    // }

    // std::random_shuffle(indexes.begin(), indexes.end());
    Index beam = Index(N / 4);

    // for (Index i = N-1; i >= N-beam; --i) {
    //     right.push_back(indexes[i]);
    //     indexes.pop_back();
    // }

    // for (Index i = N-beam-1; i >= Index(N/2); --i) {
    //     left.push_back(indexes[i]);
    //     indexes.pop_back();
    // }

    for (Index i = 0; i < Index(N / 2); i++) {
        neutro.push_back(i);
    }

    for (Index i = Index(N / 2); i < N-beam; i++) {
        right.push_back(i);
    }

    for (Index i = N-beam; i < N; i++) {
        left.push_back(i);
    }

    for (auto i : right) {
        parts.push_back(Particle(pos[i], {0.1, 0.0, 0.0}, n, -1.0, true));
    }

    for (auto i : left) {
        parts.push_back(Particle(pos[i], {-0.1, 0.0, 0.0}, n, -1.0, true));
    }

    for (auto i : neutro) {
        parts.push_back(Particle(pos[i], {0.0, 0.0, 0.0}, n, 1.0, false));
    }

    return parts;
}

std::vector<Particle> random_particles() {
    const Index N = params::N;
    const Real vd = params::vd;
    const Index gridPoints = params::gp;
    const Real dr = params::dr;
    const Real L = dr * Real(gridPoints - 1);
    const Real n = N / (L * L);
    std::mt19937_64 engine;
    std::uniform_real_distribution<Real> vel(-vd, vd);
    std::uniform_real_distribution<Real> posi(0, L);
    VecVal pos;
    std::vector<Particle> parts;

    for (Index i = 0; i < N; ++i) {
        if (i < (N / 2)) {
            parts.push_back(Particle({posi(engine), posi(engine), 0.0}, {0.0, 0.0, 0.0}, n, 1.0, true));
        } else {
            parts.push_back(Particle({posi(engine), posi(engine), 0.0}, {0.0, 0.0, 0.0}, n, -1.0, true));
        }
    }

    return parts;
}

VecVal density(std::vector<Particle>& particles, const std::vector<Real>& rho_c) {
    const Index gp = params::gp;
    const Index N = params::N;
    const Real dr = params::dr;
    VecVal rho = valarraysVector(gp, gp, 0.0);

    for (Index p = 0; p < N; ++p) {
        Index i = std::floor(particles[p].position_[0] / dr);
        Index j = std::floor(particles[p].position_[1] / dr);
        Real hx = particles[p].position_[0] - (i * dr);
        Real hy = particles[p].position_[1] - (j * dr);

        rho[i][j] += rho_c[p] * (dr - hx) * (dr - hy);
        rho[i][j+1] += rho_c[p] * (dr - hx) * hy;
        rho[i+1][j] += rho_c[p] * hx * (dr - hy);
        rho[i+1][j+1] += rho_c[p] * hx * hy;
    }

    for (Index u = 0; u < gp; ++u) {
        rho[gp-1][u] = (rho[gp-1][u] + rho[0][u]) * 0.5;
        rho[0][u] = rho[gp-1][u];
    }

    for (Index u = 0; u < gp; ++u) {
        rho[u][gp-1] = (rho[u][gp-1] + rho[u][0]) * 0.5;
        rho[u][0] = rho[u][gp-1];
    }

    for (Index i = 0; i < gp; ++i) {
        rho[i] /= (dr * dr);
    }

    return rho;
}

VecVal potential(const VecVal& rho) {   
    const Index gp = params::gp;
    const Real dr = params::dr;
    Complex rho_k[gp][gp];

    for (Index i = 0; i < gp; ++i) {
           for (Index j = 0; j < gp; ++j) {
               rho_k[i][j] = rho[i][j];
           }
       }   

    CArray f(gp);

    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            f[j] = rho_k[i][j];
        }

        fft(f);
        for (Index j = 0; j < gp; ++j) {
            rho_k[i][j] = f[j];
        }
    }

    for (Index j = 0; j < gp; ++j)
    {
        for (Index i = 0; i < gp; ++i) {    
            f[i] = rho_k[i][j];
        }
        
        fft(f);
        for (Index i = 0; i < gp; ++i) {
            rho_k[i][j] = f[i];
        }
    }

    Complex phi_k[gp][gp];

    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            phi_k[i][j] = rho_k[i][j];
        }
    }

    Complex i(0.0, 1.0);
    Complex W = std::exp(2.0 * constants::pi * i / Real(gp));
    Complex Wm = 1.0;
    Complex Wn = 1.0;

    for (Index m = 0; m < gp; ++m) {
        for (Index n = 0; n < gp; ++n) {
            Complex denom = 4.0;
            denom -= Wm + 1.0/Wm + Wn + 1.0/Wn;
            if (denom != 0.0) {
                phi_k[m][n] *= (dr * dr) / denom;
            }
            Wn *= W;
        }
        Wm *= W;
    }

    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            f[j] = phi_k[i][j];
        }

        ifft(f);
        for (Index j = 0; j < gp; ++j) {
            phi_k[i][j] = f[j];
        }
    }

    for (Index j = 0; j < gp; ++j) {
        for (Index i = 0; i < gp; ++i) {
            f[i] = phi_k[i][j];
        }
        
        ifft(f);
        for (Index i = 0; i < gp; ++i) {
            phi_k[i][j] = f[i];
        }
    }

    VecVal phi = valarraysVector(gp, gp, 0.0);
    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            phi[i][j] = std::real(phi_k[i][j]);
        }
    }

    return phi;
}

VecVecVal EField_GP(const VecVal& phi) {
    const Index gp = params::gp;
    const Real dr = params::dr;
    VecVecVal E;
    Array zeros(0.0, 3);
    for (Index i = 0; i < gp; ++i) {   
        VecVal rows;
        for (Index j = 0; j < gp; ++j) {
            rows.push_back(zeros);
        }
        E.push_back(rows);  
    }

    for (Index j = 0; j < gp; ++j) {
        for (Index i = 0; i < gp; ++i) {
            Index nxt_i = (i < gp - 1) ? (i + 1) : 0;
            Index prv_i = (i > 0) ? (i - 1) : (gp - 1);

            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (2 * dr);
        }
    }

    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            Index nxt_j = (j < gp - 1) ? (j + 1) : 0;
            Index prv_j = (j > 0) ? (j - 1) : (gp - 1);

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (2 * dr);
        }
    }

    return E;
}

VecVal EField_P(const VecVecVal& field, const std::vector<Particle>& particles) {
    const Index N = params::N;
    const Real dr = params::dr;
    VecVal E = valarraysVector(N, 3, 0.0);

    for (Index p = 0; p < N; ++p) {
        if (particles[p].move_) {
            Index i = std::floor(particles[p].position_[0] / dr);
            Index j = std::floor(particles[p].position_[1] / dr);
            Real hx = particles[p].position_[0] - (i * dr);
            Real hy = particles[p].position_[1] - (j * dr);

            Real A = (dr - hx) * (dr - hy);
            Real B = (dr - hx) * hy;
            Real C = hx * (dr - hy);
            Real D = hx * hy;

            E[p][0] = field[i][j][0] * A + field[i][j+1][0] * B + field[i+1][j][0] * C + field[i+1][j+1][0] * D;
            E[p][1] = field[i][j][1] * A + field[i][j+1][1] * B + field[i+1][j][1] * C + field[i+1][j+1][1] * D;
        }
    }

    for (Index i = 0; i < N; ++i) {
        E[i] /= (dr * dr);
    }

    return E;
}

Real norm(const Array& A) {
    return std::sqrt((A * A).sum());
}

inline Array cross(const Array& A, const Array& B) {
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

inline Real mod(const Real& a, const Real& b) {
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}

void Boris(const VecVal& E, const Array& extE, const Array& B, std::vector<Particle>& particles) {
    const Real L = params::dr * Real(params::gp - 1);
    const Real dt = params::dt;
    const Index N = params::N;
    Array t;
    Real t_2;
    Array s;
    Array v_minus;
    Array v_prime;
    Array v_plus;

    for (Index p = 0; p < N; ++p) {
        if (particles[p].move_) {
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

void outPhase(const Real& direction, const VecVal& E, const Array& extE, const Array& B, std::vector<Particle>& particles) {
    const Real dt = direction * 0.5 * params::dt;
    const Index N = params::N;
    Array t;
    Real t_2;
    Array s;
    Array v_minus;
    Array v_prime;
    Array v_plus;

    for (Index p = 0; p < N; ++p) {
        if (particles[p].move_) {
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


void writeData(std::string filename, std::string sname, const VecVec& data) {
    hid_t   file_id, dataset_id, dataspace_id, memspace_id, dcpl_id;
    hsize_t dims[2], dimsm[1], chunk_dims[2];
    hsize_t offset[2];
    hsize_t count[2];
    hsize_t stride[2];
    hsize_t block[2];
    Index   NX = data.size();
    Index   NY = data.at(0).size();
    Index   szip_options_mask, pixels_per_block;

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    dims[0] = NX;
    dims[1] = NY;

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    szip_options_mask = H5_SZIP_NN_OPTION_MASK;
    pixels_per_block = 10;
    H5Pset_szip(dcpl_id, szip_options_mask, pixels_per_block);  

    if (sname == "mesh" || sname == "energy") {
        chunk_dims[0] = Index(NX / 1);
        chunk_dims[1] = 1;
    } else {
        chunk_dims[0] = 1;
        chunk_dims[1] = Index(NY / 1);
    }

    H5Pset_chunk(dcpl_id, 2, chunk_dims);
    dataset_id = H5Dcreate2(file_id, sname.c_str(), H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    dimsm[0] = NY;
    memspace_id = H5Screate_simple(1, dimsm, NULL);

    count[0] = 1;
    count[1] = NY;

    stride[0] = 1;
    stride[1] = 1;

    block[0] = 1;
    block[1] = 1;

    for (Index i = 0; i < NX; i++) {
        offset[0] = i;
        offset[1] = 0;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, data.at(i).data());
    }

    // dataset_id = H5Dcreate2(file_id, sname.c_str(), H5T_IEEE_F64LE, dataspace_id,
    //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // for (Index i = 0; i < NX; i++) {
    //     H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.at(i).data());
    // }    

    std::cout << sname << " written" << std::endl;
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(dcpl_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
}