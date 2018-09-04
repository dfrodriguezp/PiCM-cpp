#include "particle.h"

const Real pi = 3.14159265358979323846264338328;

VecArr arraysVector(const Index& rows, const Index& cols, const Real& value) {
    VecArr f;
    Array A(value, cols);
    for (Index i = 0; i < rows; ++i) {
        f.push_back(A);
    }
    return f;
}

VecVecArr arraysVectorVector(const Index& rows, const Index& cols, const Real& value) {
    VecVecArr f;
    Array A(value, 3);
    for (Index i = 0; i < rows; ++i) {
        VecArr lines;
        for (Index j = 0; j < cols; ++j) {
            lines.push_back(A);
        }
        f.push_back(lines);
    }
    return f;
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

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
void fft(CArray& x) {
    // DFT
    Index N = x.size(), k = N, n;
    Real thetaT = 3.14159265358979323846264338328L / N;
    Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (Index l = 0; l < k; l++) {
            for (Index a = l; a < N; a += n) {
                Index b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    Index m = (Index)log2(N);
    for (Index a = 0; a < N; a++) {
        Index b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a) {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}
 
// inverse fft (in-place)
void ifft(CArray& x) {
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft(x);
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}

VecArr update_density(std::vector<Particle>& particles) {
    VecArr rho = arraysVector(gp, gp, 0.0);

    for (Index p = 0; p < N; ++p) {
        Index i = std::floor(particles.at(p).position_[0] / dr);
        Index j = std::floor(particles.at(p).position_[1] / dr);
        Real hx = particles.at(p).position_[0] - (i * dr);
        Real hy = particles.at(p).position_[1] - (j * dr);

        rho[i][j] += particles.at(p).rho_c_ * (dr - hx) * (dr - hy);
        rho[i][j+1] += particles.at(p).rho_c_ * (dr - hx) * hy;
        rho[i+1][j] += particles.at(p).rho_c_ * hx * (dr - hy);
        rho[i+1][j+1] += particles.at(p).rho_c_ * hx * hy;
    }

    // for (Index u = 0; u < gp; ++u) {
    //     rho[gp - 1][u] = (rho[gp - 1][u] + rho[0][u]) * 0.5;
    //     rho[0][u] = rho[gp - 1][u];
    // }

    // for (Index u = 0; u < gp; ++u) {
    //     rho[u][gp - 1] = (rho[u][gp - 1] + rho[u][0]) * 0.5;
    //     rho[u][0] = rho[u][gp - 1];
    // }

    for (Index i = 0; i < gp; ++i) {
        rho[i] /= (dr * dr);
    }

    return rho;
}

VecArr update_potential(const VecArr& rho) {
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

    for (Index j = 0; j < gp; ++j) {
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
    Complex W = std::exp(2.0 * pi * i / Real(gp));
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

    VecArr phi = arraysVector(gp, gp, 0.0);
    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            phi[i][j] = std::real(phi_k[i][j]);
        }
    }

    return phi;
}


VecVecArr field_n(const VecArr& phi) {
    VecVecArr E = arraysVectorVector(gp, gp, 0.0);

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

VecArr field_p(const VecVecArr& field, const std::vector<Particle>& particles) {
    VecArr E = arraysVector(N, 3, 0.0);

    for (Index p = 0; p < N; ++p) {
        if (particles.at(p).move_) {
            Index i = std::floor(particles.at(p).position_[0] / dr);
            Index j = std::floor(particles.at(p).position_[1] / dr);
            Real hx = particles.at(p).position_[0] - (i * dr);
            Real hy = particles.at(p).position_[1] - (j * dr);

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

void Boris(const VecArr& E, const Array& B, std::vector<Particle>& particles) {
    Array a;
    Real a_2;
    Array b;
    Array v_minus;
    Array v_prime;
    Array v_plus;

    for (Index p = 0; p < N; ++p) {
        if (particles.at(p).move_) {
            a = 0.5 * (particles.at(p).qm_) * B * dt;
            a_2 = norm(a) * norm(a);
            b = (2.0 * a) / (1.0 + a_2);
            v_minus = particles.at(p).velocity_ + 0.5 * (particles.at(p).qm_) * E[p] * dt;
            v_prime = v_minus + cross(v_minus, a);
            v_plus = v_minus + cross(v_prime, b);
            particles.at(p).velocity_ = v_plus + 0.5 * (particles.at(p).qm_) * E[p] * dt;

            particles.at(p).position_ += particles.at(p).velocity_ * dt;
            particles.at(p).position_[0] = mod(particles.at(p).position_[0], L);
            particles.at(p).position_[1] = mod(particles.at(p).position_[1], L);
        }
    }
}

void outphase(const Real& direction, const VecArr& E, const Array& B, std::vector<Particle>& particles) {
    const Real dT = direction * 0.5 * dt;
    Array a;
    Real a_2;
    Array b;
    Array v_minus;
    Array v_prime;
    Array v_plus;

    for (Index p = 0; p < N; ++p) {
        if (particles.at(p).move_) {
            a = 0.5 * (particles.at(p).qm_) * B * dT;
            a_2 = norm(a) * norm(a);
            b = (2.0 * a) / (1.0 + a_2);
            v_minus = particles.at(p).velocity_ + 0.5 * (particles.at(p).qm_) * E[p] * dT;
            v_prime = v_minus + cross(v_minus, a);
            v_plus = v_minus + cross(v_prime, b);
            particles.at(p).velocity_ = v_plus + 0.5 * (particles.at(p).qm_) * E[p] * dT;
        }
    }
}
