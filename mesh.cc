#include "mesh.h"
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

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
void fft(CArray& x) {
    // DFT
    Index N = x.size(), k = N, n;
    Real thetaT = pi / N;
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

Mesh::Mesh() {
    this -> rho_ = arraysVector(gp, gp, 0.0);
    this -> phi_ = arraysVector(gp, gp, 0.0);
    this -> E_n_ = arraysVectorVector(gp, gp, 0.0);
}

Mesh::~Mesh() {

}

const VecArr& Mesh::getRho() const{
    return this -> rho_;
}

const VecArr& Mesh::getPhi() const{
    return this -> phi_;
}

const VecVecArr& Mesh::getEfield() const{
    return this -> E_n_;
}

void Mesh::updateDensity(std::vector<Particle>& particles) {
    this -> rho_ = arraysVector(gp, gp, 0.0);    
    for (Index p = 0; p < particles.size(); ++p) {
        Index i = std::floor(particles[p].getPosition()[0] / dr);
        Index j = std::floor(particles[p].getPosition()[1] / dr);
        Real hx = particles[p].getPosition()[0] - (i * dr);
        Real hy = particles[p].getPosition()[1] - (j * dr);

        this -> rho_[i][j] += particles[p].getRho_c() * (dr - hx) * (dr - hy);
        this -> rho_[i][j+1] += particles[p].getRho_c() * (dr - hx) * hy;
        this -> rho_[i+1][j] += particles[p].getRho_c() * hx * (dr - hy);
        this -> rho_[i+1][j+1] += particles[p].getRho_c() * hx * hy;
    }

    for (Index u = 0; u < gp; ++u) {
        this -> rho_[gp-1][u] = (this -> rho_[gp-1][u] + this -> rho_[0][u]) * 0.5;
        this -> rho_[0][u] = this -> rho_[gp-1][u];
    }   

    for (Index u = 0; u < gp; ++u) {
        this -> rho_[u][gp-1] = (this -> rho_[u][gp-1] + this -> rho_[0][u]) * 0.5;
        this -> rho_[u][0] = this -> rho_[u][gp-1];
    }   

    for (Index i = 0; i < gp; ++i) {
        this -> rho_[i] /= (dr * dr);
    }
}

void Mesh::updatePotential() {
    Complex rho_k[gp][gp];

    for (Index i = 0; i < gp; ++i) {
           for (Index j = 0; j < gp; ++j) {
               rho_k[i][j] = this -> rho_[i][j];
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

    this -> phi_ = arraysVector(gp, gp, 0.0);
    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            this -> phi_[i][j] = std::real(phi_k[i][j]);
        }
    }
}

void Mesh::updateEfield() {
    this -> E_n_ = arraysVectorVector(gp, gp, 0.0);

    for (Index j = 0; j < gp; ++j) {
        for (Index i = 0; i < gp; ++i) {
            Index nxt_i = (i < gp - 1) ? (i + 1) : 0;
            Index prv_i = (i > 0) ? (i - 1) : (gp - 1);

            this -> E_n_[i][j][0] = (this -> phi_[prv_i][j] - this -> phi_[nxt_i][j]) / (2 * dr);
        }
    }

    for (Index i = 0; i < gp; ++i) {
        for (Index j = 0; j < gp; ++j) {
            Index nxt_j = (j < gp - 1) ? (j + 1) : 0;
            Index prv_j = (j > 0) ? (j - 1) : (gp - 1);

            this -> E_n_[i][j][1] = (this -> phi_[i][prv_j] - this -> phi_[i][nxt_j]) / (2 * dr);
        }
    }
}