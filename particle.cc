#include "particle.h"

Real norm(const Array& A) {
    return std::sqrt((A * A).sum());
}

inline Array cross(const Array& A, const Array& B) {
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

Real mod(const Real& a, const Real& b) {
    return (a < 0) ? std::fmod((a + (std::floor(-a / b) + 1) * b), b) : (std::fmod(a, b));
}

Particle::Particle(Array pos, Array vel, Real n_r, Real qm, Index move) {
    this -> position_ = pos;
    this -> velocity_ = vel;
    this -> qm_ = qm;
    this -> charge_ = (1 / qm) * (1 / n_r);
    this -> mass_ = this -> charge_ / qm;
    this -> move_ = move;
    this -> rho_c_ = this -> charge_ / (dr * dr);
    this -> Efield_ = {0.0, 0.0, 0.0};
}

Particle::~Particle() {

}

const Array& Particle::getPosition() const {
    return this -> position_;
}

const Array& Particle::getVelocity() const {
    return this -> velocity_;
}

const Real& Particle::getRho_c() const {
    return this -> rho_c_;
}

const Real& Particle::getMass() const {
    return this -> mass_;
}

const Index& Particle::Move() const {
    return this -> move_;
}

void Particle::fieldInfluence(const VecVecArr& meshField) {
    this -> Efield_ = {0.0, 0.0, 0.0};
    if (this -> move_ == 1) {
        Index i = std::floor(this -> position_[0] / dr);
        Index j = std::floor(this -> position_[1] / dr);
        Real hx = this -> position_[0] - (i * dr);
        Real hy = this -> position_[1] - (j * dr);

        Real A = (dr - hx) * (dr - hy);
        Real B = (dr - hx) * hy;
        Real C = hx * (dr - hy);
        Real D = hx * hy;

        this -> Efield_[0] = meshField[i][j][0] * A + meshField[i][j+1][0] * B + meshField[i+1][j][0] * C + meshField[i+1][j+1][0] * D;
        this -> Efield_[1] = meshField[i][j][1] * A + meshField[i][j+1][1] * B + meshField[i+1][j][1] * C + meshField[i+1][j+1][1] * D;
    }

    this -> Efield_ /= (dr * dr);
}

void Particle::update(const Array& B) {
    Array a;
    Array b;
    Array v_minus;
    Array v_prime;
    Array v_plus;
    Real a_2;

    if (this -> move_ == 1) {
        a = 0.5 * this -> qm_ * B * dt;
        a_2 = norm(a) * norm(a);
        b = (2.0 * a) / (1.0 + a_2);
        v_minus = this -> velocity_ + 0.5 * this -> qm_ * this -> Efield_ * dt;
        v_prime = v_minus + cross(v_minus, a);
        v_plus = v_minus + cross(v_prime, b);
        this -> velocity_ = v_plus + 0.5 * this -> qm_ * this -> Efield_ * dt;

        this -> position_ += this -> velocity_ * dt;
        this -> position_[0] = mod(this -> position_[0], L);
        this -> position_[1] = mod(this -> position_[1], L);
    }   
}

void Particle::outPhase(const Index& direction, const Array& B) {
    const Real dT = direction * 0.5 * dt;
    Array a;
    Array b;
    Array v_minus;
    Array v_prime;
    Array v_plus;
    Real a_2;

    if (this -> move_ == 1) {
        a = 0.5 * this -> qm_ * B * dT;
        a_2 = norm(a) * norm(a);
        b = (2.0 * a) / (1.0 + a_2);
        v_minus = this -> velocity_ + 0.5 * this -> qm_ * this -> Efield_ * dT;
        v_prime = v_minus + cross(v_minus, a);
        v_plus = v_minus + cross(v_prime, b);
        this -> velocity_ = v_plus + 0.5 * this -> qm_ * this -> Efield_ * dT;  
    }

}
