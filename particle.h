#ifndef PARTICLE_H
#define PARTICLE_H

#include "initial.h"

Real norm(const Array& A);
Real mod(const Real& a, const Real& b);

class Particle
{
public:
    Particle(Array pos, Array vel, Real n_r, Real qm, Index move);
    ~Particle();

    const Array& getPosition() const;
    const Array& getVelocity() const;
    const Real& getRho_c() const;
    const Real& getMass() const;
    const Index& Move() const;
    void fieldInfluence(const VecVecArr& meshField);
    void update(const Array& B);
    void outPhase(const Index& direction, const Array& B);

private:
    Array position_;
    Array velocity_;
    Real  qm_;
    Real  charge_;
    Real  mass_;
    Index move_;
    Real  rho_c_;
    Array  Efield_;
};

#endif