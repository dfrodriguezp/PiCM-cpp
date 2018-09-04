#ifndef PARTICLE_H
#define PARTICLE_H

#include "initial.h"


class Particle
{
public:
    Particle(Array pos, Array vel, Real n_r, Real qm, Index move);
    ~Particle();

    Array position_;
    Array velocity_;
    Real  qm_;
    Real  charge_;
    Real  mass_;
    Index move_;
    Real  rho_c_;
    Array  Efield_;
};

VecArr update_density(std::vector<Particle>& particles);
VecArr field_p(const VecVecArr& field, const std::vector<Particle>& particles);
void Boris(const VecArr& E, const Array& B, std::vector<Particle>& particles);
void outphase(const Real& direction, const VecArr& E, const Array& B, std::vector<Particle>& particles);

#endif