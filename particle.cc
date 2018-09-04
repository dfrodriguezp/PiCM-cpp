#include "particle.h"

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
