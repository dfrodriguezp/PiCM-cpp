#ifndef MESH_H
#define MESH_H

#include "particle.h"

VecVecArr arraysVectorVector(const Index& rows, const Index& cols, const Real& value);

class Mesh
{
public:
    Mesh();
    ~Mesh();

    const VecArr&    getRho() const;
    const VecArr&    getPhi() const;
    const VecVecArr& getEfield() const;

    void updateDensity(std::vector<Particle>& particles);
    void updatePotential();
    void updateEfield();

private:
    VecArr    rho_;
    VecArr    phi_;
    VecVecArr E_n_;
};

#endif