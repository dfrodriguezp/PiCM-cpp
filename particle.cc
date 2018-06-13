#ifndef PARTICLE
#define PARTICLE

#include <valarray>
#include <string>
struct Particle {
    std::valarray<double> position_;
    std::valarray<double> velocity_;
    double qm_;
    double q_;
    double m_;
    bool move_;

    Particle(std::valarray<double> position,
        std::valarray<double> velocity,
        double electronicDensity, double qm, bool move
        ) {
        position_ = position;
        velocity_ = velocity;
        qm_ = qm;
        q_ = (1 / qm) * (1 / electronicDensity);
        m_ = q_ / qm_;
        move_ = move;
    }
};

#endif