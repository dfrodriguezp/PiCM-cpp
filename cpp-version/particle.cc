#include <valarray>

struct Particle
{

    std::valarray<double> position_;
    std::valarray<double> velocity_;
    double mass_;
    double electronicDensity_;
    bool move_;

    Particle(std::valarray<double> position,
        std::valarray<double> velocity,
        double mass, double electronicDensity, bool move
        ) : position_(position), velocity_(velocity),
            mass_(mass), electronicDensity_(electronicDensity), move_(move)

    {

    }
};