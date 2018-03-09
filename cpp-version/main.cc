#include "functions.h"

void printVector(std::vector<double> vector)
{
    std::cout << "[";
    for (auto i = vector.begin(); i != vector.end(); ++i)
    {   
        std::cout << *i << ", ";
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char const *argv[])
{
    srand(parameters::seed);
    const int N = parameters::N;
    const double L = parameters::L;
    const double vt = L * parameters::p_vt;
    const double vd = L * parameters::p_vd;
    const double Bx = parameters::bx;
    const double By = parameters::by;
    const int gridPoints = parameters::gp;
    const double dr = L / double(gridPoints - 1);
    const int steps = parameters::steps;
    const double dt = parameters::dt;
    std::valarray<double> B = {Bx, By};

    if (int(std::sqrt(N)) * int(std::sqrt(N)) != N) 
    {
        std::cout << "Something went wrong: square root of N must be an integer." << std::endl;
        return 1;
    }

    double n = N / (L * L);
    double margin = dr / 10.0;

    std::mt19937_64 engine;
    std::normal_distribution<double> maxwell_right(-vd, vt);
    std::normal_distribution<double> maxwell_left(vd, vt);
    std::vector<std::valarray<double>> pos;

    double deltaL = (L - 2.0 * margin) / (std::sqrt(N) - 1);

    for (double x = margin; x <= (L - margin); x += deltaL)
    {
        for (double y = margin; y <= (L - margin); y += deltaL)
        {
            // std::cout << x << " " << y << std::endl;
            pos.push_back({x, y});
        }
    }

    std::vector<double> indexes;
    std::vector<double> right;
    std::vector<double> left;
    std::vector<Particle> parts;

    for (int i = 0; i < N; ++i)
    {
        indexes.push_back(i);
    }

    // std::random_shuffle(indexes.begin(), indexes.end());

    for (int i = N-1; i >= int(3*N/4); --i)
    {
        left.push_back(indexes[i]);
        indexes.pop_back();
    }

    for (int i = int(3*N/4)-1; i >= int(N/2); --i)
    {
        right.push_back(indexes[i]);
        indexes.pop_back();
    }

    for (auto i = indexes.begin(); i != indexes.end(); ++i)
    {
        parts.push_back(Particle(pos[*i], {0.0, 0.0}, n, 1.0, false));
    }
    for (auto i = left.begin(); i != left.end(); ++i)
    {
        // parts.push_back(Particle(pos[*i], {maxwell_left(engine), 0.0}, n, -1.0, true));
        parts.push_back(Particle(pos[*i], {0.02, 0.0}, n, -1.0, true));
    }    
    for (auto i = right.begin(); i != right.end(); ++i)
    {
        // parts.push_back(Particle(pos[*i], {maxwell_right(engine), 0.0}, n, -1.0, true));
        parts.push_back(Particle(pos[*i], {-0.02, 0.0}, n, -1.0, true));
    }


    std::vector<double> rho_c;

    for (int i = 0; i < N; ++i)
    {
        rho_c.push_back(parts[i].q_ / (dr * dr));
    }

    // for (int p = 0; p < N; ++p)
    // {
    //     std::cout << parts[p].velocity_[0] << " " << parts[p].velocity_[1] << std::endl;
    // }
    std::vector<Particle> finalParts;

    for (int i = 0; i < steps; ++i)
    {
        VecVal RHO = density(parts, rho_c, dr);
        VecVal PHI = potential(RHO, dr);
        VecVecVal EFIELDn = EField_GP(PHI, dr);
        VecVal EFIELDp = EField_P(EFIELDn, parts, dr);

        if (i == 0)
        {
            rewind(-1.0, EFIELDp, B, parts);
        }

        Boris(EFIELDp, B, parts);
        finalParts = parts;
        rewind(1.0, EFIELDp, B, finalParts);

    }

    return 0;
}