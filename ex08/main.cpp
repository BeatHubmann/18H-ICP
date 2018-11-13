#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#include "latticeview.h"

bool Coalescence(const std::vector<int>& lattice, const int L, const int x, const int y)
{
    return ((x - 1 > -1 && lattice[x - 1 + L *  y     ] > 0) ||
            (x + 1 <  L && lattice[x + 1 + L *  y     ] > 0) ||
            (y - 1 > -1 && lattice[x     + L * (y - 1)] > 0) ||
            (y + 1 <  L && lattice[x     + L * (y + 1)] > 0));
}

void RandomWalk(std::vector<int>& lattice, const int L, const int seed, const double init_speed,
           const int expiry, const int treshold, const int num_particles)
{
    for (auto p= 0; p < num_particles; p++)
    {
        std::mt19937 rd(seed + p);
        std::uniform_int_distribution<int> coord(0, 4 * L - 1); // particles may launch all around the box
        std::uniform_real_distribution<> angle(0, 2 * M_PI);
        std::uniform_real_distribution<> speed(1.0 * init_speed, 1.5 * init_speed);
        const int start{coord(rd)};
        double x_double{0};
        double y_double{0};
        switch(start / L) // determine on which border particles launches:
        {
            case 0:                      y_double= start % L; break; // west
            case 1: x_double= start % L; y_double= L - 1    ; break; // north
            case 2: x_double= L - 1    ; y_double= start % L; break; // east
            case 3: x_double= start % L;                      break; // south
        }
        double spd{init_speed}; // particle velocity
        double ang{angle(rd)}; // particle angle
        int num_steps{0};
        bool diffusing{true};
        while (diffusing && num_steps < expiry)
        {
            ang += angle(rd);
            spd= speed(rd);
            x_double += std::cos(ang) * spd;
            if (x_double > (double)L - 1)
                x_double -= L; // wrap around edges
            if (x_double < 0)
                x_double += L;
            y_double += std::sin(ang) * spd;
            if (y_double > (double)L - 1)
                y_double -= L; // wrap around edges
            if (y_double < 0)
                y_double += L;
            const int x{(int)x_double};
            const int y{(int)y_double};

            if (Coalescence(lattice, L, x, y))
            {
                lattice[x + L * y] -= 1;
                if (-lattice[x + L * y] > treshold)
                {
                    lattice[x + L * y]= (p / (num_particles / 7)) % 7 + 1; // got 7 colors
                    diffusing= false;
                }
            } 
            num_steps++;    
        }
    }
}



int main(int argc, char* argv[])
{
    if (argc < 7 || argc > 7 || argv[1] == "-h") // check command line arguments and give some help
    {  
        std::cerr << "Usage: " << argv[0]
                  << " L(int): lattice side length"
                  << " seed(int): Initial RNG seed"
                  << " speed(double): Initial particle speed"
                  << " expiry(int): Particle lifetime if no coalescence"
                  << " treshold(int): How many times site needs visiting before coalescence"
                  << " num_particles(int): How many particles to be thrown in"
                  << std::endl << std::endl;
        return 1;
    }
    static const int L{atoi(argv[1])}; // lattice side length  
    const int seed{atoi(argv[2])}; // initial RNG seed
    const double init_speed{atof(argv[3])}; // initial particle speed
    const int expiry{atoi(argv[4])}; // how many steps until particle expires when not making contact
    const int treshold{atoi(argv[5])}; // how many times a site has to be visited before coalescence
    const int num_particles{atoi(argv[6])}; // how many particles are thrown into the system

    std::vector<int> lattice(L * L, 0);
    for (auto i= (L / 2) - 1; i < (L / 2) + 1; i++)
        for (auto j= (L / 2) - 1; j < (L / 2) + 1; j++)
            lattice[L * i + j]= 1; // seed is center cluster of 4 sites

    RandomWalk(lattice, L, seed, init_speed, expiry, treshold, num_particles);

    const int image_size{std::max(L, 1024)};
    Print_lattice(lattice, L, L, image_size, image_size);

    return 0;
}