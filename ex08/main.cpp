#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <cassert>

double SquareDistance(const std::array<double, 3> a,
                      const std::array<double, 3> b = std::array<double, 3>{0.0, 0.0, 0.0})
{
    double square_dist{0.0};
    for (auto i= 0; i < a.size(); i++)
    {
        const double comp_dist{a[i] - b[i]};
        square_dist += comp_dist * comp_dist;
    }
    return square_dist;
}

double AbsDistance(const std::array<double, 3> a,
                   const std::array<double, 3> b=std::array<double, 3>{0.0, 0.0, 0.0})
{
    return std::sqrt(SquareDistance(a, b));
}

bool Overlap(const std::array<double, 3> a, const std::array<double, 3> b,
             const double r)
{
    return (AbsDistance(a, b) < r);
}

bool OverlapChain(const std::vector<std::array<double, 3>>& chain,
                  const std::array<double, 3>& pos_new,
                  const double r)
{
    bool overlap{false};
    for (auto& element : chain)
    {
        if (Overlap(element, pos_new, r))
        {
            overlap= true;
            break;
        }
    }
    return overlap;
}

double RandomWalk(const int N, const int seed, const bool self_avoid=false,
                  const bool spherical=false,          // walk in 3D i.s.o. 2D
                  const double segment_length=1.0,     // polar coordinates
                  const double r= 0.0,                 // particle radius
                  const int max_attempts= 100)         // when to give up if stuck
{
    std::mt19937 rd(seed);
    std::uniform_real_distribution<> angle(0, M_PI);
    std::vector<std::array<double, 3>> chain; 
    chain.reserve(N); // reserve memory space for efficiency
    std::array<double, 3> pos_new{0.0, 0.0, 0.0}; // place holder for walk steps
    chain.push_back(pos_new); // use origin as starting position
    double phi{0.0}; // polar coordinates
    double theta{M_PI_2}; // polar coordinates
    double sin_theta{1.0}, cos_theta{0.0}; // precalculate for 2D and/or efficiency
    int attempts{0}; // for tracking attempts when getting stuck

    while (chain.size() < N)
    {
        phi += 2 * angle(rd);
        if (spherical) // in case we're in 3D (as in task 2)
        {
            theta += angle(rd);
            sin_theta= std::sin(theta);
            cos_theta= std::cos(theta);
        }
        pos_new[0]= chain.back()[0] + segment_length * sin_theta * std::cos(phi);
        pos_new[1]= chain.back()[1] + segment_length * sin_theta * std::sin(phi);
        pos_new[2]= chain.back()[2] + segment_length * cos_theta;
        if (!self_avoid || !OverlapChain(chain, pos_new, r)) // either no collision or don't care
        {
            chain.push_back(pos_new); // add particle to chain
            attempts= 0; // reset attempts in case we were probing for new position
        }
        else // we are self-avoiding and have a collision
        {
            if (attempts++ == max_attempts) break; // seems we're trapped - abort this walk
        }
    }
    if (chain.size() == N) // successful walk
        return SquareDistance(chain.back()); // return length
    else // aborted walk because we got stuck 
        return -1.0; // error value
}

void RunExp1(const int N, const int M, const int seed,
             const double segment_length = 1.0,
             const bool self_avoid = false,
             const bool spherical = false,
             const double r = 1.0)
{
    double sum_R_squared{0.0};
    double sum_R_squared_squared{0.0};
    int k{1};
    int seed_shift{0};
    while (k < M)
    {
        const double ThisWalk{RandomWalk(N, seed + seed_shift,
                                         self_avoid,
                                         spherical,
                                         segment_length,
                                         r)};
        seed_shift++;
        if (ThisWalk > 0)
        {
            sum_R_squared += ThisWalk;
            sum_R_squared_squared += ThisWalk * ThisWalk;
            k++;
            if (k % 10 == 0) // output every tenth step
            {
                const double ens_avg_R_squared{sum_R_squared / k};
                std::cout << k << "\t"
                        << ens_avg_R_squared << "\t"
                        << std::sqrt(((sum_R_squared_squared / k)
                                        - ens_avg_R_squared * ens_avg_R_squared) / k)
                        << std::endl;
            }
        }
    }
}

void RunExp2(const int N, const int M, const int seed,
             const double segment_length = 1.0,
             const bool self_avoid = false,
             const bool spherical = false,
             const double r = 1.0)
{
    for (auto i= 8; i <= N; i *=2)
    {
        double sum_R_squared{0.0};
        double sum_R_squared_squared{0.0};
        int k{1};
        int seed_shift{0};
        while (k < M)
        {
            const double ThisWalk{RandomWalk(i, seed + seed_shift,
                                            self_avoid,
                                            spherical,
                                            segment_length,
                                            r)};
            seed_shift++;
            if (ThisWalk > 0)
            {
                sum_R_squared += ThisWalk;
                sum_R_squared_squared += ThisWalk * ThisWalk;
                k++;
            }
        }
        const double ens_avg_R_squared{sum_R_squared / M};
        std::cout << i << "\t"
                  << ens_avg_R_squared << "\t"
                  << std::sqrt(((sum_R_squared_squared / M)
                                 - ens_avg_R_squared * ens_avg_R_squared) / M)
                  << std::endl;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 8 || argc > 8 || argv[1] == "-h") // check cl args and give some help
    {  
        std::cerr << "Usage: " << argv[0] << "\n\t"
                  << " N(int): Length of random walk" << "\n\t"
                  << " M(int): Number of configurations" << "\n\t"
                  << " Seed(int): Initial RNG seed" << "\n\t"
                  << " Segment length(double): Random walk fixed segment length" << "\n\t"
                  << " Self avoid (0/1): N/Y on self-avoiding" << "\n\t"
                  << " 3D (0/1):N/Y on using 3D coordinates i.s.o. 2D coordinates" << "\n\t"
                  << " Particle radius(double): Particle radius for self-avoiding random walk"
                  << std::endl << std::endl;
        return 1;
    }
    const int N{atoi(argv[1])}; 
    const int M{atoi(argv[2])};
    const int seed{atoi(argv[3])};
    const double segment_length{atof(argv[4])};
    const bool self_avoid{(bool)atoi(argv[5])};
    const bool spherical{(bool)atoi(argv[6])};
    const double r{atof(argv[7])};
    assert(segment_length >= r);

    RunExp1(N, M, seed, segment_length, self_avoid, spherical, r); // subtask a

    std::cout << "\n\n"; // separate outputs blocks for gnuplot

    RunExp2(N, M, seed, segment_length, self_avoid, spherical, r); // subtask b

    return 0;
}
