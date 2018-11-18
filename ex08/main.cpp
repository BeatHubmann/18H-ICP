#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <cassert>

double SquareDistance(const std::vector<double> a,
                      const std::vector<double> b=std::vector<double>(3, 0.0))
{
    assert (a.size() == b.size());
    double square_dist{0.0};
    for (auto i= 0; i < a.size(); i++)
    {
        const double comp_dist{a[i] - b[i]};
        square_dist += comp_dist * comp_dist;
    }
    return square_dist;
}

double AbsDistance(const std::vector<double> a,
                   const std::vector<double> b=std::vector<double>(3, 0.0))
{
    return std::sqrt(SquareDistance(a, b));
}

bool Overlap(const std::vector<double> a, const std::vector<double> b,
             const double r)
{
    return (AbsDistance(a, b) < r);
}

bool OverlapChain(const std::vector<std::vector<double>>& chain,
                  const std::vector<double>& pos_new,
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
                  const bool spherical=false,
                  const double segment_length=1.0,
                  const double r= 0.0)
{
    std::mt19937 rd(seed);
    std::uniform_real_distribution<> angle(0, M_PI);
    std::vector<std::vector<double>> chain;
    chain.reserve(N);
    chain.push_back(std::vector<double>(3, 0.0));
    std::vector<double> pos_new(3, 0.0);
    double phi{0.0};
    double theta{M_PI_2};
    double sin_theta{1.0}, cos_theta{0.0};

    while (chain.size() < N)
    {
        phi= 2 * angle(rd);
        if (spherical)
        {
            theta= angle(rd);
            sin_theta= std::sin(theta);
            cos_theta= std::cos(theta);
        }
        pos_new[0]= chain.back()[0] + segment_length * sin_theta * std::cos(phi);
        pos_new[1]= chain.back()[1] + segment_length * sin_theta * std::sin(phi);
        pos_new[2]= chain.back()[2] + segment_length * cos_theta;
        if (!self_avoid || !OverlapChain(chain, pos_new, r))
        {
            // std::cout.precision(4);
            // std::cout << "Segment no. " << chain.size() << "\t\t"
            //           << "x:\t"  << pos_new[0] << "\t\t"
            //           << "y:\t"  << pos_new[1] << "\t\t"
            //           << "z:\t"  << pos_new[2] << "\t\t"
            //           << "Dst:\t"<< AbsDistance(chain.back(), pos_new) << "\t\t"
            //           << "Ttl:\t"<< AbsDistance(pos_new) << std::endl;
            chain.push_back(pos_new);
        }
    }
    return SquareDistance(chain.back());
}

void RunExp1(const int N, const int M, const int seed,
             const double segment_length = 1.0,
             const bool self_avoid = false,
             const bool spherical = false,
             const double r = 1.0)
{
    double sum_R_squared{0.0};
    double sum_R_squared_squared{0.0};
    for (auto k= 1; k <= M; k++)
    {
        const double ThisWalk{RandomWalk(N, seed + k,
                                         self_avoid,
                                         spherical,
                                         segment_length,
                                         r)};
        sum_R_squared += ThisWalk;
        sum_R_squared_squared += ThisWalk * ThisWalk;
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
        for (auto k= 1; k <= M; k++)
        {
            const double ThisWalk{RandomWalk(i, seed + k,
                                            self_avoid,
                                            spherical,
                                            segment_length,
                                            r)};
            sum_R_squared += ThisWalk;
            sum_R_squared_squared += ThisWalk * ThisWalk;
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

    // const std::vector<double> a(3, 2.0), b(3, 3.0);
    // std::cout << "Testing:\n"

    //           << "AbsDist:" << AbsDistance(a, b) << "\n"
    //           << "Overlap1:"<< Overlap(a, b, r) << std::endl;
    RunExp1(N, M, seed, segment_length, self_avoid, spherical, r);

    RunExp2(N, M, seed, segment_length, self_avoid, spherical, r);


    // for (auto k=0; k < M; k++)
    // {
    //     std::cout << RandomWalk(N, seed + k, self_avoid, spherical, segment_length, r) << std::endl;
    // }

    return 0;
}
