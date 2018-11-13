#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>

double SquareDistance(const std::vector<double> a, const std::vector<double> b)
{
    assert (a.size() == b.size());
    double square_dist{0.0};
    for (auto i= 0; i < a.size(); i++)
    {
        square_dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return square_dist;
}

double AbsDistance(const std::vector<double> a, const std::vector<double> b)
{
    return std::sqrt(SquareDistance(a, b));
}

bool Overlap(const std::vector<double> a, const std::vector<double> b,
             const double r)
{
    return (AbsDistance(a, b) < r);
}

double RandomWalk(const int N, const int seed, const bool self_avoid=false,
                  const bool spherical=false,
                  const double segment_length=1.0,
                  const double r= 0.0)
{
    std::mt19937 rd(seed);
    std::uniform_real_distribution<> angle(0, M_PI);

    std::vector<double> pos(3, 0.0), pos_new(3, 0.0);
    double x{0.0}, y{0.0}, z{0.0};
    double phi{0.0};
    double theta{M_PI_2};
    double sin_theta{1.0}, cos_theta{0.0};
    int num_segments{0};

    while (num_segments < N)
    {
        phi= 2 * angle(rd);
        if (spherical)
        {
            theta= angle(rd);
            sin_theta= std::sin(theta);
            cos_theta= std::cos(theta);
        }
        pos_new[0]= pos[0] + segment_length * sin_theta * std::cos(phi);
        pos_new[1]= pos[1] + segment_length * sin_theta * std::sin(phi);
        pos_new[2]= pos[2] + segment_length * cos_theta;
        if (!self_avoid || !Overlap(pos, pos_new, r))
        {
            pos= pos_new;
            num_segments++;
            // std::cout << "Segment no. " << num_segments << "\t"
        //               << "x:\t" << pos[0] << "\t"
        //               << "y:\t" << pos[1] << "\t"
        //               << "z:\t" << pos[2] << std::endl;
        }
    }
    return SquareDistance(std::vector<double>(3, 0.0), pos);
}

void RunExp1a(const int N, const int M, const int seed, const double segment_length)
{
    double sum_R_squared{0.0};
    double sum_R_squared_squared{0.0};
    for (auto k= 1; k <= M; k++)
    {
        const double ThisWalk{RandomWalk(N, seed + k)};
        sum_R_squared += ThisWalk;
        sum_R_squared_squared += ThisWalk * ThisWalk;
        if (k % 10 == 0) // output every tenth step
        {
            const double ens_avg_R_squared{sum_R_squared / k};
            std::cout << k << "\t"
                      << ens_avg_R_squared << "\t"
                      << std::sqrt(((sum_R_squared_squared / k) - ens_avg_R_squared * ens_avg_R_squared) / k) << std::endl;
        }
    }
}

void RunExp1b(const int N, const int M, const int seed, const double segment_length)
{
    
}

void RunExp2a(const int N, const int M, const int seed, const double segment_length,
               const bool self_avoid, const bool spherical, const double r)
{

}

void RunExp2b(const int N, const int M, const int seed, const double segment_length,
               const bool self_avoid, const bool spherical, const double r)
{

}


int main(int argc, char* argv[])
{
    if (argc < 8 || argc > 8 || argv[1] == "-h") // check command line arguments and give some help
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

    RunExp1a(N, M, seed, segment_length);

    RunExp1b(N, M, seed, segment_length);

    RunExp2a(N, M, seed, segment_length, true, true, r);

    RunExp2b(N, M, seed, segment_length, true, true, r);


    // for (auto k=0; k < M; k++)
    // {
    //     std::cout << RandomWalk(N, seed + k, self_avoid, spherical, segment_length, r) << std::endl;
    // }

    return 0;
}
