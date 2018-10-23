#include <cstdlib>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <cassert>

double Distance(const std::array<double, 3>& i, const std::array<double, 3>& j)
{   
    return std::sqrt(std::pow(i[0] - j[0], 2) +
                     std::pow(i[1] - j[1], 2) +
                     std::pow(i[2] - j[2], 2));
}

double MeanDistanceConfig(const std::vector<std::array<double, 3>>& config)
{
    const int n= config.size();
    double distance_sum{0.0};
    for (auto i= 0; i < n - 1; i++)
        for (auto j= i + 1; j < n, i++)
            distance_sum += Distance(config.at(i), config.at(j));
    return distance_sum * 2 / n / (n - 1);
}


GenerateRandomCoords()

bool CheckConfigValid(const std::vector<std::array<double, 3>>& config)
{
    const int n= config.size();
    bool valid{true};
    int i{0};
    do
    {
        int j{i + 1};
        do
        {
            valid= (Distance(config[i], config[j] > 2 * R);
            j++;
        } while (j < n && valid); /// try to include ++ here later
        i++;
    } while (i < n - 1 && valid);
    return valid;
}


const std::vector<std::array<double, 3>>& GenerateConfig(const int n, const int seed)
{
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    int attempt{0};
    std::vector<std::array<double, 3>> config;
    do
    {
        std::default_random_engine rd(seed + attempt);
        for (auto i= 0; i < n; i++)
            config.push_back(std::array<double>{uniform_dist(rd),
                                                uniform_dist(rd),
                                                uniform_dist(rd)});
    } while (attempt < max_attempts && ) 

}

const std::tuple<int, double> RunExperiment1(const int n, const int M, const int seed)
{   
    double distance_sum{0.0};
    for (auto i= 0; i < M; i++)
    {
        std::vector<std::array<double, 3>> config= GenerateConfig(n, seed);
        distance_sum += MeanDistanceConfig(config);
    }
    const double ensemble_avg_distance{distance_sum / M};
    return std::make_tuple(M, ensemble_avg_distance); 
}


RunExperiment2(const int n, const int M, const int seed)



int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 4 || argv[1] == "-h") // check command line arguments and give some help
    {  
        std::cerr << "Usage: " << argv[0]
                  << " n(int):Max number of particles   M(int):Max number of configurations  s(int):Initial RNG seed"
                  << std::endl << std::endl;
        return 1;
    }
    const int n{atoi(argv[1])}; // Max number of particles
    const int M{atoi(argv[2])}; // Max number of configurations
    static const int initial_seed{atoi(argv[3])}; // Initial RNG seed

    static const int L{1}; // Box edge length
    static const double R{0.01}; // Sphere radius
    static const int max_attempts{1000}; // Maximum attempts at generating configuration

    std::vector<int> num_spheres{};
    for (auto i= n; i > 0; i /= 2)
        num_spheres.push_back(i);

    std::vector<int> num_configs{};
    for (auto i= M; i > 0; i /= 2)
        num_configs.push_back(i);
    

    std::vector<std::tuple<int, int>> results_1{};
    for (const auto& config : num_configs)
        results_1.push_back(RunExperiment1(n, config, seed));

    // std::vector<double> tresholds{0.58, 0.592746, 0.61}; // Init vector of probs and make sure p_c is included

    // for (auto& treshold : tresholds) // Run experiments on all chosen occupation probs
    // {
    //     std::cout << "Experiment_1:Sandbox_Algorithm" << std::endl;

    //     std::vector<std::tuple<int, int>> results= RunExperiment1(size, treshold, initial_seed);

    //     std::cout << "p=" << treshold << std::endl;
       
    //     for (auto& result : results)
    //         std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
    //     std::cout << "\n\n";
    // }

    // for (auto& treshold : tresholds)
    // {
    //     std::cout << "Experiment_2:Box_Counting_Algorithm" << std::endl;

    //     std::vector<std::tuple<int, int>> results= RunExperiment2(size, treshold, initial_seed);

    //     std::cout << "p=" << treshold << std::endl;
       
    //     for (auto& result : results)
    //         std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
    //     std::cout << "\n\n";
    // }

    return 0;
}
