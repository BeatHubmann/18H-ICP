#include <cstdlib>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <cassert>

static const double L{1.0}; // Box edge length
static const int max_attempts{1000000}; // Number of MC attempts per config

// Returns Euclidean 3d distance between coordinates i and j
double Distance(const std::array<double, 3>& i, const std::array<double, 3>& j)
{   
    return std::sqrt(std::pow(i[0] - j[0], 2) +
                     std::pow(i[1] - j[1], 2) +
                     std::pow(i[2] - j[2], 2));
}

// Calculates the mean distance between all points of a given configuration
double MeanDistanceConfig(const std::vector<std::array<double, 3>>& config)
{
    const auto n{config.size()};
    double distance_sum{0.0};
    for (auto i= 0; i < n - 1; i++)
        for (auto j= i + 1; j < n; j++)
           distance_sum += Distance(config[i], config[j]);
    return distance_sum * 2 / n / (n - 1);
}

// Checks if a given configuration is valid i.e. has no overlapping spheres
bool CheckConfigValid(const std::vector<std::array<double, 3>>& config, const double radius)
{
    const auto n{config.size()};
    bool valid{true};
    int i{0};
    do
    {
        int j{i + 1};
        do
        {
            valid= (Distance(config[i], config[j]) > 2 * radius);
        } while (j++ < n && valid);
    } while (i++ < n - 1 && valid);
    return valid;
}

// Attempts to generate a valid configuration by sampling all coordinates for all spheres at once and only keeping the config if valid
bool GenerateConfig(std::vector<std::array<double, 3>>& config, const int n, const double radius, const int seed)
{
    std::uniform_real_distribution<double> uniform_dist(0.0, L);
    bool valid{false};
    int attempt{0};
    do
    {
        config.clear();
        std::mt19937 rd(seed + attempt);
        for (auto i= 0; i < n; i++)
            config.push_back(std::array<double, 3>{{uniform_dist(rd),
                                                    uniform_dist(rd),
                                                    uniform_dist(rd)}});
        valid= CheckConfigValid(config, radius);
        attempt++;
    } while (attempt < max_attempts && !valid);
    return valid;
}

// Conducts one iteration of the experiment for given number of spheres n and number of configurations M
const std::tuple<int, int, double> RunExperiment(const int n, const int M, const double radius, const int initial_seed)
{   
    double distance_sum{0.0};
    for (auto i= 0; i < M; i++)
    {
        std::vector<std::array<double, 3>> config{};
        bool success{GenerateConfig(config, n, radius, initial_seed + i * max_attempts)};
        assert(success && "No valid configuration found: Adjust sphere number and/or radius.");
        distance_sum += MeanDistanceConfig(config);
    }
    return std::make_tuple(n, M, distance_sum / M); // <ensemble size, ensemble avg distance>
}

int main(int argc, char* argv[])
{
    if (argc < 5 || argc > 5 || argv[1] == "-h") // check command line arguments and give some help
    {  
        std::cerr << "Usage: " << argv[0]
                  << " n(int): Number of particles     M(int): Max number of configurations   s(int): Initial RNG seed"
                  << " R(double): Sphere radius \n"
                  << std::endl << std::endl;
        return 1;
    }
    const int n{atoi(argv[1])}; // Max number of particles
    const int M{atoi(argv[2])}; // Max number of configurations
    const int initial_seed{atoi(argv[3])}; // Initial RNG seed
    const double R{atof(argv[4])}; // Sphere radius

    const double density{n * 4. / 3. * M_PI * R*R*R / (L*L*L)}; // For reference

    std::vector<int> num_configs{}; // Generate vector of Ms to conduct experiments for
    for (auto i= 1; i <= M; i *= 2)
        num_configs.push_back(i);
    
    std::vector<int> num_spheres{}; // Generate vector of ns to conduct experiments for
    for (auto i= 2; i < n + 1; i *= 2)
        num_spheres.push_back(i);

    for (const auto& spheres: num_spheres) // Run experiment for all ns
    {
        std::vector<std::tuple<int, int, double>> results{};
        for (const auto& config : num_configs) // For each n, run up to max ensemble size M
            results.push_back(RunExperiment(spheres, config, R, initial_seed));
        for (auto& result : results) // Print results
            std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << "\t" << std::get<2>(result)
                      << std::endl;
        std::cout << std::endl << std::endl; 
    }
    return 0;
}
