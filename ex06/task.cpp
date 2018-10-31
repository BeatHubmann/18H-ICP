#include <cstring>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <cassert>

static const int L{200}; // lattice side length
static const double J{1.0}; // exchange coupling constant

const int h(const std::vector<int> grid, const int i) // calculate h(i) w/ periodic b.c.
{
    const int row{i / L};
    const int col{i % L};
    const int west{col - 1 < 0 ? i + L-1  : i - 1};
    const int east{col + 1 < L ? i + 1    : i - (L-1)};
    const int north{row - 1 < 0 ? i + (L-1) * L : i - L};
    const int south{row + 1 < L ? i + L         : i - (L-1) * L};
    return grid[west] + 
           grid[east] +
           grid[north] +
           grid[south];
}

void PrintGrid(const std::vector<int> grid) // only for debugging w/ L <= 20
{
    for (auto i= 0; i < L; i++)
    {
        for (auto j= 0; j < L; j++)
            std::cout << (grid[i*L + j] == 1 ? (char)46 : (char)35) << "\t";
        std::cout << "\n";
    }
    std::cout <<"\n\n\n";
}


void FillLookup(double exp_dE[], const double T) // precalculate values for MC step
{
    for (auto s_i= -1; s_i < 2; s_i += 2) // s_i = -1, +1
        for (auto h_i= -4; h_i < 5; h_i += 2) // h_i = -4, -2, 0, +2, +4
            exp_dE[(s_i + 1) / 2 * 5 + (h_i + 4) / 2]= std::exp(-2 * (double)s_i * (double)h_i * J / T); 
}

double GetExpDeltaE(const int s_i, const int h_i, const double exp_dE[]) // lookup precalculated values in MC step
{
    return exp_dE[(s_i + 1) / 2 * 5 + (h_i + 4) / 2];
}

// core function:
const std::tuple<double, double, std::vector<double>, double, double, std::vector<double>> 
    RunExperiment(std::vector<int>& grid, double E, double M, const double T, const int steps, const int seed, const bool cold_start = false)
{
    if (cold_start)
    {
        grid.clear(); // reset grid to cold ground state
        for (auto i= 0; i < L * L; i++)
            grid.push_back(1);
        M= L * L; // reset M
        E= -J * 4 * L * L / 2; // reset E
    }

    double exp_dE[10]; // fixed size: s = -1, +1; h = -4, -2, 0, 2, 4
    FillLookup(exp_dE, T); // precalculate values for MC step

    std::mt19937 rd(seed); // random engine
    std::uniform_int_distribution<int> choice(0, L * L - 1); // for choosing site
    std::uniform_real_distribution<double> accept(0.0, 1.0); // for Metropolis MC

    std::vector<double> magnet_steps{M}; // for analysis/plotting
    std::vector<double> energy_steps{E};

    double sum_M{0}; // for equilibrium averages
    double sum_E{0};

    for (auto s= 0; s < steps; s++) // perform Metropolis MC steps
    {
        const int i{choice(rd)}; // choose site i
        const int h_i{h(grid, i)}; // get h_i
        const double delta_energy{2 * grid[i] * h_i * J}; // calculate dE
        if (delta_energy < 0) // accept
        {
            grid[i] *= -1; // flip site i
            M += 2 * grid[i]; // update M after flip
            E += delta_energy; // update E after flip
        }
        else // Monte Carlo step
        {
            if (accept(rd) < GetExpDeltaE(grid[i], h_i, exp_dE)) // lucky accept
            {
                grid[i] *= -1; // as accept above...
                M += 2 * grid[i];
                E += delta_energy;
            }
        }
        magnet_steps.push_back(M / (double)L / (double)L); // record site magnetization
        energy_steps.push_back(E / (double)L / (double)L); // record site energy
        sum_M += M; // update w.r.t. equilibrium average
        sum_E += E;
    }
    std::cout << "J/T = " << J/T << std::endl; // progress report
    // PrintGrid(grid);
    return std::make_tuple(E,
                           sum_E / (double)steps / (double)L / (double)L,
                           energy_steps,
                           M,
                           sum_M / (double)steps / (double)L / (double)L,
                           magnet_steps);
}

int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 4 || argv[1] == "-h") // check command line arguments and give some help
    {  
        std::cerr << "Usage: " << argv[0]
                  << " steps(int): Max number of MC steps   seed(int): Initial RNG seed    cold_starts(0/1): Use cold starts"
                  << std::endl << std::endl;
        return 1;
    }
    const int steps{atoi(argv[1])}; // iax number of MC steps  
    const int seed{atoi(argv[2])}; // initial RNG seed
    const bool cold_start{atoi(argv[3])}; // 1 =: start each T with fresh all +1 spin grid - not recommended

    std::vector<double> temperatures{}; // generate vector of Ts to conduct experiments for
    for (auto t= 1.0; t < 5.1; t += 0.1)
        temperatures.push_back(t);

    std::vector<int> grid(L * L, 1); // set up grid
    double E{-J * 4 * L * L / 2}; // initial energy for all sites 
    double M{L * L}; // initial magnetization for all sites

    std::vector<std::tuple<double, double, std::vector<double>, double, std::vector<double>>> results{};
    for (const auto& T : temperatures)
    { 
        const auto result= RunExperiment(grid, E, M, T, steps, seed + T, cold_start); // run experiment for T
        E= std::get<0>(result); // to feed into next T's experiment
        const auto E_avg{std::get<1>(result)}; // for output
        const auto E_steps{std::get<2>(result)}; // for output
        M= std::get<3>(result); // to feed into next T's experiment
        const auto M_avg{std::get<4>(result)}; //  for output
        const auto M_steps{std::get<5>(result)}; // for output
        results.push_back(std::make_tuple(T, E_avg, E_steps, M_avg, M_steps)); // bundle into vector
    }

    for (auto& result : results) // Print results
        std::cout << std::get<0>(result) << "\t" <<
                        std::get<1>(result) << "\t" <<
                        std::get<3>(result) << "\t" << std::endl;
    std::cout << std::endl << std::endl; 

    return 0;
}
