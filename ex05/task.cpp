#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <tuple>



int main(int argc, char* argv[])
{
    if (argc < 3 || argc > 3 || argv[1] == "-h") // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0]
                  << " int_lattice_side_length int_initial_seed \n" << std::endl;
        return 1;
     }
    const int size= atoi(argv[1]);
    const int initial_seed= atoi(argv[2]);

    std::vector<double> tresholds{0.58, 0.592746, 0.61}; // Init vector of probs and make sure p_c is included

    for (auto& treshold : tresholds) // Run experiments on all chosen occupation probs
    {
        std::cout << "Experiment_1:Sandbox_Algorithm" << std::endl;

        std::vector<std::tuple<int, int>> results= RunExperiment1(size, treshold, initial_seed);

        std::cout << "p=" << treshold << std::endl;
       
        for (auto& result : results)
            std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
        std::cout << "\n\n";
    }

    for (auto& treshold : tresholds)
    {
        std::cout << "Experiment_2:Box_Counting_Algorithm" << std::endl;

        std::vector<std::tuple<int, int>> results= RunExperiment2(size, treshold, initial_seed);

        std::cout << "p=" << treshold << std::endl;
       
        for (auto& result : results)
            std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
        std::cout << "\n\n";
    }

    return 0;
}
