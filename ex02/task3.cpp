#include "latticeview_v2.h"
#include <stdio.h> // for sprintf()

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <vector>
#include <tuple>
#include <limits>
#include <omp.h>

#define ImageWidth 1000  //image width
#define ImageHeight 1000 //image height

void PopulateLattice(int lattice[], const int size, const double treshold, const int seed)
{
    std::mt19937 mt_rand(seed);
    std::uniform_real_distribution<double> uniform_zero_one(0.0, 1.0);

    for (int ij= 0; ij < size*size; ij++)
    {
        if (uniform_zero_one(mt_rand) < treshold) lattice[ij]= 1; // set to green
    }
}

void GetNeighbors(int lattice[], const int size, const int position, std::vector<int>& neighbors)
{
    const int row= position / size;
    const int col= position % size;
    auto isValid= [&lattice, &size](int r, int c){ return !(
                                                    (r == -1) ||
                                                    (r == size) ||
                                                    (c == -1) ||
                                                    (c == size)
                                                 ) &&
                                                 (
                                                     lattice[r*size + c] == 1
                                                 );
                                       };
    if (isValid(row-1, col)) neighbors.push_back((row-1) * size + col); // north
    if (isValid(row, col+1)) neighbors.push_back((row) * size + col+1); // east
    if (isValid(row+1, col)) neighbors.push_back((row+1) * size + col); // south
    if (isValid(row, col-1)) neighbors.push_back((row) * size + col-1); // west
}

bool StartFire(int lattice[], const int size)
{
    bool fire_alive= false;
    for (int j= 0; j < size; j++) // first row only
    {
        if (lattice[j] == 1) 
        {
            lattice[j]= 2; 
            fire_alive= true;
        }
    }
    return fire_alive;
}

void OutputLattice(int lattice[], const int size, const int time_step)
{
    char cmd[160], filename[160], filename_png[160];

    sprintf(filename, "task3_%03d.ppm", time_step);
    sprintf(filename_png, "task3_%03d.png", time_step);
    sprintf(cmd, "convert %s %s; rm -f %s", filename, filename_png, filename);  

    Print_lattice (lattice, size, size, ImageWidth, ImageHeight, filename);

    system(cmd);
}

std::tuple<bool, bool> TimeStep(int lattice[], const int size, const int time_step)
{
    bool fire_alive= false;
    bool reached_opposite= false;

    #pragma omp parallel for 
    for (int ij= 0; ij < size*size; ij++)
    {
        if (lattice[ij] == time_step-1)
        {
            std::vector<int> neighbors;
            GetNeighbors(lattice, size, ij, neighbors);
            if (neighbors.size() != 0)
                fire_alive= true;
            for (int neighbor : neighbors)
            {
                if (neighbor / size == (size-1)) reached_opposite= true;
                lattice[neighbor]= time_step;
            }
        }
    }
    return std::make_tuple(fire_alive, reached_opposite); 
}


auto RunExperiment(const int size, const double treshold, const int seed)
{
    //std::cout << "Running exp " << seed << "\t" << treshold << std::endl;
    int lattice[size * size]{0};

    PopulateLattice(lattice, size, treshold, seed);

    bool fire_alive= false;

    int time_step= 2;
    int shortest_path= 0;
    bool spanning_cluster= false;

    fire_alive= StartFire(lattice, size);

    while (fire_alive)
    {
        time_step++;
        auto step_result= TimeStep(lattice, size, time_step);
        fire_alive= std::get<0>(step_result);
        if (!(spanning_cluster) && std::get<1>(step_result))
        {
            shortest_path= time_step -1;
            spanning_cluster= true;
        }
    }
    const int fire_duration= time_step - 2;

    return std::make_tuple((int)spanning_cluster, fire_duration, shortest_path);
}




int main(int argc, char* argv[])
{
    if (argc < 5) // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0]
                  << " int_occupation_prob_steps int_lattice_side_length int_initial_seed int_num_samples \n" << std::endl;
        return 1;
     }
    const int treshold_steps= atoi(argv[1]);
    const int size= atoi(argv[2]);
    const int initial_seed= atoi(argv[3]);
    const int num_samples= atoi(argv[4]);

    const double treshold_step_size= 1.0 / treshold_steps;
    std::vector<double> tresholds{0.0};
    for (auto i= 1; i <= treshold_steps; i++)
        tresholds.push_back(i * treshold_step_size);
    
    for (auto treshold : tresholds)
    {
        int num_spanning_clusters= 0;
        long int total_fire_duration= 0;
        long int total_shortest_path= 0;
        for (auto n= 0; n < num_samples; n++)
        {
            auto result= RunExperiment(size, treshold, initial_seed + n);
            num_spanning_clusters += std::get<0>(result);
            total_fire_duration += std::get<1>(result);
            total_shortest_path += std::get<2>(result);
        }
        const double percolation_rate= num_spanning_clusters / (float)num_samples;
        const double avg_fire_duration= total_fire_duration / (float)num_samples;
        const double avg_shortest_path= (num_spanning_clusters < 1 ? 0.0 :
                                         total_shortest_path / (float)num_spanning_clusters);

        std::cout << treshold << "\t"
                  << percolation_rate << "\t"
                  << avg_fire_duration << "\t"
                  << avg_shortest_path << std::endl;
    }
    return 0;
}


