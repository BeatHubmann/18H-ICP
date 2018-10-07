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

    sprintf(filename, "task2_%03d.ppm", time_step);
    sprintf(filename_png, "task2_%03d.png", time_step);
    sprintf(cmd, "convert %s %s; rm -f %s", filename, filename_png, filename);  

    Print_lattice (lattice, size, size, ImageWidth, ImageHeight, filename);

    system(cmd);
}

std::tuple<bool, bool> TimeStep(int lattice[], const int size, const int time_step, const bool generate_images)
{
    bool fire_alive= false;
    bool reached_opposite= false;
    if (time_step % 10 == 0) std::cout << "..." <<std::endl;
    int image_lattice[size * size]= {0};

    #pragma omp parallel for 
    for (int ij= 0; ij < size*size; ij++)
    {
        if (lattice[ij] == time_step-1)
        {
            image_lattice[ij]= 3;

            std::vector<int> neighbors;
            GetNeighbors(lattice, size, ij, neighbors);
            if (neighbors.size() != 0)
                fire_alive= true;
            for (int neighbor : neighbors)
            {
                if (neighbor / size == (size-1)) reached_opposite= true;
                lattice[neighbor]= time_step;
                image_lattice[neighbor]= 2;
            }
        }
        else if (lattice[ij] == 1) image_lattice[ij]= 1;
        else if (lattice[ij] < time_step-1 && lattice[ij] > 1) image_lattice[ij]= 3;
    }

    if (generate_images) OutputLattice(image_lattice, size, time_step);

    return std::make_tuple(fire_alive, reached_opposite); 
}

int main(int argc, char* argv[])
{
    if (argc < 5) // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0] << " float_site_occupation_probabilty int_lattice_side_length int_random_seed bool_print_images\n" << std::endl;
        return 1;
     }

    const double treshold= atof(argv[1]);
    const int size= atoi(argv[2]);
    const int seed= atoi(argv[3]);
    const bool generate_images= atoi(argv[4]);
    int lattice[size * size]= {0};

    PopulateLattice(lattice, size, treshold, seed);

    bool fire_alive= false;
    bool reached_opposite= false;

    int time_step= 2;
    int shortest_path= std::numeric_limits<int>::max();
    std::cout << std::left 
              << std::setw(30) << "Lattice side length: " << std::setw(15) << size << std::endl
              << std::setw(30) << "Occupation probability: " << std::setw(15) << treshold << std::endl
              << std::setw(30) << "Random seed: " << std::setw(15) << seed << std::endl
              << std::setw(30) << "Calculating..." << std::endl;

    fire_alive= StartFire(lattice, size);
    if (generate_images) Print_lattice (lattice, size, size, ImageWidth, ImageHeight, "task2_002.png");    

    while (fire_alive)
    {
        time_step++;
        auto step_result= TimeStep(lattice, size, time_step, generate_images);
        fire_alive= std::get<0>(step_result);
        reached_opposite= std::get<1>(step_result);
        if (reached_opposite) shortest_path= std::min(time_step-1, shortest_path);
    }

    std::cout << std::left
              << std::setw(30) << "Fire duration:" << std::setw(15) << time_step-2 << std::endl
              << std::setw(30) << "Shortest path:" << std::setw(15) << (shortest_path < std::numeric_limits<int>::max() ? shortest_path : -1) << std::endl;
    return 0;
}