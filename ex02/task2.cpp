#include "latticeview.h"
#include <stdio.h> // for sprintf()

#include <iostream>
#include <random>
#include <vector>

#define N 100            //Lateral number of cells
#define ImageWidth 1000  //image width
#define ImageHeight 1000 //image height
#define SEED 42


void PopulateLattice(int lattice[], const int size, const double treshold)
{
    mt19937 mt_rand(SEED);
    std::uniform_real_distribution<double> uniform_zero_one(0.0, 1.0);

    for (int ij= 0; ij < size*size; ij++)
    {
        if (uniform_zero_one(mt_rand) < treshold) lattice[ij]= 1; // set to green
        else lattice[ij]= 0; // set to white      
    }
}

void GetNeighbors(int lattice[], const int size, const int position, std::vector<int>& neighbors)
{
    const int row= position / size;
    const int col= position % size;
    auto isValid= [&size](int r, int c){ return !(
                                                    (r == -1) ||
                                                    (r == size) ||
                                                    (c == -1) ||
                                                    (c == size)
                                                 );
                                       };
    if (isValid(row-1, col)) neighbors.push_back((row-1) * N + col); // north
    if (isValid(row, col+1)) neighbors.push_back((row) * N + col+1); // east
    if (isValid(row+1, col)) neighbors.push_back((row+1) * N + col); // south
    if (isValid(row, col-1)) neighbors.push_back((row) * N + col-1); // west
}

bool TimeStep(int lattice[], const int size, const int time_step)
{
    std::cout << "Stepping " << time_step << std::endl;
    bool fire_alive= false;
    if (time_step == 1)
    {
        for (int j= 0; j < N; j++)
        {
            if (lattice[j] == 1) // if currently green == combustible
            {
                lattice[j]= 2; // set to red == on fire
                fire_alive= true;
            }
        }
    }
    else
    {
        for (int ij= 0; ij < size*size; ij++)
        {
            if (lattice[ij] == 2) // if currently red == on fire
            {
                lattice[ij]= 3; // set to black == burned out
                std::vector<int> neighbors;
                GetNeighbors(lattice, N, ij, neighbors);
                if (neighbors.size() != 0) fire_alive= true;
                for (int neighbor : neighbors)
                {
                    if (lattice[neighbor] == 1) lattice[neighbor]= 2; // set to red == on fire
                }
            }
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

    Print_lattice (lattice, N, N, ImageWidth, ImageHeight, filename);

    system(cmd);
}

int main()
{
    const double treshold= 0.5;
    int lattice[N * N];
    int time_step= 1;
    bool fire_alive= false;

    PopulateLattice(lattice, N, treshold);

    do 
    {
        fire_alive= TimeStep(lattice, N, time_step);
        OutputLattice(lattice, N, time_step);
        time_step++;
    } while (fire_alive);

    std::cout << "Time step reached: " << time_step << std::endl;

    return 0;
}