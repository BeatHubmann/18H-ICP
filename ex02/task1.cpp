#include "latticeview.h"
#include <stdio.h> // for sprintf()

#include <iostream>
#include <random>

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
        if (uniform_zero_one(mt_rand) < treshold) lattice[ij]= 1;
        else lattice[ij]= 0;        
    }
}

int main()
{
    const double treshold= 0.5;
    int lattice[N * N];
    char cmd[160], filename[160], filename_png[160];

    PopulateLattice(lattice, N, treshold);
    
    sprintf(filename, "task1.ppm");
    sprintf(filename_png, "task1.png");
    sprintf(cmd, "convert %s %s; rm -f %s", filename, filename_png, filename); 

    Print_lattice (lattice, N, N, ImageWidth, ImageHeight, filename);

    system(cmd);

    return 0;
}