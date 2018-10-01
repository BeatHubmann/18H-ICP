#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cassert>

int CongruentialRNG(int c, int p, int x_0)
{
    assert(x_0 < p && x_0 > 0); // make sure the input is smaller than p and positive
    return (c * x_0) % p;
}

int main(int argc, char* argv[])
{
     if (argc < 5) // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0] << " c p x_0 num_of_RN \n" << std::endl;
        return 1;
     }
    int c, p, n; // assign all variables
    c= atoi(argv[1]);
    p= atoi(argv[2]);
    n= atoi(argv[4]);
    int x[n+1] = {};
    x[0]= atoi(argv[3]);

    for (int i= 0; i < n; i++) // generate desired amount of RN
    {
        x[i+1]= CongruentialRNG(c, p, x[i]);
    }

    // All below: Output to file
    std::ostringstream fileNameStream_2D("rand_num_c", std::ios_base::ate);
    std::ostringstream fileNameStream_3D("rand_num_c", std::ios_base::ate);
    fileNameStream_2D << c << "_p" << p << "_2D.dat"; 
    fileNameStream_3D << c << "_p" << p << "_3D.dat"; 
    std::string fileName_2D = fileNameStream_2D.str();  
    std::string fileName_3D = fileNameStream_3D.str();  
    
    std::ofstream output_2D (fileName_2D.c_str());
    if (output_2D.is_open())
    {
        output_2D << std::fixed << std::setprecision(10);
        for (int i= 1; i < n - n%2; i++)
        {
            output_2D << (double)x[i]/(double)p << "\t" << (double)x[i+1]/(double)p
            << "\n";
        }
        output_2D.close();
    }

    std::ofstream output_3D (fileName_3D.c_str());
    if (output_3D.is_open())
    {
        output_3D << std::fixed << std::setprecision(10);
        for (int i= 1; i < n - n%3; i++)
        {
            output_3D << (double)x[i]/(double)p << "\t" << (double)x[i+1]/(double)p
            << "\t" << (double)x[i+2]/(double)p << "\n";
        }
        output_3D.close();
    }
    return 0;
}