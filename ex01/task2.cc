#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cassert>

int CongruentialRNG(int c, int p, int x_0)
{
    assert(x_0 < p && x_0 > 0);
    return (c * x_0) % p;
}

int main(int argc, char* argv[])
{
     if (argc < 5)
     {  
        std::cerr << "Usage: " << argv[0] << " c p x_0 num_of_RN \n" << std::endl;
        return 1;
     }
    int c, p, n;
    c= atoi(argv[1]);
    p= atoi(argv[2]);
    n= atoi(argv[4]);
    int x[2*n+1] = {};
    x[0]= atoi(argv[3]);

    for (int i= 0; i < 2*n; i++)
    {
        x[i+1]= CongruentialRNG(c, p, x[i]);
    }

    std::ostringstream fileNameStream_2D("rand_num_c", std::ios_base::ate);
    fileNameStream_2D << c << "_p" << p << "_2D_circle.dat"; 
    std::string fileName_2D = fileNameStream_2D.str();  
    
    std::ofstream output_2D (fileName_2D.c_str());
    if (output_2D.is_open())
    {
        output_2D << std::fixed << std::setprecision(10); 
        for (int i= 1; i < n + 1; i++)
        {
            double r= sqrt((double)x[i]/(double)p);
            double theta= 2 * M_PI * (double)x[i+n]/(double)p;
            output_2D << r * cos(theta) << "\t" << r * sin(theta) << "\n";
        }
        output_2D.close();
    }
    return 0;
}