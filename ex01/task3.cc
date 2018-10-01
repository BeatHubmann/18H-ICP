#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cassert>

int CongruentialRNG(int c, int p, int x_0) // all as in task1
{
    assert(x_0 < p && x_0 > 0);
    return (c * x_0) % p;
}

int main(int argc, char* argv[])
{
     if (argc < 6)
     {  
        std::cerr << "Usage: " << argv[0] << " c p x_0 num_of_RN num_of_bins\n" << std::endl;
        return 1;
     }
    int c, p, n, k;
    c= atoi(argv[1]);
    p= atoi(argv[2]);
    n= atoi(argv[4]);
    k= atoi(argv[5]);
    int x[n+1]= {};
    x[0]= atoi(argv[3]);


    const double p_i= 1. / k; // hypothesis: probability of RN to be in specific bin
    const double bin_width = (double)p / (double)k;

    std::cout << "bin width " << bin_width << std::endl;
    double chi_squared_total= 0;

    const int num_of_runs= 10;

    for (int j= 0; j < num_of_runs; j++) // run several times to average out chi score
    {
        double chi_squared= 0.;
        int N[k]= {};

        for (int i= 0; i < n; i++)
        {
            x[i+1]= CongruentialRNG(c, p, x[i]);
            int bin= x[i+1] / bin_width; // cheap way to get a histogram
            N[bin]++;
        }

        for (int i= 0; i < k; i++) // calculate chi score for current run
        {
            double enumerator= N[i] - n * p_i; 
            chi_squared += (enumerator * enumerator) / (n * p_i);
        }
        std::cout << "For c= " << c << " , p= " << p << " , x_0= " << x[0] << " and\n"
                << n << " RN in " << k << " bins:\n"
                << "chi squared= " << chi_squared << std::endl;
        chi_squared_total += chi_squared;
        x[0] = (x[0] + x[1]) % p;   
    }
    std::cout << num_of_runs << " runs: average chi squared= " << chi_squared_total / num_of_runs << std::endl;
    return 0;
}