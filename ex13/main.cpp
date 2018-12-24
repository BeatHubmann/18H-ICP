#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

using Real= double;
using Vector= std::vector<Real>;

template <typename F>
class Waveform
{
    public:
        Waveform(const F& f, const Real c= 1.0, const Real dt= 0.01, const Real dx= 1.0)
                : c(c), dt(dt), dx(dx), num_points(100/dx + 1), b(c * c * dt * dt / (dx * dx))
        {
            std::cout << "Waveform with c= " << c << ", "
                      << "dt= " << dt << ", "
                      << "dx= " << dx << " -> "
                      << "b= " << b << std::endl
                      << "Number of mesh points: " << num_points << std::endl;
            u_previous.resize(num_points + 2); // +2 ghosts for periodic boundary conditions
            u_current.resize(num_points + 2); //  +2 ghosts for periodic boundary conditions
            u_next.resize(num_points + 2); // +2 ghosts for periodic boundary conditions
            for (int i= 1; i <= num_points; i++)
            {
                const Real x{(i - 1) * dx};
                u_previous[i]= f(x, -dt);
                u_current[i]= f(x, 0);
            }
            u_previous[0]= u_previous[num_points]; // periodic boundary conditions
            u_previous[num_points + 1]= u_previous[1];  // periodic boundary conditions
            u_current[0]= u_current[num_points];  // periodic boundary conditions
            u_current[num_points + 1]= u_current[1];  // periodic boundary conditions
            this->write();
        }

        void advance()
        {
            for (int x= 1; x <= num_points; x++)
            {
                u_next[x]= 2 * (1 - b) * u_current[x]
                         + b * (u_current[x + 1] + u_current[x - 1])
                         - u_previous[x];
            }
            u_next[0]= u_next[num_points]; // periodic boundary conditions
            u_next[num_points + 1]= u_next[1]; // periodic boundary conditions
            u_previous= u_current;
            u_current= u_next;
        }

        void write(const int timestep= 0)
        {
            std::ostringstream file_name;
            file_name << "wave_" << std::setw(5) << std::setfill('0') << timestep << ".dat";
            std::ofstream file;
            file.open(file_name.str().c_str(), std::ios::out | std::ios::app);
            for (int x= 0; x < num_points; x++)
            {
                file << x * dx << "\t" << u_current[x] << std::endl;
            }
            file << std::endl;
            file.close();
        }

    private:
        Vector u_previous;
        Vector u_current;
        Vector u_next;
        const Real c;
        const Real dt;
        const Real dx;
        const Real b;
        const int num_points;
};

class ExpFunc
{
    public:
        ExpFunc(Real c= 0) : c(c) {}

        Real operator() (Real x, Real t) const
        {
            return std::exp(-(x - c * t - 10) * (x - c * t - 10));
        }

    private:
        const Real c;
};

class SinFunc
{
    public:
        SinFunc(Real c) : c(c) {}

        Real operator() (Real x, Real t) const
        {
            return std::sin(x - c * t);
        }
    private:
        const Real c;
};

int main(int argc, char* argv[])
{
     if (argc < 5 || argc > 5 || argv[1] == "-h") // check cl args and give some help
    {  
        std::cerr << "Usage: " << argv[0] << "\n\t"
                  << " dx(Real): mesh parameter: mesh spacing" << "\n\t"
                  << " dt(Real): time step size" << "\n\t"
                  << " c(Real): Wave coefficient" << "\n\t"
                  << " N(int): Number of timesteps" << "\n\t"
                  << std::endl << std::endl;
        return 1;
    }

    const Real dx{atof(argv[1])};
    const Real dt{atof(argv[2])};
    const Real c{atof(argv[3])};
    const int N{atoi(argv[4])};

    const Real b{c * c * dt * dt / (dx * dx)};

    ExpFunc f;
    // SinFunc f(c);

    Waveform wave(f, c, dt, dx);

    for (int n= 1; n <= N; n++)
    {
        wave.advance();
        wave.write(n);
    }

    return 0;
}