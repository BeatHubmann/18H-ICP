#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <array>
#include <iterator>

using state_type= std::vector<double>;

template<class state_type> // for generic containers
void resize(const state_type& in, state_type& out)
{
    using std::size;
    out.resize(size(in));
}

template<class T, std::size_t N> // specialize for arrays
void resize(const std::array<T, N>&, std::array<T, N>&) {} // catch arrays: no resizing needed

struct ContainerIterator
{
    template<class S_1,
             class S_2,
             class S_3,
             class Op>
    void ForEach_3(S_1& s_1, S_2& s_2, S_3& s_3,
                   Op op) const
    {
        using std::begin;
        using std::end;

        auto first_1{begin(s_1)};
        auto last_1{end(s_1)};
        auto first_2{begin(s_2)};
        auto first_3{begin(s_3)};

        for ( ; first_1 != last_1; )
            op(*first_1++, *first_2++, *first_3++);
    }

    template<class S_1,
             class S_2,
             class S_3,
             class S_4, 
             class S_5,
             class S_6,
             class Op>
    void ForEach_6(S_1& s_1, S_2& s_2, S_3& s_3,
                    S_4& s_4, S_5& s_5, S_6& s_6,
                    Op op) const
    {
        using std::begin;
        using std::end;

        auto first_1{begin(s_1)};
        auto last_1{end(s_1)};
        auto first_2{begin(s_2)};
        auto first_3{begin(s_3)};
        auto first_4{begin(s_4)};
        auto first_5{begin(s_5)};
        auto first_6{begin(s_6)};

        for ( ; first_1 != last_1; )
            op(*first_1++, *first_2++, *first_3++,
               *first_4++, *first_5++, *first_6++);
    }
};

struct ContainerOps
{
    template<class F_1 = double,
             class F_2 = F_1>
    struct ScaleSum_2
    {
        ScaleSum_2(F_1 a_1, F_2 a_2) : alpha_1(a_1), alpha_2(a_2) {}

        template<class T_0, class T_1, class T_2>
        void operator() (T_0& t_0, const T_1& t_1, const T_2& t_2) const
        {
            t_0= alpha_1 * t_1
               + alpha_2 * t_2;
        }
        
        using result_type= void;

        const F_1 alpha_1;
        const F_2 alpha_2;
    };

    template<class F_0 = double,
             class F_1 = F_0,
             class F_2 = F_1,
             class F_3 = F_2,
             class F_4 = F_3,
             class F_5 = F_4>
    struct ScaleSum_5
    {
        ScaleSum_5(F_1 a_1, F_2 a_2, F_3 a_3, F_4 a_4, F_5 a_5)
            : alpha_1(a_1), alpha_2(a_2), alpha_3(a_3), alpha_4(a_4), alpha_5(a_5) {}
        
        template<class T_0, class T_1, class T_2, class T_3, class T_4, class T_5>
        void operator() (T_0& t_0, const T_1& t_1, const T_2& t_2, const T_3& t_3,
                         const T_4& t_4, const T_5& t_5) const
        {
            t_0= alpha_1 * t_1
               + alpha_2 * t_2
               + alpha_3 * t_3
               + alpha_4 * t_4
               + alpha_5 * t_5;
        }

        using result_type= void;

        const F_1 alpha_1;
        const F_2 alpha_2;
        const F_3 alpha_3;
        const F_4 alpha_4;
        const F_5 alpha_5;
    };
};

template<class state_type,
         class value_type = double,
         class time_type = value_type,
         class iterator = ContainerIterator,
         class ops = ContainerOps>
class RungeKutta4
{
    public:
        template<typename System>
        void do_step(System& system,
                     state_type& x,
                     time_type t,
                     time_type dt)
        {
            AdjustSize(x);
            const value_type one{(value_type)1};
            const time_type dt_2{dt / 2},
                            dt_3{dt / 3},
                            dt_6{dt / 6};

            using ScaleSum_2= typename ops::template ScaleSum_2<value_type,
                                                                time_type> ;
            using ScaleSum_5= typename ops::template ScaleSum_5<value_type,
                                                                time_type,
                                                                time_type,
                                                                time_type,
                                                                time_type,
                                                                time_type>;  

            system(x,     k_1, t);
            my_iterator.ForEach_3(x_tmp, x, k_1, ScaleSum_2(one, dt_2));

            system(x_tmp, k_2, t + dt_2);
            my_iterator.ForEach_3(x_tmp, x, k_2, ScaleSum_2(one, dt_2));

            system(x_tmp, k_3, t + dt_2);
            my_iterator.ForEach_3(x_tmp, x, k_2, ScaleSum_2(one, dt_2));

            system(x_tmp, k_4, t + dt);
            my_iterator.ForEach_6(x, x, k_1, k_2, k_3, k_4,
                                                 ScaleSum_5(one, dt_6, dt_3,
                                                            dt_3, dt_6));
        }
    private:
        state_type x_tmp,
                   k_1,
                   k_2,
                   k_3,
                   k_4;
        
        iterator my_iterator;

        void AdjustSize(const state_type& x)
        {
            resize(x, x_tmp);
            resize(x, k_1);
            resize(x, k_2);
            resize(x, k_3);
            resize(x, k_4);
        }
};

struct launch
{
    launch(const double gamma,
           const double g) : gamma(gamma), g(g) {}

    void operator() (const state_type& x, state_type& dx_dt, double t)
    {
        dx_dt[0]=     - gamma * x[0]; // x coordinate
        dx_dt[1]= - g - gamma * x[1]; // z coordinate
    }

    const double gamma, g;
};

struct track
{
    track(const double v_0,
          const double v_1) : v_0(v_0), v_1(v_1) {}

    void operator() (const state_type& x, state_type& dx_dt, double t)
    {
        dx_dt[0]= v_0; // x coordinate
        dx_dt[1]= v_1; // z coordinate
    }

    const double v_0, v_1;
};

int main(int argc, char* argv[])
{
    using rk4_type= RungeKutta4<state_type>;

    if (argc < 7 || argc > 7 || argv[1] == "-h") // check cl args and give some help
    {  
        std::cerr << "Usage: " << argv[0] << "\n\t"
                  << " alpha(double): initial angle alpha in degrees" << "\n\t"
                  << " v_init(double): initial velocity v" << "\n\t"
                  << " gamma(double): drag coefficient gamma" << "\n\t"
                  << " g(double): gravitational accel g" << "\n\t"
                  << " dt(double): length of time step delta t" << "\n\t"
                  << " num_steps(int): number of time steps num_steps"
                  << std::endl << std::endl;
        return 1;
    }

    const double alpha{atof(argv[1])};
    const double v_init{atof(argv[2])};
    const double gamma{atof(argv[3])};
    const double g{atof(argv[4])};
    const double dt{atof(argv[5])};
    const int num_steps{atoi(argv[6])};

    std::vector<double> alphas, gammas, vees;

    for (auto a= 10; a < alpha; a += 10)
        alphas.push_back(a);
    alphas.push_back(alpha);
    alphas.push_back(45);

    for (auto gamma= 0.25; gamma <= 10.0; gamma += 0.25)
        gammas.push_back(gamma);

    for (auto v= 5; v < v_init; v *= 2)
        vees.push_back(v);
    vees.push_back(v_init);

    // Task 1:
    for (auto vee : vees)
        for (auto alpha : alphas)
        {
            std::cout << "alpha=" << alpha << ", "
                      << "v=" << vee << ", "
                      << "gamma=" << gamma << ", "
                      << "g=" << g << std::endl;

            state_type v{std::cos(alpha / 180 * M_PI) * vee,
                         std::sin(alpha / 180 * M_PI) * vee}; // IC velocity
            state_type x{0.0, 0.0}; // IC position
            
            rk4_type stepper;
            launch projectile_velocity{gamma, g};

            for (auto n= 0; n <= num_steps; n++)
            {
                std::cout << n * dt << "\t"
                        << v[0] << "\t"
                        << v[1] << "\t"
                        << x[0] << "\t"
                        << x[1] << std::endl;
                stepper.do_step(projectile_velocity, v, n * dt, dt);
                track projectile_position(v[0], v[1]);
                stepper.do_step(projectile_position, x, n * dt, dt);
                if (x[1] < 0) break;
            }
            std::cout << std::endl << std::endl;
        }
    
    // Task 2:
    for (auto gamma : gammas)
    {   
        double alpha_max{0}, x_max{0.0};

        for (double a= 1.0; a < 90.0; a += 0.5)
        {
            state_type v{std::cos(a / 180 * M_PI) * v_init,
                         std::sin(a / 180 * M_PI) * v_init}; // IC velocity
            state_type x{0.0, 0.0}; // IC position
            
            rk4_type stepper;
            launch projectile_velocity{gamma, g};

            for (auto n= 0; n <= num_steps; n++)
            {
                stepper.do_step(projectile_velocity, v, n * dt, dt);
                track projectile_position(v[0], v[1]);
                stepper.do_step(projectile_position, x, n * dt, dt);
                if (x[1] < 0) break; // consider z=0 as ground level
            }
            if (x[0] > x_max)
            {
                x_max= x[0];
                alpha_max= a;
            }
        }
        std::cout << gamma << "\t" << alpha_max << "\t" << x_max << std::endl;
    }
    return 0;
}