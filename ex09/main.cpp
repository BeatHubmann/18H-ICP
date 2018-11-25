#include <iostream>
#include <cmath>
#include <limits>
#include <Eigen/Dense>

class exp_func
{
    public:
        exp_func(Eigen::Vector2d init) : x_0{init[0]}, y_0{init[1]} {}
        exp_func() {}
        
        double operator() (const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return std::exp(- std::pow(x - x_0, 2) - std::pow(y - y_0, 2));
        }

        double d_dx(const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return -2 * (x - this->x_0) *  this->operator()(x_);
        }

        double d2_dx2(const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return (4 * (x - this->x_0) * (x - this->x_0) - 2) * this->operator()(x_);
        }
        double d_dy(const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return -2 * (y - this->y_0) *  this->operator()(x_);
        }

        double d2_dy2(const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return (4 * (y - this->y_0) * (y - this->y_0) - 2) * this->operator()(x_);
        }

        double d2_dxdy(const Eigen::Vector2d x_) const
        {
            double x{x_[0]}, y{x_[1]};
            return 4 * (x - this->x_0) * (y - this->y_0) * this->operator()(x_);
        }

        double d2_dydx(const Eigen::Vector2d x_) const
        {
            return this->d2_dxdy(x_);
        }
        
    private:
        double x_0{0.0}, y_0{0.0};
};

class F_func
{
    public:
        F_func(Eigen::Vector2d init) : Fx(init), Fy(init) {}

        const Eigen::Vector2d operator() (const Eigen::Vector2d x_) const
        {
            return Eigen::Vector2d(Fx.d_dx(x_), Fy.d_dy(x_));
        }

        const Eigen::Matrix2d J(const Eigen::Vector2d x_) const
        {
            Eigen::Matrix2d J;
            J << Fx.d2_dx2 (x_),   Fx.d2_dxdy(x_),
                 Fx.d2_dydx(x_),   Fx.d2_dy2 (x_);
            return J;
        }

        const Eigen::Matrix2d J_inv(const Eigen::Vector2d x_) const
        {
            return this->J(x_).inverse();
        }

        const Eigen::Matrix2d J_sec(const Eigen::Vector2d x_) const
        {
            const double h_x{x_[0] * std::sqrt(std::numeric_limits<double>::epsilon())};
            const double h_y{x_[1] * std::sqrt(std::numeric_limits<double>::epsilon())};
            const Eigen::Vector2d e_x{1.0, 0.0};
            const Eigen::Vector2d e_y{0.0, 1.0};

            Eigen::Matrix2d J_sec;
            J_sec << (this->operator()(x_ + h_x * e_x)[0] - this->operator()(x_)[0]) / h_x,
                     (this->operator()(x_ + h_y * e_y)[0] - this->operator()(x_)[0]) / h_y,
                     (this->operator()(x_ + h_x * e_x)[1] - this->operator()(x_)[1]) / h_x,
                     (this->operator()(x_ + h_y * e_y)[1] - this->operator()(x_)[1]) / h_y;
            return J_sec;
        }

        const Eigen::Matrix2d J_sec_inv(const Eigen::Vector2d x_) const
        {
            return this->J_sec(x_).inverse();
        }

    private:
        exp_func Fx, Fy;
};

inline double dist(const Eigen::Vector2d a, 
                   const Eigen::Vector2d b = Eigen::Vector2d::Zero(2))
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]));
}

int main(int argc, char* argv[])
{
    if (argc < 5 || argc > 5 || argv[1] == "-h") // check cl args and give some help
    {  
        std::cerr << "Usage: " << argv[0] << "\n\t"
                  << " x_0(double): parameter x_0" << "\n\t"
                  << " y_0(double): parameter y_0" << "\n\t"
                  << " x_init(double): initial x" << "\n\t"
                  << " y_init(double): initial y" << "\n\t"
                  << std::endl << std::endl;
        return 1;
    }
    const double x_0{atof(argv[1])};
    const double y_0{atof(argv[2])};
    const Eigen::Vector2d init(x_0, y_0);

    const double x_start{atof(argv[3])};
    const double y_start{atof(argv[4])};
    const Eigen::Vector2d start(x_start, y_start);

    exp_func f(init);
    F_func F(init);

    // Newton-Raphson method
    Eigen::Vector2d x_n;
    Eigen::Vector2d x_np1{start - F.J_inv(start) * F(start)};

    std::cout << start[0] << "\t" << start[1] << "\t" << f(start) << std::endl;

    while (dist(F(x_np1)) > std::numeric_limits<double>::epsilon())
    {
        x_n= x_np1;
        std::cout << x_n[0] << "\t" << x_n[1] << "\t" << f(x_n) << std::endl;
        x_np1= x_n - F.J_inv(x_n) * F(x_n);
    }

    std::cout << std::endl << std::endl;

    // Secant method
    x_np1= start - F.J_sec_inv(start) * F(start); 

    std::cout << start[0] << "\t" << start[1] << "\t" << f(start) << std::endl;

    while (dist(F(x_np1)) > std::numeric_limits<double>::epsilon())
    {
        x_n= x_np1;
        std::cout << x_n[0] << "\t" << x_n[1] << "\t" << f(x_n) << std::endl;
        x_np1= x_n - F.J_sec_inv(x_n) * F(x_n);
    }

    return 0;
}