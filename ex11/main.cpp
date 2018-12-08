#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cmath>
#include <cassert>
#include <chrono>
#include "writer.hpp"

using namespace std::chrono;
using Real= double;
using SparseMatrix= Eigen::SparseMatrix<Real>;
using Triplet= Eigen::Triplet<Real>;
using Vector= Eigen::VectorXd;

template <typename T>
class Rho
{
    public:
        Rho(const T x_0, const T y_0, const T dx)
         : x_0(x_0), y_0(y_0), dx(dx) {}
        T operator() (const T x, const T y)
        {   
           if (!set && abs(x - x_0) < dx / 2 && abs(y - y_0) < dx / 2)
           {
               set= true;
               return (T)1.0 / dx / dx;
           }
           else
           {
               return (T)0.0;
           }
        }
    private:
        const T x_0, y_0, dx;
        bool set{false};
};

void CreatePoissonMatrix(SparseMatrix& A, const int N)
{
    A.resize(N * N, N * N); // resize for unrolling RHS vector
    std::vector<Triplet> triplets; // to build sparse matrix A from
    triplets.reserve(N * N * 5); // each row of A in general has 5 entries - will be somewhat too large
    for (int i= 0; i < N * N; i++)
    {
        triplets.push_back(Triplet(i, i, 4)); // main diagonal
        if (i > 0 && i % N != 0)
            triplets.push_back(Triplet(i, i - 1, -1)); // left off-diagonal: i-1
        if (i < N * N - 1 && (i + 1) % N != 0 || i == 0) 
            triplets.push_back(Triplet(i, i + 1, -1)); // right off-diagonal: i+1
        if (i < N * N - N)
            triplets.push_back(Triplet(i, i + N, -1)); // j+1
        if (i > N - 1)
            triplets.push_back(Triplet(i, i - N, -1)); // j-1
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
}

template <class F, typename T>
void CreateRHSVector(Vector& RHS, const int N, const T dx, F& rho)
{
    RHS.resize(N * N); // unroll grid into vector
    
    auto x= [&dx](T i) {return (i + 1) * dx;};    
    auto y= [&dx](T j) {return (j + 1) * dx;};    

    const T dx2{dx * dx};

    for (int j= 0; j < N; j++) // row major storage
        for (int i= 0; i < N; i++)
            RHS[j * N + i]= dx2 * rho(x(i), y(j));
}

template <typename T>
Vector SolveJacobi(const SparseMatrix& A, const Vector& RHS, const int N,
                   const T treshold, const int iterations)
{
    Vector x= Vector::Zero((N + 2) * (N + 2));
    Vector x_0= Vector::Ones(N * N);
    Vector x_1;

    int i{0};
    T l2_norm{1.0};

    high_resolution_clock::time_point t_0= high_resolution_clock::now();
    const SparseMatrix D(A.diagonal().asDiagonal());
    const SparseMatrix D_inv(D.diagonal().asDiagonal().inverse());
    const SparseMatrix R= A - D;

    const Vector D_inv_x_RHS= D_inv * RHS; // to speed up calculations during loop
    const SparseMatrix D_inv_x_R= D_inv * R;

    while (l2_norm > treshold && i++ < iterations)
    {
        x_1= D_inv_x_RHS - D_inv_x_R * x_0; // x_1= D_inv * (RHS - R * x_0)
        l2_norm= (x_1 - x_0).norm();
        x_0= x_1;
    }
    high_resolution_clock::time_point t_1= high_resolution_clock::now();
    std::cout << "Jacobi relaxation solver: Success" 
              << " after " << i << " iterations in "
              << duration_cast<milliseconds>(t_1 - t_0).count() <<" ms" << std::endl;

    for (int j= 1; j <= N; j++)
        for (int i= 1; i <= N; i++)
                x[j * (N + 2) + i]= x_1[(j - 1) * N + (i - 1)];
    return x;
}

template <typename T>
Vector SolveGaussSeidel(const SparseMatrix& A, const Vector& RHS, const int N,
                        const T treshold, const int iterations)
{
    Vector x= Vector::Zero((N + 2) * (N + 2));
    Vector x_0= Vector::Ones(N * N);
    Vector x_1;

    int i{0};
    T l2_norm{1.0};
    
    high_resolution_clock::time_point t_0= high_resolution_clock::now();
    const SparseMatrix L= A.triangularView<Eigen::Lower>();
    const SparseMatrix U= A.triangularView<Eigen::StrictlyUpper>();

    Eigen::SparseLU<SparseMatrix> solver; // to get L^{-1} 
    solver.compute(L);
    SparseMatrix I(N * N, N * N); I.setIdentity();
    const SparseMatrix L_inv= solver.solve(I);

    const Vector L_inv_x_RHS= L_inv * RHS; // to speed up calculations during loop
    const SparseMatrix L_inv_x_U= L_inv * U;

    while (l2_norm > treshold && i++ < iterations)
    {
        x_1= L_inv_x_RHS - L_inv_x_U * x_0; // x_1= L_inv * (RHS - U * x_0)
        l2_norm= (x_1 - x_0).norm();
        x_0= x_1;
    }
    high_resolution_clock::time_point t_1= high_resolution_clock::now();
    std::cout << "Gauss-Seidel solver: Success" 
              << " after " << i << " iterations in " 
              << duration_cast<milliseconds>(t_1 - t_0).count() <<" ms" << std::endl;

    for (int j= 1; j <= N; j++)
        for (int i= 1; i <= N; i++)
                x[j * (N + 2) + i]= x_1[(j - 1) * N + (i - 1)];
    return x;
}

Vector SolveCholesky(const SparseMatrix& A, const Vector& RHS, const int N)
{
    Vector x= Vector::Zero((N + 2) * (N + 2)); // to fill with solution
    
    high_resolution_clock::time_point t_0= high_resolution_clock::now();
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(A);
    assert (solver.info() == Eigen::Success);
    Vector x_inner= solver.solve(RHS); // solve inner grid
    high_resolution_clock::time_point t_1= high_resolution_clock::now();
    std::cout << "Cholesky solver: Success in "
              << duration_cast<milliseconds>(t_1 - t_0).count() <<" ms" << std::endl;
    for (int j= 1; j <= N; j++)
        for (int i= 1; i <= N; i++)
                x[j * (N + 2) + i]= x_inner[(j - 1) * N + (i - 1)];
    return x;
}

int main(int argc, char* argv[])
{
     if (argc < 4 || argc > 4 || argv[1] == "-h") // check cl args and give some help
    {  
        std::cerr << "Usage: " << argv[0] << "\n\t"
                  << " N(int): grid parameter: number of points per axis" << "\n\t"
                  << " x_0(Real): x coord of point charge in 0..1" << "\n\t"
                  << " y_0(Real): y coord of point charge in 0..1" << "\n\t"
                  << std::endl << std::endl;
        return 1;
    }

    const int N{atoi(argv[1])};
    const Real dx{(Real)1.0 / (Real)(N + 1)};
    const Real x_0{atof(argv[2])}, y_0{atof(argv[3])};

    const Real treshold{0.0001};
    const int iterations{10000};

    SparseMatrix A;
    Vector RHS;
    Vector RHS2;
    Rho rho{x_0, y_0, dx};
    Rho rho2{1.0 - x_0, 1.0 - y_0, dx};
    CreatePoissonMatrix(A, N);
    CreateRHSVector(RHS, N, dx, rho);
    CreateRHSVector(RHS2, N, dx, rho2);
    RHS = (RHS + RHS2) / 2;
    Vector x_cholesky= SolveCholesky(A, RHS, N);
    Vector x_jacobi= SolveJacobi(A, RHS, N, treshold, iterations);
    Vector x_gauss_seidel= SolveGaussSeidel(A, RHS, N, treshold, iterations);
    std::cout << "L2 difference Cholesky - Jacobi: " << (x_cholesky - x_jacobi).norm() << std::endl
              << "L2 difference Cholseky - Gauss-Seidel: " << (x_cholesky - x_jacobi).norm() << std::endl;
    writeToFile("x_cholesky.txt", x_cholesky);
    writeToFile("x_jacobi.txt", x_jacobi);
    writeToFile("x_gauss_seidel.txt", x_gauss_seidel);
    writeToFile("x_diff_chol_jaco.txt", (x_cholesky - x_jacobi));
    writeToFile("x_diff_chol_gaus.txt", (x_cholesky - x_gauss_seidel));
    return 0;
}