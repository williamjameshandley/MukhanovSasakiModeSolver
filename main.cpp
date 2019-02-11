#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main(int argc, char** argv)
{
    if (argc != 11) return -1;
    std::cout.precision(18);
    //Choose Potential and set potential_ptr which is passed to background solver
    double mass = exp(atof(argv[3])) * 1e-6;
    double lambda = atof(argv[4]);
    double custom5 = atof(argv[5]);
    double custom6 = atof(argv[6]);
    double custom7 = atof(argv[7]);
    double custom8 = atof(argv[8]);
    double custom9 = atof(argv[9]);
    double custom10 = atof(argv[10]);
    Polynomial pot(mass, lambda); // (m)
    //Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 14.5); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //error for background variables and PPS
    double err = 1e-5;
    double adaptive_k_grid_err = 5e-4;
    //double err = atof(argv[1]);
    //double adaptive_k_grid_err = atof(argv[2]);
    //minimal log (natural) bin size for k-binning
    //4.6e-3=delta_ln=delta_log10/log10(e)=2e-3/log10(e)
    //bin width: delta_log10 = log10(k2) - log10(k1) for two consecutive k1, k2
    //4.6e-3 corresponds to roughly 500 bins per decade (1/500 = 2e-3)
    double dlogk = 4.6e-3;
    
    //Background Initial Conditions
    //double N_star = 55, N_dagger = 7;
    double N_dagger = atof(argv[2]);
    double N_star = atof(argv[1]);
//    double logk_test = atof(argv[4]);

    //Solve Background Variables
    auto sols = solve_equations(err, potential_ptr, N_star, N_dagger);
    
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
//    NumericModeSolver N_ms(potential_ptr, N_star);
    //Vacuum Setting Time (no. e-folds before end of inflation)
    //sols.N_end is the total no. of e-folds from beginning of background varibale evolution till end of inflation.
    //double N_r = sols.N_end - 2;
    //double N_r = sols.N_start;
    double N_r = N_star + N_dagger;  // N_end-N_r is the starting point of the mode evolution, i.e. for N_r=N_tot at the start of inflation.
    //Set Vacuum Initial Conditions
    //BD = Bunch-Davies. (We did not add any other initial conditions yet, but these could be easily added in the ModeSolver::Initial_Conditions
    ms.Initial_Conditions(BD, N_r);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;

    //////////////////////////////////////////////////////////////////////////////
    //k range to solve (logarithmic scale)
//    double k0 = 1e-6, k1 = 1;
//    std::vector<double> kplot(500);
//    for(size_t n = 0; n < kplot.size(); n++)
//        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
//
//    for(auto k : kplot)
//    {
//        std::cout<< k << "  " << ms.Find_PPS_Scalar(k) << std::endl;
//    }


//    std::cout << "k  = " << exp(logk_test) << std::endl;
//    std::cout << ms.Find_PPS_Scalar(exp(logk_test)) << std::endl;
//    std::cout << ms.Find_PPS_Tensor(exp(logk_test)) << std::endl;


//    std::cout << "k  = " << 1e-6 << std::endl;
    double kmin = 1e-6, kmax = 1e0;
    std::map<double,double> PPS_Scalar;
    std::map<double,double> PPS_Tensor;
    PPS_Scalar[log(kmin)] = log(ms.Find_PPS_Scalar(kmin));
    PPS_Scalar[log(kmax)] = log(ms.Find_PPS_Scalar(kmax));
    PPS_Tensor[log(kmin)] = log(ms.Find_PPS_Tensor(kmin));
    PPS_Tensor[log(kmax)] = log(ms.Find_PPS_Tensor(kmax));
    auto iter_s = PPS_Scalar.begin();
    auto iter_t = PPS_Tensor.begin();
    while (iter_s != std::prev(PPS_Scalar.end()) and iter_t != std::prev(PPS_Tensor.end()))
    {
        auto x0 = iter_s->first;
        auto y0 = iter_s->second;
        auto z0 = iter_t->second;
        auto iter_s_ = std::next(iter_s);
        auto iter_t_ = std::next(iter_t);
        auto x3 = iter_s_->first;
        auto y3 = iter_s_->second;
        auto z3 = iter_t_->second;

        auto x1 = 2.0*x0/3.0 + x3/3.0;
        auto x2 = x0/3.0 + 2.0*x3/3.0;
//        std::cout << "    logk0 " << x0 << "    k0 " << exp(x0) << std::endl;
//        std::cout << "    logk1 " << x1 << "    k1 " << exp(x1) << std::endl;
        auto y1 = log(ms.Find_PPS_Scalar(exp(x1))); 
//        std::cout << "    logk2 " << x2 << "    k2 " << exp(x2) << std::endl;
        auto y2 = log(ms.Find_PPS_Scalar(exp(x2))); 
//        std::cout << "    Tensor1" << std::endl;
        auto z1 = log(ms.Find_PPS_Tensor(exp(x1)));
//        std::cout << "    Tensor2" << std::endl;
        auto z2 = log(ms.Find_PPS_Tensor(exp(x2)));

        auto y1_ = 2.0*y0/3.0 + y3/3.0;
        auto y2_ = y0/3.0 + 2.0*y3/3.0;
        auto z1_ = 2.0*z0/3.0 + z3/3.0;
        auto z2_ = z0/3.0 + 2.0*z3/3.0;

        PPS_Scalar[x1]=y1;
        PPS_Scalar[x2]=y2;
        PPS_Tensor[x1]=z1;
        PPS_Tensor[x2]=z2;

        if ( std::abs(x2-x1) < dlogk or
                ((std::abs(y1-y1_) < adaptive_k_grid_err and std::abs(y2-y2_) < adaptive_k_grid_err) and
                (std::abs(z1-z1_) < adaptive_k_grid_err and std::abs(z2-z2_) < adaptive_k_grid_err)))
        {
            iter_s = iter_s_;
            iter_t = iter_t_;
            std::cout << exp(x0) << " " << exp(y0) << " " << exp(z0) << std::endl;
            int num;
            double xi, yi, zi;
            num = floor((x1 - x0) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x0 + i * dlogk;
                yi = ((x1 - xi) * y0 + (xi - x0) * y1) / (x1 - x0);
                zi = ((x1 - xi) * z0 + (xi - x0) * z1) / (x1 - x0);
                std::cout << exp(xi) << " " << exp(yi) << " " << exp(zi) << std::endl;
            }
            std::cout << exp(x1) << " " << exp(y1) << " " << exp(z1) << std::endl;
            num = floor((x2 - x1) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x1 + i * dlogk;
                yi = ((x2 - xi) * y1 + (xi - x1) * y2) / (x2 - x1);
                zi = ((x2 - xi) * z1 + (xi - x1) * z2) / (x2 - x1);
                std::cout << exp(xi) << " " << exp(yi) << " " << exp(zi) << std::endl;
            }
            std::cout << exp(x2) << " " << exp(y2) << " " << exp(z2) << std::endl;
            num = floor((x3 - x2) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x2 + i * dlogk;
                yi = ((x3 - xi) * y2 + (xi - x2) * y3) / (x3 - x2);
                zi = ((x3 - xi) * z2 + (xi - x2) * z3) / (x3 - x2);
                std::cout << exp(xi) << " " << exp(yi) << " " << exp(zi) << std::endl;
            }
        }
        if (iter_s == std::prev(PPS_Scalar.end()) and iter_t == std::prev(PPS_Tensor.end()))
        {
            std::cout << exp(x3) << " " << exp(y3) << " " << exp(z3) << std::endl;
        }
    }

    return 0;
}
