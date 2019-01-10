#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main(int argc, char** argv)
{
    if (argc != 4) return -1;
    std::cout.precision(18);
    //Choose Potential and set potential_ptr which is passed to background solver
    double mass = exp(atof(argv[1])) * 1e-6;
    Polynomial pot(mass); // (m)
    //Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 14.5); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //error for background variables and PPS
    double err = 5e-4;
    double PPS_err = 5e-4;
    //double err = atof(argv[1]);
    //double PPS_err = atof(argv[2]);
    //minimal log (natural) bin size for k-binning
    //4.6e-3=delta_ln=delta_log10/log10(e)=2e-3/log10(e)
    //4.6e-3 corresponds to roughly 500 bins per decade
    double dlogk = 4.6e-3;
    
    //Background Initial Conditions
    //double N_star = 55, N_dagger = 7;
    double N_dagger = atof(argv[2]);
    double N_star = atof(argv[3]);

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
    double N_r = N_star + N_dagger;
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

    double kmin = 1e-6, kmax = 1e0;
    std::map<double,double> PPS_Scalar;
    PPS_Scalar[log(kmin)] = log(ms.Find_PPS_Scalar(kmin));
    PPS_Scalar[log(kmax)] = log(ms.Find_PPS_Scalar(kmax));
    auto iter = PPS_Scalar.begin();
    while (iter != std::prev(PPS_Scalar.end()))
    {
        auto x0 = iter->first;
        auto y0 = iter->second;
        auto iter_ = std::next(iter);
        auto x3 = iter_->first;
        auto y3 = iter_->second;

        auto x1 = 2.0*x0/3.0 + x3/3.0;
        auto x2 = x0/3.0 + 2.0*x3/3.0;
        auto y1 = log(ms.Find_PPS_Scalar(exp(x1))); 
        auto y2 = log(ms.Find_PPS_Scalar(exp(x2))); 

        auto y1_ = 2.0*y0/3.0 + y3/3.0;
        auto y2_ = y0/3.0 + 2.0*y3/3.0;

        PPS_Scalar[x1]=y1;
        PPS_Scalar[x2]=y2;

        if ( std::abs(x2-x1) < dlogk or
                (std::abs(y1-y1_) < PPS_err and std::abs(y2-y2_) < PPS_err))
        {
            iter = iter_;
            std::cout << exp(x0) << " " << exp(y0) << std::endl;
            int num;
            double xi, yi;
            num = floor((x1 - x0) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x0 + i * dlogk;
                yi = ((x1 - xi) * y0 + (xi - x0) * y1) / (x1 - x0);
                std::cout << exp(xi) << " " << exp(yi) << std::endl;
            }
            std::cout << exp(x1) << " " << exp(y1) << std::endl;
            num = floor((x2 - x1) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x1 + i * dlogk;
                yi = ((x2 - xi) * y1 + (xi - x1) * y2) / (x2 - x1);
                std::cout << exp(xi) << " " << exp(yi) << std::endl;
            }
            std::cout << exp(x2) << " " << exp(y2) << std::endl;
            num = floor((x3 - x2) / dlogk);
            for ( int i = 1; i <= num; i += 1)
            {
                xi = x2 + i * dlogk;
                yi = ((x3 - xi) * y2 + (xi - x2) * y3) / (x3 - x2);
                std::cout << exp(xi) << " " << exp(yi) << std::endl;
            }
        }
        if (iter == std::prev(PPS_Scalar.end()))
        {
            std::cout << exp(x3) << " " << exp(y3) << std::endl;
        }
    }

    return 0;
}
