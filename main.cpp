#include <iostream>
#include <fstream>
#include <ctime>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main(int argc, char** argv)
{
    if (argc != 3) return -1;
    std::cout.precision(18);
    //Choose Potential and set potential_ptr which is passed to background solver
    Polynomial pot(6e-6); // (m)
    //Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 14.5); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    
    //error for background variables and PPS
    //double err = 1e-5;
    double err = atof(argv[1]);
    double PPS_err = atof(argv[2]);
    double dlogk = 1e-3;

    //Solve Background Variables
    auto sols = solve_equations(err, potential_ptr, N_star, N_dagger);
    
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    //Vacuum Setting Time (no. e-folds before end of inflation)
    //sols.N_end is the total no. of e-folds from beginning of background varibale evolution till end of inflation.
    double N_r = sols.N_end - 2;
    //Set Vacuum Initial Conditions
    //BD = Bunch-Davies. (We did not add any other initial conditions yet, but these could be easily added in the ModeSolver::Initial_Conditions
    ms.Initial_Conditions(BD, N_r);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;

    //////////////////////////////////////////////////////////////////////////////
    //k range to solve (logarithmic scale)
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
                (std::abs(y1-y1_) < PPS_err and std::abs(y2-y2_) < PPS_err
           ))
        {
            iter = iter_;
            std::cout << exp(x0) << " " << exp(y0) << std::endl;
            std::cout << exp(x1) << " " << exp(y1) << std::endl;
            std::cout << exp(x2) << " " << exp(y2) << std::endl;
        }
        if (iter == std::prev(PPS_Scalar.end()))
        {
            std::cout << exp(x3) << " " << exp(y3) << std::endl;
        }
    }

    return 0;
}
