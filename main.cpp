#include <iostream>
#include <fstream>
#include <ctime>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    //Choose Potential and set potential_ptr which is passed to background solver
    //Polynomial pot(6e-6); // (m)
    Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 15.3); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    
    //error for background variables and PPS
    double err = 1e-5;
    //Solve Background Variables
    auto sols = solve_equations(err*1e-1, potential_ptr, N_star, N_dagger);
    
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    //Vacuum Setting Time (no. e-folds before end of inflation)
    //sols.N_end is the total no. of e-folds from beginning of background varibale evolution till end of inflation.
    double N_r = N_star + N_dagger + 1 ;
    if(N_r > sols.N_end)
        throw std::runtime_error("Vacuum Setting Time Invalid");
    //Set Vacuum Initial Conditions
    //BD = Bunch-Davies.
    //AV = Alpha Vacua. (3rd arg is alpha, leave blank for other initial conditions)
    ms.Initial_Conditions(AV, N_r, 0.3);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;

    //////////////////////////////////////////////////////////////////////////////
    double k0 = 1e-6, k1 = 1;
    //k range to solve (logarithmic scale)
    std::vector<double> kplot_log(1000);
    for(size_t n = 0; n < kplot_log.size(); n++)
        kplot_log[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot_log.size()));
    
    //k range to solve (logarithmic scale)
    std::vector<double> kplot_lin(1000);
    for(size_t n = 0; n < kplot_lin.size(); n++)
        kplot_lin[n] = static_cast<double>(n) * (k1 - k0) / static_cast<double>(kplot_lin.size()) + k0;

    std::ofstream fout("output/PPS.txt");
    
    for(auto k : kplot_log)
    {
        fout<< k << "  " << ms.Find_PPS_Scalar(k)<<std::endl;
        std::cout<<k<<std::endl;
    }

    return 0;
}
