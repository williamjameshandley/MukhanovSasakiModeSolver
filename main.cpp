#include <iostream>
#include <fstream>
#include <ctime>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    //Choose Potential and set potential_ptr which is passed to background solver
    Polynomial pot(6e-6); // (m)
    //Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 14.5); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    
    //error for background variables and PPS
    double err = 1e-5;
    //Solve Background Variables
    auto phi_p = 23.5;
    auto t0 = 1.;
    auto ti = 1e4;
    auto sols = solve_equations(err*1e-1, potential_ptr, t0, ti, phi_p);

    auto N_star = 55;
    sols.aH_star = sols.aH(sols.N_end - N_star);
    
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    //Vacuum Setting Time (no. e-folds before end of inflation)
    //sols.N_end is the total no. of e-folds from beginning of background varibale evolution till end of inflation.
    //Set Vacuum Initial Conditions
    //BD = Bunch-Davies. (We did not add any other initial conditions yet, but these could be easily added in the ModeSolver::Initial_Conditions
    ms.Initial_Conditions(BD, 0);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;

    //////////////////////////////////////////////////////////////////////////////
    //k range to solve (logarithmic scale)
    double k0 = 1e-6, k1 = 1;
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));

    for(auto k : kplot)
    {
        std::cout<< k << "  " << ms.Find_PPS_Scalar(k) << "  " << ms.Find_PPS_Tensor(k) << std::endl;
    }

    return 0;
}
