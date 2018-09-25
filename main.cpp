#include <iostream>
#include <fstream>
#include <ctime>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    //Set Potential and ptr
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    Poly_Step pot(6.48757e-6, 0, 5e-2, 14.5);
    //Starobinsky pot(1.2e-5);
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    
    double err = 1e-5;
    //Solve Background Variables
    auto sols = solve_equations(err*1e-1, potential_ptr, N_star, N_dagger);
    
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    //Vacuum Setting Time (no. e-folds before end of inflation)
    double N_r = sols.N_end - 2;
    //Set Vacuum Initial Conditions
    ms.Initial_Conditions(BD, N_r);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;

    //////////////////////////////////////////////////////////////////////////////
    //k range
    double k0 = 1e-6, k1 = 1;
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));

    std::ofstream mout{"output/PPS.txt"};
    mout.precision(20);
    for(auto k : kplot)
    {
        //Plot
        mout << k << "  " << ms.Find_PPS_Scalar(k) << "  " << ms.Find_PPS_Tensor(k) << std::endl;
    }
    mout.close();

    return 0;
}
