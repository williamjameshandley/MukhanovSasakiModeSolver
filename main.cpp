#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main(int argc, char** argv)
{
    std::cout.precision(18);
    Polynomial pot(2.9e-6); // (m)
    //Poly_Step pot(6.48757e-6, 2e-3, 5e-2, 14.5); // (m, A, Delta, phi_0)
    //Starobinsky pot(1.2e-5); // (m)
    //AxionMonodromy pot(1.4e-5, 4./3, 12.38, 0.01, 0.05, -1./3, 0); //(m, p, phi0, f0, b, p_f, gamma0)
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //error for background variables and PPS
    double err = 1e-5;
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 4.7;

    //Solve Background Variables
    auto sols = solve_equations(err, potential_ptr, N_star, N_dagger);
    //////////////////////////////////////////////////////////////////////////////
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    double N_r = sols.N_end;  // N_end-N_r is the starting point of the  mode evolution, i.e. for N_r=N_tot at the start of   inflation.
    //Set Vacuum Initial Conditions
    ms.Initial_Conditions(discrete, N_r);
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = err;
    //////////////////////////////////////////////////////////////////////////////
    //k range to solve (logarithmic scale)
    double k0 = 1e-4, k1 = 1.1e-1;
    double N = 10000;
    std::vector<double> kplotlog(N);
    for(size_t n = 0; n < kplotlog.size(); n++)
        kplotlog[n] = k0 * exp(static_cast<double>(n) * 1.0 *   (log(k1)   - log(k0)) /   static_cast<double>(kplotlog.size()));
    std::vector<double> kplotlin(N);
    for(size_t n = 0; n < kplotlin.size(); n++)
        kplotlin[n] = k0 + (k1-k0) * (1.0 * n) / kplotlin.size();
    
    std::ofstream file;
    file.open ("output/PPS.txt");
   
    for(auto k : kplotlog)
    {
        double PPS = ms.Find_PPS_Scalar(k);
        file << k << "   " << PPS<<std::endl;
    }

    return 0;
}
