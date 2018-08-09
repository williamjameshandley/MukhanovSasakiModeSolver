#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <odepack/dlsodar.hpp>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/Transitions.hpp"

int main()
{

    std::cout<<"Solving for Background..."<<std::endl;
    //Set Potential
    Poly_Step pot(6.48757e-6, 0.01, 0.05, 14);
    auto pot_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double phi_p = 23, dphi_p = -sqrt(2./3);//17, dphi_p = - pot.dV(phi_p) / (3 * sqrt(pot.V(phi_p) / 3));
    
    //Solve Background Variables
    auto sols = solve_equations(pot_ptr, phi_p, dphi_p);
    
    std::cout<<"Solving for Transitions..."<<std::endl;
    double k = 1e-3;
    double N_i = 2, N_f = 60;
    Transitions T(N_i, N_f, sols);
    auto Tsol = T.Find(k, 1e-1);
    
    std::ofstream fout{"output/trans.txt"};
    fout.precision(20);
    
    std::ofstream mout{"output/var.txt"};
    mout.precision(20);
    
    for(size_t n = 0; n < sols.N.size(); n++)
    {
        mout<<sols.N[n]<<"   "<<(sols.omega_2[n] + k*k / pow(exp(sols.N[n]) * sols.H[n], 2))<<"   "<<sols.d_omega_2[n]<<std::endl;
    }
    mout.close();
    
    for(size_t n = 0; n < Tsol.N_step.size()-1; n++)
    {
        fout<<Tsol.N_step[n]<<"   "<<Tsol.a[n] + Tsol.b[n] * Tsol.N_step[n]<<std::endl;
    }
    fout.close();
    
    std::cout<<Tsol.N_step.size()<<std::endl;
    
    return 0;
    
    /*
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    
    double k0 = 1e-6, k1 = 1;
    double eta_r = 0.1 * sols.eta.back();
    
    ModeSolver ms(sols, 0.002);
    ms.Initial_Conditions(BD, eta_r);
    ms.Construct_PPS(k0, k1);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(10000);
    
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot) mout << k << " " << ms.PPS(k) << std::endl;
    mout.close();
*/
    return 0;
}


