#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <odepack/dlsodar.hpp>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    
    std::cout<<"Solving for Background..."<<std::endl;
    //Set Potential
    Poly_Step pot(6.48757e-6, 5e-4, 5e-3, 15.5);
    auto pot_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 3;
    
    //Solve Background Variables
    auto sols = solve_equations(pot_ptr, N_star);
    /*
    std::vector<double> N(100000);
    for(size_t n = 0; n < N.size(); n++)
        N[n] = 1.0 * n * sols.N_end / N.size();
    
    double k = 890;
    std::ofstream dout{"output/var.txt"};
    for(size_t n = 0; n < N.size(); n++)
    {
        dout << N[n] << " " << sols.omega_2(N[n]) + k * k * exp(-2 * sols.log_aH(N[n]))<< std::endl;
    }
    dout.close();
    
    Transitions T(0, sols.N_end, sols);
    auto Tsol = T.Find(k, 1e-4);
    
    std::ofstream fout{"output/var1.txt"};
    fout.precision(20);
    for(size_t n = 0; n < Tsol.log_a.size(); n++)
    {
        fout << Tsol.log_N_step[n]<< " " << exp(Tsol.log_a[n] + Tsol.log_b[n] * Tsol.log_N_step[n])<< std::endl;
    }
    fout << Tsol.log_N_step.back()<< " " << exp(Tsol.log_a.back() + Tsol.log_b.back() * Tsol.log_N_step.back())<< std::endl;
    fout.close();
    
    std::ofstream rout{"output/var2.txt"};
    rout.precision(20);
    for(size_t n = 0; n < Tsol.lin_a.size(); n++)
    {
        rout << Tsol.lin_N_step[n]<< " " << (Tsol.lin_a[n] + Tsol.lin_b[n] * Tsol.lin_N_step[n])<< std::endl;
    }
    rout.close();
    return 0;
    */
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    double k0 = 1e-6, k1 = 1;
    
    //Vacuum Setting Time
    double N_r = 1;
    
    ModeSolver ms(sols);
    ms.Initial_Conditions(BD, N_r);
    //ms.Construct_PPS(k0, k1, 1e-3);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));

    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot)
    {
        //std::cout<<k<<std::endl;
        mout << k << " " << ms.Find_PPS(k)<< std::endl;
    }
    mout.close();
    
    return 0;
    
}
