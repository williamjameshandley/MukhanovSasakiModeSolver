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
    Poly_Step pot(6.48757e-6, 1e-4, 5e-3, 15.5);
    auto pot_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 3;
    
    //Solve Background Variables
    auto sols = solve_equations(pot_ptr, N_star);
    
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
        mout << k << " " << ms.Find_PPS(k)<< std::endl;
    }
    mout.close();
    
    return 0;
    
}
