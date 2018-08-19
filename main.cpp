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
    Polynomial pot(6.48757e-6);//, 0.001, 0.005, 15.5);
    auto pot_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 3;
    
    //Solve Background Variables
    auto sols = solve_equations(pot_ptr, N_star, N_dagger);
    std::cout<<sols.N.size()<<std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    double k0 = 1e-6, k1 = 1;
    double N_r = 2;
    
    ModeSolver ms(sols);
    ms.Initial_Conditions(BD, N_r);
    //ms.Construct_PPS(k0, k1, 0.5e-2);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    FILE *PPS;
    PPS = fopen("output/PPS_.txt", "r");
    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot)
    {
        double x, y;
        fscanf(PPS,"%lf",&x);
        fscanf(PPS,"%lf",&y);
        
        mout << x << " " << (ms.Find_PPS(x) - y) / y<< std::endl;
    }
    mout.close();

    return 0;
}


