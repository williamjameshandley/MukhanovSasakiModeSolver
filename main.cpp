#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <odepack/dlsodar.hpp>
#include "src/BackgroundSolver.hpp"
#include "src/ModeSolver.hpp"
#include "src/Potential.hpp"


int main()
{

    std::cout<<"Solving for Background..."<<std::endl;
    //Background Initial Conditions
    double phi_p = 23.08546, dphi_p = -10 * sqrt(2.0/3.0);
    
    //Set Potential
    Polynomial pot(6.48757e-6);
    auto pot_ptr = static_cast<Potential*> (&pot);
    
    //Solve Background Variables
    auto background_sols = solve_equations(pot_ptr, phi_p, dphi_p);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    
    double k0 = 1e-6, k1 = 1;
    double eta_r = 0.1 * background_sols.eta.back();
    
    ModeSolver ms(background_sols, 0);
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
    /*
    FILE *test;
    test = fopen("output/test.txt", "r");
    
    double x = 0;
    for(int n = 0; n < 10000; n++)
    {
        fscanf(test,"%lf",&x);
        double k = x;
        fscanf(test,"%lf",&x);
        double PPS_true = x;
        
        mout<<k<<"   "<<(PPS_true - ms.PPS(k)) / PPS_true<<std::endl;
        if(n == 0)
        {
            std::cout<<(PPS_true - ms.PPS(k)) / PPS_true<<std::endl;
        }
    }
    */
    return 0;
}


