#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    
    std::cout<<"Solving for Background..."<<std::endl;
    
    Poly_Step pot(6.48757e-6, 0, 5e-3, 15.5);                       //Set Potential
    auto potential_ptr = static_cast<Potential*> (&pot);            //Set Potential ptr
    double N_star = 55, N_dagger = 7;                               //Background Initial Conditions
    auto sols = solve_equations(potential_ptr, N_star, N_dagger);   //Solve Background Variables

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double k0 = 1e-6, k1 = 1;                                       //k range
    double N_r = sols.N_end - 2;                                    //Set Vacuum (no. e-folds before end of inflation)
    
    std::cout<<"Finding PPS Approximately..."<<std::endl;
    ModeSolver ms(sols);                                            //Initialize ModeSolver with background solutions
    ms.Initial_Conditions(BD, N_r);                                 //Set Initial Conditions
    //ms.Construct_PPS(k0, k1, 5e-3);                               //Construct PPS linear interpolation
    
    std::cout<<"Finding PPS Numerically..."<<std::endl;
    NumericModeSolver N_ms(potential_ptr, N_star, N_dagger, N_r);   //Initialize Numeric Solver
    //N_ms.Construct_PPS(k0, k1, 5e-3);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot)                                         //ms.k_plot are k points of linear intepolation
    {
        mout << k << " " << N_ms.Find_PPS(k)<<"   "<<ms.Find_PPS(k)<<std::endl;    //Call PPS linear interpolation for plotting
    }
    mout.close();
    
    return 0;
    
}
