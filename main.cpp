#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    
    std::cout<<"Solving for Background..."<<std::endl;
    
    Poly_Step pot(6.48757e-6, 1e-4, 5e-3, 15.5);                       //Set Potential
    auto potential_ptr = static_cast<Potential*> (&pot);            //Set Potential ptr
    double N_star = 55, N_dagger = 3;                               //Background Initial Conditions
    auto sols = solve_equations(potential_ptr, N_star, N_dagger);   //Solve Background Variables

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    
    double k0 = 1e-6, k1 = 1;                                       //k range
    double N_r = sols.N_end - 2;                                    //Set Vacuum (no. e-folds before end of inflation)
    //ModeSolver ms(sols);                                            //Initialize ModeSolver with background solutions
    //ms.Initial_Conditions(BD, N_r);                                 //Set Initial Conditions
    //ms.Construct_PPS(k0, k1, 1e-2);                                 //Construct PPS linear interpolation
    
    
    NumericModeSolver Nms(potential_ptr, N_star, 1);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(100);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot)                                         //ms.k_plot are k points of linear intepolation
    {
        std::cout<<k<<std::endl;
        mout << k << " " << Nms.Find_PPS(k)<< std::endl;                  //Call PPS linear interpolation for plotting
    }
    mout.close();
    
    return 0;
    
}
