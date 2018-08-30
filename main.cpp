#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    
    std::cout<<"Solving for Background..."<<std::endl;
    
    //Set Potential
    Poly_Step pot(6.48757e-6, 0, 3e-2, 14);
    //Set Potential ptr
    auto potential_ptr = static_cast<Potential*> (&pot);
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    //Solve Background Variables
    auto sols = solve_equations(potential_ptr, N_star, N_dagger);
    
    //////////////////////////////////////////////////////////////////////////////
    //k range
    double k0 = 1e-6, k1 = 1;
    //Set Vacuum (no. e-folds before end of inflation)
    double N_r = sols.N_end - 3;
    
    std::cout<<"Finding PPS Approximately..."<<std::endl;
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    //Set Initial Conditions
    ms.Initial_Conditions(BD, N_r);
    //Construct PPS linear interpolation
    ms.Construct_PPS_Scalar(k0, k1, 5e-3);
    
    
    std::cout<<"Finding PPS Numerically..."<<std::endl;
    //Initialize Numeric Solver
    NumericModeSolver N_ms(potential_ptr, N_star, N_dagger, N_r);
    //Construct PPS linear interpolation
    N_ms.Construct_PPS_Scalar(k0, k1, 5e-3);
    
    //////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::ofstream mout{"output/PPS.txt"};
    
    //ms.k_plot_Scalar are k points of linear intepolation
    for(auto k : ms.k_plot_Scalar)
    {
        //Call PPS linear interpolation for plotting
        mout << k << " " << N_ms.PPS_Scalar(k)<<"  "<<ms.PPS_Scalar(k)<<std::endl;
    }
    mout.close();
    
    return 0;
    
}
