#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"

int main()
{
    std::cout<<"Solving for Background..."<<std::endl;
    
    //Set Potential and ptr
    Poly_Step pot(6.48757e-6, 0, 5e-3, 15.5);
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    
    //Solve Background Variables
    auto sols = solve_equations(potential_ptr, N_star, N_dagger);
   
    //////////////////////////////////////////////////////////////////////////////
    //k range
    double k0 = 1e-6, k1 = 0.2;
    //Vacuum Setting Time (no. e-folds before end of inflation)
    double N_r = sols.N_end - 2;
    
    std::cout<<"Finding PPS Approximately..."<<std::endl;
    //Initialize ModeSolver with background solutions
    ModeSolver ms(sols);
    
    //Set Vacuum Initial Conditions
    ms.Initial_Conditions(BD, N_r);
    
    //Choose error tolerance (By default set to 5e-3)
    ms.PPS_error = 1e-6;
    
    //Construct PPS linear interpolation
    //ms.Construct_PPS_Tensor(k0, k1, 3e-3);
    
    std::cout<<"Finding PPS Numerically..."<<std::endl;
    //Initialize Numeric Solver
    NumericModeSolver N_ms(potential_ptr, N_star);
    
    //Construct PPS linear interpolation
    //N_ms.Construct_PPS_Scalar(k0, k1, 3e-3);

    //////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::ofstream mout{"output/PPS.txt"};
    mout.precision(20);
    
    std::vector<double> kplot(1000);
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    //ms.k_plot_Scalar are k points of linear intepolation
    for(auto k : kplot)
    {
        double True = 0;//N_ms.Find_PPS(k);
        double Approx = ms.Find_PPS_Scalar(k);
        std::cout<<k<<std::endl;
        //Plot
        mout << k <<"  "<< Approx <<"  "<<(Approx - True) / True<<std::endl;
    }
    mout.close();
    
    return 0;
}
