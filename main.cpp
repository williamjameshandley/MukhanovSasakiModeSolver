#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.hpp"
#include "src/Potential.hpp"
#include "src/ModeSolver.hpp"
#include <map>

struct Vars
{
    double a, b;
    int i;
    Vars(double _a, double _b, int _i): a{_a}, b{_b}, i{_i} {};
};

int main()
{
    std::map<double,Vars> A;

    A.insert(std::pair<double,Vars> (5.5, Vars(1, 2, 3)));
    A.insert(std::pair<double,Vars> (6.5, Vars(3, 2, 1)));
    A.insert(std::pair<double,Vars> (7.5, Vars(2, 1, 3)));

    A.at(5.5) = Vars(7, 6, 5);

    for (auto iter = A.begin(); iter != A.end(); ++iter)
        std::cout << iter->second.a <<"  "<< iter->second.b<<"  " << iter->second.i<< std::endl;

    return 0;
    
    std::cout<<"Solving for Background..."<<std::endl;
    
    //Set Potential and ptr
    Poly_Step pot(6.48757e-6, 1e-3, 5e-3, 15.5);
    auto potential_ptr = static_cast<Potential*> (&pot);
    
    //Background Initial Conditions
    double N_star = 55, N_dagger = 7;
    
    //Solve Background Variables
    auto sols = solve_equations(potential_ptr, N_star);
   
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
    ms.PPS_error = 5e-3;
    
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
        double True = N_ms.Find_PPS(k);
        double Approx = ms.Find_PPS_Scalar(k);
        
        //Plot
        //mout << k <<"  "<< Approx <<"  "<<(Approx - True) / True<<std::endl;
    }
    mout.close();
    
    return 0;
}
