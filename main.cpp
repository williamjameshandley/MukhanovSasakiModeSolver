#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/ModeSolver.cpp"
#include "src/Transitions.cpp"
#include "src/Potential.hpp"
#include "src/Special_Functions.hpp"
#include "src/linear_interpolation.hpp"
#include <utility>

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    std::cout<<"Solving for Background"<<std::endl;
    //Background Solver Constructor
    BackgroundSolver solver;
    
    //Background Initial Conditions
    double t0 = 6.48757e-6, t1 = 20.0, phi_p = 23.08546, dphi_p = -sqrt(2.0/3.0) / t0;
    
    //Set Potential
    Polynomial pot(1.0);
    
    //Set Integrator
    double abs_err = 1.0e-5, rel_err = 1e-3;
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);
    
    //Solve Background Variables
    auto background_sols = solver.Solve(integrator, pot, t0, t1, phi_p, dphi_p);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS: ";
    double k0 = 1, k1 = 1e6;
    double eta_r = 0.1 * background_sols.eta.back();
    ModeSolver ms(background_sols);
    ms.Initial_Conditions(BD, eta_r);
    ms.Construct_PPS(k0, k1);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(10000);
    
    for(size_t n = 0; n < kplot.size(); n++)
    {
        kplot[n] = exp(n * 1.0 * (log(k1) - log(k0)) / kplot.size());
    }
    
    std::ofstream mout;
    mout.open("output/PPS.txt");
    for(size_t n = 0; n < kplot.size(); n++)
    {
        mout<<kplot[n]<<"   "<<ms.PPS(kplot[n])<<std::endl;
    }
    mout.close();
    
    
    return 0;
}
