#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <utility>
#include "src/BackgroundSolver.hpp"
#include "src/ModeSolver.hpp"
#include "src/Transitions.hpp"
#include "src/Potential.hpp"
#include "src/Special_Functions.hpp"
#include "src/linear_interpolation.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    std::cout<<"Solving for Background..."<<std::endl;
    //Background Solver Constructor
    BackgroundSolver solver;
    
    //Set Integrator
    double abs_err = 1.0e-5, rel_err = 1e-3;
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);
    
    //Background Initial Conditions
    double t0 = 0.0000124037, t1 = 72.1795, phi_p = 14.86452, dphi_p = -sqrt(2.0/3.0) / t0;
    
    //Set Potential
    Starobinsky pot(1.9119);
    
    //Solve Background Variables
    auto background_sols = solver.Solve(integrator, pot, t0, t1, phi_p, dphi_p);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    
    double k0 = 1, k1 = 1e6;
    double eta_r = 0.2 * background_sols.eta.back();
    
    ModeSolver ms(background_sols);
    ms.Initial_Conditions(RST, eta_r);
    ms.Construct_PPS(k0, k1);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(10000);
    
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    std::ofstream mout{"output/PPS.txt"};
    for(auto k : kplot) mout << k << " " << ms.PPS(k) << std::endl;
    mout.close();
    
    return 0;
}
