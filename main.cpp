#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include "src/BackgroundSolver.hpp"
#include "src/ModeSolver.hpp"
#include "src/Potential.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    std::cout<<"Solving for Background..."<<std::endl;
    //Background Solver Constructor
    BackgroundSolver solver;
    
    //Set Integrator
    double abs_err = 1.0e-5, rel_err = 1e-4;
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);
    
    //Background Initial Conditions
    double t0 = 6.48757e-6, t1 = 20.0, phi_p = 23.08546, dphi_p = -sqrt(2.0/3.0) / t0;
    
    //Set Potential
    Polynomial pot(1);
    
    //Solve Background Variables
    auto background_sols = solver.Solve(integrator, pot, t0, t1, phi_p, dphi_p);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS..."<<std::endl;
    
    double k0 = 1, k1 = 1e6;
    double eta_r = 0.1 * background_sols.eta.back();
    
    ModeSolver ms(background_sols, 0.02);
    ms.Initial_Conditions(BD, eta_r);
    ms.Construct_PPS(k0, k1);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    std::vector<double> kplot(10000);
    
    for(size_t n = 0; n < kplot.size(); n++)
        kplot[n] = k0 * exp(static_cast<double>(n) * 1.0 * (log(k1) - log(k0)) / static_cast<double>(kplot.size()));
    
    std::ofstream mout{"output/PPS.txt"};
   // for(auto k : kplot) mout << k << " " << ms.PPS(k) << std::endl;
  //  mout.close();
    
    FILE *test;
    test = fopen("output/test.txt", "r");
    
    double x = 0;
    for(int n = 0; n < 10000; n++)
    {
        fscanf(test,"%lf",&x);
        double k = x;// - 0.01 / (2 * x);
        fscanf(test,"%lf",&x);
        double PPS_true = x;
        
        mout<<k<<"   "<<(PPS_true - ms.PPS(k)) / PPS_true<<std::endl;
        if(n == 0)
        {
            std::cout<<(PPS_true - ms.PPS(k)) / PPS_true<<std::endl;
        }
    }
    return 0;
}
