#include <iostream>
#include <fstream>
#include <odepack/dlsodar.hpp>
#include <cmath>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include "src/BackgroundSolver.hpp"
#include "src/ModeSolver.hpp"
#include "src/Potential.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

void simple_pendulum_field(double qdot[], const double t, const double y[], void* params)
{
    qdot[0] = y[1];
    qdot[1] = - sin(y[0]);
}

int main()
{
    // Two variables, one 'event'
    dlsodar desolver(2, 1);

    // Initialise at a angle very close to pi
    double y[2] = {M_PI * 99999/100000.0, 0.};
    double t = 0;
    double alpha = 1.0;
    double theta0 = 0.0;
    double params[2] = {alpha, theta0};

    // now output a few oscillations to file
    std::ofstream fout{"simple_pendulum_trajectory.txt"};
    while(y[0] > theta0){
        desolver.integrate(t, t + 0.1, y, simple_pendulum_field, params);
        fout << t << " " << y[0] << " " << y[1] << std::endl;
    }

    return 0; 

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
    for(auto k : kplot) mout << k << " " << ms.PPS(k) << std::endl;
    mout.close();
    
    return 0;
}


