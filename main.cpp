#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/ModeSolver.cpp"
#include "src/Transitions.cpp"
#include "src/Potential.hpp"
#include "src/Special_Functions.hpp"

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
    std::cout<<"Finding Transitions: ";
    //Transitions Conditions
    double eta_end = background_sols.eta.back();
    double eta_i = 0.02 * eta_end;
    double eta_f = 0.95 * eta_end;
    
    //Transitions Constructor
    Transitions T(eta_i, eta_f, background_sols);
    
    //Find Transitions
    double error = 1e-4;
    auto transitions_sols = T.Find(error);
    
    std::cout<<transitions_sols.eta_step.size()<<std::endl;
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS"<<std::endl;
    size_t N = 1000;
    std::vector<double> k(N);
    double k0 = 1.0, k1 = 1.0e6;
    
    for(size_t i = 0; i < N; i++)
        k[i] = (exp(i * (log(k1) - log(k0)) / N)); //logspace
    
    ModeSolver ms(background_sols, transitions_sols);
    
    ms.Initial_Conditions("BD", 0.1 * eta_end);
    
    //Output
    std::ofstream pout;
    pout.open ("output/PPS.txt");
    
    for(size_t i = 0; i < k.size(); i++)
    {
        pout<<k[i]<<"   "<<ms.PPS(k[i])<<std::endl;
    }

    // for(auto k_i : k)
    //     pout<<k_i<<"   "<<ms.PPS(k_i)<<std::endl;
    
    pout.close();
    
    return 0;
}
