#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/Transitions.cpp"
#include "src/Potential.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    double t0 = 1e-6, t1 = 20.0, phi_p = 23.0, dphi_p = -sqrt(2.0/3.0) / t0; //Background Initial Conditions
    
    Poly pot;                                   //Set Potential
    pot.m = 1.0, pot.lambda = 0;
    
    BackgroundSolver variables(t0, t1, phi_p, dphi_p, pot);    //Background Solver Constructor
    
    double abs_err = 1.0e-5, rel_err = 1.0e-3;
    
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);   //Set Integrator
    
    std::vector<double> ddz, eta;
    
    std::tie(ddz, eta) = variables.Solve(integrator);      //Solve Background Variables
    
    //Transitions
    double eta_end = eta.back();
    double eta_i = 0.02 * eta_end;
    double eta_f = 0.95 * eta_end;
    std::vector<double> a, b, eta_step;
    double delta;
    
    Transitions T(eta_i, eta_f, eta_end, ddz, eta);     //Transitions Constructor
    
    std::tie(a, b, delta, eta_step) = T.Find();        //Find Transitions
    
    
    
    
    
    
    
    
    std::ofstream fout;
    fout.open ("bin/output/ddz.txt");
    
    for(size_t i = 0; i < ddz.size(); i++)
    {
        fout<<eta[i] / eta.back()<<"   "<<ddz[i]<<std::endl;
    }
    
    fout.close();
    
    std::ofstream mout;
    mout.open ("bin/output/segments.txt");
    
    for(size_t i = 0; i < eta_step.size(); i++)
    {
        mout<<eta_step[i] / eta.back()<<"   "<<ddz[T.search(eta, eta_step[i])]<<std::endl;
    }
    
    mout.close();
    
    return 0;
}
