#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/ModeSolver.cpp"
#include "src/Transitions.cpp"
#include "src/Potential.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    double t0 = 6.48757e-6, t1 = 20.0, phi_p = 23.08546, dphi_p = -sqrt(2.0/3.0) / t0; //Background Initial Conditions
    
    Poly pot;                                   //Set Potential
    pot.m = 1.0, pot.lambda = 0;
    
    BackgroundSolver variables(t0, t1, phi_p, dphi_p, pot);    //Background Solver Constructor
    
    double abs_err = 1.0e-5, rel_err = 1e-3;
    
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);   //Set Integrator
    
    std::vector<double> ddz, eta;
    
    std::tie(ddz, eta) = variables.Solve(integrator);      //Solve Background Variables
    
    std::cout<<ddz.size()<<std::endl;
    
    //Transitions
    double eta_end = eta.back();
    double eta_i = 0.02 * eta_end;
    double eta_f = 0.95 * eta_end;
    std::vector<double> a, b, eta_step;
    double delta;
    
    Transitions T(eta_i, eta_f, eta_end, ddz, eta);     //Transitions Constructor
    
    double error = 1e-4;

    std::tie(a, b, delta, eta_step) = T.Find(error);        //Find Transitions
    
    std::cout<<eta_step.size()<<std::endl;
    
    std::vector<double> k;
    double k0 = 1.0, k1 = 1.0e4;
    size_t N = 100;
    
    for(size_t i = 0; i < N; i++)
    {
        k.push_back(exp(i * (log(k1) - log(k0)) / N)); //logspace
    }
    
    ModeSolver ms(k, eta_step, a, b, delta, eta_end, ddz, ddz, ddz, eta);
    
    ms.Find_Mat();
    
    ms.Initial_Conditions("BD", 0.2 * eta_end);
    
    std::vector<double> PPS;
    
    PPS = ms.PPS();
    
    std::ofstream fout;
    fout.open ("bin/output/PPS.txt");
    
    for(size_t i = 0; i < k.size(); i++)
    {
        fout<<k[i]<<"   "<<PPS[i]<<std::endl;
    }
    
    fout.close();
    
    return 0;
}
