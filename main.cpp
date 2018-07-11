#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/Potential.hpp"

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    double t0 = 1e-6, t1 = 20.0, phi_p = 23.0, dphi_p = -sqrt(2.0/3.0) / t0;
    
    Poly pot;
    pot.m = 1.0, pot.lambda = 0;
    
    BackgroundSolver variables(t0, t1, phi_p, dphi_p, pot);
    
    double abs_err = 1.0e-5, rel_err = 1.0e-2;
    
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);
    
    std::vector<double> ddz, eta;
    
    std::tie(ddz, eta) = variables.Solve(integrator);
    
    std::ofstream fout;
    fout.open ("bin/output/ddz.txt");
    
    for(size_t i = 0; i < ddz.size(); i++)
    {
        fout<<eta[i]<<"   "<<ddz[i]<<std::endl;
    }
    
    std::cout<<ddz.size()<<std::endl;
    
    fout.close();
    
    return 0;
}
