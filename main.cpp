#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/Potential.hpp"


int main()
{
    double t0 = 1e-6;
    double t1 = 20.0;
    double phi_p = 23.0;
    double dphi_p = -sqrt(2.0/3.0) / t0;
    
    Poly pot;
    pot.m = 1.0;
    pot.lambda = 0;
    
    BackgroundSolver BackgroundVar(t0, t1, phi_p, dphi_p, pot);
    
    std::vector<double> ddz;
    std::vector<double> eta;
    
    boost::numeric::odeint::runge_kutta4<std::vector<double>> integrator;
    
    std::tie(ddz, eta) = BackgroundVar.Solve(integrator);
    
    std::ofstream fout;
    fout.open ("bin/output/ddz.txt");
    
    for(std::vector<double>::size_type i = 0; i < ddz.size(); i++)
    {
        fout<<eta[i]<<"   "<<ddz[i]<<std::endl;
    }
    
    fout.close();
    
    return 0;
}
