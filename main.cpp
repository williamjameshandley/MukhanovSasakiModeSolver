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
    
    BackgroundSolver Solver;
    Solver.t0 = t0;
    Solver.t1 = t1;
    Solver.phi_p = phi_p;
    Solver.dphi_p = dphi_p;
    Solver.m = 1.0;
    Solver.lambda = 0;
    
    std::vector<double> ddz;
    std::vector<double> eta;
    asc::RK4 integrator;
    
    std::tie(ddz, eta) = Solver.Solve(integrator);
    
    std::ofstream fout;
    fout.open ("bin/output/ddz.txt");
    
    for(std::vector<double>::size_type i = 0; i < ddz.size(); i++)
    {
        fout<<eta[i]<<"   "<<ddz[i]<<std::endl;
    }
    
    fout.close();
    
    return 0;
}
