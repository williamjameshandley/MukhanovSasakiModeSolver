#include <iostream>
#include "src/BackgroundSolver.cpp"
#include "src/Potential.cpp"


int main()
{
    double t0 = 1e-6;
    double t1 = 20.0;
    double phi_p = 23.0, dphi_p = -sqrt(2.0/3.0) / t0;
    
    Poly pot;
    pot.m = 1.0;
    pot.lambda = 0;
    
    BackgroundSolver(m, t0, t1, phi_p, dphi_p, pot);
    
    return 0;
}
