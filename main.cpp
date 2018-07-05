#include <iostream>
#include "src/BackgroundSolver.cpp"


int main()
{
    double t0 = 1e-6;
    double t1 = 20.0;
    double phi_p = 23.0, dphi_p = -sqrt(2.0/3.0) / t0;
    
    BackgroundSolver(t0, t1, phi_p, dphi_p);
    
    return 0;
}
