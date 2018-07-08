#pragma once
#include <iostream>
#include "Potential.hpp"
#include <vector>

class BackgroundSolver
{
    public:
        double t0;
        double t1;
        double phi_p;
        double dphi_p;
        double m;
        double lambda;
        Poly pot;
    
        void Solve();
        
        class Equations
        {
            public:
                double H(double phi, double dphi);
                void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);
        };
    
};
