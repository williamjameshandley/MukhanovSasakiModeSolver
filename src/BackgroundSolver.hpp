#pragma once
#include "ascent/Ascent.h"
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
    
        template<class Integrator>
        std::tuple<std::vector<double>, std::vector<double>> Solve(Integrator);

        double H(double phi, double dphi);
        double ddz(double phi, double dphi, double n);
        void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);
    
};
