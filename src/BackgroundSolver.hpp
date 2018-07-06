#pragma once
#include <iostream>
#include "Potential.hpp"
#include <vector>


double H(double phi, double dphi);
void BackgroundSolver(double t0, double t1, double phi_p, double dphi_p, Poly pot);
class Equations
{
    public:
        Poly pot;
        double H(double phi, double dphi);
        void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);
};
