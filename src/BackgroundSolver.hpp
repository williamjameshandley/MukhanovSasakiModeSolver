#pragma once
#include <iostream>
#include "Potential.hpp"


double H(double phi, double dphi);
void Equations(const std::vector<double>& x, std::vector<double>& dx_dt, const double);
void BackgroundSolver(double t0, double t1, double phi_p, double dphi_p, Poly pot);
