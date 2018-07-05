#pragma once
#include <iostream>


double H(double phi, double dphi);
void Equations(const std::vector<double>& x, std::vector<double>& dx_dt, const double);
void BackgroundSolver(double m, double t0, double t1, double phi_p, double dphi_p, class pot);
