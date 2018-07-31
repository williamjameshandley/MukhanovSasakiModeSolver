#include <math.h>
#include "BackgroundSolver.hpp"

double BackgroundSolver::H(double phi, double dphi)
{
    return sqrt((0.5 * dphi * dphi + V(phi)) / 3.0);
}

double BackgroundSolver::ddz(double phi, double dphi, double n)
{
    return exp(2 * n) * (2 * H(phi, dphi) * H(phi, dphi)  - (7.0/2) * (dphi * dphi) - ddV(phi) - 2 * dphi * dV(phi) / H(phi, dphi) + 0.5 * (dphi * dphi * dphi * dphi) / (H(phi, dphi) * H(phi, dphi)));
}

void BackgroundSolver::Equations(std::vector<double>& dx_dt, const double, const std::vector<double>& x, double params[])
{
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1]) * x[1] + dV(x[0]));
    dx_dt[2] = H(x[0], x[1]);
    dx_dt[3] = exp(-x[2]);
}

