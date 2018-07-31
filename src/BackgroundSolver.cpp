#include <math.h>
#include "BackgroundSolver.hpp"

double H(double phi, double dphi, Potential* pot)
{
    return sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double ddz(double phi, double dphi, double n, Potential * pot)
{
    return exp(2 * n) * (2 * H(phi, dphi) * H(phi, dphi)  - (7.0/2) * (dphi * dphi) - pot->ddV(phi) - 2 * dphi * pot->dV(phi) / H(phi, dphi) + 0.5 * (dphi * dphi * dphi * dphi) / (H(phi, dphi) * H(phi, dphi)));
}

void solve_equations(Potential* pot, double t0, double t1, double phi_p, double dphi_p)
{
    dlsodar desolver(4, 1);
    double x[4] = {phi_p, dphi_p, 1.5*t0,0};
    desolver.integrate(t0, t1, x, equations, pot);

}

void equations(double dx_dt[], const double t, const double x[], void* data)
{
    auto pot = static_cast<Potential*>(data);
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1]) * x[1] + pot->dV(x[0]));
    dx_dt[2] = H(x[0], x[1]);
    dx_dt[3] = exp(-x[2]);
}

