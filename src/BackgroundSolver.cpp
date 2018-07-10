#include <math.h>
#include "BackgroundSolver.hpp"

template<class Integrator>

std::tuple<std::vector<double>, std::vector<double>> BackgroundSolver::Solve(Integrator integrator)
{
    double t = t0;
    double dt = (t1 - t0) / 1e7;
    double eta0 = 1.5 * t0, n0 = 0;
    std::vector<double> x = {phi_p, dphi_p, n0, eta0};
    std::vector<double> DDZ;
    std::vector<double> ETA;
    
    int count = 0;
    
    while (t < t1)
    {
        if(count%1000 == 0)
        {
            DDZ.push_back(ddz(x[0], x[1], x[2]));
            ETA.push_back(x[3]);
        }
        count += 1;
        integrator(*this, x, t, dt);
    }

    return std::make_tuple(DDZ, ETA);
}


double BackgroundSolver::H(double phi, double dphi)
{
    double H = sqrt((0.5 * dphi * dphi + pot.V(phi)) / 3.0);
    return H;
}

double BackgroundSolver::ddz(double phi, double dphi, double n)
{
    double ddz = exp(2 * n) * (2 * H(phi, dphi) * H(phi, dphi)  - (7.0/2) * (dphi * dphi) - pot.ddV(phi) - 2 * dphi * pot.dV(phi) / H(phi, dphi) + 0.5 * (dphi * dphi * dphi * dphi) / (H(phi, dphi) * H(phi, dphi)));
    return ddz;
}

void BackgroundSolver::operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double)
{
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1]) * x[1] + pot.dV(x[0]));
    dx_dt[2] = H(x[0], x[1]);
    dx_dt[3] = exp(-x[2]);
}

BackgroundSolver::BackgroundSolver(double a, double b, double c, double d, Poly potential) {
    t0 = a;
    t1 = b;
    phi_p = c;
    dphi_p = d;
    pot = potential;
}
