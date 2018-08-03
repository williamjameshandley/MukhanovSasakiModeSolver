#include <math.h>
#include "BackgroundSolver.hpp"

double H(double phi, double dphi, Potential* pot)
{
    return sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double ddz(double phi, double dphi, double n, Potential* pot)
{
    return exp(2 * n) * (2 * H(phi, dphi, pot) * H(phi, dphi, pot)  - (7.0/2) * (dphi * dphi) - pot->ddV(phi) - 2 * dphi * pot->dV(phi) / H(phi, dphi, pot) + 0.5 * (dphi * dphi * dphi * dphi) / (H(phi, dphi, pot) * H(phi, dphi, pot)));
}

void equations(double dx_dt[], const double t, const double x[], void* data)
{
    auto pot = static_cast<Potential*>(data);
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1], pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = H(x[0], x[1], pot);
    dx_dt[3] = exp(-x[2]);
}

BackgroundSolution solve_equations(Potential* pot, double phi_p, double dphi_p)
{
    double t0 = 6.48757e-6 / sqrt(2 * pot->V(1)), dt = 1e-3;
    std::vector<double> x = {phi_p, dphi_p, 0, 1.5 * t0};
    
    Solutions sol;
    sol(x, t0);
    dlsodar desolver(4, 1);
    while(x[0] > 0)
    {
        desolver.integrate(t0, t0 + dt, &x[0], equations, pot);
        sol(x, t0);
        dt *= 1.0001;
    }
    
    std::vector<double> Z, DZ, DDZ, ETA;
    double a_end = exp(sol.x[sol.t.size() - 1][2]);
    
    for(size_t i = 0; i < sol.t.size(); i++)
    {
        Z.push_back(exp(sol.x[i][2]) * sol.x[i][1] / H(sol.x[i][0], sol.x[i][1], pot));
        DZ.push_back(exp(2 * sol.x[i][2]) * (-2 * sol.x[i][1] - pot->dV(sol.x[i][0]) / H(sol.x[i][0], sol.x[i][1], pot) + pow(sol.x[i][1], 3) / (2 * H(sol.x[i][0], sol.x[i][1], pot) * H(sol.x[i][0], sol.x[i][1], pot))));
        DDZ.push_back(ddz(sol.x[i][0], sol.x[i][1], sol.x[i][2], pot) / pow(a_end,2));
        ETA.push_back(sol.x[i][3] * a_end);
    }

    return BackgroundSolution(a_end, Z, DZ, DDZ, ETA);
}

