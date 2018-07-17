#include <math.h>
#include "BackgroundSolver.hpp"

template<class Integrator>
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> BackgroundSolver::Solve(Integrator integrator)
{
    double dt = (t1 - t0) / 1e7; //used only once at the start of integration
    double eta0 = 1.5 * t0, n0 = 0;
    
    std::vector<double> x0 = {phi_p, dphi_p, n0, eta0};

    Solutions sol;
    
    size_t steps = boost::numeric::odeint::integrate_adaptive(integrator, *this , x0, t0, t1, dt, boost::ref(sol));

    std::vector<double> Z, DZ, DDZ, ETA;
    
    for(size_t i=0; i<=steps; i++)
    {
        Z.push_back(exp(sol.x[i][2]) * sol.x[i][1] / H(sol.x[i][0], sol.x[i][1]));
        DZ.push_back(exp(2 * sol.x[i][2]) * (-2 * sol.x[i][1] - pot.dV(sol.x[i][0]) / H(sol.x[i][0], sol.x[i][1]) + std::pow(sol.x[i][1], 3) / (2 * H(sol.x[i][0], sol.x[i][1]) * H(sol.x[i][0], sol.x[i][1]))));
        DDZ.push_back(ddz(sol.x[i][0], sol.x[i][1], sol.x[i][2]));
        ETA.push_back(sol.x[i][3]);
    }
    
    return std::make_tuple(Z, DZ, DDZ, ETA);
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

