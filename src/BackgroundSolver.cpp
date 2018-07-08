#include <fstream>
#include <math.h>
#include "ascent/Ascent.h"
#include "Potential.hpp"
#include "BackgroundSolver.hpp"

void BackgroundSolver::Solve()
{
    double t = t0;
    double dt = (t1 - t0) / 1e7;
    double eta0 = 1.5 * t0, n0 = 0;
    std::vector<double> x = {phi_p, dphi_p, eta0, n0};
    
    asc::RK4 integrator;
    Equations eqs;
    pot.m = m;
    pot.lambda = lambda;
    
    std::ofstream fout;
    fout.open ("bin/output/phi.txt");
    
    int count = 0;
    while (t < t1)
    {
        count += 1;
        if(count % 1000 == 0)
        {
            fout<<t<<"   "<<x[0]<<std::endl;
        }
        integrator(eqs, x, t, dt);
    }
    fout.close();

}


double BackgroundSolver::Equations::H(double phi, double dphi)
{
    double H = sqrt((0.5 * dphi * dphi + pot.V(phi)) / 3.0);
    return H;
}


void BackgroundSolver::Equations::operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double)
{
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1]) * x[1] + pot.dV(x[0]));
    dx_dt[2] = H(x[0], x[1]);
    dx_dt[3] = exp(-x[2]);
}
