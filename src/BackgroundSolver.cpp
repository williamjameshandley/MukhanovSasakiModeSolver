#include <math.h>
#include "BackgroundSolver.hpp"

double H(double phi, double dphi, Potential* pot)
{
    return sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double omega_2(double phi, double dphi, Potential* pot)
{
    return (pot->ddV(phi) + 1.5 * pot->dV(phi) * dphi / H(phi, dphi, pot)) / pow(H(phi, dphi, pot), 2) - (pow(dphi/H(phi, dphi, pot), 2) - 6) * (5 * pow(dphi/H(phi, dphi, pot), 2) - 6) / 16.0;
}

double d_omega_2(double phi, double dphi, Potential* pot)
{
    return -3 * pow(pot->dV(phi) / H(phi, dphi, pot) / H(phi, dphi, pot), 2) - 9 * pot->dV(phi) * dphi / pow(H(phi, dphi, pot), 3) + (7.0/2) * pot->dV(phi) * pow(dphi, 3) / pow(H(phi, dphi, pot), 5) + (5.0/2) * pot->ddV(phi) * pow(dphi / H(phi, dphi, pot) / H(phi, dphi, pot), 2) - (27.0/2) * pow(dphi / H(phi, dphi, pot), 2) + 6 * pow(dphi / H(phi, dphi, pot), 4) - (5.0/8) * pow(dphi / H(phi, dphi, pot), 6) + pot->dddV(phi) * dphi / pow(H(phi, dphi, pot), 3);
}

void equations(double dx_dt[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x[0], x[1], pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = H(x[0], x[1], pot);
    dx_dt[3] = exp(-3 * x[2]) * H(x[0], x[1], pot) * H(x[0], x[1], pot) / x[1] / x[1];
}

void end(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    
    auto p = -(2.5 * x[1] * x[1] + x[1] * pot->dV(x[0]) / H(x[0], x[1], pot) - pow(0.5 * x[1] * x[1] / H(x[0], x[1], pot), 2));
    if(p > 0)
        g[0] = (-H(x[0], x[1], pot) + 0.5 * x[1] * x[1] / H(x[0], x[1], pot));
    else if(p < 0)
        g[0] = p;
}

void N_star_check(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto params = static_cast<double*> (ptr[1]);

    g[0] = (params[0] - params[1]) - x[2];
}

BackgroundSolution solve_equations(Potential* pot, double phi_p, double dphi_p)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    
    double t0 = 1;
    
    //Find N_end
    double t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0, 0};
    dlsodar desolver(4, 1, 10000);
    desolver.integrate(t, 1e10, &x[0], equations, end, static_cast<void*> (ptrs));
    double N_end = x[2];
    
    std::cout<<N_end<<std::endl;
    
    //Find aH_star
    double params[2];
    params[0] = N_end;
    params[1] = 55;
    ptrs[1] = static_cast<void*> (params);
    t = t0;
    x = {phi_p, dphi_p, 0, 0};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], equations, N_star_check, static_cast<void*> (ptrs));
    double aH_star = exp(x[2]) * H(x[0], x[1], pot);

    //Run Solver
    double dt = 1e-2;
    t = t0;
    x = {phi_p, dphi_p, 0, 0};
    Solutions sol;
    sol(x, t0);
    desolver = dlsodar(4, 1, 100000);
    while(x[2] < N_end)
    {
        desolver.integrate(t0, t0 + dt, &x[0], equations, end, static_cast<void*> (ptrs));
        sol(x, t0);
        dt *= 1.0001;
    }
    
    std::vector<double> _phi, _dphi, _N, _H, _z, _omega_2, _d_omega_2;
    
    for(size_t i = 0; i < sol.t.size(); i++)
    {
        _phi.push_back(sol.x[i][0]);
        _dphi.push_back(sol.x[i][1]);
        _N.push_back(sol.x[i][2]);
        _z.push_back(exp(sol.x[i][2]) * sol.x[i][1] / H(sol.x[i][0], sol.x[i][1], pot));
        _omega_2.push_back(omega_2(sol.x[i][0], sol.x[i][1], pot));
        _d_omega_2.push_back(d_omega_2(sol.x[i][0], sol.x[i][1], pot));
        _H.push_back(H(sol.x[i][0], sol.x[i][1], pot));
        
    }
    

    return BackgroundSolution(sol.t, _phi, _dphi, _N, _H, _z, _omega_2, _d_omega_2, aH_star);
}

