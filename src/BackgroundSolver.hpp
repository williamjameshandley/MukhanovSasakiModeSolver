#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include <boost/numeric/odeint.hpp>
#include "Potential.hpp"

void equations(double dx_dt[], const double t, const double x[], void* data);


double H(double phi, double dphi);
double ddz(double phi, double dphi, double n);

struct Solutions
{
    std::vector< std::vector<double> > x;
    std::vector< double > t;

    Solutions() : x{}, t{} {}
  
    void operator()( const std::vector<double> &x_ , double t_ )
    {
        x.push_back( x_ );
        t.push_back( t_ );
    }
};


//template<class Integrator>
//BackgroundSolution BackgroundSolver::Solve(Integrator integrator, Potential& potential, double t0, double t1, double phi_p, double dphi_p)
//{
//    pot = &potential;
//    double dt = (t1 - t0) / 1e6; //used only once at the start of integration
//    double eta0 = 1.5 * t0, n0 = 0;
//    
//    std::vector<double> x = {phi_p, dphi_p, n0, eta0};
//
//    Solutions sol;
//    
//    dlsodar desolver(4, 1);
//    
//    while(x[0] > 0){
//        desolver.integrate(t0, t0 + dt, &x[0], Equations, nullptr);
//        sol(x, t0);
//    }
//    
//    //size_t steps = boost::numeric::odeint::integrate_adaptive(integrator, *this , x0, t0, t1, dt, boost::ref(sol));
//
//    std::vector<double> Z, DZ, DDZ, ETA;
//    
//    for(size_t i = 0; i <= sol.t.size(); i++)
//    {
//        Z.push_back(exp(sol.x[i][2]) * sol.x[i][1] / H(sol.x[i][0], sol.x[i][1]));
//        DZ.push_back(exp(2 * sol.x[i][2]) * (-2 * sol.x[i][1] - dV(sol.x[i][0]) / H(sol.x[i][0], sol.x[i][1]) + std::pow(sol.x[i][1], 3) / (2 * H(sol.x[i][0], sol.x[i][1]) * H(sol.x[i][0], sol.x[i][1]))));
//        DDZ.push_back(ddz(sol.x[i][0], sol.x[i][1], sol.x[i][2]));
//        ETA.push_back(sol.x[i][3]);
//    }
//    pot = nullptr;
//    
//    return BackgroundSolution(Z, DZ, DDZ, ETA);
//}
