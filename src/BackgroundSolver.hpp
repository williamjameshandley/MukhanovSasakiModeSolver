#pragma once
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include "Potential.hpp"
#include <vector>

struct BackgroundSolution
{
    BackgroundSolution(std::vector<double>_z,std::vector<double>_dz,std::vector<double>_ddz,std::vector<double>_eta) :
        z{_z}, dz{_dz}, ddz{_ddz}, eta{_eta} {}

    std::vector<double> z; 
    std::vector<double> dz; 
    std::vector<double> ddz; 
    std::vector<double> eta;
};

struct BackgroundSolver
{
    double t0, t1, phi_p, dphi_p;
    Poly pot;

    BackgroundSolver(double a, double b, double c, double d, Poly potential): 
        t0(a), t1(b), phi_p(c), dphi_p(d), pot(potential) {}

    template<class Integrator>
    BackgroundSolution Solve(Integrator);

    double H(double phi, double dphi);
    double ddz(double phi, double dphi, double n);
    void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);
};

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


