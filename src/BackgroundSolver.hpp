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
    BackgroundSolver():  pot{nullptr} {}

    template<class Integrator> BackgroundSolution Solve(Integrator, Potential&, double t0, double t1, double phi_p, double dphi_p);

    void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);


    private:
        double H(double phi, double dphi);
        double ddz(double phi, double dphi, double n);
        double V(double phi)   {return pot->V(phi);}
        double dV(double phi)  {return pot->dV(phi);}
        double ddV(double phi) {return pot->ddV(phi);}
        Potential* pot;
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


