#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include "Potential.hpp"

void equations(double dx_dt[], const double t, const double x[], void* data);
double H(double phi, double dphi, Potential* pot);
double ddz(double phi, double dphi, double n, Potential* pot);

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

struct BackgroundSolution
{
    BackgroundSolution(std::vector<double> _z, std::vector<double> _dz,std::vector<double>_ddz, std::vector<double> _eta) :
    z{_z}, dz{_dz}, ddz{_ddz}, eta{_eta} {}

    std::vector<double> z;
    std::vector<double> dz;
    std::vector<double> ddz;
    std::vector<double> eta;
};

BackgroundSolution solve_equations(Potential* pot, double t0, double phi_p, double dphi_p);
