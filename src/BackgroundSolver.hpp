#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include "Potential.hpp"

void equations(double dx_dt[], const double t, const double x[], void* data);
double H(double phi, double dphi, Potential* pot);
double omega_2(double phi, double dphi, Potential* pot);
double d_omega_2(double phi, double dphi, Potential* pot);
void end(double g[], const double, const double x[], void* data);

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
    BackgroundSolution(std::vector<double> _t, std::vector<double> _phi, std::vector<double> _dphi,  std::vector<double> _N, std::vector<double> _H, std::vector<double> _z, std::vector<double> _omega_2, std::vector<double> _d_omega_2, double _aH_star) :
    t{_t}, phi{_phi}, dphi{_dphi}, N{_N}, H{_H}, z{_z}, omega_2{_omega_2}, d_omega_2{_d_omega_2}, aH_star{_aH_star} {}

    std::vector<double> t;
    std::vector<double> phi;
    std::vector<double> dphi;
    std::vector<double> N;
    std::vector<double> H;
    std::vector<double> z;
    std::vector<double> omega_2;
    std::vector<double> d_omega_2;
    double aH_star;
    
};

BackgroundSolution solve_equations(Potential* pot, double phi_p, double dphi_p);
