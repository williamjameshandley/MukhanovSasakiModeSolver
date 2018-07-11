#pragma once
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include "Potential.hpp"
#include <vector>

class BackgroundSolver
{
    public:
        double t0;
        double t1;
        double phi_p;
        double dphi_p;
        double m;
        double lambda;
        Poly pot;
    
        BackgroundSolver(double a, double b, double c, double d, Poly potential);
    
        template<class Integrator>
        std::tuple<std::vector<double>, std::vector<double>> Solve(Integrator);

        double H(double phi, double dphi);
        double ddz(double phi, double dphi, double n);
        void operator() (const std::vector<double>& x, std::vector<double>& dx_dt, const double);
    
};

struct push_back_state_and_time
{
    std::vector< std::vector<double> >& m_states;
    std::vector< double >& m_times;
    
    push_back_state_and_time( std::vector< std::vector<double> > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }
    
    void operator()( const std::vector<double> &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
