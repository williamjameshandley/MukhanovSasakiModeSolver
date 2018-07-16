#pragma once
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "linear_interpolation.hpp"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <complex>
#include <cmath>

class ModeSolver
{
    public:
    
        std::vector<double> k, eta_step, a, b;
        double delta, eta_end;
        std::vector<double> z, dz, ddz, eta_sol;
    
        std::vector<std::vector<Eigen::MatrixXd>> Mat;
        std::vector<std::complex<double>> C, D;
    
        ModeSolver(std::vector<double> kk, std::vector<double> ee, std::vector<double> aa, std::vector<double> bb, double dd, double endend, std::vector<double> zz, std::vector<double> ddzz, std::vector<double> ddddzz, std::vector<double> esol):
            k(kk), eta_step(ee), a(aa), b(bb), delta(dd), eta_end(endend), z(zz), dz(ddzz), ddz(ddddzz), eta_sol(esol) {}
    
        void Find_Mat();
        void Initial_Conditions(std::string Vacuum, double eta_r);
    
        std::vector<double> PPS();
    
};
