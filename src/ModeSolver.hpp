#pragma once
#include "Special_Functions.hpp"
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
    
        std::vector<std::vector<Eigen::Matrix2d>> Mat;
        std::vector<std::complex<double>> c, d;
    
        ModeSolver(std::vector<double> kk, std::vector<double> ee, std::vector<double> aa, std::vector<double> bb, double dd, double endend, std::vector<double> zz, std::vector<double> ddzz, std::vector<double> ddddzz, std::vector<double> esol);
    
        void Find_Mat();
        void Initial_Conditions(std::string Vacuum, double eta_r);
    
        std::vector<double> PPS();
    
};
