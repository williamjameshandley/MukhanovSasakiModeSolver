#pragma once
#include "Special_Functions.hpp"
#include "linear_interpolation.hpp"
#include "BackgroundSolver.hpp"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <complex>
#include <cmath>


class ModeSolver
{
    public:
    
        std::vector<double> eta_step, a, b;
        double delta, eta_end, eta_r;
        BackgroundSolution sol;
        std::string Vacuum;
        LinearInterpolator<double, double> DDZ, DZ, Z;
    
        Eigen::Matrix2d Mat;
        size_t initial_index;
        std::complex<double> c, d;
    
        ModeSolver(std::vector<double> ee, std::vector<double> aa, std::vector<double> bb, double dd, double endend, BackgroundSolution sol);
    
        void Find_Mat(double k);
        void Initial_Conditions(std::string Vacuum, double eta_r);
        void Match(double k);
    
        double PPS(double k);
    
};
