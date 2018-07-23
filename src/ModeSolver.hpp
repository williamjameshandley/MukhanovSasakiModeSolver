#pragma once
#include "Special_Functions.hpp"
#include "linear_interpolation.hpp"
#include "BackgroundSolver.hpp"
#include "Transitions.hpp"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <complex>
#include <cmath>


class ModeSolver
{
    public:
        BackgroundSolution Bsol;
        TransitionsSolution Tsol;
    
        double eta_r;
        std::string Vacuum;
        size_t initial_index;
        LinearInterpolator<double, double> DDZ, DZ, Z;
    
        Eigen::Matrix2d Mat;
        std::complex<double> c, d;
    
        ModeSolver(BackgroundSolution _Bsol, TransitionsSolution _Tsol):
            Bsol{_Bsol}, Tsol{_Tsol}, Mat{}, c{}, d{} {}
    
        void Find_Mat(double k);
        void Initial_Conditions(std::string Vacuum, double eta_r);
        void Match(double k);
    
        double PPS(double k);
    
};
