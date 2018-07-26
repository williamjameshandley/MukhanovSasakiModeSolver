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

enum VacuumChoice { BD, HD, RST };

class ModeSolver
{
    public:
        BackgroundSolution Bsol;
        TransitionsSolution Tsol;
    
        double eta_r;
        VacuumChoice vacuum;
        size_t initial_index;
        LinearInterpolator<double, double> DDZ, DZ, Z;
    
        ModeSolver(BackgroundSolution _Bsol, TransitionsSolution _Tsol);
    
        Eigen::Matrix2d Mat(double k);
        Eigen::Vector2cd Match(double k);
        void Initial_Conditions(VacuumChoice vacuum, double eta_r);

        Eigen::Matrix2d A(double x, double p);
        Eigen::Matrix2cd H(double eta, double k);
        Eigen::Matrix2cd H(double x, double k, double v);
        Eigen::Matrix2cd a(std::complex<double> x, double v);
    
        double PPS(double k);
    
};
