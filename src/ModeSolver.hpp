#pragma once
#include <iostream>
#include <math.h>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "Transitions.hpp"
#include "BackgroundSolver.hpp"
#include "Special_Functions.hpp"
#include "linear_interpolation.hpp"

enum VacuumChoice { BD, HD, RST };

class ModeSolver
{
    public:
        BackgroundSolution Bsol;    
        TransitionsSolution Tsol;    
    
        double eta_r, PPS_error;
        VacuumChoice vacuum;
        size_t initial_index;
        LinearInterpolator<double, double> DDZ, DZ, Z, PPS;
    
        ModeSolver(BackgroundSolution _Bsol, double PPS_error);
    
        Eigen::Matrix2d Mat(double k);
        Eigen::Vector2cd Match(double k);
        void Initial_Conditions(VacuumChoice vacuum, double eta_r);

        Eigen::Matrix2d A(double x, double p);
        Eigen::Matrix2cd H(double eta, double k);
        Eigen::Matrix2cd H(double x, double k, double v);
        Eigen::Matrix2cd z_Mat(std::complex<double> x, double v);
    
        double Find_PPS(double k);
        void Construct_PPS(double k0, double k1);
    
};
