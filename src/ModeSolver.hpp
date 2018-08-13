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
    
        double N_r, PPS_error;
        VacuumChoice vacuum;
        LinearInterpolator<double, double> Z, H, DPHI, H_aH, PPS;
    
        ModeSolver(BackgroundSolution _Bsol);
    
        Eigen::Matrix2d Mat();
        Eigen::Vector2cd Match(double k);
        void Initial_Conditions(VacuumChoice _vacuum, double _N_r);

        Eigen::Matrix2d Airy_mat(double x0, double x1, double p);
    
        double Find_PPS(double k);
        void Construct_PPS(double k0, double k1);
    
};
