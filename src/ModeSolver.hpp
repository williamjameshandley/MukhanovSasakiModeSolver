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
        std::vector<double> k_plot; 
        LinearInterpolator<double, double> OMEGA_2, Z, H, DPHI, PPS;
    
        ModeSolver(BackgroundSolution _Bsol);
    
        Eigen::Matrix2d Airy_Mat();
        Eigen::Matrix2cd Bessel_Mat();
        Eigen::Vector2cd Match(double k);
        void Initial_Conditions(VacuumChoice _vacuum, double _N_r);

        Eigen::Matrix2d Airy_gen(double p, double x1, double x0);
        Eigen::Matrix2cd Bessel_gen(double p, double x1, double x0);
    
        double Find_PPS(double k);
        void Construct_PPS(double k0, double k1, double error);
    
};
