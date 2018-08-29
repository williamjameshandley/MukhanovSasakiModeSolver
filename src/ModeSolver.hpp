#pragma once
#include <math.h>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "Transitions.hpp"
#include "BackgroundSolver.hpp"
#include "Special_Functions.hpp"
#include "linear_interpolation.hpp"

enum VacuumChoice { BD, HD, RST };

class BasicModeSolver
{
    public:
        LinearInterpolator<double, double> PPS;
        std::vector<double> k_plot; 

        BasicModeSolver() : PPS{}, k_plot{} {};

        void Construct_PPS(double k0, double k1, double error);
        virtual double Find_PPS(double k) =0;
};


class ModeSolver : public BasicModeSolver
{
    public:
        BackgroundSolution Bsol;    
        TransitionsSolution Tsol;    
    
        double N_r, PPS_error;
        VacuumChoice vacuum;
        LinearInterpolator<double, double> OMEGA_2, Z, H, DPHI;
    
        ModeSolver(BackgroundSolution _Bsol);
    
        Eigen::Matrix2d Airy_Mat();
        Eigen::Matrix2cd Bessel_Mat();
        Eigen::Vector2cd Match(double k);
        void Initial_Conditions(VacuumChoice _vacuum, double _N_r);

        Eigen::Matrix2d Airy_gen(double p, double x1, double x0);
        Eigen::Matrix2cd Bessel_gen(double p, double x1, double x0);
    
        virtual double Find_PPS(double k) override;
};

class NumericModeSolver : public BasicModeSolver
{
    public:
        Potential* pot;
        double N_star, N_r;
        NumericModeSolver(Potential* _pot, double _N_star, double _N_r) : pot{_pot}, N_star{_N_star}, N_r{_N_r} {}
    
        virtual double Find_PPS(double k) override;
};
