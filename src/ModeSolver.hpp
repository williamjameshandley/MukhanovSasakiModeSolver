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
        LinearInterpolator<double, double> PPS_Scalar, PPS_Tensor;
        std::vector<double> k_plot_Scalar, k_plot_Tensor; 

        BasicModeSolver() : PPS_Scalar{}, PPS_Tensor{}, k_plot_Scalar{}, k_plot_Tensor{} {};

        void Construct_PPS_Scalar(double k0, double k1, double error);
        void Construct_PPS_Tensor(double k0, double k1, double error);
        virtual double Find_PPS_Scalar(double k) =0;
        virtual double Find_PPS_Tensor(double k) =0;
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
    
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;
};

class NumericModeSolver : public BasicModeSolver
{
    public:
        Potential* pot;
        double N_star, N_dagger, N_r, phi_p, dphi_p, log_aH_star, phi_IC, dphi_IC, n_IC, N_end;
        NumericModeSolver(Potential* _pot, double _N_star);
        NumericModeSolver(Potential* _pot, double _N_star, double _N_r);
        NumericModeSolver(Potential* _pot, double _N_star, double N_dagger, double _N_r);
    
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;
};
