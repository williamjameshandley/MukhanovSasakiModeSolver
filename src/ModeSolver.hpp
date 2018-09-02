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
    
        ModeSolver(BackgroundSolution _Bsol);
    
        void Initial_Conditions(VacuumChoice _vacuum, double _N_r);
        Eigen::Vector2cd Match(double k);
    
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;
    
        Eigen::Matrix2d Airy_Mat(double a, double b, double N0, double N1);
        Eigen::Matrix2d Bessel_Mat(double a, double b, double N0, double N1);
        Eigen::Matrix2d Modified_Bessel_Mat(double a, double b, double N0, double N1);
    
        Eigen::Matrix2d Airy_gen(double p, double x1, double x0);
        Eigen::Matrix2d Bessel_gen(double p, double x1, double x0);
        Eigen::Matrix2d Modified_Bessel_gen(double p, double x1, double x0);
};

class NumericModeSolver : public BasicModeSolver
{
    public:
        Potential* pot;
        double N_star, N_dagger, N_r, phi_p, dphi_p, log_aH_star, phi_IC, dphi_IC, n_IC, N_end;
        NumericModeSolver(Potential* _pot, double _N_star): pot{_pot}, N_star{_N_star} {}
        NumericModeSolver(Potential* _pot, double _N_star, double _N_r);
        NumericModeSolver(Potential* _pot, double _N_star, double N_dagger, double _N_r);
    
        virtual double Find_PPS(double k);
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;
};
