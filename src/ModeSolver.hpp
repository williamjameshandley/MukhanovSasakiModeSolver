#pragma once
#include <math.h>
#include <cmath>
#include <complex>
#include <memory>
#include <Eigen/Dense>
#include "BackgroundSolver.hpp"
#include "interpolation.hpp"
#include "cephes.hpp"

const std::complex<double> I(0, 1);

enum VacuumChoice { BD, HD, RST, discrete};
enum TransitionChoice {neg_exp, lin, pos_exp};
enum PSChoice {scalar, tensor};

class BasicModeSolver
{
    public:
        LinearInterpolator<double, double> PPS_Scalar, PPS_Tensor;
        std::vector<double> k_plot_Scalar, k_plot_Tensor; 

        BasicModeSolver() : PPS_Scalar{}, PPS_Tensor{}, k_plot_Scalar{}, k_plot_Tensor{} {};

        void Construct_PPS_Scalar(double k0, double k1, double error);
        void Construct_PPS_Tensor(double k0, double k1, double error);

        virtual ~BasicModeSolver() {};
        virtual double Find_PPS_Scalar(double k) =0;
        virtual double Find_PPS_Tensor(double k) =0;
};


class ModeSolver : public BasicModeSolver
{
    public:
        BackgroundSolution Bsol;
        
        double N_r, PPS_error;
        VacuumChoice vacuum;
        
        ModeSolver(BackgroundSolution _Bsol);
        
        void Initial_Conditions(VacuumChoice _vacuum, double _N_r);
        Eigen::Vector2cd Initial_Q(double k);
        Eigen::Vector2cd Evolve(Eigen::Vector2cd Q, PSChoice _PSChoice, double k, double N_i, double& N_f);
        
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;
        
        double w_2(double N, double k, PSChoice _PSChoice);
        Eigen::Matrix2d Trig_Mat(double w2, double N0, double N1);
        Eigen::Matrix2d Airy_Mat(double a, double b, double N0, double N1);
        Eigen::Matrix2d Bessel_Mat(double a, double b, double N0, double N1);
        Eigen::Matrix2d Modified_Bessel_Mat(double a, double b, double N0, double N1);
        
        Eigen::Matrix2d Airy_gen(double p, double x1, double x0);
        Eigen::Matrix2d Bessel_gen(double p, double x1, double x0);
        Eigen::Matrix2d Modified_Bessel_gen(double p, double x1, double x0);

        Eigen::MatrixXd cst_step(double w2, double N_i, double N_f);
        Eigen::MatrixXd lin_step(double w_2_i, double w_2_f, double N_i, double N_f);
        Eigen::MatrixXd pos_exp_step(double w_2_i, double w_2_f, double N_i, double N_f);
        Eigen::MatrixXd neg_exp_step(double w_2_i, double w_2_f, double N_i, double N_f);
    
        struct Transition {
            enum type {cst, lin, pos, neg};
            double w2;
            Eigen::Matrix2d M_cst, M_lin, M_pos, M_neg;
            type i;
        };

};

double _H(double phi, double dphi, Potential* pot);
double dz_z(double phi, double dphi, Potential* pot);
void __equations(double dx_dt[], const double, const double x[], void* data);
void _equations(double dx_dt[], const double, const double x[], void* data);
void _equations_tensor(double dx_dt[], const double, const double x[], void* data);
void _inflation_end(double g[], const double, const double x[], void* data);
void _inflation_begin(double g[], const double, const double x[], void* data);
void _Find_N(double g[], const double, const double x[], void* data);
void _start(double g[], const double, const double x[], void* data);
void _finish(double g[], const double, const double x[], void* data);

class NumericModeSolver : public BasicModeSolver
{
    public:
        Potential* pot;
        double N_star, N_dagger, N_r, phi_p, dphi_p, log_aH_star, phi_IC, dphi_IC, n_IC, N_end;
        NumericModeSolver(Potential* _pot, double _N_star): 
            pot{_pot}, N_star{_N_star}, N_dagger{}, N_r{}, phi_p{}, dphi_p{}, log_aH_star{}, phi_IC{}, dphi_IC{}, n_IC{}, N_end{}
        {}
        NumericModeSolver(Potential* _pot, double _N_star, double _N_r);
        NumericModeSolver(Potential* _pot, double _N_star, double N_dagger, double _N_r);
        virtual ~NumericModeSolver() {};
    
        double Find_PPS(double k);
        virtual double Find_PPS_Scalar(double k) override;
        virtual double Find_PPS_Tensor(double k) override;

};

inline double Q_err(Eigen::Vector2cd Q,Eigen::Vector2cd Q0)
{
    auto abserr = std::max( abs(Q0[0]-Q[0]), abs(Q0[1]-Q[1]));
    auto relerr = std::max( abs(Q0[0]-Q[0])/abs(Q0[0]), abs(Q0[1]-Q[1])/abs(Q0[1]));
    //auto moderr = std::max( abs((abs(Q0[0])-abs(Q[0]))/abs(Q0[0])), abs((abs(Q0[1])-abs(Q[1]))/abs(Q0[1])));
    return std::max( abserr, relerr);
}
//{return std::max( abs(Q0[0]-Q[0])/abs(Q0[0]), abs(Q0[1]-Q[1])/abs(Q0[1]));}
//{return std::max( abs((abs(Q0[0])-abs(Q[0]))/abs(Q0[0])), abs((abs(Q0[1])-abs(Q[1]))/abs(Q0[1])));}
