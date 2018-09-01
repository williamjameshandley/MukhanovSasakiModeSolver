#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include <math.h>
#include <stddef.h>
#include "Potential.hpp"
#include "linear_interpolation.hpp"

void equations(double dx_dt[], const double t, const double x[], void* data);
double dphi_H(const double x[], Potential* pot);
double H(const double x[], Potential* pot);
double omega_2(const double x[], Potential* pot);
double d_omega_2(const double x[], Potential* pot);

void inflation_end(double g[], const double, const double x[], void* data);
void inflation_begin(double g[], const double, const double x[], void* data);
void Find_N(double g[], const double, const double x[], void* data);
void Extrema_Scalar(double g[], const double, const double x[], void* data);
void Extrema_Tensor(double g[], const double, const double x[], void* data);

struct BackgroundSolution
{
    BackgroundSolution(LinearInterpolator<double, double> _omega_2, LinearInterpolator<double, double> _omega_2_tensor, LinearInterpolator<double, double> _log_aH, LinearInterpolator<double, double> _dphi_H, std::vector<double> _N_extrema, std::vector<double> _N_extrema_tensor, double _aH_star, double _N_end) :
    omega_2{_omega_2}, omega_2_tensor{_omega_2_tensor}, log_aH{_log_aH}, dphi_H{_dphi_H}, N_extrema{_N_extrema}, N_extrema_tensor{_N_extrema_tensor}, aH_star{_aH_star}, N_end{_N_end} {}
    
    LinearInterpolator<double, double> omega_2;
    LinearInterpolator<double, double> omega_2_tensor;
    LinearInterpolator<double, double> log_aH;
    LinearInterpolator<double, double> dphi_H;
    std::vector<double> N_extrema;
    std::vector<double> N_extrema_tensor;
    double aH_star;
    double N_end;
    
};


std::vector<double> Solve_N(double t0, std::vector<double> x0, void* ptrs[]);
BackgroundSolution solve_equations(Potential* pot, double N_star);
BackgroundSolution solve_equations(Potential* pot, double N_star, double N_dagger);
LinearInterpolator<double, double> Solve_Variable(double t0, std::vector<double> x0, std::function<double(const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double error);
