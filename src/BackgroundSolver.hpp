#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include <math.h>
#include <stddef.h>
#include <limits>
#include "Potential.hpp"
#include "interpolation.hpp"

void equations(double dx_dt[], const double t, const double x[], void* data);
double H(const double n, const double x[], Potential* pot);
double dphi_H(const double n, const double x[], Potential* pot);
double log_aH(const double n, const double x[], Potential* pot);
double dlog_aH(const double n, const double x[], Potential* pot);
double omega_2(const double n, const double x[], Potential* pot);
double d_omega_2(const double n, const double x[], Potential* pot);
double omega_2_tensor(const double n, const double x[], Potential* pot);
double d_omega_2_tensor(const double n, const double x[], Potential* pot);

void inflation_end(double g[], const double, const double x[], void* data);
void inflation_begin(double g[], const double, const double x[], void* data);
void Find_N(double g[], const double, const double x[], void* data);
void Extrema_Scalar(double g[], const double, const double x[], void* data);
void Extrema_Tensor(double g[], const double, const double x[], void* data);

struct BackgroundSolution
{
    LinearInterpolator<double, double> omega_2;
    LinearInterpolator<double, double> omega_2_tensor;
    LinearInterpolator<double, double> log_aH;
    LinearInterpolator<double, double> dphi_H;
    std::vector<double> N_extrema;
    std::vector<double> N_extrema_tensor;
    double aH_star;
    double N_end;
};


BackgroundSolution solve_equations(double lim, Potential* pot, double N_star, double N_dagger=20);

LinearInterpolator<double, double> Solve_Variable(double t0, double N_end, std::vector<double> x0, std::function<double(const double n, const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double error);
