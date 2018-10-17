#pragma once
#include <iostream>
#include <vector>
#include <odepack/dlsodar.hpp>
#include <math.h>
#include <stddef.h>
#include <limits>
#include "Potential.hpp"
#include "interpolation.hpp"

void equations_n(double dx_dn[], const double n, const double x[], void* data);
void equations_t(double dx_dn[], const double n, const double x[], void* data);
double H(const double n, const double x[], Potential* pot);
double phi_dot_H(const double n, const double x[], Potential* pot);
double aH(const double n, const double x[], Potential* pot);
double dlog_aH(const double n, const double x[], Potential* pot);
double omega_2(const double n, const double x[], Potential* pot);
double d_omega_2(const double n, const double x[], Potential* pot);
double omega_2_tensor(const double n, const double x[], Potential* pot);
double d_omega_2_tensor(const double n, const double x[], Potential* pot);

void inflating(double g[], const double, const double x[], void* data);
void inflating_t(double g[], const double, const double x[], void* data);
void Extrema_Scalar(double g[], const double, const double x[], void* data);
void Extrema_Tensor(double g[], const double, const double x[], void* data);

struct BackgroundSolution
{
    SemiLogInterpolator<double, double> omega_2;
    SemiLogInterpolator<double, double> omega_2_tensor;
    SemiLogInterpolator<double, double> aH;
    SemiLogInterpolator<double, double> dphi_H;
    std::vector<double> N_extrema;
    std::vector<double> N_extrema_tensor;
    double aH_star;
    double N_end;
};


BackgroundSolution solve_equations(double lim, Potential* pot, double, double, double);

SemiLogInterpolator<double, double> Solve_Variable(double t0, double N_end, std::vector<double> x0, std::function<double(const double n, const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double error);
