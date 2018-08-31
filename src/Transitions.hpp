#pragma once
#include "BackgroundSolver.hpp"
#include "linear_interpolation.hpp"

struct TransitionsSolution
{
    TransitionsSolution() :
    a{}, b{}, seg_control{} {}
    
    std::vector<double> a;
    std::vector<double> b;
    std::vector<std::pair<double, int>> seg_control;
};

struct Transitions
{
    double N_initial, N_final;
    LinearInterpolator<double, double> omega_2, log_aH;
    std::vector<double> N_extrema;
    
    Transitions(double _N_initial, double _N_final, LinearInterpolator<double, double> _omega_2, LinearInterpolator<double, double> _log_aH, std::vector<double> _N_extrema):
    N_initial(_N_initial), N_final(_N_final), omega_2{_omega_2}, log_aH{_log_aH}, N_extrema{_N_extrema} { }
    
    TransitionsSolution Find(double k, double error);
    std::vector<std::pair<double, int>> N_Distribution(double k, double N_i, double N_f, double lim);
    double True(double N,double k);
};

