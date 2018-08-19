#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>
#include "BackgroundSolver.hpp"
#include "linear_interpolation.hpp"

struct TransitionsSolution
{
    TransitionsSolution() :
    lin_a{}, lin_b{}, lin_N_step{}, log_a{}, log_b{}, log_N_step{} {}
    
    std::vector<double> lin_a;
    std::vector<double> lin_b;
    std::vector<double> lin_N_step;
    std::vector<double> log_a;
    std::vector<double> log_b;
    std::vector<double> log_N_step;
};

struct Transitions
{
   
    double N_i, N_f;
    BackgroundSolution Bsol;

    Transitions(double _N_i, double _N_f, BackgroundSolution _Bsol):
        N_i(_N_i), N_f(_N_f), Bsol{_Bsol} { }

    TransitionsSolution Find(double k, double error);
    std::vector<double> Linear(double N_initial, double N_final, LinearInterpolator<double, double> True, double lim);
    std::vector<double> Log(double N_initial, double N_final, LinearInterpolator<double, double> True, double lim);
    
};

