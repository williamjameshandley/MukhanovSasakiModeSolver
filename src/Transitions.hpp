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
    
    double N_initial, N_final;
    BackgroundSolution Bsol;
    
    Transitions(double _N_initial, double _N_final, BackgroundSolution _Bsol):
    N_initial(_N_initial), N_final(_N_final), Bsol{_Bsol} { }
    
    double Find_N_log_end(double k, double N_i);
    TransitionsSolution Find(double k, double error);
    std::vector<double> Linear(double k, double N_i, double N_f, double lim);
    std::vector<double> Log(double k, double N_i, double N_f, double lim);
    
    double True(double N,double k);
    
};

