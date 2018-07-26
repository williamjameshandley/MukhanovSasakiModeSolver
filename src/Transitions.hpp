#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>
#include "BackgroundSolver.hpp"

struct TransitionsSolution
{
    TransitionsSolution() = default ;
    
    double delta;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> eta_step;
};

struct Transitions
{
   
    double eta_i, eta_f;
    BackgroundSolution Bsol;

    Transitions(double _eta_i, double _eta_f, BackgroundSolution _Bsol):
        eta_i(_eta_i), eta_f(_eta_f), Bsol{_Bsol} { }

    double integral(double a, double b);

    TransitionsSolution Find(double error);
    
};

