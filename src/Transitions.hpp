#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>
#include "BackgroundSolver.hpp"

struct TransitionsSolution
{
    TransitionsSolution() :
    a{}, b{}, N_step{} {}
    
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> N_step;
};

struct Transitions
{
   
    double N_i, N_f;
    BackgroundSolution Bsol;

    Transitions(double _N_i, double _N_f, BackgroundSolution _Bsol):
        N_i(_N_i), N_f(_N_f), Bsol{_Bsol} { }

    double integral(double a, double b, double k);

    TransitionsSolution Find(double k, double error);
    
};

