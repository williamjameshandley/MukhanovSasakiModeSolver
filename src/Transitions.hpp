#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>

class Transitions
{
    public:
        double eta_i, eta_f, eta_end;
        std::vector<double> ddz, eta_sol;
    
        Transitions(double a, double b, double c, std::vector<double> d, std::vector<double> e):
            eta_i(a), eta_f(b), eta_end(c), ddz(d), eta_sol(e) { }

        double integral(double a, double b);
    
        std::tuple<std::vector<double>, std::vector<double>, double,  std::vector<double>> Find(double error);
    
};

