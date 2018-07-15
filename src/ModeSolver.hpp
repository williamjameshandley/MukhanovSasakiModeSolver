#pragma once
#include <boost/math/special_functions/airy.hpp>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>

class ModeSolver
{
    public:
    
        std::vector<double> k, eta_step, a, b;
        double delta, eta_end;
        std::vector<double> z, dz, ddz;
    
    std::vector<std::vector<Eigen::MatrixXd>> Mats;
    
        ModeSolver(std::vector<double> kk, std::vector<double> ee, std::vector<double> aa, std::vector<double> bb, double dd, double endend, std::vector<double> zz, std::vector<double> ddzz, std::vector<double> ddddzz):
            k(kk), eta_step(ee), a(aa), b(bb), delta(dd), eta_end(endend), z(zz), dz(ddzz), ddz(ddddzz) {}
    
        void Find_Mat();
        void initial_conditions(double eta_r, char Vacuum);
    
        std::vector<double> Solver();
    
};
