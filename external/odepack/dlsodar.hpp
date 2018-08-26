#pragma once
#include <vector>
#include <functional>

using Field    = void(*)(double qdot[], const double t, const double q[], void *data);
using Jacobian = void(*)(double dfdq[], const double t, const double q[], void *data);
using Root     = void(*)(double root[], const double t, const double q[], void *data);

class dlsodar
{
    private:
        int neq;
        int itol;
        std::vector<double> rtol;
        std::vector<double> atol;
        int itask;
        int istate;
        int iopt;
        std::vector<double> rwork;
        int lrw;
        std::vector<int> iwork;
        int liw;
        int jt;
        int ng;
        std::vector<int> jroot;
        int max_steps;

        void _integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data);
    public:
        dlsodar(int, int, int=0); 

        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data);
        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, void *data);
        void integrate(double &t, double tout, double q[], Field f_func, void *data);

};

