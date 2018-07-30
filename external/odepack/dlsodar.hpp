#pragma once
#include <vector>
#include <functional>

using Field    = std::function<void(double qdot[], const double t, const double q[], double data[])>;
using Jacobian = std::function<void(double dfdq[], const double t, const double q[], double data[])>;
using Root     = std::function<void(double root[], const double t, const double q[], double data[])>;

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

        void _integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, double data[]);
    public:
        dlsodar(int, int); 

        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, double data[]);
        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, double data[]);
        void integrate(double &t, double tout, double q[], Field f_func, double data[]);

};

