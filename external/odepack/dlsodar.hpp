#pragma once
#include <vector>
#include <functional>

using Field    = void(*)(double qdot[], const double t, const double q[], void *data);
using Jacobian = void(*)(double dfdq[], const double t, const double q[], void *data);
using Root     = void(*)(double root[], const double t, const double q[], void *data);

static Field f_func_;
static Jacobian j_func_;
static Jacobian g_func_;
static void* data_;

void f_(const int *, const double *t_, const double *y, double *ydot);
void j_(const int *, const double *t_, const double *y, const int *, const int *, double *dfdy, const int *);
void g_(const int *, const double *t_, const double *y, const int *, double *gout);

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
        int ml;
        int mu;

        void _integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data);
    public:
        dlsodar(int, int, int=0); 

        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data);
        void integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, void *data);
        void integrate(double &t, double tout, double q[], Field f_func, void *data);

        void reset() {istate=1;}
        dlsodar& set_atol(double atol_) {atol={atol_}; return *this;}
        dlsodar& set_rtol(double rtol_) {rtol={rtol_}; return *this;}
        dlsodar& set_tol(double atol_, double rtol_) {return set_atol(atol_).set_rtol(rtol_);}

};

