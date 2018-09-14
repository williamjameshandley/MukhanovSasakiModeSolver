#include "opkdmain.h"
#include "dlsodar.hpp"
#include "pointers.hpp"
#include <algorithm>

dlsodar::dlsodar(int neq_, int ng_, int max_steps_): 
        neq{neq_},
        itol{1},
        rtol{1e-15},
        atol{1e-15},
        itask{1},
        istate{1},
        iopt{1},
        rwork{},
        lrw{},
        iwork{},
        liw{},
        jt{2},
        ng{ng_},
        jroot(static_cast<size_t>(ng)),
        ml{0},
        mu{0}
{

    lrw = std::max({ 20 + 16*neq + 3*ng,
                     22 + 9*neq + neq_*neq + 3*ng,
                     22 + 10*neq + (2*ml+mu)*neq + 3*ng});
    rwork.resize(static_cast<size_t>(lrw),0);

    liw = 20 + neq;
    iwork.resize(static_cast<size_t>(liw),0);
    iwork[0] = ml;
    iwork[1] = mu;
    iwork[5] = max_steps_;
}

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, void *data)
{ jt = 2; int ng_ = ng; ng = 0; _integrate(t, tout, q, f_func, nullptr, nullptr, data); ng = ng_; }

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, Root g_func, void *data)
{ jt = 2; _integrate(t, tout, q, f_func, nullptr, g_func, data); }

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data)
{ jt = 1; _integrate(t, tout, q, f_func, j_func, g_func, data); }

void f_(const int *, const double *t_, const double *y, double *ydot)                                           {f_func_(ydot, *t_, y, data_); }
void j_(const int *, const double *t_, const double *y, const int *, const int *, double *dfdy, const int *)    {j_func_(dfdy, *t_, y, data_);}
void g_(const int *, const double *t_, const double *y, const int *, double *gout)                              {g_func_(gout, *t_, y, data_);}


void dlsodar::_integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data){

    data_ = data;
    f_func_ = f_func;
    j_func_ = j_func;
    g_func_ = g_func;

    dlsodar_(
            &f_,
            &neq, q, &t, &tout,
            &itol, &rtol[0], &atol[0],
            &itask, &istate, &iopt,
            &rwork[0], &lrw, &iwork[0], &liw,
            &j_, &jt,
            &g_, &ng,
            &jroot[0]
            );
    // @todo istate processing
}
