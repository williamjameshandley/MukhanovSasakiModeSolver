#include "opkdmain.h"
#include "dlsodar.hpp"
#include "pointers.hpp"

dlsodar::dlsodar(int neq_, int ng_): 
        neq{neq_},
        itol{1},
        rtol{1e-8},
        atol{1e-8},
        itask{1},
        istate{1},
        iopt{0},
        rwork(static_cast<size_t>(22+neq*std::max(16,neq+9)+3*ng_)),
        lrw{static_cast<int>(rwork.size())},
        iwork(static_cast<size_t>(20+neq)),
        liw{static_cast<int>(iwork.size())},
        jt{2},
        ng{ng_},
        jroot(static_cast<size_t>(ng))
{}

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, void *data)
{ jt = 2; int ng_ = ng; ng = 0; _integrate(t, tout, q, f_func, nullptr, nullptr, data); ng = ng_; }

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, Root g_func, void *data)
{ jt = 2; _integrate(t, tout, q, f_func, nullptr, g_func, data); }

void dlsodar::integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data)
{ jt = 1; _integrate(t, tout, q, f_func, j_func, g_func, data); }

void dlsodar::_integrate(double &t, double tout, double q[], Field f_func, Jacobian j_func, Root g_func, void *data){

    auto f = [&](const int *, const double *t_, const double *y, double *ydot)                                        -> void {f_func(ydot, *t_, y, data);};
    auto j = [&](const int *, const double *t_, const double *y, const int *, const int *, double *dfdy, const int *) -> void {j_func(dfdy, *t_, y, data);};
    auto g = [&](const int *, const double *t_, const double *y, const int *, double *gout)                           -> void {g_func(gout, *t_, y, data);};

    dlsodar_(
            function_ptr(f),
            &neq, q, &t, &tout,
            &itol, &rtol[0], &atol[0],
            &itask, &istate, &iopt,
            &rwork[0], &lrw, &iwork[0], &liw,
            function_ptr(j), &jt,
            function_ptr(g), &ng,
            &jroot[0]
            );
    // @todo istate processing
}
