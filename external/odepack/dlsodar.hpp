#pragma once
#include <vector>
#include <functional>

using odepack_field_func = std::function<void(double *qdot, const double t, const double *q, void *data)>;
using odepack_jacobian_func = std::function<void(double *dfdq, const double t, const double *q, void *data)>;
using odepack_root_func = std::function<void(double *g, const double t, const double *q, void *data)>;

class dlsodar
{
    private:
        int neq, ng;

        int jac_type, jac_lower_bandwidth, jac_upper_bandwidth;

        int itol;

        std::vector<double> rtol, atol;

        int itask;
        double t_critical;

        double init_step_size, max_step_size, min_step_size;
        int max_steps, max_order_stiff, max_order_non_stiff, print_at_switch, max_num_messages; 

        int iopt, istate, liw, lrw, jt;
        std::vector<int> iwork, jroot;
        std::vector<double> rwork;

    public:
        dlsodar(int, int, int, int, std::vector<double>, std::vector<double>); 
        void integrate(double t, double *t0, double *q, odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func, void *data);

};

extern "C" void dlsodar_(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
		const int *neq, double *y, double *t, const double *tout,
		const int *itol, const double *rtol, const double *atol,
		const int *itask, int *istate, const int *iopt,
		double *rwork, const int *lrw, int *iwork, const int *liw,
		void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd), const int *jt,
		void (*g)(const int *neq, const double *t, const double *y, const int *ng, double *gout), const int *ng,
		int *jroot);
