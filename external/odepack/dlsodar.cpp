#include "dlsodar.hpp"
#include "pointers.hpp"

dlsodar::dlsodar(int neq_, int ng_, int jac_type_, int max_steps_, std::vector<double> atol_, std::vector<double> rtol_): 
    neq{neq_}, ng{ng_}, jac_type{jac_type_}, max_steps{max_steps_}, atol{atol_}, rtol{rtol_},
    jroot(ng), 
    lrw(22+neq*std::max(16,neq+9)+3*ng),
    rwork(lrw),
    liw(20+neq),
    iwork(liw),
    iopt(1)
{
    rwork[4] = init_step_size;
    rwork[5] = max_step_size;
    rwork[6] = min_step_size;

    iwork[4] = print_at_switch;
    iwork[5] = max_steps;
    iwork[6] = max_num_messages;
    iwork[7] = max_order_non_stiff;
    iwork[8] = max_order_stiff;

    if((itol < 2) || (itol > 4)) itol = 1;
    if((itask < 2) || (itask > 5)) itask = 1;
    if(itask > 3) rwork[0] = t_critical;

    istate = 1;
    jt = jac_type;
    if((jt == 4) || (jt == 5)){
        iwork[0] = jac_lower_bandwidth;
        iwork[1] = jac_upper_bandwidth;
    }

}

void dlsodar::integrate(double t,
        double *t0, double *q,
        odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
        void *data){

    auto dlsodar_field_compat = [&](const int *neq, const double *t_, const double *y, double *ydot) -> void { f_func(ydot, *t_, y, data);};
    auto dlsodar_jacobian_compat = [&](const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd) -> void {j_func(pd, *t_, y, data); };
    auto dlsodar_constraint_compat = [&] (const int *neq, const double *t_, const double *y, const int *ng, double *gout) -> void{ c_func(gout, *t_, y, data); };

    dlsodar_(
            function_ptr(dlsodar_field_compat),
            &neq, q, t0, &t,
            &itol, &rtol[0], &atol[0],
            &itask, &istate, &iopt,
            &rwork[0], &lrw, &iwork[0], &liw,
            function_ptr(dlsodar_jacobian_compat), &jt,
            function_ptr(dlsodar_constraint_compat), &ng,
            &jroot[0]
            );
}
