#pragma once

struct dlsodar_problem;

typedef void (*odepack_field_func)(double *qdot,
				   const double t, const double *q, void *data);

typedef void (*odepack_jacobian_func)(double *dfdq,
				      const double t, const double *q, void *data);

typedef void (*odepack_root_func)(double *g,
				  const double t, const double *q, void *data);

extern "C" dlsodar_problem* dlsodar_problem_create (int neq, int ng, int jac_type, int max_steps, double *atol, double *rtol);
extern "C" int dlsodar_problem_init(dlsodar_problem *dls);
extern "C" void dlsodar_problem_close(dlsodar_problem *dls);

extern "C" int dlsodar_integrate
(double t,
 double *t0, double *q,
 odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
 void *data, dlsodar_problem *dls);
