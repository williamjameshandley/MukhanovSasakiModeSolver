/*
  Copyright (C) 2011  Akshay Srinivasan 
  <akshay@ncbs.res.in>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "opkdmain.h"
#include "c-odepack.h"
#include "utility.h"

dlsode_problem* dlsode_problem_create
(int neq,
 int step_method, int iter_method,
 int max_steps, double *atol, double *rtol)
{
  dlsode_problem *ret;
  ret = calloc(1, sizeof(dlsode_problem));
  if(ret == NULL)
    goto mem_error;
  
  ret->neq = neq;
  ret->step_method = step_method;
  ret->iter_method = iter_method;
  ret->max_steps = max_steps;
  ret->atol = atol;
  ret->rtol = rtol;

  return ret;

mem_error:
  fprintf(stderr, "dlsode_problem_create: Cannot allocate memory.\n");
  return NULL;
}

int dlsode_problem_init(dlsode_problem *dls)
{
  int i, j, k;
  /*-----------------------------------------------------------------------*/
  dls->mf = dls->step_method * 10 + dls->iter_method;
  /*-----------------------------------------------------------------------*/
  int max_ord;

  if(dls->step_method == ADAMS_IMPLICIT){
    if(dls->max_order == 0)
      max_ord = 12;
    else
      max_ord = MIN(dls->max_order, 12);
  }
  else if(dls->step_method == BDF){
    if(dls->max_order == 0)
      max_ord = 5;
    else
      max_ord = MIN(dls->max_order, 5);
  }
  else
    goto unknown_step_error;
  /*-----------------------------------------------------------------------*/
  int lwm;
  
  if(dls->iter_method == 0){
    lwm = 0;
    dls->liw = 20;
  }
  else if((dls->iter_method == 1) || (dls->iter_method == 2)){
    lwm = dls->neq * dls->neq + 2;
    dls->liw = 20 + dls->neq;
  }
  else if(dls->iter_method == 3){
    lwm = dls->neq + 2;
    dls->liw = 20;
  }
  else if((dls->iter_method == 4) || (dls->iter_method == 5)){
    lwm = (2 * dls->jac_lower_bandwidth + dls->jac_upper_bandwidth + 1) * dls->neq + 2;
    dls->liw = 20 + dls->neq;
  }
  else
    goto unknown_iter_error;

  dls->lrw = 20 + dls->neq * (max_ord + 1) + 3 * dls->neq + lwm;
  /*--------------------------------------------------------------*/
  dls->iopt = 1;

  dls->rwork = calloc(dls->lrw, sizeof(double));
  if(dls->rwork == NULL)
    goto mem_error;
  
  dls->iwork = calloc(dls->liw, sizeof(int));
  if(dls->iwork == NULL)
    goto mem_error;

  dls->rwork[4] = dls->init_step_size;
  dls->rwork[5] = dls->max_step_size;
  dls->rwork[6] = dls->min_step_size;
  
  dls->iwork[4] = dls->max_order;
  dls->iwork[5] = dls->max_steps;
  /*--------------------------------------------------------------*/
  if((dls->itol < 2) || (dls->itol > 4)){
    /*Error bounds are same for all dimensions.*/
    dls->itol = 1;
  }
  /*--------------------------------------------------------------*/
  if((dls->itask < 2) || (dls->itask > 5))
    dls->itask = 1;
  else{
    dls->itask = dls->itask;
    if(dls->itask > 3)
      dls->rwork[0] = dls->t_critical;
  }
  /*--------------------------------------------------------------*/
  dls->istate = 1;
  /*--------------------------------------------------------------*/
  if((dls->iter_method == 4) || (dls->iter_method == 5)){
    dls->iwork[0] = dls->jac_lower_bandwidth;
    dls->iwork[1] = dls->jac_upper_bandwidth;
  }
  /*--------------------------------------------------------------*/
  
  return C_ODEPACK_SUCCESS;

mem_error:
  fprintf(stderr, "dlsode_problem_init: Cannot allocate memory.\n");
  return C_ODEPACK_MEM_ERROR;
unknown_step_error:
  fprintf(stderr, "dlsode_problem_init: Unknown step method. METH = %d\n", dls->step_method);
  return C_ODEPACK_UNKNOWN_OPTION;
unknown_iter_error:
  fprintf(stderr, "dlsode_problem_init: Unknown iteration method. MITER = %d\n", dls->iter_method);
  return C_ODEPACK_UNKNOWN_OPTION;
}

void dlsode_problem_close(dlsode_problem *dls){
  if(dls == NULL)
    return;
  
  free(dls->rwork);
  free(dls->iwork);
  free(dls);
}

int dlsode_integrate
(double t,
 double *t0, double *q,
 odepack_field_func f_func, odepack_jacobian_func j_func,
 void *data, dlsode_problem *dls)
{
  /*Not ANSI, don't know if this is thread safe.*/
  void dlsode_field_compat(const int *neq, const double *t_, const double *y, double *ydot){
    f_func(ydot, *t_, y, data);
    return;
  }

  void dlsode_jacobian_compat(const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd){
    j_func(pd, *t_, q, data);
    return;
  }

  dlsode_(dlsode_field_compat,
	  &dls->neq, q, t0, &t,
	  &dls->itol, dls->rtol, dls->atol,
	  &dls->itask, &dls->istate, &dls->iopt,
	  dls->rwork, &dls->lrw, dls->iwork, &dls->liw,
	  dlsode_jacobian_compat,
	  &dls->mf);

  return dls->istate;
}
