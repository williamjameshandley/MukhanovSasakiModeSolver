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

dlsodar_problem* dlsodar_problem_create
(int neq, int ng,
 int jac_type, int max_steps,
 double *atol, double *rtol){
  dlsodar_problem *dls;
  dls = calloc(1, sizeof(dlsodar_problem));
  if(dls == NULL)
    goto mem_error;

  dls->neq = neq;
  dls->ng = ng;
  dls->max_steps = max_steps;
  dls->jac_type = jac_type;
  dls->atol = atol;
  dls->rtol = rtol;

  return dls;

mem_error:
  fprintf(stderr, "dlsodar_problem_create: Cannot allocate memory.\n");
  return NULL;
}

int dlsodar_problem_init(dlsodar_problem *dls){
  int i, j, k;
  /*-----------------------------------------------------------------------*/
  dls->jroot = calloc(dls->ng, sizeof(int));
  if((dls->jroot == NULL) && (dls->ng > 0))
    goto mem_error;
  /*-----------------------------------------------------------------------*/
  dls->lrw = 22 + dls->neq * MAX(16, dls->neq + 9) + 3 * dls->ng;
  dls->rwork = calloc(dls->lrw, sizeof(double));
  if(dls->rwork == NULL)
    goto mem_error;
  /*-----------------------------------------------------------------------*/
  dls->liw = 20 + dls->neq;
  dls->iwork = calloc(dls->liw, sizeof(int));
  if(dls->iwork == NULL)
    goto mem_error;
  /*-----------------------------------------------------------------------*/
  dls->iopt = 1;

  dls->rwork[4] = dls->init_step_size;
  dls->rwork[5] = dls->max_step_size;
  dls->rwork[6] = dls->min_step_size;

  dls->iwork[4] = dls->print_at_switch;
  dls->iwork[5] = dls->max_steps;
  dls->iwork[6] = dls->max_num_messages;
  dls->iwork[7] = dls->max_order_non_stiff;
  dls->iwork[8] = dls->max_order_stiff;
/*--------------------------------------------------------------*/
  if((dls->itol < 2) || (dls->itol > 4)){
    /*Error bounds are same for all dimensions.*/
    dls->itol = 1;
  }

/*--------------------------------------------------------------*/
  if((dls->itask < 2) || (dls->itask > 5))
    dls->itask = 1;

  if(dls->itask > 3)
    dls->rwork[0] = dls->t_critical;
/*--------------------------------------------------------------*/
  dls->istate = 1;
/*--------------------------------------------------------------*/
  dls->jt = dls->jac_type;
  if((dls->jt == 4) || (dls->jt == 5)){
    dls->iwork[0] = dls->jac_lower_bandwidth;
    dls->iwork[1] = dls->jac_upper_bandwidth;
  }
/*--------------------------------------------------------------*/  
  return C_ODEPACK_SUCCESS;

mem_error:
  fprintf(stderr, "dlsodar_problem_init: Cannot allocate memory.\n");
  return C_ODEPACK_MEM_ERROR;
}


void dlsodar_problem_close(dlsodar_problem *dls){
  if(dls == NULL)
    return;
  
  free(dls->rwork);
  free(dls->iwork);
  free(dls->jroot);
  free(dls);  
}

int dlsodar_integrate(double t,
		      double *t0, double *q,
		      odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
		      void *data, dlsodar_problem *dls){

  /*Not ANSI, probably not thread safe, unless C has
    closures or something like that :)
  */
  void dlsodar_field_compat(const int *neq, const double *t_, const double *y, double *ydot){
    f_func(ydot, *t_, y, data);
    return;
  }
  
  void dlsodar_jacobian_compat(const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd){
    j_func(pd, *t_, y, data);
    return;
  }

  void dlsodar_constraint_compat(const int *neq, const double *t_, const double *y, const int *ng, double *gout){
    c_func(gout, *t_, y, data);
  }
  
  dlsodar_(dlsodar_field_compat,
	   &dls->neq, q, t0, &t,
	   &dls->itol, dls->rtol, dls->atol,
	   &dls->itask, &dls->istate, &dls->iopt,
	   dls->rwork, &dls->lrw, dls->iwork, &dls->liw,
	   dlsodar_jacobian_compat, &dls->jt,
	   dlsodar_constraint_compat, &dls->ng,
	   dls->jroot);

  return dls->istate;
}
