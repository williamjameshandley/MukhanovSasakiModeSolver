/*
  Copyright (C) 2012  Akshay Srinivasan 
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

#ifndef __C_ODEPACK_DLSODAR__
#define __C_ODEPACK_DLSODAR__

/* Quoting from opkdmain.f
C ISTATE = 2 or 3  if DLSODAR was successful, negative otherwise.
C           2 means no root was found, and TOUT was reached as desired.
C           3 means a root was found prior to reaching TOUT.
C          -1 means excess work done on this call (perhaps wrong JT).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of JT or tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means work space insufficient to finish (see messages).
*/

#define DLSODAR_NO_ROOT_FOUND 2
#define DLSODAR_ROOT_FOUND 3
#define DLSODAR_EXCESS_WORK_DONE -1
#define DLSODAR_EXCESS_ACCURACY_REQUESTED -2
#define DLSODAR_ILLEGAL_INPUT -3
#define DLSODAR_REPEATED_FAILURE -4
#define DLSODAR_REPEATED_CONVERGENCE_FAILURE -5
#define DLSODAR_ERROR_IS_ZERO -6
#define DLSODAR_INSUFFICIENT_WORKSPACE -7

/*---------------JT----------------------*/
#define USR_FULL_JAC 1
#define INT_FULL_JAC 2

#define USR_BAND_JAC 4
#define INT_BAND_JAC 5
/*-----------------------------------------------------------------------*/

typedef struct{
  int neq, ng;
  
  int jac_type;
  /* i - lower_bandwidth <= j <= i + upper_bandwidth */
  int jac_lower_bandwidth, jac_upper_bandwidth;

  int itol;
  double *rtol, *atol;
/*Quoting from opkdmain.f
C              ITOL    RTOL      ATOL      EWT(i)
C              ----    ------    ------    -----------------------------
C              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
C              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
C              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
C              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
*/   

  int itask;
  double t_critical;
/*Quoting from opkdmain.f
C     ITASK    An index specifying the task to be performed.  Input
C              only.  ITASK has the following values and meanings:
C              1   Normal computation of output values of y(t) at
C                  t = TOUT (by overshooting and interpolating).
C              2   Take one step only and return.
C              3   Stop at the first internal mesh point at or beyond
C                  t = TOUT and return.
C              4   Normal computation of output values of y(t) at
C                  t = TOUT but without overshooting t = TCRIT.  TCRIT
C                  must be input as RWORK(1).  TCRIT may be equal to or
C                  beyond TOUT, but not behind it in the direction of
C                  integration.  This option is useful if the problem
C                  has a singularity at or beyond t = TCRIT.
C              5   Take one step, without passing TCRIT, and return.
C                  TCRIT must be input as RWORK(1).
C
C              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C              (within roundoff), it will return T = TCRIT (exactly) to
C              indicate this (unless ITASK = 4 and TOUT comes before
C              TCRIT, in which case answers at T = TOUT are returned
C              first).
*/

  double init_step_size, max_step_size, min_step_size;
  int max_steps,\
    max_order_stiff, max_order_non_stiff,\
    print_at_switch, max_num_messages; 
/*
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C IXPR    IWORK(5)  flag to generate extra printing at method switches.
C                   IXPR = 0 means no extra printing (the default).
C                   IXPR = 1 means print data on each switch.
C                   T, H, and NST will be printed on the same logical
C                   unit as used for error messages.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
C                   (Adams) method.  The default value is 12.
C                   If MXORDN exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDN is held constant during the problem.
C
C MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
C                   (BDF) method.  The default value is 5.
C                   If MXORDS exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDS is held constant during the problem.
  */
  /*-----------------------------------------------------------------------
    Do not tweak these variables, if you don't know what you are doing.
    These are used in the actual dlsodar call
    -----------------------------------------------------------------------*/
  int iopt, istate, *iwork, liw, lrw, *jroot, jt;
  double *rwork;
} dlsodar_problem;

dlsodar_problem* dlsodar_problem_create
(int neq, int ng,
 int jac_type, int max_steps,
 double *atol, double *rtol);

int dlsodar_problem_init(dlsodar_problem *dls);
void dlsodar_problem_close(dlsodar_problem *dls);

int dlsodar_integrate
(double t,
 double *t0, double *q,
 odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
 void *data, dlsodar_problem *dls);
#endif
