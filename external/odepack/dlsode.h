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

#ifndef __C_ODEPACK_DLSODE__
#define __C_ODEPACK_DLSODE__
/*
  Quoting from opkdmain.f
C     ISTATE:INOUT  Index used for input and output to specify the state
C                   of the calculation.
C                   Input:
C                    1   This is the first call for a problem.
C                    2   This is a subsequent call.
C                   Output:
C                    1   Nothing was done, because TOUT was equal to T.
C                    2   DLSODE was successful (otherwise, negative).
C                        Note that ISTATE need not be modified after a
C                        successful return.
C                   -1   Excess work done on this call (perhaps wrong
C                        MF).
C                   -2   Excess accuracy requested (tolerances too
C                        small).
C                   -3   Illegal input detected (see printed message).
C                   -4   Repeated error test failures (check all
C                        inputs).
C                   -5   Repeated convergence failures (perhaps bad
C                        Jacobian supplied or wrong choice of MF or
C                        tolerances).
C                   -6   Error weight became zero during problem
C                        (solution component i vanished, and ATOL or
C                        ATOL(i) = 0.).
*/
#define DLSODE_FIRST_CALL 1
#define DLSODE_SUBSEQ_CALL 2
#define DLSODE_EXCESS_WORK_DONE -1
#define DLSODE_EXCESS_ACCURACY_REQUESTED -2
#define DLSODE_ILLEGAL_INPUT -3
#define DLSODE_REPEATED_FAILURE -4
#define DLSODE_REPEATED_CONVERGENCE_FAILURE -5
#define DLSODE_ERROR_IS_ZERO -6

/* Quoting from opkdmain.f
C              METH indicates the basic linear multistep method:
C              1   Implicit Adams method.
C              2   Method based on backward differentiation formulas
C                  (BDF's).
C
C              MITER indicates the corrector iteration method:
C              0   Functional iteration (no Jacobian matrix is
C                  involved).
C              1   Chord iteration with a user-supplied full (NEQ by
C                  NEQ) Jacobian.
C              2   Chord iteration with an internally generated
C                  (difference quotient) full Jacobian (using NEQ
C                  extra calls to F per df/dy value).
C              3   Chord iteration with an internally generated
C                  diagonal Jacobian approximation (using one extra call
C                  to F per df/dy evaluation).
C              4   Chord iteration with a user-supplied banded Jacobian.
C              5   Chord iteration with an internally generated banded
C                  Jacobian (using ML + MU + 1 extra calls to F per
C                  df/dy evaluation).
C
C              If MITER = 1 or 4, the user must supply a subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For other values of MITER, a dummy argument can be used.
*/

/*--------------METH-----------------*/
#define ADAMS_IMPLICIT 1
#define BDF 2

/*-------------MITER-----------------*/

#define FUNC_ITER 0

#define CHORD_ITER_USR_FULL_JAC 1
#define CHORD_ITER_INT_FULL_JAC 2

#define CHORD_ITER_INT_DIAG_JAC 3

#define CHORD_ITER_USR_BAND_JAC 4
#define CHORD_ITER_INT_BAND_JAC 5


typedef struct{
  int neq;
  int step_method, iter_method;
  
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
  
  int jac_lower_bandwidth, jac_upper_bandwidth;
  /*
     i - lower_bandwidth <= j <= i + upper_bandwidth
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
  int max_steps, max_order;
/*Quoting from opkdmain.f
  C     Name    Location   Meaning and default value
C     ------  ---------  -----------------------------------------------
C     H0      RWORK(5)   Step size to be attempted on the first step.
C                        The default value is determined by the solver.
C     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
C                        default value is infinite.
C     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
C                        default value is 0.  (This lower bound is not
C                        enforced on the final step before reaching
C                        TCRIT when ITASK = 4 or 5.)
C     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
C                        is 12 if METH = 1, and 5 if METH = 2. (See the
C                        MF description above for METH.)  If MAXORD
C                        exceeds the default value, it will be reduced
C                        to the default value.  If MAXORD is changed
C                        during the problem, it may cause the current
C                        order to be reduced.
C     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
C                        allowed during one call to the solver.  The
C                        default value is 500.
*/
  /*-----------------------------------------------------------------------
    Do not tweak these variables, if you don't know what you are doing.
    These are used in the actual dlsode call
    -----------------------------------------------------------------------*/
  int iopt, istate, *iwork, liw, lrw, mf;
  double *rwork;
} dlsode_problem;

dlsode_problem* dlsode_problem_create
(int neq,
 int step_method, int iter_method,
 int max_steps, double *atol, double *rtol);

int dlsode_problem_init(dlsode_problem *sess);
void dlsode_problem_close(dlsode_problem *sess);

int dlsode_integrate
(double t,
 double *t0, double *q,
 odepack_field_func func, odepack_jacobian_func jac_func,
 void *data, dlsode_problem *sess);
#endif
