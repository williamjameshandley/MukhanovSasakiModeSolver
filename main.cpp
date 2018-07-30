#include <iostream>
#include <odepack/c-odepack.h>
#include <cmath>

void simple_pendulum_field(double *qdot,
			   const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;

  double alpha;
  alpha = parms[0];

  qdot[0] = q[1];
  qdot[1] = - alpha * sin(q[0]);

  return;
}

void simple_pendulum_jacobian(double *dfdq,
			      const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;

  double alpha;
  alpha = parms[0];
  
/*Column ordered; ODEPACK is in FORTRAN 77.*/
  dfdq[0] = 0.0; dfdq[2] = 1.0;
  
  dfdq[1] = - alpha * cos(q[0]); dfdq[3] = 0.0;

  return;
}

void simple_pendulum_cycle(double *g,
			   const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;
  
  double alpha;
  alpha = parms[0];

  double theta_0;
  theta_0 = parms[1];
/*-----------------------------------------------------------------------*/
  g[0] = (q[0] - theta_0) * q[1];

  /*Make sure the initial value is not detected!*/
  if(t < 1e-2)
    g[0] += 0.1;

  return;
}

int main(){
/*-----------------------------------------------------------------------*/  
  double opkd_rtol = 0.0, opkd_atol = 1e-12;
  dlsodar_problem *opkd;

  opkd = dlsodar_problem_create(2, 1,
				USR_FULL_JAC, 10000,
				&opkd_atol, &opkd_rtol);

  if(opkd == NULL)
    return 1;
  
  if(dlsodar_problem_init(opkd) != C_ODEPACK_SUCCESS)
    return 1;
/*-----------------------------------------------------------------------*/
  double q[2], t0, tf, dt, t;

  t0 = 0.0;
  tf = 100.0;
  dt = 0.1;

  q[0] = M_PI * 999/1000.0;
  q[1] = 0.;

  double alpha = 1.0;
  double parms[2] = {alpha, q[0]};
/*-----------------------------------------------------------------------*/  
  FILE *data;
  data = fopen("simple_pendulum_trajectory", "w");

  int i = 0, cyc = 0;
  t = t0;
  
  while(t < tf){    
    dlsodar_integrate(t + dt,
		      &t, q,
		      &simple_pendulum_field, NULL// &simple_pendulum_jacobian
		      , &simple_pendulum_cycle,
		      &parms, opkd);
    fprintf(data, "%.15lf\t%.15lf\t%.15lf\n", t + dt, q[0], q[1]);
    
    if(opkd->jroot[0] == 1)
      cyc += 1;
    if(cyc == 2)
      break;
  }

  fclose(data);
/*-----------------------------------------------------------------------*/  
  dlsodar_problem_close(opkd);

  return 0;
}
