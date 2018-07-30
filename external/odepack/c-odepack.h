
#ifndef __C_ODEPACK__
#define __C_ODEPACK__

#define C_ODEPACK_SUCCESS 0
#define C_ODEPACK_MEM_ERROR -20
#define C_ODEPACK_UNKNOWN_OPTION -40

typedef void (*odepack_field_func)(double *qdot,
				   const double t, const double *q, void *data);

typedef void (*odepack_jacobian_func)(double *dfdq,
				      const double t, const double *q, void *data);

typedef void (*odepack_root_func)(double *g,
				  const double t, const double *q, void *data);

#include "dlsode.h"
#include "dlsodar.h"

#endif
