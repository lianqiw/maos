#ifndef AOS_LIB_FRACTAL_H
#define AOS_LIB_FRACTAL_H
#include "random.h"
#include "loc.h"
#define ARGS double *p0, long nx, long ny, double dx, double r0, double L0, long ninit
void fractal(ARGS);
void fractal_inv(ARGS);
void fractal_trans(ARGS);
void fractal_inv_trans(ARGS);
void fractal_vkcov_free(void);
#endif
