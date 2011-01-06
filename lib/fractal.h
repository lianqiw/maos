#ifndef AOS_LIB_FRACTAL_H
#define AOS_LIB_FRACTAL_H
#include "random.h"
#include "loc.h"
void fractal(double *p0, long nx, long ny, double dx, double r0, double L0);
void fractal_inv(double *p0, long nx, long ny, double dx, double r0, double L0);
void fractal_trans(double *p0, long nx, long ny, double dx, double r0, double L0);
void fractal_inv_trans(double *p0, long nx, long ny, double dx, double r0, double L0);
#endif
