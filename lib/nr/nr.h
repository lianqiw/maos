#ifndef AOS_LIB_NR_H
#define AOS_LIB_NR_H
#include "../common.h"
#define nrerror error
double chebev(double a, double b, double c[], int m, double x);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
#endif
