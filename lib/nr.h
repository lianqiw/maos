#ifndef AOS_LIB_NR_H
#define AOS_LIB_NR_H
/*
  Contains routines borrowed from numerical recipes
 */
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
void amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(double [], void *data), void *data, int *nfunk);
#endif
