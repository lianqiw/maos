#ifndef AOS_LIB_KALMAN_H
#define AOS_LIB_KALMAN_H

#include "dmat.h"
#include "cmat.h"
#include "fft.h"
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, double min, double max);
typedef struct{
    dmat *Ad; /*discrete state propagation at dT*/
    dmat *Fd; /*From discrete state to averaged mode dT*/
    dmat *Cd; /*From discrete state to WFS measurement*/
    dmat *AdM;/*discrete state propagation at dthi*/
    dmat *FdM;/*From discrete state to averaged mode for dthi*/
    dmat *M;  /*M is innovation gain.*/
    dmat *P;  /*Error covariance matrix*/
    double dthi;
    double dtrat;
    dmat *Gwfs;
    dmat *Rwfs;
    dmat *xhat;
    dmat *xhat2;
    dmat *xhat3;
}kalman_t;
kalman_t* sde_kalman(dmat *coeff, double dthi, int dtrat, dmat *Gwfs, dmat *Rwfs, dmat *Proj);
void kalman_free(kalman_t *kalman);
dmat *kalman_test(kalman_t *kalman, dmat *input);
void kalman_init(kalman_t *kalman);
void kalman_update(kalman_t *kalman, dmat *meas);
void kalman_output(kalman_t *kalman, dmat **out, double alpha, double beta);
#endif
