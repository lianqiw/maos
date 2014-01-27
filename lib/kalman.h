#ifndef AOS_LIB_KALMAN_H
#define AOS_LIB_KALMAN_H

#include "dmat.h"
#include "cmat.h"
#include "fft.h"
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, double min, double max, double df);
typedef struct{
    dmat *Ad; /*discrete state propagation at dT*/
    dcell *Cd; /*From discrete state to WFS measurement*/
    dmat *AdM;/*discrete state propagation at dthi*/
    dmat *FdM;/*From discrete state to averaged mode for dthi*/
    dmat *Qn;
    dcell *M;  /*M is innovation gain.*/
    dcell *P;  /*Error covariance matrix*/
    double dthi;
    dmat *dtrat;
    dcell *Gwfs;
    dcell *Rwfs;
    dcell *Rn;
    dmat *xhat;
    dmat *xhat2;
    dmat *xhat3;
}kalman_t;
void reccati(dmat **Mout, dmat **Pout, 
	     const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn);
kalman_t* sde_kalman(dmat *coeff, double dthi, dmat* dtrat, dcell *Gwfs, dcell *Rwfs, dmat *Proj);
void kalman_free(kalman_t *kalman);
dmat *kalman_test(kalman_t *kalman, dmat *input);
void kalman_init(kalman_t *kalman);
void kalman_update(kalman_t *kalman, dmat *meas, int ik);
void kalman_output(kalman_t *kalman, dmat **out, double alpha, double beta);
#endif
