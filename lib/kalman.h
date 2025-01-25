/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file kalman.h
*/
#ifndef AOS_LIB_KALMAN_H
#define AOS_LIB_KALMAN_H

#include "../math/mathdef.h"
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, real tmax_fit, int vibid);
typedef struct{
    dmat *Ad;  /**<discrete state propagation at dT*/
    dcell *Cd; /**<From discrete state to WFS measurement*/
    dmat *AdM; /**<discrete state propagation at dthi*/
    dmat *FdM; /**<From discrete state to averaged mode for dthi*/
    dmat *Qn;  /**<Covariance of the process noise term e_k*/
    dcell *M;  /**<M is innovation gain.*/
    dmat *P;   /**<Error covariance matrix*/
    real dthi;/**<Sampling period of control loop*/
    lmat *dtrat;/**<WFS sampling period over dthi*/
    dcell *Gwfs;/**<WFS measurement from modes. Can be identity.*/
    dcell *Rwfs;/**<WFS measurement noise covariance due to photon and RoN.*/
    dcell *Rn;  /**<Total WFS measurement error due to signal evolution and Rwfs. */
    dmat *xhat; /**<Temporary state*/
    dmat *xhat2;/**<Temporary state*/
    dmat *xhat3;/**<Temporary state*/
}kalman_t;
dmat* reccati(dmat **Pout, const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn);
dcell* reccati_cell(dmat **Pout, const dmat *A, const dmat *Qn, const dcell *C, const dcell *Rn);
kalman_t* sde_kalman(const dmat *coeff, const real dthi, const lmat* dtrat, 
		     const dcell *Gwfs, const dcell *Rwfs, const dmat *Proj);
void kalman_free(kalman_t *kalman);
dmat *kalman_test(kalman_t *kalman, dmat *input);
void kalman_init(kalman_t *kalman);
void kalman_update(kalman_t *kalman, dmat *meas, int ik);
void kalman_output(kalman_t *kalman, dmat **out, real alpha, real beta);
void kalman_write(kalman_t *kalman, const char *format, ...) CHECK_ARG(2);
void sde_psd(dmat **psd, const dmat *f, const real *coeff, int ncoeff, int nmod);
dmat *sde_psd2(const dmat *ff, const dmat *coeff);
#endif
