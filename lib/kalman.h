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
real sde_fit(dmat **pcoeff, const dmat *psdin, real tmax_fit, int vibid);
real sde_fit_auto(dmat **pcoeff, const_anyarray psdin, real tfit);
void sde_psd(dmat **psd, const dmat *f, const real *coeff, int ncoeff, int nmod);
dmat *sde_psd2(const dmat *ff, const dmat *coeff);
typedef struct{
    dmat *AdM; /**<discrete state propagation at dthi*/
    dmat *BM; /**<From discrete state to averaged mode for dthi*/
	dmat *Qn;  /**<Covariance of the process noise term e_k*/

    lmat *dtrat_wfs;/**<WFS sampling period over dthi*/
	lmat *dtrats;/**<uniq dtrats */
	lmat *mdirect;/**<Modes that are directly controlled by the slow loop */

	//The following are per dtrat
	dcell *Ad;  /**<discrete state propagation at dT*/
	dcell *Cd; /**<From discrete state to WFS measurement*/
    dcell *Gwfs;/**<WFS measurement from modes. Can be identity.*/
    dcell *Cnn;/**<WFS measurement noise covariance due to photon and RoN.*/
    dcell *Rn;  /**<Total WFS measurement error due to signal evolution and Rwfs. */
	dcell *M;  /**<M is innovation gain.*/
    dcell *P;  /**<Error covariance matrix*/
	dccell *Rlsq; /**<least squares reconstructor */

    real dthi;/**<Sampling period of control loop*/

    dcell *xhat; /**<stores x_k+1|k*/
    dcell *xhat2;/**<stores x_k+(1+j)/M|k*/
    dmat  *xhat3;/**<temporary*/
	dcell *xout;/**<output */
	dcell *psol;/**<current correction for psol */
}kalman_t;
real reccati(dmat **Kinf, dmat **Pout, const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn);
//dcell* reccati_mr(dmat **Pout, const dmat *A, const dmat *Qn, const dcell *C, const dcell *Rn);
kalman_t* sde_kalman(const dmat *coeff, const real dthi, const lmat* dtrat_wfs, const lmat *mdirect, 
	const dcell *Gwfs, const dcell *Rwfs, const dmat *Proj);
void kalman_free(kalman_t *kalman);
dmat *kalman_test(kalman_t *kalman, const dcell *Gwfs2, const dmat *input, int flag);
dmat* kalman_test2(const dmat *coeff, const real dthi, const lmat* dtrat, const lmat* mdirect,
	const dcell *Gwfs, const dcell *Rwfs, const dmat *Proj, const dcell *Gwfs2, const dmat *input, int flag);
dmat* kalman_test3(const dmat* input, int flag, const char* fn);
void kalman_init(kalman_t *kalman);
void kalman_update(kalman_t *kalman, dcell *meas, int idtrat);
void kalman_output(kalman_t *kalman, dmat **out, real alpha, real beta, int idtrat);
void kalman_write(kalman_t *kalman, const char *format, ...) CHECK_ARG(2);
kalman_t* kalman_read(const char* format, ...)  CHECK_ARG(1);
#endif
