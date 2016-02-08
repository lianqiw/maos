/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_LIB_KALMAN_H
#define AOS_LIB_KALMAN_H

#include "../math/mathdef.h"
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit);
typedef struct{
    dmat *Ad; /*discrete state propagation at dT*/
    dcell *Cd; /*From discrete state to WFS measurement*/
    dmat *AdM;/*discrete state propagation at dthi*/
    dmat *FdM;/*From discrete state to averaged mode for dthi*/
    dmat *Qn;
    dcell *M;  /*M is innovation gain.*/
    dmat *P;  /*Error covariance matrix*/
    double dthi;
    dmat *dtrat;
    dcell *Gwfs;
    dcell *Rwfs;
    dcell *Rn;
    dmat *xhat;
    dmat *xhat2;
    dmat *xhat3;
}kalman_t;
dmat* reccati(dmat **Pout, const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn);
dcell* reccati_cell(dmat **Pout, const dmat *A, const dmat *Qn, const dcell *C, const dcell *Rn);
kalman_t* sde_kalman(const dmat *coeff, double dthi, const dmat* dtrat, 
		     const dcell *Gwfs, const dcell *Rwfs, const dmat *Proj);
void kalman_free(kalman_t *kalman);
dmat *kalman_test(kalman_t *kalman, dmat *input);
void kalman_init(kalman_t *kalman);
void kalman_update(kalman_t *kalman, dmat *meas, int ik);
void kalman_output(kalman_t *kalman, dmat **out, double alpha, double beta);
#endif
