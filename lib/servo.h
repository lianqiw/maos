/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef SKYC_SERVO_H
#define SKYC_SERVO_H

#include "dmat.h"
#include "cmat.h"
#include "fft.h"
#include "loc.h"
/**
   \file servo.h
   Routines for servo optimization, filtering, etc.
*/

/**
   Struct for servo filtering
*/
typedef struct SERVO_T{
    dcell *mlead;       /**<lead filter temporary storage*/
    dcell *merrlast;    /**<recorded errro signal from last step*/
    dcell *mpreint;     /**<first integrator or other value.*/
    dcell **mint;       /**<second integrator. It is array to accomodate multiple ap's*/
    int nmint;          /**<number of cells in mint*/
    int initialized;   /**<is this data initialized*/
}SERVO_T;
dcell* servo_optim(const dmat *psdin, double dt, long dtrat,  double pmargin, 
		   const dmat* sigman, int servo_type);
cmat *servo_Hol(const dmat *nu, double dt, double dtrat, const dmat *gain);
double servo_residual(double *noise_amp, const dmat *psdin, double dt, long dtrat, const dmat *gain, int servo_type);
SERVO_T *servo_new(dcell *merr, const dmat *gain);
void servo_filter(SERVO_T *st, dcell *merr, double dtngs, const dmat *gain);
void servo_shift(SERVO_T *st, dmat *ap);
dmat *psd2temp(dmat *psdin, double dt, double N, rand_t* rstat);
dmat* servo_test(dmat *mideal, double dtngs, int dtrat, double sigma2n, dmat *gain);
void servo_free(SERVO_T *st);
cmat *servo_typeII_Hol(const dmat *gain, double fs, double lgsdt);
double psd_inte(double *nu, double *psd, long n);
double psd_inte2(dmat *psdin);
dmat* psd2time(dmat *psdin, rand_t *rstat, double dt, int nstep);
#endif
