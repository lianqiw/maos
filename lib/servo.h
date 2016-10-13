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

#ifndef SKYC_SERVO_H
#define SKYC_SERVO_H

#include "../math/mathdef.h"
#include "../math/mathdef.h"
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
    dccell *merrhist;   /**<Keep a short history of merr*/
    dccell *mint;       /**<second integrator. It is array to accomodate multiple ap's*/
    int initialized;    /**<is this data initialized*/
    int al;             /**<Additional latency*/
    /*Servo parameters.*/
    dmat *ap;
    dmat *ep;
    double dt;
}SERVO_T;
dcell* servo_optim(const dmat *psdin, double dt, long dtrat,  double pmargin, 
		   const dmat* sigma2n, int servo_type);
dmat *servo_rej2ol(const dmat *psdcl, double dt, long dtrat, double gain, double sigma2n);
cmat *servo_Hol(const dmat *nu, double dt, double dtrat, const dmat *gain);
double servo_residual(double *noise_amp, const dmat *psdin, double dt, long dtrat, const dmat *gain, int servo_type);
void servo_update(SERVO_T *st, const dmat *ep);
SERVO_T *servo_new(dcell *merr, const dmat *ap, int al, double dt, const dmat *ep);
int servo_filter(SERVO_T *st, const dcell *merr);
dmat* servo_test(dmat *mideal, double dtngs, int dtrat, dmat* sigma2n, dmat *gain);
void servo_reset(SERVO_T *st);
void servo_free(SERVO_T *st);
cmat *servo_typeII_Hol(const dmat *gain, double fs, double lgsdt);
/**
   struct for state space modeling of a second order harnomic oscillator with
   resonance frequency of f0 and damping ratio of zeta.
*/
typedef struct SHO_T{
    double dt;//minimum dt for advancing state.
    double c1;//2*zeta*omega0; omega0=2*pi*f0;
    double c2;//omega0^2;
    double x1;//status (derivative)
    double x2;//status (position)
}SHO_T;
SHO_T *sho_new(double f0, double zeta);
double sho_step(SHO_T *sho, double xi, double dt);
void sho_reset(SHO_T *sho);
dmat *sho_filter(const dmat *xi, double dt, double f0, double zeta);
#endif
