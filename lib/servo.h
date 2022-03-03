/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
typedef struct sho_t sho_t;
/**
   Struct for servo filtering
*/
typedef struct servo_t{
    dcell *mlead;       /**<lead filter temporary storage*/
    dcell *merrlast;    /**<recorded errro signal from last step*/
    dcell *mpreint;     /**<first integrator or other value.*/
    dccell *merrhist;   /**<Keep a short history of merr*/
    dccell *mint;       /**<second integrator. It is array to accomodate multiple ap's*/
    int initialized;    /**<is this data initialized*/
    int alint;          /**<Integral part of latency*/
    real alfrac;        /**<Fractional latency*/
    /*Servo parameters.*/
    dmat *ap;
    dmat *ep;
    real dt;
    /*if the corrector is a simple harmonic oscillator*/
    sho_t *sho;
}servo_t;
dcell *servo_optim(real dt, long dtrat, real al, real pmargin, real f0, real zeta, int servo_type, 
   const dmat *psdin,  const dmat* sigma2n);
real servo_optim_margin(real dt, long dtrat, real al, real pmargin, real f0, real zeta);
dmat *servo_rej2ol(const dmat *psdcl, real dt, long dtrat, real al, real gain, real sigma2n);
//cmat *servo_Hol(const dmat *nu, real dt, real dtrat, real al, const dmat *gain);
real servo_residual(real *noise_amp, const dmat *psdin, real dt, long dtrat, real al, const dmat *gain, int servo_type);
void servo_update(servo_t *st, const dmat *ep);
servo_t *servo_new(dcell *merr, const dmat *ap, real al, real dt, const dmat *ep);
servo_t *servo_new_sho(dcell *merr, const dmat *ap, real al, real dt, const dmat *ep, real f0, real zeta);
int servo_filter(servo_t *st, const dcell *merr);
void servo_add(servo_t *st, const dcell *madj, real alpha);
void servo_output(const servo_t *st, dcell **out);
dmat* servo_test(dmat *mideal, real dt, int dtrat, dmat* sigma2n, dmat *gain);
void servo_reset(servo_t *st);
void servo_free(servo_t *st);
//cmat *servo_typeII_Hol(const dmat *gain, real fs, real lgsdt);
/**
   struct for state space modeling of a second order harnomic oscillator with
   resonance frequency of f0 and damping ratio of zeta.
*/
struct sho_t{
    real dt;//minimum dt for advancing state.
    real c1;//2*zeta*omega0; omega0=2*pi*f0;
    real c2;//omega0^2;
    //real x1;//status (speed)
    //real x2;//status (position)
    dcell *dx; //status (speed)
    dcell *x; //status (position)
    dcell *xtmp; //temp
    dcell *ytmp; //temp
};
sho_t *sho_new(real f0, real zeta);
void sho_free(sho_t *sho);
//real sho_step(sho_t *sho, real xi, real dt);
//real sho_step_avg(sho_t *sho, real xi, real dt);
void sho_step(dcell **xout, sho_t *sho, dcell *xi, real dt);
void sho_step_avg(dcell **xout, sho_t *sho, dcell *xi, real dt);
void sho_reset(sho_t *sho);
dmat *sho_filter(const dmat *xi, real dt, real f0, real zeta);
#endif
