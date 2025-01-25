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
	union{
		cell *mint;       /**<second integrator. It is array to accomodate multiple ap's*/
		dccell *mintc;
		dcell  *mintd;
	};
	union{
		cell *mpreint;     /**<first integrator or other value.*/
		dcell *mpreintc;
		dmat  *mpreintd;
	};
	cell *mlead;       /**<lead filter temporary storage*/
	cell *merrlast;    /**<recorded errro signal from last step*/
	cell *merrhist;   /**<Keep a short history of merr*/
	
    int initialized;    /**<is this data initialized*/
    int alint;          /**<Integral part of latency*/
    real alfrac;        /**<Fractional latency*/
    /*Servo parameters.*/
    dmat *ap;
    dmat *ep;
    real dt;
	sho_t *sho;/*if the corrector is a simple harmonic oscillator*/
}servo_t;
//servo optimization related routines
dcell *servo_optim(real dt, long dtrat, real al, real pmargin, real f0, real zeta, int servo_type, 
   const dmat *psdin,  const dmat* sigma2n);
real servo_optim_margin(real dt, long dtrat, real al, real pmargin, real f0, real zeta);
dmat *servo_cl2ol(const dmat *psdcl, real dt, long dtrat, real al, real gain, real sigma2n);
//cmat *servo_Hol(const dmat *nu, real dt, real dtrat, real al, const dmat *gain);
real servo_residual(real *noise_amp, const dmat *psdin, real dt, long dtrat, real al, const dmat *gain, int servo_type);

//servo time domain filtering routines

servo_t *servo_new(anyarray merr, const dmat *ap, real al, real dt, const dmat *ep);
servo_t *servo_new_scalar(anyarray merr, real ap, real al, real dt, real ep);
servo_t *servo_new_sho(anyarray merr, const dmat *ap, real al, real dt, const dmat *ep, real f0, real zeta);
int servo_filter(servo_t *st, const anyarray merr);
void servo_add(servo_t *st, const anyarray madj, real alpha);
void servo_output(const servo_t *st, panyarray out);
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
    cell *dx; //status (speed, dx/dt)
    cell *x; //status (position)
    cell *ddx; //ddx (acceleration)
    cell *ytmp; //temp
};
sho_t *sho_new(real f0, real zeta);
void sho_free(sho_t *sho);
void sho_step(panyarray xout, sho_t *sho, anyarray xi, real dt, int average);
void sho_reset(sho_t *sho);
dmat *sho_filter(const dmat *xi, real dt, real f0, real zeta, int average);
#endif
