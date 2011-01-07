/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_RECON_H
#define AOS_RECON_H
#include "maos.h"

void TomoR(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha);
void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha);

void FitL(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha);
void FitR(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha);

void tomo(dcell **opdr, const PARMS_T *parms, 
	  const RECON_T *recon, const dcell *grad, int maxit);
void fit(dcell **adm, const PARMS_T *parms, 
		const RECON_T *recon, const dcell *opdr);
void focus_tracking(SIM_T*simu);
void tomofit(SIM_T *simu);
void lsr(SIM_T *simu);
void reconstruct(SIM_T *simu);
#endif
