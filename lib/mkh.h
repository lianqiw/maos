/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_MKH_H
#define AOS_MKH_H
#include "loc.h"
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"
dsp * mkhb(LOC_T *locin, LOC_T *locout, const double *ampout,
	   double displacex, double displacey,double scale,
	   int cubic, double cubic_iac);
dsp * mkh(LOC_T *locin, LOC_T *locout, const double *ampout,
	  double displacex, double displacey,double scale,
	  int cubic, double cubic_iac);
/*
dsp * mkhb_cubic(LOC_T *locin, LOC_T *locout, const double *ampout,
		 double displacex, double displacey,
		 double scale,double cubic);
dsp * mkh_cubic(LOC_T *locin, LOC_T *locout, const double *ampout,
		double displacex, double displacey,
		double scale,double cubic);*/
dsp *mkhbin1d(dmat *xin, dmat *xout);		   
#endif
