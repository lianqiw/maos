/*
  Copyright 2009-2026 Lianqi Wang
  
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
#ifndef AOS_LIB_ACCPHI_H
#define AOS_LIB_ACCPHI_H
#include "../math/mathdef.h"

/**
   \file lib/accphi.h

   Contains ray tracing routines optimized for different input/output
   formats. Notice that the OPDs are accumulated.
*/
/**
   Unified data structure for automatic selection of propagation.
*/
typedef struct propdata_t{
    /*Input.  */
    const map_t *mapin; /*Or */ loc_t *locin;
    const real *phiin; /*If not null, use it instead of mapin->p */
    
    /*Output */
    map_t *mapout;/*or */const pts_t *ptsout;/*Or */const loc_t *locout;/*Or*/const locstat_t *ostat;
    real *phiout; /*If not null, use it instead of mapout->p */
	real alpha;/*scale of value: */
	real hs;   /*range of source*/
    real thetax, thetay;/*ray angle*/
    real misregx, misregy;/*misregistration of the output grid*/
	real shiftx, shifty; /*shift of the input grid (due to wind, etc.) */
	real rot;  /*rotation*/

    int wrap;
	int transpose;
    int nooptim;/*disable optim. */
}propdata_t;
void* prop_thread(thread_t *data);/*A unified wrapper */
void prop_range(propdata_t* propdata, long start, long end);
void prop(propdata_t* propdata);
#endif
