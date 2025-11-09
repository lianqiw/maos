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
    const map_t *mapin;
    /*Or */
    loc_t *locin;
    const real *phiin;
    
    /*Output */
    map_t *mapout;
    /*Or */
    real *phiout;
    /*Combined with */
    const pts_t *ptsout;
    /*Or */
    const loc_t *locout; 
    /*Or  */
    const locstat_t *ostat;
	real alpha;/*scale of value: */
    real displacex, displacey;/*Constant displacement */
    real displacex2, displacey2;/*Time step dependent displacement */
    real scale;/*scale of coordinate */
	real rot;  /*rotation*/

    int wrap;
	int transpose;
    int nooptim;/*disable optim. */
    int index;
}propdata_t;
void* prop_thread(thread_t *data);/*A unified wrapper */
void prop(propdata_t* propdata, long start, long end);
#endif
