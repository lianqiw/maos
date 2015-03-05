/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "common.h"
#include "sim.h"
/**
   \file maos/cachedm.c Contains routines that prepare and carry out DM shape
   caching on fine sampled grid to speed up ray tracing.  


   Because the DM has coarse sampling, and cubic influence functions, directly
   interpolate from the DM to the WFS and science OPD takes a lot of
   computation. We can speed up this process by first propagating the DM
   commands to a square grid that matches the sampling of the WFS or science
   OPD, and then propagate from the square grid to different WFS or science
   directions. For each sampling and each DM, we only need one fine sampled
   square grid that will cover all the different directions in the FoV of the
   WFSs or science. For the ground DM, we need a plane of 1/64m sampling that
   matches the WFS and science OPD. For the upper DM, we need to two planes, one
   at 1/64m for the LGS WFS and science OPD, and another one at 1/64 m
   *(1-11.2/90) for the LGS due to cone effect. The square grid are on axis.
 */
/**
   Prepares the DM caching structs.
*/
void prep_cachedm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->ndm || !parms->sim.cachedm){
	warning("No caching is needed\n");
	return;
    }
    if(!simu->cachedm){
	simu->cachedm=cellnew(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
	    double dx=parms->dm[idm].dx/16;
	    create_metapupil(&simu->cachedm->p[idm], 0, 0, parms->dirs, parms->aper.d,
			     parms->dm[idm].ht+parms->dm[idm].vmisreg, dx, dx,
			     0, 2, 0,0,0,0);
	}
    }
    //cachedm_ha doesn't help because it is not much faster than ray tracing and
    //is not parallelized as ray tracing.
    /*new scheme for ray tracing */
    simu->cachedm_prop=calloc(parms->ndm, sizeof(thread_t*));
    simu->cachedm_propdata=calloc(parms->ndm, sizeof(PROPDATA_T));
    PROPDATA_T *cpropdata=simu->cachedm_propdata;
    for(int idm=0; idm<parms->ndm; idm++){
	simu->cachedm_prop[idm]=calloc(NTHREAD, sizeof(thread_t));
	if(simu->dmrealsq){
	    cpropdata[idm].mapin=simu->dmrealsq->p[idm];
	}else{
	    cpropdata[idm].locin=simu->recon->aloc->p[idm];
	    cpropdata[idm].phiin=simu->dmreal->p[idm]->p;
	}
	cpropdata[idm].mapout=simu->cachedm->p[idm];
	cpropdata[idm].alpha=1;
	cpropdata[idm].displacex0=0;
	cpropdata[idm].displacey0=0;
	cpropdata[idm].displacex1=0;
	cpropdata[idm].displacey1=0;
	cpropdata[idm].scale=1;
	cpropdata[idm].cubic=simu->parms->dm[idm].cubic;
	cpropdata[idm].cubic_iac=simu->parms->dm[idm].iac;
	prop_index(&cpropdata[idm]);
	thread_prep(simu->cachedm_prop[idm], 0, cpropdata[idm].mapout->ny, 
		    NTHREAD, prop, (void*)&cpropdata[idm]);
    }
}

/**
   Partition the ray tracing by DM/Destination combinations, as well as
   segments in each combination to maximum efficiency.
*/
void calc_cachedm(SIM_T *simu){
    double tk_start=myclockd();
    if(simu->parms->sim.cachedm){
	long group=0;
	/*zero out the data. */
	for(int idm=0; idm<simu->parms->ndm; idm++){
	    dzero((dmat*)simu->cachedm->p[idm]);
	    /*do the multi-threaded ray tracing */
	    QUEUE_THREAD(group,(simu->cachedm_prop[idm]), 1);
	}
	WAIT_THREAD(group);
    }
    simu->tk_cache=myclockd()-tk_start;
}
