/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   Prepares the DM caching structs that caches DM on fine sampled grid to speed
   up ray tracing with linear interpolation.

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
void prep_cachedm(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(!parms->ndm||!(parms->sim.cachedm||parms->plot.run)){
		warning("DM cache is not needed\n");
		return;
	} else{
		info("DM cache with grid");
	}
	if(!simu->cachedm){
		simu->cachedm=mapcellnew(parms->ndm, 1);
		for(int idm=0; idm<parms->ndm; idm++){
			real dx=parms->dm[idm].dx/(parms->sim.cachedm>3?parms->sim.cachedm:4);
			info(" 1/%gm", 1./dx);
			create_metapupil(&P(simu->cachedm,idm), 0, 0, parms->dirs, parms->aper.d,
				parms->dm[idm].ht+parms->dm[idm].vmisreg, dx, dx,
				0, dx, 0, 0, 0, 0);
		}
	}
	info("\n");
	//cachedm_ha doesn't help because it is not much faster than ray tracing and
	//is not parallelized as ray tracing.
	/*new scheme for ray tracing */
	simu->cachedm_prop=mycalloc(parms->ndm, thread_t*);
	simu->cachedm_propdata=mycalloc(parms->ndm, propdata_t);
	propdata_t* cpropdata=simu->cachedm_propdata;
	const int nthread=2;//operation is fast. use the same for both DMs.
	for(int idm=0; idm<parms->ndm; idm++){
		if(simu->dmrealsq){
			cpropdata[idm].mapin=P(simu->dmrealsq,idm);
		} else{
			cpropdata[idm].locin=P(simu->recon->aloc,idm);
			cpropdata[idm].phiin=P(simu->dmreal,idm);
		}
		cpropdata[idm].mapout=P(simu->cachedm,idm);
		cpropdata[idm].alpha=1;
		cpropdata[idm].displacex0=0;
		cpropdata[idm].displacey0=0;
		cpropdata[idm].displacex1=0;
		cpropdata[idm].displacey1=0;
		cpropdata[idm].scale=1;
		
		simu->cachedm_prop[idm]=thread_prep(0, NY(cpropdata[idm].mapout),
			nthread, prop, simu->cachedm_propdata+idm);
	}
}

/**
   Partition the ray tracing by DM/Destination combinations, as well as
   segments in each combination to maximum efficiency.
*/
void calc_cachedm(sim_t* simu){
	if(simu->cachedm){
		real tk_start=myclockd();
		dcellzero((dcell*)simu->cachedm);
		CALL_THREAD_ARR(simu->cachedm_prop, simu->parms->ndm, 1);
		simu->tk_cache=myclockd()-tk_start;
	}
}
