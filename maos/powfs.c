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


#include "common.h"
#include "powfs.h"
#include "powfs_utils.h"
#include "pywfs.h"
#include "recon.h"
#include "recon_utils.h"
#undef GREEN
#define GREEN BLACK
/**
   \file powfs.h

   Setting up WFS geometry. like the subaperture location, subaperture grid
   points, physical optics detection transfer function, LGS elongation transfer
   function, etc.

   \todo isolate DTF, ETF, routines, make then generic interface and relocate to the lib folder.

   Do not use sparse interpolation to replace ray tracing for fine sampled
   destination grid, especially cubic splines. The interpolation marix takes too
   much space.

   TODO: This routine and powfs_t should only contain information about the
   simulation, not about any model used during reconstruction (RTC) to avoid
   leaking information from the "real world (simulation)" to our knowledge (RTC).
*/

#define MOVES(p,i,j) p[i]=p[j]
#define MOVED(p,n,i,j) memcpy(p+n*i, p+n*j, sizeof(real)*n)
#define MOVEPTS(pts,count,isa)			\
    if(pts){					\
	MOVES(pts->origx, count,isa);		\
	MOVES(pts->origy, count,isa);		\
    }
#define MOVELOC(loc,count,isa)			\
    if(loc){					\
	MOVES(loc->locx,count,isa);		\
	MOVES(loc->locy,count,isa);		\
    }
#define MOVEDLOC(loc,nxsa,count,isa)		\
    if(loc){					\
	MOVED(loc->locx,npsa,count,isa);	\
	MOVED(loc->locy,npsa,count,isa);	\
    }

/**
   Free the powfs geometric parameters
*/
static void
free_powfs_geom(powfs_t* powfs, int ipowfs){
	ptsfree(powfs[ipowfs].pts);
	dfree(powfs[ipowfs].saa);
	locfree(powfs[ipowfs].saloc);
	locfree(powfs[ipowfs].loc);
	dfree(powfs[ipowfs].amp);
	cellfree(powfs[ipowfs].loc_tel);
	dcellfree(powfs[ipowfs].amp_tel);
	dcellfree(powfs[ipowfs].saa_tel);
	cellfree(powfs[ipowfs].loc_dm);
	dcellfree(powfs[ipowfs].realamp);
	dcellfree(powfs[ipowfs].realsaa);
	dfree(powfs[ipowfs].sumamp);
	dfree(powfs[ipowfs].sumamp2);
	locfft_free(powfs[ipowfs].fieldstop);
}
/**
   Convert amplitude map of subapertures to normalized subaperture area. Used only in setup_powfs_geom*/
static dmat* wfsamp2saa(dmat* wfsamp, long npsa){
	real areanormfactor=1./(real)npsa;
	long nsa=NX(wfsamp)/npsa;
	dmat* saa=dnew(nsa, 1);
	for(int isa=0; isa<nsa; isa++){
		real* amp=P(wfsamp)+isa*npsa;
		real area=0;
		for(int iamp=0; iamp<npsa; iamp++){
			area+=amp[iamp]*amp[iamp];
		}
		P(saa, isa)=area*areanormfactor;
	}
	return saa;
}
/**
   Creates WFS pupil mask.
 */
void wfspupmask(const parms_t* parms, loc_t* loc, dmat* amp, int iwfs){
	long nloc=loc->nloc;
	dmat* ampmask=dnew(nloc, 1);
	for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
		int jpowfs=parms->wfs[jwfs].powfs;
		if(parms->powfs[jpowfs].lo) continue;
		real hs=parms->wfs[jwfs].hs;
		real ht=parms->atm.hmax*0.7-parms->wfs[jwfs].hc;
		real r=parms->aper.d*0.5*(1.-ht/hs)/(1.-ht/hs);
		real sx=(parms->wfs[jwfs].thetax-parms->wfs[iwfs].thetax)*ht;
		real sy=(parms->wfs[jwfs].thetay-parms->wfs[iwfs].thetay)*ht;
		loc_circle_add(ampmask, loc, sx, sy, r, 0, 1);
	}
	for(int i=0; i<nloc; i++){
		if(P(ampmask, i)<0.5) P(amp, i)=0;
	}
	dfree(ampmask);
}
static void
sa_reduce(powfs_t* powfs, int ipowfs, real thresarea){
	dmat* saa=NULL;//create a temporary sa area array
	if(powfs[ipowfs].saa_tel){
	//We expact subaperture to union of all wfs
		long nsa=P(powfs[ipowfs].saa_tel, 0)->nx;
		saa=dnew(nsa, 1);
		for(long iwfs=0; iwfs<NX(powfs[ipowfs].saa_tel); iwfs++){
			for(long isa=0; isa<nsa; isa++){
				if(P(saa, isa)<P(P(powfs[ipowfs].saa_tel, iwfs), isa)){
					P(saa, isa)=P(P(powfs[ipowfs].saa_tel, iwfs), isa);
				}
			}
		}
	} else{
		saa=ddup(powfs[ipowfs].saa);
	}
	if(dmax(saa)>1.01){
		warning("The sa area maxes to %g, which should be leq 1 (misregistration can cause this).\n",
			dmax(saa));
	}

	if(powfs[ipowfs].pts->nsa>4){
		loc_t* ptsloc=LOC(powfs[ipowfs].pts);
		loc_create_map(ptsloc);
		real dx1=1./ptsloc->dx;
		real dy1=1./ptsloc->dy;
		int changed=0;
		const real thresarea2=0.1;//secondary threshold to enable being surrounded subaperturs.
		do{//validata subapertures that have enough neighbors. Disable isolated subapertures
			changed=0;
			for(long isa=0; isa<ptsloc->nloc; isa++){
				long ix=round((ptsloc->locx[isa]-ptsloc->map->ox)*dx1);
				long iy=round((ptsloc->locy[isa]-ptsloc->map->oy)*dy1);
				int nedge=0;
				int ncorner=0;
				int nself=0;
				if(P(saa, isa)>=thresarea){
					nself=1;
				}
				for(int jx=-1; jx<2; jx++){
					for(int jy=-1; jy<2; jy++){
						if(jx==0&&jy==0){
							continue;
						}
						long jsa=loc_map_get(ptsloc->map, ix+jx, iy+jy);
						if(jsa&&P(saa, jsa-1)>=thresarea){
							if(abs(jx+jy)==1){//edge
								nedge++;
							} else{//corner
								ncorner++;
							}
						}
					}
				}
				if(nself){//disable isolated valid subaperture
					if(ncorner+nedge<=2){
						P(saa, isa)=0;
						changed++;
					}
				} else if(0){//enable isolated in-valid subaperture
					if(nedge+ncorner>=6&&P(saa, isa)>=thresarea2){
						P(saa, isa)=1;
						changed++;
					}
				}
			}
		} while(changed);
		loc_free_map(ptsloc);
	}

	int count=0;
	const int npsa=powfs[ipowfs].pts->nxsa*powfs[ipowfs].pts->nysa;
	for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
		if(P(saa, isa)>=thresarea){
			/*Area is above threshold, keep.  Shift pts, ptsm, loc, locm, amp,
			  ampm, saloc area is already normalized that maxes to 1. The MOVE*
			  are defined in the beginining of this file.*/
			if(count!=isa){
				MOVEPTS(powfs[ipowfs].pts, count, isa);
				MOVES(P(powfs[ipowfs].saa), count, isa);
				MOVELOC(powfs[ipowfs].saloc, count, isa);
				MOVEDLOC(powfs[ipowfs].loc, npsa, count, isa);
				MOVED(P(powfs[ipowfs].amp), npsa, count, isa);
				if(powfs[ipowfs].saa_tel){
					for(int jwfs=0; jwfs<NX(powfs[ipowfs].saa_tel); jwfs++){
						MOVES(P(P(powfs[ipowfs].saa_tel, jwfs)), count, isa);
						MOVEDLOC(P(powfs[ipowfs].loc_tel, jwfs), npsa, count, isa);
						MOVED(P(P(powfs[ipowfs].amp_tel, jwfs)), npsa, count, isa);
					}
				}
			}
			count++;
		}
	}
	if(count==0){
		error("there are no subapertures above threshold.\n");
	}
	dfree(saa);

	ptsresize(powfs[ipowfs].pts, count);
	locresize(powfs[ipowfs].saloc, count);
	locresize(powfs[ipowfs].loc, count*npsa);
	dresize(powfs[ipowfs].saa, count, 1);
	dresize(powfs[ipowfs].amp, count*npsa, 1);
	powfs[ipowfs].saasum=dsum(powfs[ipowfs].saa);
	if(powfs[ipowfs].loc_tel){
		for(int jwfs=0; jwfs<NX(powfs[ipowfs].saa_tel); jwfs++){
			locresize(P(powfs[ipowfs].loc_tel, jwfs), count*npsa);
			dresize(P(powfs[ipowfs].amp_tel, jwfs), count*npsa, 1);
			dresize(P(powfs[ipowfs].saa_tel, jwfs), count, 1);
		}
	}
}
/**
 * wrap val to between low and high
 * */
real wrap2range(real val, real low, real high){
	if(low>=high) {
		warning("Invalid usage: low=%g should be smaller than high %g\n", low ,high);
		return val;
	}
	return remainder(val-low, high-low)+low;
}
/**
   setting up subaperture geometry.

   - powfs->pts: the location of the lower left grid coordinate in this subaperture >=saloc.

   - powfs->saloc: the lower left corner of the subaperture

   - powfs->loc: the x,y coordinates of all the points, ordered acoording to subaperture.

   2011-02-11: In the matched filter and geometric gradient computation, the
   amplitude should use the real (misregistered) amplitude map. The geometric
   optics operator should just the normal loc but misregistered amplitude map
   (on locm).  The subaperture selector should, however, use assumed amp
   (non-misregistered).

*/
static void
setup_shwfs_geom(powfs_t* powfs, const parms_t* parms,
	aper_t* aper, int ipowfs){
	free_powfs_geom(powfs, ipowfs);
	int nwfsp=parms->powfs[ipowfs].nwfs;
	/*order of the system. 60 for TMT */
	const int order=parms->powfs[ipowfs].order;
	/*Subaperture lateral length. 0.5 for TMT. */
	const real dsa=parms->powfs[ipowfs].dsa;
	/*The subaperture area that are considered as full. */
	real areafulli;
	if(order>2){
		areafulli=1.;
	} else{
		areafulli=4./M_PI;
	}

	/*Number of OPD pixels in 1 direction in each
	  subaperture. Make it even to do fft.*/
	int nxsa=2*(int)round(0.5*dsa/parms->powfs[ipowfs].dx);
	const real dx=dsa/nxsa;/*adjust dx. */
	const real dxoffset=dx*0.5;//Always keep points inside subaperture for simulation.
	if(fabs(parms->powfs[ipowfs].dx-dx)>EPS)
		info("Adjusting dx from %g to %g\n", parms->powfs[ipowfs].dx, dx);
	if(fabs(dsa-nxsa*dx)>EPS){
		warning("nx=%d,dsa=%f,dx=%f not agree\n", nxsa, dsa, dx);
	}
	info("There are %d points in each subaperture of %gm.\n", nxsa, dsa);
	const int npsa=nxsa*nxsa;/*Total Number of OPD points. */
	if(parms->powfs[ipowfs].saloc){
		powfs[ipowfs].saloc=locread("%s", parms->powfs[ipowfs].saloc);
		if(fabs(powfs[ipowfs].saloc->dx-dsa)>1e-6*fabs(dsa)){
			error("loaded saloc has dx=%g, while powfs.dsa=%g\n",
				powfs[ipowfs].saloc->dx, dsa);
		}
	} else{
		/*The coordinate of the subaperture (lower left coordinate) */
		const int orderpad=order+(order>2?2:0); //to avoid missing subapertures when there is saoffx/saoffy
		powfs[ipowfs].saloc=locnew(orderpad*orderpad, dsa, dsa);
		int count=0;
		/*Offset of the coordinate of the center most subaperture from the center. */
		real offsetx=wrap2range(parms->powfs[ipowfs].misregx,-0.5,0.5)+((order&1)?-0.5:0);
		real offsety=wrap2range(parms->powfs[ipowfs].misregy,-0.5,0.5)+((order&1)?-0.5:0);
		dbg("offsetx=%g, offsety=%g, order=%d, orderpad=%d\n", offsetx, offsety, order, orderpad);
		/*r2max: Maximum distance^2 from the center to keep a subaperture */
		real r2max=pow(order*0.5, 2);
		real r2min=dsa<parms->aper.din?pow(parms->aper.din/dsa/2, 2):-1;
		/*the lower left *grid* coordinate of the subaperture */

		/*Collect all the subapertures that are within the allowed radius*/
		for(int j=-orderpad/2; j<=(orderpad-1)/2; j++){
			for(int i=-orderpad/2; i<=(orderpad-1)/2; i++){
			//Normalized coordinate in uniq of sa size
				real xc=((real)i+offsetx);
				real yc=((real)j+offsety);
				//Radius of four corners.
				real r1=pow(xc+1, 2)+pow(yc+1, 2);
				real r2=pow(xc+1, 2)+pow(yc, 2);
				real r3=pow(xc, 2)+pow(yc+1, 2);
				real r4=pow(xc, 2)+pow(yc, 2);
				if(r1<r2max||r2<r2max||r3<r2max||r4<r2max||
					r1>r2min||r2>r2min||r3>r2min||r4>r2min){
					powfs[ipowfs].saloc->locx[count]=xc*dsa;
					powfs[ipowfs].saloc->locy[count]=yc*dsa;
					count++;
				}
			}
		}
		//info("count=%d\n", count);
		powfs[ipowfs].saloc->nloc=count;
	}
	/*convert saloc to pts*/
	powfs[ipowfs].pts=ptsnew(powfs[ipowfs].saloc->nloc, dsa, dsa, nxsa, nxsa, dx, dx);
	for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
		powfs[ipowfs].pts->origx[isa]=powfs[ipowfs].saloc->locx[isa]+dxoffset;
		powfs[ipowfs].pts->origy[isa]=powfs[ipowfs].saloc->locy[isa]+dxoffset;
	}
	/*Calculate the amplitude for each subaperture OPD point by ray tracing from
	  the pupil amplitude map. Pupil distortion is accounted for.*/
	powfs[ipowfs].loc=pts2loc(powfs[ipowfs].pts);
	/*The assumed amp. */
	if(parms->powfs[ipowfs].amp){
		powfs[ipowfs].amp=dread("%s", parms->powfs[ipowfs].amp);
		if(NX(powfs[ipowfs].amp)!=powfs[ipowfs].loc->nloc){
			error("%s is in wrong format. Need %ld, has %ld.\n",
				parms->powfs[ipowfs].amp, powfs[ipowfs].loc->nloc, NX(powfs[ipowfs].amp));
		}
	} else{
		powfs[ipowfs].amp=mkamp(powfs[ipowfs].loc, aper->ampground,
			-P(parms->aper.misreg,0), -P(parms->aper.misreg,1),
			parms->aper.d, parms->aper.din);
	}
	/*The threashold for normalized area (by areafulli) to keep subaperture. */
	real thresarea=parms->powfs[ipowfs].saat;
	dmat* ampi=NULL;
	if(parms->powfs[ipowfs].safill2d<1){
	/*subaperture amplitude map to simulate lenslet fill factor*/
		int nedge=1;
		while((nxsa-2*nedge)*(nxsa-2*nedge)>nxsa*nxsa*parms->powfs[ipowfs].safill2d){
			nedge++;
		}
		ampi=dnew(nxsa, nxsa);
		real alpha=(nxsa*nxsa*parms->powfs[ipowfs].safill2d-(nxsa-2*nedge)*(nxsa-2*nedge))
			/((nxsa-2*(nedge-1))*(nxsa-2*(nedge-1))-(nxsa-2*nedge)*(nxsa-2*nedge));
		//real tot=0;
		for(int iy=nedge-1; iy<nxsa-nedge+1; iy++){
			for(int ix=nedge-1; ix<nxsa-nedge+1; ix++){
				P(ampi, ix, iy)=1;
				if(ix==nedge-1||ix==nxsa-nedge||iy==nedge-1||iy==nxsa-nedge){
					P(ampi, ix, iy)=alpha;
				}
				//tot+=P(ampi, ix, iy);
			}
		}
		if(parms->save.setup){
			writebin(ampi, "powfs%d_ampi", ipowfs);
		}
		for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
			for(int i=0; i<nxsa*nxsa; i++){
				P(powfs[ipowfs].amp, nxsa*nxsa*isa+i)*=P(ampi, i);
			}
		}
		/*do not multiply to siglev. Already handled automatically*/
		thresarea*=parms->powfs[ipowfs].safill2d;
	}

	powfs[ipowfs].saa=wfsamp2saa(powfs[ipowfs].amp, npsa);


	setup_powfs_misreg_tel(powfs, parms, aper, ipowfs);
	/*Go over all the subapertures, calculate the normalized
	  subaperture illumination area and remove all that are below
	  the are threshold*/

	/*
	  About physical optics imaging: The subaperture area (saa) has been
	  normalized against a square subaperture to scale the physical optics
	  i0. After FFT, the PSF sum to 1 for full square subapertures, but sum to
	  PI/4 for full circular TT or full quadrant circular TTF subapertures. The
	  areascale variable (1 for square subaperture and 4/M_PI for (quad)circular
	  sa) is therefore used to normalize the PSF so that it sums to 1 for either
	  full square or full (quadrant) circular subaperture. Used in wfsints.

	  The subaperture area (saa) is subsequently scaled so that the max is 1 for
	  full square or (quadrant) circular subaperture. saa is then used in
	  setup_recon to scale geometric sanea.
	*/

	powfs[ipowfs].areascale=areafulli;
	if(fabs(areafulli-1)>EPS){
		dscale(powfs[ipowfs].saa, areafulli);
		dcellscale(powfs[ipowfs].saa_tel, areafulli);
	}
	if(!parms->powfs[ipowfs].saloc&&!parms->powfs[ipowfs].amp){
		sa_reduce(powfs, ipowfs, thresarea);
	}
	info2("There are %ld valid subaperture.\n", powfs[ipowfs].saloc->nloc);
	setup_powfs_misreg_dm(powfs, parms, aper, ipowfs);
	powfs[ipowfs].realsaa=dcellnew(nwfsp, 1);
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
		if(powfs[ipowfs].loc_tel){
			P(powfs[ipowfs].realsaa, jwfs)=dref(P(powfs[ipowfs].saa_tel, jwfs));
		} else{
			P(powfs[ipowfs].realsaa, jwfs)=dref(powfs[ipowfs].saa);
		}
	}
	powfs[ipowfs].realamp=dcellnew(nwfsp, 1);
	powfs[ipowfs].sumamp=dnew(nwfsp, 1);
	powfs[ipowfs].sumamp2=dnew(nwfsp, 1);
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
		dmat* realamp;
		if(powfs[ipowfs].loc_tel){
			realamp=P(powfs[ipowfs].amp_tel, jwfs);
		} else{
			realamp=powfs[ipowfs].amp;
		}
		real sumamp2=0;
		real sumamp=0;
		for(long i=0; i<powfs[ipowfs].loc->nloc; i++){
			sumamp2+=P(realamp, i)*P(realamp, i);
			sumamp+=P(realamp, i);
		}
		P(powfs[ipowfs].sumamp2, jwfs)=sumamp2;
		P(powfs[ipowfs].sumamp, jwfs)=sumamp;
		P(powfs[ipowfs].realamp, jwfs)=dref(realamp);
	}
	if(parms->powfs[ipowfs].fieldstop){
		warning("powfs%d: generating field stop \n", ipowfs);
		if(parms->powfs[ipowfs].nwvl>1){
			error("Not implemented yet. need to do phase unwrap in wfsgrad.\n");
		}
		powfs[ipowfs].fieldstop=locfft_init(powfs[ipowfs].loc, powfs[ipowfs].amp,
			parms->powfs[ipowfs].wvl, NULL, 2,
			parms->powfs[ipowfs].fieldstop);
	}
	dfree(ampi);

	if(parms->save.setup){
		locwrite(powfs[ipowfs].pts, "powfs%d_pts", ipowfs);
		locwrite(powfs[ipowfs].saloc, "powfs%d_saloc", ipowfs);
		writebin(powfs[ipowfs].saa, "powfs%d_saa", ipowfs);
		locwrite(powfs[ipowfs].loc, "powfs%d_loc", ipowfs);
		writebin(powfs[ipowfs].amp, "powfs%d_amp", ipowfs);
		if(powfs[ipowfs].loc_tel){
			writebin(powfs[ipowfs].saa_tel, "powfs%d_saa_tel", ipowfs);
			writebin(powfs[ipowfs].amp_tel, "powfs%d_amp_tel", ipowfs);
			writebin(powfs[ipowfs].loc_tel, "powfs%d_loc_tel", ipowfs);
		}
		if(powfs[ipowfs].loc_dm){
			for(int idm=0; idm<parms->ndm; idm++){
				for(int jwfs=0; jwfs<nwfsp; jwfs++){
					writebin(powfs[ipowfs].loc_dm, "powfs%d_loc_dm", ipowfs);
				}
			}
		}
	}
}
/**
   Setup telescope to WFS pupil misregistration.

   It is tricky to implement the misregistration of the deformable mirror and the WFS.

   Misregistration of the deformable mirror.

   Denote the grid of actuators as saloc, and with misregistration, it becomes
   salocm.

   - We should use saloc for DM fitting and pseudo open loop gradient
   computations because we are not aware of the misregistration.

   - We should use salocm for raytracing to WFS and performance evaluation
   because this is what happens in the system, no matter whether we know the
   misregistration or not.

   Misregistration of the WFS.

   Denote the grid of the WFS entry pupil as loc, and with misregistration, it
   becomes locm. Project the telescope amplitude map onto loc, we get amp, which
   is what we know. Project the telescope amplitude map onto locm, we get ampm,
   which is what happens but we don't know.

   - We should use loc and amp to compute the reconstructor.

   - We should use locm and ampm for tracing to WFS and performance evaluation.

   - We should use loc and ampm for matched filter, averaging gradient, or
   zernike tip/tilt operator computation, since this is the process that happens
   in the system, and the matched filter is updated through dithering to the
   true amplitude map. We use loc because it is our model within the WFS.

   - We should use ptsm->area (real area) to select subapertures for use in
   reconstruction, because in real system, the subapertures are selected by
   their illumination level, and is affected by the real amplitude.


 */
void
setup_powfs_misreg_tel(powfs_t* powfs, const parms_t* parms, aper_t* aper, int ipowfs){
	int nwfsp=parms->powfs[ipowfs].nwfs;
	if(parms->distortion.tel2wfs){
		TIC;tic;
		/*
		  Misregistration/distortion from Telescope pupil to WFS pupil. The
		  amplitude map after misregistration/distortion will be used for
		  wavefront sensing.
		  They are not used for wavefront reconstruction
		*/
		int isset=0;
		powfs[ipowfs].loc_tel=(loccell*)cellnew(nwfsp, 1);
		powfs[ipowfs].amp_tel=dcellnew(nwfsp, 1);
		if(parms->powfs[ipowfs].type==WFS_SH){
			powfs[ipowfs].saa_tel=dcellnew(nwfsp, 1);
		}
OMP_TASK_FOR(4)
		for(int jwfs=0; jwfs<nwfsp; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
			if(parms->distortion.tel2wfs[iwfs]){
				isset=1;
				P(powfs[ipowfs].loc_tel, jwfs)
					=loctransform(powfs[ipowfs].loc, parms->distortion.tel2wfs[iwfs]);
				P(powfs[ipowfs].amp_tel, jwfs)
					=mkamp(P(powfs[ipowfs].loc_tel, jwfs), aper->ampground,
						-P(parms->aper.misreg,0), -P(parms->aper.misreg,1),
						parms->aper.d, parms->aper.din);
				if(parms->powfs[ipowfs].type==WFS_SH){
					const int npsa=powfs[ipowfs].pts->nxsa*powfs[ipowfs].pts->nysa;
					P(powfs[ipowfs].saa_tel, jwfs)=wfsamp2saa(P(powfs[ipowfs].amp_tel, jwfs), npsa);
				}
			}
		}
		if(!isset){
			cellfree(powfs[ipowfs].loc_tel);
			dcellfree(powfs[ipowfs].saa_tel);
			dcellfree(powfs[ipowfs].amp_tel);
		} else{
			toc2("distortion.tel2wfs");
		}
	}/*if misreg */
}
/**
   setup DM to WFS misregistration.
*/
void
setup_powfs_misreg_dm(powfs_t* powfs, const parms_t* parms, aper_t* aper, int ipowfs){
	(void)aper;
	int nwfsp=parms->powfs[ipowfs].nwfs;
	if(parms->distortion.dm2wfs){
		TIC;tic;
		/*
		  Misregistration/distortion from DM to WFS pupil. Not about telescope
		  pupil. The distorted grid are used for ray tracing from DM to WFS.
		  They are not used for wavefront reconstruction
		*/
		powfs[ipowfs].loc_dm=loccellnew(nwfsp,parms->ndm);
		int isset=0;
OMP_TASK_FOR_COLLAPSE(2, NTHREAD)
		for(int idm=0; idm<parms->ndm; idm++){
			for(int jwfs=0; jwfs<nwfsp; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				if(parms->distortion.dm2wfs[iwfs+idm*parms->nwfs]){
					P(powfs[ipowfs].loc_dm, jwfs, idm)
						=loctransform(powfs[ipowfs].loc, parms->distortion.dm2wfs[iwfs+idm*parms->nwfs]);
					isset=1;
				}
			}
		}
		if(!isset){
			cellfree(powfs[ipowfs].loc_dm);
		} else{
			toc2("distortion.dm2wfs");
		}
	}
}
/**
   Creating geometric wavefront gradient operator GS0 and ZA0 from WFS OPD to
   subaperture grads. Use loc, and ampm. Do not locm. locm is only used to do
   ray tracing.  This is simulation, not RTC, so use real system distortion,
   misregistration, etc. */
static void
setup_shwfs_grad(powfs_t* powfs, const parms_t* parms, int ipowfs){
	if(parms->powfs[ipowfs].gtype_recon==GTYPE_G||parms->powfs[ipowfs].gtype_sim==GTYPE_G){
		dspcellfree(powfs[ipowfs].GS0);
		/*Setting up every gradient tilt (g-ztilt) */
		if(parms->load.GS0){
			powfs[ipowfs].GS0=dspcellread("powfs%d_GS0", ipowfs);
			if(powfs[ipowfs].amp_tel){
				assert(NX(powfs[ipowfs].GS0)==parms->powfs[ipowfs].nwfs);
			} else{
				assert(NX(powfs[ipowfs].GS0)==1);
			}
		} else{
			/*This mkg takes about 5 seconds. */
			if(powfs[ipowfs].amp_tel){
				powfs[ipowfs].GS0=dspcellnew(parms->powfs[ipowfs].nwfs, 1);
			} else{
				powfs[ipowfs].GS0=dspcellnew(1, 1);
			}
			for(int iwfs=0; iwfs<NX(powfs[ipowfs].GS0); iwfs++){
				P(powfs[ipowfs].GS0, iwfs)=mkg(powfs[ipowfs].loc,
					powfs[ipowfs].loc,
					P(powfs[ipowfs].realamp, iwfs),
					powfs[ipowfs].saloc,
					1, 0, 0, 1);
			}
			if(parms->save.setup&&powfs[ipowfs].GS0){
				writebin(powfs[ipowfs].GS0, "powfs%d_GS0", ipowfs);
			}
		}
	}
	if(parms->powfs[ipowfs].gtype_recon==GTYPE_Z||parms->powfs[ipowfs].gtype_sim==GTYPE_Z){
	/*setting up zernike best fit (ztilt) inv(M'*W*M). good for NGS. */
		if(parms->powfs[ipowfs].order>4){
			warning("Ztilt for high order powfs %d is not good\n", ipowfs);
		}
		powfs[ipowfs].nsaimcc=MAX(1, (powfs[ipowfs].loc_tel?parms->powfs[ipowfs].nwfs:1));
		int nsaimcc=powfs[ipowfs].nsaimcc;
		cellfree(powfs[ipowfs].saimcc);
		powfs[ipowfs].saimcc=(dccell*)cellnew(nsaimcc, 1);
		for(int imcc=0; imcc<nsaimcc; imcc++){
			dcell* mcc=pts_mcc_ptt(powfs[ipowfs].pts, P(P(powfs[ipowfs].realamp, imcc)));
			P(powfs[ipowfs].saimcc, imcc)=dcellinvspd_each(mcc);
			dcellfree(mcc);
		}
	}
}
void setup_powfs_neasim(const parms_t* parms, powfs_t* powfs){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		const long nsa=powfs[ipowfs].saloc->nloc;
		const int ng=parms->powfs[ipowfs].ng;
		dcell* nea=0;
		//if(parms->powfs[ipowfs].neaphy || parms->powfs[ipowfs].phystep>-1){
		if(powfs[ipowfs].sanea){
			info("Use sanea to derive neasim\n");
			nea=dcelldup(powfs[ipowfs].sanea);
			for(int ii=0; ii<NX(nea); ii++){
				nea_chol(&P(nea, ii), P(nea, ii),ng);
			}
		} else{
			//neasimfile/neareconfile is saved by skyc in radian.
			char* file=parms->powfs[ipowfs].neasimfile;
			if(!file&&parms->powfs[ipowfs].neasim==-1){
				file=parms->powfs[ipowfs].neareconfile;
			}
			if(file){
				nea=dcellread_prefix(file, parms, ipowfs);
			}
		}
		if(nea){
			for(int ii=0; ii<NX(nea)*NY(nea); ii++){
				nea_check(P(nea, ii), nsa, ng);
			}
		} else{
			int nnea=powfs[ipowfs].loc_tel?parms->powfs[ipowfs].nwfs:1;
			nea=dcellnew(nnea, 1);
			for(int jwfs=0; jwfs<nnea; jwfs++){
				real nea_rad;
				if(parms->powfs[ipowfs].neasim<0){
					nea_rad=parms->powfs[ipowfs].nearecon;//in mas
				} else{
					nea_rad=parms->powfs[ipowfs].neasim;//in mas
				}
				nea_rad=nea_rad*MAS2RAD/sqrt(parms->powfs[ipowfs].dtrat*parms->sim.dt/parms->sim.dtref);//in rad
				real* saa=powfs[ipowfs].realsaa?P(P(powfs[ipowfs].realsaa, jwfs)):P(powfs[ipowfs].saa);
				dmat* nea_each=P(nea, jwfs)=dnew(nsa, 3);
				for(int isa=0; isa<nsa; isa++){
					P(nea_each, isa, 0)=P(nea_each, isa, 1)=nea_rad/sqrt(saa[isa]);
				}
			}
		}
		powfs[ipowfs].neasim=nea; nea=0;
		if(parms->save.setup){
			writebin(powfs[ipowfs].neasim, "powfs%d_neasim", ipowfs);
		}
	}
}
/**
   Prepare the parameters for physical optics setup.

   \todo Make LGS number of pixels depend on subapertue elongation in radial
   coordinate CCD mode. Make notfx, notfy dependent on number of pixels for radial ccd.

*/
static void
setup_powfs_prep_phy(powfs_t* powfs, const parms_t* parms, int ipowfs){
	const real pixthetax=parms->powfs[ipowfs].radpixtheta;
	const real pixthetay=parms->powfs[ipowfs].pixtheta;
	const int pixpsay=parms->powfs[ipowfs].pixpsa;
	const int radpix=parms->powfs[ipowfs].radpix;
	const real dsa=powfs[ipowfs].pts->dsa;
	const int nsa=powfs[ipowfs].saloc->nloc;
	if(parms->powfs[ipowfs].llt){
		const int nllt=parms->powfs[ipowfs].llt->n;
		real rsa2, rsa2max=0;
		dcellfree(powfs[ipowfs].srot);
		dcellfree(powfs[ipowfs].srsa);
		cellfree(powfs[ipowfs].sprint);

		powfs[ipowfs].srot=dcellnew(nllt, 1);
		powfs[ipowfs].srsa=dcellnew(nllt, 1);
		powfs[ipowfs].srsamax=dnew(nllt, 1);
		powfs[ipowfs].sprint=lcellnew(nllt, 1);

		for(int illt=0;illt<nllt;illt++){
			/*adjusted llt center because saloc->locx/y is corner */
			real ox2=P(parms->powfs[ipowfs].llt->ox, illt)-dsa*0.5;
			real oy2=P(parms->powfs[ipowfs].llt->oy, illt)-dsa*0.5;
			P(powfs[ipowfs].srot, illt)=dnew(nsa, 1);
			P(powfs[ipowfs].srsa, illt)=dnew(nsa, 1);

			for(int isa=0; isa<nsa; isa++){
				real ddx=(powfs[ipowfs].saloc->locx[isa]-ox2);
				real ddy=(powfs[ipowfs].saloc->locy[isa]-oy2);
				rsa2=pow(ddx, 2)+pow(ddy, 2);
				P(P(powfs[ipowfs].srot, illt), isa)=atan2(ddy, ddx);
				P(P(powfs[ipowfs].srsa, illt), isa)=sqrt(rsa2);
				if(rsa2>rsa2max) rsa2max=rsa2;
			}
			P(powfs[ipowfs].srsamax, illt)=sqrt(rsa2max);
			real dprint=MAX(parms->aper.d*0.1, dsa);
			int pnsa=(int)ceil(sqrt(rsa2max)/dprint)+1;

			real prot[pnsa];
			P(powfs[ipowfs].sprint, illt)=lnew(pnsa, 1);
			long* pp=P(P(powfs[ipowfs].sprint, illt));
			for(int ind=0; ind<pnsa; ind++){
				prot[ind]=INFINITY;
				pp[ind]=-1;
			}
			real desrot=0;//LLT polar angle
			if(fabs(P(parms->powfs[ipowfs].llt->ox, illt))<dsa
				&&fabs(P(parms->powfs[ipowfs].llt->oy, illt))<dsa){
				desrot=0;
			} else{
				real ddx=(0-P(parms->powfs[ipowfs].llt->ox, illt));
				real ddy=(0-P(parms->powfs[ipowfs].llt->oy, illt));
				desrot=atan2(ddy, ddx);
			}
			dmat* saa=PR(powfs[ipowfs].realsaa, illt, 0);
			for(int isa=0; isa<nsa; isa++){
				int ind=(int)round(P(P(powfs[ipowfs].srsa, illt), isa)/dprint);
				real irot=fabs(P(P(powfs[ipowfs].srot, illt), isa)-desrot);
				if(ind>=pnsa){
					error("ind=%d>=pnsa=%d\n", ind, pnsa);
				}
				if(irot<prot[ind]&&P(saa, isa)>0.9){
					prot[ind]=irot;
					pp[ind]=isa;
				}
			}
		}
		if(parms->save.setup){
			writebin(powfs[ipowfs].sprint, "powfs%d_sprint", ipowfs);
			writebin(powfs[ipowfs].srot, "powfs%d_srot", ipowfs);
			writebin(powfs[ipowfs].srsa, "powfs%d_srsa", ipowfs);
		}
	}
	int pixpsax=pixpsay;
	if(radpix<0){/*determine radpix adaptively */
	/*use rsa2max to determine radpix */
		error("Not emplemented yet\n");
	}
	pixpsax=(radpix)?radpix:pixpsay;

	/*record # of detector pixels. */
	powfs[ipowfs].pixpsax=pixpsax;
	powfs[ipowfs].pixpsay=pixpsay;

	{
	/*to avoid aliasing, fft warping. usually 2. */
		const real embfac=parms->powfs[ipowfs].embfac;
		real safov=MAX(pixpsax*pixthetax, pixpsay*pixthetay);
		real wvlmin=P(parms->powfs[ipowfs].wvl, 0);
		real dtheta=wvlmin/(dsa*embfac);/*PSF sampling. */

		/*size required to form detector image. */
		int notf;
		const char* okind;
		if(parms->powfs[ipowfs].notf){
			notf=(parms->powfs[ipowfs].notf+1)/2*2;
			okind="input";
		} else{
			notf=MAX(powfs[ipowfs].pts->nxsa*embfac, ceil(safov/dtheta));
			/*
			if(notf>8 && notf<16){
			notf=16;
			}else if(notf>16 && notf<32){
			notf=32;
			}else if(notf>64 && notf<67){
			notf=64;
			}else if(notf<128 && notf>120){
			notf=128;
			}*/
			notf=nextfftsize(notf);
			okind="calculated";
		}
		powfs[ipowfs].notfx=notf;
		powfs[ipowfs].notfy=notf;
		info("OTF dimension is %dx%d (%s)\n", powfs[ipowfs].notfx, powfs[ipowfs].notfy, okind);
		if(safov>dtheta*notf){
			warning("Subaperture PSF size (%.3f\") is smaller than detector FoV (%.3f\").\n",
				dtheta*notf*RAD2AS, safov*RAD2AS);
		}
	}/*notf */


	if(parms->powfs[ipowfs].bkgrndfn){
		char* fn=parms->powfs[ipowfs].bkgrndfn;
		info("Loading sky background/rayleigh backscatter from %s\n", fn);
		dcellfree(powfs[ipowfs].bkgrnd);
		powfs[ipowfs].bkgrnd=dcellread("%s", fn);
	}
	if(parms->powfs[ipowfs].bkgrndfnc){
		char* fn=parms->powfs[ipowfs].bkgrndfnc;
		info("Loading sky background/rayleigh backscatter correction from %s\n", fn);
		dcellfree(powfs[ipowfs].bkgrndc);
		powfs[ipowfs].bkgrndc=dcellread("%s", fn);
	}
	if(parms->powfs[ipowfs].bkgrndfn||parms->powfs[ipowfs].bkgrndfnc){
		real bkscale=(parms->sim.dt/parms->sim.dtref)*parms->powfs[ipowfs].dtrat;
		if(fabs(bkscale-1)>1.e-20){
			dcellscale(powfs[ipowfs].bkgrnd, bkscale);
			dcellscale(powfs[ipowfs].bkgrndc, bkscale);
			warning("Scaling bkgrndfn by %g", bkscale);
		}
		if(parms->powfs[ipowfs].bkgrndfn&&
			(NX(powfs[ipowfs].bkgrnd)!=powfs[ipowfs].saloc->nloc
			 ||(NY(powfs[ipowfs].bkgrnd)!=1
				&&NY(powfs[ipowfs].bkgrnd)!=parms->powfs[ipowfs].nwfs))){
			error("powfs%d: bkgrnd is of dimension %ld x %ld, "
				"but should be %ld x 1 or %d\n",
				ipowfs, NX(powfs[ipowfs].bkgrnd), NY(powfs[ipowfs].bkgrnd),
				powfs[ipowfs].saloc->nloc, parms->powfs[ipowfs].nwfs);
		}
		if(parms->powfs[ipowfs].bkgrndfnc&&
			(NX(powfs[ipowfs].bkgrndc)!=powfs[ipowfs].saloc->nloc
			 ||(NY(powfs[ipowfs].bkgrndc)!=1
				&&NY(powfs[ipowfs].bkgrndc)!=parms->powfs[ipowfs].nwfs))){
			error("powfs%d: bkgrndc is of dimension %ld x %ld, "
				"but should be %ld x 1 or %d\n",
				ipowfs, NX(powfs[ipowfs].bkgrndc), NY(powfs[ipowfs].bkgrndc),
				powfs[ipowfs].saloc->nloc, parms->powfs[ipowfs].nwfs);
		}
	}
}

/**
   Setting up Detector transfer function used to translate PSFs from FFTs to
   detector pixel images. Integration over pixel size and leaking between pixels
   are considered.

   We are going to sample the PSF \f$\phi(\theta)\f$ onto detectors. For a detector
   pixel \f$i\f$ at \f$\theta_{i}\f$, oriented along angle \f$\alpha\f$, we have
   \f[
   I(\theta_{i})=\int\textrm{d}\theta\phi(\theta)S(\theta-\theta_{i})=\int\textrm{d}\theta\phi(\theta-\theta_{i})S(\theta)\f]
   where \f$S(\theta)\f$ is the pixel weighting function.  Rewrite
   the equation in Fourier domain,

   \f{eqnarray*}{ I(\theta_{i}) & = &
   \int\textrm{d}\theta\int\textrm{d} u\hat{\phi}(u)e^{2\pi i(\theta-\theta_{i})u}S(\theta)\\ & =
   & \int\textrm{d} u\hat{\phi}(u)\hat{S}(-u)e^{2\pi i\theta_{i}u}\\ & = &
   \mathcal{F}^{-1}[\hat{\phi}(u)\hat{S}(-u)](\theta_{i})\f}

   For a polar coordinate CCD with pixel edges alonging , we have

   The detector pixels can be modeled as a square, convolved with a Gaussian to
   model the pixel blurring (charge leaking)
   \f[
   S(\theta)=S^{\prime}(\theta_{ra})=\sqcap(\theta_{r}/\Delta)\sqcap(\theta_{a}/\Delta)*\frac{1}{2\pi
   \theta_b^2}\exp(-\frac{\theta_r^2+\theta_u^2}{2\theta_b^2})
   \f], where \f$R\f$ is the rotation matrix and
   \f$\theta_{ra}=R^{T}\theta\f$ is the coordinate in the local radial-azimuthal
   coordinate system for a Polar coordinate CCD, \f$\theta_b\f$ is a measure of the blurring of the pixel. For a pixel aligned along \f$x/y\f$ directions, \f$R\f$ is simply identity.

   The Fourier transform of the detector pixel function \f$S(\theta)\f$ is then

   \f{eqnarray*}{
   \hat{S}(u) & = & \int\textrm{d}\theta S(\theta)e^{-2\pi i\theta^{T}u}\\
   & = & \int\textrm{d}\theta
   S^{\prime}(R^{T}\theta)e^{-2\pi i\theta^{T}u}\\
   & = &\int\textrm{d}(R^{T}\theta)S^{\prime}(R^{T}\theta)e^{-2\pi i(R^{T}\theta)^{T}(R^{T}u)}\\
   & = & \hat{S^{\prime}}(R^{T}u)\f}

   The Fourier transform of \f$\hat{S^{\prime}}\f$ is
   \f[
   \hat{S^{\prime}}(u_{ra})=\Delta^{2}\textrm{sinc}(u_{r}\Delta)\textrm{sinc}(u_{\theta}\Delta)
   \exp[-2\pi^2\theta_b^2(u_r^2+u_a^2)]
   \f]

   Finally, we have \f[
   I(\theta_{i})=\mathcal{F}^{-1}[\hat{\phi}(u)\Delta^{2}\textrm{sinc}(\Delta(u_{x}\cos\theta+u_{y}\sin\theta))\textrm{sinc}(\Delta(-u_{x}\sin\theta+u_{y}\cos\theta))\exp[-2\pi^2\theta_b^2(u_r^2+u_a^2)](\theta_{i})\f]
   the last \f$\theta_{i}\f$ means evaluate the function at \f$\theta_{i}\f$ by
   interpolation. When there is leaking between actuators. We can multiply the
   sinc functions with a Gaussian blurring function.
*/
/*
   \section WFS gradient pixel offset

   Use powfs.pixoffx and powfs.pixoffy to set WFs gradient pix offset

   - If abs(pixoffx) is less than 1. All subapertures have uniform pixel offset along x/r: pixoffx and along y/a: pixoffy (pixel).

   - If pixoffx==1 : There is a rotational offset with maximum value of pixoffy (pixel) at the edge.

   - If pixoffx==2: There is a global offset along x and y of pixoffy (pixel) that accounts for PCCCD geometry

*/
static void
setup_powfs_dtf(powfs_t* powfs, const parms_t* parms, int ipowfs){
	dmat* pixoffx=0;
	dmat* pixoffy=0;
	if(parms->powfs[ipowfs].pixoffx||parms->powfs[ipowfs].pixoffy){
		const int nsa=powfs[ipowfs].pts->nsa;
		if(fabs(parms->powfs[ipowfs].pixoffx)<1){
			info("powfs%d: uniform pixel offset\n", ipowfs);
			//both pixoff within 1 denotes constant offset in unit of pixel.
			pixoffx=dnew(nsa, 1); dset(pixoffx, parms->powfs[ipowfs].pixoffx);
			pixoffy=dnew(nsa, 1); dset(pixoffy, parms->powfs[ipowfs].pixoffy);
		} else{
			const int nwfs=1; //parms->powfs[ipowfs].nwfs;
			pixoffx=dnew(nsa, nwfs);
			pixoffy=dnew(nsa, nwfs);
			if(parms->powfs[ipowfs].pixoffx==1){
				info("powfs%d: CCD rotation wrt lenslet @ %g pixel.\n", ipowfs, parms->powfs[ipowfs].pixoffy);
				//pixoffx==1 denotes clocking effect with pixoffy denotes maximum amount in unit of pixel CCW.
				real pixtheta=0.5*(parms->powfs[ipowfs].pixtheta+parms->powfs[ipowfs].radpixtheta);
				const real dsah=powfs[ipowfs].pts->dsa*0.5;
				const real angle=parms->powfs[ipowfs].pixoffy*pixtheta/(parms->aper.d*0.5-dsah);
				const real ct=cos(angle);
				const real st=sin(angle);

				for(int jwfs=0; jwfs<nwfs; jwfs++){
					for(int isa=0; isa<nsa; isa++){
						real ox=powfs[ipowfs].saloc->locx[isa]+dsah;
						real oy=powfs[ipowfs].saloc->locy[isa]+dsah;
						real dx=ox*ct-oy*st-ox;
						real dy=ox*st+oy*ct-oy;
						if(parms->powfs[ipowfs].radpix){//radial coordinate.
							dy=sqrt(dx*dx+dy*dy);
							dx=0;
						}
						P(pixoffx, isa, jwfs)=dx/pixtheta;
						P(pixoffy, isa, jwfs)=dy/pixtheta;
					}
				}

			} else if(parms->powfs[ipowfs].pixoffx==2){
				info("powfs%d: CCD global shift wrt lenslet @ %g pixel.\n", ipowfs, parms->powfs[ipowfs].pixoffy);
				for(int jwfs=0; jwfs<nwfs; jwfs++){
					for(int isa=0; isa<nsa; isa++){
						real gx=parms->powfs[ipowfs].pixoffy;
						real gy=0;
						if(parms->powfs[ipowfs].radpix){
							real angle=PR(PR(powfs[ipowfs].srot, jwfs, 0), isa, 0);
							real ct=cos(angle);
							real st=sin(angle);
							real gx2=gx*ct+gy*st;//XY->RA: CW
							gy=-gx*st+gy*ct;
							gx=gx2;
						}
						P(pixoffx, isa, jwfs)=gx;
						P(pixoffy, isa, jwfs)=gy;
					}
				}
			} else{
				error("Invalid input in pixoffx\n");
			}
		}
		if(parms->save.setup>1){
			writebin(pixoffx, "powfs%d_pixoffx", ipowfs);
			writebin(pixoffy, "powfs%d_pixoffy", ipowfs);
		}
	}

	powfs[ipowfs].pixoffx=pixoffx;
	powfs[ipowfs].pixoffy=pixoffy;

	powfs[ipowfs].dtf=mkdtf(parms->powfs[ipowfs].wvl,
		powfs[ipowfs].pts->dsa,
		parms->powfs[ipowfs].embfac,
		powfs[ipowfs].notfx,
		powfs[ipowfs].notfy,
		powfs[ipowfs].pixpsax,
		powfs[ipowfs].pixpsay,
		parms->powfs[ipowfs].radpixtheta,
		parms->powfs[ipowfs].pixtheta,
		powfs[ipowfs].pixoffx,
		powfs[ipowfs].pixoffy,
		parms->powfs[ipowfs].pixblur,
		parms->powfs[ipowfs].radpix?powfs[ipowfs].srot:NULL
	);

	int nwvl=parms->powfs[ipowfs].nwvl;
	powfs[ipowfs].dtheta=dnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		P(powfs[ipowfs].dtheta, iwvl)=powfs[ipowfs].dtf[iwvl].dtheta;
		if(parms->save.setup>1){
			writebin(powfs[ipowfs].dtf[iwvl].nominal,
				"powfs%d_dtf%d_nominal", ipowfs, iwvl);
			writebin(powfs[ipowfs].dtf[iwvl].si,
				"powfs%d_dtf%d_si", ipowfs, iwvl);
		}
	}
}
/**
   setup the range to sodium layer as an additional parameter.
*/
static void setup_powfs_focus(powfs_t* powfs, const parms_t* parms, int ipowfs){
	if(!parms->powfs[ipowfs].llt||!parms->powfs[ipowfs].llt->fnrange) return;
	char* fnrange=parms->powfs[ipowfs].llt->fnrange;
	warning("loading sodium range variation from %s\n", fnrange);
	if(powfs[ipowfs].focus) dfree(powfs[ipowfs].focus);
	powfs[ipowfs].focus=dread("%s", fnrange);
	if(NY(powfs[ipowfs].focus)!=1&&NY(powfs[ipowfs].focus)!=parms->powfs[ipowfs].nwfs){
		error("fnrange has wrong format. Must be column vectors of 1 or %d columns\n",
			parms->powfs[ipowfs].nwfs);
	}
	/*(D/h)^2/(16*sqrt(3)) convert from range to WFE in m but here we want
	  focus mode, so just do 1/(2*h^2).*/
	/*1./cos() is for zenith angle adjustment of the range.*/
	real range2focus=0.5*pow(1./parms->powfs[ipowfs].hs, 2)*(1./cos(parms->sim.za));
	dscale(powfs[ipowfs].focus, range2focus);
	if(parms->save.setup){
		writebin(powfs[ipowfs].focus, "powfs%d_focus", ipowfs);
	}
}

/**
   Load and smooth out sodium profile and adjust for telescope altitude and
   zenith angle. We preserve the sum of the sodium profile, which represents a
   scaling factor of the signal level.

*/
static void setup_powfs_sodium(powfs_t* powfs, const parms_t* parms, int ipowfs){
	char* fnprof=parms->powfs[ipowfs].llt->fnprof;
	dcell* Nains=dcellread("%s", fnprof);
	int nprof=PN(Nains);
	if(nprof!=1&&nprof!=parms->powfs[ipowfs].nwfs){
		error("The sodium profile input %s is in wrong fromat\n", fnprof);
	}
	real dxnew=0;
	if(parms->dbg.na_smooth){/*resampling the sodium profile by binning. */
		/*Make new sampling: */
		const real rsamax=dmax(powfs[ipowfs].srsamax);
		const real dthetamin=dmin(powfs[ipowfs].dtheta);
		/*minimum sampling required. */
		dxnew=pow(parms->powfs[ipowfs].hs, 2)/rsamax*dthetamin;
	}
	powfs[ipowfs].sodium=smooth_cell(Nains, dxnew);
	dcellfree(Nains);
	if(parms->save.setup){
		writebin(powfs[ipowfs].sodium, "powfs%d_sodium", ipowfs);
	}
	
	//profile used for i0 and matched filter (if set)
	fnprof=parms->powfs[ipowfs].llt->fnprep;
	if(fnprof){
		Nains=dcellread("%s", fnprof);
		nprof=PN(Nains);
		if(nprof!=1&&nprof!=parms->powfs[ipowfs].nwfs){
			error("The sodium profile input %s is in wrong fromat\n", fnprof);
		}
		powfs[ipowfs].sodiumprep=smooth_cell(Nains, dxnew);
		dcellfree(Nains);

		if(parms->save.setup){
			writebin(powfs[ipowfs].sodiumprep, "powfs%d_sodiumprep", ipowfs);
		}
	}
}
static int etf_match(etf_t* etf, int icol, real hs, real thresh){
	return (etf&&etf->icol==icol&&fabs(etf->hs-hs)<thresh);
}
/**
   Compute Elongation Transfer function.
   - mode=0: for preparation of matched filter.
   - mode=1: for simulation.
   - mode=2: for simulation, next profile (linear interpolation.)
*/
void setup_powfs_etf(powfs_t *powfs, const parms_t *parms, int ipowfs, int mode, int icol, real deltah, real thresh){
	if(!parms->powfs[ipowfs].llt) return;
    etf_t** petf=0;
	const cell* sodium=CELL(powfs[ipowfs].sodium);
	if(mode==0){/*preparation. */
		if(powfs[ipowfs].etfprep&&powfs[ipowfs].etfsim!=powfs[ipowfs].etfprep){
			etf_free(powfs[ipowfs].etfprep);
		}
		if(powfs[ipowfs].sodiumprep){
			dbg("Using sodiumprep for etfprep\n");
			sodium=CELL(powfs[ipowfs].sodiumprep);
			petf=&powfs[ipowfs].etfprep;
		} else if(etf_match(powfs[ipowfs].etfsim, icol, parms->powfs[ipowfs].hs+deltah, thresh)){
			powfs[ipowfs].etfprep=powfs[ipowfs].etfsim;
			dbg("Using etfsim as etfprep\n");
		}else{
			petf=&powfs[ipowfs].etfprep;
		}
	} else if(mode==1){/*first pair for interpolation*/
		if(etf_match(powfs[ipowfs].etfsim, icol, parms->powfs[ipowfs].hs+deltah, thresh)){
			dbg("No need to update etfsim\n");
		} else{
			if(powfs[ipowfs].etfsim!=powfs[ipowfs].etfprep){
				etf_free(powfs[ipowfs].etfsim);
			}
			powfs[ipowfs].etfsim=0;
			if(etf_match(powfs[ipowfs].etfsim2, icol, parms->powfs[ipowfs].hs+deltah, thresh)){
				//reuse etfsim2 as etfsim
				dbg("reuse etfsim2 as etfsim.\n");
				powfs[ipowfs].etfsim=powfs[ipowfs].etfsim2;
			} else{
				petf=&powfs[ipowfs].etfsim;
			}
		}
	} else if(mode==2){/*second pair for interpolation*/
		if(etf_match(powfs[ipowfs].etfsim2, icol, parms->powfs[ipowfs].hs+deltah, thresh)){
			dbg("No need to update etfsim2\n");
		} else{
			if(powfs[ipowfs].etfsim2!=powfs[ipowfs].etfsim){
				etf_free(powfs[ipowfs].etfsim2);
			}
			powfs[ipowfs].etfsim2=0;
			petf=&powfs[ipowfs].etfsim2;
		}
	} else{
		error("Invalid mode=%d\n", mode);
	}
	if(petf){
		dbg("mketf: powfs %d using column %d with dh=%g\n", ipowfs, icol, deltah);
		*petf=mketf(powfs[ipowfs].dtf, sodium, icol, 
					powfs[ipowfs].srot, powfs[ipowfs].srsa,
					parms->powfs[ipowfs].hs+deltah, 
					parms->sim.htel, parms->sim.za,	!parms->dbg.na_interp);
	}
}

/**
   setting up uplink pts/amp lotf
*/
static void
setup_powfs_llt(powfs_t* powfs, const parms_t* parms, int ipowfs){
	if(!parms->powfs[ipowfs].llt) return;
	const int nwvl=parms->powfs[ipowfs].nwvl;
	real wvl0=P(parms->powfs[ipowfs].wvl, 0);
	llt_t* llt=powfs[ipowfs].llt=mycalloc(1, llt_t);
	const llt_cfg_t* lltcfg=parms->powfs[ipowfs].llt;

	real lltd=lltcfg->d;
	const real dx=powfs[ipowfs].pts->dx;//so that OTF has the same sampling
	real lltdsa=MAX(lltd, powfs[ipowfs].pts->dsa);
	int nxsa=round(lltdsa/dx); lltdsa=dx*nxsa;
	pts_t* lpts=llt->pts=ptsnew(1, lltdsa, lltdsa, nxsa, nxsa, dx, dx);
	llt->pts->origx[0]=(dx-lltdsa)*0.5;
	llt->pts->origy[0]=(dx-lltdsa)*0.5;

	real sumamp2=0;
	llt->amp=dnew(nxsa*nxsa, 1);
	if(lltcfg->fnamp){
		map_t* lltamp=mapread("%s", lltcfg->fnamp);
		prop_grid_pts(lltamp, llt->pts, P(llt->amp), 1, 0, 0, 1, 1, 0, 0);
		sumamp2=dinn(llt->amp, llt->amp);
		mapfree(lltamp);
	} else{
		real* amps=P(llt->amp);
		const real l2max=pow(lltd*0.5, 2);
		/*the waist is defined as the radius where amplitude
		  drop to 1/e or intensity to 1/e^2.*/
		real r2waist=pow(lltd*0.5*parms->powfs[ipowfs].llt->widthp, 2);
		//misreg is treated as pupil centering error. It clips the laser
		const real misregx=P(parms->powfs[ipowfs].llt->misreg, 0);
		const real misregy=P(parms->powfs[ipowfs].llt->misreg, 1);
		const real ox=llt->pts->origx[0];
		const real oy=llt->pts->origy[0];
		for(int iy=0; iy<nxsa; iy++){
			real yy=iy*dx+oy;
			for(int ix=0; ix<nxsa; ix++){
				real xx=ix*dx+ox;
				real r2=xx*xx+yy*yy;
				real r2m=pow(xx-misregx, 2)+pow(yy-misregy, 2);
				real amp=exp(-r2m/r2waist);
				sumamp2+=pow(amp, 2);
				if(r2<=l2max){
					amps[iy*nxsa+ix]=amp;
				}
			}
		}
	}
	/*normalized so that max(otf)=1 for unclipped beam; */
	sumamp2=1./(sqrt(sumamp2));
	dscale(llt->amp, sumamp2);
	if(lltcfg->fnsurf){
		mapcell* ncpa=genscreen_str(lltcfg->fnsurf);
		int nlotf=NX(ncpa)*NY(ncpa);
		assert(nlotf==1||nlotf==parms->powfs[ipowfs].nwfs);
		llt->ncpa=dcellnew(nlotf, 1);
		for(int ilotf=0; ilotf<nlotf; ilotf++){
			P(llt->ncpa, ilotf)=dnew(nxsa, nxsa);
			prop_grid_pts(P(ncpa, ilotf), llt->pts, P(P(llt->ncpa, ilotf)), 1, 0, 0, 1, 0, 0, 0);
		}
		cellfree(ncpa);
	}
	/*find negative values in llt->amp and transfer then to surface error.*/
	for(int i=0; i<NX(llt->amp)*NY(llt->amp); i++){
		if(P(llt->amp, i)<0){
			if(!llt->ncpa){
				llt->ncpa=dcellnew(1, 1);
				P(llt->ncpa, 0)=dnew(nxsa, nxsa);
			}
			for(int ic=0; ic<NX(llt->ncpa); ic++){
				if(nwvl>1) error("Please implement\n");
				P(P(llt->ncpa, ic), i)+=wvl0/2;
			}
			P(llt->amp, i)=-P(llt->amp, i);
		}
	}

	llt->loc=mksqloc(nxsa, nxsa, dx, dx, lpts->origx[0], lpts->origy[0]);
	llt->mcc=pts_mcc_ptt(llt->pts, P(llt->amp));
	llt->imcc=dcellinvspd_each(llt->mcc);
	if(lltcfg->focus){
		if(!llt->ncpa){
			llt->ncpa=dcellnew(1, 1);
			P(llt->ncpa, 0)=dnew(nxsa, nxsa);
		}
		real nm2rad=1e-9*2*sqrt(3)*pow(lltcfg->d, -2);
		dmat* focus=dnew(nxsa, nxsa);
		loc_add_focus(focus, llt->loc, lltcfg->focus*nm2rad);
		real var=0, piston=0;
		long count=0;
		for(long i=0; i<llt->loc->nloc; i++){
			if(P(llt->amp, i)>0){
				count++;
				var+=P(focus, i)*P(focus, i);
				piston+=P(focus, i);
			}
		}
		var/=count;
		piston/=count;
		dadds(focus, -piston);
		var=sqrt(var-piston*piston);
		for(int ic=0; ic<NX(llt->ncpa); ic++){
			dadd(&P(llt->ncpa, ic), 1, focus, lltcfg->focus*1e-9/var);
		}
		cellfree(focus);
		info("Adding focus %g nm (unweighted) to LLT ncpa\n", lltcfg->focus);
	}
	/*Remove tip/tilt/focus from ncpa */
	if(llt->ncpa&&parms->powfs[ipowfs].llt->ttfr){
		int nmod=parms->powfs[ipowfs].llt->ttfr==2?4:3;
		dmat* pttfa=zernike(llt->loc, parms->powfs[ipowfs].llt->d, 0, 2, 0);
		dmat* pttf=drefcols(pttfa, 0, nmod);
		dmat* proj=dpinv(pttf, CELL(llt->amp));
		dmat* res=dnew(nmod, 1);
		for(int ilotf=0; ilotf<NX(llt->ncpa)*NY(llt->ncpa); ilotf++){
			dzero(res);
			dmulvec(P(res), proj, P(P(llt->ncpa, ilotf)), 1);
			dmulvec(P(P(llt->ncpa, ilotf)), pttf, P(res), -1);
		}
		if(parms->save.setup){
			writebin(pttf, "powfs%d_llt_pttf", ipowfs);
			writebin(proj, "powfs%d_llt_proj", ipowfs);
		}
		dfree(pttfa);
		dfree(pttf);
		dfree(proj);
		dfree(res);
	}
	if(parms->save.setup){
		locwrite(llt->loc, "powfs%d_llt_loc", ipowfs);
		writebin(llt->amp, "powfs%d_llt_amp", ipowfs);
		writebin(llt->imcc, "powfs%d_llt_imcc", ipowfs);
		if(llt->ncpa){
			writebin(llt->ncpa, "powfs%d_llt_ncpa", ipowfs);
		}
	}
	if(parms->save.setup>1){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			if(powfs[ipowfs].etfprep&&powfs[ipowfs].etfprep[iwvl].etf){
				writebin(powfs[ipowfs].etfprep[iwvl].etf,
					"powfs%d_etfprep%d", ipowfs, iwvl);
			}
		}
		if(powfs[ipowfs].etfsim&&powfs[ipowfs].etfsim!=powfs[ipowfs].etfprep){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				if(powfs[ipowfs].etfsim[iwvl].etf){
					writebin(powfs[ipowfs].etfsim[iwvl].etf,
						"powfs%d_etfsim%d", ipowfs, iwvl);
				}
			}
		}
	}
}
/**
   Setup CoG NEA for reconstruction.
*/
static void
setup_powfs_cog_nea(const parms_t* parms, powfs_t* powfs, int ipowfs){
	TIC;tic;
	const int nwfs=parms->powfs[ipowfs].nwfs;
	const int nsa=powfs[ipowfs].saloc->nloc;
	const int ntry=500;
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const real pixthetax=parms->powfs[ipowfs].radpixtheta;
	const real pixthetay=parms->powfs[ipowfs].pixtheta;
	const real rne=parms->powfs[ipowfs].rne;
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
	const real cogthres=parms->powfs[ipowfs].cogthres;
	const real cogoff=parms->powfs[ipowfs].cogoff;
	intstat_t* intstat=powfs[ipowfs].intstat;
	rand_t rstat; seed_rand(&rstat, 1);
	if(!intstat||!intstat->i0){
		if(!parms->powfs[ipowfs].phyusenea){
			error("powfs[%d].i0 is not available, please enable phyusenea.\n", ipowfs);
		}
	}
	info("Compute NEA in CoG using Monte Carlo simulation with %d trials\n", ntry);
	dcellfree(powfs[ipowfs].sanea);
	powfs[ipowfs].sanea=dcellnew(NY(intstat->i0), 1);

	for(int jwfs=0; jwfs<nwfs; jwfs++){
		//int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
		if(jwfs==0||(intstat&&NY(intstat->i0)>1)){
			real* srot=NULL;
			if(parms->powfs[ipowfs].radpix){
				srot=P(PR(powfs[ipowfs].srot, jwfs, 0));
			}
			dmat** bkgrnd2=NULL;
			dmat** bkgrnd2c=NULL;
			dmat* psanea=P(powfs[ipowfs].sanea, jwfs)=dnew(nsa, 3);

			if(powfs[ipowfs].bkgrnd){
				bkgrnd2=PCOLR(powfs[ipowfs].bkgrnd, jwfs);
			}
			if(powfs[ipowfs].bkgrndc){
				bkgrnd2c=PCOLR(powfs[ipowfs].bkgrndc, jwfs);
			}
			dmat* nea=dnew(2, 2);
			for(int isa=0; isa<nsa; isa++){
				dmat* bkgrnd2i=bkgrnd2?bkgrnd2[isa]:NULL;
				dmat* bkgrnd2ic=bkgrnd2c?bkgrnd2c[isa]:NULL;
				dzero(nea);
				dmat* ints=P(intstat->i0, isa, jwfs);/*equivalent noise*/
				cog_nea(P(nea), ints, cogthres, cogoff, ntry, &rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
				P(nea, 0)=P(nea, 0)*pixthetax*pixthetax;
				P(nea, 3)=P(nea, 3)*pixthetay*pixthetay;
				P(nea, 1)=P(nea, 1)*pixthetax*pixthetay;
				P(nea, 2)=P(nea, 1);

				if(srot){
					dmat* nea2=0;
					drotvecnn(&nea2, nea, srot[isa]);
					dfree(nea); nea=nea2; nea2=0;
				}
				P(psanea, isa, 0)=P(nea, 0);//xx
				P(psanea, isa, 1)=P(nea, 3);//yy
				P(psanea, isa, 2)=P(nea, 1);//xy
			}
			dfree(nea);
		}
	}//for jwfs

	toc2("setup_powfs_cog_nea");
}


/**
   Setup the (matched filter or CoG) pixel processing parameters for physical optics wfs.
*/
static void
setup_powfs_phygrad(powfs_t* powfs, const parms_t* parms, int ipowfs){
	long nsa=powfs[ipowfs].saloc->nloc;
	if(powfs[ipowfs].intstat){
		error("Should only be called once\n");
	}
	if(parms->powfs[ipowfs].phytype_recon==PTYPE_MF||parms->powfs[ipowfs].phytype_sim==PTYPE_MF||!parms->powfs[ipowfs].phyusenea
		||(powfs[ipowfs].opdbias&&parms->powfs[ipowfs].ncpa_method==NCPA_I0)
		||parms->powfs[ipowfs].phytype_sim==PTYPE_CORR
		){
		intstat_t* intstat=powfs[ipowfs].intstat=mycalloc(1, intstat_t);
		if(parms->powfs[ipowfs].i0load){
			info("Loading i0, gx, gy\n");
			if(zfexist("%s/powfs%d_i0", parms->powfs[ipowfs].i0load, ipowfs)){
				intstat->i0=dcellread("%s/powfs%d_i0", parms->powfs[ipowfs].i0load, ipowfs);
			} else{
			//convert runtime saved ints.
				intstat->i0=dcellnew(nsa, parms->powfs[ipowfs].nwfs);
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dcell* i0=dcellread("%s/ints_1_wfs%d", parms->powfs[ipowfs].i0load, iwfs);
					for(int isa=0; isa<nsa; isa++){
						P(intstat->i0, isa, jwfs)=dref(P(i0, isa));
					}
					dcellfree(i0);
				}
			}
			if((parms->powfs[ipowfs].phytype_recon==PTYPE_MF||parms->powfs[ipowfs].phytype_sim==PTYPE_MF)
				&&!parms->powfs[ipowfs].mtchfft){
				intstat->gx=dcellread("%s/powfs%d_gx", parms->powfs[ipowfs].i0load, ipowfs);
				intstat->gy=dcellread("%s/powfs%d_gy", parms->powfs[ipowfs].i0load, ipowfs);
			}
			if(intstat->i0->keywords){
				real dt=search_keyword_num(intstat->i0->keywords, "dt");
				if(!isinf(dt)){
					real ratio=parms->sim.dt*parms->powfs[ipowfs].dtrat/dt;
					info("Scale loaded i0/gx/gy by %g\n", ratio);
					dcellscale(intstat->i0, ratio);
					dcellscale(intstat->gx, ratio);
					dcellscale(intstat->gy, ratio);
				}
			}
		} else{
			if(parms->powfs[ipowfs].piinfile){
			/*load psf. 1 for each wavefront sensor. */
				info("Using 1 sepsf for each wfs when loading sepsf\n");
				intstat->sepsf=dccellnew(parms->powfs[ipowfs].nwfs, 1);
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dcell* sepsf=P(intstat->sepsf, jwfs)=dcellread("%s_wfs%d", parms->powfs[ipowfs].piinfile, iwfs);
					real pmax=dmax(P(sepsf, 0));
					if(P(P(sepsf, 0), 0)>pmax*0.5){
						error("wfs %d:  psf must have peak at center, not corner.\n", iwfs);
					}
					if(NX(sepsf)!=nsa){
						error("piinfile doesn't match\n");
					}
				}
			} else if(parms->powfs[ipowfs].lo&&parms->powfs[ipowfs].order<=2){
				error("Please specify piinfile for lo order phy wfs\n");
			} else{
				/*Gen short exposure OTFs due to atmosphere*/
				cccell* otf=0;
				cccell* lotf=0;

				dcell* opdbias=parms->powfs[ipowfs].ncpa_method==NCPA_I0?powfs[ipowfs].opdbias:0;
				otf=genseotf(powfs[ipowfs].pts, CELL(powfs[ipowfs].realamp),
					opdbias, CELL(powfs[ipowfs].realsaa), parms->powfs[ipowfs].wvl,
					parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
					parms->powfs[ipowfs].embfac);
				const int print_psf=2; //1: uplink 2: downlink
				if(parms->powfs[ipowfs].llt){
					lotf=genseotf(powfs[ipowfs].llt->pts, CELL(powfs[ipowfs].llt->amp),
						powfs[ipowfs].llt->ncpa, NULL, parms->powfs[ipowfs].wvl,
						parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
						parms->powfs[ipowfs].embfac);

					if(print_psf==1 || 1){
						//print uplink PSF width
						dccell* lltpsf=0;
						const dmat* wvl=parms->powfs[ipowfs].wvl;
						gensepsf(&lltpsf, lotf, NULL, NULL, wvl, 0, 0);
						for(int iwvl=0; iwvl<PN(wvl); iwvl++){
							for(int illt=0; illt<NX(lotf); illt++){
								dmat* psf=P(P(lltpsf, illt, iwvl), 0);
								const real dpsf=P(wvl, 0)/(NX(psf)*powfs[ipowfs].llt->pts->dx)*RAD2AS;
								real fwhm=dfwhm_gauss(psf)*dpsf;
								info("Uplink FWHM (illt %d, iwvl %d) is %g\"\n", illt, iwvl, fwhm);
							}
						}
						if(parms->save.setup){//Save uplink PSF.
							writebin(lltpsf, "powfs%d_llt_psf", ipowfs);
						}
						cellfree(lltpsf);
					}
				}

				/*Generating short exposure psfs for both uplink and downlink
				turbulence effect. */
				
				gensepsf(&intstat->sepsf, otf, lotf, CELL(powfs[ipowfs].realsaa),
					parms->powfs[ipowfs].wvl, powfs[ipowfs].notfx, powfs[ipowfs].notfy);
				if(print_psf==2){
					const dmat *wvl=parms->powfs[ipowfs].wvl;
					for(int iwvl=0; iwvl<PN(wvl); iwvl++){
						for(int illt=0; illt<NX(lotf); illt++){
							dcell *psfs=P(intstat->sepsf, illt, iwvl);
							for(int isa=0; isa<PN(psfs); isa++){
								if(P(P(powfs[ipowfs].realsaa, illt), isa)>0.95){
									dmat *psf=P(psfs, isa);
									const real dpsf=P(wvl, 0)/(NX(psf)*powfs[ipowfs].llt->pts->dx)*RAD2AS;
									real fwhm=dfwhm_gauss(psf)*dpsf;
									info("Downlink FWHM (illt %d, iwvl %d) is %g\"\n", illt, iwvl, fwhm);
									break;
								}
							}
						}
					}
				}

				if(parms->save.setup>1&&intstat){
					writebin(P(intstat->sepsf, 0), "powfs%d_sepsf", ipowfs);
				}
				/*Free short exposure otf. */
				ccellfree(lotf);
				cellfree(otf);
			}
			/*generate short exposure i0,gx,gy from psf. */
			{	
				if(parms->powfs[ipowfs].llt){
					setup_powfs_etf(powfs, parms, ipowfs, 0, parms->powfs[ipowfs].llt->colprep, 0, 0);
				}
				cccell** pfotf=(parms->powfs[ipowfs].phytype_sim==PTYPE_MAP
					||(parms->dbg.wfslinearity!=-1&&parms->wfs[parms->dbg.wfslinearity].powfs==ipowfs))?&intstat->fotf:0;
				gensei(&intstat->i0, &intstat->gx, &intstat->gy, pfotf,
					intstat->sepsf, powfs[ipowfs].dtf, powfs[ipowfs].etfprep, powfs[ipowfs].realsaa,
					parms->powfs[ipowfs].radgx?powfs[ipowfs].srot:NULL,
					parms->powfs[ipowfs].siglevs, parms->powfs[ipowfs].wvlwts, NULL,
					parms->powfs[ipowfs].i0scale, parms->powfs[ipowfs].mtchstc);
				//gensei(parms, powfs, ipowfs);
			}
			if(parms->save.setup){
				writebin(intstat->i0, "powfs%d_i0", ipowfs);
				writebin(intstat->gx, "powfs%d_gx", ipowfs);
				writebin(intstat->gy, "powfs%d_gy", ipowfs);
			}
			{
				//Test sodium fitting
				//First, recreate sepsf without opdbias
				int SODIUM_FIT=0;
				READ_ENV_INT(SODIUM_FIT, 0, 1);
				if(SODIUM_FIT){
					sodium_fit_wrap(NULL, NULL, &intstat->i0, &intstat->gx, &intstat->gy, intstat->i0, parms, powfs, ipowfs, 
					parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 1, 0);
				}
			}
		}

	}
	/*Generating Matched filter */
	if(parms->powfs[ipowfs].phytype_recon==PTYPE_MF||parms->powfs[ipowfs].phytype_sim==PTYPE_MF){
		genmtch(parms, powfs, ipowfs);
		if(parms->save.setup){
			writebin(powfs[ipowfs].intstat->mtche, "powfs%d_mtche", ipowfs);
		}
	}
	if(parms->powfs[ipowfs].phytype_recon==PTYPE_COG&&parms->powfs[ipowfs].skip!=3&&!parms->powfs[ipowfs].phyusenea){
		setup_powfs_cog_nea(parms, powfs, ipowfs);
	}
	if(parms->save.setup){
		writebin(powfs[ipowfs].sanea, "powfs%d_sanea", ipowfs);
	}
}
/*
  Setup gradient offset for calibration. opdadd is the wavefront aberration in
  WFS due to optics, without DM correction. opdbias is the wavefront aberration
  in WFS after DM system flat is applied.
*/
void setup_powfs_calib(const parms_t* parms, powfs_t* powfs){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	//if opdadd is null, but dm_ncpa is not, there will still be opdbias.
		if(powfs[ipowfs].opdbias){
			//Always compute gradncpa. May not use it in CMF non updated case.
			if(!powfs[ipowfs].gradncpa){
				powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs, 1);
			} else{
				warning("gradncpa already exists\n");
				dcellzero(powfs[ipowfs].gradncpa);
			}
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				if(P(powfs[ipowfs].opdbias, jwfs)){
					real* realamp=P(P(powfs[ipowfs].realamp, jwfs));
					if(parms->powfs[ipowfs].type==WFS_PY){//pywfs
						dmat* ints=0;
						pywfs_ints(&ints, powfs[ipowfs].pywfs, P(powfs[ipowfs].opdbias, jwfs), parms->powfs[ipowfs].siglev);
						//writebin(P(powfs[ipowfs].opdbias,jwfs), "opdbias\n");exit(0);
						pywfs_grad(&P(powfs[ipowfs].gradncpa, jwfs), powfs[ipowfs].pywfs, ints);
						dfree(ints);
					} else if(parms->powfs[ipowfs].gtype_sim==GTYPE_Z){//Ztilt
						pts_ztilt(&P(powfs[ipowfs].gradncpa, jwfs), powfs[ipowfs].pts,
							PR(powfs[ipowfs].saimcc, jwfs, 0),
							realamp, P(P(powfs[ipowfs].opdbias, jwfs)));
					} else{//Gtilt
						if(parms->powfs[ipowfs].ncpa_method==NCPA_G 
							|| (parms->powfs[ipowfs].ncpa_method==NCPA_I0&&parms->powfs[ipowfs].dither)){//GS0*opd
							dspmm(&P(powfs[ipowfs].gradncpa, jwfs), PR(powfs[ipowfs].GS0, jwfs, 0),
								P(powfs[ipowfs].opdbias, jwfs), "nn", 1);
							//MF drift control moved to wfsgrad.c
							//need gradncpa for sodium profile fit.
						}
					}
				}
			}
		}
		if(powfs[ipowfs].pixoffx){
			if(!powfs[ipowfs].gradncpa){
				powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs, 1);
			}
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int nsa=powfs[ipowfs].pts->nsa;
				if(!P(powfs[ipowfs].gradncpa, jwfs)){
					P(powfs[ipowfs].gradncpa, jwfs)=dnew(nsa*2, 1);
				}

				real pixthetax=parms->powfs[ipowfs].pixtheta;
				real pixthetay=parms->powfs[ipowfs].radpixtheta;
				for(int isa=0; isa<nsa; isa++){
					real gx=pixthetax*PR(powfs[ipowfs].pixoffx, isa, jwfs);
					real gy=pixthetay*PR(powfs[ipowfs].pixoffy, isa, jwfs);
					if(parms->powfs[ipowfs].radpix){
						real angle=PR(PR(powfs[ipowfs].srot, jwfs, 0), isa, 1);
						real ct=cos(angle);
						real st=sin(angle);
						real gx2=gx*ct-gy*st;//RA->XY; CCW
						gy=gx*st+gy*ct;
						gx=gx2;
					}
					P(P(powfs[ipowfs].gradncpa, jwfs), isa)+=gx;
					P(P(powfs[ipowfs].gradncpa, jwfs), isa+nsa)+=gy;
				}
			}
		}
		if(parms->save.setup){
			writebin(powfs[ipowfs].gradncpa, "powfs%d_gradncpa", ipowfs);
		}
	}
}

/**
   Setup the powfs struct based on parms and aper. Everything about wfs are
   setup here.  \callgraph */
powfs_t* setup_powfs_init(const parms_t* parms, aper_t* aper){
	TIC;tic;
	powfs_t* powfs=mycalloc(parms->npowfs, powfs_t);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs==0) continue;
		if(parms->powfs[ipowfs].type==WFS_SH){
			info2("\n%sSetting up powfs %d geom%s\n\n", GREEN, ipowfs, BLACK);
			setup_shwfs_geom(powfs, parms, aper, ipowfs);
			setup_shwfs_grad(powfs, parms, ipowfs);
		} else if(parms->powfs[ipowfs].type==WFS_PY){
			info2("\n%sSetting up powfs %d in Pyramid mode%s\n\n", GREEN, ipowfs, BLACK);
			pywfs_setup(parms->powfs[ipowfs].pycfg, powfs, parms, aper, ipowfs);
		} else{
			error("powfs %d: invalid wfstype=%d\n", ipowfs, parms->powfs[ipowfs].type);
		}
	}
	toc2("setup_powfs_init");
	return powfs;
}
/**
   Setup physical optics parameters for SHWFS, such as DTF, ETF, LLT, pixel processing.
*/
void setup_powfs_phy(const parms_t* parms, powfs_t* powfs){
	TIC;tic;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs&&parms->powfs[ipowfs].type==WFS_SH
			&&(parms->powfs[ipowfs].usephy
			   ||parms->powfs[ipowfs].psfout
			   ||parms->powfs[ipowfs].pistatout
			   ||parms->powfs[ipowfs].neaphy)){
			info2("\n%sSetting up powfs %d PO WFS%s\n\n", GREEN, ipowfs, BLACK);
			/*We have physical optics. setup necessary struct */
			setup_powfs_prep_phy(powfs, parms, ipowfs);
			setup_powfs_dtf(powfs, parms, ipowfs);
			if(parms->powfs[ipowfs].llt){
				//load and smooth sodium profile
				setup_powfs_sodium(powfs, parms, ipowfs);/*read sodium profile and smooth it */
				//etf for first time step. etfprep is moved to before i0 generation
				setup_powfs_etf(powfs, parms, ipowfs, 1, parms->powfs[ipowfs].llt->colsim, 0., 0.);
				setup_powfs_llt(powfs, parms, ipowfs);
			}

			if(parms->powfs[ipowfs].llt){
				/*If there is LLT, setup the extra focus term if needed. */
				setup_powfs_focus(powfs, parms, ipowfs);
			}
			if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy){
				setup_powfs_phygrad(powfs, parms, ipowfs);
			}
		}
	}/*ipowfs */
	toc2("setup_powfs_phy");
}
/**
   free unused parameters before simulation starts
*/
void free_powfs_unused(const parms_t* parms, powfs_t* powfs){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(powfs[ipowfs].intstat){
			cellfree(powfs[ipowfs].intstat->sepsf);
			cellfree(powfs[ipowfs].intstat->gx);
			cellfree(powfs[ipowfs].intstat->gy);
		}
		if(!parms->powfs[ipowfs].needGS0&&powfs[ipowfs].GS0){
			dspcellfree(powfs[ipowfs].GS0);
			powfs[ipowfs].GS0=NULL;
		}
	}
}

static void free_powfs_shwfs(powfs_t* powfs, int ipowfs){
	dtf_free(powfs[ipowfs].dtf);
	dspcellfree(powfs[ipowfs].GS0);
	dcellfree(powfs[ipowfs].neasim);
	dcellfree(powfs[ipowfs].sanea);
	if(powfs[ipowfs].intstat){
		intstat_t* intstat=powfs[ipowfs].intstat;
		cellfree(intstat->fotf);
		cellfree(intstat->mtche);
		cellfree(powfs[ipowfs].intstat->i0);
		dfree(intstat->i0sum);
		dfree(intstat->i0sumsum);
		free(intstat);
		powfs[ipowfs].intstat=NULL;
	}

	dcellfree(powfs[ipowfs].srot);
	dcellfree(powfs[ipowfs].srsa);
	dfree(powfs[ipowfs].srsamax);
	cellfree(powfs[ipowfs].sprint);
	dfree(powfs[ipowfs].pixoffx);
	dfree(powfs[ipowfs].pixoffy);

	cellfree(powfs[ipowfs].saimcc);
	if(powfs[ipowfs].llt){
		ptsfree(powfs[ipowfs].llt->pts);
		dfree(powfs[ipowfs].llt->amp);
		locfree(powfs[ipowfs].llt->loc);
		dcellfree(powfs[ipowfs].llt->mcc);
		dcellfree(powfs[ipowfs].llt->imcc);
		dcellfree(powfs[ipowfs].llt->ncpa);
		free(powfs[ipowfs].llt);
	}
	dcellfree(powfs[ipowfs].sodium);
	dcellfree(powfs[ipowfs].sodiumprep);
	if(powfs[ipowfs].etfprep!=powfs[ipowfs].etfsim){
		etf_free(powfs[ipowfs].etfprep);
	}
	etf_free(powfs[ipowfs].etfsim);
	etf_free(powfs[ipowfs].etfsim2);
	dcellfree(powfs[ipowfs].opdadd);
	dcellfree(powfs[ipowfs].opdbias);
	dcellfree(powfs[ipowfs].gradncpa);
	dfree(powfs[ipowfs].dtheta);
}
/**
   Free all parameters of powfs at the end of simulation.
*/
void free_powfs(const parms_t* parms, powfs_t* powfs){
	free_powfs_unused(parms, powfs);
	free_powfs_fit(powfs, parms);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		free_powfs_geom(powfs, ipowfs);
		free_powfs_shwfs(powfs, ipowfs);
		pywfs_free(powfs[ipowfs].pywfs);
	}
	free(powfs);
}
