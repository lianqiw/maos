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

#include "common.h"
#include "setup_powfs.h"
#include "mtch.h"
#include "genseotf.h"
#include "pywfs.h"

#define MOVES(p,i,j) p[i]=p[j]
#define MOVED(p,n,i,j) memcpy(p+n*i, p+n*j, sizeof(double)*n)
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
	MOVED(loc->locx,nxsa,count,isa);	\
	MOVED(loc->locy,nxsa,count,isa);	\
    }

/**
   \file maos/setup_powfs.c

   Setting up WFS geometry. like the subaperture location, subaperture grid
   points, physical optics detection transfer function, LGS elongation transfer
   function, etc.

   \todo isolate DTF, ETF, routines, make then generic interface and relocate to the lib folder.

   Do not use sparse interpolation to replace ray tracing for fine sampled
   destination grid, especially cubic splines. The interpolation marix takes too
   much space.  

   TODO: This routine and POWFS_T should only contain information about the
   simulation, not about any model used during reconstruction (RTC) to avoid
   leaking information from the "real world (simulation)" to our knowledge (RTC).
*/



/**
   Free the powfs geometric parameters
*/
static void
free_powfs_geom(POWFS_T *powfs, int ipowfs){
    if(!powfs[ipowfs].pts){
	return;
    }
    ptsfree(powfs[ipowfs].pts);
    dfree(powfs[ipowfs].saa);
    locfree(powfs[ipowfs].saloc);
    locfree(powfs[ipowfs].loc);
    dfree(powfs[ipowfs].amp);
    if(powfs[ipowfs].loc_tel){
	cellfree(powfs[ipowfs].loc_tel);
	dcellfree(powfs[ipowfs].amp_tel);
	dcellfree(powfs[ipowfs].saa_tel);
    }
    if(powfs[ipowfs].loc_dm){
	cellfree(powfs[ipowfs].loc_dm);
    }
    dcellfree(powfs[ipowfs].realamp);
    dcellfree(powfs[ipowfs].realsaa);
    dfree(powfs[ipowfs].sumamp);
    dfree(powfs[ipowfs].sumamp2);
    locfft_free(powfs[ipowfs].fieldstop);
}
/**
   Convert amplitude map of subapertures to normalized subaperture area. Used only in setup_powfs_geom*/
static dmat *wfsamp2saa(dmat *wfsamp, long nxsa){
    double areanormfactor=1./(double)nxsa;
    long nsa=wfsamp->nx/nxsa;
    dmat* saa=dnew(nsa, 1);
    for(int isa=0; isa<nsa; isa++){
	double *amp=wfsamp->p+isa*nxsa;
	double area=0;
	for(int iamp=0; iamp<nxsa; iamp++){
	    area+=amp[iamp]*amp[iamp];
	}
	saa->p[isa]=area*areanormfactor;
    }
    return saa;
}

void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, int iwfs){
    long nloc=loc->nloc;
    dmat *ampmask=dnew(nloc, 1);
    double ht=parms->atm.hmax*0.7;
    for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
	int jpowfs=parms->wfs[jwfs].powfs;
	if(parms->powfs[jpowfs].lo) continue;
	double hs=parms->wfs[iwfs].hs;
	double r=parms->aper.d*0.5*(1.-ht/hs)/(1.-ht/hs);
	double sx=(parms->wfs[jwfs].thetax-parms->wfs[iwfs].thetax)*ht;
	double sy=(parms->wfs[jwfs].thetay-parms->wfs[iwfs].thetay)*ht;
	loccircle(ampmask->p, loc, sx,sy, r, 1);
    }
    for(int i=0; i<nloc; i++){
	if(ampmask->p[i]<0.5) amp->p[i]=0;
    }
    dfree(ampmask);
}
static void
sa_reduce(POWFS_T *powfs, int ipowfs, double thresarea, int reduce_isolated){
    dmat *saa=NULL;
    if(powfs[ipowfs].saa_tel){
	warning_once("Todo: Improve to allow different sa for same wfs type\n");
	const double scale=1./(double)powfs[ipowfs].saa_tel->nx;
	for(int i=0; i<powfs[ipowfs].saa_tel->nx; i++){
	    dadd(&saa, 1, powfs[ipowfs].saa_tel->p[i], scale);
	}
    }else{
	saa=ddup(powfs[ipowfs].saa);
    }
    if(dmax(saa)>1.01){
	warning("The sa area maxes to %g, which should be leq 1 (misregistration can cause this).\n", 
		dmax(saa));
    }
    if(reduce_isolated){
	//Disable subapertures that doesn't have at least two close neighbors on its side
	loc_t *ptsloc=(loc_t*)powfs[ipowfs].pts;
	loc_create_map(ptsloc);
	double dx1=1./ptsloc->dx;
	double dy1=1./ptsloc->dy;
	int disabled=0;
	do{
	    disabled=0;
	    for(long isa=0; isa<ptsloc->nloc; isa++){
		if(saa->p[isa]>thresarea){
		    long ix=round((ptsloc->locx[isa]-ptsloc->map->ox)*dx1);
		    long iy=round((ptsloc->locy[isa]-ptsloc->map->oy)*dy1);
		    int ncount=0;
		    for(int jx=-1; jx<2; jx++){
			for(int jy=-1; jy<2; jy++){
			    if(abs(jx+jy)==1){
				long jsa=loc_map_get(ptsloc->map, ix+jx, iy+jy);
				if(jsa && saa->p[jsa-1]>thresarea){
				    ncount++;
				}
			    }
			}
		    }
		    if(ncount<=1){
			saa->p[isa]=0;
			disabled++;
		    }
		}
	    }
	}while(disabled);
	loc_free_map(ptsloc);
    }
    int count=0;
    const int nxsa=powfs[ipowfs].pts->nx * powfs[ipowfs].pts->nx;
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	if(saa->p[isa]>=thresarea){
	    /*Area is above threshold, keep.  Shift pts, ptsm, loc, locm, amp,
	      ampm, saloc area is already normalized that maxes to 1. The MOVE*
	      are defined in the beginining of this file.*/
	    if(count!=isa){
		MOVEPTS(powfs[ipowfs].pts,  count, isa);
		MOVES(powfs[ipowfs].saa->p, count, isa);
		MOVELOC(powfs[ipowfs].saloc, count, isa);
		MOVEDLOC(powfs[ipowfs].loc, nxsa, count, isa);
		MOVED(powfs[ipowfs].amp->p, nxsa, count, isa);
		if(powfs[ipowfs].saa_tel){
		    for(int jwfs=0; jwfs<powfs[ipowfs].saa_tel->nx; jwfs++){
			MOVES(powfs[ipowfs].saa_tel->p[jwfs]->p, count, isa);
			MOVEDLOC(powfs[ipowfs].loc_tel->p[jwfs], nxsa, count, isa);
			MOVED(powfs[ipowfs].amp_tel->p[jwfs]->p,nxsa, count, isa);
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
    locresize(powfs[ipowfs].loc, count*nxsa);
    dresize(powfs[ipowfs].saa, count, 1);
    dresize(powfs[ipowfs].amp, count*nxsa, 1);
    powfs[ipowfs].saasum=dsum(powfs[ipowfs].saa);
    if(powfs[ipowfs].loc_tel){
	for(int jwfs=0; jwfs<powfs[ipowfs].saa_tel->nx; jwfs++){
	    locresize(powfs[ipowfs].loc_tel->p[jwfs], count*nxsa);
	    dresize(powfs[ipowfs].amp_tel->p[jwfs], count*nxsa, 1);
	    dresize(powfs[ipowfs].saa_tel->p[jwfs], count, 1);
	}
    }	
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
setup_powfs_geom(POWFS_T *powfs, const PARMS_T *parms, 
		 APER_T *aper, int ipowfs){
    free_powfs_geom(powfs,ipowfs);
    int nwfsp=parms->powfs[ipowfs].nwfs;
    /*order of the system. 60 for TMT */
    const int order  = parms->powfs[ipowfs].order;
    /*Subaperture lateral length. 0.5 for TMT. */
    const double dsa = parms->powfs[ipowfs].dsa;
    /*The subaperture area that are considered as full. */
    double areafulli;
    if(order>2){
	areafulli=1.;
    }else{
	areafulli=4./M_PI;
    }

    /*Number of OPD pixels in 1 direction in each
      subaperture. Make it even to do fft.*/
    int nx = 2*(int)round(0.5*dsa/parms->powfs[ipowfs].dx);
    const double dx=dsa/nx;/*adjust dx. */
    const double dxoffset=dx*0.5;//Always keep points inside subaperture for simulation.
    if(fabs(parms->powfs[ipowfs].dx-dx)>EPS)
	warning("Adjusting dx from %g to %g\n", parms->powfs[ipowfs].dx,dx);
    if(fabs(dsa - nx * dx)>EPS){
	warning("nx=%d,dsa=%f,dx=%f not agree\n", nx, dsa, dx);
    }
    info2("There are %d points in each subaperture of %gm.\n", nx, dsa);
    const int nxsa=nx*nx;/*Total Number of OPD points. */
    if(parms->powfs[ipowfs].saloc){
	powfs[ipowfs].saloc=locread("%s", parms->powfs[ipowfs].saloc);
	if(fabs(powfs[ipowfs].saloc->dx-dsa)>1e-6*fabs(dsa)){
	    error("loaded saloc has dx=%g, while powfs.dsa=%g\n",
		  powfs[ipowfs].saloc->dx, dsa);
	}
    }else{
	/*The coordinate of the subaperture (lower left coordinate) */
	powfs[ipowfs].saloc=locnew(order*order, dsa, dsa);
	int count = 0;
	/*Offset of the coordinate of the center most subaperture from the center. */
	double offset;
	if(order & 1){/*odd */
	    offset = -0.5;
	}else{
	    offset = 0.0;
	}
	/*r2max: Maximum distance^2 from the center to keep a subaperture */
	double r2max=pow(order*0.5, 2);
	double r2min=dsa<parms->aper.din?pow(parms->aper.din/dsa/2,2):-1;
	/*the lower left *grid* coordinate of the subaperture */

	/*Collect all the subapertures that are within the allowed radius*/
	for(int j=-order/2; j<=(order-1)/2; j++){
	    for(int i=-order/2; i<=(order-1)/2; i++){
		//Normalized coordinate in uniq of sa size
		double xc=((double)i+offset);
		double yc=((double)j+offset);
		//Radius of four corners.
		double r1=pow(xc+1,2)+pow(yc+1,2);
		double r2=pow(xc+1,2)+pow(yc,2);
		double r3=pow(xc,2)+pow(yc+1,2);
		double r4=pow(xc,2)+pow(yc,2);
		if(r1<r2max || r2<r2max || r3<r2max || r4<r2max||
		   r1>r2min || r2>r2min || r3>r2min || r4>r2min){
		    powfs[ipowfs].saloc->locx[count]=xc*dsa;
		    powfs[ipowfs].saloc->locy[count]=yc*dsa;
		    count++;
		}
	    }
	}
	powfs[ipowfs].saloc->nloc=count;
    }
    /*convert saloc to pts*/
    powfs[ipowfs].pts=ptsnew(powfs[ipowfs].saloc->nloc, dsa, dsa, nx, dx, dx);
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	powfs[ipowfs].pts->origx[isa]=powfs[ipowfs].saloc->locx[isa]+dxoffset;
	powfs[ipowfs].pts->origy[isa]=powfs[ipowfs].saloc->locy[isa]+dxoffset;
    }
    /*Calculate the amplitude for each subaperture OPD point by ray tracing from
      the pupil amplitude map. Pupil distortion is accounted for.*/
    powfs[ipowfs].loc=pts2loc(powfs[ipowfs].pts);
    /*The assumed amp. */
    powfs[ipowfs].amp=mkamp(powfs[ipowfs].loc, aper->ampground, 
			    -parms->misreg.pupil->p[0],-parms->misreg.pupil->p[1], 
			    parms->aper.d, parms->aper.din);
    /*The threashold for normalized area (by areafulli) to keep subaperture. */
    double thresarea=parms->powfs[ipowfs].saat;
    dmat *ampi=NULL;
    if(parms->powfs[ipowfs].safill2d<1){
	/*subaperture amplitude map to simulate lenslet fill factor*/
	int nedge=1;
	while ((nx-2*nedge)*(nx-2*nedge)>nx*nx*parms->powfs[ipowfs].safill2d){
	    nedge++;
	}
	ampi=dnew(nx,nx);
	dmat*  pampi=ampi/*PDMAT*/;
	double alpha=(nx*nx*parms->powfs[ipowfs].safill2d-(nx-2*nedge)*(nx-2*nedge))
	    /((nx-2*(nedge-1))*(nx-2*(nedge-1))-(nx-2*nedge)*(nx-2*nedge));
	double tot=0;
	for(int iy=nedge-1; iy<nx-nedge+1; iy++){
	    for(int ix=nedge-1; ix<nx-nedge+1; ix++){
		IND(pampi,ix,iy)=1;
		if(ix==nedge-1 || ix==nx-nedge || iy==nedge-1 || iy==nx-nedge){
		    IND(pampi,ix,iy)=alpha;
		}
		tot+=IND(pampi,ix,iy);
	    }
	}
	if(parms->save.setup){
	    writebin(ampi,"powfs%d_ampi", ipowfs);
	}
	for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	    for(int i=0; i<nx*nx; i++){
		powfs[ipowfs].amp->p[nx*nx*isa+i]*=ampi->p[i];
	    }
	}
	/*do not multiply to siglev. Already handled automatically*/
	thresarea*=parms->powfs[ipowfs].safill2d;
    }

    powfs[ipowfs].saa=wfsamp2saa(powfs[ipowfs].amp, nxsa);

   
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
    if(!parms->powfs[ipowfs].saloc){
	sa_reduce(powfs, ipowfs, thresarea, 0);
    }
    setup_powfs_misreg_dm(powfs, parms, aper, ipowfs);
    powfs[ipowfs].realsaa=dcellnew(nwfsp, 1);
    for(int jwfs=0; jwfs<nwfsp; jwfs++){
	if(powfs[ipowfs].loc_tel){
	    powfs[ipowfs].realsaa->p[jwfs]=dref(powfs[ipowfs].saa_tel->p[jwfs]);
	}else{
	    powfs[ipowfs].realsaa->p[jwfs]=dref(powfs[ipowfs].saa);
	}
    }
    powfs[ipowfs].realamp=dcellnew(nwfsp, 1);
    powfs[ipowfs].sumamp=dnew(nwfsp, 1);
    powfs[ipowfs].sumamp2=dnew(nwfsp, 1);
    for(int jwfs=0; jwfs<nwfsp; jwfs++){
	dmat *realamp;
	if(powfs[ipowfs].loc_tel){
	    realamp=powfs[ipowfs].amp_tel->p[jwfs];
	}else{
	    realamp=powfs[ipowfs].amp;
	}
	double sumamp2=0;
	double sumamp=0;
	for(long i=0; i<powfs[ipowfs].loc->nloc; i++){
	    sumamp2+=realamp->p[i]*realamp->p[i];
	    sumamp+=realamp->p[i];
	}
	powfs[ipowfs].sumamp2->p[jwfs]=sumamp2;
	powfs[ipowfs].sumamp->p[jwfs]=sumamp;
	powfs[ipowfs].realamp->p[jwfs]=dref(realamp);
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
	locwrite((loc_t*)powfs[ipowfs].pts,"powfs%d_pts",ipowfs);
	locwrite(powfs[ipowfs].saloc,"powfs%d_saloc",ipowfs); 
	writebin(powfs[ipowfs].saa,"powfs%d_saa",ipowfs);
	locwrite(powfs[ipowfs].loc,"powfs%d_loc",ipowfs);
	writebin(powfs[ipowfs].amp,"powfs%d_amp",ipowfs);
	if(powfs[ipowfs].loc_tel){
	    writebin(powfs[ipowfs].saa_tel,"powfs%d_saa_tel",ipowfs);
	    writebin(powfs[ipowfs].amp_tel,"powfs%d_amp_tel",ipowfs);
	    writebin(powfs[ipowfs].loc_tel,"powfs%d_loc_tel",ipowfs);
	}
	if(powfs[ipowfs].loc_dm){
	    for(int idm=0; idm<parms->ndm; idm++){
		for(int jwfs=0; jwfs<nwfsp; jwfs++){
		    writebin(powfs[ipowfs].loc_dm,"powfs%d_loc_dm", ipowfs);
		}
	    }
	}
    }
}
void 
setup_powfs_misreg_tel(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    int nwfsp=parms->powfs[ipowfs].nwfs;
    if(parms->misreg.tel2wfs){
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
	if(parms->powfs[ipowfs].type==0){
	    powfs[ipowfs].saa_tel=dcellnew(nwfsp, 1);
	}
#pragma omp parallel for shared(isset)
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    if(parms->misreg.tel2wfs[iwfs]){
		isset=1;
		powfs[ipowfs].loc_tel->p[jwfs]
		    =loctransform(powfs[ipowfs].loc, parms->misreg.tel2wfs[iwfs]);
		powfs[ipowfs].amp_tel->p[jwfs]
		    =mkamp(powfs[ipowfs].loc_tel->p[jwfs], aper->ampground,
			   -parms->misreg.pupil->p[0], -parms->misreg.pupil->p[1], 
			   parms->aper.d, parms->aper.din);
		if(parms->powfs[ipowfs].type==0){
		    const int nxsa=powfs[ipowfs].pts->nx * powfs[ipowfs].pts->nx;
		    powfs[ipowfs].saa_tel->p[jwfs]=wfsamp2saa(powfs[ipowfs].amp_tel->p[jwfs], nxsa);
		}
	    }
	}
	if(!isset){
	    cellfree(powfs[ipowfs].loc_tel);
	    dcellfree(powfs[ipowfs].saa_tel); 
	    dcellfree(powfs[ipowfs].amp_tel);
	}else{
	    toc2("misreg.tel2wfs");
	}
    }/*if misreg */ 
}
void 
setup_powfs_misreg_dm(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    (void)aper;
    int nwfsp=parms->powfs[ipowfs].nwfs;
    if(parms->misreg.dm2wfs){
	TIC;tic;
	/*
	  Misregistration/distortion from DM to WFS pupil. Not about telescope
	  pupil. The distorted grid are used for ray tracing from DM to WFS.
	  They are not used for wavefront reconstruction
	*/
	powfs[ipowfs].loc_dm=(loccell*)cellnew(nwfsp*parms->ndm, 1);
	int isset=0;
#pragma omp parallel for collapse(2) shared(isset)
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int jwfs=0; jwfs<nwfsp; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		if(parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]){
		    powfs[ipowfs].loc_dm->p[jwfs+nwfsp*idm]
			=loctransform(powfs[ipowfs].loc, parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]);
		    isset=1;
		}
	    }
	}
	if(!isset){
	    cellfree(powfs[ipowfs].loc_dm);
	}else{
	    toc2("misreg.dm2wfs");
	}
    }
}
/**
   Creating geometric wavefront gradient operator GS0 and ZA0 from WFS OPD to
   subaperture grads. Use loc, and ampm. Do not locm. locm is only used to do
   ray tracing.  This is simulation, not RTC, so use real system distortion,
   misregistration, etc. */
static void 
setup_powfs_grad(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(parms->powfs[ipowfs].gtype_recon==0 ||parms->powfs[ipowfs].gtype_sim==0){
	dspcellfree(powfs[ipowfs].GS0);
	/*Setting up every gradient tilt (g-ztilt) */
	if(parms->load.GS0){
	    powfs[ipowfs].GS0=dspcellread("powfs%d_GS0",ipowfs);
	    if(powfs[ipowfs].amp_tel){
		assert(powfs[ipowfs].GS0->nx==parms->powfs[ipowfs].nwfs);
	    }else{
		assert(powfs[ipowfs].GS0->nx==1);
	    }
	}else{
	    /*This mkg takes about 5 seconds. */
	    if(powfs[ipowfs].amp_tel){
		powfs[ipowfs].GS0=dspcellnew(parms->powfs[ipowfs].nwfs, 1);
	    }else{
		powfs[ipowfs].GS0=dspcellnew(1, 1);
	    }
	    for(int iwfs=0; iwfs<powfs[ipowfs].GS0->nx; iwfs++){
		powfs[ipowfs].GS0->p[iwfs]=mkg(powfs[ipowfs].loc, 
					       powfs[ipowfs].loc,
					       powfs[ipowfs].realamp->p[iwfs],
					       powfs[ipowfs].saloc,
					       1, 0, 0, 1);
	    }
	}
    }
    if(parms->powfs[ipowfs].gtype_recon==1 ||parms->powfs[ipowfs].gtype_sim==1){
	/*setting up zernike best fit (ztilt) inv(M'*W*M). good for NGS. */
	if(parms->powfs[ipowfs].order>4) 
	    warning("Ztilt for high order wfs is not good");
	powfs[ipowfs].nsaimcc=MAX(1,(powfs[ipowfs].loc_tel?parms->powfs[ipowfs].nwfs:1));
	int nsaimcc=powfs[ipowfs].nsaimcc;
	cellfree(powfs[ipowfs].saimcc);
	powfs[ipowfs].saimcc=(dccell*)cellnew(nsaimcc, 1);
	for(int imcc=0; imcc<nsaimcc; imcc++){
	    dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].realamp->p[imcc]->p);
	    powfs[ipowfs].saimcc->p[imcc]=dcellinvspd_each(mcc);
	    dcellfree(mcc);
	}
    }
 
    if(!parms->powfs[ipowfs].neaphy){
	powfs[ipowfs].neasim=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    const long nsa=powfs[ipowfs].saloc->nloc;
	    dmat *nea=NULL;
	    if(parms->powfs[ipowfs].neasimfile){
		if(parms->powfs[ipowfs].neasim!=-1){
		    error("neasimfile and neasim can not both be supplied\n");
		}
		nea=dread("%s_wfs%d",parms->powfs[ipowfs].neasimfile,iwfs);/*rad */
		if(nea->nx!=nsa || (nea->ny!=2 && nea->ny!=3)){
		    error("wfs %d: NEA read from %s_wfs%d has incorrect format. We want %ldx2(3) array.\n",
			  iwfs, parms->powfs[ipowfs].neasimfile, iwfs,powfs[ipowfs].saloc->nloc);
		}
	    }else{
		nea=dnew(nsa, 2);
		if(parms->powfs[ipowfs].neasim<0){
		    dset(nea, parms->powfs[ipowfs].nearecon);
		}else{
		    dset(nea, parms->powfs[ipowfs].neasim);
		}
		dscale(nea, 1./206265000);/*convert from mas to rad. */
	    }
	    if(nea->ny==2){
		dresize(nea,nsa,3);
		memset(nea->p+nsa*2,0,nsa*sizeof(double));
	    }
	    /*Sanity check */
	    double neamax=dmax(nea);
	    if(neamax>0){
		if(neamax<parms->powfs[ipowfs].pixtheta*1e-5){
		    warning("wfs %d: NEA(%g mas)<pixtheta(%g mas)*1e-6. Unit error?\n",
			    iwfs, neamax*206265000, parms->powfs[ipowfs].pixtheta*206265000);
		}else if(neamax>parms->powfs[ipowfs].pixtheta){
		    warning("wfs %d: NEA(%g mas)>pixtheta(%g mas). Unit error?\n",
			    iwfs, neamax*206265000, parms->powfs[ipowfs].pixtheta*206265000);
		}
	    }
	    /*Scale by dtrat */
	    dscale(nea, 1./sqrt(parms->powfs[ipowfs].dtrat));
	    /*Scale by normalized subaperture area. */
	    double *saa=powfs[ipowfs].realsaa->p[jwfs]->p;
	    for(long isa=0; isa<nsa; isa++){
		/*scale nea by sqrt(1/area). (seeing limited) */
		/*scale nea by 1/area if diffraction limited (NGS) */
		const int dl=0; /*Assume not diffraction limited. minor point*/
		double scale=dl?(1./saa[isa]):(1./sqrt(saa[isa]));
		nea->p[isa]*=scale;
		nea->p[isa+nsa]*=scale;
		nea->p[isa+nsa*2]*=scale;/*cross term. */
	    }
	    powfs[ipowfs].neasim->p[jwfs]=nea;
	}
	if(parms->save.setup){
	    writebin(powfs[ipowfs].neasim,"powfs%d_neasim", ipowfs);
	}
    }
    if(parms->save.setup && powfs[ipowfs].GS0){
	writebin(powfs[ipowfs].GS0,"powfs%d_GS0",ipowfs);
    }
}
/**
   Prepare the parameters for physical optics setup.

   \todo Make LGS number of pixels depend on subapertue elongation in radial
   coordinate CCD mode. Make ncompx, ncompy dependent on number of pixels for radial ccd.
   

*/
static void 
setup_powfs_prep_phy(POWFS_T *powfs,const PARMS_T *parms,int ipowfs){
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int pixpsay=parms->powfs[ipowfs].pixpsa;
    const int radpix=parms->powfs[ipowfs].radpix;
    const double dsa=powfs[ipowfs].pts->dsa;
    const int nsa=powfs[ipowfs].saloc->nloc;
    if(parms->powfs[ipowfs].llt){
	const int nllt=parms->powfs[ipowfs].llt->n;
	double rsa2, rsa2max=0;
	dcellfree(powfs[ipowfs].srot);
	dcellfree(powfs[ipowfs].srsa);
	dcellfree(powfs[ipowfs].sprint);

	powfs[ipowfs].srot=dcellnew(nllt,1);
	powfs[ipowfs].srsa=dcellnew(nllt,1);
	powfs[ipowfs].srsamax=dnew(nllt,1);
	powfs[ipowfs].sprint=dcellnew(nllt,1);

	for(int illt=0;illt<nllt;illt++){
	    /*adjusted llt center because pts->orig is corner */
	    double ox2=parms->powfs[ipowfs].llt->ox->p[illt]-dsa*0.5;
	    double oy2=parms->powfs[ipowfs].llt->oy->p[illt]-dsa*0.5;
	    powfs[ipowfs].srot->p[illt]=dnew(nsa,1);
	    powfs[ipowfs].srsa->p[illt]=dnew(nsa,1);

	    for(int isa=0; isa<nsa; isa++){
		double ddx=(powfs[ipowfs].saloc->locx[isa]-ox2);
		double ddy=(powfs[ipowfs].saloc->locy[isa]-oy2);
		rsa2=pow(ddx,2)+pow(ddy,2);
		powfs[ipowfs].srot->p[illt]->p[isa]=atan2(ddy,ddx);
		powfs[ipowfs].srsa->p[illt]->p[isa]=sqrt(rsa2);
		if(rsa2>rsa2max) rsa2max=rsa2;
	    }
	    powfs[ipowfs].srsamax->p[illt]=sqrt(rsa2max);
	    double dprint=parms->aper.d*0.1;
	    int pnsa=(int)ceil(sqrt(rsa2max)/dsa)+1;
	    
	    double prot[pnsa];
	    powfs[ipowfs].sprint->p[illt]=dnew(pnsa,1);
	    double *pp=powfs[ipowfs].sprint->p[illt]->p;
	    for(int ind=0; ind<pnsa; ind++){
		prot[ind]=INFINITY;
		pp[ind]=-1;
	    }
	    double desrot=0;
	    if(fabs(parms->powfs[ipowfs].llt->ox->p[illt])<dsa 
	       && fabs(parms->powfs[ipowfs].llt->oy->p[illt])<dsa){
		desrot=0;
	    }else{
		double ddx=(0-parms->powfs[ipowfs].llt->ox->p[illt]);
		double ddy=(0-parms->powfs[ipowfs].llt->oy->p[illt]);
		desrot=atan2(ddy,ddx);
	    }
	    for(int isa=0; isa<nsa; isa++){
		int ind=(int)round(powfs[ipowfs].srsa->p[illt]->p[isa]/dprint);
		double irot=fabs(powfs[ipowfs].srot->p[illt]->p[isa]-desrot);
		if(irot<prot[ind]){
		    prot[ind]=irot;
		    pp[ind]=(double)isa;
		}
	    }
	}
	if(parms->save.setup){
	    writebin(powfs[ipowfs].sprint,"powfs%d_sprint",ipowfs);
	    writebin(powfs[ipowfs].srot,"powfs%d_srot",ipowfs);
	    writebin(powfs[ipowfs].srsa,"powfs%d_srsa",ipowfs);
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

    /*to avoid aliasing, fft warping. usually 2. */
    const double embfac=parms->powfs[ipowfs].embfac;
 
    /*size required to form detector image. */
    int ncompx, ncompy;
    if(parms->powfs[ipowfs].ncomp){
	ncompx=parms->powfs[ipowfs].ncomp;
	ncompy=parms->powfs[ipowfs].ncomp;
	warning2("ncomp is specified in input file to %dx%d\n", ncompx,ncompy);
    }else{
	/*
	  compute required otf size to cover the detector FoV
	  Need to revise this part: when ncomp is less than the
	  size of the full psf, may need padding.  study aliasing
	  extensively with full, cropped psf and detector
	  transfer function.
	*/
	double wvlmin=parms->powfs[ipowfs].wvl->p[0];
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(wvlmin>parms->powfs[ipowfs].wvl->p[iwvl])
		wvlmin=parms->powfs[ipowfs].wvl->p[iwvl];
	}
	double dtheta=wvlmin/(dsa*embfac);/*Min PSF sampling. */
	ncompx=ceil(pixpsax*pixthetax/dtheta);
	ncompy=ceil(pixpsay*pixthetay/dtheta);
	
	/*
	  Found that: Must set ncompx==ncompy even for rotationg either psf or
	  otf. reduce aliasing and scattering of image intensities.
	*/
	ncompx=ncompy=MAX(ncompx, ncompy);
	/*A few manual optimizations. */
	if(ncompx==ncompy){
	    if(ncompx>8 && ncompx<16){
		ncompx=ncompy=16;
	    }else if(ncompx>16 && ncompx<32){
		ncompx=ncompy=32;
	    }else if(ncompx>64 && ncompx<67){
		ncompx=ncompy=64;
	    }else if(ncompx<128 && ncompx>120){
		ncompx=ncompy=128;
	    }
	}
	ncompx=ncompy=nextfftsize(ncompx);
	info2("Subaperture DTF is %dx%d\n", ncompx,ncompy);
    }/*ncomp */
    powfs[ipowfs].ncompx=ncompx;
    powfs[ipowfs].ncompy=ncompy;

    if(parms->powfs[ipowfs].bkgrndfn){
	char *fn=parms->powfs[ipowfs].bkgrndfn;
	info2("Loading sky background/rayleigh backscatter from %s\n",fn);
	dcellfree(powfs[ipowfs].bkgrnd);
	powfs[ipowfs].bkgrnd=dcellread("%s",fn);
    }
    if(parms->powfs[ipowfs].bkgrndfnc){
	char *fn=parms->powfs[ipowfs].bkgrndfnc;
	info2("Loading sky background/rayleigh backscatter correction from %s\n",fn);
	dcellfree(powfs[ipowfs].bkgrndc);
	powfs[ipowfs].bkgrndc=dcellread("%s",fn);
    }
    if(parms->powfs[ipowfs].bkgrndfn || parms->powfs[ipowfs].bkgrndfnc){
	double bkscale = parms->sim.dt*800*parms->powfs[ipowfs].dtrat;
	if(fabs(bkscale-1)>1.e-20){
	    dcellscale(powfs[ipowfs].bkgrnd, bkscale);
	    dcellscale(powfs[ipowfs].bkgrndc, bkscale);
	    warning("Scaling bkgrnd by %g", bkscale);
	}
	if(parms->powfs[ipowfs].bkgrndfn && 
	   (powfs[ipowfs].bkgrnd->nx!=powfs[ipowfs].saloc->nloc
	    ||(powfs[ipowfs].bkgrnd->ny!=1
	       && powfs[ipowfs].bkgrnd->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrnd is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrnd->nx, powfs[ipowfs].bkgrnd->ny,
		  powfs[ipowfs].saloc->nloc, parms->powfs[ipowfs].nwfs);
	}
	if(parms->powfs[ipowfs].bkgrndfnc && 
	   (powfs[ipowfs].bkgrndc->nx!=powfs[ipowfs].saloc->nloc
	    ||(powfs[ipowfs].bkgrndc->ny!=1
	       && powfs[ipowfs].bkgrndc->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrndc is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrndc->nx, powfs[ipowfs].bkgrndc->ny,
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
static void 
setup_powfs_dtf(POWFS_T *powfs,const PARMS_T *parms,int ipowfs){
    powfs[ipowfs].dtf=mkdtf(parms->powfs[ipowfs].wvl, 
			    powfs[ipowfs].pts->dsa,
			    parms->powfs[ipowfs].embfac,
			    powfs[ipowfs].ncompx,
			    powfs[ipowfs].ncompy,
			    powfs[ipowfs].pixpsax,
			    powfs[ipowfs].pixpsay,
			    parms->powfs[ipowfs].radpixtheta,
			    parms->powfs[ipowfs].pixtheta,
			    parms->powfs[ipowfs].pixoffx,
			    parms->powfs[ipowfs].pixoffy,
			    parms->powfs[ipowfs].pixblur,
			    powfs[ipowfs].srot,
			    parms->powfs[ipowfs].radpix,
			    parms->powfs[ipowfs].radrot);
    if(parms->powfs[ipowfs].radrot){
	info2("Rotating PSF for Polar CCD\n");/*Used mainly for on-axis launch */
    }
 
    int nwvl=parms->powfs[ipowfs].nwvl;
    powfs[ipowfs].dtheta=dnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	powfs[ipowfs].dtheta->p[iwvl]=powfs[ipowfs].dtf[iwvl].dtheta;
	if(parms->save.setup>1){
	    writebin(powfs[ipowfs].dtf[iwvl].nominal,
		     "powfs%d_dtf%d_nominal",ipowfs,iwvl);
	    writebin(powfs[ipowfs].dtf[iwvl].si,
		     "powfs%d_dtf%d_si",ipowfs,iwvl);
	}
    }
}
/**
   setup the range to sodium layer as an additional parameter.
*/
static void setup_powfs_focus(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt || !parms->powfs[ipowfs].llt->fnrange) return;
    char *fnrange=parms->powfs[ipowfs].llt->fnrange;
    warning("loading sodium range from %s\n", fnrange);
    if(powfs[ipowfs].focus) dfree(powfs[ipowfs].focus);
    powfs[ipowfs].focus=dread("%s",fnrange);
    if(powfs[ipowfs].focus->ny!=1 &&powfs[ipowfs].focus->ny!=parms->powfs[ipowfs].nwfs){
	error("fnrange has wrong format. Must be column vectors of 1 or %d columns\n",
	      parms->powfs[ipowfs].nwfs);
    }
    /*(D/h)^2/(16*sqrt(3)) convert it from range to WFE in m but here we want
      focus mode, so just do 1/(2*h^2).*/
    /*1./cos() is for zenith angle adjustment of the range.*/
    double range2focus=0.5*pow(1./parms->powfs[ipowfs].hs,2)*(1./cos(parms->sim.za));
    dscale(powfs[ipowfs].focus, range2focus);
    writebin(powfs[ipowfs].focus, "powfs%d_focus", ipowfs);
}

/**
   Load and smooth out sodium profile. We preserve the sum of the sodium profile,
   which represents a scaling factor of the signal level.  */
static void setup_powfs_sodium(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    char *fnprof=parms->powfs[ipowfs].llt->fnprof;
    dcell *Nains=dcellread("%s",fnprof);
    int nprof=Nains->nx*Nains->ny;
    if(nprof!=1 && nprof!=parms->powfs[ipowfs].nwfs){
	error("The sodium profile input %s is in wrong fromat\n", fnprof);
    }
    powfs[ipowfs].sodium=dcellnew(nprof, 1);
    for(int i=0; i<Nains->nx*Nains->ny; i++){
	dmat *Nain=Nains->p[i];
	if(Nain->ny<2 || Nain->nx!=Nains->p[0]->nx){
	    error("The sodium profile input %s is in wrong fromat\n", fnprof);
	}
	if(parms->dbg.na_smooth){/*resampling the sodium profile by binning. */
	    /*Make new sampling: */
	    const double rsamax=dmax(powfs[ipowfs].srsamax);
	    const double dthetamin=dmin(powfs[ipowfs].dtheta);
	    /*minimum sampling required. */
	    const double dxnew=pow(parms->powfs[ipowfs].hs,2)/rsamax*dthetamin;
	    powfs[ipowfs].sodium->p[i]=smooth(Nain, dxnew);
	}else{
	    info2("Not smoothing sodium profile\n");
	    powfs[ipowfs].sodium->p[i]=dref(Nain);
	}
    }
    dcellfree(Nains);
    if(parms->save.setup){
	writebin(powfs[ipowfs].sodium,"powfs%d_sodium",ipowfs);
    }
}
typedef struct{
    DTF_T *dtfs;  /**<The dtfs*/
    double hs;    /**<Guide star focus range*/
    dcell *sodium;/**<The sodium profile. First column is coordinate.*/
    int icol;     /**<Which sodium profile to use*/
    int nwvl;     /**<Number of wavelength*/
    dcell *srot;  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
    dcell *srsa;  /**<Subaperture to LLT distance*/
    double za;    /**<Zenith angle*/
    int no_interp;/**<Use direct sum instead of interpolation + FFT. Slower */
    int free;     /**<Free this array after using?*/
}mketf_t;
ETF_T* mketf_wrap(mketf_t *data){
    ETF_T*result=mketf(data->dtfs, data->hs, data->sodium, data->icol, data->nwvl,
		       data->srot, data->srsa, data->za, data->no_interp);
    if(data->free) free(data);
    return result;
}
/**
   Compute Elongation Transfer function.
   - mode=0: for preparation.
   - mode=1: for simulation.
*/
void setup_powfs_etf(POWFS_T *powfs, const PARMS_T *parms, int ipowfs, int mode, int istep){
    if(!parms->powfs[ipowfs].llt) return;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    mketf_t etfdata={powfs[ipowfs].dtf, 
		     parms->powfs[ipowfs].hs,
		     powfs[ipowfs].sodium,
		     istep,
		     parms->powfs[ipowfs].nwvl, 
		     powfs[ipowfs].srot,
		     powfs[ipowfs].srsa,
		     parms->sim.za,
		     !parms->dbg.na_interp, 0};
    if(mode==0){/*preparation. */
	if(powfs[ipowfs].etfprep && powfs[ipowfs].etfsim!=powfs[ipowfs].etfprep){
	    etf_free(powfs[ipowfs].etfprep, nwvl);
	}
	powfs[ipowfs].etfprep=mketf_wrap(&etfdata);
    }else{/*simulation*/
	if(mode==1){/*first pair for interpolation*/
	    if(powfs[ipowfs].etfsim==powfs[ipowfs].etfprep){
		powfs[ipowfs].etfsim=0;
	    }
	    etf_free(powfs[ipowfs].etfsim, nwvl);
	    if(istep!=0 && powfs[ipowfs].etfsim2){//reuse etfsim2 as etfsim
		powfs[ipowfs].etfsim=powfs[ipowfs].etfsim2;
		powfs[ipowfs].etfsim2=0;
	    }else{
		powfs[ipowfs].etfsim=mketf_wrap(&etfdata);
	    }
	}else if(mode==2){/*second pair for interpolation*/
	    static pthread_t etfthread=0;
	    if(powfs[ipowfs].etfsim2){
		error("etsim2 should be null\n");
		etf_free(powfs[ipowfs].etfsim2, nwvl);
	    }
	    if(etfthread){//preparation already running in a thread
		pthread_join(etfthread, (void**)(void*)&powfs[ipowfs].etfsim2);
		etfthread=0;
	    }else{
		powfs[ipowfs].etfsim2=mketf_wrap(&etfdata);
	    }
	    //asynchronously preparing for next update.
	    //Copy data to heap so that they don't disappear during thread execution
	    mketf_t *etfdata2=mycalloc(1,mketf_t);//freed by mketf_wrap.
	    memcpy(etfdata2, &etfdata, sizeof(mketf_t));
	    etfdata2->icol++;
	    etfdata2->free=1;
	    if(pthread_create(&etfthread, NULL, (void*(*)(void *))mketf_wrap, etfdata2)){
		warning("Thread creation failed\n");
		free(etfdata2);
		etfthread=0;
	    }
	}else{
	    error("Invalid mode=%d\n", mode);
	}
    }
}

/**
   setting up uplink pts/amp lotf
*/
static void 
setup_powfs_llt(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    double wvl0=parms->powfs[ipowfs].wvl->p[0];
    LLT_T *llt=powfs[ipowfs].llt=mycalloc(1,LLT_T);
    LLT_CFG_T *lltcfg=parms->powfs[ipowfs].llt;
    pts_t *lpts=llt->pts=mycalloc(1,pts_t);
    lpts->nsa=1;
    double lltd=lltcfg->d;
    lpts->dsay=lpts->dsa=MAX(lltd, powfs[ipowfs].pts->dsa);
    int notf=MAX(powfs[ipowfs].ncompx, powfs[ipowfs].ncompy);
    /*The otf would be dx/lambda. Make it equal to embfac*pts->dsa/lambda/notf)*/
    const double dx=lpts->dx=parms->powfs[ipowfs].embfac*powfs[ipowfs].pts->dsa/notf;
    lpts->dy=lpts->dx;
    const int nx=lpts->nx=round(lpts->dsa/lpts->dx);
    lpts->dsay=lpts->dsa=lpts->dx*lpts->nx;
    
    lpts->origx=mycalloc(1,double);
    lpts->origy=mycalloc(1,double);

    double oy=lpts->origx[0]=(dx-lpts->dsa)*0.5;
    double ox=lpts->origy[0]=(dx-lpts->dsa)*0.5;
    double sumamp2=0;
    llt->amp=dnew(nx*nx, 1);
    if(lltcfg->fnamp){
	map_t *lltamp=mapread("%s", lltcfg->fnamp);
	prop_grid_pts(lltamp, llt->pts, llt->amp->p, 1, 0, 0, 1, 1, 0, 0);
	sumamp2=dinn(llt->amp, llt->amp);
	mapfree(lltamp);
    }else{
	double *amps=llt->amp->p;
	double l2max =pow(lltd*0.5,2);
	/*the waist is defined as the radius where amplitude
	  drop to 1/e or intensity to 1/e^2.*/
	double r2waist=pow(lltd*0.5*parms->powfs[ipowfs].llt->widthp,2);
	for(int iy=0; iy<nx; iy++){
	    double yy=iy*dx+oy;
	    yy*=yy;
	    for(int ix=0; ix<nx; ix++){
		double xx=ix*dx+ox;
		xx*=xx;
		double r2=xx+yy;
		if(r2<=l2max){
		    amps[iy*nx+ix]=exp(-r2/r2waist);
		    sumamp2+=pow(amps[iy*nx+ix],2);
		}
	    }
	}
    }
    /*normalized so that max(otf)=1; */
    sumamp2=1./(sqrt(sumamp2));
    dscale(llt->amp, sumamp2);
    if(parms->powfs[ipowfs].llt->fnsurf){
	mapcell *ncpa=mapcellread("%s",parms->powfs[ipowfs].llt->fnsurf);
	int nlotf=ncpa->nx*ncpa->ny;
	assert(nlotf==1 || nlotf==parms->powfs[ipowfs].nwfs);
	llt->ncpa=dcellnew(nlotf, 1);
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    llt->ncpa->p[ilotf]=dnew(nx,nx);
	    prop_grid_pts(ncpa->p[ilotf], llt->pts, llt->ncpa->p[ilotf]->p, 1, 0, 0, 1, 0, 0, 0);
	}
	cellfree(ncpa);
    }
    /*find negative values in llt->amp and transfer then to surface error.*/
    for(int i=0; i<llt->amp->nx*llt->amp->ny; i++){
	if(llt->amp->p[i]<0){
	    if(!llt->ncpa){
		llt->ncpa=dcellnew(1,1);
		llt->ncpa->p[0]=dnew(nx,nx);
	    }
	    for(int ic=0; ic<llt->ncpa->nx; ic++){
		if(nwvl>1) error("Please implement\n");
		llt->ncpa->p[ic]->p[i]+=wvl0/2;
	    }
	    llt->amp->p[i]=-llt->amp->p[i];
	}
    }
  
    llt->loc=mksqloc(nx,nx,dx,dx, lpts->origx[0], lpts->origy[0]);
    llt->mcc =pts_mcc_ptt(llt->pts, llt->amp->p);
    llt->imcc =dcellinvspd_each(llt->mcc);
    /*Remove tip/tilt/focus from ncpa */
    if(llt->ncpa && parms->powfs[ipowfs].llt->ttfr){
	int nmod=4;//pick first four modes
	dmat *pttfa=zernike(llt->loc, parms->powfs[ipowfs].llt->d, 0, 2, 0);
	dmat *pttf=drefcols(pttfa, 0, nmod);
	dmat *proj=dpinv(pttf, llt->amp);
	dmat *res=dnew(nmod, 1);
	for(int ilotf=0; ilotf<llt->ncpa->nx*llt->ncpa->ny; ilotf++){
	    dzero(res);
	    dmulvec(res->p, proj, llt->ncpa->p[ilotf]->p, 1);
	    dmulvec(llt->ncpa->p[ilotf]->p, pttf, res->p, -1);
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
	locwrite(llt->loc,"powfs%d_llt_loc",ipowfs);
	writebin(llt->amp,"powfs%d_llt_amp",ipowfs);
	writebin(llt->imcc,"powfs%d_llt_imcc",ipowfs);
	if(llt->ncpa){
	    writebin(llt->ncpa,"powfs%d_llt_ncpa", ipowfs);
	}
    }
    if(parms->save.setup>1){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(powfs[ipowfs].etfprep[iwvl].p1){
		writebin(powfs[ipowfs].etfprep[iwvl].p1, 
			 "powfs%d_etfprep%d_1d",ipowfs,iwvl);
	    }
	    if(powfs[ipowfs].etfprep[iwvl].p2){
		writebin(powfs[ipowfs].etfprep[iwvl].p2,
			   "powfs%d_etfprep%d_2d",ipowfs,iwvl);
	    }
	}
	if(powfs[ipowfs].etfsim != powfs[ipowfs].etfprep){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(powfs[ipowfs].etfsim[iwvl].p1){
		    writebin(powfs[ipowfs].etfsim[iwvl].p1, 
			       "powfs%d_etfsim%d_1d",ipowfs,iwvl);
		}
		if(powfs[ipowfs].etfsim[iwvl].p2){
		    writebin(powfs[ipowfs].etfsim[iwvl].p2,
			       "powfs%d_etfsim%d_2d",ipowfs,iwvl);
		}
	    }
	}
    }
}
/*compute cog NEA using Monte Carlo realizations of noise*/
static void cog_nea(double *nea, dmat *ints, double cogthres, double cogoff, int ntry, 
		    rand_t *rstat, double bkgrnd, double bkgrndc, dmat *bkgrnd2i, dmat *bkgrnd2ic, double rne
		    ){
    dmat *ints2=dnew(ints->nx, ints->ny);
    double gnf[2]={0,0};
    double gny[2]={0,0};
    dcog(gnf, ints, 0, 0, cogthres, cogoff, 0);
    seed_rand(rstat, 1);/*reset the seed each time to make dminsearch work.*/
    nea[0]=0; nea[1]=0; nea[2]=0; nea[3]=0;
    for(int i=0; i<ntry; i++){
	dcp(&ints2, ints);
	addnoise(ints2, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, 0, rne);
	dcog(gny, ints2, 0, 0, cogthres, cogoff, 0);
	double errx=gny[0]-gnf[0];
	double erry=gny[1]-gnf[1];
	nea[0]+=errx*errx;
	nea[1]+=errx*erry;
	nea[3]+=erry*erry;
    }
    dfree(ints2);
    double stry=1./ntry;
    nea[0]=nea[0]*stry;
    nea[3]=nea[3]*stry;
    nea[1]=nea[1]*stry;
    nea[2]=nea[1];
}
typedef struct {
    dmat *ints;
    double bkgrnd;
    double bkgrndc;
    dmat *bkgrnd2i;
    dmat *bkgrnd2ic;
    double rne;
    rand_t *rstat;
    int ntry;
}cogdata_t;
static double
cogfun(double *x, cogdata_t *info){
    double nea[4];
    cog_nea(nea, info->ints, x[0], x[1], info->ntry, info->rstat, info->bkgrnd, info->bkgrndc, info->bkgrnd2i, info->bkgrnd2ic, info->rne);
    return nea[0]+nea[3];
}
/**
   Setup CoG gradient offset for simulation and NEA for reconstruction.
*/
static void 
setup_powfs_cog(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    TIC;tic;
    const int nwfs=parms->powfs[ipowfs].nwfs;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int ntry=500;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
    const double bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    int do_nea=0;
    rand_t rstat;
    double neaspeckle2=0;
    dcell *sanea=NULL;
    if(!intstat || !intstat->i0){
	if(!parms->powfs[ipowfs].phyusenea){
	    error("powfs[%d].i0 is not available, please enable phyusenea.\n", ipowfs);
	}
	if(parms->powfs[ipowfs].cogthres<0 && parms->powfs[ipowfs].cogoff<0){
	     error("i0 is not available, please specify powfs.cogthres and powfs.cogoff to non-negative numbers.\n");
	}
    }
    if(parms->powfs[ipowfs].phytype==2 && parms->powfs[ipowfs].skip!=3 && !parms->powfs[ipowfs].phyusenea){/*need nea in reconstruction*/
	do_nea=1;
	powfs[ipowfs].saneaxy=dcellnew(nsa, intstat->i0->ny);
	sanea=dcellnew(intstat->i0->ny, 1);
	seed_rand(&rstat, 1);
	double neaspeckle=parms->powfs[ipowfs].neaspeckle/206265000.;
	if(neaspeckle>pixthetax){
	    error("parms->powfs[%d].neaspeckle=%g is bigger than pixel size\n",
		  ipowfs, neaspeckle);
	}
	if(neaspeckle>0){
	    warning2("powfs%d: Adding speckle noise of %.2f mas\n", ipowfs, neaspeckle*206265000);
	}
	neaspeckle2=pow(neaspeckle,2);
    }
    powfs[ipowfs].cogcoeff=dcellnew(nwfs,1);
    
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(iwfs==0 || (intstat && intstat->i0->ny>1)){
	    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	    double *srot=NULL;
	    if(parms->powfs[ipowfs].radpix){
		srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p;
	    }
	    dmat **bkgrnd2=NULL;
	    dmat **bkgrnd2c=NULL;
	    if(do_nea){
		sanea->p[iwfs]=dnew(nsa,2);
	    }
	    if(powfs[ipowfs].bkgrnd){
		if(powfs[ipowfs].bkgrnd->ny==1){
		    bkgrnd2=powfs[ipowfs].bkgrnd->p;
		}else{
		    bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*wfsind;
		}
	    }
	    if(powfs[ipowfs].bkgrndc){
		if(powfs[ipowfs].bkgrndc->ny==1){
		    bkgrnd2c=powfs[ipowfs].bkgrndc->p;
		}else{
		    bkgrnd2c=powfs[ipowfs].bkgrndc->p+nsa*wfsind;
		}
	    }

	    powfs[ipowfs].cogcoeff->p[iwfs]=dnew(2,nsa);
	    dmat*  cogcoeff=powfs[ipowfs].cogcoeff->p[iwfs]/*PDMAT*/;
	    for(int isa=0; isa<nsa; isa++){
		dmat *bkgrnd2i=bkgrnd2?bkgrnd2[isa]:NULL;
		dmat *bkgrnd2ic=bkgrnd2c?bkgrnd2c[isa]:NULL;

		if(parms->powfs[ipowfs].cogthres<0 && parms->powfs[ipowfs].cogoff<0){
		    /*automatically figure out the optimal thres/offset for each subaperture.*/
		    double coeff[2]={1., 1.};
		    dmat *ints=intstat->i0->p[isa+iwfs*nsa];/*equivalent noise*/
		    cogdata_t data={ints, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne, &rstat, 100};
		    double ftol=0.0001;/*tolerance of nea^2 in pixel*/
		    int ncall=dminsearch(coeff, 2, ftol, 5000, (dminsearch_fun)cogfun, &data);
		    IND(cogcoeff,0,isa)=coeff[0];
		    IND(cogcoeff,1,isa)=coeff[1];
		    info("isa %d:, ncall=%3d, coeff=%g %g\n", isa, ncall, coeff[0], coeff[1]);
		}else{
		    IND(cogcoeff,0,isa)=parms->powfs[ipowfs].cogthres;
		    IND(cogcoeff,1,isa)=parms->powfs[ipowfs].cogoff;
		}
		if(do_nea){
		    dmat *nea=dnew(2,2);
		    dmat *ints=intstat->i0->p[isa+iwfs*nsa];/*equivalent noise*/
		    cog_nea(nea->p, ints, IND(cogcoeff,0,isa), IND(cogcoeff,1,isa), ntry, &rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
		    nea->p[0]=nea->p[0]*pixthetax*pixthetax+neaspeckle2;
		    nea->p[3]=nea->p[3]*pixthetay*pixthetay+neaspeckle2;
		    nea->p[1]=nea->p[1]*pixthetax*pixthetay;
		    nea->p[2]=nea->p[1];

		    sanea->p[iwfs]->p[isa]=nea->p[0];
		    sanea->p[iwfs]->p[isa+nsa]=nea->p[3];
		    if(srot){
			drotvecnn(&powfs[ipowfs].saneaxy->p[isa+nsa*iwfs], nea, srot[isa]);
		    }else{
			powfs[ipowfs].saneaxy->p[isa+nsa*iwfs]=dref(nea);
		    }
		    dfree(nea);
		}
	    }
	}else{
	    powfs[ipowfs].cogcoeff->p[iwfs]=dref(powfs[ipowfs].cogcoeff->p[0]);
	}
    }
    if(parms->save.setup){
	if(sanea){
	    writebin(sanea,"powfs%d_sanea", ipowfs);
	}
	writebin(powfs[ipowfs].cogcoeff,"powfs%d_cogcoeff", ipowfs);
    }
    dcellfree(sanea);
    toc2("setup_powfs_cog");
}


/**
   Setup the matched filter pixel processing parameters for physical optics wfs.
*/
static void 
setup_powfs_phygrad(POWFS_T *powfs,const PARMS_T *parms, int ipowfs){
    int disable_save_save=disable_save;
    disable_save=0;//temporarily disable this feature.
    long nsa=powfs[ipowfs].saloc->nloc;
    if(powfs[ipowfs].intstat){
	error("Should only be called once\n");
    }
    if(parms->powfs[ipowfs].phytype==1 || parms->powfs[ipowfs].phytypesim==1 || !parms->powfs[ipowfs].phyusenea){
	INTSTAT_T *intstat=powfs[ipowfs].intstat=mycalloc(1,INTSTAT_T);
	if(parms->load.i0){
	    warning("Loading i0, gx, gy\n");
	    intstat->i0=dcellread("%s/powfs%d_i0",parms->load.i0, ipowfs);
	    intstat->gx=dcellread("%s/powfs%d_gx",parms->load.i0, ipowfs);
	    intstat->gy=dcellread("%s/powfs%d_gy",parms->load.i0, ipowfs);
	}else{
	    if(parms->powfs[ipowfs].piinfile){
		/*load psf. 1 for each wavefront sensor. */
		info2("Using 1 sepsf for each wfs when loading sepsf\n");
		intstat->nsepsf=parms->powfs[ipowfs].nwfs;
		intstat->sepsf=dccellnew(parms->powfs[ipowfs].nwfs, 1);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    dcell *sepsf=intstat->sepsf->p[jwfs]=dcellread("%s_wfs%d",parms->powfs[ipowfs].piinfile,iwfs);
		    double pmax=dmax(sepsf->p[0]);
		    if(sepsf->p[0]->p[0]>pmax*0.5){
			error("wfs %d:  psf must have peak at center, not corner.\n", iwfs);
		    }
		    if(sepsf->nx!=nsa){
			error("piinfile doesn't match\n");
		    }
		}
	    }else if(parms->powfs[ipowfs].lo && parms->powfs[ipowfs].order<=2){
		error("Please specify piinfile for lo order phy wfs\n");
	    }else{
		/*Gen short exposure OTFs due to atmosphere*/
		genseotf(parms, powfs, ipowfs);
		if(parms->powfs[ipowfs].llt){
		    genselotf(parms, powfs, ipowfs);
		}
		/*Generating short exposure psfs for both uplink and downlink turbulence effect. */
		gensepsf(parms,powfs,ipowfs);
		if(parms->save.setup>1 && intstat){
		    writebin(intstat->sepsf->p[0],"powfs%d_sepsf",ipowfs);
		}
		/*Free short exposure otf. */
		ccellfree(intstat->lotf);
		cellfree(intstat->otf);
	    }
	    /*generate short exposure i0,gx,gy from psf. */
	    gensei(parms,powfs,ipowfs);
	    if(parms->save.setup){
		writebin(intstat->i0,"powfs%d_i0",ipowfs);
		writebin(intstat->gx,"powfs%d_gx",ipowfs);
		writebin(intstat->gy,"powfs%d_gy",ipowfs);
	    }
   	}
    }
    /*Generating Matched filter */
    if(parms->powfs[ipowfs].phytype==1 || parms->powfs[ipowfs].phytypesim==1){
	genmtch(parms,powfs,ipowfs);
	if(parms->save.setup){
	    writebin(powfs[ipowfs].intstat->mtche,"powfs%d_mtche",ipowfs);
	}
    }
    if(parms->powfs[ipowfs].phytype==2 || parms->powfs[ipowfs].phytypesim==2 || parms->powfs[ipowfs].dither){
	setup_powfs_cog(parms, powfs, ipowfs);
    }
    if(parms->save.setup){
	writebin(powfs[ipowfs].saneaxy,"powfs%d_saneaxy",ipowfs);
    }
    if(parms->powfs[ipowfs].neaphy){/*use physical optics nea for geom grad*/
	powfs[ipowfs].neasim=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    dmat **sanea=NULL;
	    if(powfs[ipowfs].saneaxy->ny==1){
		sanea=powfs[ipowfs].saneaxy->p;
	    }else{
		sanea=powfs[ipowfs].saneaxy->p+jwfs;
	    }
	    dmat *nea=dnew(nsa,3);
	    dmat*  pnea=nea/*PDMAT*/;
	    for(int isa=0; isa<nsa; isa++){
		dmat *neal=dchol(sanea[isa]);
		IND(pnea,isa,0)=neal->p[0];
		IND(pnea,isa,1)=neal->p[3];
		IND(pnea,isa,2)=neal->p[1];
		dfree(neal);
	    }
	    powfs[ipowfs].neasim->p[jwfs]=nea;
	}
	if(parms->save.setup){
	    writebin(powfs[ipowfs].neasim,"powfs%d_neasim", ipowfs);
	}
    }
  
    disable_save=disable_save_save;//put it back.
}
/*
  Setup gradient offset for calibration. opdadd is the wavefront aberration in
  WFS due to optics, without DM correction. opdbias is the wavefront aberration
  in WFS after DM system flat is applied.
*/
void setup_powfs_calib(const PARMS_T *parms, POWFS_T *powfs, loccell *aloc, dcell *dm_ncpa){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	dcellcp(&powfs[ipowfs].opdbias, powfs[ipowfs].opdadd);
	if(aloc && dm_ncpa){//apply DM figure to powfs surface
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		const double hs=parms->wfs[iwfs].hs;
		const double hc=parms->powfs[ipowfs].hc;
		double thetax=parms->wfs[iwfs].thetax;
		double thetay=parms->wfs[iwfs].thetay;
		if(!powfs[ipowfs].opdbias){
		    powfs[ipowfs].opdbias=dcellnew(parms->powfs[ipowfs].nwfs, 1);
		}
		if(!powfs[ipowfs].opdbias->p[jwfs]){
		    powfs[ipowfs].opdbias->p[jwfs]=dnew(powfs[ipowfs].loc->nloc, 1);
		}
		for(int idm=0; idm<parms->ndm; idm++){
		    if(!dm_ncpa->p[idm] || dm_ncpa->p[idm]->nx==0) continue;
		    double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg-hc;
		    double scale=1.-ht/hs;
		    double dispx=ht*thetax;
		    double dispy=ht*thetay;
		    prop_nongrid(aloc->p[idm], dm_ncpa->p[idm]->p, 
				 powfs[ipowfs].loc, powfs[ipowfs].opdbias->p[jwfs]->p, 
				 -1, dispx, dispy, scale, 0, 0);
		}
	    }
	}
	//if opdadd is null, but dm_ncpa is not, there will still be opdbias.
	if(powfs[ipowfs].opdbias){
	    if(parms->sim.ncpa_ttr){
		/*remove average tilt from opdbias and same amount from
		  opdadd. Does not need to be very accurate.*/
		dmat *mcc=loc_mcc_ptt(powfs[ipowfs].loc, powfs[ipowfs].amp->p);
		dinvspd_inplace(mcc);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    double ptt[3]={0,0,0};
		    loc_calc_ptt(NULL, ptt, powfs[ipowfs].loc, 1./mcc->p[0], mcc, 
				 powfs[ipowfs].amp->p, powfs[ipowfs].opdbias->p[jwfs]->p);
		    loc_remove_ptt(powfs[ipowfs].opdbias->p[jwfs]->p, ptt, powfs[ipowfs].loc);
		    loc_remove_ptt(powfs[ipowfs].opdadd->p[jwfs]->p, ptt, powfs[ipowfs].loc);
		}
		dfree(mcc);
	    }
	    if(parms->save.setup){
		writebin(powfs[ipowfs].opdbias, "powfs%d_opdbias", ipowfs);
	    }
	    if(parms->powfs[ipowfs].ncpa_method==1){//gradient offset
		if(!powfs[ipowfs].gradncpa){
		    powfs[ipowfs].gradncpa=dcellnew(parms->powfs[ipowfs].nwfs,1);
		}
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    if(powfs[ipowfs].opdbias->p[jwfs]){
			double *realamp=powfs[ipowfs].realamp->p[jwfs]->p;
			if(parms->powfs[ipowfs].type==1){//pywfs
			    dmat *ints=0;
			    pywfs_fft(&ints, powfs[ipowfs].pywfs, powfs[ipowfs].opdbias->p[jwfs]);
			    //writebin(powfs[ipowfs].opdbias->p[jwfs], "opdbias\n");exit(0);
			    pywfs_grad(&powfs[ipowfs].gradncpa->p[jwfs], powfs[ipowfs].pywfs, ints);
			    dfree(ints);
			}else if(parms->powfs[ipowfs].gtype_sim==1){
			    pts_ztilt(&powfs[ipowfs].gradncpa->p[jwfs], powfs[ipowfs].pts,
				      powfs[ipowfs].saimcc->p[powfs[ipowfs].nsaimcc>1?jwfs:0], 
				      realamp, powfs[ipowfs].opdbias->p[jwfs]->p);
			}else{
			    dspmm(&powfs[ipowfs].gradncpa->p[jwfs],adpind(powfs[ipowfs].GS0, jwfs),
				  powfs[ipowfs].opdbias->p[jwfs],"nn",1);
			}
		    }
		}
		if(powfs[ipowfs].gradncpa && parms->save.setup){
		    writebin(powfs[ipowfs].gradncpa, "powfs%d_gradncpa", ipowfs);
		}
	    }
	}
    }
}

/**
   Setup the powfs struct based on parms and aper. Everything about wfs are
   setup here.  \callgraph */
POWFS_T * setup_powfs_init(const PARMS_T *parms, APER_T *aper){
    TIC;tic;
    POWFS_T *powfs=mycalloc(parms->npowfs,POWFS_T);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].type==0){
	    info2("\n\033[0;32mSetting up powfs %d geom\033[0;0m\n\n", ipowfs);
	    setup_powfs_geom(powfs,parms,aper,ipowfs);
	    setup_powfs_grad(powfs,parms,ipowfs);
	}else if(parms->powfs[ipowfs].type==1){
	    info2("\n\033[0;32mSetting up powfs %d in Pyramid mode\033[0;0m\n\n", ipowfs);
	    pywfs_setup(powfs, parms, aper, ipowfs);
	}else{
	    error("powfs %d: invalid wfstype=%d\n", ipowfs, parms->powfs[ipowfs].type);
	}
    }
    toc2("setup_powfs_init");
    return powfs;
}
//this is only for shwfs.
void setup_powfs_phy(const PARMS_T *parms, POWFS_T *powfs){
    TIC;tic;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].type!=0) continue;
	info2("\n\033[0;32mSetting up powfs %d PO WFS\033[0;0m\n\n", ipowfs);
	if(parms->powfs[ipowfs].usephy
	   ||parms->powfs[ipowfs].psfout
	   ||parms->powfs[ipowfs].pistatout
	   ||parms->powfs[ipowfs].neaphy){
	    /*We have physical optics. setup necessary struct */
	    setup_powfs_prep_phy(powfs,parms,ipowfs);
	    setup_powfs_dtf(powfs,parms,ipowfs);
	    if(parms->powfs[ipowfs].llt){
		/*prepare Laser launch telescope. */
		setup_powfs_sodium(powfs,parms,ipowfs);/*read sodium profile and smooth it */
		setup_powfs_etf(powfs,parms,ipowfs,0,parms->powfs[ipowfs].llt->colprep);/*etf for prep */
		if(!parms->powfs[ipowfs].llt->colsimdtrat){/*const etf for sim */
		    if(parms->powfs[ipowfs].llt->colprep==parms->powfs[ipowfs].llt->colsim){
			powfs[ipowfs].etfsim=powfs[ipowfs].etfprep;
		    }else{
			setup_powfs_etf(powfs,parms,ipowfs,1,parms->powfs[ipowfs].llt->colsim);
		    }
		}
		setup_powfs_llt(powfs,parms,ipowfs);
	    }
	}
	if(parms->powfs[ipowfs].llt){
	    /*If there is LLT, setup the extra focus term if needed. */
	    setup_powfs_focus(powfs,parms,ipowfs);
	}
	if(parms->powfs[ipowfs].usephy || parms->powfs[ipowfs].neaphy){
	    setup_powfs_phygrad(powfs,parms,ipowfs);
	}
    }/*ipowfs */
    toc2("setup_powfs_phy");
}
void free_powfs_unused(const PARMS_T *parms, POWFS_T *powfs){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(powfs[ipowfs].intstat){
	    cellfree(powfs[ipowfs].intstat->sepsf);
	    cellfree(powfs[ipowfs].intstat->i0);
	    cellfree(powfs[ipowfs].intstat->gx);
	    cellfree(powfs[ipowfs].intstat->gy);
	}
    }
}
void free_powfs_shwfs(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    free_powfs_geom(powfs,ipowfs);
    dtf_free(powfs[ipowfs].dtf, parms->powfs[ipowfs].nwvl);
    dspcellfree(powfs[ipowfs].GS0);
    dcellfree(powfs[ipowfs].neasim);
    dcellfree(powfs[ipowfs].saneaxy);
    cellfree(powfs[ipowfs].cogcoeff);
    if(powfs[ipowfs].intstat){
	INTSTAT_T *intstat=powfs[ipowfs].intstat;
	cellfree(intstat->fotf);
	cellfree(intstat->potf);
	dcellfree(intstat->mtche);

	dfree(intstat->i0sum);
	free(intstat);
	powfs[ipowfs].intstat=NULL;
    }
    if(parms->powfs[ipowfs].llt){
	dcellfree(powfs[ipowfs].srot);
	dcellfree(powfs[ipowfs].srsa);
	dfree(powfs[ipowfs].srsamax);
	dcellfree(powfs[ipowfs].sprint);
    }

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
    if(powfs[ipowfs].etfprep!=powfs[ipowfs].etfsim){
	etf_free(powfs[ipowfs].etfprep, parms->powfs[ipowfs].nwvl);
    }
    etf_free(powfs[ipowfs].etfsim, parms->powfs[ipowfs].nwvl);
    etf_free(powfs[ipowfs].etfsim2, parms->powfs[ipowfs].nwvl);
    dcellfree(powfs[ipowfs].opdadd);
    dcellfree(powfs[ipowfs].opdbias);
    dcellfree(powfs[ipowfs].gradncpa);
    dfree(powfs[ipowfs].dtheta);
}
/**
   Free all parameters of powfs at the end of simulation.
*/
void free_powfs(const PARMS_T *parms, POWFS_T *powfs){
    free_powfs_unused(parms, powfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].type==0){
	    free_powfs_shwfs(parms, powfs, ipowfs);
	}else if(parms->powfs[ipowfs].type==1){
	    pywfs_free(powfs[ipowfs].pywfs);
	}
    }
    free(powfs);
}
