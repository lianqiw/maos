/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include "maos.h"
#include "setup_powfs.h"
#include "mtch.h"
#include "genseotf.h"
#define TEST_POWFS 0
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
free_powfs_geom(POWFS_T *powfs,  const PARMS_T *parms, int ipowfs){
    if(!powfs[ipowfs].pts){
	return;
    }
    ptsfree(powfs[ipowfs].pts);
    dfree(powfs[ipowfs].saa);
    locfree(powfs[ipowfs].saloc);
    locfree(powfs[ipowfs].loc);
    dfree(powfs[ipowfs].amp);
    if(powfs[ipowfs].loc_tel){
	locarrfree(powfs[ipowfs].loc_tel, parms->powfs[ipowfs].nwfs);
	dcellfree(powfs[ipowfs].amp_tel);
	dcellfree(powfs[ipowfs].saa_tel);
    }
    if(powfs[ipowfs].loc_dm){
	locarrfree(powfs[ipowfs].loc_dm, parms->powfs[ipowfs].nwfs*parms->ndm);
    }
    dcellfree(powfs[ipowfs].realamp);
    dcellfree(powfs[ipowfs].realsaa);
    dfree(powfs[ipowfs].sumamp);
    dfree(powfs[ipowfs].sumamp2);
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
/**
   Create WFS amplitude map from coordinate, masked with annular defined by (D,Din). 
*/
static dmat *mkwfsamp(loc_t *loc, map_t *ampground, double misregx, double misregy, double D, double Din){
    dmat *amp=dnew(loc->nloc, 1);
    if(ampground){
	prop_grid(ampground, loc, NULL, amp->p, 1, misregx, misregy,1,0,0,0);
    }else{
	locannular(amp->p, loc, -misregx, -misregy, D*0.5, Din*0.5, 1);
    }
    return amp;
}

static void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, int iwfs){
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
    free_powfs_geom(powfs, parms, ipowfs);
    int nwfsp=parms->powfs[ipowfs].nwfs;
    /*order of the system. 60 for TMT */
    const int order  = parms->powfs[ipowfs].order;
    /*Subaperture lateral length. 0.5 for TMT. */
    const double dxsa = parms->aper.d/(double)order;
    /*The subaperture area that are considered as full. */
    double areafulli;
    if(order>2){
	areafulli=1.;
    }else{
	areafulli=4./M_PI;
    }

    /*Number of OPD pixels in 1 direction in each
      subaperture. Make it even to do fft.*/
    int nx = 2*(int)round(0.5*dxsa/parms->powfs[ipowfs].dx);
    const double dx=dxsa/nx;/*adjust dx. */
    const double dxoffset=dx*0.5;//Always keep points inside subaperture for simulation.
    if(fabs(parms->powfs[ipowfs].dx-dx)>EPS)
	warning("Adjusting dx from %g to %g\n", parms->powfs[ipowfs].dx,dx);
    if(fabs(dxsa - nx * dx)>EPS){
	warning("nx=%d,dsa=%f,dx=%f not agree\n", nx, dxsa, dx);
    }
    info2("There are %d points in each subaperture of %gm.\n", nx, dxsa);
    const int nxsa=nx*nx;/*Total Number of OPD points. */
    if(parms->powfs[ipowfs].saloc){
	powfs[ipowfs].saloc=locread("%s", parms->powfs[ipowfs].saloc);
	if(fabs(powfs[ipowfs].saloc->dx-dxsa)>0.1*(powfs[ipowfs].saloc->dx+dxsa)){
	    error("loaded saloc has dx=%g, while powfs.order implies %g\n",
		  powfs[ipowfs].saloc->dx, dxsa);
	}
    }else{
	/*The coordinate of the subaperture (lower left coordinate) */
	powfs[ipowfs].saloc=locnew(order*order, dxsa, dxsa);
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
	double r2min=dxsa<parms->aper.din?pow(parms->aper.din/dxsa/2,2):-1;
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
		    powfs[ipowfs].saloc->locx[count]=xc*dxsa;
		    powfs[ipowfs].saloc->locy[count]=yc*dxsa;
		    count++;
		}
	    }
	}
	powfs[ipowfs].saloc->nloc=count;
    }
    /*convert saloc to pts*/
    powfs[ipowfs].pts=ptsnew(powfs[ipowfs].saloc->nloc, dxsa, dxsa, nx, dx, dx);
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	powfs[ipowfs].pts->origx[isa]=powfs[ipowfs].saloc->locx[isa]+dxoffset;
	powfs[ipowfs].pts->origy[isa]=powfs[ipowfs].saloc->locy[isa]+dxoffset;
    }
    /*Calculate the amplitude for each subaperture OPD point by ray tracing from
      the pupil amplitude map.*/
    powfs[ipowfs].loc=pts2loc(powfs[ipowfs].pts);
    /*The assumed amp. */
    powfs[ipowfs].amp=mkwfsamp(powfs[ipowfs].loc, aper->ampground, 
			       parms->misreg.pupil[0],parms->misreg.pupil[1], 
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
	PDMAT(ampi, pampi);
	double alpha=(nx*nx*parms->powfs[ipowfs].safill2d-(nx-2*nedge)*(nx-2*nedge))
	    /((nx-2*(nedge-1))*(nx-2*(nedge-1))-(nx-2*nedge)*(nx-2*nedge));
	double tot=0;
	for(int iy=nedge-1; iy<nx-nedge+1; iy++){
	    for(int ix=nedge-1; ix<nx-nedge+1; ix++){
		pampi[iy][ix]=1;
		if(ix==nedge-1 || ix==nx-nedge || iy==nedge-1 || iy==nx-nedge){
		    pampi[iy][ix]=alpha;
		}
		tot+=pampi[iy][ix];
	    }
	}
	if(parms->save.setup){
	    dwrite(ampi, "%s/powfs%d_ampi", dirsetup, ipowfs);
	}
	for(int isa=0; isa<powfs[ipowfs].pts->nsa; isa++){
	    for(int i=0; i<nx*nx; i++){
		powfs[ipowfs].amp->p[nx*nx*isa+i]*=ampi->p[i];
	    }
	}
	/*do not multiply to siglev. Already handled automatically*/
	thresarea*=parms->powfs[ipowfs].safill2d;
    }
 
    powfs[ipowfs].saa=wfsamp2saa(powfs[ipowfs].amp, nxsa);
    
    //Create another set of loc/amp that can be used to build GP. It has points on edge of subapertures
    map_t *map=create_metapupil_wrap(parms, 0, dx, dx, 0, 0, 0, 0, 0, 0);
    powfs[ipowfs].gloc=map2loc(map); mapfree(map);
    //do not use misregistration since this is the model
    powfs[ipowfs].gamp=mkwfsamp(powfs[ipowfs].gloc, aper->ampground, 
				0,0, parms->aper.d, parms->aper.din);
    loc_reduce(powfs[ipowfs].gloc, powfs[ipowfs].gamp, 1, NULL);

    if(parms->dbg.pupmask && parms->powfs[ipowfs].lo){//for NGS WFS only.
	if(nwfsp>1){
	    error("dbg.pupmask=1, powfs can only have 1 wfs.\n");
	}
	int iwfs=parms->powfs[ipowfs].wfs[0];
	wfspupmask(parms, powfs[ipowfs].loc, powfs[ipowfs].amp, iwfs);
	wfspupmask(parms, powfs[ipowfs].gloc, powfs[ipowfs].gamp, iwfs);
    }

    if(parms->misreg.tel2wfs){
	TIC;tic;
	/*
	  Misregistration/distortion from Telescope pupil to WFS pupil. The
	  amplitude map after misregistration/distortion will be used for
	  wavefront sensing.
	  They are not used for wavefront reconstruction
	*/
	int isset=0;
	powfs[ipowfs].loc_tel=calloc(nwfsp, sizeof(loc_t*));
	powfs[ipowfs].saa_tel=dcellnew(nwfsp, 1);
	powfs[ipowfs].amp_tel=dcellnew(nwfsp, 1);
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs[jwfs];
	    if(parms->misreg.tel2wfs[iwfs])
#pragma omp task shared(isset)
	    {
		isset=1;
		powfs[ipowfs].loc_tel[jwfs]
		    =loctransform(powfs[ipowfs].loc, parms->misreg.tel2wfs[iwfs]);
		powfs[ipowfs].amp_tel->p[jwfs]
		    =mkwfsamp(powfs[ipowfs].loc_tel[jwfs], aper->ampground,
			      -parms->misreg.pupil[0], -parms->misreg.pupil[1], 
			      parms->aper.d, parms->aper.din);
		if(parms->dbg.pupmask && parms->powfs[ipowfs].lo){
		    if(nwfsp>1){
			error("dbg.pupmask=1, powfs can only have 1 wfs.\n");
		    }
		    wfspupmask(parms, powfs[ipowfs].loc_tel[jwfs],
			       powfs[ipowfs].amp_tel->p[jwfs], iwfs);	
		}
		powfs[ipowfs].saa_tel->p[jwfs]=wfsamp2saa(powfs[ipowfs].amp_tel->p[jwfs], nxsa);
	    }
	}
#pragma omp taskwait
	if(!isset){
	    free(powfs[ipowfs].loc_tel); powfs[ipowfs].loc_tel=0;
	    dcellfree(powfs[ipowfs].saa_tel); 
	    dcellfree(powfs[ipowfs].amp_tel);
	}else{
	    toc("misreg.tel2wfs");
	}
    }/*if misreg */
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
    dmat *saa=NULL;
    if(powfs[ipowfs].saa_tel){
	warning_once("Todo: Improve to allow different sa for same wfs type\n");
	const double scale=1./(double) nwfsp;
	for(int i=0; i<nwfsp; i++){
	    dadd(&saa, 1, powfs[ipowfs].saa_tel->p[0], scale);
	}
    }else{
	saa=dref(powfs[ipowfs].saa);
    }
    if(dmax(saa)>1.01){
	warning("The sa area maxes to %g, which should be leq 1 (misregistration can cause this).\n", 
	      dmax(saa));
    }
    int count=0;
    for(int isa=0; isa<powfs[ipowfs].pts->nsa; isa++){
	if(saa->p[isa]>thresarea){
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
		    for(int jwfs=0; jwfs<nwfsp; jwfs++){
			MOVES(powfs[ipowfs].saa_tel->p[jwfs]->p, count, isa);
			MOVEDLOC(powfs[ipowfs].loc_tel[jwfs], nxsa, count, isa);
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

    powfs[ipowfs].npts = count*nxsa;
    powfs[ipowfs].nthread=count<parms->sim.nthread?count:parms->sim.nthread;
    ptsresize(powfs[ipowfs].pts, count);
    locresize(powfs[ipowfs].saloc, count);
    locresize(powfs[ipowfs].loc, count*nxsa);
    dresize(powfs[ipowfs].saa, count, 1);
    dresize(powfs[ipowfs].amp, count*nxsa, 1);
    if(powfs[ipowfs].loc_tel){
	for(int jwfs=0; jwfs<nwfsp; jwfs++){
	    locresize(powfs[ipowfs].loc_tel[jwfs], count*nxsa);
	    dresize(powfs[ipowfs].amp_tel->p[jwfs], count*nxsa, 1);
	    dresize(powfs[ipowfs].saa_tel->p[jwfs], count, 1);
	}
    }	
    powfs[ipowfs].realamp=dcellnew(nwfsp, 1);
    powfs[ipowfs].realsaa=dcellnew(nwfsp, 1);
    powfs[ipowfs].sumamp=dnew(nwfsp, 1);
    powfs[ipowfs].sumamp2=dnew(nwfsp, 1);
    for(int jwfs=0; jwfs<nwfsp; jwfs++){
	dmat *realamp, *realsaa;
	if(powfs[ipowfs].loc_tel){
	    realamp=powfs[ipowfs].amp_tel->p[jwfs];
	    realsaa=powfs[ipowfs].saa_tel->p[jwfs];
	}else{
	    realamp=powfs[ipowfs].amp;
	    realsaa=powfs[ipowfs].saa;
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
	powfs[ipowfs].realsaa->p[jwfs]=dref(realsaa);
    }
    if(parms->powfs[ipowfs].fieldstop){
	warning("powfs%d: generating field stop \n", ipowfs);
	if(parms->powfs[ipowfs].nwvl>1){
	    error("Not implemented yet. need to do phase unwrap in wfsgrad.\n");
	}
	powfs[ipowfs].embed=loc_create_embed(&powfs[ipowfs].nembed, powfs[ipowfs].loc, 1);
	long nembed=powfs[ipowfs].nembed;
	powfs[ipowfs].fieldstop=dnew(nembed, nembed);
	double dtheta=parms->powfs[ipowfs].wvl[0]/(powfs[ipowfs].loc->dx*nembed); 
	double radius=parms->powfs[ipowfs].fieldstop/dtheta/2;
	dcircle(powfs[ipowfs].fieldstop, nembed/2+1, nembed/2+1, 1, 1, radius, 1);
	dfftshift(powfs[ipowfs].fieldstop);
    }
    dfree(ampi);
    if(parms->misreg.dm2wfs){
	TIC;tic;
	/*
	  Misregistration/distortion from DM to WFS pupil. Not about telescope
	  pupil. The distorted grid are used for ray tracing from DM to WFS.
	  They are not used for wavefront reconstruction
	*/
	powfs[ipowfs].loc_dm=calloc(nwfsp*parms->ndm, sizeof(loc_t*));
	int isset=0;
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int jwfs=0; jwfs<nwfsp; jwfs++)
#pragma omp task shared(isset)
	    {
		int iwfs=parms->powfs[ipowfs].wfs[jwfs];
		if(parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]){
		    powfs[ipowfs].loc_dm[jwfs+nwfsp*idm]
			=loctransform(powfs[ipowfs].loc, parms->misreg.dm2wfs[iwfs+idm*parms->nwfs]);
		    isset=1;
		}
	    }
	}
#pragma omp taskwait
	if(!isset){
	    free(powfs[ipowfs].loc_dm); powfs[ipowfs].loc_dm=0;
	}else{
	    toc("misreg.dm2wfs");
	}
    }
    if(parms->save.setup){
	locwrite((loc_t*)powfs[ipowfs].pts, "%s/powfs%d_pts",dirsetup,ipowfs);
	locwrite(powfs[ipowfs].saloc, "%s/powfs%d_saloc",dirsetup,ipowfs); 
	locwrite(powfs[ipowfs].gloc,"%s/powfs%d_gloc", dirsetup, ipowfs);
	dwrite(powfs[ipowfs].gamp,"%s/powfs%d_gamp", dirsetup, ipowfs);
	dwrite(powfs[ipowfs].saa,"%s/powfs%d_saa", dirsetup,ipowfs);
	locwrite(powfs[ipowfs].loc,"%s/powfs%d_loc",dirsetup,ipowfs);
	dwrite(powfs[ipowfs].amp, "%s/powfs%d_amp", dirsetup,ipowfs);
	if(powfs[ipowfs].loc_tel){
	    dcellwrite(powfs[ipowfs].saa_tel,"%s/powfs%d_saa_tel", dirsetup,ipowfs);
	    dcellwrite(powfs[ipowfs].amp_tel,"%s/powfs%d_amp_tel", dirsetup,ipowfs);
	    locarrwrite(powfs[ipowfs].loc_tel, nwfsp, "%s/powfs%d_loc_tel", dirsetup,ipowfs);
	}
	if(powfs[ipowfs].loc_dm){
	    for(int idm=0; idm<parms->ndm; idm++){
		for(int jwfs=0; jwfs<nwfsp; jwfs++){
		    locarrwrite(powfs[ipowfs].loc_dm, parms->ndm*nwfsp,
				"%s/powfs%d_loc_dm", dirsetup, ipowfs);
		}
	    }
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
	spcellfree(powfs[ipowfs].GS0);
	/*Setting up every gradient tilt (g-ztilt) */
	if(parms->load.GS0){
	    powfs[ipowfs].GS0=spcellread("powfs%d_GS0",ipowfs);
	    if(powfs[ipowfs].amp_tel){
		assert(powfs[ipowfs].GS0->nx==powfs[ipowfs].nwfs);
	    }else{
		assert(powfs[ipowfs].GS0->nx==1);
	    }
	}else{
	    double displace[2]={0,0};
	    /*This mkg takes about 5 seconds. */
	    if(powfs[ipowfs].amp_tel){
		powfs[ipowfs].GS0=spcellnew(powfs[ipowfs].nwfs, 1);
	    }else{
		powfs[ipowfs].GS0=spcellnew(1, 1);
	    }
	    for(int iwfs=0; iwfs<powfs[ipowfs].GS0->nx; iwfs++){
		powfs[ipowfs].GS0->p[iwfs]=mkg(powfs[ipowfs].loc, 
					       powfs[ipowfs].loc,
					       powfs[ipowfs].realamp->p[iwfs]->p,
					       powfs[ipowfs].saloc,
					       1, 1, displace, 1);
	    }
	}
    }
    if(parms->powfs[ipowfs].gtype_recon==1 ||parms->powfs[ipowfs].gtype_sim==1){
	/*setting up zernike best fit (ztilt) inv(M'*W*M). good for NGS. */
	if(parms->powfs[ipowfs].order>4) 
	    warning("Ztilt for high order wfs is not good");
	powfs[ipowfs].nsaimcc=MAX(1,(powfs[ipowfs].loc_tel?powfs[ipowfs].nwfs:1));
	int nsaimcc=powfs[ipowfs].nsaimcc;
	dcellfreearr(powfs[ipowfs].saimcc, nsaimcc);
	powfs[ipowfs].saimcc=calloc(nsaimcc, sizeof(dcell*));
	for(int imcc=0; imcc<nsaimcc; imcc++){
	    dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].realamp->p[imcc]->p);
	    powfs[ipowfs].saimcc[imcc]=dcellinvspd_each(mcc);
	    dcellfree(mcc);
	}
    }
 
    if(!parms->powfs[ipowfs].neaphy){
	powfs[ipowfs].neasim=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs[jwfs];
	    const long nsa=powfs[ipowfs].pts->nsa;
	    dmat *nea=NULL;
	    if(parms->powfs[ipowfs].neasimfile){
		if(parms->powfs[ipowfs].neasim!=-1){
		    error("neasimfile and neasim can not both be supplied\n");
		}
		nea=dread("%s_wfs%d",parms->powfs[ipowfs].neasimfile,iwfs);/*rad */
		if(nea->nx!=nsa || (nea->ny!=2 && nea->ny!=3)){
		    error("wfs %d: NEA read from %s_wfs%d has incorrect format. We want %ldx2(3) array.\n",
			  iwfs, parms->powfs[ipowfs].neasimfile, iwfs,powfs[ipowfs].pts->nsa);
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
		    warning("wfs %d: NEA=%g mas, too small. Unit error?\n",iwfs, neamax*206265000);
		}else if(neamax>parms->powfs[ipowfs].pixtheta){
		    warning("wfs %d: NEA=%g mas, too big. Unit error?\n",iwfs, neamax*206265000);
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
	    dcellwrite(powfs[ipowfs].neasim,"%s/powfs%d_neasim", dirsetup, ipowfs);
	}
    }
    if(parms->save.setup && powfs[ipowfs].GS0){
	spcellwrite(powfs[ipowfs].GS0,"%s/powfs%d_GS0",dirsetup,ipowfs);
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
    const double dxsa=powfs[ipowfs].pts->dsa;
    const int nsa=powfs[ipowfs].pts->nsa;
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
	    double ox2=parms->powfs[ipowfs].llt->ox[illt]-dxsa*0.5;
	    double oy2=parms->powfs[ipowfs].llt->oy[illt]-dxsa*0.5;
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
	    int pnsa=(int)ceil(sqrt(rsa2max)/dxsa)+1;
	    
	    double prot[pnsa];
	    powfs[ipowfs].sprint->p[illt]=dnew(pnsa,1);
	    double *pp=powfs[ipowfs].sprint->p[illt]->p;
	    for(int ind=0; ind<pnsa; ind++){
		prot[ind]=INFINITY;
		pp[ind]=-1;
	    }
	    double desrot=0;
	    if(fabs(parms->powfs[ipowfs].llt->ox[illt])<dxsa 
	       && fabs(parms->powfs[ipowfs].llt->oy[illt])<dxsa){
		desrot=0;
	    }else{
		double ddx=(0-parms->powfs[ipowfs].llt->ox[illt]);
		double ddy=(0-parms->powfs[ipowfs].llt->oy[illt]);
		desrot=atan2(ddy,ddx);
	    }
	    for(int isa=0; isa<nsa; isa++){
		int ind=(int)round(powfs[ipowfs].srsa->p[illt]->p[isa]/dxsa);
		double irot=fabs(powfs[ipowfs].srot->p[illt]->p[isa]-desrot);
		if(irot<prot[ind]){
		    prot[ind]=irot;
		    pp[ind]=(double)isa;
		}
	    }
	}
	if(parms->save.setup){
	    dcellwrite(powfs[ipowfs].sprint,"%s/powfs%d_sprint", dirsetup,ipowfs);
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
	double wvlmin=parms->powfs[ipowfs].wvl[0];
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(wvlmin>parms->powfs[ipowfs].wvl[iwvl])
		wvlmin=parms->powfs[ipowfs].wvl[iwvl];
	}
	double dtheta=wvlmin/(dxsa*embfac);/*Min PSF sampling. */
	ncompx=4*(int)round(0.25*pixpsax*pixthetax/dtheta);
	ncompy=4*(int)round(0.25*pixpsay*pixthetay/dtheta);
	
	if(!parms->powfs[ipowfs].radrot){
	    /*Follow laos method, need square  */
	    ncompx=ncompx>ncompy?ncompx:ncompy;
	    ncompy=ncompx;
	}else{
	    /*Found that: Must set ncompx==ncompy for
	      rotationg either psf or otf. reduce aliasing
	      and scattering of image intensities.
	    */
	    ncompx=ncompx>ncompy?ncompx:ncompy;
	    ncompy=ncompx;
	}
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
	   (powfs[ipowfs].bkgrnd->nx!=powfs[ipowfs].pts->nsa
	    ||(powfs[ipowfs].bkgrnd->ny!=1
	       && powfs[ipowfs].bkgrnd->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrnd is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrnd->nx, powfs[ipowfs].bkgrnd->ny,
		  powfs[ipowfs].pts->nsa, parms->powfs[ipowfs].nwfs);
	}
	if(parms->powfs[ipowfs].bkgrndfnc && 
	   (powfs[ipowfs].bkgrndc->nx!=powfs[ipowfs].pts->nsa
	    ||(powfs[ipowfs].bkgrndc->ny!=1
	       && powfs[ipowfs].bkgrndc->ny!=parms->powfs[ipowfs].nwfs))){
	    error("powfs%d: bkgrndc is of dimension %ld x %ld, "
		  "but should be %ld x 1 or %d\n",
		  ipowfs, powfs[ipowfs].bkgrndc->nx, powfs[ipowfs].bkgrndc->ny,
		  powfs[ipowfs].pts->nsa, parms->powfs[ipowfs].nwfs);
	}
    }
}
/**
   Free the detector transfer function.
*/
static void 
free_powfs_dtf(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(powfs[ipowfs].dtf){
	for(int iwvl=0;iwvl<parms->powfs[ipowfs].nwvl;iwvl++){
	    ccellfree(powfs[ipowfs].dtf[iwvl].nominal);
	    spcellfree(powfs[ipowfs].dtf[iwvl].si);
	    cfree(powfs[ipowfs].dtf[iwvl].Ux);
	    cfree(powfs[ipowfs].dtf[iwvl].Uy);
	}
	free(powfs[ipowfs].dtf);
	dfree(powfs[ipowfs].dtheta);
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
    /*wvl independent parameters */
    const double dxsa=powfs[ipowfs].pts->dsa;
    const int nsa=powfs[ipowfs].pts->nsa;
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const double blurx=parms->powfs[ipowfs].pixblur*pixthetax;
    const double blury=parms->powfs[ipowfs].pixblur*pixthetay;/*Is this right?*/
    const double e0x=-2*M_PI*M_PI*blurx*blurx;/*blurring factors */
    const double e0y=-2*M_PI*M_PI*blury*blury;
    const int do_blur=fabs(blurx)>EPS && fabs(blury)>EPS;
    const int ncompx=powfs[ipowfs].ncompx;
    const int ncompy=powfs[ipowfs].ncompy;
    const int ncompx2=ncompx>>1;
    const int ncompy2=ncompy>>1;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int pixpsax=powfs[ipowfs].pixpsax;
    const int pixpsay=powfs[ipowfs].pixpsay;
    const double embfac=parms->powfs[ipowfs].embfac;
    const double pxo=-(pixpsax*0.5-0.5+parms->powfs[ipowfs].pixoffx)*pixthetax;
    const double pyo=-(pixpsay*0.5-0.5+parms->powfs[ipowfs].pixoffy)*pixthetay;
    int ndtf;
    int nllt;
    int multi_dtf=0;
    if(parms->powfs[ipowfs].llt&&!parms->powfs[ipowfs].radrot
       &&parms->powfs[ipowfs].radpix){
	/*When we have llt, there is elongation, radial pix but
	  not use rotating psf/otf Need to create nominal/si for
	  each subaperture. DTF is on x/y coordinate, so pixels
	  and pixel coordinates need to rotated to r/a
	  direction. */
	ndtf=nsa;
	nllt=parms->powfs[ipowfs].llt->n;
	multi_dtf=1;
    }else{
	/*We only need a single DTF. */
	ndtf=1;
	nllt=1;
	multi_dtf=0;
    }
    powfs[ipowfs].dtf=calloc(nwvl,sizeof(DTF_T));
    powfs[ipowfs].dtheta=dnew(nwvl,1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double wvl=parms->powfs[ipowfs].wvl[iwvl];
	const double dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
	const double dux=1./(dtheta*ncompx);
	const double duy=1./(dtheta*ncompy);
	const double dux2=dux*dux;
	const double duy2=duy*duy;
	const double pdtheta=pixthetax*pixthetay/(dtheta*dtheta);
	const double duxp=dux*pixthetax;
	const double duyp=duy*pixthetay;
	powfs[ipowfs].dtheta->p[iwvl]=dtheta;
	powfs[ipowfs].dtf[iwvl].nominal=ccellnew(ndtf,nllt);
	powfs[ipowfs].dtf[iwvl].si=spcellnew(ndtf,nllt);
	cmat*(*nominals)[ndtf]=
	    (void*)powfs[ipowfs].dtf[iwvl].nominal->p;
	/*both nominal and si depends on wavelength.*/
	dsp*(*sis)[ndtf]=(void*)powfs[ipowfs].dtf[iwvl].si->p;
	cmat *nominal=cnew(ncompx,ncompy);
	cfft2plan(nominal,-1);
	cfft2plan(nominal,1);
	PCMAT(nominal,pn);
	loc_t *loc_psf=mksqloc(ncompx,ncompy,dtheta,dtheta,-ncompx2*dtheta, -ncompy2*dtheta);
	double theta=0;
	double ct=cos(theta);
	double st=sin(theta);
	for(int illt=0; illt<nllt; illt++){
	    for(int isa=0; isa<ndtf; isa++){
		if(multi_dtf){
		    theta=powfs[ipowfs].srot->p[illt]->p[isa];
		    ct=cos(theta);
		    st=sin(theta);
		}
		
		for(int iy=0; iy<ncompy; iy++){
		    int jy=iy-ncompy2;
		    for(int ix=0; ix<ncompx; ix++){
			int jx=ix-ncompx2;
			double ir=ct*jx+st*jy;
			double ia=-st*jx+ct*jy;
			pn[iy][ix]=sinc(ir*duyp)*sinc(ia*duxp)*pdtheta;
			if(do_blur){
			    pn[iy][ix]*=exp(e0x*(ir*ir*duy2)+e0y*(ia*ia*dux2));
			}
		    }
		}
		/*put peak in corner. pretreat nominal so
		  that we avoid fftshift when forming i0.*/
		
		cfftshift(nominal);
		cfft2(nominal,-1);
		cfftshift(nominal);
		cfft2(nominal,1);
		cscale(nominal,1./(double)(nominal->nx*nominal->ny));
		ccp(&nominals[illt][isa], nominal);
		loc_t *loc_ccd=mksqlocrot(pixpsax,pixpsay, pixthetax,pixthetay,pxo,pyo,theta);
		sis[illt][isa]=mkh(loc_psf,loc_ccd,NULL,0,0,1,0,0);
		locfree(loc_ccd);
	    }/*isa */
	}/*illt */
	cfree(nominal);
	locfree(loc_psf);
	if(parms->save.setup>1){
	    ccellwrite(powfs[ipowfs].dtf[iwvl].nominal,
		       "%s/powfs%d_dtf%d_nominal",dirsetup,ipowfs,iwvl);
	    spcellwrite(powfs[ipowfs].dtf[iwvl].si,
			"%s/powfs%d_dtf%d_si",dirsetup,ipowfs,iwvl);
	}
	/*Create an excessive high frequency in nominal so that
	  we don't have to do fftshift later.*/
	powfs[ipowfs].dtf[iwvl].Ux=cnew(ncompx,1);
	powfs[ipowfs].dtf[iwvl].Uy=cnew(ncompy,1);
	dcomplex *Ux=powfs[ipowfs].dtf[iwvl].Ux->p;
	dcomplex *Uy=powfs[ipowfs].dtf[iwvl].Uy->p;

	/*The following is used in genseotf to compute shifted
	  i0.  not used currently*/
	for(int ix=0; ix<ncompx; ix++){
	    int jx=ix<ncompx2?ix:(ix-ncompx);
	    Ux[ix]=-2.*I*M_PI*jx*dux;
	}
	for(int iy=0; iy<ncompy; iy++){
	    int jy=iy<ncompy2?iy:(iy-ncompy);
	    Uy[iy]=-2.*I*M_PI*jy*duy;
	}

    }/*iwvl */
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
    dwrite(powfs[ipowfs].focus, "powfs%d_focus", ipowfs);
}

/**
   Load and smooth out sodium profile. We preserve the sum of the sodium profile,
   which represents a scaling factor of the signal level.  */
static void setup_powfs_sodium(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    char *fnprof=parms->powfs[ipowfs].llt->fnprof;
    dmat *Nain=dread("%s",fnprof);
    if(Nain->ny<2){
	error("The sodium profile input %s is in wrong fromat\n", fnprof);
    }

    if(parms->dbg.na_smooth){/*resampling the sodium profile by binning. */
	/*Make new sampling: */
	double rsamax=dmax(powfs[ipowfs].srsamax);
	double dthetamin=dmin(powfs[ipowfs].dtheta);
	const double x0in=Nain->p[0];
	const long nxin=Nain->nx;
	double dxin=(Nain->p[nxin-1]-x0in)/(nxin-1);
	/*minimum sampling required. */
	double dxnew=pow(parms->powfs[ipowfs].hs,2)/rsamax*dthetamin;
	if(dxnew > dxin * 2){
	    info2("Smoothing sodium profile\n");
	    const long nxnew=ceil((Nain->p[nxin-1]-x0in)/dxnew);
	    loc_t *loc_in=mk1dloc_vec(Nain->p, nxin);
	    loc_t *loc_out=mk1dloc(x0in, dxnew, nxnew);
	    dsp *ht=mkhb(loc_out, loc_in, NULL, 0, 0, 1,0,0);
	    powfs[ipowfs].sodium=dnew(nxnew, Nain->ny);
	    memcpy(powfs[ipowfs].sodium->p, loc_out->locx, sizeof(double)*nxnew);
			     
	    for(long icol=1; icol<Nain->ny; icol++){
		/*input profile */
		double *pin=Nain->p + nxin * icol;
		/*output profile */
		double *pout=powfs[ipowfs].sodium->p + nxnew*icol;
		/*preserve sum of input profile */
		double Nasum=dblsum(pin, nxin);
		spmulvec(pout, ht, pin, 1);
		normalize_sum(pout, nxnew, Nasum);
	    }
	    spfree(ht);
	    locfree(loc_in);
	    locfree(loc_out);
	    dfree(Nain);
	}else{
	    info2("Not smoothing sodium profile\n");
	    powfs[ipowfs].sodium=Nain;
	}
    }else{
	info2("Not smoothing sodium profile\n");
	powfs[ipowfs].sodium=Nain;
    }
    if(parms->save.setup){
	dwrite(powfs[ipowfs].sodium, "%s/powfs%d_sodium",dirsetup,ipowfs);
    }
}
/**
   Compute Elongation Transfer function.
   - mode=0: for preparation.
   - mode=1: for simulation.
*/
void setup_powfs_etf(POWFS_T *powfs, const PARMS_T *parms, int ipowfs, int mode, int istep){
    if(!parms->powfs[ipowfs].llt) return;
    const double dxsa=powfs[ipowfs].pts->dsa;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int nllt=parms->powfs[ipowfs].llt->n;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int ncompx=powfs[ipowfs].ncompx;
    const int ncompy=powfs[ipowfs].ncompy;
    const int ncompx2=ncompx>>1;
    const int ncompy2=ncompy>>1;
    ETF_T *powfsetf=NULL;
    int colstart, colskip;
    if(mode==0){/*preparation. */
	powfs[ipowfs].etfprep=calloc(nwvl, sizeof(ETF_T));
	powfsetf=powfs[ipowfs].etfprep;
	colskip=parms->powfs[ipowfs].llt->colprep;
    }else{
	if(parms->powfs[ipowfs].llt->colprep==parms->powfs[ipowfs].llt->colsim
	   && parms->powfs[ipowfs].llt->colsimdtrat==0){
	    if(powfs[ipowfs].etfsim){
		error("Etfsim should be empty\n");
	    }
	    if(!powfs[ipowfs].etfprep){
		error("Please call setup_powfs_etf with mode = 0 first\n");
	    }
	    powfs[ipowfs].etfsim=powfs[ipowfs].etfprep;
	    info2("Simulation and reconstruction using same etf\n");
	    return;
	}
	if(powfs[ipowfs].etfsim){
	    info2("Free previous etfsim\n");
	    for(int iwvl=0;iwvl<parms->powfs[ipowfs].nwvl;iwvl++){
		ccellfree(powfs[ipowfs].etfsim[iwvl].p1);
		ccellfree(powfs[ipowfs].etfsim[iwvl].p2);
	    }
	    free(powfs[ipowfs].etfsim);
	}
	powfs[ipowfs].etfsim=calloc(nwvl, sizeof(ETF_T));
	powfsetf=powfs[ipowfs].etfsim;
	colskip=parms->powfs[ipowfs].llt->colsim;
    }

    if(!powfs[ipowfs].sodium){
	error("Sodium profile is NULL\n");
    }
    dmat *sodium=powfs[ipowfs].sodium;
    const int nhp=sodium->nx; 
    const double hs=parms->powfs[ipowfs].hs;
    //adjusting for the zenith angle;
    double hpmin=0, dhp1=0;
    if(parms->dbg.na_interp){
	hpmin=sodium->p[0]/cos(parms->sim.za);
	dhp1=cos(parms->sim.za)/(sodium->p[1]-sodium->p[0]);
	/*assume linear spacing. check the assumption valid */
	if(fabs(sodium->p[nhp-1]-sodium->p[0]-(nhp-1)*(sodium->p[1]-sodium->p[0]))>1.e-7){
	    error("llt profile is not evenly spaced\n");
	}
    }

    if(sodium->ny-1<colskip){
	error("Invalid configuration. colprep or colsim is too big\n");
    }
    colstart=colskip+istep%(sodium->ny-1-colskip);
    info2("Na using column %d in %s\n",colstart, mode==0?"preparation":"simulation");
    /*points to the effective sodium profile. */
    double* pp=sodium->p+nhp*(1+colstart);
    const double* px=sodium->p;
    /*the sum of pp determines the scaling of the pixel intensity. */
    const double i0scale=dblsum(pp, nhp);
    if(fabs(i0scale-1)>0.01){
	warning("Siglev is scaled by %g by sodium profile\n", i0scale);
    }
    assert(nwvl==1);/*sanity check. need to double check the code is nwvl!=1. */
    const double embfac=parms->powfs[ipowfs].embfac;
    /*setup elongation along radial direction. don't care azimuthal. */
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double wvl=parms->powfs[ipowfs].wvl[iwvl];
	const double dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
	const double dux=1./(dtheta*ncompx);
	const double duy=1./(dtheta*ncompy);
	/**
	   fixme: should do integrate instead of just interpolating closest point. 
	*/
	cmat *(*petf)[nsa]=NULL;
	const int ndtf=powfs[ipowfs].dtf[iwvl].nominal->nx;
	cmat *(*pnominal)[ndtf]=(void*)powfs[ipowfs].dtf[iwvl].nominal->p;
	int mnominal;/*multiply with isa to get index into pnominal. */
	if(ndtf==1)
	    mnominal=0;
	else if(ndtf==nsa)
	    mnominal=1;
	else{
	    mnominal=-1;
	    error("Invalid nominal.\n");
	}
	int use1d;

	if(parms->powfs[ipowfs].radrot){
	    if(!parms->powfs[ipowfs].radpix){
		error("radrot can only be used with radpix\n");
	    }
	    /*
	      Use rotating psf/otf method to do radial
	      pixel. ETF is 1D only. Good for off-axis
	      launch, otherwise, ETF and DTF takes a lot of
	      space for 6 LGS in NFIRAOS setup.
	    */
	    warning("Rotate PSF to do radial format detector (preferred)\n");
	    powfsetf[iwvl].p1=ccellnew(nsa,nllt);
	    petf=(void*)powfsetf[iwvl].p1->p;
	    use1d=1;
	}else{
	    /*
	      Applied to 2 cases: 
	      1) radial ccd, not rotating PSF or OTF
	      2) Non-radial ccd
	      2010-01-04: Fuse dtf nominal into etf for this case.
	    */
	    if(parms->powfs[ipowfs].radpix){
		info2("2D ETF for Radial CCD\n");
	    }else{
		info2("Non-Radial CCD\n");
	    }
	    powfsetf[iwvl].p2=ccellnew(nsa,nllt);
	    petf=(void*)powfsetf[iwvl].p2->p;
	    use1d=0;
	}

	const int fuse_etf=(use1d==0);
	TIC;tic;
	if(!parms->dbg.na_interp){
	    /*
	      The ETF is computed as DFT
	      ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*{\theta}_j)P({\theta}_j)
	      Where \theta_j is the radial coord of pixel j, and
	      P({\theta}_j}=sodium(h_j) with h_j=rsa*hs/(rsa+hs*{\theta}_j)
	      where hs is GS distance and rsa is subaperture to LLT distance
	      We replace the above interpolation (and FFT) by doing summation directly
	      ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*(rsa/h_j-rsa/hs))P(h_j)
	    */
	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    /*1d ETF along radius. */
		    double rsa=powfs[ipowfs].srsa->p[illt]->p[isa];
		    double rsa_za=rsa*cos(parms->sim.za);
		    /* No interpolation, no fft */
		    if(use1d){
			petf[illt][isa]=cnew(ncompx,1);
			dcomplex *etf1d=petf[illt][isa]->p;
			for(int icompx=0; icompx<ncompx; icompx++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
			{
			    const double kr=dux*(icompx>=ncompx2?(icompx-ncompx):icompx);
			    for(int ih=0; ih<nhp; ih++){
				const double tmp=(-2*M_PI*(kr*(rsa_za/sodium->p[ih]-rsa/hs)));
				etf1d[icompx]+=pp[ih]*cos(tmp)+pp[ih]*sin(tmp)*I;
			    }
			}
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
		    }else{
			const double theta=powfs[ipowfs].srot->p[illt]->p[isa];
			const double ct=cos(theta);
			const double st=sin(theta);
			petf[illt][isa]=cnew(ncompx,ncompy);
			dcomplex (*etf2d)[ncompx]=(void*)petf[illt][isa]->p; 
			for(int icompy=0; icompy<ncompy; icompy++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
			{
			    const double ky=duy*(icompy>=ncompy2?(icompy-ncompy):icompy);
			    for(int icompx=0; icompx<ncompx; icompx++){
				const double kx=dux*(icompx>=ncompx2?(icompx-ncompx):icompx);
				double kr=(ct*kx+st*ky);/*along radial*/
				for(int ih=0; ih<nhp; ih++){
				    const double tmp=(-2*M_PI*(kr*(rsa_za/px[ih]-rsa/hs)));
				    etf2d[icompy][icompx]+=pp[ih]*cos(tmp)+pp[ih]*sin(tmp)*I;
				}
			    }
			}
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
		    }
		}//isa
	    }//illt
	}else{
	    const int npad=2;/*zero padding to reduce aliasing? */
	    const int nover=2;/*enough size for rotation*/
	    const int netf=ncompx*nover*npad;
	    const double dtetf=dtheta/nover;
	    const double dusc=(netf*dtetf)/(dtheta*ncompx);
	    cmat *etf=cnew(netf,1);
	    double *thetas=calloc(netf, sizeof(double));
	    const int netf2=netf>>1;
	    /*Only interpolating the center part. the rest is padding. */
	    const int etf0=netf2-(int)round(ncompx2*(dtheta/dtetf));
	    const int etf1=etf0+(int)round(ncompx*(dtheta/dtetf));
	    if(etf0<0) error("Invalid configuration\n");
	    for(int it=etf0; it<etf1; it++){
		thetas[it]=(it-netf2)*dtetf;
	    }
	    cfft2plan(etf, -1);

	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    /*1d ETF along radius. */
		    double rsa=powfs[ipowfs].srsa->p[illt]->p[isa];
		    double etf2sum=0;
		    czero(etf);
		    for(int icomp=etf0; icomp<etf1; icomp++){
			/*peak in center */
			const double itheta=thetas[icomp];
			/*non linear mapping. changed from - to + on 2014-05-07.*/
			const double ih=hs*rsa/(rsa+itheta*hs);
			/*interpolating to get Na profile strenght. */
			/*this is bilinear interpolation. not good. need to do averaging */
			const double iih=(ih-hpmin)*dhp1;
			const int iihf=ifloor(iih);
			const double iihw=iih-iihf;
			if(iihf<0 || iihf>nhp-2){
			    etf->p[icomp]=0.;
			}else{
			    double tmp=pp[iihf]*(1.-iihw)+pp[iihf+1]*iihw;
			    /*neglected rsa1 due to renormalization. */
			    etf->p[icomp]=tmp;
			    etf2sum+=tmp;
			}
		    }
		    if(fabs(etf2sum)>1.e-20){
			/*2010-11-09:

			  We used to normalize the etf before fft so that after fft
			  it max to 1. The strength of original profile doesn't
			  matter.
		    
			  Changed: We no longer normalize the etf, so we can model
			  the variation of the intensity and meteor trails.
		      
			*/
			cscale(etf,i0scale/etf2sum);
			cfftshift(etf);/*put peak in corner; */
			cfft2(etf, -1);
			if(use1d){
			    if(npad==1 && nover==1){
				ccp(&petf[illt][isa],etf);
			    }else{
				cfftshift(etf);
				petf[illt][isa]=cnew(ncompx,1);
				dcomplex *etf1d=petf[illt][isa]->p;
				for(int icompx=0; icompx<ncompx; icompx++){
				    double ir=dusc*(icompx-ncompx2)+netf2;
				    int iir=ifloor(ir);
				    ir=ir-iir;
				    if(iir>=0 && iir<netf-1){
					etf1d[icompx]=etf->p[iir]*(1.-ir)
					    +etf->p[iir+1]*ir;
				    }/*else{etf1d[icompx]=0;}*/
				}
				cfftshift(petf[illt][isa]);
			    }
			}else{
			    /*Rotate the ETF. */
			    /*need to put peak in center. */
			    cfftshift(etf);
			    double theta=powfs[ipowfs].srot->p[illt]->p[isa];
			    double ct=cos(theta);
			    double st=sin(theta);
			    petf[illt][isa]=cnew(ncompx,ncompy);
			    dcomplex (*etf2d)[ncompx]=(void*)petf[illt][isa]->p;
			    for(int icompy=0; icompy<ncompy; icompy++){
				double iy=(icompy-ncompy2);
				for(int icompx=0; icompx<ncompx; icompx++){
				    double ix=(icompx-ncompx2);
				    double ir=(dusc*(ct*ix+st*iy))+netf2;/*index in etf */
				    int iir=ifloor(ir);
				    ir=ir-iir;
				    if(iir>=0 && iir<netf-1){
					/*bilinear interpolation. */
					etf2d[icompy][icompx]=etf->p[iir]*(1.-ir)
					    +etf->p[iir+1]*ir;
				    }/*else{etf2d[icompy][icompx]=0;}*/
				}
			    }
			    cfftshift(petf[illt][isa]);/*peak in corner; */
			}
		    }else{
			warning_once("Wrong focus!\n");
			if(use1d){
			    petf[illt][isa]=cnew(ncompx,1);
			}else{
			    petf[illt][isa]=cnew(ncompx,ncompy);
			}
			cset(petf[illt][isa],1);
		    }
		}//for isa
	    }//for illt.
	    cfree(etf);
	    free(thetas);
	}//if na_interp
	toc2("ETF");
	if(fuse_etf){
	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    ccwm(petf[illt][isa], pnominal[illt][isa*mnominal]);
		}
	    }
	    powfs[ipowfs].dtf[iwvl].fused=1;
	    if(parms->powfs[ipowfs].llt->colprep==parms->powfs[ipowfs].llt->colsim
	       && parms->powfs[ipowfs].llt->colsimdtrat==0){
		ccellfree(powfs[ipowfs].dtf[iwvl].nominal);
		info2("DTF nominal is fused to ETF and freed\n");
	    }else{
		info2("DTF nominal is fused to ETF but kept\n");
	    }
	}	    
    }//for iwvl
}

/**
   setting up uplink pts/amp lotf
*/
static void 
setup_powfs_llt(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    double wvl0=parms->powfs[ipowfs].wvl[0];
    LLT_T *llt=powfs[ipowfs].llt=calloc(1, sizeof(LLT_T));
    LLT_CFG_T *lltcfg=parms->powfs[ipowfs].llt;
    pts_t *lpts=llt->pts=calloc(1, sizeof(pts_t));
    lpts->nsa=1;
    double lltd=lltcfg->d;
    lpts->dsay=lpts->dsa=MAX(lltd, powfs[ipowfs].pts->dsa);
    int notf=MAX(powfs[ipowfs].ncompx, powfs[ipowfs].ncompy);
    /*The otf would be dx/lambda. Make it equal to pts->dsa/lambda/notf)*/
    const double dx=lpts->dx=parms->powfs[ipowfs].embfac*powfs[ipowfs].pts->dsa/notf;
    lpts->dy=lpts->dx;
    const int nx=lpts->nx=round(lpts->dsa/lpts->dx);
    lpts->dsay=lpts->dsa=lpts->dx*lpts->nx;
    
    lpts->origx=calloc(1, sizeof(double));
    lpts->origy=calloc(1, sizeof(double));

    double oy=lpts->origx[0]=(dx-lpts->dsa)*0.5;
    double ox=lpts->origy[0]=(dx-lpts->dsa)*0.5;
    double sumamp2=0;
    llt->amp=dnew(nx,nx);
    if(lltcfg->fnamp){
	map_t *lltamp=mapread("%s", lltcfg->fnamp);
	prop_grid_pts(lltamp, llt->pts, NULL, llt->amp->p, 1, 0, 0, 1, 1, 0, 0);
	sumamp2=dinn(llt->amp, llt->amp);
	mapfree(lltamp);
    }else{
	PDMAT(llt->amp, amps);
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
		    amps[iy][ix]=exp(-r2/r2waist);
		    sumamp2+=pow(amps[iy][ix],2);
		}
	    }
	}
    }
    /*normalized so that max(otf)=1; */
    sumamp2=1./(sqrt(sumamp2));
    dscale(llt->amp, sumamp2);
    if(parms->powfs[ipowfs].llt->fnsurf){
	int nlotf;
	map_t **ncpa=maparrread(&nlotf, "%s",parms->powfs[ipowfs].llt->fnsurf);
	assert(nlotf==1 || nlotf==parms->powfs[ipowfs].nwfs);
	llt->ncpa=dcellnew(nlotf, 1);
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    llt->ncpa->p[ilotf]=dnew(nx,nx);
	    prop_grid_pts(ncpa[ilotf], llt->pts, NULL, llt->ncpa->p[ilotf]->p, 1, 0, 0, 1, 0, 0, 0);
	}
	maparrfree(ncpa, nlotf);
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
    if(parms->save.setup){
	locwrite(llt->loc, "%s/powfs%d_llt_loc",dirsetup,ipowfs);
	dwrite(llt->amp, "%s/powfs%d_llt_amp", dirsetup,ipowfs);
	dcellwrite(llt->imcc, "%s/powfs%d_llt_imcc",dirsetup,ipowfs);
	if(llt->ncpa){
	    dcellwrite(llt->ncpa, "%s/powfs%d_llt_ncpa", dirsetup, ipowfs);
	}
	dcellwrite(powfs[ipowfs].srot, "%s/powfs%d_srot",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].srsa, "%s/powfs%d_srsa",dirsetup,ipowfs);
    }
    if(parms->save.setup>1){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(powfs[ipowfs].etfprep[iwvl].p1){
		ccellwrite(powfs[ipowfs].etfprep[iwvl].p1, 
			   "%s/powfs%d_etfprep%d_1d",dirsetup,ipowfs,iwvl);
	    }
	    if(powfs[ipowfs].etfprep[iwvl].p2){
		ccellwrite(powfs[ipowfs].etfprep[iwvl].p2,
			   "%s/powfs%d_etfprep%d_2d",dirsetup,ipowfs,iwvl);
	    }
	}
	if(powfs[ipowfs].etfsim != powfs[ipowfs].etfsim){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(powfs[ipowfs].etfsim[iwvl].p1){
		    ccellwrite(powfs[ipowfs].etfsim[iwvl].p1, 
			       "%s/powfs%d_etfsim%d_1d",dirsetup,ipowfs,iwvl);
		}
		if(powfs[ipowfs].etfsim[iwvl].p2){
		    ccellwrite(powfs[ipowfs].etfsim[iwvl].p2,
			       "%s/powfs%d_etfsim%d_2d",dirsetup,ipowfs,iwvl);
		}
	    }
	}
    }
}
/*compute cog NEA using Monte Carlo realizations of noise*/
static void cog_nea(double *nea, dmat *ints, double cogthres, double cogoff, int ntry, 
		    rand_t *rstat, double bkgrnd, double bkgrndc, double *bkgrnd2i, double *bkgrnd2ic, double rne
		    ){
    dmat *ints2=dnew(ints->nx, ints->ny);
    double gnf[2]={0,0};
    double gny[2]={0,0};
    dcog(gnf, ints, 0, 0, cogthres, cogoff);
    seed_rand(rstat, 1);/*reset the seed each time to make dminsearch work.*/
    nea[0]=0; nea[1]=0; nea[2]=0; nea[3]=0;
    for(int i=0; i<ntry; i++){
	dcp(&ints2, ints);
	addnoise(ints2, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
	dcog(gny, ints2, 0, 0, cogthres, cogoff);
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
    double *bkgrnd2i;
    double *bkgrnd2ic;
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
    const int nsa=powfs[ipowfs].pts->nsa;
    const int ntry=500;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
    const double bkgrndc=parms->powfs[ipowfs].bkgrndc;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    int do_nea=0;
    if(parms->powfs[ipowfs].phytypesim==2){;
	powfs[ipowfs].gradphyoff=dcellnew(nwfs, 1);
    }
    rand_t rstat;
    double neaspeckle2=0;
    dcell *sanea=NULL;
    if(parms->powfs[ipowfs].phytype==2){/*need nea in reconstruction*/
	do_nea=1;
	intstat->saneaxy=dcellnew(nsa, intstat->i0->ny);
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
    intstat->cogcoeff=dcellnew(nwfs,1);
    
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(iwfs==0 || intstat->i0->ny>1){
	    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
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

	    double *restrict gx=NULL;
	    double *restrict gy=NULL;
	    if(powfs[ipowfs].gradphyoff){
		powfs[ipowfs].gradphyoff->p[iwfs]=dnew(nsa*2,1);
		gx=powfs[ipowfs].gradphyoff->p[iwfs]->p;
		gy=gx+nsa;
	    }
	    intstat->cogcoeff->p[iwfs]=dnew(2,nsa);
	    PDMAT(intstat->cogcoeff->p[iwfs], cogcoeff);
	    double g[2]={0,0};
	    for(int isa=0; isa<nsa; isa++){
		dmat *ints=intstat->i0->p[isa+iwfs*nsa];/*equivalent noise*/
		double *bkgrnd2i=(bkgrnd2 && bkgrnd2[isa])?bkgrnd2[isa]->p:NULL;
		double *bkgrnd2ic=(bkgrnd2c && bkgrnd2c[isa])?bkgrnd2c[isa]->p:NULL;

		if(parms->powfs[ipowfs].cogthres<0 && parms->powfs[ipowfs].cogoff<0){
		    /*automatically figure out the optimal thres/offset for each subaperture.*/
		    double coeff[2]={1., 1.};
		    double scale[2]={0.1,0.1};
		    cogdata_t data={ints, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne, &rstat, 100};
		    double ftol=0.0001;/*tolerance of nea^2 in pixel*/
		    int ncall=dminsearch(coeff, scale, 2, ftol, (dminsearch_fun)cogfun, &data);
		    cogcoeff[isa][0]=coeff[0];
		    cogcoeff[isa][1]=coeff[1];
		    info("isa %d:, ncall=%3d, coeff=%g %g\n", isa, ncall, coeff[0], coeff[1]);
		}else{
		    cogcoeff[isa][0]=parms->powfs[ipowfs].cogthres;
		    cogcoeff[isa][1]=parms->powfs[ipowfs].cogoff;
		}
		if(do_nea){
		    dmat *nea=dnew(2,2);
		    cog_nea(nea->p, ints, cogcoeff[isa][0], cogcoeff[isa][1], ntry, &rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
		    nea->p[0]=nea->p[0]*pixthetax*pixthetax+neaspeckle2;
		    nea->p[3]=nea->p[3]*pixthetay*pixthetay+neaspeckle2;
		    nea->p[1]=nea->p[1]*pixthetax*pixthetay;
		    nea->p[2]=nea->p[1];

		    sanea->p[iwfs]->p[isa]=nea->p[0];
		    sanea->p[iwfs]->p[isa+nsa]=nea->p[3];
		    if(srot){
			drotvecnn(&intstat->saneaxy->p[isa+nsa*iwfs], nea, srot[isa]);
		    }else{
			intstat->saneaxy->p[isa+nsa*iwfs]=dref(nea);
		    }
		    dfree(nea);
		}
		if(gx){/*gradient offset*/
		    dcog(g, ints, 0, 0, cogcoeff[isa][0], cogcoeff[isa][1]);
		    g[0]*=pixthetax;
		    g[1]*=pixthetay;
		    if(srot){
			double theta=srot[isa];
			double cx=cos(theta);
			double sx=sin(theta);
			double tmp=g[0]*cx-g[1]*sx;
			g[1]=g[0]*sx+g[1]*cx;
			g[0]=tmp;
		    }
		    gx[isa]=g[0];
		    gy[isa]=g[1];
		}
	    }
	}else{
	    if(powfs[ipowfs].gradphyoff){
		powfs[ipowfs].gradphyoff->p[iwfs]=dref(powfs[ipowfs].gradphyoff->p[0]);
	    }
	    powfs[ipowfs].intstat->cogcoeff->p[iwfs]=dref(powfs[ipowfs].intstat->cogcoeff->p[0]);
	}
    }
    if(parms->save.setup){
	if(powfs[ipowfs].gradphyoff){
	    dcellwrite(powfs[ipowfs].gradphyoff, "%s/powfs%d_gradphyoff", dirsetup, ipowfs);
	}
	if(sanea){
	    dcellwrite(sanea, "%s/powfs%d_sanea", dirsetup, ipowfs);
	}
	dcellwrite(powfs[ipowfs].intstat->cogcoeff, "%s/powfs%d_cogcoeff", dirsetup, ipowfs);
    }
    dcellfree(sanea);
    toc2("setup_powfs_cog");
}


/**
   Setup the matched filter pixel processing parameters for physical optics wfs.
*/
static void 
setup_powfs_mtch(POWFS_T *powfs,const PARMS_T *parms, int ipowfs){
    int disable_save_save=disable_save;
    disable_save=0;//temporarily disable this feature.
    long nsa=powfs[ipowfs].pts->nsa;
    if(powfs[ipowfs].intstat){
	error("Should only be called once\n");
    }
    INTSTAT_T *intstat=powfs[ipowfs].intstat=calloc(1, sizeof(INTSTAT_T));
    if(parms->load.i0){
	warning("Loading i0, gx, gy\n");
	intstat->i0=dcellread("powfs%d_i0",ipowfs);
	intstat->gx=dcellread("powfs%d_gx",ipowfs);
	intstat->gy=dcellread("powfs%d_gy",ipowfs);
    }else{
	if(parms->powfs[ipowfs].piinfile){
	    /*load psf. 1 for each wavefront sensor. */
	    info2("Using 1 sepsf for each wfs when loading sepsf\n");
	    intstat->nsepsf=parms->powfs[ipowfs].nwfs;
	    intstat->sepsf=calloc(intstat->nsepsf, sizeof(dcell*));
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs[jwfs];
		intstat->sepsf[jwfs]=dcellread("%s_wfs%d",parms->powfs[ipowfs].piinfile,iwfs);
		double pmax=dmax(intstat->sepsf[jwfs]->p[0]);
		if(intstat->sepsf[jwfs]->p[0]->p[0]>pmax*0.5){
		    error("wfs %d:  psf must have peak at center, not corner.\n", iwfs);
		}
		if(intstat->sepsf[jwfs]->nx!=nsa){
		    error("piinfile doesn't match\n");
		}
	    }
	}else if(parms->powfs[ipowfs].lo && parms->powfs[ipowfs].order<=2){
	    error("Please specify piinfile for lo order phy wfs\n");
	}else{
	    /*LGS */
	    {
		int npsfx=powfs[ipowfs].pts->nx *parms->powfs[ipowfs].embfac;
		int npsfy=npsfx;
		char fnprefix[PATH_MAX]; fnprefix[0]='\0';
		uint32_t key=0;
		if(parms->aper.fnamp){
		    char *tmp=mybasename(parms->aper.fnamp);
		    strncat(fnprefix, tmp, PATH_MAX-strlen(fnprefix)-1);
		    free(tmp);
		}else{
		    strncat(fnprefix, "noamp", PATH_MAX-strlen(fnprefix)-1);
		}
		if(powfs[ipowfs].amp_tel){
		    for(int iamp=0; iamp<powfs[ipowfs].nwfs; iamp++){
			key=dhash(powfs[ipowfs].amp_tel->p[iamp], key);
		    }
		}else{
		    key=dhash(powfs[ipowfs].amp, key);
		}
		info("powfs %d: ncpa_method=%d, opdbias=%p\n",
		     ipowfs, parms->powfs[ipowfs].ncpa_method, powfs[ipowfs].opdbias);
		if(powfs[ipowfs].opdbias && parms->powfs[ipowfs].ncpa_method==2){
		    info("Puting opdbias to key\n");
		    for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
			key=dhash(powfs[ipowfs].opdbias->p[iwfs],key);
		    }
		}
		if(key!=0){
		    char tmp2[80];
		    snprintf(tmp2,80,"_%ud",key);
		    strncat(fnprefix, tmp2, PATH_MAX-strlen(fnprefix)-1);
		}
	    
		char fnotf[PATH_MAX];
		char fnlock[PATH_MAX];
		snprintf(fnotf,PATH_MAX,"%s/.aos/otfc/",HOME);
		if(!exist(fnotf)) 
		    mymkdir("%s",fnotf);
	
		snprintf(fnotf,PATH_MAX,"%s/.aos/otfc/%s_D%g_%g_"
			 "r0_%g_L0%g_dsa%g_nsa%ld_dx1_%g_"
			 "nwvl%d_%g_embfac%d_%dx%d_SEOTF_v2",
			 HOME, fnprefix,
			 parms->aper.d,parms->aper.din, 
			 parms->atm.r0, parms->atm.l0, 
			 powfs[ipowfs].pts->dsa,nsa,
			 1./powfs[ipowfs].pts->dx, 
			 parms->powfs[ipowfs].nwvl,
			 parms->powfs[ipowfs].wvl[0]*1.e6,
			 parms->powfs[ipowfs].embfac,npsfx,npsfy);
		snprintf(fnlock, PATH_MAX, "%s.lock", fnotf);
	    retry:
		if(exist(fnlock) || !zfexist(fnotf)){/*need to create data */
		    int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
		    if(fd>=0){/*succeed */
			info2("Generating WFS OTF for %s...", fnotf);
			TIC;tic;
			genseotf(parms,powfs,ipowfs);
			toc2("done");
			ccellwritearr(intstat->otf, intstat->notf, 1, "%s", fnotf);
			close(fd);
			remove(fnlock);
		    }else{
			fd=lock_file(fnlock,1,0);/*blocking exclusive lock */
			close(fd);
			remove(fnlock);
			goto retry;
		    }
		}else{
		    long nx, ny;
		    info2("Reading WFS OTF from %s\n", fnotf);
		    intstat->otf=ccellreadarr(&nx, &ny, "%s",fnotf);
		    intstat->notf=nx*ny;
		    zftouch(fnotf);
		}
	    }
	    if(parms->powfs[ipowfs].llt){
		char fnprefix[80];
		uint32_t key=0;
		key=dhash(powfs[ipowfs].llt->amp, key);
		if(powfs[ipowfs].llt->ncpa){
		    dcell *ncpa=powfs[ipowfs].llt->ncpa;
		    long nlotf=ncpa->nx*ncpa->ny;
		    for(long ilotf=0; ilotf<nlotf; ilotf++){
			key=dhash(ncpa->p[ilotf], key);
		    }
		}
		snprintf(fnprefix,80,"SELOTF_%0x",key);
		char fnlotf[PATH_MAX];
		snprintf(fnlotf,PATH_MAX,"%s/.aos/otfc/%s_"
			 "r0_%g_L0%g_lltd%g_dx1_%g_W%g_"
			 "nwvl%d_%g_embfac%d_v2", 
			 HOME, fnprefix,
			 parms->atm.r0, parms->atm.l0, 
			 powfs[ipowfs].llt->pts->dsa,
			 1./powfs[ipowfs].llt->pts->dx,
			 parms->powfs[ipowfs].llt->widthp,
			 parms->powfs[ipowfs].nwvl,
			 parms->powfs[ipowfs].wvl[0]*1.e6,
			 parms->powfs[ipowfs].embfac);
		char fnllock[PATH_MAX];
		snprintf(fnllock, PATH_MAX, "%s.lock", fnlotf);
	    retry2:
		if(exist(fnllock) || !zfexist(fnlotf)){/*need to create data */
		    int fd2=lock_file(fnllock, 0, 0);/*nonblocking exclusive lock */
		    if(fd2>=0){/*succeed */
			info2("Generating WFS LLT OTF for %s\n", fnlotf);
			genselotf(parms,powfs,ipowfs);
			ccellwrite(intstat->lotf, "%s",fnlotf);
			close(fd2);
			remove(fnllock);
		    }else{
			fd2=lock_file(fnllock, 1, 0);/*blocking, exclusive */
			close(fd2);
			remove(fnllock);
			goto retry2;
		    }
		    intstat->lotf=ccellread("%s",fnlotf);
		    zftouch(fnlotf);
		}else{
		    intstat->lotf=ccellread("%s",fnlotf);
		    zftouch(fnlotf);
		    info2("Reading WFS LLT OTF from %s\n", fnlotf);
		}
		int nwvl=intstat->lotf->nx;
		PCCELL(intstat->lotf, lotf);
		int nlpsf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
		cmat *psfhat=cnew(nlpsf, nlpsf);
		dmat *psf=dnew(nlpsf, nlpsf);
		cellarr *lltpsfsave=NULL;
		if(parms->save.setup){
		    lltpsfsave=cellarr_init(nwvl, intstat->lotf->ny, "%s/powfs%d_llt_psf", dirsetup, ipowfs);
		}
		cfft2plan(psfhat, 1);
		for(int illt=0; illt<intstat->lotf->ny; illt++){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			const double dx=powfs[ipowfs].llt->pts->dx;
			const double wvl=parms->powfs[ipowfs].wvl[iwvl];
			const double dpsf=wvl/(nlpsf*dx)*206265.;
			ccp(&psfhat, lotf[illt][iwvl]);
			cfftshift(psfhat);
			cfft2i(psfhat, 1);
			cfftshift(psfhat);
			creal2d(&psf, 0, psfhat, 1);
			info2("illt %d, iwvl %d has FWHM of %g\"\n",
			      illt, iwvl, sqrt(4.*(double)dfwhm(psf)/M_PI)*dpsf);
			if(lltpsfsave) cellarr_dmat(lltpsfsave, illt*nwvl+iwvl, psf);
		    }
		}
		if(lltpsfsave) cellarr_close(lltpsfsave);
		cfree(psfhat);
		dfree(psf);
	    }
	    /*Generating short exposure images. */
	    gensepsf(parms,powfs,ipowfs);
	    if(parms->save.setup>1 && intstat){
		dcellwrite(intstat->sepsf[0], "%s/powfs%d_sepsf",dirsetup,ipowfs);
	    }
	    /*Free short exposure otf. */
	    if(!parms->sim.ncpa_calib){
		ccellfree(intstat->lotf);
	    }
	    ccellfreearr(intstat->otf, intstat->notf);
	}
	/*generate short exposure i0,gx,gy from psf. */
	gensei(parms,powfs,ipowfs);
	dcellfreearr(intstat->sepsf,intstat->nsepsf);
	if(parms->save.setup){
	    dcellwrite(intstat->i0,"%s/powfs%d_i0",dirsetup,ipowfs);
	    dcellwrite(intstat->gx,"%s/powfs%d_gx",dirsetup,ipowfs);
	    dcellwrite(intstat->gy,"%s/powfs%d_gy",dirsetup,ipowfs);
	}
    
	/*Remove OTFs that are older than 30 days. */
	char *dirotf=stradd(HOME, "/.aos/otfc", NULL);
	remove_file_older(dirotf, 30*24*3600);
	free(dirotf);
    }
    /*Generating Matched filter */
    if(parms->powfs[ipowfs].phytype==1 || parms->powfs[ipowfs].phytypesim==1){
	genmtch(parms,powfs,ipowfs);
	if(parms->save.setup){
	    dcellwrite(powfs[ipowfs].intstat->mtche, "%s/powfs%d_mtche",dirsetup,ipowfs);
	}
    }
    if(parms->powfs[ipowfs].phytype==2 || parms->powfs[ipowfs].phytypesim==2){
	setup_powfs_cog(parms, powfs, ipowfs);
    }
    intstat->saneaxyl=dcellnew(intstat->saneaxy->nx, intstat->saneaxy->ny);//cholesky decomposition
    intstat->saneaixy=dcellnew(intstat->saneaxy->nx, intstat->saneaxy->ny);//inverse
    for(int i=0; i<intstat->saneaxy->nx*intstat->saneaxy->ny; i++){
	intstat->saneaxyl->p[i]=dchol(intstat->saneaxy->p[i]);
	intstat->saneaixy->p[i]=dinvspd(intstat->saneaxy->p[i]);
    }
    if(parms->save.setup){
	dcellwrite(intstat->saneaxy, "%s/powfs%d_saneaxy",dirsetup,ipowfs);
    }
    if(parms->powfs[ipowfs].neaphy){/*use physical optics nea for geom grad*/
	powfs[ipowfs].neasim=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    dmat **sanea=NULL;
	    /*saneaxyl is the cholesky decomposition of the Cnn. */
	    if(powfs[ipowfs].intstat->saneaxyl->ny==1){
		sanea=powfs[ipowfs].intstat->saneaxyl->p;
	    }else{
		sanea=powfs[ipowfs].intstat->saneaxyl->p+jwfs;
	    }
	    dmat *nea=dnew(nsa,3);
	    PDMAT(nea, pnea);
	    for(int isa=0; isa<nsa; isa++){
		pnea[0][isa]=sanea[isa]->p[0];
		pnea[1][isa]=sanea[isa]->p[3];
		pnea[2][isa]=sanea[isa]->p[1];
	    }
	    powfs[ipowfs].neasim->p[jwfs]=nea;
	}
	if(parms->save.setup){
	    dcellwrite(powfs[ipowfs].neasim,"%s/powfs%d_neasim", dirsetup, ipowfs);
	}
    }
    dcellfree(powfs[ipowfs].intstat->i0);
    dcellfree(powfs[ipowfs].intstat->gx);
    dcellfree(powfs[ipowfs].intstat->gy);
    disable_save=disable_save_save;//put it back.
}
/*
  Setup gradient offset for calibration.
*/
void setup_powfs_calib(const PARMS_T *parms, POWFS_T *powfs, loc_t **aloc, dcell *dm_ncpa){
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	dcellcp(&powfs[ipowfs].opdbias, powfs[ipowfs].opdadd);
	if(aloc && dm_ncpa){
	    for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
		int iwfs0=parms->powfs[ipowfs].wfs[iwfs];
		double hs=parms->wfs[iwfs].hs;
		double thetax=parms->wfs[iwfs0].thetax;
		double thetay=parms->wfs[iwfs0].thetay;
		if(!powfs[ipowfs].opdbias){
		    powfs[ipowfs].opdbias=dcellnew(parms->powfs[ipowfs].nwfs, 1);
		}
		if(!powfs[ipowfs].opdbias->p[iwfs]){
		    powfs[ipowfs].opdbias->p[iwfs]=dnew(powfs[ipowfs].npts, 1);
		}
		for(int idm=0; idm<parms->ndm; idm++){
		    if(!dm_ncpa->p[idm]) continue;
		    double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
		    double scale=1.-ht/hs;
		    double dispx=ht*thetax;
		    double dispy=ht*thetay;
		    if(parms->dm[idm].cubic){
			prop_nongrid_pts_cubic(aloc[idm], dm_ncpa->p[idm]->p, 
					       powfs[ipowfs].pts, NULL, powfs[ipowfs].opdbias->p[iwfs]->p, 
					       -1, dispx, dispy, scale, parms->dm[idm].iac, 0, 0);
		    }else{
			prop_nongrid_pts(aloc[idm], dm_ncpa->p[idm]->p, 
					 powfs[ipowfs].pts, NULL, powfs[ipowfs].opdbias->p[iwfs]->p, 
					 -1, dispx, dispy, scale, 0, 0);
		    }
		}
	    }
	}
	//if opdadd is null, but dm_ncpa is not, there will be opdbias.
	if(powfs[ipowfs].opdbias){
	    if(parms->sim.ncpa_ttr){
		/*remove average tilt from opdbias and same amount from
		  opdadd. Does not need to be very accurate.*/
		dmat *mcc=loc_mcc_ptt(powfs[ipowfs].loc, powfs[ipowfs].amp->p);
		dinvspd_inplace(mcc);
		for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
		    double ptt[3]={0,0,0};
		    loc_calc_ptt(NULL, ptt, powfs[ipowfs].loc, 1./mcc->p[0], mcc, 
				 powfs[ipowfs].amp->p, powfs[ipowfs].opdbias->p[iwfs]->p);
		    loc_remove_ptt(powfs[ipowfs].opdbias->p[iwfs]->p, ptt, powfs[ipowfs].loc);
		    loc_remove_ptt(powfs[ipowfs].opdadd->p[iwfs]->p, ptt, powfs[ipowfs].loc);
		}
		dfree(mcc);
	    }
	    if(parms->powfs[ipowfs].ncpa_method==1){
		if(!powfs[ipowfs].gradoff){
		    powfs[ipowfs].gradoff=dcellnew(parms->powfs[ipowfs].nwfs,1);
		}
		for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
		    if(powfs[ipowfs].opdbias->p[iwfs]){
			double *realamp=powfs[ipowfs].realamp->p[iwfs]->p;
			if(parms->powfs[ipowfs].gtype_sim==1){
			    pts_ztilt(&powfs[ipowfs].gradoff->p[iwfs], powfs[ipowfs].pts,
				      powfs[ipowfs].saimcc[powfs[ipowfs].nsaimcc>1?iwfs:0], 
				      realamp, powfs[ipowfs].opdbias->p[iwfs]->p);
			}else{
			    spmulmat(&powfs[ipowfs].gradoff->p[iwfs],adpind(powfs[ipowfs].GS0, iwfs),
				     powfs[ipowfs].opdbias->p[iwfs],1);
			}
		    }
		}
		if(powfs[ipowfs].gradoff){
		    dcellwrite(powfs[ipowfs].gradoff, "powfs%d_gradoff", ipowfs);
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
    POWFS_T *powfs=calloc(parms->npowfs, sizeof(POWFS_T));
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	info2("\n\033[0;32mSetting up powfs %d geom\033[0;0m\n\n", ipowfs);
	setup_powfs_geom(powfs,parms,aper,ipowfs);
        setup_powfs_grad(powfs,parms,ipowfs);
    }
    toc("setup_powfs_init");
    return powfs;
}
void setup_powfs_phy(const PARMS_T *parms, POWFS_T *powfs){
    TIC;tic;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	info2("\n\033[0;32mSetting up powfs %d PO WFS\033[0;0m\n\n", ipowfs);
	if(TEST_POWFS||parms->powfs[ipowfs].usephy
	   ||parms->powfs[ipowfs].psfout
	   ||parms->powfs[ipowfs].pistatout
	   ||parms->powfs[ipowfs].neaphy){
	    /*We have physical optics. setup necessary struct */
	    setup_powfs_prep_phy(powfs,parms,ipowfs);
	    setup_powfs_dtf(powfs,parms,ipowfs);
	    if(parms->powfs[ipowfs].llt){
		/*prepare Laser launch telescope. */
		setup_powfs_sodium(powfs,parms,ipowfs);/*read sodium profile and smooth it */
		setup_powfs_etf(powfs,parms,ipowfs,0,0);/*etf for prep */
		setup_powfs_etf(powfs,parms,ipowfs,1,0);/*etf for sim */
		setup_powfs_llt(powfs,parms,ipowfs);
	    }
	}
	if(parms->powfs[ipowfs].llt){
	    /*If there is LLT, setup the extra focus term if needed. */
	    setup_powfs_focus(powfs,parms,ipowfs);
	}
	if(TEST_POWFS||parms->powfs[ipowfs].usephy || parms->powfs[ipowfs].neaphy){
	    setup_powfs_mtch(powfs,parms,ipowfs);
	}
    }/*ipowfs */
    toc("setup_powfs_phy");
}
/**
   Free all parameters of powfs at the end of simulation.
*/
void free_powfs(const PARMS_T *parms, POWFS_T *powfs){
    int ipowfs;
    for(ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	free_powfs_geom(powfs, parms, ipowfs);
	free_powfs_dtf(powfs, parms, ipowfs);
	spcellfree(powfs[ipowfs].GS0);
	dcellfree(powfs[ipowfs].neasim);
	if(powfs[ipowfs].intstat){
	    INTSTAT_T *intstat=powfs[ipowfs].intstat;
	    ccellfreearr(intstat->fotf, intstat->nsepsf);
	    dcellfree(intstat->mtche);
	    dcellfree(intstat->saneaxyl);
	    dcellfree(intstat->saneaxy);
	    dcellfree(intstat->saneaixy);
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

	dcellfreearr(powfs[ipowfs].saimcc, powfs[ipowfs].nsaimcc);
	if(powfs[ipowfs].llt){
	    ptsfree(powfs[ipowfs].llt->pts);
	    dfree(powfs[ipowfs].llt->amp);
	    locfree(powfs[ipowfs].llt->loc);
	    dcellfree(powfs[ipowfs].llt->mcc);
	    dcellfree(powfs[ipowfs].llt->imcc);
	    dcellfree(powfs[ipowfs].llt->ncpa);
	    free(powfs[ipowfs].llt);
	}
	dfree(powfs[ipowfs].sodium);
	if(powfs[ipowfs].etfprep){
	    if(powfs[ipowfs].etfprep!=powfs[ipowfs].etfsim){
		for(int iwvl=0;iwvl<parms->powfs[ipowfs].nwvl;iwvl++){
		    ccellfree(powfs[ipowfs].etfsim[iwvl].p1);
		    ccellfree(powfs[ipowfs].etfsim[iwvl].p2);
		}
		free(powfs[ipowfs].etfsim);
	    }
	    for(int iwvl=0;iwvl<parms->powfs[ipowfs].nwvl;iwvl++){
		ccellfree(powfs[ipowfs].etfprep[iwvl].p1);
		ccellfree(powfs[ipowfs].etfprep[iwvl].p2);
	    }
	    free(powfs[ipowfs].etfprep);
	}
	dcellfree(powfs[ipowfs].opdadd);
	dcellfree(powfs[ipowfs].opdbias);
	dcellfree(powfs[ipowfs].gradoff);
    }
    free(powfs);
}
