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

#include <math.h>
#include <string.h>
#include <alloca.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include "maos.h"
TIC;
#include "setup_powfs.h"
#include "mtch.h"
#include "genseotf.h"
#define TEST_POWFS 0
#define MOVES(p,i,j) p[i]=p[j]
#define MOVED(p,n,i,j) memcpy(p+n*i, p+n*j, sizeof(double)*n)
#define MOVEPTS(pts,count,isa)			\
    if(pts){\
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
   much space.  */

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
    if(powfs[ipowfs].nlocm){
	for(int ilocm=0; ilocm<powfs[ipowfs].nlocm; ilocm++){
	    locfree(powfs[ipowfs].locm[ilocm]);
	}
	free(powfs[ipowfs].locm);powfs[ipowfs].locm=NULL;
	dcellfree(powfs[ipowfs].ampm);
	dcellfree(powfs[ipowfs].saam);
    }
    free(powfs[ipowfs].realamp);
    free(powfs[ipowfs].realsaa);
    dfree(powfs[ipowfs].sumamp);
    dfree(powfs[ipowfs].sumamp2);
    free(powfs[ipowfs].misreg);
    dcellfree(powfs[ipowfs].saam);
    dcellfree(powfs[ipowfs].ampm);
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
static void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, double misreg[2], int iwfs){
    long nloc=loc->nloc;
    dmat *ampmask=dnew(nloc, 1);
    double ht=parms->atm.hmax*0.7;
    int ipowfs=parms->wfs[iwfs].powfs;
    for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
	int jpowfs=parms->wfs[jwfs].powfs;
	if(parms->powfs[jpowfs].lo) continue;
	double r=parms->aper.d*0.5*(1.-ht/parms->powfs[jpowfs].hs)/(1.-ht/parms->powfs[ipowfs].hs);
	double sx=(parms->wfs[jwfs].thetax-parms->wfs[iwfs].thetax)*ht-misreg[0];
	double sy=(parms->wfs[jwfs].thetay-parms->wfs[iwfs].thetay)*ht-misreg[1];
	loccircle(ampmask->p, loc, sx,sy, r, 1);
    }
    for(int i=0; i<nloc; i++){
	if(ampmask->p[i]<0.5) amp->p[i]=0;
    }
    dfree(ampmask);
}
static void 
setup_powfs_geom(POWFS_T *powfs, const PARMS_T *parms, 
		 APER_T *aper, int ipowfs){
    free_powfs_geom(powfs, parms, ipowfs);
    double offset,dxoffset;
    int count;
 
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
    /*The threashold for normalized area (by areafulli) to keep subaperture. */
    const double thresarea=parms->powfs[ipowfs].saat;
    /*Offset of the coordinate of the center most subaperture from the center. */
    if(order & 1){/*odd */
	offset = -0.5;
    }else{
	offset = 0.0;
    }
    /*r2max: Maximum distance^2 from the center to keep a subaperture */
    double r2max=pow(order/2+0.5, 2);
    double r2min=dxsa<parms->aper.din?pow(parms->aper.din/dxsa/2-0.5,2):-1;
    /*the lower left *grid* coordinate of the subaperture */
    
    /*The coordinate of the subaperture (lower left coordinate) */
    powfs[ipowfs].saloc=locnew(order*order, dxsa);
    /*Number of OPD pixels in 1 direction in each
      subaperture. Make it even to do fft.*/
    int nx = 2*(int)round(0.5*dxsa/parms->powfs[ipowfs].dx);
    const double dx=dxsa/nx;/*adjust dx. */
    powfs[ipowfs].pts=ptsnew(order*order, dxsa, nx, dx);/*calloc(1, sizeof(pts_t)); */

    if(fabs(parms->powfs[ipowfs].dx-powfs[ipowfs].pts->dx)>EPS)
	warning("Adjusting dx from %g to %g\n",
		parms->powfs[ipowfs].dx,powfs[ipowfs].pts->dx);
    if(fabs(dxsa - nx * powfs[ipowfs].pts->dx)>EPS){
	warning("nx=%d,dsa=%f,dx=%f for not agree\n", nx, dxsa, dx);
    }
    dxoffset=dx*0.5;//Always keep points inside subaperture for simulation.

    /*if(parms->dbg.dxonedge){
	warning2("Put points on the edge of subapertures\n");
	dxoffset=0;
	}else{*/
    //}
    info2("There are %d points in each subapeture of %gm.\n", nx, dxsa);
    const int nxsa=nx*nx;/*Total Number of OPD points. */
    count = 0;
    /*Collect all the subapertures that are within the allowed radius*/
    for(int j=-order/2; j<=(order-1)/2; j++){
	for(int i=-order/2; i<=(order-1)/2; i++){
	    /*subaperture distance from the center*/
	    double r2=pow(i+0.5+offset, 2)+pow(j+0.5+offset, 2);
	    if(r2 <= r2max && r2 >=r2min){
		powfs[ipowfs].pts->origx[count]=
		    ((double)i+offset)*dxsa+dxoffset;
		powfs[ipowfs].pts->origy[count]=
		    ((double)j+offset)*dxsa+dxoffset;
		powfs[ipowfs].saloc->locx[count]=
		    ((double)i+offset)*dxsa;
		powfs[ipowfs].saloc->locy[count]=
		    ((double)j+offset)*dxsa;
		count++;
	    }
	}
    }
    /*Calculate the amplitude for each subaperture OPD point by ray tracing from
      the pupil amplitude map.*/
    powfs[ipowfs].pts->nsa=count;
    powfs[ipowfs].loc=pts2loc(powfs[ipowfs].pts);
    /*The assumed amp. */
    powfs[ipowfs].amp=mkwfsamp(powfs[ipowfs].loc, aper->ampground, 
			       0,0, parms->aper.d, parms->aper.din);
    powfs[ipowfs].saa=wfsamp2saa(powfs[ipowfs].amp, nxsa);
    if(parms->dbg.dxonedge){
	//Create another set of loc/amp that can be used to build GP.
	map_t *map=create_metapupil_wrap(parms, 0, dx, 0, 0, 0, 0, 0, 0);
	powfs[ipowfs].gloc=map2loc(map); mapfree(map);
	powfs[ipowfs].gamp=mkwfsamp(powfs[ipowfs].gloc, aper->ampground, 
				    0,0, parms->aper.d, parms->aper.din);
    }else{
	powfs[ipowfs].gloc=powfs[ipowfs].loc;
	powfs[ipowfs].gamp=powfs[ipowfs].amp;
    }
    if(parms->dbg.pupmask && parms->powfs[ipowfs].lo){
	double misreg[2]={0,0};
	if(parms->powfs[ipowfs].nwfs>1){
	    error("dbg.pupmask=1, only powfs can only have 1 wfs.\n");
	}
	int iwfs=parms->powfs[ipowfs].wfs[0];
	wfspupmask(parms, powfs[ipowfs].loc, powfs[ipowfs].amp, misreg, iwfs);
	if(parms->dbg.dxonedge){
	    wfspupmask(parms, powfs[ipowfs].gloc, powfs[ipowfs].gamp, misreg, iwfs);
	}
    }
    powfs[ipowfs].misreg=calloc(2*parms->powfs[ipowfs].nwfs, sizeof(double));
    if(parms->powfs[ipowfs].misreg){
	/*Create misregistered coordinate and amplitude map to do physical
	  optics wavefront sensing. the un-misregistered coordinate and
	  amplitude map is still used to do the reconstructor because the RTC is
	  not aware of the misregistration. The misregistered coordinates are
	  only used for raytracing that is what happens in reality that we are
	  not aware of. */
	warning("Creating misregistration for WFS.\n");
	dcell *misreg=dcellread("%s",parms->powfs[ipowfs].misreg); 
	int nlocm=0;
	if(misreg->nx!=2)
	    error("%s is in wrong format\n",parms->powfs[ipowfs].misreg);
	if(misreg->ny==1){
	    nlocm=1;
	}else if(misreg->ny==parms->powfs[ipowfs].nwfs){
	    nlocm=parms->powfs[ipowfs].nwfs;
	}else{
	    error("%s is in wrong format\n",parms->powfs[ipowfs].misreg);
	}
	powfs[ipowfs].nlocm=nlocm;
	powfs[ipowfs].locm=calloc(nlocm, sizeof(loc_t*));
	powfs[ipowfs].saam=dcellnew(nlocm, 1);
	powfs[ipowfs].ampm=dcellnew(nlocm, 1);
	PDCELL(misreg,pmisreg);
	for(int ilocm=0; ilocm<nlocm; ilocm++){
	    /*Transforms loc using misregistration information. The do ray
	      tracing to determine the misregistered ampltidue map*/
	    double *shiftxy=NULL;
	    /*may return NUL; if it is purely transform.*/
	    powfs[ipowfs].locm[ilocm]=loctransform(powfs[ipowfs].loc, &shiftxy, pmisreg[ilocm]);
	    powfs[ipowfs].misreg[ilocm][0]=shiftxy[0];
	    powfs[ipowfs].misreg[ilocm][1]=shiftxy[1];
	    powfs[ipowfs].saam->p[ilocm]=dnew(powfs[ipowfs].pts->nsa, 1);
	    powfs[ipowfs].ampm->p[ilocm]=dnew(powfs[ipowfs].loc->nloc, 1);
	    loc_t *loc=powfs[ipowfs].locm[ilocm]?powfs[ipowfs].locm[ilocm]:powfs[ipowfs].loc;
	    powfs[ipowfs].ampm->p[ilocm]=mkwfsamp(loc, aper->ampground, 
						  shiftxy[0]-parms->aper.misreg[0], 
						  shiftxy[1]-parms->aper.misreg[1],
						  parms->aper.d, parms->aper.din);
	    if(parms->dbg.pupmask && parms->powfs[ipowfs].lo){
		//don't use aper.misreg since the pupil mask is in wfs.
		if(parms->powfs[ipowfs].nwfs>1){
		    error("dbg.pupmask=1, only powfs can only have 1 wfs.\n");
		}
		int iwfs=parms->powfs[ipowfs].wfs[0];
		wfspupmask(parms, powfs[ipowfs].locm[ilocm], powfs[ipowfs].ampm->p[ilocm], shiftxy, iwfs);
	    }
	    powfs[ipowfs].saam->p[ilocm]=wfsamp2saa(powfs[ipowfs].ampm->p[ilocm], nxsa);
	}
        for(int ilocm=nlocm; nlocm==1 && ilocm<parms->powfs[ipowfs].nwfs; ilocm++){
	    powfs[ipowfs].misreg[ilocm][0]=powfs[ipowfs].misreg[0][0];
	    powfs[ipowfs].misreg[ilocm][1]=powfs[ipowfs].misreg[0][1];
	}
    }else if(parms->aper.ismisreg){/*do not need powfs[ipowfs].locm. */
	powfs[ipowfs].nlocm=1;
	powfs[ipowfs].locm=calloc(1, sizeof(loc_t*));
	powfs[ipowfs].ampm=dcellnew(1,1);
	powfs[ipowfs].saam=dcellnew(1,1);
	powfs[ipowfs].ampm->p[0]=mkwfsamp(powfs[ipowfs].loc, aper->ampground, 
					  -parms->aper.misreg[0], -parms->aper.misreg[1], 
					  parms->aper.d, parms->aper.din);
	powfs[ipowfs].saam->p[0]=wfsamp2saa(powfs[ipowfs].ampm->p[0], nxsa);
    }/*if misreg */
    count=0;
    /*Go over all the subapertures, calculate the normalized
      subaperture illumination area and remove all that are below
      the are threshold*/
    const int nlocm=powfs[ipowfs].nlocm;
    dmat *saa=NULL;
    switch(nlocm){
    case 0:
	saa=ddup(powfs[ipowfs].saa);break;
    case 1:
	saa=ddup(powfs[ipowfs].saam->p[0]);break;
    default:
	{/*We take the average now to avoid different geometry in different LGS WFS.  */
	    const double scale=1./(double)nlocm;
	    saa=dnew(powfs[ipowfs].pts->nsa, 1);
	    for(int ilocm=0; ilocm<nlocm; ilocm++){
		dadd(&saa, 1, powfs[ipowfs].saam->p[ilocm], scale);
	    }
	}
    }
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
		for(int ilocm=0; ilocm<nlocm; ilocm++){
		    MOVES(powfs[ipowfs].saam->p[ilocm]->p, count, isa);
		    MOVEDLOC(powfs[ipowfs].locm[ilocm], nxsa, count, isa);
		    MOVED(powfs[ipowfs].ampm->p[ilocm]->p,nxsa, count, isa);
		}
	    }
	    count++;
	}
    }
    /*area was already normalized against square subaperture.  to scale the
      physical optics i0. After FFT, the PSF sum to 1 for square subapertures,
      but sum to PI/4 for circular TT or TTF subapertures. maxarea is 1 for
      square subapertures, and PI/4 for TT and TTF WFS due to cut off by
      pupil. for square subapertures, areascale=1, no effect. for TT or TTF WFS,
      areascale is 4/PI, which scales the PSF from sum to PI/4 to sum to 1
      again. Used in wfsints.

      We normalize area against 1 for square or PI/4 for arc subapertures. The
      area is not used in physical optics. It is used to scale geometric
      sanea.*/
    double maxarea=dmax(saa);
    if(maxarea>1+EPS){
	error("The area maxes to %g, which should be leq 1\n", maxarea);
    }
    dfree(saa);
    powfs[ipowfs].areascale=1./maxarea;
    if(fabs(areafulli-1)>EPS){
	dscale(powfs[ipowfs].saa, areafulli);
	dcellscale(powfs[ipowfs].saam, areafulli);
    }
    if(count==0){
	error("there are no subapertures above threshold.\n");
    }
    powfs[ipowfs].npts = count*nxsa;
    powfs[ipowfs].nthread=count<parms->sim.nthread?count:parms->sim.nthread;
    ptsresize(powfs[ipowfs].pts, count);
    locresize(powfs[ipowfs].saloc, count);
    locresize(powfs[ipowfs].loc, count*nxsa);
    dresize(powfs[ipowfs].saa, count, 1);
    dresize(powfs[ipowfs].amp, count*nxsa, 1);
    for(int ilocm=0; ilocm<nlocm; ilocm++){
	locresize(powfs[ipowfs].locm[ilocm], count*nxsa);
	dresize(powfs[ipowfs].ampm->p[ilocm], count*nxsa, 1);
	dresize(powfs[ipowfs].saam->p[ilocm], count, 1);
    }	
    powfs[ipowfs].realamp=calloc(parms->powfs[ipowfs].nwfs, sizeof(double*));
    powfs[ipowfs].realsaa=calloc(parms->powfs[ipowfs].nwfs, sizeof(double*));
    powfs[ipowfs].sumamp=dnew(parms->powfs[ipowfs].nwfs, 1);
    powfs[ipowfs].sumamp2=dnew(parms->powfs[ipowfs].nwfs, 1);
    for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
	double *realamp, *realsaa;
	if(powfs[ipowfs].locm){
	    int ilocm=(powfs[ipowfs].nlocm>1)?iwfs:0;
	    realamp=powfs[ipowfs].ampm->p[ilocm]->p;
	    realsaa=powfs[ipowfs].saam->p[ilocm]->p;
	}else{
	    realamp=powfs[ipowfs].amp->p;
	    realsaa=powfs[ipowfs].saa->p;
	}
	double sumamp2=0;
	double sumamp=0;
	for(long i=0; i<powfs[ipowfs].loc->nloc; i++){
	    sumamp2+=realamp[i]*realamp[i];
	    sumamp+=realamp[i];
	}
	powfs[ipowfs].sumamp2->p[iwfs]=sumamp2;
	powfs[ipowfs].sumamp->p[iwfs]=sumamp;
	powfs[ipowfs].realamp[iwfs]=realamp;
	powfs[ipowfs].realsaa[iwfs]=realsaa;
    }
   
    if(parms->plot.setup){
	drawopd("amp", powfs[ipowfs].loc, powfs[ipowfs].amp->p,NULL,
		"WFS Amplitude Map","x (m)","y (m)","powfs %d", ipowfs);
	for(int ilocm=0; ilocm<nlocm; ilocm++){
	    drawopd("ampm", powfs[ipowfs].loc, powfs[ipowfs].ampm->p[ilocm]->p,NULL,
		    "WFS Amplitude Map","x (m)","y (m)","powfs %d", ipowfs);
	}
    }
  
    if(parms->save.setup){
	locwrite((loc_t*)powfs[ipowfs].pts, "%s/powfs%d_pts",dirsetup,ipowfs);
	locwrite(powfs[ipowfs].saloc, "%s/powfs%d_saloc",dirsetup,ipowfs); 
	dwrite(powfs[ipowfs].saa,"%s/powfs%d_saa", dirsetup,ipowfs);
	locwrite(powfs[ipowfs].loc,"%s/powfs%d_loc",dirsetup,ipowfs);
	dwrite(powfs[ipowfs].amp, "%s/powfs%d_amp", dirsetup,ipowfs);
	if(nlocm){
	    dcellwrite(powfs[ipowfs].saam,"%s/powfs%d_saam", dirsetup,ipowfs);
	    dcellwrite(powfs[ipowfs].ampm,"%s/powfs%d_ampm", dirsetup,ipowfs);
	}
	if(parms->powfs[ipowfs].misreg){
	    for(int ilocm=0; ilocm<powfs[ipowfs].nlocm; ilocm++){
		locwrite(powfs[ipowfs].locm[ilocm],"%s/powfs%d_locm%d", dirsetup,ipowfs,ilocm);
	    }
	}
    }
}
/**
   Creating geometric wavefront gradient operator GS0 and ZA0 from WFS OPD to
   subaperture grads. Use loc, and ampm. Do not locm. locm is only used to do
   ray tracing.  */
static void 
setup_powfs_grad(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(parms->powfs[ipowfs].gtype_recon==0 ||parms->powfs[ipowfs].gtype_sim==0){
	spcellfree(powfs[ipowfs].GS0);
	/*Setting up every gradient tilt (g-ztilt) */
	if(parms->load.GS0){
	    powfs[ipowfs].GS0=spcellread("powfs%d_GS0",ipowfs);
	    if(powfs[ipowfs].ampm && powfs[ipowfs].ampm->nx>1){
		assert(powfs[ipowfs].GS0->nx==powfs[ipowfs].ampm->nx);
	    }else{
		assert(powfs[ipowfs].GS0->nx==1);
	    }
	}else{
	    double displace[2]={0,0};
	    /*This mkg takes about 5 seconds. */
	    if(powfs[ipowfs].ampm && powfs[ipowfs].ampm->nx>1){
		powfs[ipowfs].GS0=spcellnew(powfs[ipowfs].ampm->nx, 1);
	    }else{
		powfs[ipowfs].GS0=spcellnew(1, 1);
	    }
	    for(int iwfs=0; iwfs<powfs[ipowfs].GS0->nx; iwfs++){
		powfs[ipowfs].GS0->p[iwfs]=mkg(powfs[ipowfs].loc, 
					       powfs[ipowfs].loc,
					       powfs[ipowfs].realamp[iwfs],
					       powfs[ipowfs].saloc,
					       1, 1, displace, 1);
	    }
	}
    }
    if(parms->powfs[ipowfs].gtype_recon==1 ||parms->powfs[ipowfs].gtype_sim==1){
	/*setting up zernike best fit (ztilt) inv(M'*W*M). good for NGS. */
	if(parms->powfs[ipowfs].order>4) 
	    warning("Ztilt for high order wfs is not good");
	powfs[ipowfs].nsaimcc=MAX(1,powfs[ipowfs].nlocm);
	int nsaimcc=powfs[ipowfs].nsaimcc;
	dcellfreearr(powfs[ipowfs].saimcc, nsaimcc);
	powfs[ipowfs].saimcc=calloc(nsaimcc, sizeof(dcell*));
	for(int imcc=0; imcc<nsaimcc; imcc++){
	    dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].realamp[imcc]);
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
		    warning("wfs %d: NEA=%g mas, read from file is too small. Unit error?\n",iwfs, neamax*206265000);
		}else if(neamax>parms->powfs[ipowfs].pixtheta){
		    warning("wfs %d: NEA=%g mas, read from file is too big. Unit error?\n",iwfs, neamax*206265000);
		}
	    }
	    /*Scale by dtrat */
	    dscale(nea, 1./sqrt(parms->powfs[ipowfs].dtrat));
	    /*Scale by normalized subaperture area. */
	    double *saa=powfs[ipowfs].realsaa[jwfs];
	    const int dl=parms->powfs[ipowfs].dl; /*diffraction limited. */
	    for(long isa=0; isa<nsa; isa++){
		/*scale nea by sqrt(1/area). (seeing limited) */
		/*scale nea by 1/area if diffraction limited (NGS) */
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
	info2("ncompx=%d, ncompy=%d\n", ncompx,ncompy);
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
   Set up NCPA maps for each WFS. Load NCPA maps and ray trace to grid of loc */
static void
setup_powfs_ncpa(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    int do_rot=(fabs(rot)>1.e-10);
    const char *fn_ncpa = parms->powfs[ipowfs].ncpa;
    if(fn_ncpa){/*Always make nwfs NCPA's */
	int nncpa_in;
	map_t **ncpa=maparrread(&nncpa_in, "%s",fn_ncpa);
	const int nwfs=parms->powfs[ipowfs].nwfs;
	if(nncpa_in != 1 && nncpa_in != nwfs){
	    error("powfs[%d].ncpa is in wrong format. nncpa_in=%d, nwfs=%d\n",
		  ipowfs,nncpa_in, nwfs);
	}

	int incpa_mul=0;
	if(nncpa_in==nwfs){
	    incpa_mul=1;
	}
	info2("Creating NCPA\n");
	dcellfree(powfs[ipowfs].ncpa);
	powfs[ipowfs].ncpa=dcellnew(nwfs,1);
	
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ilocm=0;
	    if(powfs[ipowfs].locm && powfs[ipowfs].nlocm>1){/*misregistration. */
		ilocm=parms->powfs[ipowfs].wfsind[iwfs];
	    }
	    int incpa_in=incpa_mul*iwfs;
	    info2("iwfs=%d, incpa_in=%d\n", iwfs, incpa_in);

	    if(!powfs[ipowfs].ncpa->p[iwfs]){
		powfs[ipowfs].ncpa->p[iwfs]=dnew(powfs[ipowfs].npts,1);
	    }
	    loc_t *locwfsin, *locwfs;
	    if(powfs[ipowfs].locm){
		locwfsin=powfs[ipowfs].locm[ilocm];
	    }else{
		locwfsin=powfs[ipowfs].loc;
	    }
	    if(do_rot){
		info2("Rotating telescope pupil\n");
		locwfs=locdup(locwfsin);
		locrot(locwfs,rot);
	    }else{
		locwfs=locwfsin;
	    }
	    double hl=ncpa[incpa_in]->h;
	    double hs=parms->powfs[ipowfs].hs;
	    const double scale=1.-hl/hs;
	    const double displacex=parms->wfs[iwfs].thetax*hl;
	    const double displacey=parms->wfs[iwfs].thetay*hl;
	    if(ncpa[incpa_in]->nx==1 || ncpa[incpa_in]->ny==1){
		error("ncpa is a vector: Invalid format\n");
	    }
	    prop_grid(ncpa[incpa_in], locwfs, NULL, powfs[ipowfs].ncpa->p[iwfs]->p,
		      1, displacex, displacey, scale, 1, 0, 0);
	    if(do_rot){
		locfree(locwfs);
	    }
	}/*iwfs; */
	maparrfree(ncpa,nncpa_in);
	if(parms->powfs[ipowfs].ncpa_method<0 || parms->powfs[ipowfs].ncpa_method>2){
	    error("Invalid ncpa_method=%d\n", parms->powfs[ipowfs].ncpa_method);
	}
	/*
	  ncpa_method:
	  0: do not calibrate
	  1: use gradient offset
	  2: use offset in matched filter.
	 */
	if(parms->powfs[ipowfs].ncpa_method==1){
	    info2("Calculating gradient offset due to NCPA\n");
	    powfs[ipowfs].ncpa_grad=dcellnew(nwfs,1);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		double *realamp=powfs[ipowfs].realamp[iwfs];
		if(parms->powfs[ipowfs].gtype_sim==1){
		    pts_ztilt(&powfs[ipowfs].ncpa_grad->p[iwfs], powfs[ipowfs].pts,
			      powfs[ipowfs].nsaimcc>1?powfs[ipowfs].saimcc[iwfs]:powfs[ipowfs].saimcc[0], 
			      realamp, powfs[ipowfs].ncpa->p[iwfs]->p);
		}else{
		    spmulmat(&powfs[ipowfs].ncpa_grad->p[iwfs],adpind(powfs[ipowfs].GS0, iwfs),
			     powfs[ipowfs].ncpa->p[iwfs],1);
		}
	    }
	    if(parms->save.setup){
		dcellwrite(powfs[ipowfs].ncpa_grad, "%s/powfs%d_ncpa_grad",dirsetup,ipowfs);
	    }
	}
	if(parms->save.setup){
	    dcellwrite(powfs[ipowfs].ncpa, "%s/powfs%d_ncpa",dirsetup,ipowfs);

	}
    }/*if(fn_ncpa) */
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
	loc_t *loc_psf=mksqloc(ncompx,ncompy,dtheta, -ncompx2*dtheta, -ncompy2*dtheta);
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
			pn[iy][ix]=sinc(ir*duyp)*sinc(ia*duxp)
			    *exp(e0x*(ir*ir*duy2)+e0y*(ia*ia*dux2))
			    *pdtheta;
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
	if(parms->save.setup){
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
    double *p=powfs[ipowfs].focus->p;
    /*Convert range to focus. */
    double range2focus=pow(parms->aper.d/parms->powfs[ipowfs].hs,2)/(16.*sqrt(3));
    for(int ii=0; ii<powfs[ipowfs].focus->nx*powfs[ipowfs].focus->ny; ii++){
	p[ii]*=range2focus;
    }
}

/**
   Load and smooth out sodium profile. We preserve the sum of the sodium profile,
   which represents a scaling factor of the signal level.  */
static void setup_powfs_sodium(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    char *lltfn=parms->powfs[ipowfs].llt->fn;
    dmat *Nain=dread("%s",lltfn);
    if(Nain->ny<2){
	error("The sodium profile input %s is in wrong fromat\n", lltfn);
    }

    if(parms->powfs[0].llt->smooth){/*resampling the sodium profile by binning. */
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
		normalize(pout, nxnew, Nasum);
	    }
	    spfree(ht);
	    locfree(loc_in);
	    locfree(loc_out);
	    dfree(Nain);
	}else{
	    warning("Not smoothing sodium profile\n");
	    powfs[ipowfs].sodium=Nain;
	}
    }else{
	warning("Not smoothing sodium profile\n");
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
    /*find maximum elongation.  */
    /*Then test if # of pixels large enough. */
 
    assert(nwvl==1);/*sanity check. need to double check the code is nwvl!=1. */
    const double embfac=parms->powfs[ipowfs].embfac;
    /*setup elongation along radial direction. don't care azimuthal. */
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double wvl=parms->powfs[ipowfs].wvl[iwvl];
	const double dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
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
	      space for 6 LGS in NFIRAOS setup
	      warning("Using rot psf/otf method to do radial
	      pixel detector\n");
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
	cmat *etf=NULL;
	int netf;
	double dusc;
	double dtetf;
	int npad;
	int nover;
	if(use1d){
	    /*
	      2009-12-18:Thinking of padding here also to reduce aliasing. 
	      result: helps little, not much.
	    */
	    npad=2;/*zero padding to reduce aliasing? */
	    nover=2;
	    /*padding by x2. */
	    netf=ncompx*nover*npad;
	    dtetf=dtheta/nover;
	}else{
	    /*
	      We padd the array to incrase the sampling of
	      du. Don't shrink dtheta.  in use1d==0 case, du in
	      etf is smaller.  so interpolation will be more
	      accurate.
	    */
	    npad=2;
	    nover=2;/*doesn't change a lot. reduce aliasing */
	    /*padding by x2. */
	    netf=ncompx*nover*npad;
	    /*make du array bigger so that rotation is covered */
	    dtetf=dtheta/nover;
	}
	dusc=(netf*dtetf)/(dtheta*ncompx);
	etf=cnew(netf,1);
	double *thetas=calloc(netf, sizeof(double));
	int netf2=netf>>1;
	/*Only interpolating the center part. the rest is padding. */
	int etf0=netf2-(int)round(ncompx2*(dtheta/dtetf));
	int etf1=etf0+(int)round(ncompx*(dtheta/dtetf));
	if(etf0<0) error("Invalid configuration\n");
	for(int it=etf0; it<etf1; it++){
	    thetas[it]=(it-netf2)*dtetf;
	}
	cfft2plan(etf, -1);
	if(!powfs[ipowfs].sodium){
	    error("Sodium profile is NULL\n");
	}
	dmat *sodium=powfs[ipowfs].sodium;
       	int nhp=sodium->nx; 

	const double hs=parms->powfs[ipowfs].hs;
	double hpmin=sodium->p[0]/cos(parms->sim.za);
	double dhp1=1./(sodium->p[1]-sodium->p[0]);
	/*assume linear spacing. check the assumption valid */
	if(fabs(sodium->p[nhp-1]-sodium->p[0]-(nhp-1)/dhp1)>1.e-7){
	    error("llt profile is not evenly spaced:%g\n",
		  fabs(sodium->p[nhp-1]-sodium->p[0]-(nhp-1)/dhp1));
	}
	dhp1=dhp1*cos(parms->sim.za);
	if(sodium->ny-1<colskip){
	    error("Invalid configuration. colprep or colsim is too big\n");
	}
	colstart=colskip+istep%(sodium->ny-1-colskip);
	info2("Na using column %d in %s\n",colstart, mode==0?"preparation":"simulation");
	/*points to the effective sodium profile. */
	double* pp=sodium->p+nhp*(1+colstart);
	/*the sum of pp determines the scaling of the pixel intensity. */
	double i0scale=dblsum(pp, nhp);
	if(fabs(i0scale-1)>0.01){
	    warning("Siglev is scaled by %g by sodium profile\n", i0scale);
	}
	if(i0scale<0.8 || i0scale>2){
	    error("Check whether this is valid. Relax the restriction if desired.\n");
	}
	for(int illt=0; illt<nllt; illt++){
	    for(int isa=0; isa<nsa; isa++){
		/*1d ETF along radius. */
		double rsa=powfs[ipowfs].srsa->p[illt]->p[isa];
		double etf2sum=0;
		czero(etf);
		for(int icomp=etf0; icomp<etf1; icomp++){
		    /*peak in center */
		    const double itheta=thetas[icomp];
		    /*non linear mapping. */
		    const double ih=hs*rsa/(rsa-itheta*hs);
		    /*interpolating to get Na profile strenght. */
		    /*this is bilinear interpolation. not good.  */
		    /*need to do averaging */
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
		    warning("Wrong focus!\n");
		    if(use1d){
			petf[illt][isa]=cnew(ncompx,1);
		    }else{
			petf[illt][isa]=cnew(ncompx,ncompy);
		    }
		    cset(petf[illt][isa],1);
		}
		if(fuse_etf){
		    /*Fuse dtf into this 2d etf. */
		    ccwm(petf[illt][isa], pnominal[illt][isa*mnominal]);
		}
	    }
	}
	if(fuse_etf){
	    powfs[ipowfs].dtf[iwvl].fused=1;
	    if(parms->powfs[ipowfs].llt->colprep==parms->powfs[ipowfs].llt->colsim
	       && parms->powfs[ipowfs].llt->colsimdtrat==0){
		ccellfree(powfs[ipowfs].dtf[iwvl].nominal);
		info2("DTF nominal is fused to ETF and freed\n");
	    }else{
		info2("DTF nominal is fused to ETF but kept\n");
	    }
	}	    
	
	cfree(etf);
	free(thetas);
    }
}

/**
   setting up uplink pts/amp lotf
*/
static void 
setup_powfs_llt(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    LLT_T *llt=powfs[ipowfs].llt=calloc(1, sizeof(LLT_T));
    pts_t *lpts=llt->pts=calloc(1, sizeof(pts_t));
    lpts->nsa=1;
    double lltd=parms->powfs[ipowfs].llt->d;
    lpts->dsa=MAX(lltd, powfs[ipowfs].pts->dsa);
    int notf=MAX(powfs[ipowfs].ncompx, powfs[ipowfs].ncompy);
    /*The otf would be dx/lambda. Make it equal to pts->dsa/lambda/notf)*/
    const double dx=lpts->dx=parms->powfs[ipowfs].embfac*powfs[ipowfs].pts->dsa/notf;
    info("llt dx=%g\n", dx);
    const int nx=lpts->nx=round(lpts->dsa/lpts->dx);
    lpts->dsa=lpts->dx*lpts->nx;
    
    lpts->origx=calloc(1, sizeof(double));
    lpts->origy=calloc(1, sizeof(double));

    double oy=lpts->origx[0]=(dx-lpts->dsa)*0.5;
    double ox=lpts->origy[0]=(dx-lpts->dsa)*0.5;

    llt->amp=dnew(nx,nx);
    PDMAT(llt->amp, amps);
    double l2max =pow(lltd*0.5,2);
    double r2eff=l2max*pow(parms->powfs[ipowfs].llt->widthp,2);
    double sumamp2=0;
    for(int iy=0; iy<nx; iy++){
	double yy=iy*dx+oy;
	yy*=yy;
	for(int ix=0; ix<nx; ix++){
	    double xx=ix*dx+ox;
	    xx*=xx;
	    double r2=xx+yy;
	    if(r2<=l2max){
		amps[iy][ix]=exp(-r2/r2eff);
		sumamp2+=pow(amps[iy][ix],2);
	    }
	}
    }
    /*normalized so that max(otf)=1; */
    sumamp2=1./(sqrt(sumamp2));
    dscale(llt->amp, sumamp2);
    llt->loc=mksqloc(nx,nx,dx, lpts->origx[0], lpts->origy[0]);
    llt->mcc =pts_mcc_ptt(llt->pts, llt->amp->p);
    llt->imcc =dcellinvspd_each(llt->mcc);
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
    if(parms->save.setup){
	locwrite(llt->loc, "%s/powfs%d_llt_loc",dirsetup,ipowfs);
	dwrite(llt->amp, "%s/powfs%d_llt_amp", dirsetup,ipowfs);
	dcellwrite(llt->imcc, "%s/powfs%d_llt_imcc",dirsetup,ipowfs);
	if(llt->ncpa){
	    dcellwrite(llt->ncpa, "%s/powfs%d_llt_ncpa", dirsetup, ipowfs);
	}
	dcellwrite(powfs[ipowfs].srot, "%s/powfs%d_srot",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].srsa, "%s/powfs%d_srsa",dirsetup,ipowfs);

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
/**
   Setup the matched filter pixel processing parameters for physical optics wfs.
*/
static void 
setup_powfs_mtch(POWFS_T *powfs,const PARMS_T *parms, int ipowfs){
    if(parms->powfs[ipowfs].phytype!=1 ){/*not matched filter */
	warning("This is only intended for matched filter || \n");
    }
    long nsa=powfs[ipowfs].pts->nsa;
    if(powfs[ipowfs].intstat){
	dcellfree(powfs[ipowfs].intstat->mtche);
	dcellfree(powfs[ipowfs].intstat->mtchera);
	dcellfree(powfs[ipowfs].intstat->sanea);
	dcellfree(powfs[ipowfs].intstat->saneaxy);
	dcellfree(powfs[ipowfs].intstat->saneaixy);
	dfree(powfs[ipowfs].intstat->i0sum);
	free(powfs[ipowfs].intstat);
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
	}else if(parms->powfs[ipowfs].lo){
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
		if(powfs[ipowfs].ampm){
		    for(int iamp=0; iamp<powfs[ipowfs].ampm->nx; iamp++){
			key=dhash(powfs[ipowfs].ampm->p[iamp], key);
		    }
		}else{
		    key=dhash(powfs[ipowfs].amp, key);
		}
		if(parms->powfs[ipowfs].ncpa && parms->powfs[ipowfs].ncpa_method==2){
		    for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
			key=dhash(powfs[ipowfs].ncpa->p[iwfs],key);
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
			info2("Generating WFS OTF for %s...", fnotf);tic;
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
		cfft2plan(psfhat, 1);
		for(int illt=0; illt<intstat->lotf->ny; illt++){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			const double dx=powfs[ipowfs].llt->pts->dx;
			const double wvl=parms->powfs[ipowfs].wvl[iwvl];
			const double dpsf=wvl/(nlpsf*dx)*206265.;
			ccp(&psfhat, lotf[illt][iwvl]);
			cfftshift(psfhat);
			cifft2(psfhat, 1);
			cfftshift(psfhat);
			creal2d(&psf, 0, psfhat, 1);
			info2("illt %d, iwvl %d has FWHM of %g\"\n",
			      illt, iwvl, sqrt(4.*(double)dfwhm(psf)/M_PI)*dpsf);
		    }
		}
		cfree(psfhat);
		dfree(psf);
	    }
	    /*Generating short exposure images. */
	    gensepsf(parms,powfs,ipowfs);
	    if(parms->save.setup && intstat){
		dcellwrite(intstat->sepsf[0],
			   "%s/powfs%d_sepsf",dirsetup,ipowfs);
	    }
	    /*Free short exposure otf. */
	    ccellfree(intstat->lotf);
	    for(int iotf=0; iotf<intstat->notf; iotf++){
		ccellfree(intstat->otf[iotf]);
	    }
	    free(intstat->otf);
	}
	/*generate short exposure i0,gx,gy from psf. */
	gensei(parms,powfs,ipowfs);
	dcellfreearr(intstat->sepsf,intstat->nsepsf);
	if(parms->save.setup && intstat){
	    dcellwrite(intstat->i0,"%s/powfs%d_i0",dirsetup,ipowfs);
	    dcellwrite(intstat->gx,"%s/powfs%d_gx",dirsetup,ipowfs);
	    dcellwrite(intstat->gy,"%s/powfs%d_gy",dirsetup,ipowfs);
	}
    }
    /*Generating Matched filter */
    genmtch(parms,powfs,ipowfs);
 
    if(parms->powfs[ipowfs].neaphy){
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
    /*Remove OTFs that are older than 30 days. */
    char *dirotf=stradd(HOME, "/.aos/otfc", NULL);
    remove_file_older(dirotf, 30*24*3600);
    free(dirotf);
}

static void 
setup_powfs_cog(POWFS_T *powfs, const PARMS_T *parms, int ipowfs){
    if(parms->powfs[ipowfs].phytypesim!=2) return;
    const int nwfs=parms->powfs[ipowfs].nwfs;
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    powfs[ipowfs].gradphyoff=dcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(iwfs==0 || powfs[ipowfs].intstat->i0->ny>1){
	    int nsa=powfs[ipowfs].pts->nsa;
	    powfs[ipowfs].gradphyoff->p[iwfs]=dnew(nsa*2,1);
	    double *restrict gx=powfs[ipowfs].gradphyoff->p[iwfs]->p;
	    double *restrict gy=gx+nsa;
	    double g[2]={0,0};
	    for(int isa=0; isa<nsa; isa++){
		dmat *ints=powfs[ipowfs].intstat->i0->p[isa+iwfs*nsa];
		double maxi0=dmax(ints);
		dcog(g, ints, 0, 0, parms->powfs[ipowfs].cogthres*maxi0, parms->powfs[ipowfs].cogoff*maxi0);
		gx[isa]=g[0]*pixthetax;
		gy[isa]=g[1]*pixthetay;
	    }
	}else{
	    powfs[ipowfs].gradphyoff->p[iwfs]=dref(powfs[ipowfs].gradphyoff->p[0]);
	}
    }
}

/**
   Setup the powfs struct based on parms and aper. Everything about wfs are
   setup here.  \callgraph */
POWFS_T * setup_powfs(const PARMS_T *parms, APER_T *aper){
    POWFS_T *powfs=calloc(parms->npowfs, sizeof(POWFS_T));
    int ipowfs;
    tic;
    for(ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	info2("\n\033[0;32mSetting up powfs %d\033[0;0m\n\n", ipowfs);
	setup_powfs_geom(powfs,parms,aper,ipowfs);
        setup_powfs_grad(powfs,parms,ipowfs);
	setup_powfs_ncpa(powfs,parms,ipowfs);
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
	    if(parms->powfs[ipowfs].phytype==1){/*matched filter */
		setup_powfs_mtch(powfs,parms,ipowfs);
	    }else{
		warning("Please fill in this part. Use matched filter to estimate noise for the moment.\n");
		setup_powfs_mtch(powfs,parms,ipowfs);
	    }
	    if(parms->powfs[ipowfs].usephy && parms->powfs[ipowfs].phytypesim==2){
		setup_powfs_cog(powfs,parms,ipowfs);
	    }
	    dcellfree(powfs[ipowfs].intstat->i0);
	    dcellfree(powfs[ipowfs].intstat->gx);
	    dcellfree(powfs[ipowfs].intstat->gy);
	}
    }/*ipowfs */
    toc("setup_powfs");
    return powfs;
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
	    dcellfree(intstat->mtchera);
	    dcellfree(intstat->sanea);
	    dcellfree(intstat->saneara);
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
	    free(powfs[ipowfs].llt);
	    dcellfree(powfs[ipowfs].llt->ncpa);
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
    }
    free(powfs);
}
