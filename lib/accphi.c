/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
#include "thread.h"
#include "accphi.h"
/*
   2010-01-02: use a=(int)b instead of a=(int)floor(b).  call to floor is slow
   and unnecessary.  This only works if b is positive.  */
/**
   \file accphi.c Contains ray tracing routines optimized for different
   input/output formats. Notice that the OPDs are accumulated.
 */
#undef  EPS
#define EPS 1.e-12 /**<A threashold*/
/*
  The myfma() function computes x * y + z. without rounding, may be slower than x*y+z*/
/*
  Carefull: when design any algorithm, make sure you consider all possible cases
  and handle them properly. */

#define ONLY_FULL 0 /**<set to 1 to be compatible with other props.*/
#define USE_OPTIM 1 /**<Use the optimized version.*/
/**
   A wrapper prop routine that handles all the different cases by calling the
   different routines. Handles threading.
 */
void prop(thread_t *data){
    PROPDATA_T *propdata=data->data;
    const double displacex=propdata->displacex0+propdata->displacex1;
    const double displacey=propdata->displacey0+propdata->displacey1;
    switch(propdata->index){
    case 0:
	prop_grid_map(propdata->mapin, propdata->mapout,
		      propdata->alpha, displacex, displacey,
		      propdata->scale, propdata->wrap, data->start, data->end);
	break;
    case 1:
	prop_grid_pts(propdata->mapin, propdata->ptsout, 
		      propdata->ampout, propdata->phiout,
		      propdata->alpha, displacex, displacey,
		      propdata->scale, propdata->wrap, 
		      data->start, data->end);
	break;
    case 2:
	prop_grid(propdata->mapin, propdata->locout, 
		  propdata->ampout, propdata->phiout,
		  propdata->alpha, displacex, displacey,
		  propdata->scale, propdata->wrap, 
		  data->start, data->end);
	break;
    case 3:
	prop_grid_stat(propdata->mapin, propdata->ostat, propdata->phiout,
		       propdata->alpha, displacex, displacey,
		       propdata->scale, propdata->wrap, 
		       data->start, data->end);
	break;
    case 4:
	prop_nongrid_cubic(propdata->locin, propdata->phiin,
			   propdata->locout, propdata->ampout, propdata->phiout,
			   propdata->alpha, displacex, displacey, 
			   propdata->scale, propdata->cubic_iac, 
			   data->start, data->end);
	break;
    case 5:
	prop_nongrid_map_cubic(propdata->locin, propdata->phiin,
			       propdata->mapout, 
			       propdata->alpha, displacex, displacey, 
			       propdata->scale, propdata->cubic_iac, 
			       data->start, data->end);
	break;
    case 6:
	prop_nongrid_pts_cubic(propdata->locin, propdata->phiin,
			       propdata->ptsout, propdata->ampout, propdata->phiout,
			       propdata->alpha, displacex, displacey, 
			       propdata->scale, propdata->cubic_iac, 
			       data->start, data->end);
	break;
    case 7:
	prop_nongrid(propdata->locin, propdata->phiin,
		     propdata->locout, propdata->ampout, propdata->phiout,
		     propdata->alpha, displacex, displacey, 
		     propdata->scale,
		     data->start, data->end);
	break;
    case 8:
	prop_nongrid_map(propdata->locin, propdata->phiin,
			 propdata->mapout, 
			 propdata->alpha, displacex, displacey, 
			 propdata->scale, 
			 data->start, data->end);
	break;
    case 9:
	prop_nongrid_pts(propdata->locin, propdata->phiin,
			 propdata->ptsout, propdata->ampout, propdata->phiout,
			 propdata->alpha, displacex, displacey, 
			 propdata->scale, 
			 data->start, data->end);
	break;
    case 10:
	prop_grid_map_cubic(propdata->mapin, propdata->mapout,
			    propdata->alpha, displacex, displacey,
			    propdata->scale, propdata->cubic_iac,
			    data->start, data->end);
	break;
    case 11:
	prop_grid_pts_cubic(propdata->mapin, propdata->ptsout, 
			    propdata->ampout, propdata->phiout,
			    propdata->alpha, displacex, displacey,
			    propdata->scale, propdata->cubic_iac, 
			    data->start, data->end);
	break;
    case 12:
	prop_grid_cubic(propdata->mapin, propdata->locout, 
			propdata->ampout, propdata->phiout,
			propdata->alpha, displacex, displacey,
			propdata->scale, propdata->cubic_iac, 
			data->start, data->end);
	break;
    case 13:
	error("Invalid\n");
	break;
    default:
	error("Invalid\n");
    }
}
/**
   Identify the index of this ray tracing.
*/
void prop_index(PROPDATA_T *propdata){
    int done=0;
    if(propdata->mapin){
	if(propdata->cubic){
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		/*case 10 */
		propdata->index=10;
		done=1;
	    }
	    if(propdata->phiout){
		if(propdata->ptsout){
		    if(done) error("Invalid\n");
		    /*case 11 */
		    propdata->index=11;
		    done=1;
		}
		if(propdata->locout){
		    if(done) error("Invalid\n");
		    /*case 12 */
		    propdata->index=12;
		    done=1;
		}
		if(propdata->ostat){
		    error("Invalid\n");
		    /*case 13 */
		    propdata->index=13;
		    done=1;
		}
	    }
	}else{
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		/*case 0 */
		propdata->index=0;
		done=1;
	    }
	    if(propdata->phiout){
		if(propdata->ptsout){
		    if(done) error("Invalid\n");
		    /*case 1 */
		    propdata->index=1;
		    done=1;
		}
		if(propdata->locout){
		    if(done) error("Invalid\n");
		    /*case 2 */
		    propdata->index=2;
		    done=1;
		}
		if(propdata->ostat){
		    if(done) error("Invalid\n");
		    /*case 3 */
		    propdata->index=3;
		    done=1;
		}
	    }
	}
    }
    if(propdata->locin){
	if(propdata->cubic){
	    if(propdata->locout){
		if(done) error("Invalid\n");
		/*case 4 */
		propdata->index=4;
		done=1;
	    }
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		/*case 5 */
		propdata->index=5;
		done=1;
	    }
	    if(propdata->ptsout){
		if(done) error("Invalid\n");
		/*case 6 */
		propdata->index=6;
		done=1;
	    }
	}else{
	    if(propdata->locout){
		if(done) error("Invalid\n");
		/*case 7 */
		propdata->index=7;
		done=1;
	    }
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		/*case 8 */
		propdata->index=8;
		done=1;
	    }
	    if(propdata->ptsout){
		if(done) error("Invalid\n");
		/*case 9 */
		propdata->index=9;
		done=1;
	    }
	}
    }
    if(done==0) error("Invalid\n");
}

#define PREPIN_NONGRID(npad, nskip)			\
    /*padding to avoid test boundary*/			\
    loc_create_map_npad(locin, npad);			\
    const double dx_in1 = 1./locin->dx;			\
    const double dx_in2 = scale*dx_in1;			\
    displacex = (displacex-locin->map->ox)*dx_in1;	\
    displacey = (displacey-locin->map->oy)*dx_in1;	\
    long (*map)[locin->map->nx]=(void*)(locin->map->p);	\
    const int nxmax=locin->map->nx-nskip;		\
    const int nymax=locin->map->ny-nskip;		\
    /*-1 because we count from 1 in the map.*/		\
    const double *phiin0=phiin-1;			

#define PREPIN_GRID(nskip)				\
    const double dx_in1 = 1./mapin->dx;			\
    const double dx_in2 = scale*dx_in1;			\
    displacex = (displacex-mapin->ox)*dx_in1;		\
    displacey = (displacey-mapin->oy)*dx_in1;		\
    const int nxmax  = mapin->nx-nskip;			\
    const int nymax  = mapin->ny-nskip;			\
    double (*phiin)[mapin->nx]=(void*)(mapin->p);	

#define PREPOUT_LOC				\
    if(!locout) error("locout is NULL!");	\
    const double *px=locout->locx;		\
    const double *py=locout->locy;		\
    if(!end) end=locout->nloc;			

#define PREPOUT_PTS				\
    const double dxout=pts->dx;			\
    if(!end) end=pts->nsa;

#define PREPOUT_MAP				\
    const double dxout=mapout->dx;		\
    const double ox=mapout->ox;			\
    const double oy=mapout->oy;			\
    double *phiout=mapout->p;			\
    const int nxout=mapout->nx;			\
    if(!end) end=mapout->ny;

#if ONLY_FULL == 1
#define RUNTIME_LINEAR				\
    double dplocx, dplocy;			\
    int nplocx, nplocy, nplocx1, nplocy1;	\
    int missing=0;				\
    long iphi1,iphi2,iphi3,iphi4;			
#else
#define RUNTIME_LINEAR				\
    double dplocx, dplocy;			\
    int nplocx, nplocy, nplocx1, nplocy1;	\
    int missing=0;				\
    long iphi=0;					
#endif

#define MAKE_CUBIC_PARAM					\
    double fx[4],fy[4];						\
    const double cubicn=1./(1.+2.*cubic_iac);			\
    const double c0=1.*cubicn;					\
    const double c1=(4.*cubic_iac-2.5)*cubicn;			\
    const double c2=(1.5-3.*cubic_iac)*cubicn;			\
    const double c3=(2.*cubic_iac-0.5)*cubicn;			\
    const double c4=(0.5-cubic_iac)*cubicn; 

#define WARN_MISSING\
    ({if(missing>0) {warning("%d points not covered by input screen\n", missing); \
    print_backtrace(0);\
	 }})

#if ONLY_FULL == 1 
#define LINEAR_ADD_NONGRID						\
    iphi1=map[nplocy][nplocx];						\
    iphi2=map[nplocy][nplocx1];						\
    iphi3=map[nplocy1][nplocx];						\
    iphi4=map[nplocy1][nplocx1];					\
    if(iphi1 && iphi2 && iphi3 && iphi4){				\
	phiout[iloc]+=alpha*(phiin0[iphi1]*(1.-dplocx)*(1.-dplocy));	\
	phiout[iloc]+=alpha*(phiin0[iphi2]*(dplocx)*(1.-dplocy));	\
	phiout[iloc]+=alpha*(phiin0[iphi3]*(1.-dplocx)*(dplocy));	\
	phiout[iloc]+=alpha*(phiin0[iphi4]*(dplocx)*(dplocy));		\
    }
#else
#define LINEAR_ADD_NONGRID						\
    if((iphi=map[nplocy][nplocx]))					\
	phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(1.-dplocy));	\
    if((iphi=map[nplocy][nplocx1]))					\
	phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(1.-dplocy));	\
    if((iphi=map[nplocy1][nplocx]))					\
	phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(dplocy));	\
    if((iphi=map[nplocy1][nplocx1]))					\
	phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(dplocy));
#endif

#define RUNTIME_CUBIC					\
    register double dplocx, dplocy, dplocx0, dplocy0;	\
    int nplocx, nplocy;					\
    int missing=0;					\
    MAKE_CUBIC_PARAM;					

#define CUBIC_ADD_GRID						\
    fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);			\
    fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);			\
    fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);			\
    fx[3]=dplocx*dplocx*(c3+c4*dplocx);				\
								\
    fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);			\
    fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);			\
    fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);			\
    fy[3]=dplocy*dplocy*(c3+c4*dplocy);				\
    register double sum=0;					\
    for(int ky=-1; ky<3; ky++){					\
	for(int kx=-1; kx<3; kx++){				\
	    sum+=fx[kx+1]*fy[ky+1]*phiin[ky+nplocy][kx+nplocx];	\
	}							\
    }								\
    phiout[iloc]+=sum*alpha;

#define CUBIC_ADD_NONGRID						\
    fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);				\
    fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);				\
    fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);				\
    fx[3]=dplocx*dplocx*(c3+c4*dplocx);					\
									\
    fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);				\
    fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);				\
    fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);				\
    fy[3]=dplocy*dplocy*(c3+c4*dplocy);					\
									\
    for(int jy=-1; jy<3; jy++){						\
	for(int jx=-1; jx<3; jx++){					\
	    long iphi=map[jy+nplocy][jx+nplocx];			\
	    if(iphi){							\
		phiout[iloc]+=alpha*fx[jx+1]*fy[jy+1]*phiin0[iphi];	\
	    }								\
	}								\
    }


#include "prop_grid_pts.c"
#define TRANSPOSE 0
#include "prop_grid_stat.c"
#include "prop_grid_map.c"
#undef  TRANSPOSE

#define TRANSPOSE 1
#include "prop_grid_stat.c"
#include "prop_grid_map.c"
#undef  TRANSPOSE
/**
   Propagate OPD defines on grid mapin to coordinate locout.  alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect.*/
void prop_grid(ARGIN_GRID,
	       ARGOUT_LOC,
	       ARG_PROP,
	       int wrap,           /**<[in] wrap input OPD or not*/
	       long start,         /**<[in] First point to do*/
	       long end            /**<[in] Last point to do*/
	       ){
    PREPIN_GRID(1);
    PREPOUT_LOC;
    RUNTIME_LINEAR;
    (void)iphi;
    const int nx = mapin->nx;
    const int ny = mapin->ny;

    for(long iloc=start; iloc<end; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;/*skip points that has zero amplitude */
	dplocx=myfma(px[iloc],dx_in2,displacex);
	dplocy=myfma(py[iloc],dx_in2,displacey);
	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>=nxmax||nplocy<0||nplocy>=nymax){
	    if(wrap){
		while(nplocx<0)
		    nplocx+=nx;
		while(nplocx>nxmax)
		    nplocx-=nx;
		while(nplocy<0)
		    nplocy+=ny;
		while(nplocy>nymax)
		    nplocy-=ny;
		
		nplocx1=(nplocx==nxmax?0:nplocx+1);
		nplocy1=(nplocy==nymax?0:nplocy+1);
	    }else{
		missing++;
		continue;
	    }
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}
	phiout[iloc]+=alpha*
	    (+(phiin[nplocy][nplocx]*(1.-dplocx)
	       +phiin[nplocy][nplocx1]*dplocx)*(1.-dplocy)
	     +(phiin[nplocy1][nplocx]*(1.-dplocx)
	       +phiin[nplocy1][nplocx1]*dplocx)*dplocy);
    }
    WARN_MISSING;
}

/**
   Propagate OPD defines on coordinate locin to coordinate locout.  This is the
   <em>slowest</em>. alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. See prop_grid() for definition of other parameters.*/
void prop_nongrid(ARGIN_NONGRID,
		  ARGOUT_LOC,
		  ARG_PROP,
		  long start,          /**<[in] First point to do*/
		  long end             /**<[in] Last point to do*/
		  ){
    PREPIN_NONGRID(1,2);
    PREPOUT_LOC;
    RUNTIME_LINEAR;
    for(long iloc=start; iloc<end; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;/*skip points that has zero amplitude */
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>nxmax||nplocy<0||nplocy>nymax){
	    missing++;
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}
	LINEAR_ADD_NONGRID;
    }
    WARN_MISSING;
}
/**
   Propagate OPD defines on coordinate locin to grid mapout. alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect. */
void prop_nongrid_map(ARGIN_NONGRID,
		      ARGOUT_MAP,
		      ARG_PROP,
		      long start,       /**<[in] First point to do*/
		      long end          /**<[in] Last point to do*/
		      ){
    PREPIN_NONGRID(1,2);/*Do we need 2 in first? */
    PREPOUT_MAP;
    RUNTIME_LINEAR ;
    for(int iy=start; iy<end; iy++){
	dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<0||nplocy>nymax){
	    missing++;
	    continue;
	}else{
	    nplocy1=nplocy+1;
	}
	for(int ix=0; ix<nxout; ix++){
	    int iloc=ix+iy*nxout;
	    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
	    SPLIT(dplocx,dplocx,nplocx);
	    if(nplocx<0||nplocx>nxmax){
		missing++;
		continue;
	    }
	    nplocx1=nplocx+1;
	    LINEAR_ADD_NONGRID;
	}
    }
    WARN_MISSING;
}
/**
   Propagate OPD defines on coordinate locin to subapertures pts. alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect. See prop_groid().*/
void prop_nongrid_pts(ARGIN_NONGRID,
		      ARGOUT_PTS,
		      ARG_PROP,
		      long start,           /**<[in] First point to do*/
		      long end              /**<[in] Last point to do*/
		      ){
    PREPIN_NONGRID(1,2);
    PREPOUT_PTS;
    RUNTIME_LINEAR;
 
    for(int isa=start; isa<end; isa++){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(nplocy<0||nplocy>nymax){
		continue;
	    }
	    nplocy1=nplocy+1;
	    for(int ix=0; ix<pts->nx; ix++){
		iloc++;
		if(ampout && fabs(ampout[iloc])<EPS)
		    continue;/*skip points that has zero amplitude */
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		SPLIT(dplocx,dplocx,nplocx);
		if(nplocx<0||nplocx>nxmax){
		    missing++;
		    continue;
		}
		nplocx1=nplocx+1;
		LINEAR_ADD_NONGRID;
	    }    
	}
    }
    WARN_MISSING;
}


/**
   Propagate OPD defines on grid mapin to coordinate locout with cubic influence
   functions.  alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. The input grid must cover the output loc completely (we aim for best
   efficiency).*/
void prop_grid_cubic(ARGIN_GRID,
		     ARGOUT_LOC,
		     ARG_PROP,
		     double cubic_iac,   /**<[in] inter-actuator-coupling for cubic*/
		     long start,         /**<[in] First point to do*/
		     long end            /**<[in] Last point to do*/
		     ){
    PREPIN_GRID(3);
    PREPOUT_LOC;
    RUNTIME_CUBIC;

    for(long iloc=start; iloc<end; iloc++){
	dplocx=myfma(px[iloc],dx_in2,displacex);
	dplocy=myfma(py[iloc],dx_in2,displacey);
	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<1||nplocx>nxmax||nplocy<1||nplocy>nymax){
	    missing++;
	    continue;
	}
	dplocy0=1.-dplocy;
	dplocx0=1.-dplocx;
	CUBIC_ADD_GRID;
    }
    WARN_MISSING;
}
/**
   Propagate OPD defines on grid mapin to coordinate locout with cubic influence
   functions.  alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. The input grid must cover the output loc completely (we aim for best
   efficiency).*/
void prop_grid_pts_cubic(ARGIN_GRID,
			 ARGOUT_PTS,
			 ARG_PROP,
			 double cubic_iac,   /**<[in] inter-actuator-coupling for cubic*/
			 long start,         /**<[in] First point to do*/
			 long end            /**<[in] Last point to do*/
			 ){
    PREPIN_GRID(3);
    PREPOUT_PTS;
    RUNTIME_CUBIC;
    for(int isa=start; isa<end; isa++){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(nplocy<1||nplocy>nymax){
		continue;
	    }
	    dplocy0=1.-dplocy;
	    for(int ix=0; ix<pts->nx; ix++){
		iloc++;
		if(ampout && fabs(ampout[iloc])<EPS)
		    continue;/*skip points that has zero amplitude */
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		SPLIT(dplocx,dplocx,nplocx);
		if(nplocx<1||nplocx>nxmax){
		    continue;
		}
		dplocx0=1.-dplocx;
		CUBIC_ADD_GRID;
	    }
	}
    }
    WARN_MISSING;
}
/**
   like prop_grid_map() but with cubic influence functions. cubic_iac is the
inter-actuator-coupling. */
void prop_grid_map_cubic(ARGIN_GRID,
			 ARGOUT_MAP,
			 ARG_PROP,
			 double cubic_iac,
			 long start, long end){
    PREPIN_GRID(3);
    PREPOUT_MAP;
    RUNTIME_CUBIC;
    for(int iy=start; iy<end; iy++){
	dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<1||nplocy>nymax){
	    continue;
	}
	dplocy0=1.-dplocy;
	for(int ix=0; ix<nxout; ix++){
	    int iloc=ix+iy*nxout;
	    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
	    SPLIT(dplocx,dplocx,nplocx);
	    if(nplocx<1||nplocx>nxmax){
		continue;
	    }
	    dplocx0=1.-dplocx;
	    CUBIC_ADD_GRID;
	}
    }
    WARN_MISSING;
}
/**
   like prop_nongrid() but with cubic influence functions. cubic_iac is the
   inter-actuator coupling. Consider embed the input into a map_t and call
   prop_grid_cubic instead. That one is must faster.*/
void prop_nongrid_cubic(ARGIN_NONGRID,
			ARGOUT_LOC,
			ARG_PROP,
			double cubic_iac,
			long start, long end){
    PREPIN_NONGRID(2,3);
    PREPOUT_LOC;
    RUNTIME_CUBIC;

    for(long iloc=start; iloc<end; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;/*skip points that has zero amplitude */
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<1 || nplocy>nymax || nplocx<1 || nplocx>nxmax){
	    missing++;
	    continue;
	}
	dplocy0=1.-dplocy;
	dplocx0=1.-dplocx;
	CUBIC_ADD_NONGRID;
    }
    WARN_MISSING;
}
/**
   like prop_nongrid_pts() but with cubic influence functions. cubic_iac is the
inter-actuator-coupling.  */
void prop_nongrid_pts_cubic(ARGIN_NONGRID,
			    ARGOUT_PTS,
			    ARG_PROP,
			    double cubic_iac, 
			    long start, long end){
    PREPIN_NONGRID(2,3);
    PREPOUT_PTS;
    RUNTIME_CUBIC;
    for(int isa=start; isa<end; isa++){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(nplocy<1||nplocy>nymax){
		continue;
	    }
	    dplocy0=1.-dplocy;
	    for(int ix=0; ix<pts->nx; ix++){
		iloc++;
		if(ampout && fabs(ampout[iloc])<EPS)
		    continue;/*skip points that has zero amplitude */
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		SPLIT(dplocx,dplocx,nplocx);
		if(nplocx<1||nplocx>nxmax){
		    continue;
		}
		dplocx0=1.-dplocx;
		CUBIC_ADD_NONGRID;
	    }
	}
    }
    WARN_MISSING;
}
/**
   like prop_nongrid_map() but with cubic influence functions. cubic_iac is the
inter-actuator-coupling. */
void prop_nongrid_map_cubic(ARGIN_NONGRID,
			    ARGOUT_MAP,
			    ARG_PROP,
			    double cubic_iac,
			    long start, long end){
    PREPIN_NONGRID(2,3);
    PREPOUT_MAP;
    RUNTIME_CUBIC;
    for(int iy=start; iy<end; iy++){
	dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<1||nplocy>nymax){
	    continue;
	}
	dplocy0=1.-dplocy;
	for(int ix=0; ix<nxout; ix++){
	    int iloc=ix+iy*nxout;
	    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
	    SPLIT(dplocx,dplocx,nplocx);
	    if(nplocx<1||nplocx>nxmax){
		continue;
	    }
	    dplocx0=1.-dplocx;
	    CUBIC_ADD_NONGRID;
	}
    }
    WARN_MISSING;
}

/**
  the following routine is used to do down sampling by doing binning ray tracing
  using reverse interplation.  locout is coarse sampling, locin is fine
  sampling. phiout is the destination OPD. The weightings are obtained by
  interpolating from locout to locin, but the OPD are reversed computed. Simply
  replace prop_nongrid by prop_nongrid_reverse without changing arguments,
  except removing start, end, will do the same ray tracing using reverse
  interpolation (binning). ampout is not used.

  2011-04-27: Revised usage of alpha, displacex/y so that the result agrees with
  prop_nongrid when used in the same situation. Report missing does not make
  sense here since locin is usually bigger than locout.
  */
void prop_nongrid_bin(const loc_t *locin,
		      const double* phiin,
		      loc_t *locout, 
		      const double *ampout,
		      double* phiout,
		      ARG_PROP){
    if(locout->dx<locin->dx) {
	error("This routine is designed for down sampling.\n");
    }
    (void) ampout;
    loc_create_map_npad(locout,1);/*will only do once and save in locout. */
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    const int wrapx1 = locout->map->nx;
    const int wrapy1 = locout->map->ny;
    const int nxmax = wrapx1-1;
    const int nymax = wrapy1-1;
    const double dx_out1 = 1./locout->dx;
    /*notice inverse of scale. */
    const double dx_out2 = (1./scale)*dx_out1;
    /*notice negative sign in displacex/y. */
    displacex = (-displacex/scale-locout->map->ox)*dx_out1;
    displacey = (-displacey/scale-locout->map->oy)*dx_out1;
    const double *px=locin->locx;
    const double *py=locin->locy;
    /*Scale alpha to cancel out scaling */
    double alpha2 = alpha*pow(locin->dx/locout->dx/scale,2);
#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    long (*map)[locout->map->nx]
	=(long(*)[locout->map->nx])(locout->map->p);
    /*-1 because we count from 1 in the map. */
    double *phiout0=phiout-1;
    for(long iloc=0; iloc<locin->nloc; iloc++){
	dplocy=myfma(py[iloc],dx_out2,displacey);
	dplocx=myfma(px[iloc],dx_out2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>nxmax||nplocy<0||nplocy>nymax){
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}

#if ONLY_FULL == 1 /*only proceed if all four points exist. */
	iphi1=map[nplocy][nplocx];
	iphi2=map[nplocy][nplocx1];
	iphi3=map[nplocy1][nplocx];
	iphi4=map[nplocy1][nplocx1];
	if(iphi1 && iphi2 && iphi3 && iphi4){
	    phiout0[iphi1]+=alpha2*(phiin[iloc]*(1.-dplocx)*(1.-dplocy));
	    phiout0[iphi2]+=alpha2*(phiin[iloc]*(dplocx)*(1.-dplocy));
	    phiout0[iphi3]+=alpha2*(phiin[iloc]*(1.-dplocx)*(dplocy));
	    phiout0[iphi4]+=alpha2*(phiin[iloc]*(dplocx)*(dplocy));
	}
#else	
	if((iphi=map[nplocy][nplocx])) 
	    phiout0[iphi]+=alpha2*(phiin[iloc]*(1.-dplocx)*(1.-dplocy));
	if((iphi=map[nplocy][nplocx1]))
	    phiout0[iphi]+=alpha2*(phiin[iloc]*(dplocx)*(1.-dplocy));
	if((iphi=map[nplocy1][nplocx]))
	    phiout0[iphi]+=alpha2*(phiin[iloc]*(1.-dplocx)*(dplocy));
	if((iphi=map[nplocy1][nplocx1]))
	    phiout0[iphi]+=alpha2*(phiin[iloc]*(dplocx)*(dplocy));
#endif
    }
}
