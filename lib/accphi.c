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
#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-variable"
#elif defined __GNUC__
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include "accphi.h"
#undef  EPS
#define EPS 1.e-12 /**<A threashold*/

/**
   Identify the index of this ray tracing.
*/
static void prop_index(PROPDATA_T *propdata){
    int done=0;
    if(propdata->mapin){
	if(propdata->mapin->iac){
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
		    if(done) error("Invalid\n");
		    /*case 13 */
		    propdata->index=14;
		    done=1;
		}
	    }
	}else{
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		/*case 0 */
		propdata->index=13;
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
	if(propdata->locin->iac){
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

/**
   A wrapper prop routine that handles all the different cases by calling the
   different routines. Handles threading.
*/
void prop(thread_t *data){
    PROPDATA_T *propdata=data->data;
    if(!propdata->index){
	prop_index(propdata);
    }
    const double displacex=propdata->displacex0+propdata->displacex1;
    const double displacey=propdata->displacey0+propdata->displacey1;
    switch(propdata->index){
    case 1:
	prop_grid_pts(propdata->mapin, propdata->ptsout, propdata->phiout,
		      propdata->alpha, displacex, displacey,
		      propdata->scale, propdata->wrap, 
		      data->start, data->end);
	break;
    case 2:
	prop_grid(propdata->mapin, propdata->locout, propdata->phiout,
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
			   propdata->locout, propdata->phiout,
			   propdata->alpha, displacex, displacey, 
			   propdata->scale, propdata->locin->iac, 
			   data->start, data->end);
	break;
    case 5:
	prop_nongrid_map_cubic(propdata->locin, propdata->phiin,
			       propdata->mapout, 
			       propdata->alpha, displacex, displacey, 
			       propdata->scale, propdata->locin->iac, 
			       data->start, data->end);
	break;
    case 6:
	prop_nongrid_pts_cubic(propdata->locin, propdata->phiin,
			       propdata->ptsout, propdata->phiout,
			       propdata->alpha, displacex, displacey, 
			       propdata->scale, propdata->locin->iac, 
			       data->start, data->end);
	break;
    case 7:
	prop_nongrid(propdata->locin, propdata->phiin,
		     propdata->locout, propdata->phiout,
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
			 propdata->ptsout, propdata->phiout,
			 propdata->alpha, displacex, displacey, 
			 propdata->scale, 
			 data->start, data->end);
	break;
    case 10:
	prop_grid_map_cubic(propdata->mapin, propdata->mapout,
			    propdata->alpha, displacex, displacey,
			    propdata->scale, propdata->mapin->iac,
			    data->start, data->end);
	break;
    case 11:
	prop_grid_pts_cubic(propdata->mapin, propdata->ptsout, propdata->phiout,
			    propdata->alpha, displacex, displacey,
			    propdata->scale, propdata->mapin->iac, 
			    data->start, data->end);
	break;
    case 12:
	prop_grid_cubic(propdata->mapin, propdata->locout, propdata->phiout,
			propdata->alpha, displacex, displacey,
			propdata->scale, propdata->mapin->iac, 
			data->start, data->end);
	break;
    case 13:
	prop_grid_map(propdata->mapin, propdata->mapout,
		      propdata->alpha, displacex, displacey,
		      propdata->scale, propdata->wrap, data->start, data->end);
	break;
    case 14:
	prop_grid_stat_cubic(propdata->mapin, propdata->ostat, propdata->phiout,
			    propdata->alpha, displacex, displacey,
			    propdata->scale, propdata->mapin->iac,
			    data->start, data->end);
	break;
    default:
	error("Invalid\n");
    }
}


#define PREPIN_NONGRID(nskip)			\
    /*padding to avoid test boundary*/			\
    if(!locin->map) loc_create_map_npad(locin, nskip,0,0);	\
    const double dx_in1 = 1./locin->dx;			\
    const double dx_in2 = scale*dx_in1;			\
    const double dy_in1 = 1./locin->dy;			\
    const double dy_in2 = scale*dy_in1;			\
    displacex = (displacex-locin->map->ox)*dx_in1;	\
    displacey = (displacey-locin->map->oy)*dy_in1;	\
    map_t*  map=locin->map;				\
    const int nxmin=MAX(nskip,locin->npad);		\
    const int nymin=nxmin;				\
    const int nxmax=locin->map->nx-nxmin-1;		\
    const int nymax=locin->map->ny-nymin-1;		\
    /*-1 because we count from 1 in the map.*/		\
    const double *phiin0=phiin-1;			

#define PREPIN_GRID(nskip)				\
    const double dx_in1 = 1./mapin->dx;			\
    const double dx_in2 = scale*dx_in1;			\
    const double dy_in1 = 1./mapin->dy;			\
    const double dy_in2 = scale*dy_in1;			\
    displacex = (displacex-mapin->ox)*dx_in1;		\
    displacey = (displacey-mapin->oy)*dy_in1;		\
    const int nxmin=nskip;				\
    const int nymin=nskip;				\
    const int nxmax  = mapin->nx-nskip-1;		\
    const int nymax  = mapin->ny-nskip-1;		\

#define PREPOUT_LOC				\
    if(!locout) error("locout is NULL!");	\
    const double *px=locout->locx;		\
    const double *py=locout->locy;		\
    if(!end) end=locout->nloc;			\

#define PREPOUT_PTS				\
    const double dxout=pts->dx;			\
    const double dyout=pts->dy;			\
    if(!end) end=pts->nsa;

#define PREPOUT_MAP				\
    const double dxout=mapout->dx;		\
    const double dyout=mapout->dy;		\
    const double ox=mapout->ox;			\
    const double oy=mapout->oy;			\
    double *phiout=mapout->p;			\
    const int nxout=mapout->nx;			\
    if(!end) end=mapout->ny;

#define PREPOUT_STAT				\
    const double dxout  = ostat->dx;		\
    if(colend==0) colend = ostat->ncol;

#define RUNTIME_LINEAR				\
    double dplocx, dplocy;			\
    int nplocx, nplocy, nplocx1, nplocy1;	\
    int missing=0;				

#define MAKE_CUBIC_PARAM			\
    const double cubicn=1./(1.+2.*cubic_iac);	\
    const double c0=1.*cubicn;			\
    const double c1=(4.*cubic_iac-2.5)*cubicn;	\
    const double c2=(1.5-3.*cubic_iac)*cubicn;	\
    const double c3=(2.*cubic_iac-0.5)*cubicn;	\
    const double c4=(0.5-cubic_iac)*cubicn; 

#define WARN_MISSING							\
    ({static int printed=0; if(missing>0 && !printed) {printed=1; warning("%d points not covered by input screen\n", missing); \
	    print_backtrace(); }})

#define LINEAR_ADD_NONGRID						\
    long iphi; double tmp=0; double wt=0; double wtsum=0;		\
    wt=(1.-dplocx)*(1.-dplocy);						\
    if(wt>EPS){/*this test fixed to top/right boundary defect*/		\
	if((iphi=fabs(IND(map,nplocx,nplocy)))){				\
	    tmp+=(phiin0[iphi]*wt);					\
	    wtsum+=wt;							\
	}								\
    }									\
    wt=(dplocx)*(1.-dplocy);						\
    if(wt>EPS){								\
	if((iphi=fabs(IND(map,nplocx1,nplocy)))){				\
	    tmp+=(phiin0[iphi]*wt);					\
	    wtsum+=wt;							\
	}								\
    }									\
    wt=(1.-dplocx)*(dplocy);						\
    if(wt>EPS){								\
	if((iphi=fabs(IND(map,nplocx,nplocy1)))){				\
	    tmp+=(phiin0[iphi]*wt);					\
	    wtsum+=wt;							\
	}								\
    }									\
    wt=(dplocx)*(dplocy);						\
    if(wt>EPS){								\
	if((iphi=fabs(IND(map,nplocx1,nplocy1)))){				\
	    tmp+=(phiin0[iphi]*wt);					\
	    wtsum+=wt;							\
	}								\
    }									\
    /*Wt use sum(opd*wt)/sum(wt) to avoid edge roll off*/		\
    if(wtsum>EPS) phiout[iloc]+=alpha*tmp/wtsum; 


#define RUNTIME_CUBIC					\
    register double dplocx, dplocy, dplocx0, dplocy0;	\
    int nplocx, nplocy;					\
    int missing=0;					\
    MAKE_CUBIC_PARAM;					

#define MAKE_CUBIC_COEFF			\
    double fx[4],fy[4];				\
    fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);	\
    fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);	\
    fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);	\
    fx[3]=dplocx*dplocx*(c3+c4*dplocx);		\
						\
    fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);	\
    fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);	\
    fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);	\
    fy[3]=dplocy*dplocy*(c3+c4*dplocy);			

#define CUBIC_ADD_GRID					\
    register double sum=0, sumwt=0;			\
    for(int ky=-1; ky<3; ky++){				\
	for(int kx=-1; kx<3; kx++){			\
	    double wt=fx[kx+1]*fy[ky+1];		\
	    if(wt>EPS){					\
		double tmp=IND(mapin,kx+nplocx,ky+nplocy);	\
		sum+=wt*tmp;sumwt+=wt;			\
	    }						\
	}						\
    }							\
    if(sumwt>EPS && sum==sum) phiout[iloc]+=alpha*sum/sumwt;

#define CUBIC_ADD_NONGRID					\
    register double sum=0,sumwt=0;				\
    for(int jy=-1; jy<3; jy++){					\
	for(int jx=-1; jx<3; jx++){				\
	    long iphi;						\
	    double wt=fx[jx+1]*fy[jy+1];			\
	    if(wt>EPS){						\
		if((iphi=fabs(IND(map,jx+nplocx,jy+nplocy)))){	\
		    sum+=wt*phiin0[iphi];			\
		    sumwt+=wt;					\
		}						\
	    }							\
	}							\
    }								\
    if(sumwt>EPS) phiout[iloc]+=alpha*sum/sumwt;

/**
   A building block for ray tracing onto a rectangular block of points
*/

DEF_ENV_FLAG(PROP_GRID_MAP_OPTIM, 1);

#define TRANSPOSE 0
#include "prop_grid.c"
#undef  TRANSPOSE

#define TRANSPOSE 1
#include "prop_grid.c"
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
    if(mapin->iac){
	if(wrap) error("wrap=1 is invalid for cubic\n");
	prop_grid_cubic(ARGIN_GRID2, ARGOUT_LOC2, ARG_PROP2, mapin->iac, start, end);
	return;
    }
    PREPIN_GRID(0);
    PREPOUT_LOC;
    RUNTIME_LINEAR;
    (void)nxmin; (void)nymin; 
    const int nx = mapin->nx;
    const int ny = mapin->ny;
    OMPTASK_FOR(iloc, start, end, private(nplocx,nplocy,nplocx1,nplocy1,dplocx,dplocy)){
	// The myfma() function computes x * y + z. without rounding, may be slower than x*y+z
	dplocx=myfma(px[iloc],dx_in2,displacex);
	dplocy=myfma(py[iloc],dy_in2,displacey);
	if(dplocx<0||dplocx>nxmax||dplocy<0||dplocy>nymax){
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    //Handle top/right boundary correctly
	    if(wrap){
		while(nplocx<0)
		    nplocx+=nx;
		while(nplocx>nxmax)
		    nplocx-=nx;
		while(nplocy<0)
		    nplocy+=ny;
		while(nplocy>nymax)
		    nplocy-=ny;
		
	    }else{
		missing++;
		continue;
	    }
	}else{
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	}
	nplocx1=(nplocx==nxmax?0:nplocx+1);
	nplocy1=(nplocy==nymax?0:nplocy+1);
	double tmp=(+(IND(mapin,nplocx,nplocy)*(1.-dplocx)
		      +IND(mapin,nplocx1,nplocy)*dplocx)*(1.-dplocy)
		    +(IND(mapin,nplocx,nplocy1)*(1.-dplocx)
		      +IND(mapin,nplocx1,nplocy1)*dplocx)*dplocy);
	add_valid(phiout[iloc],alpha, tmp);
    }
    OMPTASK_END;
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
    if(locin->iac){
	prop_nongrid_cubic(ARGIN_NONGRID2, ARGOUT_LOC2, ARG_PROP2, locin->iac, start, end);
	return;
    }
    PREPIN_NONGRID(0);
    PREPOUT_LOC;
    RUNTIME_LINEAR;
    OMPTASK_FOR(iloc, start, end, private(nplocx,nplocy,nplocx1,nplocy1,dplocx,dplocy)){
	dplocx=myfma(px[iloc],dx_in2,displacex);
	dplocy=myfma(py[iloc],dy_in2,displacey);
	if(dplocy>=nymin && dplocy<=nymax && dplocx>=nxmin && dplocx<=nxmax){
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	    LINEAR_ADD_NONGRID;
	}else{
	    missing++;
	    continue;
	}
    }
    OMPTASK_END;
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
    if(locin->iac){
	prop_nongrid_map_cubic(ARGIN_NONGRID2, ARGOUT_MAP2, ARG_PROP2, locin->iac, start, end);
	return;
    }
    PREPIN_NONGRID(0);
    PREPOUT_MAP;
    RUNTIME_LINEAR ;
    OMPTASK_FOR(iy, start, end, private(nplocy, dplocy, nplocy1, nplocx,dplocx,nplocx1)){
	dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	if(dplocy>=nymin && dplocy<=nymax){
	    SPLIT(dplocy,dplocy,nplocy);
	    nplocy1=nplocy+1;
	    for(int ix=0; ix<nxout; ix++){
		int iloc=ix+iy*nxout;
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		if(dplocx>=nxmin && dplocx<=nxmax){
		    SPLIT(dplocx,dplocx,nplocx);
		    nplocx1=nplocx+1;
		    LINEAR_ADD_NONGRID;
		}else{
		    missing++;
		}
	    }
	}else{
	    missing++;
	}
    }
    OMPTASK_END;
    //WARN_MISSING;//will always be triggered.
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
    if(locin->iac){
	prop_nongrid_pts_cubic(ARGIN_NONGRID2, ARGOUT_PTS2, ARG_PROP2, locin->iac, start, end);
	return;
    }
    PREPIN_NONGRID(0);
    PREPOUT_PTS;
    RUNTIME_LINEAR;
    
    OMPTASK_FOR(isa, start, end, private(nplocx,nplocy,dplocx,dplocy,nplocx1,nplocy1)){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	    if(dplocy>=nymin && dplocy<=nymax){
		SPLIT(dplocy,dplocy,nplocy);
		nplocy1=nplocy+1;
		for(int ix=0; ix<pts->nx; ix++){
		    iloc++;
		    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		    if(dplocx>=nxmin && dplocx<=nxmax){
			SPLIT(dplocx,dplocx,nplocx);
			nplocx1=nplocx+1;
			LINEAR_ADD_NONGRID;
		    }else{
			missing++;
		    }
		} 
	    }else{
		missing++;
	    }
	}
    }
    OMPTASK_END;
    WARN_MISSING;
}


/**
   Propagate OPD defines on grid mapin to coordinate locout with cubic influence
   functions.  alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. The input grid must cover the output loc completely (we aim for best
   efficiency), with at least one full guard ring.
*/
void prop_grid_cubic(ARGIN_GRID,
		     ARGOUT_LOC,
		     ARG_PROP,
		     double cubic_iac,   /**<[in] inter-actuator-coupling for cubic*/
		     long start,         /**<[in] First point to do*/
		     long end            /**<[in] Last point to do*/
    ){
    PREPIN_GRID(1);
    (void)nymin;
    PREPOUT_LOC;
    RUNTIME_CUBIC;

    OMPTASK_FOR(iloc, start, end, private(dplocx,dplocy,nplocx,nplocy,dplocx0,dplocy0)){
	dplocx=myfma(px[iloc],dx_in2,displacex);
	dplocy=myfma(py[iloc],dy_in2,displacey);
	if(dplocx>=nxmin&&dplocx<=nxmax&&dplocy>=nxmin&&dplocy<=nymax){
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    dplocy0=1.-dplocy;
	    dplocx0=1.-dplocx;
	    MAKE_CUBIC_COEFF;
	    CUBIC_ADD_GRID;
	}else{
	    missing++;
	}
    }
   
    OMPTASK_END;
    WARN_MISSING;
}
/**
   Propagate OPD defines on grid mapin to coordinate locout with cubic influence
   functions.  alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. The input grid must cover the output loc completely (we aim for best
   efficiency), with at least one full guard ring.
*/
void prop_grid_pts_cubic(ARGIN_GRID,
			 ARGOUT_PTS,
			 ARG_PROP,
			 double cubic_iac,   /**<[in] inter-actuator-coupling for cubic*/
			 long start,         /**<[in] First point to do*/
			 long end            /**<[in] Last point to do*/
    ){
    PREPIN_GRID(1);
    PREPOUT_PTS;
    RUNTIME_CUBIC;
    OMPTASK_FOR(isa, start, end, private(nplocx,nplocy,dplocx0,dplocy0,dplocx,dplocy)){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];

	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	    if(dplocy>=nymin && dplocy<=nymax){
		SPLIT(dplocy,dplocy,nplocy);
		dplocy0=1.-dplocy;
		for(int ix=0; ix<pts->nx; ix++){
		    iloc++;
		    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		    if(dplocx>=nxmin && dplocx<=nxmax){
			SPLIT(dplocx,dplocx,nplocx);
			dplocx0=1.-dplocx;
			MAKE_CUBIC_COEFF;
			CUBIC_ADD_GRID;
		    }else{
			missing++;
		    }
		}
	    }else{
		missing++;
	    }
	}
    }
    OMPTASK_END;
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
    PREPIN_GRID(1);
    PREPOUT_MAP;
    RUNTIME_CUBIC;
    OMPTASK_FOR(iy, start, end, private(dplocx,nplocx,dplocx0,nplocy,dplocy0,dplocy)){
	dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	if(dplocy>=nymin && dplocy<=nymax){
	    SPLIT(dplocy,dplocy,nplocy);
	    dplocy0=1.-dplocy;
	    for(int ix=0; ix<nxout; ix++){
		int iloc=ix+iy*nxout;//output index
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		if(dplocx>=nxmin && dplocx<=nxmax){
		    SPLIT(dplocx,dplocx,nplocx);
		    dplocx0=1.-dplocx;
		    MAKE_CUBIC_COEFF;
		    CUBIC_ADD_GRID;
		}else{
		    missing++;
		}
	    }
	}else{
	    missing++;
	}
    }
    OMPTASK_END;
    //WARN_MISSING;
}
/**
   like prop_grid_stat() but with cubic influence function.
 */
void prop_grid_stat_cubic(ARGIN_GRID,
			  ARGOUT_STAT,
			  ARG_PROP,
			  double cubic_iac,
			  long colstart, long colend){
    PREPIN_GRID(1);
    PREPOUT_STAT;
    RUNTIME_CUBIC;
    OMPTASK_FOR(icol, colstart, colend, private(dplocx,nplocx,dplocx0,nplocy,dplocy0,dplocy)){	
	dplocy=myfma(ostat->cols[icol].ystart,dy_in2,displacey);
	if(dplocy>=nymin && dplocy<=nymax){
	    SPLIT(dplocy,dplocy,nplocy);
	    dplocy0=1.-dplocy;
	    const int nxout=ostat->cols[icol+1].pos-ostat->cols[icol].pos;
	    for(int ix=0; ix<nxout; ix++){
		int iloc=ostat->cols[icol].pos+ix;//output index.
		dplocx=myfma(ostat->cols[icol].xstart+ix*dxout,dx_in2,displacex); 
		if(dplocx>=nxmin && dplocx<=nxmax){
		    SPLIT(dplocx,dplocx,nplocx);
		    dplocx0=1.-dplocx;
		    MAKE_CUBIC_COEFF;
		    CUBIC_ADD_GRID;
		}else{
		    missing++;
		}
	    }
	}else{
	    missing++;
	}
	OMPTASK_END;
	WARN_MISSING;
    }
}
/**
   like prop_nongrid() but with cubic influence functions. cubic_iac is the
   inter-actuator coupling. Consider embed the input into a map_t and call
   prop_grid_cubic instead, which is must faster.
*/
void prop_nongrid_cubic(ARGIN_NONGRID,
			ARGOUT_LOC,
			ARG_PROP,
			double cubic_iac,
			long start, long end){
    PREPIN_NONGRID(1);
    PREPOUT_LOC;
    RUNTIME_CUBIC;
    OMPTASK_FOR(iloc, start, end, private(dplocx,dplocy,nplocx,nplocy,dplocx0,dplocy0)){
	dplocy=myfma(py[iloc],dy_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);
	if(dplocy>=nymin && dplocy<=nymax && dplocx>=nxmin && dplocx<=nxmax){
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    dplocy0=1.-dplocy;
	    dplocx0=1.-dplocx;
	    MAKE_CUBIC_COEFF;
	    CUBIC_ADD_NONGRID;
	}else{
	    missing++;
	}
    }
    OMPTASK_END;
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
    PREPIN_NONGRID(1);
    PREPOUT_PTS;
    RUNTIME_CUBIC;
    OMPTASK_FOR(isa, start, end, private(dplocx,dplocy,dplocx0,dplocy0,nplocx,nplocy)){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	    if(dplocy>=nymin && dplocy<=nymax){
		SPLIT(dplocy,dplocy,nplocy);
		dplocy0=1.-dplocy;
		for(int ix=0; ix<pts->nx; ix++){
		    iloc++;
		    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		    if(dplocx>=nxmin && dplocx<=nxmax){
			SPLIT(dplocx,dplocx,nplocx);
		        dplocx0=1.-dplocx;
			MAKE_CUBIC_COEFF;
		        CUBIC_ADD_NONGRID;
		    }else{
			missing++;
		    }
		}
	    }else{
		missing++;
	    }
	}
    }
    OMPTASK_END;
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
    PREPIN_NONGRID(1);
    PREPOUT_MAP;
    RUNTIME_CUBIC;
    OMPTASK_FOR(iy, start, end, private(dplocx,nplocx,dplocx0,nplocy,dplocy0,dplocy)){
	dplocy=myfma(oy+iy*dyout,dy_in2,displacey);
	if(dplocy>=nymin && dplocy<=nymax){
	    SPLIT(dplocy,dplocy,nplocy);
	    dplocy0=1.-dplocy;
	    for(int ix=0; ix<nxout; ix++){
		int iloc=ix+iy*nxout;
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		if(dplocx>=nxmin && dplocx<=nxmax){
		    SPLIT(dplocx,dplocx,nplocx);
		    dplocx0=1.-dplocx;
		    MAKE_CUBIC_COEFF;
		    CUBIC_ADD_NONGRID;
		}else{
		    missing++;
		}
	    }
	}else{
	    missing++;
	}
    }
    OMPTASK_END;
    //WARN_MISSING;//will always be triggered
}

/**
   The following routine is used to do down sampling by doing binning ray
   tracing using reverse interplation.  locout is coarse sampling, locin is fine
   sampling. phiout is the destination OPD. The weightings are obtained by
   interpolating from locout to locin, but the OPD are reversed computed. Simply
   replace prop_nongrid by prop_nongrid_bin without changing arguments, except
   removing start, end, will do the same ray tracing using reverse interpolation
   (binning). 

   2011-04-27: Revised usage of alpha, displacex/y so that the result agrees with
   prop_nongrid when used in the same situation. Report missing does not make
   sense here since locin is usually bigger than locout.

   2014-06-04: Notice that this operation is not direct transpose of
   prop_nongrid() to avoid accumulate bumps around the edge.
*/
void prop_nongrid_bin(const loc_t *locin,
		      const double* phiin,
		      loc_t *locout, 
		      double* phiout,
		      ARG_PROP){
    if(locout->dx<locin->dx) {
	error("This routine is designed for down sampling.\n");
    }
    loc_create_map_npad(locout, 1,0,0);
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    const int nxmin=locout->npad; 
    const int nymin=locout->npad; 
    const int nxmax=locout->map->nx-nxmin-1; 
    const int nymax=locout->map->ny-nymin-1; 
    /*notice inverse of scale. */
    const double dx_out1 = 1./locout->dx;
    const double dx_out2 = (1./scale)*dx_out1;
    const double dy_out1 = 1./locout->dy;
    const double dy_out2 = (1./scale)*dy_out1;
    /*notice negative sign in displacex/y. */
    displacex = (-displacex/scale-locout->map->ox)*dx_out1;
    displacey = (-displacey/scale-locout->map->oy)*dy_out1;
    const double *px=locin->locx;
    const double *py=locin->locy;
    /*Scale alpha to cancel out scaling */
    double alpha2 = alpha*(locin->dx/locout->dx/scale)*(locin->dy/locout->dy/scale);
    long iphi1,iphi2,iphi3,iphi4;
    map_t* map=locout->map;
    /*-1 because we count from 1 in the map. */
    double *phiout0=phiout-1;
    for(long iloc=0; iloc<locin->nloc; iloc++){
	dplocy=myfma(py[iloc],dy_out2,displacey);
	dplocx=myfma(px[iloc],dx_out2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy>=nymin && nplocy<nymax && nplocx>=nxmin && nplocx<nxmax){
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	
	    /*only proceed if all four points exist, otherwise we run the risk
	     * of accumulating edge results beyond the map to within the map. */
	    iphi1=IND(map,nplocx,nplocy);
	    iphi2=IND(map,nplocx1,nplocy);
	    iphi3=IND(map,nplocx,nplocy1);
	    iphi4=IND(map,nplocx1,nplocy1);
	    if(iphi1>0 && iphi2>0 && iphi3>0 && iphi4>0){
		phiout0[iphi1]+=alpha2*(phiin[iloc]*(1.-dplocx)*(1.-dplocy));
		phiout0[iphi2]+=alpha2*(phiin[iloc]*(dplocx)*(1.-dplocy));
		phiout0[iphi3]+=alpha2*(phiin[iloc]*(1.-dplocx)*(dplocy));
		phiout0[iphi4]+=alpha2*(phiin[iloc]*(dplocx)*(dplocy));
	    }
	}
    }
}
