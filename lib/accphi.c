/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
	prop_grid_grid(propdata->mapin, propdata->mapout,
		       propdata->alpha, displacex, displacey,
		       propdata->scale, propdata->wrap);
	break;
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
	if(propdata->mapout){
	    if(done) error("Invalid\n");
	    //case 0
	    propdata->index=0;
	    done=1;
	}
	if(propdata->phiout){
	    if(propdata->ptsout){
		if(done) error("Invalid\n");
		//case 1
		propdata->index=1;
		done=1;
	    }
	    if(propdata->locout){
		if(done) error("Invalid\n");
		//case 2
		propdata->index=2;
		done=1;
	    }
	    if(propdata->ostat){
		if(done) error("Invalid\n");
		//case 3
		propdata->index=3;
		done=1;
	    }
	}
    }
    if(propdata->locin){
	if(propdata->cubic){
	    if(propdata->locout){
		if(done) error("Invalid\n");
		//case 4
		propdata->index=4;
		done=1;
	    }
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		//case 5
		propdata->index=5;
		done=1;
	    }
	    if(propdata->ptsout){
		if(done) error("Invalid\n");
		//case 6
		propdata->index=6;
		done=1;
	    }
	}else{
	    if(propdata->locout){
		if(done) error("Invalid\n");
		//case 7
		propdata->index=7;
		done=1;
	    }
	    if(propdata->mapout){
		if(done) error("Invalid\n");
		//case 8
		propdata->index=8;
		done=1;
	    }
	    if(propdata->ptsout){
		if(done) error("Invalid\n");
		//case 9
		propdata->index=9;
		done=1;
	    }
	}
    }
    if(done==0) error("Invalid\n");
}
/**
   Propagate OPD defines on grid mapin to grid mapout.  alpha is the scaling of
   data. displacex, displacy is the displacement of the center of the beam on
   the input grid. scale is the cone effect.*/
void prop_grid_grid(const map_t *mapin, /**<[in] OPD defind on a square grid*/
		    map_t *mapout,      /**<[in,out] OPD defind on a square grid*/
		    double alpha,       /**<[in] scaling of OPD*/
		    double displacex,   /**<[in] displacement of the ray */
		    double displacey,   /**<[in] displacement of the ray*/
		    double scale,       /**<[in] scaling of the beam diameter (cone)*/
		    int wrap            /**<[in] wrap input OPD or not*/
		    ){
    //A convenient function. Not optimized
    pts_t pts;
    pts.nsa=1;
    pts.origx=&mapout->ox;
    pts.origy=&mapout->oy;
    pts.dx=mapout->dx;
    pts.nx=mapout->nx;
    prop_grid_pts(mapin, &pts, mapout->p, alpha, displacex, displacey,
		  scale, wrap, 0, 1);
}
/**
   Propagate OPD defines on grid mapin to subapertures.  alpha is the scaling of
   data. displacex, displacy is the displacement of the center of the beam on
   the input grid.  scale is the cone effect.*/
void prop_grid_pts(const map_t *mapin, /**<[in] OPD defind on a square grid*/
		   const pts_t *pts,   /**<[in] defining each subaperture*/
		   double *phiout0,    /**<[in,out] OPDs for each subaperture*/
		   double alpha,       /**<[in] scaling of OPD*/
		   double displacex,   /**<[in] displacement of the beam */
		   double displacey,   /**<[in] displacement of the beam */
		   double scale,       /**<[in] scaling of the beam diameter (cone)*/
		   int wrap,           /**<[in] wrap input OPD or not*/
		   long sastart,       /**<[in] The starting subaperture to trace ray*/
		   long saend         /**<[in] The last (exclusive) subaperture to trace ray*/
		   ){
    /*
       2010-01-02: Improved by saving x interpolations for each subaperture
       during non-matched case. original function renamed to prop_grid_pts_old
       time cut by almost 3 fold for each LGS WFS.  0.40s->0.16s
     */
    int isa, ipix,jpix;
    int mx,my;
    int nplocx2;
    const int ninx   = mapin->nx;
    const int niny   = mapin->ny;
    const int nx     = pts->nx;
    const double dx_in1 = 1/mapin->dx;
    const double dxout  = pts->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dx_in1;
    if(!saend) saend=pts->nsa;
    if(USE_OPTIM && fabs(dx_in2*dxout-1)<EPS){
	double dplocx, dplocy;
	int nplocx, nplocy;
        for(isa=sastart; isa<saend; isa++){
	    double (*phiout)[nx]=(double(*)[nx])(phiout0+nx*nx*isa);
	    double (*phiin)[ninx]=(double(*)[ninx])(mapin->p);
	    const double origx=pts->origx[isa];
	    const double origy=pts->origy[isa];
	    double w11,w10,w01,w00;
	    
	    dplocy = myfma(origy,dx_in2,displacey);
	    dplocx = myfma(origx,dx_in2,displacex);
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(!wrap){
		int sx, sy;
		if(nplocy<0){
		    sy=-nplocy;
		}else
		    sy=0;
		if(nplocx<0){
		    sx=-nplocx;
		}else
		    sx=0;
		
		my=niny-nplocy-1;//remaining possible points.
		if(my>nx){
		    my=nx;
		}
		if(my<=sy){
		    continue;
		}
		mx=ninx-nplocx-1;
		if(mx>nx){
		    mx=nx;
		}
		if(mx<=sx){
		    continue;
		}

		if((dplocx)<EPS && (dplocy)<EPS){
		    /*aligned perfectly.*/
		    for(jpix=sy; jpix<my; jpix++){
			//nplocy2=nplocy+jpix;
			double *restrict phiin2=phiin[nplocy+jpix];
			double *restrict phiout2=phiout[jpix];
			for(ipix=sx; ipix<mx; ipix++){
			    phiout2[ipix]+=alpha*phiin2[nplocx+ipix];
			}
		    }
		}else{	
		    w11=dplocx*dplocy;
		    w10=(1.-dplocx)*dplocy;
		    w01=dplocx*(1.-dplocy);
		    w00=(1.-dplocx)*(1.-dplocy);	
		    for(jpix=sy; jpix<my; jpix++){
			//nplocy2=nplocy+jpix;
			double *restrict phiin2=phiin[nplocy+jpix];
			double *restrict phiin3=phiin[nplocy+jpix+1];
			double *restrict phiout2=phiout[jpix];
			for(ipix=sx; ipix<mx; ipix++){
			    phiout2[ipix]+=
				alpha*(+phiin2[nplocx+ipix]*w00
				       +phiin2[nplocx+ipix+1]*w01
				       +phiin3[nplocx+ipix]*w10
				       +phiin3[nplocx+ipix+1]*w11);
			}
		    }
		}
	    }else{//wraping
		double *phiin_1, *phiin_2;

		if(ninx < nx || niny < nx){
		    error("Input map is too small. wraps more than once\n");
		}
		w11=dplocx*dplocy;
		w10=(1.-dplocx)*dplocy;
		w01=dplocx*(1.-dplocy);
		w00=(1.-dplocx)*(1.-dplocy);

		nplocy-=niny*(nplocy/niny);
		if(nplocy<0)
		    nplocy+=niny;
		nplocx-=ninx*(nplocx/ninx);
		if(nplocx<0)
		    nplocx+=ninx;
		my=niny-nplocy-1;//remaining possible points.
		mx=ninx-nplocx-1;
		if(my>nx) my=nx;
		if(mx>nx) mx=nx;

		for(jpix=0; jpix<nx; jpix++){
		    if(jpix<my){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1];
		    }else if(jpix==my){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1-niny];
		    }else{
			phiin_1=phiin[nplocy-niny];
			phiin_2=phiin[nplocy+1-niny];
		    }
		    nplocx2=nplocx;
		    for(ipix=0; ipix<mx; ipix++){
			phiout[jpix][ipix]+=alpha*
			    (phiin_1[nplocx2]*w00
			     +phiin_1[nplocx2+1]*w01
			     +phiin_2[nplocx2]*w10
			     +phiin_2[nplocx2+1]*w11);
			nplocx2++;
		    }
		    if(mx<nx){
			ipix=mx;
			phiout[jpix][ipix]+=alpha*
			    (phiin_1[nplocx2]*w00
			     +phiin_1[nplocx2+1-ninx]*w01
			     +phiin_2[nplocx2]*w10
			     +phiin_2[nplocx2+1-ninx]*w11);
			nplocx2++;

			nplocx2-=ninx;
			for(ipix=mx+1; ipix<nx; ipix++){
			    phiout[jpix][ipix]+=alpha*
				(phiin_1[nplocx2]*w00
				 +phiin_1[nplocx2+1]*w01
				 +phiin_2[nplocx2]*w10
				 +phiin_2[nplocx2+1]*w11);
			    nplocx2++;
			}
		    }
		    nplocy++;
		}
	    }	    
	}
    }else{ /*different spacing. or non optim*/
	double dplocx, dplocy;
	double ratio;
	const double (*phiin)[ninx]=(const double(*)[ninx])(mapin->p);
	int nplocx, nplocy;
	ratio = dxout*dx_in2;
	for(isa=sastart; isa<saend; isa++){
	    double (*phiout)[nx]=(double(*)[nx])(phiout0+nx*nx*isa);
	    const double origx=pts->origx[isa];
	    const double origy=pts->origy[isa];
	    double dplocx0, dplocy0;
	    dplocy0 = myfma(origy,dx_in2,displacey);
	    dplocx0 = myfma(origx,dx_in2,displacex);
	    if(!wrap){
		int sx, sy;
		if(dplocx0<0){
		    sx=iceil(-dplocx0/ratio);
		}else
		    sx=0;
		if(dplocy0<0){
		    sy=iceil(-dplocy0/ratio);
		}else
		    sy=0;

		my=iceil((niny-1-dplocy0)/ratio);
		if(my>nx){
		    my=nx;
		}
		mx=iceil((ninx-1-dplocx0)/ratio);
		if(mx>nx){
		    mx=nx;
		}
		
		int nplocxs[mx];
		double dplocxs[mx];

		dplocy0  = myfma(origy+(double)sy*dxout,dx_in2,displacey);
		dplocx0  = myfma(origx+(double)sx*dxout,dx_in2,displacex);
		
		for(ipix=sx; ipix<mx; ipix++){
		    SPLIT(dplocx0,dplocx,nplocx);
		    nplocxs[ipix]=nplocx;
		    dplocxs[ipix]=dplocx;
		    dplocx0+=ratio;
		}

		for(jpix = sy; jpix < my; jpix++){
		    SPLIT(dplocy0,dplocy,nplocy);
		    const double dplocy1=1.-dplocy;
		    const double *phiin2=phiin[nplocy];
		    const double *phiin3=phiin[nplocy+1];
		    double *phiout2=phiout[jpix];
		    for(ipix=sx; ipix<mx; ipix++){
			nplocx=nplocxs[ipix];
			dplocx=dplocxs[ipix];
			phiout2[ipix]+=alpha*
			    (+(phiin2[nplocx]
			       +(phiin2[nplocx+1]-phiin2[nplocx])*dplocx)
			     *dplocy1
			     +(phiin3[nplocx]
			       +(phiin3[nplocx+1]-phiin3[nplocx])*dplocx)
			     *dplocy);
		    }
		    dplocy0+=ratio;
		}
	    }else{
		const double *phiin_1, *phiin_2;
		double dplocy1;
		int mx0;
		if(ninx < nx*ratio || niny < nx*ratio){
		    info("nx=%d, ratio=%g, ninx=%d\n",nx,ratio,ninx);
		    error("Input map is too small. wraps more than once\n");
		}
		dplocy0-=niny*ifloor(dplocy0/niny);
		dplocx0-=ninx*ifloor(dplocx0/ninx);
		mx0=iceil((ninx-1.-dplocx0)/ratio);
		if(mx0>nx) mx0=nx;
		if(mx0<0) mx0=0;
		
		int nplocxs[nx];
		int nplocxs2[nx];
		double dplocxs[nx];
		mx=mx0;
		for(ipix=0; ipix<mx; ipix++){
		    SPLIT(dplocx0,dplocx,nplocx);
		    nplocxs[ipix]=nplocx;
		    nplocxs2[ipix]=nplocx+1;
		    dplocxs[ipix]=dplocx;
		    dplocx0+=ratio;
		}
		if(mx<nx){
		    while(dplocx0<ninx && mx<nx){//falls on the edge
			SPLIT(dplocx0,dplocx,nplocx);
			nplocxs[mx]=nplocx;
			nplocxs2[mx]=0;
			dplocxs[mx]=dplocx;
			dplocx0+=ratio;
			mx++;
		    }
		    dplocx0-=ninx;
		    for(ipix=mx; ipix<nx; ipix++){
			SPLIT(dplocx0,dplocx,nplocx);
			nplocxs[ipix]=nplocx;
			nplocxs2[ipix]=nplocx+1;
			dplocxs[ipix]=dplocx;
			dplocx0+=ratio;
		    }
		}
		for(jpix=0; jpix<nx; jpix++){
		    mx=mx0;
		    SPLIT(dplocy0,dplocy,nplocy);
		    dplocy1=1.-dplocy;
		    if(nplocy<niny-1){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1];
		    }else if(nplocy==niny-1){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[0];
		    }else{
			phiin_1=phiin[nplocy-niny];
			phiin_2=phiin[nplocy+1-niny];
		    }
		    double *phiout2=phiout[jpix];
		    for(ipix=0; ipix<nx; ipix++){
			nplocx=nplocxs[ipix];
			nplocx2=nplocxs2[ipix];
			dplocx=dplocxs[ipix];
			phiout2[ipix]+=alpha*
			    (+(phiin_1[nplocx]
			       +(phiin_1[nplocx2]-phiin_1[nplocx])*dplocx)
			     *dplocy1
			     +(phiin_2[nplocx]
			       +(phiin_2[nplocx2]-phiin_2[nplocx])*dplocx)
			     *dplocy);
			
		    }
		    dplocy0+=ratio;
		}
	    }/*wrap*/
	}/*end isa*/
    }
}
/**
   Propagate OPD defines on grid mapin to locstat_t that stores the starting
   point of each column.  alpha is the scaling of data. displacex, displacy is
   the displacement of the center of the beam on the input grid.  scale is the
   cone effect.*/
void prop_grid_stat(const map_t *mapin, /**<[in] OPD defind on a square grid*/
		    const locstat_t *ostat, /**<[in] information about each clumn of a output loc grid*/
		    double *phiout0,    /**<[in,out] OPD defined on ostat*/
		    double alpha,       /**<[in] scaling of OPD*/
		    double displacex,   /**<[in] displacement of the ray */
		    double displacey,   /**<[in] displacement of the ray */
		    double scale,       /**<[in] scaling of the beam diameter (cone)*/
		    int wrap,           /**<[in] wrap input OPD or not*/
		    long colstart,      /**<[in] First column to do ray tracing*/
		    long colend         /**<[in] Last column (exclusive) to do ray tracing*/
		    ){
    double *phiout;
    double dplocx, dplocy;
    int nplocx, nplocy;
    int icol,irow;
    long offset;
    int collen;
    int missing=0;
    if(colend==0) colend = ostat->ncol;

    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx  = wrapx1 - 1;
    const int wrapy  = wrapy1 - 1;
    
    const double dx_in1 = 1./mapin->dx;
    const double dx_in2 = scale*dx_in1;
    const double dxout  = ostat->dx;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dx_in1;
    const double *phiin  = mapin->p;
    const double ratio  = dxout*dx_in2;
#if USE_OPTIM == 1
	if(fabs(ratio-1)<EPS){
	    /*grid size of loc_in and loc_out agree*/
	    double bl, br, tl, tr;
	    const double *phicol, *phicol2;
	    int rowdiv,rowdiv2;
	    int irows;
	    /*loc_out and loc_in has the same grid sampling.*/
	    for(icol=colstart; icol<colend; icol++){
		/*starting address of that col*/
		offset=ostat->cols[icol].pos;
		collen=ostat->cols[icol+1].pos-offset;/*exclusive*/

		phiout=phiout0+offset;
		dplocy=ostat->cols[icol].ystart*dx_in2+displacey;
		if(wrap){
		    dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
		}
		SPLIT(dplocy,dplocy,nplocy);
		dplocx=ostat->cols[icol].xstart*dx_in2+displacex;
		if(wrap){
		    dplocx-=wrapx1*floor(dplocx/(double)wrapx1);
		    irows=0;
		}else{
		    if(dplocx<0)
			irows=iceil(-dplocx);
		    else
			irows=0;
		}
		SPLIT(dplocx,dplocx,nplocx);

		if(wrap){
		    phicol  = phiin+nplocy*wrapx1+nplocx;
		    phicol2 = phiin+(nplocy==wrapy?0:nplocy+1)*wrapx1+nplocx;
		}else{
		    if(nplocy<wrapy && nplocy>=0){
			phicol=phiin+nplocy*wrapx1+nplocx;
			phicol2=phiin+(nplocy+1)*wrapx1+nplocx;
		    }else{
			missing+=collen;
			continue;
		    }
		}
 
	    
		bl=dplocx*dplocy;
		br=(1.-dplocx)*dplocy;
		tl=dplocx*(1-dplocy);
		tr=(1.-dplocx)*(1-dplocy);
	    
		rowdiv=wrapx-nplocx;
		/*max number of rows possible before wraping*/
		if(rowdiv>collen) rowdiv=collen;
		if(rowdiv<0)rowdiv=0;
		/*reduce to number of needed*/

		for(irow=irows; irow<rowdiv; irow++){
		    phiout[irow]+=alpha*
			(bl*phicol2[irow+1]
			 +br*phicol2[irow]
			 +tl*phicol[irow+1]
			 +tr*phicol[irow]);
		}

		if(wrap){
		    while(rowdiv < collen){/*can wrap several times*/ 
			irow=rowdiv;
			phiout[irow]+=alpha*
			    (bl*phicol2[irow-wrapx]
			     +br*phicol2[irow]
			     +tl*phicol[irow-wrapx]
			     +tr*phicol[irow]);
			rowdiv++;
			
			rowdiv2=rowdiv+wrapx;
			if(rowdiv2>collen) rowdiv2=collen;

			phicol-=wrapx1;/*wrap the second part*/
			phicol2-=wrapx1;

			for(irow=rowdiv; irow<rowdiv2; irow++){
			    phiout[irow]+=alpha*
				(+bl*phicol2[irow+1]
				 +br*phicol2[irow]
				 +tl*phicol[irow+1]
				 +tr*phicol[irow]);
			}
			rowdiv=rowdiv2;
		    }
		}
	    }/*end for icol*/
	}else{
	    /*grid size of loc_in and loc_out doesn't agree*/
	    const double *phicol, *phicol2;
	    double dplocx0;
	    int nplocx0;
	    int rowdiv,rowdiv2;
	    int irows;
	    for(icol=colstart; icol<colend; icol++){
		/*starting address of that col*/
		offset=ostat->cols[icol].pos;
		collen=ostat->cols[icol+1].pos-offset;/*exclusive*/

		phiout=phiout0+offset;
		dplocy=ostat->cols[icol].ystart*dx_in2+displacey;
		if(wrap){
		    dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
		}
		SPLIT(dplocy,dplocy,nplocy);
		const double dplocy1=1.-dplocy;
		dplocx=ostat->cols[icol].xstart*dx_in2+displacex;
		if(wrap){
		    dplocx-=wrapx1*floor(dplocx/(double)wrapx1);
		    irows=0;
		}else{
		    if(dplocx<0){
			irows=iceil(-dplocx/ratio);
		    }else
			irows=0;
		}

		if(wrap){
		    phicol  = phiin+nplocy*wrapx1;
		    phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*wrapx1;
		}else{
		    if(nplocy<wrapy && nplocy>=0){
			phicol=phiin+nplocy*wrapx1;
			phicol2=phiin+(nplocy+1)*wrapx1;
		    }else{
			missing+=collen;
			continue;
		    }
		}
 	  
		/*last row to do if not wrap*/
		/*figure out which row we can go in this segment. */
		  
		rowdiv = iceil((wrapx-dplocx)/ratio);
		if(rowdiv<0) rowdiv=0;//rowdiv may be -1 if dplocx=wrapx+0.*
		if(rowdiv>collen) rowdiv=collen;
		dplocx = dplocx + ratio*irows;
		for(irow=irows; irow<rowdiv; irow++){/*no wrap*/
		    SPLIT(dplocx,dplocx0,nplocx0);
		    phiout[irow]+=alpha*
			(+dplocy1*((1-dplocx0)*phicol[nplocx0]
				   +dplocx0*phicol[nplocx0+1])
			 +dplocy *((1-dplocx0)*phicol2[nplocx0]
				   +dplocx0*phicol2[nplocx0+1]));
		    dplocx+=ratio;
		}
		if(wrap){
		    while(rowdiv < collen){/*can wrap several times*/ 
			//falls on the edge.
			while(dplocx<wrapx1 && rowdiv<collen){
			    SPLIT(dplocx,dplocx0,nplocx0);
			    phiout[rowdiv]+=alpha*
				(+dplocy1*((1-dplocx0)*phicol[nplocx0]
					   +dplocx0*phicol[nplocx0-wrapx])
				 +dplocy *((1-dplocx0)*phicol2[nplocx0]
					   +dplocx0*phicol2[nplocx0-wrapx]));
			    dplocx+=ratio;
			    rowdiv++;
			}
			dplocx=dplocx-wrapx1;			
			rowdiv2=iceil((wrapx-dplocx)/ratio);
			if(rowdiv2>collen) rowdiv2=collen;
			for(irow=rowdiv; irow<rowdiv2; irow++){
			    SPLIT(dplocx,dplocx0,nplocx0);
			    phiout[irow]+=alpha*
				(+dplocy1*((1-dplocx0)*phicol[nplocx0]
					   +dplocx0*phicol[nplocx0+1])
				 +dplocy *((1-dplocx0)*phicol2[nplocx0]
					   +dplocx0*phicol2[nplocx0+1]));
			    dplocx+=ratio;
			}
			rowdiv=rowdiv2;
		    }
		}
	    }/*end for icol*/
	}/*fabs(dx_in2*dxout-1)<EPS*/
#else
	double dplocx0;
	int nplocy1, nplocx1;

	/*non optimized case. slower, but hopefully accurate*/
	for(icol=colstart; icol<colend; icol++){
	    offset=ostat->cols[icol].pos;
	    collen=ostat->cols[icol+1].pos-offset;
	    phiout=phiout0+offset;
	    dplocy=ostat->cols[icol].ystart*dx_in2+displacey;
	    SPLIT(dplocy,dplocy,nplocy);
	    const double dplocy1=1.-dplocy;
	    if(nplocy<0||nplocy>=wrapy){
		if(wrap){
		    while(nplocy<0)
			nplocy+=wrapy1;
		    while(nplocy>wrapy)
			nplocy-=wrapy1;
		    nplocy1=(nplocy==wrapx?0:nplocy+1);
		}else{
		    missing+=collen;
		    continue;
		}
	    }else{
		nplocy1=nplocy+1;
	    }
	    dplocx0=(ostat->cols[icol].xstart)*dx_in2+displacex;
	    for(irow=0; irow<collen; irow++){
		SPLIT(dplocx0,dplocx,nplocx);
		dplocx0+=ratio;
		if(nplocx>=0 && nplocx<wrapx){
		    nplocx1=nplocx+1;
		}else{
		    if(wrap){
			while(nplocx<0)
			    nplocx+=wrapx1;
			while(nplocx>wrapx)
			    nplocx-=wrapx1;
			nplocx1=(nplocx==wrapx?0:nplocx+1);
		    }else{
			missing++;
			continue;
		    }
		}
		
		phiout[irow]+= alpha*
		    (dplocx * (dplocy * phiin[(nplocx1) + (nplocy1)*wrapx1]
			       +dplocy1 * phiin[(nplocx1) + nplocy*wrapx1])
		     + (1-dplocx) * (dplocy * phiin[nplocx + (nplocy1)*wrapx1]
				     +dplocy1 * phiin[nplocx + nplocy*wrapx1]));
		
	    }/*for irow*/
	}/*for icol*/
#endif
    if(missing>0){
	warning("%d points not covered by input screen\n", missing);
    }
}/*function*/
/**
   Propagate OPD defines on grid mapin to coordinate locout.  alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect.*/
void prop_grid(const map_t *mapin, /**<[in] OPD defind on a square grid*/
	       const loc_t *locout,/**<[in] coordinate of irregular destination grid*/
	       double *phiout,     /**<[in,out] OPD defined on locout*/
	       double alpha,       /**<[in] scaling of OPD*/
	       double displacex,   /**<[in] displacement of the ray */
	       double displacey,   /**<[in] displacement of the ray */
	       double scale,       /**<[in] wrap input OPD or not*/
	       int wrap,           /**<[in] wrap input OPD or not*/
	       long start,         /**<[in] First point to do*/
	       long end            /**<[in] Last point to do*/
	       ){
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    int iloc;
    int missing=0;
    if(!locout) error("locout is NULL!");
    if(end==0){
	end=locout->nloc;
    }
    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx=wrapx1-1;
    const int wrapy=wrapy1-1;
    const double dx_in1 = 1./mapin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;

    
    double (*phiin)[wrapx1]=(double(*)[wrapx1])(mapin->p);
    for(iloc=start; iloc<end; iloc++){
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);
	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    if(wrap){
		while(nplocx<0)
		    nplocx+=wrapx1;
		while(nplocx>wrapx)
		    nplocx-=wrapx1;
		while(nplocy<0)
		    nplocy+=wrapy1;
		while(nplocy>wrapy)
		    nplocy-=wrapy1;
		
		nplocx1=(nplocx==wrapx?0:nplocx+1);
		nplocy1=(nplocy==wrapy?0:nplocy+1);
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
    if(missing>0){
	warning("%d points not covered by input screen\n", missing);
    }
}
/**
   Propagate OPD defines on coordinate locin to coordinate locout.  This is the
   <em>slowest</em>. alpha is the scaling of data. displacex, displacy is the
   displacement of the center of the beam on the input grid.  scale is the cone
   effect. See prop_grid() for definition of other parameters.*/
void prop_nongrid(loc_t *locin,        /**<[in] Coordinate of iregular source grid*/
		  const double* phiin, /**<[in] Input OPD defined in locin*/
		  const loc_t *locout, /**<[in] Coordinate of irregular output grid*/
		  const double *ampout,/**<[in] Amplitude defined on locout. skip point of amp is 0*/
		  double* phiout,      /**<[in,out] Output OPD defined in locout*/
		  double alpha,        /**<[in] scaling of OPD*/
		  double displacex,    /**<[in] displacement of the ray */
		  double displacey,    /**<[in] displacement of the ray */
		  double scale,        /**<[in] wrap input OPD or not*/
		  long start,          /**<[in] First point to do*/
		  long end             /**<[in] Last point to do*/
		  ){
    loc_create_map_npad(locin,1);
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    long iloc;
    int missing=0;
 
    const int wrapx1 = locin->map->nx;
    const int wrapy1 = locin->map->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;

#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    //-1 because we count from 1 in the map.
    const double *phiin0=phiin-1;
    if(!end) end=locout->nloc;
    for(iloc=start; iloc<end; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;//skip points that has zero amplitude
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    missing++;
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}

#if ONLY_FULL == 1 //only proceed if all four points exist.
	iphi1=map[nplocy][nplocx];
	iphi2=map[nplocy][nplocx1];
	iphi3=map[nplocy1][nplocx];
	iphi4=map[nplocy1][nplocx1];
	if(iphi1 && iphi2 && iphi3 && iphi4){
	    phiout[iloc]+=alpha*(phiin0[iphi1]*(1.-dplocx)*(1.-dplocy));
	    phiout[iloc]+=alpha*(phiin0[iphi2]*(dplocx)*(1.-dplocy));
	    phiout[iloc]+=alpha*(phiin0[iphi3]*(1.-dplocx)*(dplocy));
	    phiout[iloc]+=alpha*(phiin0[iphi4]*(dplocx)*(dplocy));
	}else{
	    missing++;
	}
#else	
	if((iphi=map[nplocy][nplocx]))
	    phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(1.-dplocy));
	if((iphi=map[nplocy][nplocx1]))
	    phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(1.-dplocy));
	if((iphi=map[nplocy1][nplocx]))
	    phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(dplocy));
	if((iphi=map[nplocy1][nplocx1]))
	    phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(dplocy));
#endif
    }
    if(missing>0){
	warning("%d points not covered by input screen\n", missing);
    }
}
/**
   Propagate OPD defines on coordinate locin to grid mapout. alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect. */
void prop_nongrid_map(loc_t *locin,     /**<[in] Coordinate of iregular source grid*/
		      const double *phiin,/**<[in] Input OPD defined in locin*/
		      map_t *mapout,    /**<[in,out] Output OPD defined in a square grid*/
		      double alpha,     /**<[in] scaling of OPD*/
		      double displacex, /**<[in] displacement of the ray */
		      double displacey, /**<[in] displacement of the ray */
		      double scale,     /**<[in] wrap input OPD or not*/
		      long start,       /**<[in] First point to do*/
		      long end          /**<[in] Last point to do*/
		      ){
    //propagate to a square map.
    loc_create_map_npad(locin,2);//will only do once and save in locin.

    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    long iloc;

    const int wrapx1 = locin->map->nx;
    const int wrapy1 = locin->map->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double dxout=mapout->dx;
    const double ox=mapout->ox;
    const double oy=mapout->oy;
    double *phiout=mapout->p;
#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    //-1 because we count from 1 in the map.
    long (*map)[locin->map->nx]=(long(*)[locin->map->nx])(locin->map->p);
    const int nxout=mapout->nx;
    const double *phiin0=phiin-1;
    if(!end) end=mapout->ny;
    for(int iy=start; iy<end; iy++){
	dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<0||nplocy>=wrapy){
	    continue;
	}else{
	    nplocy1=nplocy+1;
	}
	for(int ix=0; ix<nxout; ix++){
	    iloc=ix+iy*nxout;
	    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
	    SPLIT(dplocx,dplocx,nplocx);
	    if(nplocx<0||nplocx>=wrapx){
		continue;
	    }else{
		nplocx1=nplocx+1;
	    }

#if ONLY_FULL == 1 //only proceed if all four points exist.
	    iphi1=map[nplocy][nplocx];
	    iphi2=map[nplocy][nplocx1];
	    iphi3=map[nplocy1][nplocx];
	    iphi4=map[nplocy1][nplocx1];
	    if(iphi1 && iphi2 && iphi3 && iphi4){
		phiout[iloc]+=alpha*(phiin0[iphi1]*(1.-dplocx)*(1.-dplocy));
		phiout[iloc]+=alpha*(phiin0[iphi2]*(dplocx)*(1.-dplocy));
		phiout[iloc]+=alpha*(phiin0[iphi3]*(1.-dplocx)*(dplocy));
		phiout[iloc]+=alpha*(phiin0[iphi4]*(dplocx)*(dplocy));
	    }
#else	
	    if((iphi=map[nplocy][nplocx])) 
		phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(1.-dplocy));
	    if((iphi=map[nplocy][nplocx1]))
		phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(1.-dplocy));
	    if((iphi=map[nplocy1][nplocx]))
		phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(dplocy));
	    if((iphi=map[nplocy1][nplocx1]))
		phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(dplocy));
#endif
	}
    }
}
/**
   Propagate OPD defines on coordinate locin to subapertures pts. alpha is the
   scaling of data. displacex, displacy is the displacement of the center of the
   beam on the input grid.  scale is the cone effect. See prop_groid().*/
void prop_nongrid_pts(loc_t *locin,         /**<[in] Coordinate of iregular source grid*/
		      const double *phiin,  /**<[in] Input OPD defined in locin*/
		      const pts_t *pts,     /**<[in] defining each subaperture*/
		      const double *ampout, /**<[in] Amplitude of subaps. skip point of amp is 0*/
		      double *phiout,       /**<[in,out] OPD for subaps*/
		      double alpha,         /**<[in] scaling of OPD*/
		      double displacex,     /**<[in] displacement of the ray */
		      double displacey,     /**<[in] displacement of the ray */
		      double scale,         /**<[in] wrap input OPD or not*/
		      long start,           /**<[in] First point to do*/
		      long end              /**<[in] Last point to do*/
		      ){
    //propagate to a square map.
    loc_create_map_npad(locin,1);//will only do once and save in locin.

    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;

    const int wrapx1 = locin->map->nx;
    const int wrapy1 = locin->map->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double dxout=pts->dx;
#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    //-1 because we count from 1 in the map.
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    const double *phiin0=phiin-1;
    if(!end) end=pts->nsa;
    for(int isa=start; isa<end; isa++){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(nplocy<0||nplocy>=wrapy){
		continue;
	    }else{
		nplocy1=nplocy+1;
	    }
	    for(int ix=0; ix<pts->nx; ix++){
		iloc++;
		if(ampout && fabs(ampout[iloc])<EPS)
		    continue;//skip points that has zero amplitude
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		SPLIT(dplocx,dplocx,nplocx);
		if(nplocx<0||nplocx>=wrapx){
		    continue;
		}else{
		    nplocx1=nplocx+1;
		}

#if ONLY_FULL == 1 //only proceed if all four points exist.
		iphi1=map[nplocy][nplocx];
		iphi2=map[nplocy][nplocx1];
		iphi3=map[nplocy1][nplocx];
		iphi4=map[nplocy1][nplocx1];
		if(iphi1 && iphi2 && iphi3 && iphi4){
		    phiout[iloc]+=alpha*(phiin0[iphi1]*(1.-dplocx)*(1.-dplocy));
		    phiout[iloc]+=alpha*(phiin0[iphi2]*(dplocx)*(1.-dplocy));
		    phiout[iloc]+=alpha*(phiin0[iphi3]*(1.-dplocx)*(dplocy));
		    phiout[iloc]+=alpha*(phiin0[iphi4]*(dplocx)*(dplocy));
		}
#else	
		if((iphi=map[nplocy][nplocx]))
		    phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(1.-dplocy));
		if((iphi=map[nplocy][nplocx1]))
		    phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(1.-dplocy));
		if((iphi=map[nplocy1][nplocx]))
		    phiout[iloc]+=alpha*(phiin0[iphi]*(1.-dplocx)*(dplocy));
		if((iphi=map[nplocy1][nplocx1]))
		    phiout[iloc]+=alpha*(phiin0[iphi]*(dplocx)*(dplocy));
#endif
	    }    
	}
    }
}
#define MAKE_CUBIC_PARAM					\
    double fx[4],fy[4];						\
    const double cubicn=1./(1.+2.*cubic_iac);			\
    const double c0=1.*cubicn;					\
    const double c1=(4.*cubic_iac-2.5)*cubicn;			\
    const double c2=(1.5-3.*cubic_iac)*cubicn;			\
    const double c3=(2.*cubic_iac-0.5)*cubicn;			\
    const double c4=(0.5-cubic_iac)*cubicn /**<For cubic interpolation.*/
/**
   like prop_nongrid() but with cubic influence functions. cubic_iac is the
   inter-actuator coupling.*/
void prop_nongrid_cubic(loc_t *locin, const double* phiin, 
			const loc_t *locout, const double *ampout,
			double* phiout, double alpha,
			double displacex, double displacey,
			double scale,
			double cubic_iac,
			long start, long end){
    loc_create_map_npad(locin,2);//padding to avoid test boundary
    double dplocx, dplocy;
    int nplocx, nplocy;
    int ix,iy;
    long iloc;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;

    //-1 because we count from 1 in the map.
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    const double *phiin0=phiin-1;
    //cubic
    MAKE_CUBIC_PARAM;
    double dplocx0, dplocy0;
    const int nmapx3=locin->map->nx-3;
    const int nmapy3=locin->map->ny-3;
    if(!end) end=locout->nloc;
    for(iloc=start; iloc<end; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;//skip points that has zero amplitude
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<1 || nplocy>nmapy3 || nplocx<1 || nplocx>nmapx3){
	    continue;
	}
	dplocy0=1.-dplocy;
	dplocx0=1.-dplocx;
	    
	fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);
	fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);
	fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);
	fx[3]=dplocx*dplocx*(c3+c4*dplocx);
	    
	fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);
	fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);
	fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);
	fy[3]=dplocy*dplocy*(c3+c4*dplocy);

	for(iy=nplocy-1; iy<nplocy+3; iy++){
	    for(ix=nplocx-1; ix<nplocx+3; ix++){
		long iphi=map[iy][ix];
		if(iphi){
		    phiout[iloc]+=alpha*fx[ix-nplocx+1]
			*fy[iy-nplocy+1]*phiin0[iphi];
		}
	    }
	}
    }
}
/**
   like prop_nongrid_pts() but with cubic influence functions. cubic_iac is the
inter-actuator-coupling.  */
void prop_nongrid_pts_cubic(loc_t *locin, const double* phiin, 
			    const pts_t *pts, const double *ampout, 
			    double* phiout, double alpha,
			    double displacex, double displacey,
			    double scale,double cubic_iac, 
			    long start, long end){
    loc_create_map_npad(locin,2);//padding to avoid test boundary
    double dplocx, dplocy;
    int nplocx, nplocy;

    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    //-1 because we count from 1 in the map.
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    const double *phiin0=phiin-1;
    //cubic
    MAKE_CUBIC_PARAM;
    double dplocx0, dplocy0;
    const int nmapx3=locin->map->nx-3;
    const int nmapy3=locin->map->ny-3;
    const double dxout=pts->dx;
    if(!end) end=pts->nsa;
    for(int isa=start; isa<end; isa++){
	const long iloc0=isa*pts->nx*pts->nx;
	const double ox=pts->origx[isa];
	const double oy=pts->origy[isa];
	for(int iy=0; iy<pts->nx; iy++){
	    long iloc=iloc0+iy*pts->nx-1;
	    dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(nplocy<1||nplocy>=nmapy3){
		continue;
	    }
	    for(int ix=0; ix<pts->nx; ix++){
		iloc++;
		if(ampout && fabs(ampout[iloc])<EPS)
		    continue;//skip points that has zero amplitude
		dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
		SPLIT(dplocx,dplocx,nplocx);
		if(nplocx<1||nplocx>=nmapx3){
		    continue;
		}
		
		dplocy0=1.-dplocy;
		dplocx0=1.-dplocx;
	    
		fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);
		fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);
		fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);
		fx[3]=dplocx*dplocx*(c3+c4*dplocx);
	    
		fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);
		fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);
		fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);
		fy[3]=dplocy*dplocy*(c3+c4*dplocy);
	
		for(int jy=nplocy-1; jy<nplocy+3; jy++){
		    for(int jx=nplocx-1; jx<nplocx+3; jx++){
			long iphi=map[jy][jx];
			if(iphi){
			    phiout[iloc]+=alpha*fx[jx-nplocx+1]
				*fy[jy-nplocy+1]*phiin0[iphi];
			}
		    }
		}
	    }
	}
    }
}
/**
   like prop_nongrid_map() but with cubic influence functions. cubic_iac is the
inter-actuator-coupling. */
void prop_nongrid_map_cubic(loc_t *locin, const double* phiin, 
			    map_t* mapout, double alpha,
			    double displacex, double displacey,
			    double scale,double cubic_iac,
			    long start, long end){
    loc_create_map_npad(locin,2);//padding to avoid test boundary
    double dplocx, dplocy;
    int nplocx, nplocy;
    long iloc;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;

    const double dxout=mapout->dx;
    const double ox=mapout->ox;
    const double oy=mapout->oy;
    double *phiout=mapout->p;
    //cubic
    MAKE_CUBIC_PARAM;
    double dplocx0, dplocy0;
    const int nmapx3=locin->map->nx-3;
    const int nmapy3=locin->map->ny-3;
    const int nxout=mapout->nx;
    //-1 because we count from 1 in the map.
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    const double *phiin0=phiin-1;
    if(!end) end=mapout->ny;
    for(int iy=start; iy<end; iy++){
	dplocy=myfma(oy+iy*dxout,dx_in2,displacey);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocy<1||nplocy>=nmapy3){
	    continue;
	}
	for(int ix=0; ix<nxout; ix++){
	    iloc=ix+iy*nxout;
	    dplocx=myfma(ox+ix*dxout,dx_in2,displacex); 
	    SPLIT(dplocx,dplocx,nplocx);
	    if(nplocx<1||nplocx>=nmapx3){
		continue;
	    }
	    dplocy0=1.-dplocy;
	    dplocx0=1.-dplocx;
	    
	    fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);
	    fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);
	    fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);
	    fx[3]=dplocx*dplocx*(c3+c4*dplocx);
	    
	    fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);
	    fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);
	    fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);
	    fy[3]=dplocy*dplocy*(c3+c4*dplocy);

	    for(int jy=nplocy-1; jy<nplocy+3; jy++){
		for(int jx=nplocx-1; jx<nplocx+3; jx++){
		    long iphi=map[jy][jx];
		    if(iphi){
			phiout[iloc]+=alpha*fx[jx-nplocx+1]
			    *fy[jy-nplocy+1]*phiin0[iphi];
		    }
		}
	    }
	}
    }
}

/**
  the following routine is used to do down sampling by doing *reverse* ray
  tracing.  locin is coarse sampling, locout is fine sampling. phiin is the
  destination OPD. The weightings are obtained by interpolating from locin to
  locout, but the OPD are reversed computed */
void prop_nongrid_reverse(loc_t *locin,       
			  double* phiin,      
			  const loc_t *locout,
			  const double *ampout,
			  const double* phiout,
			  double alpha,
			  double displacex,
			  double displacey,
			  double scale){
    if(locin->dx<locout->dx) {
	error("This routine is designed for down sampling.\n");
    }
    loc_create_map_npad(locin,1);//will only do once and save in locin.
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    long iloc;
    int missing=0;
    const int wrapx1 = locin->map->nx;
    const int wrapy1 = locin->map->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;

#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->nx])(locin->map->p);
    //-1 because we count from 1 in the map.
    double *phiin0=phiin-1;
    for(iloc=0; iloc<locout->nloc; iloc++){
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;//skip points that has zero amplitude
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    missing++;
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}

#if ONLY_FULL == 1 //only proceed if all four points exist.
	iphi1=map[nplocy][nplocx];
	iphi2=map[nplocy][nplocx1];
	iphi3=map[nplocy1][nplocx];
	iphi4=map[nplocy1][nplocx1];
	if(iphi1 && iphi2 && iphi3 && iphi4){
	    phiin0[iphi1]+=alpha*(phiout[iloc]*(1.-dplocx)*(1.-dplocy));
	    phiin0[iphi2]+=alpha*(phiout[iloc]*(dplocx)*(1.-dplocy));
	    phiin0[iphi3]+=alpha*(phiout[iloc]*(1.-dplocx)*(dplocy));
	    phiin0[iphi4]+=alpha*(phiout[iloc]*(dplocx)*(dplocy));
	}else{
	    //if(!ampout || fabs(ampout[iloc])>EPS)
	    missing++;
	}
#else	
	if((iphi=map[nplocy][nplocx])) 
	    phiin0[iphi]+=alpha*(phiout[iloc]*(1.-dplocx)*(1.-dplocy));
	if((iphi=map[nplocy][nplocx1]))
	    phiin0[iphi]+=alpha*(phiout[iloc]*(dplocx)*(1.-dplocy));
	if((iphi=map[nplocy1][nplocx]))
	    phiin0[iphi]+=alpha*(phiout[iloc]*(1.-dplocx)*(dplocy));
	if((iphi=map[nplocy1][nplocx1]))
	    phiin0[iphi]+=alpha*(phiout[iloc]*(dplocx)*(dplocy));
#endif
    }
    if(missing>0){
	warning("%d points not covered by input screen\n", missing);
    }
}
