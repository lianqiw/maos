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

/*
  2009-09-19
  forked from mkgfull.c in matlab_utils
  2009-10-11
  seems to consume a lot more memory when ploc and xloc are different.
  the memory consumption is high, not release when mkg returns. 
  2009-10-12
  changed to create G' instead of G. the memory consumption is low now.
  2009-10-12
  make mkg.c a standalone routine. gather mex related into a separate file.
*/
#define _ISOC99_SOURCE
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#ifdef MATLAB_MEX_FILE
#include "../mex/mkgmex.h"
#else
#include "misc.h"
#include "mkg.h"
#include "bin.h"
#include "common.h"
#include "dsp.h"
#include "loc.h"
#include "dmat.h"
#include "cmat.h"
#endif

#ifndef DEBUG
#define DEBUG 0
#endif
#if DEBUG == 1
#define pline fprintf(stderr,"Line: %d\n", __LINE__);
#else
#define pline
#endif
#define BOUND0(A) A=(A<0.)?0.:A;
#define BOUND1(A) A=(A>1.)?1.:A;
#define SQ(A) (A)*(A)
#define CUBIC(A) (A)*(A)*(A)
#define CEIL(A) (fabs(A-round(A))<TOL?(int)round(A):iceil(A))
#define FLOOR(A) (fabs(A-round(A))<TOL?(int)round(A):(int)floor(A))
#define DBLNZ(A) (fabs(A)>TOL)

static const double TOL=1.e-14;
/**
   \file mkg.c
   Contains function that creates average gradient operator
*/
/*
  even if we do do_partial, we still not accounting for boundary pixels (some amp2 are missing)
  in this way, we assume the wavefront is flat in the boundary.
*/
#define DO_CALC								\
    if(do_partial){							\
	BOUND0(alpha1);							\
	BOUND1(alpha2);							\
	BOUND0(beta1);							\
	BOUND1(beta2);							\
        if(alpha1>1.-TOL||alpha2<-TOL||beta1>1.-TOL||beta2<-TOL){	\
	    valid=0;							\
	}else{								\
	    valid=1;							\
	}								\
	wtalpha[0]=(SQ(alpha2)-SQ(alpha1))*0.5-(CUBIC(alpha2)-CUBIC(alpha1))/3; \
	wtalpha[1]=(CUBIC(1-alpha1)-CUBIC(1-alpha2))/3;			\
	wtalpha[2]=(SQ(alpha2)-SQ(alpha1))*0.5;				\
	wtalpha[3]=(alpha2-alpha1)-(SQ(alpha2)-SQ(alpha1))*0.5;		\
	wtbeta[0]=(SQ(beta2)-SQ(beta1))*0.5-(CUBIC(beta2)-CUBIC(beta1))/3; \
	wtbeta[1]=(CUBIC(1-beta1)-CUBIC(1-beta2))/3;			\
	wtbeta[2]=(SQ(beta2)-SQ(beta1))*0.5;				\
	wtbeta[3]=(beta2-beta1)-(SQ(beta2)-SQ(beta1))*0.5;		\
    }									\
    if (valid &&							\
	DBLNZ(amp2[indy[0]][indx[0]]) &&				\
	DBLNZ(amp2[indy[0]][indx[1]]) &&				\
	DBLNZ(amp2[indy[1]][indx[0]]) &&				\
	DBLNZ(amp2[indy[1]][indx[1]]) ){				\
	weight[0]+=((wtbeta[0]*amp2[indy[0]][indx[0]]+			\
		     wtbeta[1]*amp2[indy[1]][indx[0]])*wtalpha[2]+	\
		    (wtbeta[0]*amp2[indy[0]][indx[1]]+			\
		     wtbeta[1]*amp2[indy[1]][indx[1]])*wtalpha[3])*dp1*signx; \
	weight[1]+=((wtbeta[2]*amp2[indy[0]][indx[0]]+			\
		     wtbeta[3]*amp2[indy[1]][indx[0]])*wtalpha[0]+	\
		    (wtbeta[2]*amp2[indy[0]][indx[1]]+			\
		     wtbeta[3]*amp2[indy[1]][indx[1]])*wtalpha[1])*dp1*signy; \
	ampsum+=wtalpha[3]*wtbeta[3]*amp2[1][1];			\
    }
/*
  dp1 above comes from the gradient.
*/

#define PMAP(ipix,limx,limy)			\
    {ipix=(limx>=0 && limx<ploc->map->nx	\
	   && limy>=0 && limy<ploc->map->ny)?	\
	    pmaps[limy][limx]:0;}
#define XMAP(iphi,nplocx,nplocy)			\
    {iphi=(nplocx>=0 && nplocx<xloc->map->nx		\
	   && nplocy>=0 && nplocy<xloc->map->ny)?	\
	    xmaps[nplocy][nplocx]:0;}
typedef struct{
    long i;
    double x;
}record_t;
#define ADDWT(A,B,C,D)				\
    XMAP(iphi,(A),(C));				\
    if(iphi) wtsum+=(B)*(D);

/**
   Returns the transpose of a average gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc with ploc as the intermediate pupil plane.
 */
dsp * mkgt(loc_t* xloc,     /**<the grid on which OPDs are defined*/
	   loc_t *ploc,     /**<the grid on the aperture plan*/
	   double *amp,     /**<the amplitude on ploc*/
	   loc_t *saloc,    /**<Coordinate of the subapertures*/
	   int saorc,       /**<0: saloc is the center of the subapertures. 1: saloc is the origin of the subapertures.*/
	   double scale,    /**<cone effect*/
	   double *displace,/**<displacement due to beam angle (2 vector). similar as accphi routines*/
	   int do_partial   /**<1: use points that are outside of by close to a subaperture. 0: do not use.*/
	   )
{
    /*
      
      Create an average gradient influence matrix from XLOC to SALOC, with
      amplitude maps defined on PLOC. PLOC and SALOC are defined on the same
      plan. The offset *displace* and scaling *scale* is between PLOC/SALOC and
      XLOC, when projected onto XLOC plane,
      SALOC->SALOC*scale+displace. PLOC->PLOC*scale+dispace. There are no offset
      or scaling between PLOC and SALOC

      saorc: SALOC is subaperture origin or center. 
      1: origin (lower left corner), 
      0: center.
    */
    double *ampcopy;
    double dsa2;
    double scx, scy;
    double weight[2];
    double dx1,dx2;
    double dp1,dp2,poffset[2];
    double displacex=0, displacey=0;
   
    int nsa,isa;
    int limx, limy;
    int same;
    double dx=xloc->dx;
    double dp=ploc->dx;
    double dsa=saloc->dx;

    nsa=saloc->nloc;
    dx1=1./dx;
    dp1=1./dp;
    dx2=scale*dx1;

    /*
      if do_partial is 1, will use points that are outside of the subaperture,
      but contribute.
    */
    /*notice, this map is in matlab index*/
    loc_create_map(ploc);
    poffset[0]=ploc->map->ox;
    poffset[1]=ploc->map->oy;
    if(xloc->locx && ploc->locx != xloc->locx){
	/*
	  If XLOC is present and different from PLOC, 
	  treat the misregistration
	  as between PLOC and XLOC. SALOC and PLOC has no misregistration.
	*/
	same=0;
	loc_create_map(xloc);
	poffset[0]*=-dp1;
	poffset[1]*=-dp1;

	/*displace between PLOC(SALOC) and XLOC. Only used if same==0;*/
	displacex=(displace[0]-xloc->map->ox)*dx1;
	displacey=(displace[1]-xloc->map->oy)*dx1;	
	dp2=dp1;
	if(fabs(scale-1)>1.e-10)
	    fprintf(stderr,"This is a bug in mkg to do three plane mkg when scale!=1\n");
    }else{
	/*
	  If XLOC is missing or XLOC and PLOC are the same, treat the
	  misregistration as between SALOC and PLOC
	*/
	same=1;
	poffset[0]=(displace[0]-ploc->map->ox)*dp1;
	poffset[1]=(displace[1]-ploc->map->oy)*dp1;
	dp2=dp1*scale;
    }

    
    /*Spacing of half subaperture on pupil plane.*/
    dsa2=dsa/2.*dp2;
    /*To store the weights.*/

    poffset[0]+=dsa2*saorc;/*offset to SALOC to make it subaperture center.*/
    poffset[1]+=dsa2*saorc;
    
    long (*pmaps)[ploc->map->nx]=(long (*)[ploc->map->nx])ploc->map->p;
    long (*xmaps)[xloc->map->nx]=(long (*)[xloc->map->nx])xloc->map->p;
    double amp2[3][3];
    double wtfull[4]={0.5/3, 1./3, 0.5, 0.5};
    double *wtalpha, *wtbeta;
    double signx, signy;
    double ampsum;
    double limx1, limx2, limy1, limy2;
    double alpha1, alpha2, beta1, beta2;
    double dplocx,plocx,dplocy,plocy;
    double wt0,wtsum;
    int nplocx,nplocy;
    int ipix,ipix2,itmp,jtmp;
    int indx[2],indy[2];
    int valid=1;
    int iw;
    int iphi;
    if(!same && amp){
	/*
	  Copy and modify amplitude map to fix boundry
	  issues when propage between two planes.  This
	  generally applies to the case that XLOC is smaller
	  than PLOC.
	*/
	ampcopy=malloc(sizeof(double)*ploc->nloc);
	memcpy(ampcopy,amp,sizeof(double)*ploc->nloc);
	for(ipix=0; ipix<ploc->nloc; ipix++){
	    plocx=ploc->locx[ipix]*dx2+displacex;
	    plocy=ploc->locy[ipix]*dx2+displacey;
	    SPLIT(plocx,dplocx,nplocx);
	    SPLIT(plocy,dplocy,nplocy);
	    /*Find out whether this PLOC point has XLOC points associated,
	      if not, will scale the amplitude of this PLOC point.*/
	    wtsum=0;
	    ADDWT(nplocx,1-dplocx,nplocy,1-dplocy);
	    ADDWT(1+nplocx,dplocx,nplocy,1-dplocy);
	    ADDWT(nplocx,1-dplocx,1+nplocy,dplocy);
	    ADDWT(1+nplocx,dplocx,1+nplocy,dplocy);
	    ampcopy[ipix]*=wtsum;
	}
    }else{
	ampcopy=amp;
    }
    if(do_partial){
	wtalpha=(double*)alloca(sizeof(double)*4);
	wtbeta =(double*)alloca(sizeof(double)*4);
    }else{
	wtalpha=(double*)wtfull;
	wtbeta =(double*)wtfull;
    }
    /*if no amplitude weighting, set all to one.*/
    if(!amp){
	for(jtmp=0; jtmp<3; jtmp++){
	    for(itmp=0;itmp<3;itmp++){
		amp2[jtmp][itmp]=1;
	    }
	}
    }
    dsp *GS0t[2];
    long nzmax[2];
    spint *restrict pp[2];
    spint *restrict pi[2];
    double  *restrict px[2];
    long count[2];
    for(iw=0; iw<2; iw++){
	nzmax[iw]=nsa*(dsa/dx)*16;//four row in each side.
	GS0t[iw]=spnew(xloc->nloc, nsa, nzmax[iw]);
	pp[iw]=GS0t[iw]->p; 
	pi[iw]=GS0t[iw]->i; 
	px[iw]=GS0t[iw]->x;
	count[iw]=0;
    }
    for(isa=0; isa<nsa; isa++){
	for(iw=0; iw<2; iw++){
	    pp[iw][isa]=count[iw];
	}
	/*center of subaperture when mapped onto PLOC MAP.*/
	scx=saloc->locx[isa]*dp2+poffset[0];
	scy=saloc->locy[isa]*dp2+poffset[1];
	ampsum=0; /*Accumulated amplitude integral in this subaperture*/
	/*Boundary of PLOC points covered by this subaperture.*/
	limx1=(scx-dsa2);
	limx2=(scx+dsa2);
	limy1=(scy-dsa2);
	limy2=(scy+dsa2);
	/*
	  was the following when only account for fully illuminated subaps.
	  for (limy=FLOOR(limy1); limy<CEIL(limy2)+1; limy++){
	  for (limx=FLOOR(limx1); limx<CEIL(limx2)+1; limx++){
	  The new limits are used to include points that are outside of 
	  but adjacent to the subaperture
	*/
	    
	for (limy=CEIL(limy1)-1; limy<FLOOR(limy2)+2; limy++){
	    for (limx=CEIL(limx1)-1; limx<FLOOR(limx2)+2; limx++){
		PMAP(ipix,limx,limy);

		if (ipix--){
		    weight[0]=0;
		    weight[1]=0;
		    /*if amplitude map is defined, fetch the 3x3 map surrounding this point.*/
		    if(amp){
			for(jtmp=0; jtmp<3; jtmp++){
			    for(itmp=0; itmp<3; itmp++){
				PMAP(ipix2,limx+itmp-1,limy+jtmp-1);
				if(ipix2--){
				    amp2[jtmp][itmp]=ampcopy[ipix2];
				}else{
				    amp2[jtmp][itmp]=0;
				}
			    }
			}
		    }
		    /*
		      Integrate over the four quadrants surround the point
		      ipix, inside the subaperture
		      limx, limy are the coordinate of the point in PLOC MAP.
		      [limx1,limx2] [limy1,limy2] bound the subaperture.

		      Notice that the coordinates are positive to the left/top for integral.
		    */
		    /*left use 0 1*/
		    alpha1=limx-limx2;/*The lower bound to do integral (right)*/
		    alpha2=limx-limx1;/*The upper bound to do integral (left)*/
		    if(do_partial||(alpha1<TOL && alpha2>1.-TOL)){
			indx[0]=0; indx[1]=1;
			signx=1.;
			/*bottom use 0 1*/
			beta1=limy-limy2;
			beta2=limy-limy1;
			if(do_partial||(beta1<TOL && beta2>1.-TOL)){
			    indy[0]=0; indy[1]=1;
			    signy=1.;
			    DO_CALC;
			}
			/*top use 2 1 */
			beta1=limy1-limy;
			beta2=limy2-limy;
			if(do_partial||(beta1<TOL && beta2>1.-TOL)){
			    indy[0]=2; indy[1]=1;
			    signy=-1.;
			    DO_CALC;
			}
		    }
		    /*right  use 2 1*/
		    alpha1=limx1-limx;
		    alpha2=limx2-limx;
		    if(do_partial||(alpha1<TOL && alpha2>1.-TOL)){
			indx[0]=2; indx[1]=1;
			signx=-1.;
			/*bottom use 0 1*/
			beta1=limy-limy2;
			beta2=limy-limy1;
			if(do_partial||(beta1<TOL && beta2>1.-TOL)){
			    indy[0]=0; indy[1]=1;
			    signy=1.;
			    DO_CALC;
			}
			/*top use 2 1 */
			beta1=limy1-limy;
			beta2=limy2-limy;
			if(do_partial||(beta1<TOL && beta2>1.-TOL)){
			    indy[0]=2; indy[1]=1;
			    signy=-1.;
			    DO_CALC;
			}
		    }
		    for(iw=0; iw<2; iw++){
			if(count[iw]+4>nzmax[iw]){
			    nzmax[iw]=nzmax[iw]*2;
			    spsetnzmax(GS0t[iw], nzmax[iw]);
			    pp[iw]=GS0t[iw]->p; 
			    pi[iw]=GS0t[iw]->i; 
			    px[iw]=GS0t[iw]->x;
			}
			if(fabs(weight[iw])>TOL){
			    if(same){
				/*XLOC and PLOC are the same, 
				  and no misregistration.*/
				pi[iw][count[iw]]=ipix;
				px[iw][count[iw]++]=weight[iw];
			    }else{
				/*
				  XLOC and PLOC are not the same, find the
				  points on XLOC that couples in to ipix.
				*/

				plocx=ploc->locx[ipix]*dx2+displacex;
				plocy=ploc->locy[ipix]*dx2+displacey;
				SPLIT(plocx,dplocx,nplocx);
				SPLIT(plocy,dplocy,nplocy);
				//info2("weight[%d]=%g\n",iw,weight[iw]);
				//info2("(%g, %g, %d), (%g, %g, %d)\n",
				//    plocx,dplocx,nplocx,
				//   plocy,dplocy,nplocy);
				wtsum=0;
				ADDWT(nplocx,1-dplocx,nplocy,1-dplocy);
				ADDWT(1+nplocx,dplocx,nplocy,1-dplocy);
				ADDWT(nplocx,1-dplocx,1+nplocy,dplocy);
				ADDWT(1+nplocx,dplocx,1+nplocy,dplocy);
				//info2("wtsum=%g\n",wtsum);
				wtsum=1./wtsum;
				/*We are computing gradients on the PLOC.
				  This gradient is scaled by 1/scale in XLOC*/
				if(!same)
				    wtsum/=scale;
#define INTERP(A,B,C,D)							\
    XMAP(iphi,A,C);							\
    if(iphi--){								\
	/*Check whether we already have the weight for this point.*/	\
	if(fabs(wt0=weight[iw]*(B)*(D)*wtsum)>TOL){			\
	    /*info2("(%d,%d)=%g)\n",isa+nsa*iw,iphi,wt0);*/		\
	    for(itmp=pp[iw][isa]; itmp<count[iw]; itmp++){		\
		if(pi[iw][itmp]==iphi){					\
		    px[iw][itmp]+=wt0;					\
		    break;						\
		}							\
	    }								\
	    if(itmp == count[iw]){/*not found*/				\
		pi[iw][count[iw]]=iphi;					\
		px[iw][count[iw]++]=wt0;				\
	    }								\
	}								\
    }
				    
				INTERP(nplocx,1-dplocx,nplocy,1-dplocy);
				INTERP(1+nplocx,dplocx,nplocy,1-dplocy);
				INTERP(nplocx,1-dplocx,1+nplocy,dplocy);
				INTERP(1+nplocx,dplocx,1+nplocy,dplocy);
#undef INTERP
			    }
			}
		    }
		}/*if ipix*/
	    }/*limx*/
	}/*limy*/
	//info2("ampsum=%g\n",ampsum);
	if(ampsum>TOL){
	    ampsum=1./ampsum;
	}else{
	    ampsum=0;
	}
	for(iw=0; iw<2; iw++){
	    for(int ii=pp[iw][isa]; ii<count[iw]; ii++){
		px[iw][ii]*=ampsum;
	    }
	}
		//raise(SIGSEGV);
    }/*isa*/
    for(iw=0; iw<2; iw++){
	if(count[iw]>nzmax[iw]){
	    error("Over flow\n");//should never happen
	}
	pp[iw][nsa]=count[iw];
	spsetnzmax(GS0t[iw],count[iw]);
    }
    if(ploc->map){
	free(ploc->map->p);
	free(ploc->map);
	ploc->map=NULL;
    }
    if(xloc->map){
	free(xloc->map->p);
	free(xloc->map);
	xloc->map=NULL;
    }
    if(ampcopy!=amp)
	free(ampcopy);
    //concatenate x and y gradient operators together.
    dsp *GS0=spcat(GS0t[0], GS0t[1],1);
    spfree(GS0t[0]); GS0t[0]=NULL;
    spfree(GS0t[1]); GS0t[1]=NULL;
    loc_free_map(ploc);
    loc_free_map(xloc);
    return GS0;
}
#ifndef MATLAB_MEX_FILE
/**
   Returns the transpose of mkgt()
 */
dsp *mkg(loc_t* xloc, loc_t *ploc, double *amp, loc_t *saloc, 
	 int saorc, double scale, double *displace, int do_partial){
    dsp *GS0T=mkgt(xloc, ploc, amp, saloc, 
		   saorc, scale, displace, do_partial);
    dsp *GS0=sptrans(GS0T);
    spfree(GS0T);
    return GS0;
}
#endif

