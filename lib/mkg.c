/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
//#define _ISOC99_SOURCE

#include "../math/mathdef.h"
#include "mkg.h"

#ifndef DEBUG
#define DEBUG 0
#endif
#if DEBUG == 1
#define pline info("Line: %d\n", __LINE__);
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

static const real TOL=10*EPS;
#define OFF 1 /*0: Only use points inside subaperture. introduced on 2011-08-28. 1: original approach. */

/*
  even if we do do_partial, we still not accounting for boundary pixels (some amp2 are missing)
  in this way, we assume the wavefront is flat in the boundary.
*/
#define DO_CALC\
    if(do_partial){\
	BOUND0(alpha1);\
	BOUND1(alpha2);\
	BOUND0(beta1);\
	BOUND1(beta2);\
        if(alpha1>1.-TOL||alpha2<-TOL||beta1>1.-TOL||beta2<-TOL){\
	    valid=0;\
	}else{\
	    valid=1;\
	}\
	wtalpha[0]=(SQ(alpha2)-SQ(alpha1))*0.5-(CUBIC(alpha2)-CUBIC(alpha1))/3; \
	wtalpha[1]=(CUBIC(1-alpha1)-CUBIC(1-alpha2))/3;\
	wtalpha[2]=(SQ(alpha2)-SQ(alpha1))*0.5;\
	wtalpha[3]=(alpha2-alpha1)-(SQ(alpha2)-SQ(alpha1))*0.5;\
	wtbeta[0]=(SQ(beta2)-SQ(beta1))*0.5-(CUBIC(beta2)-CUBIC(beta1))/3; \
	wtbeta[1]=(CUBIC(1-beta1)-CUBIC(1-beta2))/3;\
	wtbeta[2]=(SQ(beta2)-SQ(beta1))*0.5;\
	wtbeta[3]=(beta2-beta1)-(SQ(beta2)-SQ(beta1))*0.5;\
    }\
    if (valid &&\
	DBLNZ(amp2[indy[0]][indx[0]]) &&\
	DBLNZ(amp2[indy[0]][indx[1]]) &&\
	DBLNZ(amp2[indy[1]][indx[0]]) &&\
	DBLNZ(amp2[indy[1]][indx[1]]) ){\
	weight[0]+=((wtbeta[0]*amp2[indy[0]][indx[0]]+\
		     wtbeta[1]*amp2[indy[1]][indx[0]])*wtalpha[2]+\
		    (wtbeta[0]*amp2[indy[0]][indx[1]]+\
		     wtbeta[1]*amp2[indy[1]][indx[1]])*wtalpha[3])*dp1*signx; \
	weight[1]+=((wtbeta[2]*amp2[indy[0]][indx[0]]+\
		     wtbeta[3]*amp2[indy[1]][indx[0]])*wtalpha[0]+\
		    (wtbeta[2]*amp2[indy[0]][indx[1]]+\
		     wtbeta[3]*amp2[indy[1]][indx[1]])*wtalpha[1])*dp1*signy; \
	ampsum+=wtalpha[3]*wtbeta[3]*amp2[1][1];\
    }
/*
  dp1 above comes from the gradient.
*/

typedef struct{
	long i;
	real x;
}record_t;
#define ADDWT(A,B,C,D)\
    iphi=loc_map_get(xloc->map, (A),(C));\
    if(iphi>0) wtsum+=(B)*(D);

/**
   Returns the transpose of a average gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc with ploc as the intermediate pupil plane.
 */
dsp* mkgt(loc_t* xloc,/**<the grid on which OPDs are defined*/
	loc_t* ploc,     /**<the grid on the aperture plane*/
	dmat* pamp,      /**<the amplitude on ploc*/
	loc_t* saloc,    /**<Lower left origin of the subapertures*/
	dmat *saa,	     /**<subaperture amplitude, normalized to max at 1.*/
	real saat,		 /**<saa threshold to enable gradient operation. */
	real scale,      /**<cone effect*/
	real dispx,      /**<displacement due to beam angle (2 vector). similar as accphi routines*/
	real dispy,      /**<displacement due to beam angle (2 vector). similar as accphi routines*/
	int do_partial   /**<1: use points that are outside of but close to a subaperture. 0: do not use.*/
){
/*

  Create an average gradient influence matrix from XLOC to SALOC, with
  amplitude maps defined on PLOC. PLOC and SALOC are defined on the same
  plan. The offset *displace* and scaling *scale* is between PLOC/SALOC and
  XLOC, when projected onto XLOC plane,
  SALOC->SALOC*scale+displace. PLOC->PLOC*scale+dispace. There are no offset
  or scaling between PLOC and SALOC

  1: origin (lower left corner),
  0: center.
*/
	real* ampcopy;
	real dsa2;
	real scx, scy;
	real weight[2];
	real dx1, dx2;
	real dp1, dp2, poffset[2];
	int nsa, isa;
	int limx, limy;
	int same;
	real dx=xloc->dx;
	real dp=ploc->dx;
	real dsa=saloc->dx;

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
	if(xloc->locx&&ploc->locx!=xloc->locx){
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
		dispx=(dispx-xloc->map->ox)*dx1;
		dispy=(dispy-xloc->map->oy)*dx1;
		dp2=dp1;
		if(fabs(scale-1)>1.e-10)
			info("This is a bug in mkg to do three plane mkg when scale!=1\n");
	} else{
	/*
	  If XLOC is missing or XLOC and PLOC are the same, treat the
	  misregistration as between SALOC and PLOC
	*/
		same=1;
		poffset[0]=(dispx-ploc->map->ox)*dp1;
		poffset[1]=(dispy-ploc->map->oy)*dp1;
		dispx=0;
		dispy=0;
		dp2=dp1*scale;
	}


	/*Spacing of half subaperture on pupil plane.*/
	dsa2=dsa/2.*dp2;
	/*To store the weights.*/
	real amp_thres=dsa2*dsa2*4*0.01;/*1% area is the lower threshold to use subaperture*/
	poffset[0]+=dsa2;/*offset to SALOC to make it subaperture center.*/
	poffset[1]+=dsa2;
	real amp2[3][3];
	real wtalpha[4]={0.5/3, 1./3, 0.5, 0.5};//if do_partial is set, the values will be updated
	real wtbeta[4]={0.5/3, 1./3, 0.5, 0.5};
	real signx, signy;
	real ampsum;
	real limx1, limx2, limy1, limy2;
	real alpha1, alpha2, beta1, beta2;
	real dplocx, plocx, dplocy, plocy;
	real wt0, wtsum;
	int nplocx, nplocy;
	int ipix, ipix2, itmp, jtmp;
	int indx[2], indy[2];
	int valid=1;
	int iw;
	int iphi;
	real* amp=pamp?P(pamp):0;
	if(!same&&amp){
		/*
		Copy and modify amplitude map to fix boundry
		issues when propagate between two planes.  This
		generally applies to the case that XLOC is smaller
		than PLOC.
		*/
		ampcopy=mymalloc(ploc->nloc, real);
		memcpy(ampcopy, amp, sizeof(real)*ploc->nloc);
		for(ipix=0; ipix<ploc->nloc; ipix++){
			plocx=ploc->locx[ipix]*dx2+dispx;
			plocy=ploc->locy[ipix]*dx2+dispy;
			SPLIT(plocx, dplocx, nplocx);
			SPLIT(plocy, dplocy, nplocy);
			/*Find out whether this PLOC point has XLOC points associated,
			  if not, will scale the amplitude of this PLOC point.*/
			wtsum=0;
			ADDWT(nplocx, 1-dplocx, nplocy, 1-dplocy);
			ADDWT(1+nplocx, dplocx, nplocy, 1-dplocy);
			ADDWT(nplocx, 1-dplocx, 1+nplocy, dplocy);
			ADDWT(1+nplocx, dplocx, 1+nplocy, dplocy);
			ampcopy[ipix]*=wtsum;
		}
	} else{
		ampcopy=amp;
	}
	
	/*if no amplitude weighting, set all to one.*/
	if(!amp){
		for(jtmp=0; jtmp<3; jtmp++){
			for(itmp=0;itmp<3;itmp++){
				amp2[jtmp][itmp]=1;
			}
		}
	}
	dsp* GS0t[2];
	long nzmax[2];
	spint* restrict pp[2];
	spint* restrict pi[2];
	real* restrict px[2];
	long count[2];
	for(iw=0; iw<2; iw++){
		nzmax[iw]=nsa*pow(dsa/dx+3,2);
		GS0t[iw]=dspnew(xloc->nloc, nsa, nzmax[iw]);
		pp[iw]=GS0t[iw]->pp;
		pi[iw]=GS0t[iw]->pi;
		px[iw]=GS0t[iw]->px;
		count[iw]=0;
	}
	for(isa=0; isa<nsa; isa++){
		for(iw=0; iw<2; iw++){
			pp[iw][isa]=count[iw];
		}
		if(saa && saat && P(saa,isa)<saat) continue;
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
		for(limy=CEIL(limy1)-OFF; limy<FLOOR(limy2)+OFF+1; limy++){
			for(limx=CEIL(limx1)-OFF; limx<FLOOR(limx2)+OFF+1; limx++){
				ipix=loc_map_get(ploc->map, limx, limy);
				if(ipix>0){
					ipix--;
					weight[0]=0;
					weight[1]=0;
					/*if amplitude map is defined, fetch the 3x3 map surrounding this point.*/
					if(amp){
						for(jtmp=0; jtmp<3; jtmp++){
							for(itmp=0; itmp<3; itmp++){
								ipix2=loc_map_get(ploc->map, limx+itmp-1, limy+jtmp-1);
								if(ipix2>0){
									amp2[jtmp][itmp]=ampcopy[ipix2-1];
								} else{
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
					if(do_partial||(alpha1<TOL&&alpha2>1.-TOL)){
						indx[0]=0; indx[1]=1;
						signx=1.;
						/*bottom use 0 1*/
						beta1=limy-limy2;
						beta2=limy-limy1;
						if(do_partial||(beta1<TOL&&beta2>1.-TOL)){
							indy[0]=0; indy[1]=1;
							signy=1.;
							DO_CALC;
						}
						/*top use 2 1 */
						beta1=limy1-limy;
						beta2=limy2-limy;
						if(do_partial||(beta1<TOL&&beta2>1.-TOL)){
							indy[0]=2; indy[1]=1;
							signy=-1.;
							DO_CALC;
						}
					}
					/*right  use 2 1*/
					alpha1=limx1-limx;
					alpha2=limx2-limx;
					if(do_partial||(alpha1<TOL&&alpha2>1.-TOL)){
						indx[0]=2; indx[1]=1;
						signx=-1.;
						/*bottom use 0 1*/
						beta1=limy-limy2;
						beta2=limy-limy1;
						if(do_partial||(beta1<TOL&&beta2>1.-TOL)){
							indy[0]=0; indy[1]=1;
							signy=1.;
							DO_CALC;
						}
						/*top use 2 1 */
						beta1=limy1-limy;
						beta2=limy2-limy;
						if(do_partial||(beta1<TOL&&beta2>1.-TOL)){
							indy[0]=2; indy[1]=1;
							signy=-1.;
							DO_CALC;
						}
					}
					for(iw=0; iw<2; iw++){
						if(count[iw]+4>nzmax[iw]){
							dbg("Increasing size of GP.\n");
							nzmax[iw]=nzmax[iw]*2;
							dspsetnzmax(GS0t[iw], nzmax[iw]);
							pp[iw]=GS0t[iw]->pp;
							pi[iw]=GS0t[iw]->pi;
							px[iw]=GS0t[iw]->px;
						}
						if(fabs(weight[iw])>TOL){
							if(same){
							/*XLOC and PLOC are the same,
							  and no misregistration.*/
								pi[iw][count[iw]]=ipix;
								px[iw][count[iw]++]=weight[iw];
							} else{
							/*
							  XLOC and PLOC are not the same, find the
							  points on XLOC that couples in to ipix.
							*/

								plocx=ploc->locx[ipix]*dx2+dispx;
								plocy=ploc->locy[ipix]*dx2+dispy;
								SPLIT(plocx, dplocx, nplocx);
								SPLIT(plocy, dplocy, nplocy);
								/*info("weight[%d]=%g\n",iw,weight[iw]); */
								/*info("(%g, %g, %d), (%g, %g, %d)\n", */
								/*    plocx,dplocx,nplocx, */
								/*   plocy,dplocy,nplocy); */
								wtsum=0;
								ADDWT(nplocx, 1-dplocx, nplocy, 1-dplocy);
								ADDWT(1+nplocx, dplocx, nplocy, 1-dplocy);
								ADDWT(nplocx, 1-dplocx, 1+nplocy, dplocy);
								ADDWT(1+nplocx, dplocx, 1+nplocy, dplocy);
								/*info("wtsum=%g\n",wtsum); */
								wtsum=1./wtsum;
								/*We are computing gradients on the PLOC.
								  This gradient is scaled by 1/scale in XLOC*/
								if(!same)
									wtsum/=scale;
								#define INTERP(A,B,C,D)\
								iphi=loc_map_get(xloc->map, A,C);\
								if(iphi>0){\
									iphi--;\
									/*Check whether we already have the weight for this point.*/\
									if(fabs(wt0=weight[iw]*(B)*(D)*wtsum)>TOL){\
										/*info("(%d,%d)=%g)\n",isa+nsa*iw,iphi,wt0);*/\
										for(itmp=pp[iw][isa]; itmp<count[iw]; itmp++){\
											if(pi[iw][itmp]==iphi){\
												px[iw][itmp]+=wt0;\
												break;\
											}\
										}\
										if(itmp == count[iw]){/*not found*/\
											pi[iw][count[iw]]=iphi;\
											px[iw][count[iw]++]=wt0;\
										}\
									}\
								}

								INTERP(nplocx, 1-dplocx, nplocy, 1-dplocy);
								INTERP(1+nplocx, dplocx, nplocy, 1-dplocy);
								INTERP(nplocx, 1-dplocx, 1+nplocy, dplocy);
								INTERP(1+nplocx, dplocx, 1+nplocy, dplocy);
								#undef INTERP
							}
						}
					}
				}/*if ipix*/
			}/*limx*/
		}/*limy*/
		/*info("ampsum=%g\n",ampsum); */
		if(ampsum>amp_thres){//2017-09-19: was TOL. Changed to amp_thres to filter weak subapertures
			ampsum=1./ampsum;
		} else{
			ampsum=0;
		}
		for(iw=0; iw<2; iw++){
			int jj=pp[iw][isa];
			for(int ii=pp[iw][isa]; ii<count[iw]; ii++){
				px[iw][ii]*=ampsum;
				if(fabs(px[iw][ii])>TOL){//<TOL after accumulation, remove.
					if(jj!=ii){
						px[iw][jj]=px[iw][ii];
						pi[iw][jj]=pi[iw][ii];
					}
					jj++;
				}
			}
			count[iw]=jj;
		}
	}/*isa*/
	for(iw=0; iw<2; iw++){
		if(count[iw]>nzmax[iw]){
			error("Over flow\n");/*should never happen */
		}
		if(count[iw]==0){
			error("count[%d]=0\n", iw);
		}
		pp[iw][nsa]=count[iw];
		dspsetnzmax(GS0t[iw], count[iw]);
	}
	if(ampcopy!=amp)
		free(ampcopy);
		/*concatenate x and y gradient operators together. */
	dsp* GS0=dspcat(GS0t[0], GS0t[1], 1);
	dspfree(GS0t[0]); GS0t[0]=NULL;
	dspfree(GS0t[1]); GS0t[1]=NULL;
	/*loc_free_map(ploc);
	  loc_free_map(xloc);*/
	dspdroptol(GS0, 1e-14);/*drop small values. */
	return GS0;
}
#ifndef MATLAB_MEX_FILE
/**
   Returns the transpose of mkgt()
 */
dsp* mkg(loc_t* xloc, loc_t* ploc, dmat* amp, loc_t* saloc, dmat *saa, real saat,
	real scale, real dispx, real dispy, int do_partial){
	dsp* GS0T=mkgt(xloc, ploc, amp, saloc, saa, saat, scale, dispx, dispy, do_partial);
	dsp* GS0=dsptrans(GS0T);
	dspfree(GS0T);
	return GS0;
}
#endif
/**
 * @brief Compute the NGS modal to t/t or t/t/f gradient interaction matrix. Use for sky coverage code for adhoc split tomography
 * 
 * @param sacent 	Amplitude weighted center coordinate of each subaperture (nsa*2)
 * @param thetax 	angular coordinate of NGS
 * @param thetay 	angular coordinate of NGS
 * @param hc 		upper DM conjugation range (for MCAO)
 * @param hs 		LGS range
 * @param indps 	Plate scale mode index (0 or 2) (for MCAO)
 * @param indastig  Astigmatism mode index (0 or 2) (for LTAO)
 * @param indfocus  Focus mode index (0 or indps+3 or indastig+2
 * @param ahstfocus =1: no focus mode in first PS mode in science OPD.
 * @return dmat 
 */
dmat* mkg_ngs(dmat *sacent, real thetax, real thetay, real hc, real hs, int indps, int indastig, int indfocus, int ahstfocus){
	if(!sacent) error("sacent must be specified as a nsa*2 matrix\n");
	const int nsa=NX(sacent);
	const int nmod=2+(indps?3:0)+(indastig?2:0)+(indfocus?1:0);
	const real scale=pow(1.-hc/hs, -2);
	const real scale1=1.-scale;
	dmat *pg=dnew(nsa*2, nmod);
	if(indps+3>nmod || indastig+2>nmod || indfocus+1>nmod){
		error("indps=%d, indastig=%d or indfocus=%d is incorrect\n", indps, indastig, indfocus);
	}
	for(long isa=0; isa<nsa; isa++){
		const real xm=P(sacent,isa,0);/*dot of x with amp. */
		const real ym=P(sacent,isa,1);

		P(pg, isa, 0)=1.;//tip
		P(pg, isa+nsa, 1)=1.;//tilt
		if(indps){
			//if ahstfocus==1, first PS mode has no focus component in science.
			P(pg, isa, indps)=(ahstfocus?0:scale1*2*xm)-2*thetax*hc*scale;
			P(pg, isa+nsa, indps)=(ahstfocus?0:scale1*2*ym)-2*thetay*hc*scale;
			
			P(pg, isa, indps+1)=(scale1*2*xm-2*thetax*hc*scale);
			P(pg, isa+nsa, indps+1)=(-scale1*2*ym+2*thetay*hc*scale);
			P(pg, isa, indps+2)=(scale1*ym-thetay*hc*scale);
			P(pg, isa+nsa, indps+2)=(scale1*xm-thetax*hc*scale);
		}
		if(indastig){
			P(pg, isa, indastig)=(2*xm);//d(x^2-y^2)/dx
			P(pg, isa+nsa, indastig)=(-2*ym);
			P(pg, isa, indastig+1)=(ym);//d(xy)/dx
			P(pg, isa+nsa, indastig+1)=(xm);
		}
		if(indfocus){
			P(pg, isa, indfocus)=xm*2;
			P(pg, isa+nsa, indfocus)=ym*2;
		}
	}
	return pg;
}

