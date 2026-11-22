/*
  Copyright 2009-2026 Lianqi Wang

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
#include "shwfs.h"
#include "mkg.h"
#include "libmisc.h"

/**
   Create the lower left corner of the subapertures of a Shack-Hartmann WFS on a
   square grid centered on the aperture. A subaperture is kept if it overlaps
   with the annular aperture enlarged by the misregistration. Subapertures with
   insufficient illumination are removed later with shwfs_sa_isolated() and the
   area threshold.

   @param[in] order Number of subapertures across the aperture.
   @param[in] dsa Subaperture size in meter.
   @param[in] din Inner diameter of the aperture (central obscuration) in meter.
   @param[in] misreg Maximum misregistration of the WFS (lenslet array) in meter.
   @param[in] misreg_aper Misregistration of the telescope aperture in meter.
   @return The subaperture grid with spacing dsa. Newly allocated; the caller frees it.
*/
loc_t* shwfs_saloc(int order, real dsa, real din, real misreg, real misreg_aper){
	/*Offset of the coordinate of the center most subaperture from the center. */
	real offsetx=((order&1)?-0.5:0);
	real offsety=((order&1)?-0.5:0);
	//make sure orderpad has the same evenness as order. orderpad is larger than order to handle misregistration.
	const int orderpad=order+(order>2?(1+misreg/dsa)*2:0);//2*MAX(fabs(offsetx), fabs(offsety)); //to avoid missing subapertures when there is saoffx/saoffy
	dbg("saloc: offsetx=%g, offsety=%g, order=%d, orderpad=%d\n", offsetx, offsety, order, orderpad);
	/*Keep a subaperture if it overlaps with the annular aperture, enlarged
	  by WFS and aperture misregistration (in unit of dsa). Subapertures
	  with insufficient illumination are removed later by saat.*/
	const real rmargin=(misreg+misreg_aper)/dsa;
	/*r2max: Maximum distance^2 from the center to keep a subaperture */
	const real r2max=pow(order*0.5+rmargin, 2);
	/*r2min: subapertures entirely within this distance^2 are obscured. */
	const real rin=din/dsa/2-rmargin;
	const real r2min=(dsa<din&&rin>0)?pow(rin, 2):-1;
	/*the lower left *grid* coordinate of the subaperture */
	loc_t* saloc=locnew(orderpad*orderpad, dsa, dsa);
	int count=0;
	/*Collect all the subapertures that are within the allowed radius*/
	for(int j=-orderpad/2; j<=(orderpad-1)/2; j++){
		for(int i=-orderpad/2; i<=(orderpad-1)/2; i++){
			//Normalized coordinate in uniq of sa size
			real xc=((real)i+offsetx);
			real yc=((real)j+offsety);
			//Minimum distance^2 from the center to the subaperture
			const real dxmin=xc>0?xc:(xc+1<0?-(xc+1):0);
			const real dymin=yc>0?yc:(yc+1<0?-(yc+1):0);
			const real rmin2=dxmin*dxmin+dymin*dymin;
			//Maximum distance^2 from the center to the subaperture (farthest corner)
			const real dxmax=MAX(fabs(xc), fabs(xc+1));
			const real dymax=MAX(fabs(yc), fabs(yc+1));
			const real rmax2=dxmax*dxmax+dymax*dymax;
			//Keep it if it is partially inside the outer edge and not entirely inside the inner obscuration.
			if(rmin2<r2max && rmax2>r2min){
				saloc->locx[count]=xc*dsa;
				saloc->locy[count]=yc*dsa;
				count++;
			}
		}
	}
	locresize(saloc, count);//also moves locy. Do not set nloc directly.
	return saloc;
}
/**
   Create the points in each subaperture. The number of points along each side
   of a subaperture is made even (for FFT) and the sampling is adjusted so that
   it divides dsa. The points are offset by half of the sampling from the lower
   left corner so that they are inside the subaperture.

   @param[in] saloc Lower left corner of each subaperture.
   @param[in] dsa Subaperture size in meter.
   @param[in] dx Requested sampling of the points in meter.
   @return The points. Its nxsa and dx are the number of points along each side and the adjusted sampling. Newly allocated; the caller frees it.
*/
pts_t* shwfs_pts(const loc_t* saloc, real dsa, real dx){
	/*Number of OPD pixels in 1 direction in each
	  subaperture. Make it even to do fft.*/
	const int nxsa=2*(int)round(0.5*dsa/dx);
	const real dx2=dsa/nxsa;/*adjust dx. */
	const real dxoffset=dx2*0.5;//Always keep points inside subaperture for simulation.
	if(fabs(dx-dx2)>EPS)
		info("Adjusting dx from %g to %g\n", dx, dx2);
	if(fabs(dsa-nxsa*dx2)>EPS){
		warning("nx=%d,dsa=%f,dx=%f not agree\n", nxsa, dsa, dx2);
	}
	dbg("There are %dx%d points in each subaperture of %gm.\n", nxsa, nxsa, dsa);
	pts_t* pts=ptsnew(saloc->nloc, dsa, dsa, nxsa, nxsa, dx2, dx2);
	for(int isa=0; isa<saloc->nloc; isa++){
		pts->origx[isa]=saloc->locx[isa]+dxoffset;
		pts->origy[isa]=saloc->locy[isa]+dxoffset;
	}
	return pts;
}
/**
   Create the amplitude map of a subaperture that simulates the lenslet fill
   factor: the points along the edge are turned off, and the outermost points
   that are kept are partially transmitting, so that the sum is
   nxsa*nxsa*safill2d.

   @param[in] nxsa Number of points along each side of the subaperture.
   @param[in] safill2d Fill factor of the lenslet in 2-d, between 0 and 1.
   @return nxsa*nxsa amplitude map. Newly allocated; the caller frees it.
*/
dmat* shwfs_safill(int nxsa, real safill2d){
	int nedge=1;
	while((nxsa-2*nedge)*(nxsa-2*nedge)>nxsa*nxsa*safill2d){
		nedge++;
	}
	dmat* ampi=dnew(nxsa, nxsa);
	real alpha=(nxsa*nxsa*safill2d-(nxsa-2*nedge)*(nxsa-2*nedge))
		/((nxsa-2*(nedge-1))*(nxsa-2*(nedge-1))-(nxsa-2*nedge)*(nxsa-2*nedge));
	for(int iy=nedge-1; iy<nxsa-nedge+1; iy++){
		for(int ix=nedge-1; ix<nxsa-nedge+1; ix++){
			P(ampi, ix, iy)=1;
			if(ix==nedge-1||ix==nxsa-nedge||iy==nedge-1||iy==nxsa-nedge){
				P(ampi, ix, iy)=alpha;
			}
		}
	}
	return ampi;
}
/**
   Convert amplitude map of subapertures to normalized subaperture area.

   @param[in] amp Amplitude map of all subapertures, ordered by subaperture.
   @param[in] npsa Number of points per subaperture.
   @return Subaperture area normalized by the full subaperture area. Newly allocated; the caller frees it.
*/
dmat* shwfs_saa(const dmat* amp, long npsa){
	real areanormfactor=1./(real)npsa;
	long nsa=NX(amp)/npsa;
	dmat* saa=dnew(nsa, 1);
	for(int isa=0; isa<nsa; isa++){
		const real* ampi=P(amp)+isa*npsa;
		real area=0;
		for(int iamp=0; iamp<npsa; iamp++){
			area+=ampi[iamp]*ampi[iamp];
		}
		P(saa, isa)=area*areanormfactor;
	}
	return saa;
}
/**
   Compute saa[isa]=max(P(P(saa, iwfs), isa)) over iwfs.

   @param[in] saa Subaperture area (or any value) of each WFS.
   @param[in] alpha 1 to compute the maximum. -1 to compute the minimum.
   @return The element wise maximum (or minimum). Newly allocated (or a reference to saa[0] if saa has a single element); the caller frees it. NULL if saa is NULL.
*/
dmat* shwfs_saa_max(dcell* saa, real alpha){
	if(!saa) return NULL;
	if(PN(saa)==1) return dref(P(saa,0));
	long nsa=P(saa, 0)->nx;
	dmat *saamax=dnew(nsa, 1);dset(saamax, -INFINITY);
	for(long iwfs=0; iwfs<NX(saa); iwfs++){
		for(long isa=0; isa<nsa; isa++){
			if(P(saamax, isa)<P(P(saa, iwfs), isa)*alpha){
				P(saamax, isa)=P(P(saa, iwfs), isa)*alpha;
			}
		}
	}
	dscale(saamax, 1./alpha);
	return saamax;
}
/**
   Turn off isolated subapertures: a subaperture with area above the threshold
   is turned off (its area set to 0) if at most two of its eight neighbors are
   above the threshold. Repeated until nothing changes. Nothing is done if there
   are 4 or fewer subapertures.

   @param[in] pts The subapertures. Its map is created and freed.
   @param[in,out] saa Normalized area of each subaperture.
   @param[in] saat Threshold of the normalized area.
*/
void shwfs_sa_isolated(pts_t* pts, dmat* saa, real saat){
	if(pts->nsa<=4) return;
	loc_t* ptsloc=LOC(pts);
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
			if(P(saa, isa)>=saat){
				nself=1;
			}
			for(int jx=-1; jx<2; jx++){
				for(int jy=-1; jy<2; jy++){
					if(jx==0&&jy==0){
						continue;
					}
					long jsa=loc_map_get(ptsloc->map, ix+jx, iy+jy);
					if(jsa&&P(saa, jsa-1)>=saat){
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
/**
   Create the geometric gradient operator (average gradient over each
   subaperture) from OPD on loc to the gradients of the subapertures for each WFS.

   @param[in] loc The points of all subapertures (pts2loc()).
   @param[in] amp The amplitude on loc for each WFS.
   @param[in] saloc Lower left corner of each subaperture.
   @param[in] saa Normalized subaperture area for each WFS.
   @param[in] saat Subapertures with saa below saat have no gradient.
   @return The operator for each WFS, which gives x gradients of all subapertures followed by y. Newly allocated; the caller frees it.
*/
dspcell* shwfs_mkg(loc_t* loc, dcell* amp, loc_t* saloc, dcell* saa, real saat){
	dspcell* GS0=dspcellnew(NX(amp), 1);
OMP_FOR(NX(amp))
	for(int jwfs=0; jwfs<(int)NX(amp); jwfs++){
		P(GS0, jwfs)=mkg(loc, loc, P(amp, jwfs), saloc, PR(saa, jwfs), saat, 1, 0, 0, 1);
	}
	return GS0;
}
/**
   Create the inverse of the piston, tip, and tilt mode covariance of each
   subaperture for each WFS, used by shwfs_grad_geom() to compute the zernike tilt
   (best fit tip/tilt) gradients.

   @param[in] pts The subapertures.
   @param[in] amp The amplitude of the points in pts for each WFS.
   @return The 3x3 inverse covariance of each subaperture for each WFS. Newly allocated; the caller frees it.
*/
dccell* shwfs_saimcc(const pts_t* pts, const dcell* amp){
	dccell* saimcc=(dccell*)cellnew(PN(amp), 1);
	for(int imcc=0; imcc<PN(amp); imcc++){
		dcell* mcc=pts_mcc_ptt(pts, P(P(amp, imcc)));
		P(saimcc, imcc)=dcellinvspd_each(mcc);
		dcellfree(mcc);
	}
	return saimcc;
}
/**
   Compute the geometric gradients of a Shack-Hartmann WFS from the OPD and add
   them to grad: the zernike tilt (best fit) if imcc is set, otherwise the
   average gradient with GS0.

   @param[in,out] grad Gradients, x of all subapertures followed by y. Allocated if NULL.
   @param[in] pts The subapertures.
   @param[in] opd OPD on the points of pts (pts2loc()).
   @param[in] amp Amplitude on the points of pts. Used with imcc.
   @param[in] GS0 Average gradient operator from shwfs_mkg(). Used if imcc is NULL.
   @param[in] imcc Inverse mode covariance of each subaperture from shwfs_saimcc(), or NULL.
*/
void shwfs_grad_geom(dmat** grad, const pts_t* pts, const dmat* opd, const dmat* amp, const dsp* GS0, const dcell* imcc){
	if(imcc){
		pts_ztilt(grad, pts, imcc, P(amp), P(opd));
	} else{
		dspmm(grad, GS0, opd, "nn", 1);
	}
}
/**
   Create the mask (weighting) for the CoG of each subaperture image: the
   image normalized by thres and clipped at 1.

   @param[in] i0 Subaperture images (nsa x nwfs).
   @param[in] thres The image value above which the weight is 1.
   @return The mask with the same dimensions as i0. Newly allocated; the caller frees it.
*/
dcell* shwfs_cog_mask(const dcell* i0, real thres){
	dcell* cogmask=dcellnew_same(NX(i0), NY(i0), NX(i0, 0), NY(i0,0));
	for(int icol=0; icol<NY(i0); icol++){
		for(int isa=0; isa<NX(i0); isa++){
			const dmat *i0i=P(i0, isa, icol);
			dmat *cogmaski=P(cogmask, isa, icol);
			for(int i=0; i<PN(i0i); i++){
				real ratio=P(i0i, i)/thres;
				if(ratio>1.) ratio=1.;
				P(cogmaski, i)=ratio;
			}
		}
	}
	return cogmask;
}
/**
   Compute the noise equivalent angle (NEA) of the CoG of each subaperture of
   a WFS with a Monte Carlo simulation of the noise with cog_nea().

   @param[in] i0 Subaperture images of the WFS (nsa).
   @param[in] cogmask Mask for the CoG of each subaperture (nsa). May be NULL.
   @param[in] bkgrnd2 Background image of each subaperture (nsa). May be NULL.
   @param[in] bkgrnd2c Background image removed by calibration of each subaperture (nsa). May be NULL.
   @param[in] srot Angle to rotate the NEA of each subaperture (for radial pixels). May be NULL.
   @param[in] nsa Number of subapertures.
   @param[in] cogthres Threshold of the CoG (see dcog()).
   @param[in] cogoff Offset subtracted in the CoG (see dcog()).
   @param[in] ntry Number of Monte Carlo trials.
   @param[in,out] rstat Random number generator.
   @param[in] bkgrnd Background level per pixel.
   @param[in] bkgrndc Background level per pixel removed by calibration.
   @param[in] rne Read out noise.
   @param[in] pixthetax Pixel size along x (or radial) in radian.
   @param[in] pixthetay Pixel size along y (or azimuthal) in radian.
   @return nsa x 3 matrix with the xx, yy, and xy NEA (radian^2) of each subaperture. Newly allocated; the caller frees it.
*/
dmat* shwfs_cog_nea(dmat** i0, dmat** cogmask, dmat** bkgrnd2, dmat** bkgrnd2c, const real* srot, long nsa,
	real cogthres, real cogoff, int ntry, rand_t* rstat, real bkgrnd, real bkgrndc, real rne, real pixthetax, real pixthetay){
	dmat* psanea=dnew(nsa, 3);
	dmat* nea=dnew(2, 2);
	for(long isa=0; isa<nsa; isa++){
		dmat* bkgrnd2i=bkgrnd2?bkgrnd2[isa]:NULL;
		dmat* bkgrnd2ic=bkgrnd2c?bkgrnd2c[isa]:NULL;
		dzero(nea);
		cog_nea(P(nea), i0[isa], cogmask?cogmask[isa]:NULL, cogthres, cogoff, ntry, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
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
	return psanea;
}
/**
   Add the CoG of the subaperture images of a WFS to its gradient offset, so
   that the CoG of these images gives zero gradient.

   @param[in,out] gradoff Gradient offset, x of all subapertures followed by y.
   @param[in] i0 Subaperture images of the WFS (nsa).
   @param[in] cogmask Mask for the CoG of each subaperture (nsa). May be NULL.
   @param[in] nsa Number of subapertures.
   @param[in] cogthres Threshold of the CoG (see dcog()).
   @param[in] cogoff Offset subtracted in the CoG (see dcog()).
   @param[in] pixthetax Pixel size along x in radian.
   @param[in] pixthetay Pixel size along y in radian.
*/
void shwfs_cog_gradoff(dmat* gradoff, dmat** i0, dmat** cogmask, long nsa, real cogthres, real cogoff, real pixthetax, real pixthetay){
	for(long isa=0; isa<nsa; isa++){
		real grad[2];
		dcog(grad, i0[isa], 0, 0, cogthres, cogoff, 0, cogmask?cogmask[isa]:NULL);
		P(gradoff, isa)+=grad[0]*pixthetax;
		P(gradoff, isa+nsa)+=grad[1]*pixthetay;
	}
}
