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
#include <unistd.h>
#include "pywfs.h"
#include "accphi.h"
#include "zernike.h"
#include "turbulence.h"
#include "libmisc.h"
#include "cure.h"
static double PYWFS_PSIZE=3; //Size of pywfs detector image. sub-pupil center to center dist is PYWFS_SIZE/2
int PYWFS_DEBUG=0;//debugging implementation
static int PYWFS_TT_DUAL=0;//average tip/tilt response along +/-
static int PYWFS_FULL=0;//make saloc entire focal plane
/**
   \file pywfs.h

   Setup pyramid WFS and do simulation.
*/
/**
 * Create detector pixel to subaperture mapping
 * nside==4: (regular 4-side pyramid)
 * 	2 3
 * 	0 1
 * nside==3: (3-sided pyramid)
 * 1  2
 *   0
 * nside==2: (roof)
 * 0 1
 * nside==1: (zernike sensor)
 * 0
 * */
static void pywfs_mksi(pywfs_t *pywfs, loc_t *loc_fft, loc_t *saloc0, real pupelong){
	dspcellfree(pywfs->si);
	cellfree(pywfs->msaloc);
	const int pyside=pywfs->cfg->nside;
	pywfs->si=dspcellnew(pyside, 1);
	if(!pywfs->sioff){
		pywfs->sioff=dnew(pyside, 2);
	}
	if(pupelong){//msaloc is also used by cuda code for direct interpolation
		pywfs->msaloc=loccellnew(pyside, 1);
	}
	const real dxp=loc_fft->dx;
	const real dsa=saloc0->dx;
	const long notf2=NX(pywfs->nominal)/2;
	for(int ind=0; ind<pyside; ind++){
		const int iy=ind/2;
		const int ix=ind%2;
		loc_t *saloc=0;
		real shx=0, shy=0;
		if(pywfs->pupilshift){
			shx=P(pywfs->pupilshift, ind, 0)*dsa;
			shy=P(pywfs->pupilshift, ind, 1)*dsa;
		}
		real offx=0, offy=0;
		switch(pyside){
		case 1:
			offx=0;
			offy=0;
			break;
		case 2:
			offx=ix-0.5;
			offy=0;
			break;
		case 3:
			if(ind==0){
				offx=0;
				offy=-0.5;
			} else{
				offx=(ind-1.5)*sqrt(3.)*0.5;
				offy=0.25;
			}
			break;
		case 4:
			offx=ix-0.5;
			offy=iy-0.5;
			break;
		default:
			error("Invalid dbg.pwfs_side=%d\n", pyside);
		}
		if(pupelong){//pupil elongation (along radial direction)
			real angle=atan2(offy, offx);
			real frac=1.-pupelong;
			saloc=P(pywfs->msaloc, ind)=locdup(saloc0);
			//squeeze the detector pixel coordinate radially to simulate pupil elongation
			locstretch(saloc, angle, frac);
		} else{
			saloc=saloc0;
		}
		P(pywfs->sioff, ind, 0)=offx;
		P(pywfs->sioff, ind, 1)=offy;
		P(pywfs->si, ind)=mkh(loc_fft, saloc, (offx*notf2)*dxp+shx, (offy*notf2)*dxp+shy, 1., 0);
	}
}
pywfs_t *pywfs_new(pywfs_cfg_t *pycfg, loc_t *loc, const dmat *amp){
	READ_ENV_INT(PYWFS_DEBUG, -10, 10);
	READ_ENV_INT(PYWFS_TT_DUAL, 0, 1);
	READ_ENV_DBL(PYWFS_PSIZE, 1, 10);
	READ_ENV_INT(PYWFS_FULL, 0, 1);
	//check input and fill in optional parameters
	if(!pycfg->order || !pycfg->siglev || !pycfg->wvl){
		error("Some necessary parameters are not provided. order=%d, siglev=%g, wvl=%p\n", pycfg->order, pycfg->siglev, pycfg->wvl);
	}
	
	if(!pycfg->dx){
		if(loc) pycfg->dx=loc->dx; 
		else pycfg->dx=1./64.;
	}else if(pycfg->dx && loc && pycfg->dx!=loc->dx){
		error("pycfg->dx and loc->dx does not match\n");
	}
	if(!pycfg->D){
		if(loc) pycfg->D=loc_diam(loc);
		else error("D must be provided if loc is not present\n");
	}
	if(!pycfg->nside) pycfg->nside=4;
	if(!pycfg->dsa) pycfg->dsa=pycfg->D/pycfg->order;
	if(!pycfg->hs) pycfg->hs=INFINITY;
	if(!pycfg->modulate) pycfg->modulate=5*dmax(pycfg->wvl)/pycfg->D;
	if(!pycfg->modulpos) pycfg->modulpos=(32/pycfg->nside)*pycfg->nside;
	if(!pycfg->modulring) pycfg->modulring=1;
	if(!pycfg->sigmatch) pycfg->sigmatch=2;
	if(!pycfg->poke) pycfg->poke=1e-7;
	if(!pycfg->wvlwts) {pycfg->wvlwts=dnew(PN(pycfg->wvl),1); dset(pycfg->wvlwts, 1./PN(pycfg->wvl));};
	if(!pycfg->pixblur) pycfg->pixblur=0.3;
	if(!pycfg->fieldstop) pycfg->fieldstop=0;
	pywfs_t *pywfs=mycalloc(1, pywfs_t);
	pywfs->cfg=pycfg; //consider copy struct if necessary.
	const int pyside=pycfg->nside;
	const int ng=pywfs_ng(pycfg);

	const real siglev=pycfg->siglev;//signal level for noise free calculation
	const real dx=pycfg->dx;
	//For convenience.
	if(!loc){
		loc=mkannloc(pycfg->D, 0, pycfg->dx, 0);
	}
	if(amp){
		pywfs->amp=dref(amp);
	}else{
		pywfs->amp=dnew(loc->nloc, 1);
		dset(pywfs->amp, 1);
	}
	pywfs->loc=locref(loc);
	int nwvl=NX(pycfg->wvl);
	pywfs->locfft=locfft_init(loc, pywfs->amp, pycfg->wvl, 0, PYWFS_PSIZE, 0);

	long npsf=P(pywfs->locfft->nembed, 0);//size of PSF at the Pyramid tip.
	pywfs->gain=1;
	//size of the part of the PSF captured by pyramid
	const long notf=npsf;//size of OTF array
	const long notf2=notf/2;
	const comp coeff=COMPLEX(0, M_PI*0.5);
	pywfs->pyramid=ccellnew(nwvl, 1);
	dmat *pyramid=dnew(notf, notf);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		P(pywfs->pyramid, iwvl)=cnew(notf, notf);
		cmat *pp=P(pywfs->pyramid, iwvl)/*PCMAT*/;
		const real dtheta=P(pycfg->wvl, iwvl)/(dx*npsf);//PSF sampling
		int nstop=notf;//needed coverage by the Pyramid
		if(pyside==1){
			nstop=round(1.063*dx*npsf/pycfg->D);
			dbg("zernike wfs: nstop=%d\n", nstop);
		} else if(pycfg->fieldstop && isfinite(pycfg->fieldstop)){//Limit fov per wvl
			nstop=ceil(pycfg->fieldstop/dtheta*0.5)*2;
			if(nstop>npsf){
				warning("field stop (%g\") is bigger than PSF size (%g\").\n", pycfg->fieldstop*RAD2AS, dtheta*RAD2AS*npsf);
			}
			dbg("fieldstop is %g\", dtheta is %g\", nstop is %d, npsf is %ld\n", pycfg->fieldstop*RAD2AS, dtheta*RAD2AS, nstop, npsf);
		}
		const long skip=notf>nstop?(notf-nstop)/2:0;
		const real radius2=nstop*nstop*0.25;
		//Make pyramid edge or vertax flat within certain range
		const real eskip=(pycfg->flate/dtheta/2);
		const real vskip=(pycfg->flatv/dtheta/2);
		const real sqrt3=sqrt(3.);
		//const real ratio2=acos(sqrt(0.5))/acos(0.5);
		for(long iy=skip; iy<notf-skip; iy++){
			for(long ix=skip; ix<notf-skip; ix++){
				real xd=labs(ix-notf2);
				real yy=iy-notf2;
				real yd=fabs(yy);
				real opd=0;
				if((xd*xd+yd*yd)<radius2){
					 //not on flat edge
					switch(pyside){
					case 1://zernike
						opd=2.;//opd*coeff=pi
						break;
					case 2://roof with slope
						opd=xd;
						break;
					case 3://3 sided pyrmaid
						//We keep one side of the pyramid the same as slope as a
						//roof in order to shift pupil by half.
						if(yy*sqrt3>xd){
							opd=yy;
						} else{
							//Rotate the coordinate by 120 and apply above surface formula.
							opd=fabs(yy-xd*sqrt3)*0.5;
						}
						break;
					case 4://4-sided pyramid
						if(!(xd<eskip||yd<eskip||(xd<vskip&&yd<vskip))){//pyramid defacts
							opd=(xd+yd);
						}
						break;
					default:
						error("Invalid pywfs.nside=%d\n", pyside);
					}
					P(pp, ix, iy)=cexp(opd*coeff);
					P(pyramid, ix, iy)=opd;//saving only.
				}
			}
		}
	}
	/*if(parms->save.setup){
		writebin(pyramid, "powfs%d_pyramid", ipowfs);
	}*/
	dfree(pyramid);
	//Detector transfer function (sampling onto pixels).
	cmat *nominal=pywfs->nominal=cnew(notf, notf);
	const long order=pycfg->order;
	//when pyside==1 (zernike wfs), dsa needs to be scaled
	const real dsa=(pyside==1?2:1)*pycfg->dsa; //size of detector pixel mapped on pupil
	const real dx2=dx*npsf/notf;//sampling of pupil after inverse fft
	const real du=1./(dx2*notf);
	const real dupix=dsa*du;
	const real pdmeter=pow(dsa/dx2, 2);
	const real pixblur=pycfg->pixblur*dsa;
	const real e0b=-2*pow(M_PI*pixblur*du, 2);
	for(int iy=0; iy<notf; iy++){
		int jy=iy-notf2;
		for(int ix=0; ix<notf; ix++){
			int jx=ix-notf2;
			P(nominal, ix, iy)=sinc(jy*dupix)*sinc(jx*dupix)*pdmeter;
			if(pixblur){
				P(nominal, ix, iy)*=exp(e0b*(jx*jx+jy*jy));
			}
		}
	}
	cfftshift(nominal);
	/*
	cfft2(nominal, -1);
	cfftshift(nominal);
	cfft2(nominal, 1);
	cscale(nominal, 1./(NX(nominal)*NY(nominal)));
	*/
	if(pycfg->psx){
		if(NX(pycfg->psx)!=4){
			error("dbg.pwfs_psx has wrong format: expected 4x1, got %ldx%ld\n",
				NX(pycfg->psx), NY(pycfg->psx));
		}
		if(NX(pycfg->psy)!=4){
			error("dbg.pwfs_psy has wrong format: expected 4x1, got %ldx%ld\n",
				NX(pycfg->psy), NY(pycfg->psy));
		}
		int redefine=1;
		pywfs->pupilshift=dnew(4, 2);
		for(int i=0; i<4; i++){
			real tmp;
			tmp=P(pycfg->psx, i);
			P(pywfs->pupilshift, i, 0)=tmp-redefine?round(tmp):0;
			tmp=P(pycfg->psy, i);
			P(pywfs->pupilshift, i, 1)=tmp-redefine?round(tmp):0;

		}

	}
	loc_t *loc_fft=mksqloc(notf, notf, dx2, dx2, (-notf2+0.5)*dx2, (-notf2+0.5)*dx2);
	const real pupelong=pycfg->pupelong*sqrt(2.)/(order*0.5);
	if(PYWFS_FULL){//sample the full image. for debugging only
		pywfs->saloc=locdup(loc_fft);
	} else{//Pad the grid to avoid missing significant pixels (subapertures).
		long order2=ceil(order)+2*MAX(0, ceil(pycfg->pupelong));
		if(order2>order){
			dbg("Elongated pupil: order %ld increased to %ld.\n", order, order2);
		}
		//Make loc_t symmetric to ensure proper sampling onto detector. Center of subaperture
		pywfs->saloc=mksqloc(order2, order2, dsa, dsa,
			(-order2*0.5+0.5)*dsa, (-order2*0.5+0.5)*dsa);
	}

	pywfs_mksi(pywfs, loc_fft, pywfs->saloc, pupelong);//for each quadrant.

	{
	//Determine subapertures area
		dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
		dmat *ints=0;
		pywfs_ints(&ints, pywfs, opd, siglev);
		//if(parms->save.setup){
		//writebin(ints, "pywfs_ints0");
		//}
		int nints=NX(ints);
		pywfs->saa=dnew(nints, 1);
		for(int i=0; i<NX(ints); i++){
			for(int j=0; j<NY(ints); j++){
				P(pywfs->saa, i)+=P(ints, i, j);
			}
		}
		cellfree(ints);
		dfree(opd);
	}
	if(pycfg->saat>0){
	//Get ride of insufficiently illuminated subapertures.
		dmat *saa=pywfs->saa;
		//real samax=dmaxabs(saa);
		real samean=dsum(saa)/PN(saa);
		loc_reduce(pywfs->saloc, saa, pycfg->saat*samean, 1, 0);
		dscale(saa, NX(saa)/dsum(saa));//saa average to one.
		pywfs_mksi(pywfs, loc_fft, pywfs->saloc, pupelong);
	}
	const int nsa=pywfs->saloc->nloc;
	dscale(pywfs->saa, NX(pywfs->saa)/dsum(pywfs->saa));//saa average to one.
	locfree(loc_fft);
	//Determine the gain and offset of PyWFS
	{
		//offset: grad of a flat wavefront
		dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
		dmat *ints=0;
		dmat *goff=0;
		pywfs_ints(&ints, pywfs, opd, siglev);
		pywfs_grad(&goff, pywfs, ints);
		/*if(parms->save.setup){
			writebin(ints, "powfs%d_piston_ints", ipowfs);
			writebin(goff, "powfs%d_piston_grad", ipowfs);
		}*/
		dadd(&pywfs->gradoff, 1, goff, 1);

		dfree(goff);
		dfree(opd);
		dfree(ints);
		//Determine optical gain.
		dmat *TT=pywfs_tt(pywfs);
		real gsum=dsum(TT);
		real gainscl=ng*nsa/gsum;
		/*
		  pywfs->gain is inverse of optical gain, to insure 1rad of input
		  tip/tilt wavefront gives 1rad of gradient output.*/
		pywfs->gain*=gainscl;
		dscale(pywfs->gradoff, gainscl);
		dscale(TT, gainscl);
		pywfs->GTT=TT;TT=NULL;
		info("pywfs_gain=%g\n", pywfs->gain);
	}
	return pywfs;
}
/**
   Perform FFT over the complex PSF with additional phases caused by the
   pyramid. FFT on each quadrant of the PSF creates diffraction effects.
   @param[in,out] ints	The intensity. Accumulate.
   @param[in] pywfs		PYWFS parameters
   @param[in] opd		The OPD
   @param[in] siglev	The signal level.
*/
void pywfs_ints(dmat **ints, const pywfs_t *pywfs, const dmat *opd, real siglev){
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	locfft_t *locfft=pywfs->locfft;
	ccell *psfs=0;
	locfft_psf(&psfs, locfft, opd, NULL, 1);//psfs.^2 sum to 1. peak in center
	/*if(global->setupdone&&global->parms->plot.run){
		draw("Ints", (plot_opts){ .cimage=P(psfs, 0), .ctype=0, .zlog=1 }, "PWFS PSF", "x", "y", "wfs %d focus", pywfs->iwfs0);
	}*/
	const int nwvl=NX(locfft->wvl);
	const real dx=locfft->loc->dx;
	const long npsf=P(locfft->nembed, 0);
	const long npsf2=npsf/2;
	const dmat *wvlwts=pycfg->wvlwts;
	//position of pyramid for modulation
	const int pos_n=pycfg->modulpos;
	const int pos_nr=pycfg->modulring;
	const real pos_r=pycfg->modulate;
	const long notf=NX(pywfs->nominal);
	const long notf2=notf/2;
	cmat *otf=cnew(notf, notf);
	dmat *pupraw=dnew(notf, notf);
	for(int ir=0; ir<pos_nr; ir++){//Radius of the current ring
		const real pos_ri=pos_r*(ir+1)/pos_nr;
		//Scale number of points by ring size to have even surface brightness
		const int pos_ni=pos_n*(ir+1)/pos_nr;
		const int ipos0=pycfg->modulpos_i>0?(pycfg->modulpos_i-1):0;
		const int ipos1=pycfg->modulpos_i>0?pycfg->modulpos_i:pos_ni;
		for(int ipos=ipos0; ipos<ipos1; ipos++){
			//whether the first point falls on the edge or not makes little difference
			const real theta=2*M_PI*((real)ipos/pos_ni);
			const real posx=cos(theta)*pos_ri;
			const real posy=sin(theta)*pos_ri;
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				real dtheta=P(locfft->wvl, iwvl)/(dx*npsf);
				const long offy=(long)round(posy/dtheta);
				const long offy2=npsf2+offy-notf2;
				const long iy0=MAX(-offy2, 0);
				const long ny2=MIN(notf, npsf-offy2)-iy0;

				const long offx=(long)round(posx/dtheta);
				const long offx2=npsf/2+offx-notf2;
				const long ix0=MAX(-offx2, 0);
				const long nx2=MIN(notf, npsf-offx2)-ix0;

				czero(otf);
				const comp *pyramid=P(P(pywfs->pyramid, iwvl));
				for(long iy=iy0; iy<ny2; iy++){
					for(long ix=ix0; ix<nx2; ix++){
						const long indin=ix+offx2+(iy+offy2)*npsf;
						const long indout=ix+iy*notf;
						P(otf, indout)=P(P(psfs, iwvl), indin)*pyramid[indout];
					}
				}
				cfft2(otf, 1);
				cabs22d(&pupraw, 1., otf, P(wvlwts, iwvl)/(notf*notf*pos_ni*pos_nr));
			}//for iwvl
		}//for ipos
	}//for ir
	/*if(global&&global->setupdone&&global->parms->plot.run){
		draw("Ints", (plot_opts){ .image=pupraw, .zlog=1 }, "PWFS Pupil", "x", "y", "wfs %d pupil", pywfs->iwfs0);
	}*/
	//writebin(pupraw, "cpu_psf"); exit(0);
	ccpd(&otf, pupraw);//pupraw sum to one.
	//writebin(otf, "cpu_wvf4");
	dfree(pupraw);
	cfft2(otf, -1);
	//writebin(otf, "cpu_wvf5");
	ccwm(otf, pywfs->nominal);
	cfft2(otf, 1);
	const int nsa=P(pywfs->si, 0)->nx;
	if(!(*ints)){
		(*ints)=dnew(nsa, NX(pywfs->si));
	}
	for(int i=0; i<NX(pywfs->si); i++){
		//normalized so that each "subaperture" sum to 1.
		dspmulcreal(PCOL(*ints, i), P(pywfs->si, i), P(otf), (real)nsa*siglev/(notf*notf));
	}
	//writebin(*ints, "cpu_ints"); exit(0);
	ccellfree(psfs);
	cfree(otf);
}
/**
   Compute gradients. It replaces the result, not accumulate.
 */
void pywfs_grad(dmat **pgrad, const pywfs_t *pywfs, const dmat *ints){
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	const long nsa=NX(ints);
	const int ng=pywfs_ng(pycfg);
	if(!*pgrad){
		*pgrad=dnew(nsa*ng, 1);
	}
	dmat *grad=*pgrad;
	real *pgx=P(grad);
	real *pgy=P(grad)+nsa;
	real gain=pywfs->gain;
	real triscalex=sqrt(3.)/2;
	real imean=0;
	if(pycfg->sigmatch==2){
		imean=dsum(ints)/nsa;
	}
	for(int isa=0; isa<nsa; isa++){
		if(pycfg->raw||pycfg->nside<3){
			real alpha2=gain/pycfg->siglev;
			for(int iside=0; iside<pycfg->nside; iside++){
				P(grad, isa+nsa*iside)=P(ints, isa, iside)*alpha2;
			}
		} else{
			real isum=0;//subaperture intensity for normalization
			switch(pycfg->sigmatch){
			case 0:
				info_once("PWFS: No siglev correction.\n");
				isum=pycfg->siglev*P(pywfs->saa, isa);
				break;
			case 1:
				info_once("PWFS: Individual siglev correction.\n");
				for(int i=0; i<pycfg->nside; i++){
					isum+=P(ints, isa, i);
				}
				break;
			case 2:
				info_once("PWFS: Global siglev correction.\n");//preferred.
				isum=imean*P(pywfs->saa, isa);
				break;
			}
			real alpha2=gain/isum;
			switch(pycfg->nside){
			case 3:
				pgx[isa]=(P(ints, isa, 1)-P(ints, isa, 2))*alpha2*triscalex;
				pgy[isa]=(P(ints, isa, 0)-0.5*(P(ints, isa, 1)+P(ints, isa, 2)))*alpha2;
				break;
			case 4:
				pgx[isa]=(P(ints, isa, 0)-P(ints, isa, 1)
					+P(ints, isa, 2)-P(ints, isa, 3))*alpha2;
				pgy[isa]=(P(ints, isa, 0)+P(ints, isa, 1)
					-P(ints, isa, 2)-P(ints, isa, 3))*alpha2;
				break;
			}
		}
	}

	if(pywfs->gradoff){
		dadd(pgrad, 1, pywfs->gradoff, -1);
	}
}
/**
   Return measurement of T/T mode, normalized for 1 unit of input.
*/
dmat *pywfs_tt(const pywfs_t *pywfs){
	TIC;tic;info("Computing pywfs_tt...");
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	const loc_t *loc=pywfs->locfft->loc;
	dmat *opd=dnew(loc->nloc, 1);
	dmat *ints=0;
	const long nsa=P(pywfs->si, 0)->nx;
	const int ng=pywfs_ng(pycfg);
	dmat *out=dnew(nsa*ng, 2);
	dmat *gradx=drefcols(out, 0, 1);
	dmat *grady=drefcols(out, 1, 1);

	real ptt[3]={0,0,0};
	const real alpha=0.005*AS2RAD;
	const real siglev=pycfg->siglev;//signal level for noise free calculation
	//+x
	ptt[1]=alpha;  ptt[2]=0;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_ints(&ints, pywfs, opd, siglev);
	pywfs_grad(&gradx, pywfs, ints);
	if(PYWFS_DEBUG){
		writebin(ints, "pywfs_ttx_ints");
		writebin(gradx, "pywfs_ttx_grad");
	}
	//+y
	ptt[1]=-alpha; ptt[2]=alpha;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_ints(&ints, pywfs, opd, siglev);
	pywfs_grad(&grady, pywfs, ints);
	if(PYWFS_DEBUG){
		writebin(ints, "pywfs_tty_ints");
		writebin(grady, "pywfs_tty_grad");
	}

	if(PYWFS_TT_DUAL){
		dmat *gradx2=dnew(nsa*2, 1);
		dmat *grady2=dnew(nsa*2, 1);
		//-x
		ptt[1]=-alpha; ptt[2]=-alpha;
		loc_add_ptt(opd, ptt, loc);
		dzero(ints);
		pywfs_ints(&ints, pywfs, opd, siglev);
		pywfs_grad(&gradx2, pywfs, ints);
		if(PYWFS_DEBUG){
			writebin(ints, "pywfs_ttx2_ints");
			writebin(gradx2, "pywfs_ttx2_grad");
		}
		//-y
		ptt[1]=+alpha; ptt[2]=-alpha;
		loc_add_ptt(opd, ptt, loc);
		dzero(ints);
		pywfs_ints(&ints, pywfs, opd, siglev);
		pywfs_grad(&grady2, pywfs, ints);
		if(PYWFS_DEBUG){
			writebin(ints, "pywfs_tty2_ints");
			writebin(grady, "pywfs_tty2_grad");
		}
		dadd(&gradx, 1, gradx2, -1);
		dadd(&grady, 1, grady2, -1);
		dscale(out, 0.5/alpha);
		dfree(gradx2);
		dfree(grady2);
	} else{
		dscale(out, 1./alpha);
	}
	dfree(gradx);
	dfree(grady);
	dfree(opd);
	dfree(ints);
	toc2("done");
	return out;
}
static uint32_t pywfs_hash(const pywfs_t *pywfs, uint32_t key){
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	key=lochash(pywfs->loc, key);
	key=dhash(pywfs->amp, key);
	key=dhash(pywfs->saa, key);
	key=dhash(pycfg->wvl, key);
	key=dhash(pycfg->wvlwts, key);
	key=chash(P(pywfs->pyramid, 0), key);
	if(pywfs->pupilshift){
		key=dhash(pywfs->pupilshift, key);
	}
	key=hashlittle(&pycfg->raw, sizeof(int), key);
	key=hashlittle(&pycfg->nside, sizeof(int), key);
	return key;
}
dmat *(*pywfs_mkg_ext)(const pywfs_t *pywfs, const loc_t *locin, const loc_t *locfft, const dmat *mod, real displacex, real displacey)=NULL;
/**
   There is no need to simulate turbulence to fake optical gain. Any optical
   gain can be used as long as it "correct", i.e., a radian of tilt produces one
   radian of tilt, which is gauranteed when pywfs->gain is computed under the
   same conditions.
 */
static dmat *pywfs_mkg_do(const pywfs_t *pywfs, const loc_t *locin, const loc_t *locfft, const dmat *mod,
	real displacex, real displacey){
#if USE_CUDA
	if(pywfs_mkg_ext){
		return pywfs_mkg_ext(pywfs, locin, locfft, mod, displacex, displacey);
	}
#endif
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	const int nsa=P(pywfs->si, 0)->nx;
	const int ng=pywfs_ng(pycfg);
	dmat *grad0=dnew(nsa*ng, 1);
	dmat *opd0;
	if(pywfs->opdadd){
		opd0=dref(pywfs->opdadd);
	} else{
		opd0=dnew(locfft->nloc, 1);
	}
	{
		dmat *ints=0;
		pywfs_ints(&ints, pywfs, opd0, pycfg->siglev);
		pywfs_grad(&grad0, pywfs, ints);
		//writebin(grad0, "grad0_cpu");
		//writebin(ints, "ints0_cpu");
		dfree(ints);
	}
	unsigned int count=0;
	const int nmod=mod?NY(mod):locin->nloc;
	dmat *ggd=dnew(nsa*ng, nmod);
	if(mod&&NX(mod)!=locin->nloc){
		error("NX(mod) must equal to %ld", locin->nloc);
	}

	const real scale=1.-locin->ht/pycfg->hs;
	TIC;tic;
	OMP_FOR(NTHREAD)
	for(int imod=0; imod<nmod; imod++){
		dmat *opdin=dnew(locin->nloc, 1);
		dmat *opdfft=ddup(opd0);
		dmat *ints=0;
		dmat *grad=drefcols(ggd, imod, 1);
		real poke=pycfg->poke;
		if(mod){
			dmat *tmp=drefcols(mod, imod, 1);
			//the equivalent radimodl order of zernike.
			//real radial=ceil((sqrt(8.*(imod+1)+1)-3)*0.5)+1;
			//real std=dstd(tmp);
			real tmax, tmin;
			dmaxmin(tmp, &tmax, &tmin);
			poke/=(tmax-tmin);//sqrt(radial);
			dadd(&opdin, 0, tmp, poke);
			dfree(tmp);
		} else{
			P(opdin, imod)=poke;
		}
		prop_nongrid((loc_t *)locin, P(opdin), locfft, P(opdfft), 1, displacex, displacey, scale, 0, 0);
		//writebin(opdfft, "phiout_cpu_%d", imod);
		pywfs_ints(&ints, pywfs, opdfft, pycfg->siglev);
		//writebin(ints, "ints_cpu_%d", imod);
		pywfs_grad(&grad, pywfs, ints);
		dadd(&grad, 1, grad0, -1);
		dscale(grad, 1./poke);
		atomic_add_fetch(&count, 1);
		if(count%((nmod+9)/10)==0){
			real ts=myclockd()-tk;
			info2("%d of %d. %.2f of %.2f seconds.\n", count, nmod, ts, ts*nmod/count);
		}
		dfree(opdfft);
		dfree(opdin);
		dfree(grad);
		dfree(ints);
	}
	info2("\n");
	dfree(grad0);
	dfree(opd0);
	return ggd;
}
/**
   locin is on pupil.
 */
dmat *pywfs_mkg(pywfs_t *pywfs, const loc_t *locin, const char *distortion, const dmat *mod, const dmat *opdadd,
	real displacex, real displacey){
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	if(opdadd){
		dfree(pywfs->opdadd);
		pywfs->opdadd=dnew(pywfs->locfft->loc->nloc, 1);
		prop_nongrid(pywfs->loc, P(opdadd), pywfs->locfft->loc, P(pywfs->opdadd), 1, 0, 0, 1, 0, 0);
	}
	loc_t *locfft=pywfs->locfft->loc;
	if(distortion){
		locfft=loctransform(locfft, distortion);
	}
	/*if(mod && NY(mod)<=6){
	return pywfs_mkg_do(pywfs, locin, locfft, mod, displacex, displacey);
	}*/
	uint32_t key=0;
	key=lochash(locin, key);
	key=pywfs_hash(pywfs, key);
	if(mod) key=dhash(mod, key);
	if(opdadd) key=dhash(opdadd, key);
	char fn[PATH_MAX-10];
	snprintf(fn, sizeof(fn), "G/G_%ld_%ld_%ld_%g_%d_%g_%g_%g_%g_%g_%u_v2.bin", 
		P(pywfs->locfft->nembed, 0), locin->nloc, NY(mod), pycfg->modulate*RAD2MAS, pycfg->modulpos,
		locin->iac, displacex, displacey, 1., pycfg->poke*1e9, key);
	dmat *gg=0;
	CACHE_FILE(gg, fn, dread, ({gg=pywfs_mkg_do(pywfs, locin, locfft, mod, displacex, displacey);}), writebin);
	if(distortion){
		locfree(locfft);
	}
	return gg;
}
/**
   frees pywfs_t
*/
void pywfs_free(pywfs_t *pywfs){
	if(!pywfs) return;
	cellfree(pywfs->amp);
	locfree(pywfs->loc);

	locfft_free(pywfs->locfft);
	cellfree(pywfs->pyramid);
	cellfree(pywfs->nominal);
	cellfree(pywfs->si);
	cellfree(pywfs->sioff);
	cellfree(pywfs->opdadd);
	cellfree(pywfs->gradoff);
	cellfree(pywfs->GTT);
	cellfree(pywfs->pupilshift);
	cellfree(pywfs->msaloc);
	cellfree(pywfs->saa);
	free(pywfs);
}
void pycfg_free(pywfs_cfg_t *pycfg){
	if(!pycfg) return;
	cellfree(pycfg->wvl);
	cellfree(pycfg->wvlwts);
	dfree(pycfg->psx);
	dfree(pycfg->psy);
	free(pycfg);
}
/**
 * A convenient wrapper to simulate pywfs from Python.
 * @param[out] ints		Subaperture images
 * @param[out] grad		Subaperture gradients
 * @param[in] pycfg		Configuration struct. Optional if order, wvl, siglev are all set
 * @param[in] order		PWFS order. optional if pycfg is set
 * @param[in] wvl		Wavelength. optional if pycfg is set
 * @param[in] siglev	Signal level of each subaperture. optional if pycfg is set
 * @param[in] loc		OPD loc. 
 * @param[in] amp		OPD amplitude map.
 * @param[in] opd		OPD.
*/
void pywfs_simu(dmat **ints, dmat **grad, pywfs_cfg_t *pycfg, int order, dmat *wvl, real siglev, loc_t *loc, const dmat *amp, const dmat *opd){
	int free_pycfg=0;
	if(!pycfg){
		pycfg=mycalloc(1, pywfs_cfg_t);
		free_pycfg=1;
	}
	if(!pycfg->order) pycfg->order=order;
	if(!pycfg->wvl) pycfg->wvl=dref(wvl);
	if(!pycfg->siglev) pycfg->siglev=siglev;
	
	pywfs_t *pywfs=pywfs_new(pycfg, loc, amp);
	dmat *ints2=NULL;
	if(!ints) ints=&ints2;
	pywfs_ints(ints, pywfs, opd, siglev);
	pywfs_grad(grad, pywfs, *ints);
	pywfs_free(pywfs);
	dfree(ints2);
	if(free_pycfg) pycfg_free(pycfg);
}
/**
 * PYWFS gain calibrate using measurements.
 * The idea is as follows.
 * With a given input OPD, the PWFS measures OPDR.
 * Feed OPDR to PWFS, it measures OPDR2. Optionally add a fitting error to the input.
 * Scale OPDR by alpha, it measures OPDR3.
 * Iterate the alpha until OPDR3 has the same magnitude as OPDR.
 * The final alpha is the gain adjustment.
 * This does not quite work because the original input OPD and OPDR has different spatial frequency content.
 * Most of OPD is not measurable.
 */
void pywfs_gain_calibrate(pywfs_t *pywfs, const dmat *grad, real r0){
	dmat *ints=NULL;
	dmat *grad2=NULL;
	
	dmat *opd2=dnew(pywfs->locfft->loc->nloc, 1);//resampled opdr onto locfft
	dmat *opdr2=NULL;//2nd reconstructed OPD
	real amp2r;
	{	//compute measured OPDR and resample to input grid
		dmat *opdr=NULL; //reconstructed OPD
		cure_loc(&opdr, grad, pywfs->saloc);
		amp2r=dsumsq(opdr);
		map_t *opdr_map=map_convert(opdr);opdr=NULL;//reference OPD into a map_t
		prop_grid(opdr_map, pywfs->locfft->loc, P(opd2), 1, 0, 0, 1, 0, 0, 0);
		mapfree(opdr_map);
	}
		
	real gaintot=1;
	dmat *opd2fit=r0?genatm_loc(pywfs->locfft->loc, r0, 1, -11./3., 3):NULL;
	dmat *opd2tot=NULL;
	real siglev=pywfs->cfg->siglev;
	for(int i=0; i<5; i++){
		//compute PYWFS gain error with opdr and fix it
		dzero(ints);
		dadd(&opd2tot, 0, opd2, gaintot);
		dadd(&opd2tot, 1, opd2fit, 1);
		pywfs_ints(&ints, pywfs, opd2tot, siglev);
		pywfs_grad(&grad2, pywfs, ints);
		cure_loc(&opdr2, grad2, pywfs->saloc);
		//writebin(opdr2, "opdr2_%d", i);
		real gain=sqrt(amp2r/dsumsq(opdr2));//make the output agree with original opdr.
		gaintot*=gain;
		dbg("gain=%.3f, total gain adjustment is %.3f\n", gain, gaintot);
	}
	pywfs->gain*=gaintot;
	dfree(ints);
	dfree(grad2);
	
	dfree(opd2);
	dfree(opdr2);
	
	dfree(opd2fit);
	dfree(opd2tot);
}
extern int PYWFS_DEBUG;
///Test PYWFS implementation
void pywfs_test(pywfs_t *pywfs){
	if(!PYWFS_DEBUG) return;
	if(!pywfs) return;
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	real siglev=pycfg->siglev;
	if(abs(PYWFS_DEBUG)==1){//Test linearity of PWFS with a zernike mode.
		dmat* ints=0;
		real wve=1e-9*20;
		dmat* opds=NULL;
		if(PYWFS_DEBUG==-1){//replace with petal modes
			opds=dnew(pywfs->locfft->loc->nloc, 7);
			real *px=pywfs->locfft->loc->locx;
			real *py=pywfs->locfft->loc->locy;
			real dang=M_PI/3;
			real ang0=M_PI;
			for(int ix=0; ix<NX(opds);ix++){
				real angle=atan2(px[ix], py[ix]);//intentionally reverse x/y to match the pupil amplitude map
				int ip=ifloor((angle+ang0)/dang)+1;
				for(int iy=1; iy<7; iy++){
					P(opds, ix, iy)=0;
				}
				P(opds, ix, ip)=1;
			}
		}else{
			opds=zernike(pywfs->locfft->loc, 0, 2, 2, -5);
			warning("Using zernike for pywfs gain testing\n");
		}
		writebin(opds, "pywfs_input_opds");
		dmat* grad=0;
		int nn=1;
		dmat *atm=NULL;
		{
			nn=40;
			real r0=0.186;
			real L0=30;
			atm=genatm_loc(pywfs->locfft->loc, r0, L0, -11./3., 1);
		}
		
		zfarr *pupsave=zfarr_init(nn,NY(opds), "pywfs_modal_ints");
		zfarr *grads=zfarr_init(nn,NY(opds), "pywfs_modal_grad");
		dmat *opd=dnew(NX(opds),1);

		for(int im=0; im<NY(opds); im++){
			dmat *opdi=drefcols(opds, im, 1);
			for(int j=0; j<nn; j++){
				info2("im=%d, j=%d\n", im, j);
				//for(int posi=0; posi<pywfs->cfg->modulpos; posi++){
					//((pywfs_cfg_t*)pywfs->cfg)->modulpos_i=posi+1;//for testing modulate subframe
					dzero(opd);
					if(j%2==1){
						dadd(&opd, 1, opdi, wve);
					}
					if(atm){
						dadd(&opd, 1, atm, (j/2)*0.05);
					}
					dzero(ints);
					pywfs_ints(&ints, pywfs, opd, siglev);
					pywfs_grad(&grad, pywfs, ints);
					zfarr_push(pupsave, 0, ints);
					zfarr_push(grads, 0, grad);
				//}
			}
			dfree(opdi);
		}
		zfarr_close(pupsave);
		zfarr_close(grads);
		cellfree(ints);
		dfree(grad);
		dfree(opds);
		dfree(opd);
		dfree(atm);
	}else if(PYWFS_DEBUG==2){//Test linearity of a zenike mode with noise
		real wve=20e-9;
		dmat* opds=zernike(pywfs->locfft->loc, 0, 0, 0, -5);
		dmat* opdi=0;
		dmat* ints=0, * grad=0;
		dadd(&opdi, 0, opds, wve);
		pywfs_ints(&ints, pywfs, opdi, siglev);
		pywfs_grad(&grad, pywfs, ints);
		dmat* reg=dpinv(grad, 0);
		writebin(opds, "pywfs_dither_opd");
		writebin(ints, "pywfs_dither_ints");
		writebin(grad, "pywfs_dither_grad");

		rand_t rstat;
		seed_rand(&rstat, 1);
		dmat* tmp=0;
		int nj=10, nn=10;
		dmat* res=dnew(nj, nn);
		dmat* ints2=0;
		for(int j=0; j<nj; j++){
			dzero(ints);
			dadd(&opdi, 0, opds, wve*(j+1));
			pywfs_ints(&ints, pywfs, opdi, siglev);
			for(int in=0; in<nn; in++){
				dadd(&ints2, 0, ints, 1);
				addnoise(ints2, &rstat, 0, 0, 0, 0, 0, in, 1);
				pywfs_grad(&grad, pywfs, ints2);
				dmm(&tmp, 0, reg, grad, "nn", 1);
				P(res, j, in)=P(tmp,0);
				info2("%d of %d, %d of %d: %g\n", j, nj, in, nn, P(tmp,0));
			}
		}

		writebin(res, "pywfs_dither_response");
		dfree(opds); dfree(opdi); dfree(ints); dfree(grad); dfree(reg); dfree(tmp); dfree(res); dfree(ints2);
	}else if(PYWFS_DEBUG==3){//Test NCPA calibration
		dmat* opdatm=dread("opdatm");
		dmat* opdbias_full=dread("opdbias_full");
		dmat* opdbias_astigx=dread("opdbias_astigx");
		dmat* opdbias_polish=dread("opdbias_polish");
		const real atmscale=1;
		for(int i=0; i<100; i++){
			info2("%d ", i);
			dmat* ints=0;
			dmat* grad=0;
			dmat* opd=0;

			dadd(&opd, 0, opdatm, (i+1)*0.02);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "grad_atm_%d", i);

			dadd(&opd, 0, opdbias_full, (i+1)*0.02);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradbias_full_%d", i);

			dadd(&opd, 0, opdbias_astigx, (i+1)*0.02);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradbias_astigx_%d", i);

			dadd(&opd, 0, opdbias_polish, (i+1)*0.02);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradbias_polish_%d", i);

			dadd(&opd, 0, opdbias_full, (i+1)*0.02);
			dadd(&opd, 1, opdatm, atmscale);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradboth_full_%d", i);


			dadd(&opd, 0, opdbias_astigx, (i+1)*0.02);
			dadd(&opd, 1, opdatm, atmscale);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradboth_astigx_%d", i);

			dadd(&opd, 0, opdbias_polish, (i+1)*0.02);
			dadd(&opd, 1, opdatm, atmscale);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradboth_polish_%d", i);
			dfree(opd);
			dfree(ints);
			dfree(grad);


			dadd(&opd, 0, opdbias_full, (i+1)*0.02);
			dadd(&opd, 1, opdatm, (i+1)*0.02);
			dzero(ints);
			pywfs_ints(&ints, pywfs, opd, siglev);
			pywfs_grad(&grad, pywfs, ints);
			writebin(grad, "gradall_%d", i);
			dfree(opd);
			dfree(ints);
			dfree(grad);

		}
	}else if(PYWFS_DEBUG==4){//test automatic gain calibration 2nd version
		dmat *ints=NULL;
		dmat *grad=NULL;
		real r0=0.2;
		real L0=2;//smaller value to simulate residual
		dmat *opd=genatm_loc(pywfs->locfft->loc, r0, L0, -11./3., 1);
		pywfs_ints(&ints, pywfs, opd, siglev);
		pywfs_grad(&grad, pywfs, ints);
		pywfs_gain_calibrate(pywfs, grad, r0);
		//writebin(opd, "opd_in");
		dfree(ints);
		dfree(grad);
		dfree(opd);
	}
	if(PYWFS_DEBUG){
		exit(0);
	}
}
