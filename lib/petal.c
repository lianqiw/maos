/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

#include "../math/mathdef.h"
#include "petal.h"
#include "phase.h"
#include "mkdtf.h"
/**
 * Create petal piston mode to OPD interaction matrix and/or Convert petal modes to petal OPD.
 * no longer included as common tip/tilt can be recreated with global tip/tilt
 * and differential petal tip/tilt creates high order modes that can be sensed.
 * @param[out]ph, 	Outupt the interaction matrix from petal modes to OPD [can be NULL]
 * @param[out] popd Output the opd from petal mode vector [can be NULL]
 * @param loc		The phase screen grid. Replacing nx, ny, cx, cy.
 * @param nx, 	    Size of the phase screen grid.
 * @param ny,   	Size of the phase screen grid.
 * @param cx,       Center of aperture 
 * @param cy,       Center of aperture 
 * @param npetal,   Number of petals.
 * @param theta, 	Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.
 * @param mode  	If set, output opd with the mode vector.
){
*/
static void petal_do(dsp **ph, dmat **popd, loc_t *loc, long nx, long ny, real cx, real cy, long npetal, real theta, const dmat *mode){
	if(!ph&&!popd){
		warning("Either ph or popd should be set. No action\n");
		return;
	}
	if (loc){
		nx=loc->nloc;
		ny=1;
	}
	long ncol=npetal;
	real dtheta=TWOPI/npetal;//central angle
	
	theta=M_PI*1.5-theta;
	while(theta<M_PI) theta+=TWOPI;
	dsp *ht=NULL;
	spint *pp=NULL, *pi=NULL;
	real *px=NULL;
	if(ph){
		ht=dspnew(ncol, nx*ny, nx*ny*ncol/(npetal-1));//transpose of the final result
		pp=ht->pp;
		pi=ht->pi;
		px=ht->px;
	}
	if(popd){
		dinit(popd, nx, ny);
		if(!mode){
			error("Mode must be set in order to output OPD\n");
		}else if(PN(mode)!=npetal){
			error("Mode must have length of %ld (is %ld)\n", npetal, PN(mode));
		}
	}
	//info("nx=%ld, ny=%ld, cx=%g, cy=%g\n", nx, ny, cx, cy);
	spint count=0;
	long icol=0;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			real thetai;
			if(loc){
				thetai=atan2(loc->locy[ix], loc->locx[ix])+theta;//guaranteed to >0
			}else{
				thetai=atan2(iy-cy, ix-cx)+theta;
			}
			thetai=fmod(thetai, TWOPI);//between 0 and TWOPI
			real fseg=thetai/dtheta;
			long iseg=ifloor(fseg);
			int onedge=0;
			if(fabs(iseg-fseg)<1e-4){
				onedge=-1;
			}else if(fabs(iseg+1-fseg)<1e-4){
				onedge=1;//skip points right at the petal gap
			}
			if(ph){
				pp[icol++]=count;
				px[count]=1; 		   		   
				pi[count]=iseg; 
				count++;
				if(onedge){
					px[count-1]=0.5;
					px[count]=0.5;
					pi[count]=(iseg+npetal+onedge)%npetal;
					count++;
				}
				if(count>ht->nzmax){
					error("hpetal matrix overflows. count=%ld, nzmax=%ld\n", count, ht->nzmax);
				}
			}
			if(popd){
				P(*popd, ix, iy)=P(mode, iseg);
				if(onedge){
					P(*popd, ix, iy)=P(mode, iseg)*0.5;
					P(*popd, ix, iy)=P(mode, (iseg+onedge+npetal)%npetal)*0.5;
				}
			}
		}
	}
	if(ph){
		pp[nx*ny]=count;
		dspsetnzmax(ht, count);
		if(*ph) dspfree(*ph);
		*ph=dsptrans(ht);
		dspfree(ht);
	}
}
/**
 * Create petal piston mode to OPD interaction matrix.
 * @return  	    Outupt the interaction matrix from petal modes to OPD
 * @param loc		The phase screen grid. 
 * @param npetal,   Number of petals.
 * @param theta, 	Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.
*/
dsp *petal_mkh_loc(loc_t *loc, long npetal, real theta){
	dsp *h=NULL;
	petal_do(&h, NULL, loc, 0, 0, 0, 0, npetal, theta, NULL);
	return h;
}
/**
 * Create petal piston mode to OPD interaction matrix.
 * @return  	    Outupt the interaction matrix from petal modes to OPD
 * @param nx, 	    Size of the phase screen grid.
 * @param ny,   	Size of the phase screen grid.
 * @param cx,       Center of aperture
 * @param cy,       Center of aperture
 * @param npetal,   Number of petals.
 * @param theta, 	Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.
 */
dsp *petal_mkh(long nx, long ny, real cx, real cy, long npetal, real theta){
	dsp *h=NULL;
	dbg("nx=%ld, ny=%ld, cx=%g, cy=%g\n", nx, ny, cx, cy);
	petal_do(&h, NULL, NULL, nx, ny, cx, cy, npetal, theta, NULL);
	return h;
}
/**
 * Convert petal modes to petal OPD.
 * @param[out] opd  Output the opd from petal mode vector
 * @param nx, 	    Size of the phase screen grid.
 * @param ny,   	Size of the phase screen grid.
 * @param cx,       Center of aperture
 * @param cy,       Center of aperture
 * @param npetal,   Number of petals.
 * @param theta, 	Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.
 * @param mode  	If set, output opd with the mode vector.
){
*/
void petal_opd(anydmat opd, real cx, real cy, long npetal, real theta, const dmat *mode){
	petal_do(NULL, &opd.dm, NULL, NX(opd.dm), NY(opd.dm), cx, cy, npetal, theta, mode);
}
/**
 * Convert petal modes to petal OPD. Alternative interface
 * @return	  		Output the opd from petal mode vector
 * @param nx, 	    Size of the phase screen grid.
 * @param ny,   	Size of the phase screen grid.
 * @param cx,       Center of aperture
 * @param cy,       Center of aperture
 * @param npetal,   Number of petals.
 * @param theta, 	Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.
 * @param mode  	If set, output opd with the mode vector.
){
*/
dmat *petal_mkopd(long nx, long ny, real cx, real cy, long npetal, real theta, const dmat *mode){
	dmat *opd=NULL;
	petal_do(NULL, &opd, NULL, nx, ny, cx, cy, npetal, theta, mode);
	return opd;
}
/**
 * Struct for petal mode reconstruction
*/
typedef struct petal_t{
	dmat *amp; /**<Pupil amplitude*/
	dsp *hpetal;/**<Petal modes*/
	dmat *mod; /**<Model matrix*/
	dmat *rmod;/**<Pseudo Inverse of mod*/
	lmat *ind; /**<Mapping from mod columns (from 2+) to petals*/
	real fembed;/**<Embedding factor for FFT. Use 2. */
	dmat *otf; /**<Pixel size and blur transfer function to deblur PSF.*/
	int nrep;  /**<Number of gerchberg_saxton iterations*/
	int npetal;/**<Number of petals*/
	int npsf;  /**<Number of pixels to use*/
	int nsa;   /**<Number of subapertures in this group.*/
	int withtt;/**<Whether mod include tip/tilt*/
}petal_t;
/**
 * Free the petal_t array
*/
void petal_free(petal_t *p, int nsa){
	if(!p) return;
	if(!nsa) nsa=p->nsa;
	if(!nsa) nsa=1;
	for(int isa=0; isa<nsa; isa++){
		dfree(p[isa].amp);
		dspfree(p[isa].hpetal);
		dfree(p[isa].mod);
		dfree(p[isa].rmod);
		lfree(p[isa].ind);
	}
	free(p);
}

/**
 * Connect petal measurements by TTF quadrants.
 * @param mout	The petal measurements for the entire pupil
 * @param mphi  The petal measurements for each quadrant
 * @param npetal The number of petals.
 * @param nsa   Number of subapertures.
*/
void petal_connect(dmat **mout, dmat *mphi, int npetal, int nsa){
	if(NX(mphi)!=npetal||NY(mphi)!=nsa){
		error("mphi has invalid size\n");
	}
	dinit(mout, npetal, 1);
	dmat *diff=*mout;
	dzero(diff);
	int inan=npetal;
	for(int ip=0; ip<npetal; ip++){
		long ct=0;
		int ip1=(ip+1)%npetal;
		for(int isa=0; isa<nsa; isa++){
			real dp=P(mphi, ip1, isa)-P(mphi, ip, isa);
			if(!isnan(dp)){
				P(diff, ip1)+=dp;
				ct++;
			}
		}
		if(ct>1){
			P(diff, ip1)/=ct;
		} else if(ct==0){
			P(diff, ip1)=NAN;
			if(inan==npetal) inan=ip;
		}
	}
	//dshow(diff, "diff1");
	for(int jp=0; jp<npetal; jp++){
		int ip=(jp+inan)%npetal;//start from behind first gap
		int ip1=(ip+1)%npetal;
		if(isnan(P(diff, ip1))){
			P(diff, ip1)=0;//reset accumulation
		} else{
			P(diff, ip1)+=P(diff, ip);
		}
	}
	//dshow(diff, "diff2");
	if(inan==npetal){
		//no invalid gap, need to connect head and tail.
		//When we started we assume the value at jp+inan is 0. 
		//After accumulation, the value at jp+inan is the excessive phase that needs to be removed
		real reduce=-P(diff, 0)/(npetal);
		P(diff, 0)=0;
		for(int ip=1; ip<npetal; ip++){
			P(diff, ip)+=reduce*ip;
		}
	}
}

/**
 * Phase retrieval setup of petal measurements from PSF for a subapertures of a WFS.
 * @param[out] petal: The petal reconstruction parameters.
 * @param amp:		The pupil plane amplitude map. Sorted by subapertures.
 * @param cx,cy:    Coordinate of the pupil center.
 * @param npetal:	Number of petals
 * @param pdtheta:  Ratio of pixtheta/dtheta
 * @param pixblur: 	Pixel blur ratio
 * @param theta: 	The pupil clocking angle (CCW) in radian.
 * @param nsa:	 	Number of subapertures within the given WFS.
 * @param withtt:	Whether tip/tilt in included in the reconstruction.
*/
void petal_setup_sa(petal_t *petal, const dmat *amp, real cx, real cy, int npetal, real pdtheta, real pixblur, real theta, int nsa, int withtt){
	const long nloc=PN(amp);
	const long nx=round(sqrt(nloc));
	const long ny=nloc/nx;
	if(!petal){
		error("petal should be pre-allocated\n");
	}
	dmat *mod=NULL;
	lmat *map=NULL;
	
	dsp* hp=petal_mkh(nx, ny, cx, cy, npetal, theta);
	if(nsa>1){
		map=dspdropemptycol(hp);
	}
	dspfull(&mod, hp, 'n', 1); 
	if(nsa==1){//modify the modes for tip/tilt with opposing petals
		for(int i=0; i<npetal/2; i++){
			for(long ix=0; ix<NX(mod); ix++){
				P(mod,ix,i)-=P(mod,ix,i+npetal/2);
			}
		}
		dresize(mod, NX(mod), npetal/2);
	}
	if(withtt){
		loc_t *tt=mksqloc(nx, ny, 2./nx, 2./ny, -1, -1);
		dmat *hmod=dcat(tt->dmat, mod, 2);
		locfree(tt);
		dfree(mod);
		mod=hmod;
	}
	
	dmat *wt=dref_reshape(amp, PN(amp),1);
	petal->amp=dref_reshape(amp, nx, ny);
	petal->hpetal=hp;
	petal->mod=mod;
	petal->rmod=dpinv(mod, wt);
	petal->ind=map;
	petal->fembed=2;
	petal->npsf=12;
	long notf=petal->npsf;//no need to use the FULL psf.
	if(pdtheta || pixblur) {
		dtf_otf(&petal->otf, notf, notf, pdtheta, pdtheta, 0, pixblur, 2);
		dbg("creating otf with %ld, %g, %g\n", notf, pdtheta, pixblur);
	}
	petal->npetal=npetal;
	petal->nsa=nsa;
	petal->withtt=withtt;
	dfree(wt);
	//info("setup_sa done\n");
}

/**
 * Phase retrieval setup of petal measurements from PSFs for all subapertures of a WFS.
 * @param saloc:    Coordinate of the lower left point (where amp is defined) of each subaperture.
 * @param dx:		Grid spatial sampling 
 * @param amp:		The pupil plane amplitude map. Sorted by subapertures.
 * @param theta: 	The pupil clocking angle (CCW) in radian.
 * @param withtt:	Whether tip/tilt in included in the reconstruction.
*/
petal_t *petal_setup(const loc_t *saloc, real dx, const dmat *amp, real pdtheta, real pixblur, real theta, int withtt){
	const long nsa=(saloc&&saloc->nloc>0)?saloc->nloc:1;
	const long nloc=PN(amp)/nsa;
	const long nx=round(sqrt(nloc));
	const long ny=nloc/nx;
	petal_t *petal=mycalloc(nsa, petal_t);
	for(int isa=0; isa<nsa; isa++){
		real ox=nsa>1?saloc->locx[isa]/dx:-(nx/2-0.5);
		real oy=nsa>1?saloc->locy[isa]/dx:-(ny/2-0.5);
		//info("isa %d: ox=%g, oy=%g, dx=%g, dy=%g\n", isa, ox, oy, saloc->dx, saloc->dy);
		dmat *ampi=nsa>1?dsub(amp, isa*nloc, nloc, 0, 1):(dmat*)amp;
		petal_setup_sa(&petal[isa], ampi, -ox, -oy, 6, pdtheta, pixblur, theta, nsa, withtt); 
		if(nsa>1) dfree(ampi);
		/*writebin(petal[isa].mod, "mode_%d", isa);
		writebin(petal[isa].rmod, "rmode_%d", isa);
		writebin(petal[isa].amp, "amp_%d", isa);*/
	}
	//info("petal_setup_wfs done\n");
	return petal;
}
/**
 * Remove blurry caused by otf.
 * @param psf	PSF with peak in center
 * @param otf	detector otf with peak in corner
*/
void deblur(dmat *psf, const dmat *otf){
	cmat *wvf=cnew(NX(otf), NY(otf));
	cembedd(wvf, psf, 0);
	cfft2(wvf, -1);
	for(long i=0; i<PN(wvf); i++){
		if(P(otf, i)) P(wvf, i)/=P(otf, i);
	}
	cfft2(wvf, 1);
	if(NX(psf)==NX(otf)){
		creal2d(&psf, 0, wvf, 1./PN(wvf));//cancel FFT scaling effect.
	} else{
		dmat *tmp=NULL;
		creal2d(&tmp, 0, wvf, 1./PN(wvf));
		dembed(psf, tmp, 0);
		dfree(tmp);
	}
	cfree(wvf);
}
/**
 * Phase retrieval of petal measurements from PSF for a subaperture
 * @param[out] phi1:  The reconstructed pupil phase in zonal space.
 * @param[out] mphi1: The reconstructed pupil phase (if phi1b is specified, the update from phi1b) in modal space.
 * @param[in] petal: The petal reconstruction parameters.
 * @param ints: 	The (subaperture) PSF.
 * @param phi1b     The initial pupil phase estimate in zonal or modal (only nsa==1) space.
 * @param nrep:	 	Number of gerchberg saxton iterations.
*/
void petal_solve_sa(dmat **phi1, dmat **mphi1, const petal_t *petal, const dmat *ints, const dmat *phi1b, int nrep){
	if(!nrep) nrep=8;
	dmat *amp2=dnew(petal->npsf, petal->npsf);
	dembed(amp2, ints, 0);//pick only central region
	if(petal->otf){
		deblur(amp2, petal->otf);
	}
	if(!petal->withtt){
		//dbg_once("shift image to center\n");
		dshift2center(amp2, .5, .5);
	}
	real amax=dmax(amp2);
	if(amax==0){
		warning("amp2 is zero. cannot proceed\n"); 
		return;
	}
	for(long i=0; i<PN(amp2); i++){
		P(amp2, i)=sqrt(MAX(0,P(amp2, i))/amax);
	}
	dmat *mphi1t=NULL;
	dmat *phi1b2=NULL;
	if(phi1b){
		if(PN(phi1b)==PN(petal->amp)){
			phi1b2=(dmat*)phi1b;
		}else if(PN(phi1b)==petal->npetal){
			dspmm(&phi1b2, petal->hpetal, phi1b, "nn", 1);
		}else{
			error("Phi1b has invalid size\n");
		}
	}
	//dshow(phi1b, "phi1b");
	//dshow(phi1b2, "phi1b2");
	gerchberg_saxton(phi1, &mphi1t, petal->amp, amp2, phi1b2, petal->mod, petal->rmod, petal->fembed, nrep);
	if(mphi1){
		dinit(mphi1, petal->npetal, 1);
		if(petal->nsa>1){
			dset(*mphi1, NAN);
			for(int i=0; i<PN(petal->ind); i++){
				P(*mphi1, P(petal->ind, i))=P(mphi1t,i+(petal->withtt?2:0));//first two values in mphi1t are tip/tilt
			}
		}else{
			for(int i=0; i<petal->npetal/2; i++){
				P(*mphi1, i)=P(mphi1t,i+(petal->withtt?2:0));
				P(*mphi1, i+petal->npetal/2)=-P(*mphi1, i);
			}
		}
	}
	dfree(mphi1t);
	dfree(amp2);
	if(phi1b2!=phi1b) dfree(phi1b2);
}
/**
 * Phase retrieval setup and reconstruction of petal measurements from PSFs for all subapertures of a WFS.
 * @param[out] phi1  The reconstructed pupil phase in zonal space, sorted by subapertures.
 * @param[out] mphi1 The reconstructed pupil phase (if zonal phi1b is specified, the update from phi1b) in modal space for the entire pupil.
 * @param petal     The petal reconstruction parameters.
 * @param ints      The (subaperture) PSFs sorted by subapertures.
 * @param phi1b     The initial pupil phase estimate in zonal or modal (only nsa==1) space.
 * @param nrep      Number of gerchberg saxton iterations.
*/
void petal_solve(dcell **phi1, dmat **mphi1, const petal_t *petal, const dcell *ints, const dmat *phi1b, int nrep){
	if(!petal || !ints){
		warning("Unable to proceed: petal=%p, ints=%p\n", petal, ints);
		return;
	}
	const long nloc=PN(petal->amp);
	const long nsa=petal->nsa;
	if(phi1b && ((nsa>1 && PN(phi1b)!=nloc*nsa) || (nsa==1 && PN(phi1b)!=petal->npetal))){
		error("phi1b should have the same size as amp or is npetal for tt wfs.\n");
	}
	if(phi1) cellinit(phi1, nsa, 1);
	dmat *mval=dnew(petal->npetal, nsa);
	//if(mphi1) dinit(mphi1, petal->npetal, nsa>1?(nsa+1):1);
	//info("solve_sa wfs. amp is %ld, nsa is %ld, nloc is %ld.\n", PN(amp), nsa, nloc);
	for(long isa=0; isa<nsa; isa++){
		dmat *phi1bi=phi1b?(nsa>1?dsub(phi1b, isa*nloc, nloc, 0, 1):(dmat*)phi1b):NULL;
		dmat *mvali=nsa>1?drefcols(mval, isa, 1):mval;
		petal_solve_sa(phi1?&P(*phi1,isa):NULL, &mvali, &petal[isa], P(ints, isa), phi1bi, nrep);
		if(nsa>1){
			dfree(phi1bi); 
			dfree(mvali);
		}
	}
	if(mphi1){
		if(nsa>1){
			petal_connect(mphi1, mval, petal->npetal, nsa);
		}else{
			if(PN(phi1b)==petal->npetal){
				//don't add phi1b to mphi1. They may be the same memory.
				//For TT WFS, if modal phi1b is provided, add it to mphi1.
				dadd(&mval, 1, phi1b, 1);
			}
			dcp(mphi1, mval);
		}
	}
	dfree(mval);
}
void petal_save(petal_t *petal, const char *format, ...){
	format2fn;
	if(!petal) return;
	for(int isa=0; isa<petal->nsa; isa++){
		writebin(petal[isa].hpetal, "%s_hpetal_%d", fn, isa);
		writebin(petal[isa].mod, "%s_mode_%d", fn, isa);
		writebin(petal[isa].amp, "%s_amp_%d", fn, isa);
		writebin(petal[isa].rmod, "%s_rmod_%d", fn, isa);
	}
}
/**
 * A convenient wrapper for petal_setup and petal_solve to be called by Python.
 * Phase retrieval setup and reconstruction of petal measurements from PSFs for all subapertures of a WFS. 
 * @param[out] phi1  The reconstructed pupil phase in zonal space, sorted by subapertures.
 * @param[out] mphi1 The reconstructed pupil phase (if phi1b is specified, the update from phi1b) in modal space for the entire pupil.
 * @param ints      The (subaperture) PSFs sorted by subapertures.
 * @param phi1b     The initial pupil phase estimate.
 * @param saloc     Coordinate of the lower left point (where amp is defined) of each subaperture.
 * @param amp       The pupil plane amplitude map. Sorted by subapertures.
 * @param theta    The pupil clocking angle (CCW) in radian.
 * @param nrep      Number of gerchberg saxton iterations.
 * @param withtt    Whether tip/tilt in included in the reconstruction.
*/
void petal_solve_wfs(dcell **phi1, dmat **mphi1, const dcell *ints, const dmat *phi1b, const loc_t *saloc, const dmat *amp, real pdtheta, real pixblur, real theta, int nrep, int withtt){
	petal_t *petal=petal_setup(saloc, 1, amp, pdtheta, pixblur, theta, withtt);
	if(petal) petal_solve(phi1, mphi1, petal, ints, phi1b, nrep);
	if(petal) petal_free(petal, petal->nsa);
}
