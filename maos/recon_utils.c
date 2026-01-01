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
#include "recon_utils.h"
#include "recon.h"
#include "ahst.h"
#include "sim.h"

/**
 * \file recon_utils.h
   Reusable utilities for wavefront reconstruction and DM fitting, supporting recon.c
*/

/**
   Apply Laplacian2 to xin and accumulates to xout.
*/
void apply_L2(dcell** xout, const dspcell* L2, const dcell* xin,
	real alpha){
	dcell* xx=NULL;
	dcellmm(&xx, L2, xin, "nn", 1.);
	dcellmm(xout, L2, xx, "tn", alpha);
	dcellfree(xx);
}
/**
   Apply turbulence invpsd to xin in Fourier space, scaled by alpha and add to xout.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks.
   \todo: this is not thread safe.
*/
void apply_invpsd(dcell** xout, const void* A, const dcell* xin, real alpha, int xb, int yb){
	if(xb!=yb) return;
	const invpsd_t* extra=(const invpsd_t*)A;
	dcell* invpsd=extra->invpsd;
	ccell* fftxopd=extra->fftxopd;
	int ips1, ips2;
	if(xb<0){/*do all cells */
		ips1=0;
		ips2=NX(xin)*NY(xin);
	} else{/*do the specified cell */
		ips1=xb;
		ips2=xb+1;
	}
	for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<NX(xin)*NY(xin); ips++){ */
		long nx=P(fftxopd, ips)->nx;
		long ny=P(fftxopd, ips)->ny;
		if(extra->square){
			dmat* xini=dref_reshape(P(xin, ips), nx, ny);
			ccpd(&P(fftxopd, ips), xini);
			dfree(xini);
		} else{
			czero(P(fftxopd, ips));
			cembed_locstat(&P(fftxopd, ips), 0, P(extra->xloc, ips), P(P(xin, ips)), alpha, 0);
		}
		cfft2(P(fftxopd, ips), -1);
		ccwmd(P(fftxopd, ips), P(invpsd, ips));
		cfft2(P(fftxopd, ips), 1);
		if(extra->square){
			dmat* xouti=NULL;
			xouti=dref_reshape(P(*xout, ips), nx, ny);
			creal2d(&xouti, 1, P(fftxopd, ips), alpha);
			dfree(xouti);
		} else{
			cembed_locstat(&P(fftxopd, ips), 1, P(extra->xloc, ips), P(P(*xout, ips)), 1, 1);
		}
	}
}

/**
   Apply fractal regularization to x, scaled by alpha.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks.
   \todo: this is not thread safe.
*/
void apply_fractal(dcell** xout, const void* A, const dcell* xin, real alpha, int xb, int yb){
	if(xb!=yb) return;
	const fractal_t* extra=(const fractal_t*)A;
	int ips1, ips2;
	if(xb<0){/*do all cells */
		ips1=0;
		ips2=NX(xin)*NY(xin);
	} else{/*do the specified cell */
		ips1=xb;
		ips2=xb+1;
	}
	for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<NX(xin)*NY(xin); ips++){ */
		dzero(P(extra->xopd, ips));
		real r0i=extra->r0*pow(extra->wt[ips], -3./5.);
		dembed_locstat(&P(extra->xopd, ips), 0, P(extra->xloc, ips), P(P(xin, ips)),
			alpha*extra->scale, 0);
		fractal_inv(P(extra->xopd, ips),
			P(extra->xloc, ips)->dx, r0i, extra->L0, extra->ninit);
		fractal_inv_trans(P(extra->xopd, ips),
			P(extra->xloc, ips)->dx, r0i, extra->L0, extra->ninit);
		dembed_locstat(&P(extra->xopd, ips), 1, P(extra->xloc, ips), P(P(*xout, ips)), 1, 1);
	}
}
/**
 * @brief Remove projection
 * x = x - M * (PM * x)
 * @param x 	The vector
 * @param M 	The mode
 * @param PM 	The mode projection matrix
 */
void remove_mode(dcell *x, const dcell *M, const dcell *PM){
	if(!x || !M || !PM) return;
	dcell* junk=NULL;
	dcellmm(&junk, PM, x, "nn", 1);
	dcellmm(&x, M, junk, "nn", -1);
	dcellfree(junk);
}

/**
   Apply weighting W0/W1 to a vector. \f$W0*x-W1*(W1'*x)\f$
*/
static dmat* applyWeach(dmat* xin, const dsp* W0, const dmat* W1, const real wt){
	if(!W0||!W1){
		warning("W0 or W1 is NULL\n");
		return dref(xin);
	}
	dmat* xout=NULL;//do not use out. it maybe the same as xin
	dmat* tmp=NULL;
	dspmm(&xout, W0, xin, "nn", wt);
	dmm(&tmp, 0, W1, xin, "tn", -1);
	dmm(&xout, 1, W1, tmp, "nn", wt);
	dfree(tmp);
	return xout;
}
/**
   apply weighting W0/W1 with weighting wt for each block.
   \f$xin<-(W0-W1*W1')*xin;\f$
*/
void applyW(dcell* xin, const dsp* W0, const dmat* W1, const real* wt){
	const int nevl=NX(xin);
	for(int iy=0; iy<NY(xin); iy++){
		for(int ievl=0; ievl<nevl; ievl++){
			dmat *xout=applyWeach(P(xin, ievl, iy), W0, W1, wt[ievl]);
			dcp(&P(xin, ievl, iy), xout);
			dfree(xout);
		}
	}
}

/**
   Compute W0/W1 weighted dot product: \f$A'*(W0-W1*W1')*B\f$
*/
dcell* calcWmcc(const dcell* A, const dcell* B, const dsp* W0,
	const dmat* W1, const dmat* wt){

	assert(NX(wt)==NX(B)&&NY(wt)==1&&NX(A)==NX(B));
	const int nevl=NX(B);
	dcell* res=dcellnew(NY(A), NY(B));
	OMP_FOR(4)
	for(int iy=0; iy<NY(B); iy++){
		for(int ievl=0; ievl<nevl; ievl++){
			dmat* xout=applyWeach(P(B, ievl, iy), W0, W1, P(wt, ievl));
			OMP_FOR(4)
			for(int ix=0; ix<NY(A); ix++){
				dmm(&P(res, ix, iy), 1, P(A, ievl, ix), xout, "tn", 1);
			}
			dfree(xout);
		}
	}
	return res;
}


/**
 * Convert column vector nsax2 or nsax3 to sparse matrix. The third column of nea is coupling.
 * ll: controls lower left
 * ur: controls upper right 
 * if ng>2, treat as diagonal matrix (for raw PWFS)
 */
dsp* nea2sp(dmat* nea, int ll, int ur, int ng){
	if(NY(nea) != ng && NY(nea)!=3){
		error("nea has wrong format:%ldx%ld. Expect 2 or 3 columns.\n", NX(nea), NY(nea));
	}
	if(NY(nea)==ng){//off-diagonal not available
		ll=0;
		ur=0;
	}
	const long nsa=NX(nea);
	dsp* sanea=dspnew(nsa*ng, nsa*ng, (ng+(ll?1:0)+(ur?1:0))*nsa);
	spint* pp=sanea->pp;
	spint* pi=sanea->pi;
	real* px=sanea->px;
	long count=0;
	for(long ig=0; ig<ng; ig++){
		for(long isa=0; isa<nsa; isa++){
			pp[isa+nsa*ig]=count;
			if(ig==1&&ur&&P(nea, isa, 2)){/*Cxy */
				pi[count]=isa;
				px[count]=P(nea, isa, 2);
				count++;
			}
			{/*Cxx */
				pi[count]=isa+nsa*ig;
				px[count]=P(nea, isa, ig);
				count++;
			}
			if(ig==0&&ll&&P(nea, isa, 2)){/*Cyx */
				pi[count]=isa+nsa;
				px[count]=P(nea, isa, 2);
				count++;
			}
		}
	}
	pp[nsa*ng]=count;
	if(count>sanea->nzmax){
		error("memory overflow\n");
	}else if(count<sanea->nzmax){
		dspsetnzmax(sanea, count);
	}
	return sanea;
}
/**
   Check the dimension of NEA and fix if possible
*/
void nea_check(dmat *nea, int nsa, int ng){
	if(!nea){
		error("nea is not defined.\n");
	} else if(NY(nea)==1){
		reshape(nea, nsa, NX(nea)/nsa);
	}
	if(NX(nea)!=nsa||(NY(nea)!=ng&&NY(nea)!=3)){
		error("nea has wrong format (%ldx%ld), should be (%dx%d).\n", NX(nea), NY(nea), nsa, 3);
	}
}
/**
 * Common routine for the following to prepare the data.
 * */
static int nea_prep(dmat **pout, const dmat *in, const int ng){
	if(NY(in)!=ng&&NY(in)!=3){
		error("in has wrong format: %ldx%ld\n", NX(in), NY(in));
	}
	int isxy=(NY(in)>ng&&ng==2)?1:0;//cross term
	int ncol=ng==2?3:ng;
	if(!*pout){
		*pout=dnew(NX(in), ncol);
	}else if(isxy&&NY((*pout))!=ncol){
		dresize(*pout, NX(in), ncol);
	}
	return isxy;
}
/**
   Apply cholesky in a 2x2 symmetric matrix packed in [a[0,0],a[1,1],a[0,1]] as row vector.
   Input and output may be the same.
   if ng>2, treat as diagonal matrix
 */
void nea_chol(dmat **pout, const dmat *in, const int ng){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
		//Use temporary variable to handle the case that out and in is the same.
			real a=sqrt(P(in, isa, 0));
			real b=(isxy && a)?(P(in, isa, 2)/a):0;
			real c=sqrt(P(in, isa, 1)-b*b);
			P(out, isa, 0)=a;
			P(out, isa, 1)=c;
			if(isxy) P(out, isa, 2)=b;
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			P(out, ig)=sqrt(P(in, ig));
		}
	}
}
/**
   Apply LL' to lower diagonal matrix packed in [a[0,0],a[1,1],a[0,1]] as row vector.
   Input and output may be the same.
   if ng>2, treat as diagonal matrix
 */
void nea_mm(dmat **pout, const dmat *in, const int ng){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
		//Use temporary variable to handle the case that out and in is the same.
			real a=P(in, isa, 0);
			real b=isxy?(P(in, isa, 2)):0;
			real c=P(in, isa, 1);
			P(out, isa, 0)=a*a;
			P(out, isa, 1)=c*c+b*b;
			if(isxy) P(out, isa, 2)=a*b;
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			P(out, ig)=P(in, ig)*P(in, ig);
		}
	}
}
/**
   Apply matrix inversion in a 2x2 symmetrisa matrix packed in row vector [a[0,0],a[1,1],a[0,1]]
   Input and output may be the same.
*/

void nea_inv(dmat **pout, const dmat *in, int ng, real scale){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
			//use double to prevent overflow of intermediate value
			double xx=P(in, isa, 0);
			double yy=P(in, isa, 1);
			double xy=isxy?P(in, isa, 2):0;
			double invdet=scale/(xx*yy-xy*xy);
			if(isfinite(invdet)){
				P(out, isa, 0)=invdet*yy;
				P(out, isa, 1)=invdet*xx;
				if(xy) P(out, isa, 2)=-invdet*xy;
			}
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			if(P(in, ig)){
				P(out, ig)=scale/P(in, ig);
			}
		}
	}
}

/**
   Compute and save PSF reconstruction telemetry.  For pseudo open loop
   estimations, like high order loop, we add opdr to, and subtract dmpsol from
   the OPD.  For closed loop estimations, like ahst low order or lsr, we add
   dmerr_lo, and dmerr to the OPD.*/
void psfr_calc(sim_t* simu, dcell* opdr, dcell* dmpsol, dcell* dmerr, dcell* dmerr_lo){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	/* The tomography estimates, opdr is pseudo open loop estimates. We need to
	   subtract the constribution of the added DM command to form closed loop
	   estimates.
	*/
	dcell* dmtmp=NULL;
	if(opdr&&parms->dbg.useopdr){
	/* The original formulation using Hx*x-Ha*a.  Changed from ploc to plocs
	   on July 18, 2011. Ploc is too low sampled. Make sure the sampling of
	   plocs is not too big.  Deprecated. Use dm space instead.
	*/
		if(dmpsol){/*Pseudo OL estimates */
			if(parms->recon.split==1){
			/* We will remove NGS modes from dmlast which is in NULL modes of
			   tomography reconstructor (is this 100% true)?  SHould we remove NGS
			   modes from final OPD, xx, instead?*/
				dcell* tmp=dcelldup(dmpsol);/*The DM command used for high order. */
				ngsmod_remove(simu, tmp);/*remove NGS modes as we do in ahst. */
				dcelladd(&dmtmp, 1, tmp, -1);
				dcellfree(tmp);
			} else{
				dcelladd(&dmtmp, 1, dmpsol, -1);
			}
		}

		loc_t* locs=simu->aper->locs;
		dmat* xx=dnew(locs->nloc, 1);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			real hs=P(parms->evl.hs, ievl);
			if(P(parms->evl.psfr, ievl)){
				dzero(xx);
				if(opdr){
					const int npsr=recon->npsr;
					/*First compute residual opd: Hx*x-Ha*a*/
					for(int ips=0; ips<npsr; ips++){
						propdata_t propdata={.phiin=P(P(opdr, ips)), .locout=locs, .phiout=P(xx), .alpha=1, .hs=hs, .thetax=P(parms->evl.thetax, ievl), .thetay=P(parms->evl.thetay, ievl)};
						if(parms->tomo.square){/*square xloc */
							propdata.mapin=P(recon->xmap, ips);
						} else{
							propdata.locin=P(recon->xloc, ips);
						}
						prop(&propdata);
					}
				}
				if(dmtmp){
					for(int idm=0; idm<parms->ndm; idm++){
						prop(&(propdata_t){.locin=P(recon->aloc, idm), .phiin=P(P(dmtmp, idm)), .locout=locs,
							.phiout=P(xx), .alpha=1, .hs=hs, .thetax=P(parms->evl.thetax, ievl), .thetay=P(parms->evl.thetay, ievl)});
					}
				}
				dmm(&P(simu->ecov, ievl), 1, xx, xx, "nt", 1);
				if(parms->dbg.ecovxx){
					zfarr_push(simu->save->ecovxx[ievl], simu->reconisim, xx);
				}
			}/*if psfr[ievl] */
		}/*ievl */
		dfree(xx);
	} else{/*Do Ha in postproc, so just do a. */
		if(dmerr){/*high order closed loop estimates. (lsr)*/
			dcelladd(&dmtmp, 1, dmerr, 1);
		}
		if(dmerr_lo){/*In AHST, dmerr_lo is CL Estimation.*/
			addlow2dm(&dmtmp, simu, dmerr_lo, 1);
		}
		dcellmm(&simu->ecov, dmtmp, dmtmp, "nt", 1);
	}
	dcellfree(dmtmp);
}


/**
   Prepare arrays for cn2 estimation. Multiple pairs can be used to do Cn2
Estimation. The result will be an average of them.  */


cn2est_t* cn2est_prepare(const parms_t* parms, const powfs_t* powfs){
	dmat* pair=parms->cn2.pair;
	int npair=NX(pair)*NY(pair);
	int ipowfs=-1;
	for(int ind=0; ind<npair; ind++){
		int iwfs=(int)P(pair, ind);
		if(ipowfs==-1){
			ipowfs=parms->wfs[iwfs].powfs;
		} else if(ipowfs!=parms->wfs[iwfs].powfs){
			error("All wfs in parms->cn2.pair do not belong to the same powfs\n");
		}
	}
	dmat* wfstheta=dnew(parms->nwfs, 2);
	dmat* ptheta=wfstheta/*PDMAT*/;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		P(ptheta, iwfs, 0)=parms->wfs[iwfs].thetax;
		P(ptheta, iwfs, 1)=parms->wfs[iwfs].thetay;
	}
	dmat* ht=0;
	if(parms->cn2.keepht){
		ht=dnew(parms->atmr.nps, 1);
		for(int iht=0; iht<NX(ht); iht++){
			P(ht, iht)=P(parms->atmr.ht, iht);
		}
	} else{
		int nht=parms->cn2.nhtomo;
		ht=dnew(nht, 1);
		real dht=parms->cn2.hmax/(real)(nht-1);
		if(dht<=0) error("Invalid hmax or nhtomo\n");
		for(int iht=0; iht<nht; iht++){
			P(ht, iht)=dht*iht;
		}
	}
	dmat* hs=dnew(parms->nwfs, 1);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		P(hs, iwfs)=parms->wfs[iwfs].hs;
	}
	cn2est_t* cn2est=cn2est_new(pair, wfstheta, powfs[ipowfs].saloc, powfs[ipowfs].saamax, parms->cn2.saat,
		hs, ht, parms->cn2.keepht, parms->atmr.L0);
	dfree(hs);
	dfree(wfstheta);
	if(cn2est){
		cn2est->os=dnew(NX(ht), 1);
		if(!parms->cn2.keepht){
			/*preserve the number of over sampled layers. */
			int osc=0;
			for(int ips=0; ips<parms->atmr.nps; ips++){
				if(P(parms->atmr.os, ips)>1){
					osc++;
					if(P(parms->atmr.os, ips)!=2){
						error("os is not 2. adept this code to it.\n");
					}
				}
			}
			for(int iht=0; iht<NX(ht); iht++){
				if(iht<osc){
					P(cn2est->os, iht)=2;
				} else{
					P(cn2est->os, iht)=1;
				}
			}
		} else{
			for(int iht=0; iht<NX(ht); iht++){
				P(cn2est->os, iht)=P(parms->atmr.os, iht);
			}
		}
		if(parms->cn2.verbose){
			cn2est->dmht=dnew(parms->ndm, 1);
			for(int idm=0; idm<parms->ndm; idm++){
				P(cn2est->dmht, idm)=parms->dm[idm].ht;
			}
		}
		if(parms->save.setup>1){
			writebin(cn2est->overlapi, "cn2_overlapi");
			writebin(cn2est->iPnk, "cn2_iPnk");
			writebin(cn2est->Pnk, "cn2_Pnk");
			writebin(cn2est->ht, "cn2_ht");
			writebin(cn2est->wtconvert, "cn2_wtconvert");
		}
	}
	dfree(ht);
	return cn2est;
}
/**
   Implemented mechanism to move height of layers.
 */
static void cn2est_moveht(recon_t* recon){
	(void)recon;
	/*cn2est_t *cn2est=recon->cn2est; */
	/*
	  Implemented mechanism to move height of layers. Need to redo HXF, GX, etc.
	*/
	error("moveht not yet implemented");
}

/**
   Wrapper of Cn2 Estimation operations in recon.c
*/
void cn2est_isim(dcell* cn2res, recon_t* recon, const parms_t* parms, const dcell* grad, int* tomo_update){
	cn2est_t* cn2est=recon->cn2est;
	cn2est_push(cn2est, grad);
	static int icn2=-1;
	if(cn2est->count%parms->cn2.step==0){
		icn2++;
		int nset=cn2est->count/parms->cn2.step;
		cn2est_est(cn2est, parms->cn2.verbose);/*do the CN2 estimation */
		if(parms->cn2.reset&&(nset%parms->cn2.reset)==0){
			cn2est_reset(cn2est);
		}

		if(cn2res){
			P(P(cn2res, 0), icn2)=cn2est->r0m;
			memcpy(PCOL(P(cn2res, 1), icn2), P(P(cn2est->wtrecon, 0)),
				NX(cn2est->htrecon)*sizeof(real));
		}
		if(parms->cn2.tomo){
			if(parms->cn2.moveht){
				cn2est_moveht(recon);
			}

			if(parms->cn2.verbose){
				info2("Updating tomography weights\n");
			}
			/*Changes recon parameters. cannot be parallel with tomofit(). */
			/*wtrecon is referenced so should be updated automaticaly. */
			if(P(recon->wt)!=P(P(cn2est->wtrecon, 0))){
				dfree(recon->wt);
				recon->wt=dref(P(cn2est->wtrecon, 0));
			}
			recon->r0=cn2est->r0m;
			recon->L0=cn2est->L0;

			*tomo_update=1;
		}
	}
}
