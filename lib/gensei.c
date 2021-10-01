/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
  2009-11-26: changed to rotate OTF instead of psf to comply with the approach
  in wfsints.  this gives slightly larger intensity because max value of OTF is
  preserved which corresponds to the sum of psf.


*/

#include "gensei.h"
#include "genotf.h"
/**
   see genseotf. Check whether cells are different.
*/
static int count_unique(const dcell* opdbias){
	if(!opdbias) return 1;
	int different=0;
	for(int iwfs=1; iwfs<PN(opdbias); iwfs++){
		real diff=ddiff(P(opdbias, 0), P(opdbias, iwfs));
		if(diff>1e-4){
			dbg("opdbias[%d] is different from opdbias[0] by %g.\n", iwfs, diff);
			different=1;
		}
	}
	return different?PN(opdbias):1;
}
static cccell* genseotf_do(const pts_t* pts,
	const cell* amp, const dcell* opdbias,
	const cell* saa, const dmat* wvl, real r0, real L0, int embfac){

	/*create a grid representing a sub-aperture. */
	loc_t* loc=mksqloc_auto(pts->nxsa, pts->nxsa, pts->dx, pts->dy);
	/*The embeding factor for embedding the aperture */
	const int npsfx=pts->nxsa*embfac;
	const int npsfy=pts->nysa*embfac;
	const int nwvl=PN(wvl);
	const int nsa=pts->nsa;

	int notf=count_unique(opdbias);
	int notf2=iscell(amp)?count_unique((dcell*)amp):1;
	if(notf==1&&notf2!=1){
		notf=notf2;
	} else if(notf2!=1&&notf!=notf2){
		error("Mismatch for notf: %d vs %d\n", notf, notf2);
	}
	cccell* otf=cccellnew(notf, nwvl);

	info("There is %s bias\n", opdbias?"NCPA":"no");
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		for(int iotf=0; iotf<notf; iotf++){
			dmat* opdi=opdbias?P(opdbias, iotf):NULL;
			real thres=opdi?1:(1-1e-10);
			const dmat* ampi=iscell(amp)?PR((dcell*)amp, iotf, 1):(dmat*)amp;
			const dmat* saai=iscell(saa)?PR((dcell*)saa, iotf, 1):(dmat*)saa;
			//OTFs are always generated with native sampling. It is upsampled at gensepsf if necessary.
			OMPTASK_SINGLE
				genotf(&P(otf, iotf, iwvl), loc, ampi, opdi, saai,
					thres, P(wvl, iwvl), NULL, r0, L0, npsfx, npsfy, nsa, 1);
		}
	}/*iwvl */
	locfree(loc);
	return otf;
}
/**
   Generates short exposure OTF by calling genotf() with p/t/t removal set.
*/
cccell* genseotf(const pts_t* pts, /**<[in]subaperture low left coordinate*/
				const cell* amp, /**<[in] amplitude map. can be dcell or dmat*/
				const dcell* opdbias, /**<[in] opd bias for otf.*/
				const cell* saa, /**<[in] list of subaperture area, can be dcell or dmat*/
				const dmat* wvl, /**<[in] list of wavelength*/
				const real r0, /**<[in] Fried parameter*/
				const real L0, /**<[in] Outer scale*/
				const int embfac/**<[in] Embedding factor, normally 2*/
){
	char fnprefix[200]; fnprefix[0]='\0';
	uint32_t key=0;
	strcat(fnprefix, "SEOTF");
	if(amp){
		if(iscell(amp)){
			key=dcellhash((dcell*)amp, key);
		} else{
			key=dhash((dmat*)amp, key);
		}
	}
	if(wvl){
		key=dhash(wvl, key);
	}
	if(opdbias){
		key=dcellhash(opdbias, key);
	}
	if(key!=0){
		char tmp2[80];
		snprintf(tmp2, 80, "_%ud", key);
		strcat(fnprefix, tmp2);
	}
	char fnotf[PATH_MAX+20];
	char fnlock[PATH_MAX+40];
	snprintf(fnotf, PATH_MAX, "%s/SEOTF/", CACHE);
	if(!exist(fnotf)){
		mymkdir("%s", fnotf);
	}

	snprintf(fnotf, sizeof(fnotf), "%s/SEOTF/%s_r0_%g_L0%g_dsa%g_nsa%ld_dx1_%g_embfac%d_v3",
		CACHE, fnprefix, r0, L0, pts->dsa, pts->nsa, 1./pts->dx, embfac);
	snprintf(fnlock, sizeof(fnlock), "%s.lock", fnotf);
	cccell* otf=0;
	while(!otf){
		if(exist(fnlock)||!zfexist("%s", fnotf)){/*need to create data */
			int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
			if(fd>=0){/*succeed */
				info("Generating OTF for %s...", fnotf);
				TIC;tic;
				otf=genseotf_do(pts, amp, opdbias, saa, wvl, r0, L0, embfac);
				toc2("done");
				writebin(otf, "%s", fnotf);
			} else{
				warning("Waiting for previous lock to release ...");
				fd=lock_file(fnlock, 1, 0);
			}
			close(fd); remove(fnlock);
		} else{
			info("Reading OTF from %s\n", fnotf);
			otf=cccellread("%s", fnotf);
		}
	}
	return otf;
}
/**
   Upsample the otf in to out while preserving the PSF.
 */
void upsample_otf(cmat* out, const cmat* in){
	if(NX(in)==NX(out)&&NY(in)==NY(out)){
		ccp(&out, in);
	} else{
		cmat* temp=0;
		ccp(&temp, in);
		cfft2(temp, -1);
		cscale(temp, 1./(NX(in)*NY(in)));
		czero(out);
		ccpcorner(out, temp, C_FULL);
		cfft2(out, 1);
		cfree(temp);
	}
}
/**
   Createing subaperture short exposure PSF from the tip/tilt removed turbulence
   OTF and uplink OTF. Not including detector or elongation characteristics.  */
void gensepsf(dccell** psepsfs, const cccell* otfs, const cccell* lotf,
	const cell* saa, const dmat* wvl, int notfx, int notfy
){
	const int nwvl=PN(wvl);
	const int nsa=NX(P(otfs, 0, 0));
	const int nlotf=lotf?NX(lotf):0;
	int nsepsf=NX(otfs);
	if(nlotf>1){
		if(nsepsf==1){
			nsepsf=nlotf;
		} else if(nsepsf!=nlotf){
			error("mismatch: notf is %d, lotf is %d\n", nsepsf, nlotf);
		}
	}
	dccell* sepsfs=NULL;
	cellinit((cell**)psepsfs, nsepsf, nwvl);

	if(!*psepsfs){
		*psepsfs=dccellnew(nsepsf, nwvl);
	}
	sepsfs=*psepsfs;
	if(!notfx) notfx=NX(P(P(otfs, 0), 0));
	if(!notfy) notfy=NX(P(P(otfs, 0), 0));
	for(int isepsf=0; isepsf<nsepsf; isepsf++){
		const dmat* saai=saa?(iscell(saa)?PR((dcell*)saa, isepsf, 1):(dmat*)saa):NULL;
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			P(sepsfs, isepsf, iwvl)=dcellnew_same(nsa, 1, notfx, notfy);
			const ccell* otf=PR(otfs, isepsf, iwvl);
			cmat* sepsf=cnew(notfx, notfy);
			cmat* lotfi=0;
			if(nlotf){
				cmat* lotfi2=P(PR(lotf, isepsf, iwvl), 0);
				if(NX(lotfi2)!=notfx||NY(lotfi2)!=notfy){
					lotfi=cnew(notfx, notfy);
					upsample_otf(lotfi, lotfi2);
				} else{
					lotfi=cref(lotfi2);
				}
			}
			for(int isa=0; isa<nsa; isa++){
				real norm=(saai?P(saai, isa):1)/((real)(notfx*notfy));
				if(P(otf, isa, 0)){
					upsample_otf(sepsf, P(otf, isa, 0));/*peak in center */
				} else{
					czero(sepsf);
				}
				if(lotfi){
					ccwm(sepsf, lotfi);
				}
				cfftshift(sepsf); /*peak now in corner. */
				cfft2(sepsf, 1);   /*turn to psf. FFT 1th */
				cfftshift(sepsf); /*psf with peak in center */
				creal2d(&P(P(sepsfs, isepsf, iwvl), isa, 0), 0, sepsf, norm);/*copy to output. */
			}
			cfree(sepsf);
			cfree(lotfi);
		}
	}
}
/**
   generate subaperture short exposure average pixel intensities sampled on
   detector from short expsoure PSF, the elongation transfer function of the
   sodium layer, and the detector transfer function. */
void gensei(dcell** pi0, dcell** pgx, dcell** pgy, cccell** pfotf,
	const dccell* sepsfs, const dtf_t* dtf, const etf_t* etf, const dcell* saa,
	const dcell* gxyrot, const dmat* siglev, const dmat* wvlwts, const dcell* goff,
	int i0scale, int shift2center
){
	if(!sepsfs||!dtf){
		error("sepsfs and dtf must be set\n");
	}
	if(shift2center&&goff){
		error("cannot specify both goff and shift2center");
	}
	const int notfx=NX(P(P(sepsfs, 0, 0), 0, 0));
	const int notfy=NY(P(P(sepsfs, 0, 0), 0, 0));
	const int nwvl=NY(sepsfs);
	const int nsa=NX(P(sepsfs, 0, 0));
	const int nllt=etf?NY(etf[0].etf):0;
	if(etf&&nsa!=NX(etf[0].etf)){
		error("mismatch: nsa=%d, etf is %ldx%ld\n", nsa, NX(etf[0].etf), NY(etf[0].etf));
	}

	/**
	   ni0 may be greater than 1 in the following two cases
	   1) multiple LLT
	   2) different signal level or wvlwts
	   3) powfs[ipowfs].bkgrnd contains rayleigh scatter bkgrnd for each wfs in this powfs.
	*/
	int ni0=MAX(NX(sepsfs), nllt);
	if(NY(wvlwts)>1){
		if(ni0==1){
			ni0=NY(wvlwts);
		} else if(ni0!=NY(wvlwts)){
			error("Mismatch: ni0=%d, wvlwts is %ldx%ld\n", ni0, NX(wvlwts), NY(wvlwts));
		}
	}
	if(PN(siglev)>1){
		if(ni0==1){
			ni0=PN(siglev);
		} else if(ni0!=PN(siglev)){
			error("Mismatch: ni0=%d, siglev is %ldx1\n", ni0, PN(siglev));
		}
	}
	if(goff&&PN(goff)>1){
		if(ni0==1){
			ni0=PN(goff);
		} else if(ni0!=PN(goff)){
			error("Mismatch: ni0=%d, goff is %ldx1\n", ni0, NX(goff));
		}
	}
	//dbg("number of i0 is %d\n", ni0);

	const int pixpsax=dtf[0].pixpsax;
	const int pixpsay=dtf[0].pixpsay;
	const real dtheta=dtf[0].dtheta;
	if(pi0){
		if(!*pi0 || NX(*pi0)!=nsa||NY(*pi0)!=ni0){
			dcellfree(*pi0);
			*pi0=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		} else{
			dcellzero(*pi0);
		}
	}
	if(pgx){
		if(!*pgx || NX(*pgx)!=nsa||NY(*pgx)!=ni0){
			dcellfree(*pgx);
			*pgx=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		} else{
			dcellzero(*pgx);
		}
	}
	if(pgy){
		if(!*pgy || NX(*pgy)!=nsa||NY(*pgy)!=ni0){
			dcellfree(*pgy);
			*pgy=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		} else{
			dcellzero(*pgy);
		}
	}
	if(pfotf){
		cellfree(*pfotf);
		*pfotf=cccellnew(ni0, 1);
		for(int i=0; i<ni0; i++){
			P(*pfotf, i)=ccellnew(nsa, nwvl);
		}
	}
	dcell* i0=pi0?*pi0:0;
	dcell* gx=pgx?*pgx:0;
	dcell* gy=pgy?*pgy:0;

	/*
	  Notice, the generation of shifted i0s are not accurate
	  because the PSF is not enough to cover the size.
	  Disable the computation.
	*/

	if(i0scale){
		warning("i0 is scaled to match sa area\n");
	}

	for(int iwvl=0; iwvl<nwvl; iwvl++){
		const comp* Ux=P(dtf[iwvl].Ux);
		const comp* Uy=P(dtf[iwvl].Uy);
		const real norm=1./(real)(notfx*notfy);
		const ccell* petf=etf?etf[iwvl].etf:0;
		for(int ii0=0; ii0<ni0; ii0++){
			real* area=P(PR(saa, ii0, 0));
			real wvlsig=PR(wvlwts, iwvl, ii0)*PR(siglev, ii0, 0);

			dcell* psepsf=PR(sepsfs, ii0, iwvl);
			real* angles=(gxyrot)?P(PR(gxyrot, ii0, 0)):0;
			ccell* se_cache=ccellnew_same(2, MAXTHREAD, notfx, notfy);

			OMP_FOR
			for(int isa=0; isa<nsa; isa++){
				int ith=0;
				/*loaded psepsf. sum to 1 for full sa. peak in center */
				dmat* sepsfi=P(psepsf, isa, 0);
#if DEBUG				
				if(dsum(sepsfi)>1.1){
					error("Short exposure PSF has wrong scaling. It should total to <=1\n");
				}
#endif				
#ifdef _OPENMP
				ith=omp_get_thread_num();
#endif
#define seotfj P(se_cache,0,ith)
#define seotfk P(se_cache,1,ith)
				/*if(!seotfk){
					seotfk=cnew(notfx, notfy);
				}*/
				cmat* nominal=NULL;
				dsp* si=NULL;
				if(!etf){
					nominal=PR(dtf[iwvl].nominal, isa, ii0);
				}
				si=PR(dtf[iwvl].si, isa, ii0);
				real pgrad[2]={0,0};

				if(goff){//convert goff in radian to sepsf pixel unit.
					pgrad[0]=-P(PR(goff, ii0, 0), isa)/dtheta;
					pgrad[1]=-P(PR(goff, ii0, 0), isa+nsa)/dtheta;
				} else if(shift2center){
					/*Forst psf to be centered. */
					real pmax=dmax(sepsfi);
					dcog(pgrad, sepsfi, 0.5, 0.5, 0.1*pmax, 0.2*pmax, 0);
				}

				/*C_ABS causes sum of PSF to increase when there are negative values. Switch to literal copy.*/
				cembedd(seotfk, sepsfi, 0);
				cfftshift(seotfk);/*PSF, peak in corner; */
				cfft2(seotfk, -1);/*turn to OTF, peak in corner, max is 1 */
				cscale(seotfk, norm);/*normalized so that after fft, psf sum to 1*/
				if(pgrad[0]||pgrad[1]){
					ctilt(seotfk, -pgrad[0], -pgrad[1], 0);
				}
				if(petf){/*elongation. */
					ccwm(seotfk, PR(petf, isa, ii0));
				}else if(nominal){//nominal is fused into etf.
					ccwm(seotfk, nominal);
				}
				ccp(&seotfj, seotfk);/*save for later use*/
				if(pfotf){
					ccp(&P(P(*pfotf, ii0), isa, iwvl), seotfk);
				}
				cfft2(seotfk, 1);/*PSF with peak in center. sum to (pixtheta/dtheta)^2 due to nominal.*/
				/*no need fftshift because nominal is pre-treated */
				if(i0) dspmulcreal(P(P(i0, isa, ii0)), si, P(seotfk), wvlsig);
				if(gx||gy){
					ccp(&seotfk, seotfj);
					if(angles){//Apply derivative in rotated coordinate
						//derivative of i0 along radial/azimuthal direction
						const real ct=cos(angles[isa]);
						const real st=sin(angles[isa]);

						for(int iy=0; iy<notfy; iy++){
							for(int ix=0; ix<notfx; ix++){
								P(seotfk, ix, iy)*=ct*Ux[ix]+st*Uy[iy];
								P(seotfj, ix, iy)*=-st*Ux[ix]+ct*Uy[iy];
							}
						}
					} else{
						for(int iy=0; iy<notfy; iy++){
							for(int ix=0; ix<notfx; ix++){
								P(seotfk, ix, iy)*=Ux[ix];
								P(seotfj, ix, iy)*=Uy[iy];
							}
						}
					}
					if(gx){
						cfft2(seotfk, 1);
						dspmulcreal(P(P(gx, isa, ii0)), si, P(seotfk), wvlsig);
					}
					if(gy){
						cfft2(seotfj, 1);
						dspmulcreal(P(P(gy, isa, ii0)), si, P(seotfj), wvlsig);
					}
				}
				if(i0&&i0scale){
					real scale=area[isa]/dsum(P(i0, isa, ii0));
					if(i0) dscale(P(i0, isa, ii0), scale);
					if(gx) dscale(P(gx, isa, ii0), scale);
					if(gy) dscale(P(gy, isa, ii0), scale);
				}
			}/*for isa */
			cellfree(se_cache);
		}/*for ii0*/
	}/*iwvl */
}
