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

#include "common.h"
#include "genseotf.h"

/**
   see genseotf
*/
int count_unique(const dcell *opdbias){
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
	const void* amp, const dcell* opdbias,
	const void *saa, const dmat *wvl, real r0, real L0, int embfac){

	/*create a grid representing a sub-aperture. */
	loc_t* loc=mksqloc_auto(pts->nx, pts->nx, pts->dx, pts->dy);
	/*The embeding factor for embedding the aperture */
	const int npsfx=pts->nx*embfac;
	const int npsfy=pts->ny*embfac;
	const int nwvl=PN(wvl);
	const int nsa=pts->nsa;

	int notf=count_unique(opdbias);
	int notf2=iscell(amp)?count_unique((dcell*)amp):1;
	if(notf==1 && notf2!=1){
		notf=notf2;
	}else if(notf2!=1 && notf!=notf2){
		error("Mismatch for notf: %d vs %d\n", notf, notf2);
	}
	cccell *otf=cccellnew(notf, nwvl);
	
	info("There is %s bias\n", opdbias?"NCPA":"no");
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		for(int iotf=0; iotf<notf; iotf++){
			dmat* opdi=opdbias?P(opdbias,iotf):NULL;
			real thres=opdi?1:(1-1e-10);
			const dmat *ampi=iscell(amp)?PR((dcell*)amp,iotf,1):(dmat*)amp;
			const dmat* saai=iscell(saa)?PR((dcell*)saa, iotf, 1):(dmat*)saa;
			//OTFs are always generated with native sampling. It is upsampled at gensepsf if necessary.
			OMPTASK_SINGLE
				genotf(&P(otf,iotf,iwvl),loc, ampi, opdi, saai,
					thres, P(wvl,iwvl), NULL, r0, L0, npsfx, npsfy, nsa, 1);
		}
	}/*iwvl */
	locfree(loc);
	return otf;
}
/**
   Generates short exposure OTF by calling genotf() with p/t/t removal set.
*/
cccell* genseotf(const pts_t* pts, /**<[in]subaperture low left coordinate*/
				const void* amp, /**<[in] amplitude map. can be dcell or dmat*/
				const dcell* opdbias, /**<[in] opd bias for otf.*/
				const void* saa, /**<[in] list of subaperture area, can be dcell or dmat*/
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
		}else{
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
		CACHE, fnprefix, r0, L0, pts->dsa, pts->nsa,1./pts->dx, embfac);
	snprintf(fnlock, sizeof(fnlock), "%s.lock", fnotf);
	cccell *otf=0;
	while(!otf){
		if(exist(fnlock)||!zfexist("%s",fnotf)){/*need to create data */
			int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
			if(fd>=0){/*succeed */
				info("Generating WFS OTF for %s...", fnotf);
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
			info("Reading WFS OTF from %s\n", fnotf);
			otf=cccellread("%s", fnotf);
		}
	}
	return otf;
}
/**
   Upsample the otf in to out while preserving the PSF.
 */
static void upsample_otf(cmat* out, const cmat* in){
	if(in->nx==out->nx&&in->ny==out->ny){
		ccp(&out, in);
	} else{
		cmat* temp=0;
		ccp(&temp, in);
		cfft2(temp, -1);
		cscale(temp, 1./(in->nx*in->ny));
		czero(out);
		ccpcorner(out, temp, C_FULL);
		cfft2(out, 1);
		cfree(temp);
	}
}
/**
   Createing subaperture short exposure PSF from the tip/tilt removed turbulence
   OTF and uplink OTF. Not including detector or elongation characteristics.  */
void gensepsf(dccell** psepsfs, const cccell* otfs, const cccell* lotf, const void* saa, dmat* wvl, int notfx, int notfy
	){
	const int nwvl=PN(wvl);
	const int nsa=NX(P(otfs,0,0));
	const int nlotf=lotf?NX(lotf):0;
	int nsepsf=NX(otfs);
	if(nlotf>1){
		if (nsepsf==1){
			nsepsf=nlotf;
		}else if(nsepsf!=nlotf){
			error("mismatch: notf is %d, lotf is %d\n", nsepsf, nlotf);
		}
	}
	dccell *sepsfs=NULL; 
	cellinit((cell**)psepsfs, nsepsf, nwvl);

	if(!*psepsfs){
		*psepsfs=dccellnew(nsepsf, nwvl);
	}
	sepsfs=*psepsfs;
	for(int isepsf=0; isepsf<nsepsf; isepsf++){
		const dmat* saai=iscell(saa)?PR((dcell*)saa, isepsf, 1):(dmat*)saa;
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			P(sepsfs, isepsf, iwvl)=dcellnew_same(nsa, 1, notfx, notfy);
			const ccell* otf=PR(otfs, isepsf, iwvl);
			cmat* sepsf=cnew(notfx, notfy);
			cmat* lotfi=0;
			if(nlotf){
				cmat* lotfi2=P(PR(lotf, isepsf, iwvl),0);
				if(lotfi2->nx!=notfx||lotfi2->ny!=notfy){
					lotfi=cnew(notfx, notfy);
					upsample_otf(lotfi, lotfi2);
				} else{
					lotfi=cref(lotfi2);
				}
			}
			for(int isa=0; isa<nsa; isa++){
				real norm=P(saai, isa)/((real)(notfx*notfy));
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
void gensei(dcell **pi0, dcell **pgx, dcell **pgy, cccell **pfotf, cccell **ppotf,
	dccell *sepsfs, dtf_t *dtf, etf_t *etf, dcell *saa, dcell *srot, dmat *siglev, dmat *wvlwts, int dtrat,
	int i0scale, int radgx, int shift2center
){
	if(!sepsfs){
		error("sepsfs must be set\n");
	}
	const int notfx=NX(P(P(sepsfs, 0, 0), 0, 0)); 
	const int notfy=NY(P(P(sepsfs, 0, 0), 0, 0));
	const int nwvl=NY(sepsfs);
	const int nsa=NX(P(sepsfs,0,0));
	const int nllt=etf?NY(etf[0].etf):0;
	if(etf && nsa!=NX(etf[0].etf)){
		error("mismatch: nsa=%d, etf is %ldx%ld\n", nsa, NX(etf[0].etf), NY(etf[0].etf));
	}
	
	/**
	   ni0 may be greater than 1 in the following two cases
	   1) multiple LLT
	   2) different signal level or wvlwts
	   3) powfs[ipowfs].bkgrnd contains rayleigh scatter bkgrnd for each wfs in this powfs.
	*/
	const int nsepsf=NX(sepsfs);
	int ni0=MAX(nsepsf, nllt);
	if(NY(wvlwts)>1){
		if(ni0==1){
			ni0=NY(wvlwts);
		}else if(ni0!=NY(wvlwts)){
			error("Mismatch: ni0=%d, wvlwts is %ldx%ld\n", ni0, NX(wvlwts), NY(wvlwts));
		}
	}
	if(PN(siglev)>1){
		if(ni0==1){
			ni0=PN(siglev);
		} else if(ni0!=PN(siglev)){
			error("Mismatch: ni0=%d, siglev is %ldx1\n", ni0, NX(siglev));
		}
	}
	
	if(ni0>1){
		info("number of i0 for matched filter is %d\n", ni0);
	}
	const int pixpsax=dtf[0].pixpsax;
	const int pixpsay=dtf[0].pixpsay;

	if(pi0) {
		if(NX(*pi0)!=nsa || NY(*pi0)!=ni0){
			dcellfree(*pi0);
			*pi0=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		}
	}
	if(pgx){ 
		if(NX(*pgx)!=nsa || NY(*pgx)!=ni0){
			dcellfree(*pgx);
			*pgx=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		}
	}
	if(pgy){
		if(NX(*pgy)!=nsa||NY(*pgy)!=ni0){
			dcellfree(*pgy);
			*pgy=dcellnew_same(nsa, ni0, pixpsax, pixpsay);
		}
	}
	if(pfotf){
		cellfree(*pfotf);
		*pfotf=cccellnew(nsepsf, 1);
		for(int i=0; i<nsepsf; i++){
		 	P(*pfotf, i)=ccellnew(nsa, nwvl);
		}
	}
	if(ppotf){
		cellfree(*ppotf);
		*ppotf=cccellnew(nsepsf, 1);
		for(int i=0; i<nsepsf; i++){
			P(*ppotf,i)=ccellnew(nsa, nwvl);
		}
	}
	dcell* i0=*pi0;
	dcell* gx=pgx?*pgx:0;
	dcell* gy=pgy?*pgy:0;

	/*
	  Notice, the generation of shifted i0s are not accurate
	  because the PSF is not enough to cover the size.
	  Disable the computation.
	*/
	/*
	for(int ii0=0; ii0<ni0; ii0++){
		for(int isa=0; isa<nsa; isa++){
			P(i0, isa, ii0)=dnew(pixpsax, pixpsay);
			P(gx, isa, ii0)=dnew(pixpsax, pixpsay);
			P(gy, isa, ii0)=dnew(pixpsax, pixpsay);
		}
	}
	*/
	if(i0scale){
		warning("i0 is scaled to match sa area\n");
	}
	
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		const comp* Ux=dtf[iwvl].Ux->p;
		const comp* Uy=dtf[iwvl].Uy->p;
		const real norm=1./(real)(notfx*notfy);
		const ccell* petf=etf?etf[iwvl].etf:0;
		for(int ii0=0; ii0<ni0; ii0++){
			real* area=PR(saa,ii0,0)->p;
			real wvlsig=PR(wvlwts,iwvl, ii0)*PR(siglev, ii0, 0)*dtrat;
				
			dcell* psepsf=PR(sepsfs,ii0, iwvl);
			real* angles=(radgx)?(PR(srot,ii0,0)->p):0;
			ccell* se_save=ccellnew(3, NTHREAD);
#ifdef _OPENMP
			if(omp_in_parallel()){
				warning("Already in parallel\n");
			}
#endif
#pragma omp parallel default(shared)
#pragma omp for 
			for(int isa=0; isa<nsa; isa++){
				int ith=0;
				dmat *sepsfi=P(psepsf, isa, 0);
#ifdef _OPENMP
				ith=omp_get_thread_num();
#endif
#define seotfj P(se_save,0,ith)
#define seotfk P(se_save,1,ith)
				if(!seotfk){
					seotfk=cnew(notfx, notfy);
				}
				cmat* nominal=NULL;
				dsp* si=NULL;
				if(!dtf[iwvl].fused){
					nominal=PR(dtf[iwvl].nominal, isa, ii0);
				}
				si=PR(dtf[iwvl].si, isa, ii0); 
				real pgrad[2];
				/*loaded psepsf. sum to 1 for full sa. peak in center */
				if(shift2center){
					/*Forst psf to be centered. */
					real pmax=dmax(sepsfi);
					dcog(pgrad, sepsfi, 0.5, 0.5, 0.1*pmax, 0.2*pmax, 0);
				}
				if(dsum(sepsfi)>1.1){
					error("Short exposure PSF has wrong scaling. It should total to <=1\n");
				}
				/*C_ABS causes sum of PSF to increase when there are negative values. Switch to literal copy.*/
				cembedd(seotfk, sepsfi, 0);
				cfftshift(seotfk);/*PSF, peak in corner; */
				cfft2(seotfk, -1);/*turn to OTF, peak in corner, max is 1 */
				if(shift2center&&fabs(pgrad[0])>EPS&&fabs(pgrad[1])>EPS){
					ctilt(seotfk, -pgrad[0], -pgrad[1], 0);
				}
				if(nominal) ccwm(seotfk, nominal);
				cscale(seotfk, norm);/*normalized so that after fft, psf sum to 1*/
				if(ppotf){
					ccp(&P(P(*ppotf,ii0), isa, iwvl), seotfk);
				}
				if(nllt){/*elongation. */
					ccwm(seotfk, PR(petf, isa, ii0));
				}
				ccp(&seotfj, seotfk);/*backup */
				if(pfotf){
					ccp(&P(P(*pfotf,ii0), isa, iwvl), seotfk);
				}
				cfft2(seotfk, 1);/*PSF with peak in center. sum to (pixtheta/dtheta)^2 due to nominal.*/
				/*no need fftshift because nominal is pre-treated */
				dspmulcreal(P(i0, isa, ii0)->p, si, seotfk->p, wvlsig);
				if(gx || gy){
					ccp(&seotfk, seotfj);
					if(radgx){//Apply derivative in rotated coordinate
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
						dspmulcreal(P(gx, isa, ii0)->p, si, seotfk->p, wvlsig);
					}
					if(gy){
						cfft2(seotfj, 1);
						dspmulcreal(P(gy, isa, ii0)->p, si, seotfj->p, wvlsig);
					}
				}
				if(i0scale){
					real scale=area[isa]/dsum(P(i0, isa, ii0));
					if(i0) dscale(P(i0, isa, ii0), scale);
					if(gx) dscale(P(gx, isa, ii0), scale);
					if(gy) dscale(P(gy, isa, ii0), scale);
				}
			}/*for isa */
			cellfree(se_save);
		}/*for ii0*/
	}/*iwvl */
}

/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs){
	intstat_t* intstat=powfs[ipowfs].intstat;
	const real pixthetax=parms->powfs[ipowfs].radpixtheta;
	const real pixthetay=parms->powfs[ipowfs].pixtheta;
	const real rne=parms->powfs[ipowfs].rne;
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
	const int radgx=parms->powfs[ipowfs].radgx;
	int ni0=intstat->i0->ny;

	if(powfs[ipowfs].bkgrnd&&powfs[ipowfs].bkgrnd->ny>ni0){
		error("Please implement the case when bkgrnd has more columns\n");
	}
	info("Generating matched filter for %d\n", ipowfs);
	if(ni0!=1&&ni0!=parms->powfs[ipowfs].nwfs){
		error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
	}
	const int nsa=powfs[ipowfs].saloc->nloc;
	//Prevent printing of NEA during recomputing of matched filter
	//const int print_nea=intstat->mtche?0:1;
		
	dcellfree(powfs[ipowfs].sanea);
	dcell* sanea=powfs[ipowfs].sanea=dcellnew_same(ni0, 1, nsa, 3);
	dfree(intstat->i0sum);
	dmat* i0sum=intstat->i0sum=dnew(nsa, ni0);
	dfree(intstat->i0sumsum);
	intstat->i0sumsum=dnew(ni0, 1);

	const dcell* i0s=intstat->i0;
	const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx;
	const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy;
	
	long npix=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
	dcellfree(intstat->mtche);
	dcell* mtche=intstat->mtche=dcellnew_same(nsa, ni0, 2, npix);

	//dcell *saneaxy=powfs[ipowfs].saneaxy;
	int nllt;
	if(parms->powfs[ipowfs].llt){
		nllt=parms->powfs[ipowfs].llt->n;
	} else{
		nllt=0;
	}
	int irot_multiplier=nllt>1?1:0;
	const int mtchadp=parms->powfs[ipowfs].mtchadp;
	real sigratio=parms->powfs[ipowfs].sigrecon>0?(parms->powfs[ipowfs].sigrecon/parms->powfs[ipowfs].siglev):1;
	real sigratior=1./sigratio;
	if(sigratio<0){
		error("sigratio cannot be negative\n");
	}

	for(int ii0=0; ii0<ni0; ii0++){
		int iwfs=P(parms->powfs[ipowfs].wfs,ii0);
		const real siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
		real i0thres=MAX(0.1*siglev, rne*10);
		real nea2thres=pixthetax*pixthetay*100;
		real* srot=NULL;
		if(powfs[ipowfs].srot){
			int irot=ii0*irot_multiplier;
			srot=P(powfs[ipowfs].srot,irot)->p;
		}
		//P(sanea,ii0)=dnew(nsa, 3);
		dmat* psanea=P(sanea,ii0)/*PDMAT*/;
		real i0sumsum=0;
		int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
		int ncrdisable=0;
		dmat* nea2=0;
		for(int isa=0; isa<nsa; isa++){
			real pixrot=0;//pixel rotation
			if(srot&&parms->powfs[ipowfs].radpix){
				pixrot=srot[isa];
			}
			if(mtchadp){
				long fwhm=dfwhm_gauss(P(i0s, isa, ii0));
				if(fwhm>4){
					crdisable=0;
				} else{
					crdisable=1;
					ncrdisable++;
				}
			}
			dmat* bkgrnd2=NULL;
			dmat* bkgrnd2c=NULL;
			if(powfs[ipowfs].bkgrnd){
				bkgrnd2=P(powfs[ipowfs].bkgrnd, isa, ii0);
			}
			if(powfs[ipowfs].bkgrndc){
				bkgrnd2c=P(powfs[ipowfs].bkgrndc, isa, ii0);
			}
			mtch(&P(mtche, isa, ii0), &nea2, P(i0s, isa, ii0),
				gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
				parms->powfs[ipowfs].qe,
				bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
			if(fabs(sigratio-1)>1e-5){
				dscale(P(i0s, isa, ii0), sigratio);
				if(gxs) dscale(P(gxs, isa, ii0), sigratio);
				if(gys) dscale(P(gys, isa, ii0), sigratio);
				mtch(NULL, &nea2, P(i0s, isa, ii0),
					gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
					parms->powfs[ipowfs].qe,
					bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
					pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
				dscale(P(i0s, isa, ii0), sigratior);
				if(gxs) dscale(P(gxs, isa, ii0), sigratior);
				if(gys) dscale(P(gys, isa, ii0), sigratior);
			}
			P(i0sum, isa, ii0)=dsum(P(i0s, isa, ii0));
			i0sumsum+=P(i0sum, isa, ii0);

			if(P(i0sum, isa, ii0)<i0thres||P(nea2,0)>nea2thres||P(nea2,3)>nea2thres){
			//Signal level too low or error to high.
				P(nea2,0)=P(nea2,3)=nea2thres;
				P(nea2,1)=P(nea2,2)=0;
				dset(P(mtche, isa, ii0), 0);
			}
			if(parms->powfs[ipowfs].mtchcpl==0
				&&(!parms->powfs[ipowfs].radpix||radgx)){
			 /*remove coupling between r/a (x/y) measurements. */
				P(nea2,1)=P(nea2,2)=0;
			}
			P(psanea, isa, 0)=P(nea2,0);
			P(psanea, isa, 1)=P(nea2,3);
			P(psanea, isa, 2)=P(nea2,1);
		}/*isa  */
		dfree(nea2);

		if(mtchadp){
			info("Mtched filter contraint are disabled for %d subaps out of %d.\n",
				ncrdisable, nsa);
		}
		P(intstat->i0sumsum,ii0)=i0sumsum;
	}/*ii0 */
	if(1 /*print_nea*/){
		info2("Matched filter sanea:\n");
		if(powfs[ipowfs].sprint){/*print nea for select subapertures.*/
			for(int ii0=0; ii0<ni0; ii0++){
				int illt=0;
				if(ni0==parms->powfs[ipowfs].llt->n){
					illt=ii0;
				} else if(ni0==parms->powfs[ipowfs].nwfs&&parms->powfs[ipowfs].llt->n==1){
					illt=0;
				} else{
					error("Invalid combination\n");
				}
				info2("ii0 %d, llt %d.\n", ii0, illt);
				info2("sa index   dist   noise equivalent angle\n");
				dmat* psanea=P(sanea,ii0)/*PDMAT*/;
				for(int ksa=0; ksa<P(powfs[ipowfs].sprint,illt)->nx; ksa++){
					int isa=(int)P(P(powfs[ipowfs].sprint,illt),ksa);
					if(isa>0){
						info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n",
							isa, P(P(powfs[ipowfs].srsa,illt),isa),
							sqrt(P(psanea, isa, 0))*206265000,
							sqrt(P(psanea, isa, 1))*206265000);
					}
				}
			}
		} else{
			real dsa=powfs[ipowfs].saloc->dx;
			real llimit=-dsa/2;
			real ulimit=dsa/2;
			info2("index: position noise equivalent angle\n");
			for(int isa=0; isa<nsa; isa++){
				real locx=powfs[ipowfs].saloc->locx[isa];
				real locy=powfs[ipowfs].saloc->locy[isa];
				if((parms->powfs[ipowfs].llt&&(nsa<10||(locx>0&&locy>llimit&&locy<ulimit)))
					||(!parms->powfs[ipowfs].llt&&locx>=0&&locx<dsa*0.6&&locy>=0&&locy<dsa*0.6)
					||nsa<=4){
					info2("sa%4d:%6.1fm", isa, locx);
					for(int ii0=0; ii0<ni0; ii0++){
						info2(" (%4.1f,%4.1f)",
							sqrt(P(P(sanea,ii0), isa, 0))*206265000,
							sqrt(P(P(sanea,ii0), isa, 1))*206265000);
					}//for ii0
					info2(" mas\n");
				}
			}/*isa  */
		}
	}

	if(parms->recon.glao&&ni0>0){
		info("Averaging saneaxy of different WFS for GLAO mode\n");
		dmat* sanea2=0;
		real scale=1./ni0;
		for(int ii0=0; ii0<ni0; ii0++){
			dadd(&sanea2, 1, P(sanea,ii0), scale);
		}
		dcellfree(powfs[ipowfs].sanea);
		powfs[ipowfs].sanea=dcellnew(1, 1);
		P(powfs[ipowfs].sanea,0)=sanea2;
	}
}
