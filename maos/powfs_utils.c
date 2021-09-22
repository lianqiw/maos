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
#include "powfs_utils.h"


void print_nea(const dcell* sanea, const dcell* sprint, const loc_t* saloc, const dcell* srsa){
	info2("Matched filter sanea:\n");
	if(sprint){/*print nea for select subapertures.*/
		for(int ii0=0; ii0<NX(sanea); ii0++){
			dmat* saindex=PR(sprint, ii0, 0);
			info2("ii0 %d\n", ii0);
			info2("sa index   dist   noise equivalent angle\n");
			dmat* psanea=P(sanea, ii0);
			for(int ksa=0; ksa<NX(saindex); ksa++){
				int isa=(int)P(saindex, ksa);
				if(isa>0){
					info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n",
						isa, P(PR(srsa, ii0, 0), isa),
						sqrt(P(psanea, isa, 0))*206265000,
						sqrt(P(psanea, isa, 1))*206265000);
				}
			}
		}
	} else{
		real dsa=saloc->dx;
		real llimit=-dsa*0.6;
		real ulimit=dsa*0.4;
		info2("index: position noise equivalent angle\n");
		for(int isa=0; isa<saloc->nloc; isa++){
			real locx=saloc->locx[isa];
			real locy=saloc->locy[isa];
			if((locx>=0&&locy>llimit&&locy<ulimit)||saloc->nloc<=4){
				info2("sa%4d:%6.1fm", isa, locx);
				for(int ii0=0; ii0<NX(sanea); ii0++){
					info2(" (%4.1f,%4.1f)",
						sqrt(P(P(sanea, ii0), isa, 0))*206265000,
						sqrt(P(P(sanea, ii0), isa, 1))*206265000);
				}//for ii0
				info2(" mas\n");
			}
		}/*isa  */
	}
}

void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs){
	info("Generating matched filter for %d\n", ipowfs);
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;

	intstat_t* intstat=powfs[ipowfs].intstat;
	const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx;
	const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy;
	real sigratio=parms->powfs[ipowfs].sigrecon>0?(parms->powfs[ipowfs].sigrecon/parms->powfs[ipowfs].siglev):1;
	const int mtchadp=parms->powfs[ipowfs].mtchadp;
	int mtchcr=mtchadp?-1:parms->powfs[ipowfs].mtchcr;
	mtch_cell(&intstat->mtche, &powfs[ipowfs].sanea, &intstat->i0sum, &intstat->i0sumsum,
		intstat->i0, gxs, gys, parms->powfs[ipowfs].qe,
		powfs[ipowfs].bkgrnd, powfs[ipowfs].bkgrndc, bkgrnd, bkgrndc,
		parms->powfs[ipowfs].rne, parms->powfs[ipowfs].radpixtheta, parms->powfs[ipowfs].pixtheta,
		parms->powfs[ipowfs].radpix?powfs[ipowfs].srot:NULL, parms->powfs[ipowfs].radgx, mtchcr, sigratio
	);
	print_nea(powfs[ipowfs].sanea, powfs[ipowfs].sprint, powfs[ipowfs].saloc, powfs[ipowfs].srsa);
}


/*compute cog NEA using Monte Carlo realizations of noise*/
void cog_nea(real* nea, const dmat* ints, real cogthres, real cogoff, int ntry,
	rand_t* rstat, real bkgrnd, real bkgrndc, const dmat* bkgrnd2i, const dmat* bkgrnd2ic, real rne
){
	dmat* ints2=dnew(NX(ints), NY(ints));
	real gnf[2]={0,0};
	real gny[2]={0,0};
	dcog(gnf, ints, 0, 0, cogthres, cogoff, 0);
	seed_rand(rstat, 1);/*reset the seed each time.*/
	nea[0]=0; nea[1]=0; nea[2]=0; nea[3]=0;
	for(int i=0; i<ntry; i++){
		dcp(&ints2, ints);
		addnoise(ints2, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, 0, rne, 1);
		dcog(gny, ints2, 0, 0, cogthres, cogoff, 0);
		real errx=gny[0]-gnf[0];
		real erry=gny[1]-gnf[1];
		nea[0]+=errx*errx;
		nea[1]+=errx*erry;
		nea[3]+=erry*erry;
	}
	dfree(ints2);
	real stry=1./ntry;
	nea[0]=nea[0]*stry;
	nea[3]=nea[3]*stry;
	nea[1]=nea[1]*stry;
	nea[2]=nea[1];
}
struct fit_cache{
	dccell *sepsf;
	etf_t **etfs;
	dccell *i0m;
	dcell *i0mv;
	int nx;
}fit_cache={0};
void fit_cache_free(){
	if(!fit_cache.nx) return;
	cellfree(fit_cache.sepsf);
	cellfree(fit_cache.i0m);
	cellfree(fit_cache.i0mv);
	for(long ix=0; ix<fit_cache.nx; ix++){
		etf_free(fit_cache.etfs[ix]);
	}
	free(fit_cache.etfs);
	fit_cache.etfs=0;
	fit_cache.nx=0;
}
/**
 * Fit i0 to sodium profile using iterative algorithm. The steps are as follows
 * 1. Create sub images for each sodium profile bin
 * 2. Fit such sub-images against i0 to determine the profile
 * 3. Determine subaperture tip/tilt comparing i0 and fitted i0 (using matched filter)
 * 4. Repeat 1-4.
 * */
void fit_sodium_profile(
	dmat** sodium, /**<The sodium profile determined by fit*/
	dcell** pgrad, /**<The gradients determined by fit.*/
	dcell** pi0,   /**<The output i0*/
	dcell** pgx,   /**<The output gx*/
	dcell** pgy,   /**<The output gy*/
	const dcell* i0i, /**<The input i0*/
	const dccell* sepsf,   /**<Short exposure PSF*/
	const dtf_t* dtf,     /**<Detector transfer function*/
	const void* saa,      /**<Subaperture area. dmat or dcell*/
	const dcell* srsa,    /**<Subaperture to LLT distance*/
	const dcell* srot,    /**<Subaperture to LLT clocking*/
	const dmat* siglev,  /**<Subaperture signal level*/
	const dmat* wvlwts,    /**<Wavelength weights*/
	const dcell* gradncpa,/**<NCPA gradient to be used for pi0,pgx,pgy output.*/
	real dh,      /**<The sodium profile sampling in meters*/
	real hs,      /**<LGS focusing height*/
	real htel,    /**<Telescope hegith*/
	real za,      /**<Telescope zenith angle*/
	real cogthres,/**<Threshold for cog*/
	real tikcr,   /**<Tikhonov regularization*/
	real svdthres, /**<SVD threshold*/
	int use_mtche, /**<Use mtche to compute gradient error*/
	int nrep,     /**<Number of iterations*/
	int save,      /**<Save results to file*/
	int use_cache  /**<Use cache*/
){
	static int count=-1; count++;
	//const double ht=25000;//total thickness
	const double hmin=80000;
	const double hmax=105000;

	long nx=(long)floor((hmax-hmin)/dh)+1;
	const int radgx=0;
	dmat *nai=(sodium&&*sodium)?*sodium:NULL;
	if(!nai){
		nai=dnew(nx, 2);
		if(sodium) *sodium=nai;
	}else if(NX(nai)!=nx || NY(nai)!=2){
		dresize(nai, nx, 2);
	}

	for(long ix=0; ix<nx; ix++){
		P(nai, ix, 0)=hmin+ix*dh;
	}
	const long nsa=NX(i0i);
	const long ni0=NY(i0i);
	real pixthetax=dtf[0].pixthetax;
	real pixthetay=dtf[0].pixthetay;
	if(*pi0&&NY(*pi0)!=ni0){
		warning("i0 has wrong dimensions %ldx%ld, recreate i0, gx, gy.\n", NX(*pi0), NY(*pi0));
		dcellfree(*pi0);
		dcellfree(*pgx);
		dcellfree(*pgy);
	}
	dcell* i0tmp=0;
	dcell* gxtmp=0;
	dcell* gytmp=0;
	dcell* gradtmp=0;
	
	//avoid overriding i0 input
	dcell **pi0tmp=(pi0 && (nrep==1 || *pi0 != i0i))?pi0:&i0tmp;
	dcell **pgxtmp=pgx?pgx:&gxtmp;//ok to override input
	dcell **pgytmp=pgy?pgy:&gytmp;
	if(!pgrad) pgrad=&gradtmp;
	else if(*pgrad && NY(*pgrad)!=ni0){
		warning("pgrad has wrong dimensions %ldx%ld, recreate\n", NX(*pgrad), NY(*pgrad));
		dcellfree(*pgrad);
	}
	if(!*pgrad){
		*pgrad=dcellnew_same(ni0, 1, nsa*2, 1);
		dbg("Initial gradient uses %s.\n", gradncpa?"gradncpa":"cog of i0");
		for(long ii0=0; ii0<ni0; ii0++){
			dmat *gradi=P(*pgrad, ii0);
			if(gradncpa){
				dcp(&gradi, PR(gradncpa, ii0, 0));
			}else{
				for(long isa=0; isa<nsa; isa++){
					real g[2]={0,0};
					dcog(g, P(i0i, isa, ii0), 0, 0, cogthres, 0, 0);
					P(gradi, isa)    =g[0]*pixthetax;
					P(gradi, isa+nsa)=g[1]*pixthetay;
				}
			}
		}
		if(save){
			writebin(*pgrad, "sodium_grad_in");
		}
	}
	print_mem("grad");
	TIC;tic;
	dcell* i02=dcellnew(1, 1);
	if(!i0i->m){
		error("i0i->m is not set\n");
	}
	P(i02, 0)=dref(i0i->m);
	dcell* mtche=0;
	dcell* res=0;
	etf_t** etfs=use_cache?fit_cache.etfs:NULL;
	dccell* i0m=use_cache?fit_cache.i0m:NULL;//sa image for each sodium bin
	dcell* i0mv=use_cache?fit_cache.i0mv:NULL;//vectorized i0m
	int skip_first=0;
	if(!etfs){
		etfs=mycalloc(nx, etf_t*);
		dmat *na2i=dnew(1,2);
		//OMP_FOR	
		for(long ix=0; ix<nx; ix++){
			P(na2i, 0, 0)=P(nai, ix, 0);
			P(na2i, 0, 1)=1;
			etfs[ix]=mketf(dtf, na2i->base, 0, srot, srsa, hs, htel, za, 1);//no need to update
		}
		dfree(na2i);
		toc2("mketf");tic;
		print_mem("mketf");
		i0mv=dcellnew(1, nx);
		i0m=dccellnew(nx, 1);
		if(use_cache){
			fit_cache.etfs=etfs;
			fit_cache.i0m=i0m;
			fit_cache.i0mv=i0mv;
			fit_cache.nx=nx;
		}
	}else{
		skip_first=1;
		dbg("reuse previous etf, i0m, i0mv\n");
	}
	dbg("svdthres=%g, tikcr=%g, nrep=%d\n", svdthres, tikcr, nrep);
	//don't try to cachc fotf. It is per WFS and uses too much storage.
	dcell* ata=0, * atb=0;
	etf_t *etfi=0;
	for(int irep=0; irep<nrep; irep++){
		dbg("repeat %d of %d\n", irep+1, nrep);
		if(irep>0 || !skip_first){
			OMP_FOR
			for(long ix=0; ix<nx; ix++){
				gensei(&P(i0m, ix), NULL, NULL, NULL, sepsf, dtf, etfs[ix], saa, radgx?srot:NULL, siglev, wvlwts, *pgrad, 0, 0);
				if(!P(i0mv, 0, ix)||P(P(i0mv, 0, ix))!=P(P(i0m, ix)->m)){
					dfree(P(i0mv, 0, ix));
					P(i0mv, 0, ix)=dref(P(i0m, ix)->m);
				}
			}
			toc2("gensei");tic;
			print_mem("gensei");
		}
		dcellzero(ata);
		dcellzero(atb);
		dcellmm(&ata, i0mv, i0mv, "tn", 1);
		dcellmm(&atb, i0mv, i02, "tn", 1);

		dcellsvd_pow(ata, -1, svdthres, tikcr);
		dcellzero(res);
		dcellmm(&res, ata, atb, "nn", 1);
		real scale=1./dcellsum(res);//make sure nai sum to 1.
		for(long ix=0; ix<nx; ix++){
			P(nai, ix, 1)=P(P(res, ix), 0)*scale;
		}
		if(use_mtche){
			if(etfi) etf_free(etfi);
			//mketf for full profile must use the same no_interp flag 
			etfi=mketf(dtf, nai->base, 0, srot, srsa, hs, htel, za, 1);
			toc2("mketf full"); tic;
			gensei(pi0tmp, pgxtmp, pgytmp, NULL, sepsf, dtf, etfi, saa, radgx?srot:NULL, siglev, wvlwts, *pgrad, 0, 0);
			toc2("gensei full"); tic;
			mtch_cell(&mtche, NULL, NULL, NULL, *pi0tmp, *pgxtmp, *pgytmp, NULL, NULL, NULL, 0, 0, 3,
				pixthetax, pixthetay, NULL, radgx, 1, 1);
			toc2("mtche create"); tic;
			OMP_FOR_COLLAPSE(2)
			for(long ii0=0; ii0<ni0; ii0++){
				for(long isa=0; isa<nsa; isa++){
					dmat* gradi=P(*pgrad, ii0);
					real g[2]={0,0};
					dmulvec(g, P(mtche, isa, ii0), P(P(i0i, isa, ii0)), 1);
					P(gradi, isa)+=g[0];
					P(gradi, isa+nsa)+=g[1];
				}
			}
			toc2("mtche apply"); tic;
		}else{
			if(*pi0tmp) dcellzero(*pi0tmp);
			for(long ix=0; ix<nx; ix++){
				dcelladd(pi0tmp, 1, P(i0m, ix), P(nai, ix, 1));
			}
			const dcell* i0o=*pi0tmp;
			OMP_FOR_COLLAPSE(2)
			for(long ii0=0; ii0<ni0; ii0++){
				for(long isa=0; isa<nsa; isa++){
					dmat* gradi=P(*pgrad, ii0);
					real g1[2],g2[2];
					dcog(g1, P(i0o, isa, ii0), 0, 0, cogthres, 0, 0);
					dcog(g2, P(i0i, isa, ii0), 0, 0, cogthres, 0, 0);
					//i0o has real gradient of pgrad.
					P(gradi, isa)+=(g2[0]-g1[0])*pixthetax;
					P(gradi, isa+nsa)+=(g2[1]-g1[1])*pixthetay;
				}
			}
			toc2("tcog diff"); tic;
		}
		if(save){
			writebin(ata, "sodium_ata_%d_%d", count, irep);
			writebin(atb, "sodium_atb_%d_%d", count, irep);
			writebin(nai, "sodium_prof_%d_%d", count, irep);
			writebin(*pgrad, "sodium_grad_%d_%d", count, irep);
			writebin(mtche, "sodium_mtche_%d_%d", count, irep);
			writebin(*pi0, "sodium_i0_%d_%d", count, irep);
		}
	}
	//output is desired. build final i0, gx, gy with the final gradient or ncpa gradient.
	if(pi0 || pgx || pgy){
		if(!etfi) etfi=mketf(dtf, nai->base, 0, srot, srsa, hs, htel, za, 1);
		const dcell *gradf=gradncpa?gradncpa:(*pgrad);
		gensei(pi0, pgx, pgy, NULL, sepsf, dtf, etfi, saa, radgx?srot:NULL, siglev, wvlwts, gradf, 0, 0);
	}
	toc2("gensei final");tic;

	if(etfi) etf_free(etfi);
	dcellfree(ata);
	dcellfree(atb);
	if(!use_cache){
		cellfree(i0m);
		cellfree(i0mv);
	
		for(long ix=0; ix<nx; ix++){
			etf_free(etfs[ix]);
		}
		free(etfs);
	}
	cellfree(mtche);
	cellfree(i02);
	cellfree(res);
	if(!sodium) dfree(nai);

	cellfree(i0tmp);
	cellfree(gxtmp);
	cellfree(gytmp);
	cellfree(gradtmp);
}
/**
 * Fit i0 to sodium profile and replace i0, gx, gy with derived parameters
 * */
void fit_sodium_profile_wrap(dmat** psodium, dcell** pgrad, const dcell* i0in, const parms_t* parms, 
	powfs_t* powfs, int ipowfs, int use_mtche, int nrep, int use_ncpa, int use_cache){
	dccell* sepsf=use_cache?fit_cache.sepsf:NULL;
	if(!sepsf){
		cccell* otf=NULL, * lotf=NULL;
		otf=genseotf(powfs[ipowfs].pts, powfs[ipowfs].realamp,
			NULL, powfs[ipowfs].realsaa, parms->powfs[ipowfs].wvl,
			parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
			parms->powfs[ipowfs].embfac);
		if(parms->powfs[ipowfs].llt){
						//genselotf(parms, powfs, ipowfs);
			lotf=genseotf(powfs[ipowfs].llt->pts, powfs[ipowfs].llt->amp,
				NULL, NULL, parms->powfs[ipowfs].wvl,
				parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
				parms->powfs[ipowfs].embfac);
		}
		gensepsf(&sepsf, otf, lotf, powfs[ipowfs].realsaa,
			parms->powfs[ipowfs].wvl, powfs[ipowfs].notfx, powfs[ipowfs].notfy);
		cellfree(otf);
		cellfree(lotf);
		if(use_cache){
			fit_cache.sepsf=sepsf;
		}
	}
	//i0, gx, gy
	if(!powfs[ipowfs].intstat){
		powfs[ipowfs].intstat=mycalloc(1, intstat_t);
	}
	intstat_t* intstat=powfs[ipowfs].intstat;
	
	info("Replacing i0, gx, gy with fitted value\n");
	
	fit_sodium_profile(psodium, pgrad, &intstat->i0, &intstat->gx, &intstat->gy, i0in,
		sepsf, powfs[ipowfs].dtf, powfs[ipowfs].realsaa,
		powfs[ipowfs].srsa, powfs[ipowfs].srot,
		parms->powfs[ipowfs].siglevs, parms->powfs[ipowfs].wvlwts, use_ncpa?powfs[ipowfs].gradncpa:NULL,
		parms->dbg.na_fit_dh, parms->powfs[ipowfs].hs, parms->sim.htel, parms->sim.za, parms->powfs[ipowfs].cogthres,
		0, parms->dbg.na_fit_svdthres, use_mtche, nrep, 0, use_cache);//parms->save.setup);
	
	if(parms->save.setup){
		writebin(intstat->i0, "powfs%d_i0_fit", ipowfs);
		writebin(intstat->gx, "powfs%d_gx_fit", ipowfs);
		writebin(intstat->gy, "powfs%d_gy_fit", ipowfs);
		if(psodium) writebin(*psodium, "powfs%d_sodium_fit", ipowfs);
	}
	if(!use_cache){
		dcellfree(sepsf);
	}else{
		static int registered=0;
		if(!registered){
			registered=1;
			register_deinit(fit_cache_free, NULL);
		}
	}
}