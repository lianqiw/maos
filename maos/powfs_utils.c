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



/*
  2009-11-26: changed to rotate OTF instead of psf to comply with the approach
  in wfsints.  this gives slightly larger intensity because max value of OTF is
  preserved which corresponds to the sum of psf.


*/

#include "common.h"
#include "powfs_utils.h"


void print_nea(const dcell* sanea, const lcell* sprint, const loc_t* saloc){
	//info2("Matched filter sanea:\n");
	if(sprint){/*print nea for select subapertures.*/
		for(int ii0=0; ii0<NX(sanea); ii0++){
			lmat* saindex=PR(sprint, ii0, 0);
			info2("sa index   location       noise equivalent angle\n");
			dmat* psanea=P(sanea, ii0);
			for(int ksa=0; ksa<NX(saindex); ksa++){
				long isa=P(saindex, ksa);
				if(isa>0){
					info2("%8ld: (%5.1f, %5.1f) m, (%6.2f, %6.2f) mas\n",
						isa, saloc->locx[isa], saloc->locy[isa],
						sqrt(P(psanea, isa, 0))*RAD2MAS,
						sqrt(P(psanea, isa, 1))*RAD2MAS);
				}
			}
		}
	} else{
		real dsa=saloc->dx;
		real llimit=-dsa*0.6;
		real ulimit=dsa*0.4;
		info2("sa index   location       noise equivalent angle\n");
		for(int isa=0; isa<saloc->nloc; isa++){
			real locx=saloc->locx[isa];
			real locy=saloc->locy[isa];
			if((locx>=0&&locy>llimit&&locy<ulimit)||saloc->nloc<=4){
				info2("%8d: (%5.1f, %5.1f) m", isa, locx, locy);
				for(int ii0=0; ii0<NX(sanea); ii0++){
					info2(" (%6.2f, %6.2f)",
						sqrt(P(P(sanea, ii0), isa, 0))*RAD2MAS,
						sqrt(P(P(sanea, ii0), isa, 1))*RAD2MAS);
				}//for ii0
				info2(" mas\n");
			}
		}/*isa  */
	}
}

void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs){
	info("Generating matched filter for powfs %d\n", ipowfs);
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
	if(powfs[ipowfs].srsa){
		print_nea(powfs[ipowfs].sanea, powfs[ipowfs].sprint, powfs[ipowfs].saloc);
	}
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
	dccell *i0m;
	dcell *grad;
	dcell *i0mv;
}fit_cache={0};
void fit_cache_free(){
	cellfree(fit_cache.sepsf);
	cellfree(fit_cache.i0m);
	cellfree(fit_cache.i0mv);
}

/**
 * Fit i0 to sodium profile using iterative algorithm. The steps are as follows
 * 1. Create sub images for each sodium profile bin
 * 2. Fit such sub-images against i0 to determine the profile
 * 3. Determine subaperture tip/tilt comparing i0 and fitted i0 (using matched filter)
 * 4. Repeat 1-4.
 * */
void sodium_fit(
	dmat** sodium, /**<The sodium profile determined by fit*/
	dcell** pgrad, /**<The gradients determined by fit.*/
	dcell** pi0,   /**<The output i0*/
	dcell** pgx,   /**<The output gx*/
	dcell** pgy,   /**<The output gy*/
	const dcell* i0i, /**<The input i0*/
	const dccell* sepsf,   /**<Short exposure PSF*/
	const dtf_t* dtf,     /**<Detector transfer function*/
	const loc_t* saloc,   /**<Saloc*/
	const dcell* saa,      /**<Subaperture area. */
	const dcell* srsa,    /**<Subaperture to LLT distance*/
	const dcell* srot,    /**<Subaperture to LLT clocking*/
	const dmat* siglev,  /**<Subaperture signal level*/
	const dmat* wvlwts,    /**<Wavelength weights*/
	const dcell* gradncpa,/**<NCPA gradient to be used for pi0,pgx,pgy output.*/
	real dh,      /**<The sodium profile sampling in meters*/
	real hs,      /**<LGS focusing height*/
	real htel,    /**<Telescope hegith*/
	real za,      /**<Telescope zenith angle*/
	real svdthres, /**<SVD threshold*/
	int nrep,     /**<Number of iterations*/
	int save,      /**<Save results to file*/
	int use_cache  /**<Use cache*/
){
	static int count=-1; count++;
	//const real ht=25000;//total thickness
	const real hmin=80000;
	const real hmax=105000;//wrapp around happens at 10000 for 10 pixel if alined along x/y.

	long nh=(long)floor((hmax-hmin)/dh)+1;
	const int radgx=0;
	dmat *nai=(sodium&&*sodium)?*sodium:NULL;
	if(!nai){
		nai=dnew(nh, 2);
		if(sodium) *sodium=nai;
	}else if(NX(nai)!=nh || NY(nai)!=2){
		dresize(nai, nh, 2);
	}

	for(long ix=0; ix<nh; ix++){
		P(nai, ix, 0)=hmin+ix*dh;
	}
	const long nsa=NX(i0i);
	const long ni0=NY(i0i);
	real pixthetax=dtf[0].pixthetax;
	real pixthetay=dtf[0].pixthetay;

	dcell* i0tmp=0;
	dcell* gxtmp=0;
	dcell* gytmp=0;
	dcell* gradtmp=0;
	
	//avoid overriding i0 input. Array for gensei full output to build temporal matched filter
	dcell **pi0tmp=(pi0 && (nrep==1 || *pi0 != i0i))?pi0:&i0tmp;
	dcell **pgxtmp=pgx?pgx:&gxtmp;//ok to override input
	dcell **pgytmp=pgy?pgy:&gytmp;
	if(!pgrad) pgrad=&gradtmp;
	if(!*pgrad && (nrep>1 || pgrad!=&gradtmp)){//need output
		*pgrad=dcellnew_same(ni0, 1, nsa*2, 1);
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
	//etf_t** etfs=use_cache?fit_cache.etfs:NULL;
	dccell* i0m=use_cache?fit_cache.i0m:NULL;//sa image for each sodium bin
	dcell* i0mv=use_cache?fit_cache.i0mv:NULL;//vectorized i0m
	dcell *grad=use_cache?fit_cache.grad:*pgrad;
	int skip_first=0;

	if(!i0m){//prepare to create the sub layer image model
		i0mv=dcellnew(1, nh);
		i0m=dccellnew(nh, 1);
		if(gradncpa){
			dcellcp(&grad, gradncpa);
		}else{
			grad=dcellnew_same(ni0, 1, nsa*2, 1);
		}
		dbg("Initial gradient uses %s.\n", gradncpa?"gradncpa":"zero");
		
		if(use_cache){
			fit_cache.i0m=i0m;
			fit_cache.i0mv=i0mv;
			fit_cache.grad=grad;
		}
		if(save){
			writebin(grad, "sodium_grad_in");
		}
	}else{
		skip_first=1;//cache available.
		dbg("reuse previous grad, i0m, i0mv\n");
	}

	dbg("svdthres=%g, nrep=%d\n", svdthres, nrep);
	//don't try to cachc fotf. It is per WFS and uses too much storage.
	dcell* ata=0, * atb=0;
	etf_t *etf_full=0;
	dcell *na2s=dcellnew_same(nh, 1, 1, 2);
	for(int irep=0; irep<nrep; irep++){
		dbg("repeat %d of %d\n", irep+1, nrep);
		if(irep>0 || !skip_first){//Compute subaperture sublayer imaging model
			if(irep>0 && &grad!=pgrad){
				dcellcp(&grad, *pgrad);
			}
			OMP_TASK_FOR(4)
			for(long ix=0; ix<nh; ix++){
				dmat *na2i=P(na2s, ix);
				P(na2i, 0, 0)=P(nai, ix, 0);
				P(na2i, 0, 1)=1;
				//ETF takes a lot of storage but is inexpensive to build. So we choose to build it on the fly
				etf_t *etf_i=mketf(dtf, na2i, 0, srot, srsa, hs, htel, za, 1);
				gensei(&P(i0m, ix), NULL, NULL, NULL, sepsf, dtf, etf_i, saa, radgx?srot:NULL, siglev, wvlwts, grad, 0, 0);
				etf_free(etf_i);
				if(!P(i0mv, 0, ix)||P(P(i0mv, 0, ix))!=P(P(i0m, ix)->m)){
					dfree(P(i0mv, 0, ix));
					P(i0mv, 0, ix)=dref(P(i0m, ix)->m);
				}

			}
			toc2("gensei each");tic;
			print_mem("gensei");
		}
		dcellzero(ata);
		dcellzero(atb);
		dcellmm(&ata, i0mv, i0mv, "tn", 1);
		dcellmm(&atb, i0mv, i02, "tn", 1);
		if(save){
			writebin(ata, "sodium_ata_%d_%d", count, irep);
			writebin(atb, "sodium_atb_%d_%d", count, irep);
			if(count==0&&irep==0) writebin(i0m, "sodium_i0m_%d_%d", count, irep);
		}
		dcellsvd_pow(ata, -1, svdthres);
		dcellzero(res);
		dcellmm(&res, ata, atb, "nn", 1);
		real scale=1./dcellsum(res);//make sure nai sum to 1.
		for(long ix=0; ix<nh; ix++){
			P(nai, ix, 1)=P(P(res, ix), 0)*scale;
		}
		if(nrep>1 || pgrad!=&gradtmp){//need to determine error in applying gradient offset
			if(etf_full) etf_free(etf_full);
			//mketf for full profile must use the same no_interp flag 
			etf_full=mketf(dtf, nai, 0, srot, srsa, hs, htel, za, 1);
			toc2("mketf full"); tic;
			gensei(pi0tmp, pgxtmp, pgytmp, NULL, sepsf, dtf, etf_full, saa, radgx?srot:NULL, siglev, wvlwts, grad, 0, 0);
			toc2("gensei full"); tic;
			mtch_cell(&mtche, NULL, NULL, NULL, *pi0tmp, *pgxtmp, *pgytmp, NULL, NULL, NULL, 0, 0, 3,
				pixthetax, pixthetay, NULL, radgx, 1, 1);
			toc2("mtche create"); tic;

			for(long ii0=0; ii0<ni0; ii0++){
				dmat* grad1=P(grad, ii0);//model
				dmat* grad2=P(*pgrad, ii0);//output
OMP_TASK_FOR(8)
				for(long isa=0; isa<nsa; isa++){
					real g[2]={0,0};
					dmulvec(g, P(mtche, isa, ii0), P(P(i0i, isa, ii0)), 1);
					P(grad2, isa)    =P(grad1, isa    )+g[0];
					P(grad2, isa+nsa)=P(grad1, isa+nsa)+g[1];
				}
			}
			toc2("mtche apply"); tic;
			
OMP_TASK_FOR(4)
			for(long ii0=0; ii0<ni0; ii0++){
				//Remove focus mode from the gradients as it degenerates with sodium profile shift.
				loc_remove_focus_grad(P(*pgrad, ii0), saloc, 1);
			}
			if(save){
				writebin(nai, "sodium_prof_%d_%d", count, irep);
				writebin(*pgrad, "sodium_grad_%d_%d", count, irep);
				//writebin(mtche, "sodium_mtche_%d_%d", count, irep);
				//writebin(*pi0tmp, "sodium_i0_%d_%d", count, irep);
			}
		}
	}

	//output is desired. build final i0, gx, gy with the final gradient or ncpa gradient.
	if(pi0 || pgx || pgy){
		dbg("Replacing i0, gx, gy with fitted value\n");
		if(!etf_full){
			etf_full=mketf(dtf, nai, 0, srot, srsa, hs, htel, za, 1);
			toc2("mketf final");tic;
		}
		const dcell *gradf=gradncpa?gradncpa:(pgrad?(*pgrad):grad);
		gensei(pi0, pgx, pgy, NULL, sepsf, dtf, etf_full, saa, radgx?srot:NULL, siglev, wvlwts, gradf, 0, 0);
		toc2("gensei final");tic;
	}

	if(etf_full) etf_free(etf_full);
	dcellfree(ata);
	dcellfree(atb);
	if(!use_cache){
		cellfree(i0m);
		cellfree(i0mv);
		cellfree(grad);
	}
	cellfree(mtche);
	cellfree(i02);
	cellfree(res);
	if(!sodium) dfree(nai);
	cellfree(na2s);
	cellfree(i0tmp);
	cellfree(gxtmp);
	cellfree(gytmp);
	cellfree(gradtmp);
	
}
/**
 * Fit i0 to sodium profile and replace i0, gx, gy with derived parameters
 * */
void sodium_fit_wrap(dmat** psodium, /**<[out] sodium profile*/
	dcell** pgrad, /**<[out] estimated actual gradient*/
	dcell** pi0,   /**<[out] The output i0*/
	dcell** pgx,   /**<[out] The output gx*/
	dcell** pgy,   /**<[out] The output gy*/
	const dcell* i0in, /**<[in]The input sa intensities. may equal to *pi0 */
	const parms_t* parms,/**<[in]parms*/
	powfs_t* powfs, /**<[in]powfs*/
	int ipowfs, /**<[in] ipowfs*/
	real r0,  /**<[in] Fried parameter*/
	real L0,  /**<[in] outer scale*/
	int nrep, /**<[in] Number of iterations. 1 for mtche, 3 for cog*/
	int use_cache /**<[in] cache intermediate results.*/
	){
	dccell* sepsf=use_cache?fit_cache.sepsf:NULL;
	if(!sepsf){
		cccell* otf=NULL, * lotf=NULL;
		otf=genseotf(powfs[ipowfs].pts, powfs[ipowfs].realamp,
			NULL, powfs[ipowfs].realsaa, parms->powfs[ipowfs].wvl, r0, L0,
			parms->powfs[ipowfs].embfac);
		if(parms->powfs[ipowfs].llt){
						//genselotf(parms, powfs, ipowfs);
			lotf=genseotf(powfs[ipowfs].llt->pts, powfs[ipowfs].llt->amp,
				NULL, NULL, parms->powfs[ipowfs].wvl, r0, L0, 
				parms->powfs[ipowfs].embfac);
		}
		gensepsf(&sepsf, otf, lotf, powfs[ipowfs].realsaa,
			parms->powfs[ipowfs].wvl, powfs[ipowfs].notfx, powfs[ipowfs].notfy);
		cellfree(otf);
		cellfree(lotf);
		if(use_cache){
			fit_cache.sepsf=sepsf;
		}
		if(parms->save.dither){
			static int count=0;
			writebin(sepsf, "sodium_sepsf_%d", count);
			count++;
		}
	}

	sodium_fit(psodium, pgrad, pi0, pgx, pgy, i0in,
		sepsf, powfs[ipowfs].dtf, powfs[ipowfs].saloc, powfs[ipowfs].realsaa,
		powfs[ipowfs].srsa, powfs[ipowfs].srot,
		parms->powfs[ipowfs].siglevs, parms->powfs[ipowfs].wvlwts, powfs[ipowfs].gradncpa,
		parms->powfs[ipowfs].llt->na_fit_dh, parms->powfs[ipowfs].hs, parms->sim.htel, parms->sim.za, 
		parms->powfs[ipowfs].llt->na_fit_svdthres, nrep, parms->save.dither>1, use_cache);//parms->save.setup);
	
	/*if(parms->save.setup){
		if(pi0) writebin(*pi0, "powfs%d_i0_fit", ipowfs);
		if(pgx) writebin(*pgx, "powfs%d_gx_fit", ipowfs);
		if(pgy) writebin(*pgy, "powfs%d_gy_fit", ipowfs);
		if(psodium) writebin(*psodium, "powfs%d_sodium_fit", ipowfs);
	}*/
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