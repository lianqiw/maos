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
	dmat* ints2=dnew(ints->nx, ints->ny);
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
	real dh,      /**<The sodium profile sampling in meters*/
	real hs,      /**<LGS focusing height*/
	real htel,    /**<Telescope hegith*/
	real za,      /**<Telescope zenith angle*/
	real tikcr,   /**<Tikhonov regularization*/
	real svdthres, /**<SVD threshold*/
	int save      /**<Save results to file*/
){
	dcell* nac=dcellnew(1, 1);
	//const double ht=25000;//total thickness
	const double hmin=80000;
	const double hmax=105000;

	long nx=(long)floor((hmax-hmin)/dh)+1;
	const int radgx=0;
	dmat* nai=P(nac, 0)=dnew(nx, 2);
	if(sodium){
		*sodium=dref(nai);
	}
	for(long ix=0; ix<nx; ix++){
		P(nai, ix, 0)=hmin+ix*dh;
	}
	const long nsa=NX(i0i);
	const long ni0=NY(i0i);
	real pixthetax=dtf[0].pixthetax;
	real pixthetay=dtf[0].pixthetay;
	dcell* i0=0;
	dcell* gx=0;
	dcell* gy=0;
	dcell* grad=0;
	if(!pi0) pi0=&i0;
	if(!pgx) pgx=&gx;
	if(!pgy) pgy=&gy;
	if(!pgrad) pgrad=&grad;
	if(!*pgrad){
		*pgrad=dcellnew_same(ni0, 1, nsa, 2);
	}
	TIC;tic;
	for(long ii0=0; ii0<ni0; ii0++){
		for(long isa=0; isa<nsa; isa++){
			real g[2]={0,0};
			dcog(g, P(i0i, isa, ii0), 0, 0, 9, 0, 0);
			P(*pgrad, isa, 0, ii0, 0)=g[0]*pixthetax;
			P(*pgrad, isa, 1, ii0, 0)=g[1]*pixthetay;
		}
	}
	toc("cog");tic;
	print_mem("grad");
	if(save){
		writebin(*pgrad, "sodium_grad_in");
	}
	dcell* i02=dcellnew(1, 1);
	P(i02, 0)=dref(i0i->m);
	dcell* mtche=0;
	dcell* res=0;
	etf_t** etfs=mycalloc(nx, etf_t*);
	for(long ix=0; ix<nx; ix++){
		P(nai, ix, 1)=1;
		etfs[ix]=mketf(dtf, nac, 0, srot, srsa, hs, htel, za, 1);//no need to update
		P(nai, ix, 1)=0;
	}
	toc("mketf");tic;
	print_mem("mketf");
	int nrep=3;
	dbg("svdthres=%g, tikcr=%g, nrep=%d\n", svdthres, tikcr, nrep);
	//don't try to cachc fotf. It is per WFS and uses too much storage.
	dcell* ata=0, * atb=0;
	dccell* i0m=dccellnew(nx, 1);
	for(int irep=0; irep<nrep; irep++){
		dbg("repeat %d of %d\n", irep+1, nrep);
		dcell* i0m2=dcellnew(1, nx);
		for(long ix=0; ix<nx; ix++){
			gensei(&P(i0m, ix), NULL, NULL, NULL, sepsf, dtf, etfs[ix], saa, radgx?srot:NULL, siglev, wvlwts, *pgrad, 0, 0);
			P(i0m2, 0, ix)=dref(P(i0m, ix)->m);
		}
		dcellzero(ata);
		dcellzero(atb);
		dcellmm(&ata, i0m2, i0m2, "tn", 1);
		dcellmm(&atb, i0m2, i02, "tn", 1);
		cellfree(i0m2);

		dcellsvd_pow(ata, -1, svdthres, tikcr);
		//dcell *pi0m=dcellpinv2(i0m2, NULL, svdthres, tikcr);
		dcellzero(res);
		dcellmm(&res, ata, atb, "nn", 1);
		real scale=1./dcellsum(res);//make sure nai sum to 1.
		for(long ix=0; ix<nx; ix++){
			P(nai, ix, 1)=P(P(res, ix), 0)*scale;
		}
		etf_t* etfi=mketf(dtf, nac, 0, srot, srsa, hs, htel, za, 1);
		gensei(pi0, pgx, pgy, NULL, sepsf, dtf, etfi, saa, radgx?srot:NULL, siglev, wvlwts, *pgrad, 0, 0);
		etf_free(etfi);
		mtch_cell(&mtche, NULL, NULL, NULL, *pi0, *pgx, *pgy, NULL, NULL, NULL, 0, 0, 3,
			pixthetax, pixthetay, NULL, radgx, 1, 1);
		real gmax=0;
		for(long ii0=0; ii0<ni0; ii0++){
			for(long isa=0; isa<nsa; isa++){
				real g[2]={0,0};
				dmulvec(g, P(mtche, isa, ii0), P(P(i0i, isa, ii0)), 1);
				P(P(*pgrad, ii0), isa)+=g[0];
				P(P(*pgrad, ii0), isa+nsa)+=g[1];
				gmax=MAX(MAX(gmax, fabs(g[0])), fabs(g[1]));
			}
		}

		info("gradient has maximum %g mas\n", gmax*206265000);
		if(save){
			writebin(ata, "sodium_ata_%d", irep);
			writebin(atb, "sodium_atb_%d", irep);
			writebin(nai, "sodium_prof_%d", irep);
			writebin(*pgrad, "sodium_grad_%d", irep);
			writebin(mtche, "sodium_mtche_%d", irep);
			writebin(*pi0, "sodium_i0_%d", irep);
		}
	}
	dcellfree(ata);
	dcellfree(atb);
	cellfree(i0m);
	for(long ix=0; ix<nx; ix++){
		etf_free(etfs[ix]);
	}
	free(etfs);
	cellfree(mtche);
	cellfree(i02);
	cellfree(res);
	cellfree(nac);

	cellfree(i0);
	cellfree(gx);
	cellfree(gy);
	cellfree(grad);
}
