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
#include "utils.h"
#include "accphi.h"
#include "cuwfs.h"
#include "cudata.h"

/**
   This file is work in progress and is not used.
 */
namespace cuda_wfs{
cuwfs_info::cuwfs_dbg(const parms_t* parms, const powfs_t* powfs, int _iwfs, int _igpu)
	:iwfs(_iwfs), igpu(_igpu){
	dbg("cuwfs_info[%d]\n", iwfs);
	int ipowfs=parms->wfs[iwfs].powfs;
	int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	loc=new culoc_t(powfs[ipowfs].loc);
	cp2gpu(&amp, powfs[ipowfs].realamp->p[wfsind]);
	cp2gpu(&embed, powfs[ipowfs].embed, powfs[ipowfs].loc->nloc, 1);
	embednx=embedny=powfs[ipowfs].nembed;
}

void cufieldstop_t::apply(cuwfs_info* wfsinfo, curmat& opd, cudaStream_t stream){
	if(nwvl>1) error("Please implemente\n");
	embed_wvf_do<<<DIM(opd->nx, 256), 0, stream>>>
		(wvf->p, opd->p, wfsinfo->amp->p, wfsinfo->embed->p, opd->nx, wvl[0]);
	CUFFT(plan, (Comp*)wvf->p, CUFFT_FORWARD);
	cwm_do<<<DIM(wvf->nx*wvf->ny, 256), 0, stream>>>
		(wvf->p, fieldmask->p, wvf->nx*wvf->ny);
	CUFFT(plan, (Comp*)wvf->p, CUFFT_INVERSE);
	unwrap_phase_do<<<DIM2(wvf->nx, wvf->ny, 16), 0, stream>>>
		(wvf->p, opd->p, wfsinfo->embed->p, opd->nx, wvl[0]);
}
__global__ static void setup_rand(curandState* rstat, int seed){
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	curand_init(seed, id, 0, &rstat[id]);
}
curand_t::curand_t(int seed, int nsa){
#define RAND_BLOCK 16
#define RAND_THREAD 32
	int ng=nsa*2;
	if(nsa<RAND_THREAD){
		nt=ng;
		nb=1;
	} else if(nsa<RAND_THREAD*RAND_BLOCK){
		nt=RAND_THREAD;
		nb=(ng+RAND_THREAD-1)/RAND_THREAD;
	} else{
		nt=RAND_THREAD;
		nb=RAND_BLOCK;
	}
	cudaMalloc(&rstat, nt*nb*sizeof(curandState));
	setup_rand<<<nb, nt>>>(rstat, seed);
}
struct wfscfg_t{
	const parms_t* parms;
	const powfs_t* powfs;
	int iwfs;
	int ipowfs;
	int wfsind;
	stream_t& stream;
	wfscfg_t(const parms_t* _parms, const powfs_t* _powfs, int _iwfs, stream_t& _stream)
		:parms(_parms), powfs(_powfs), iwfs(_iwfs), ipowfs(parms->wfs[_iwfs].powfs),
		wfsind(parms->powfs[ipowfs].wfsind->p[_iwfs]), stream(_stream){};

};
cuwfs_base::cuwfs_base(wfscfg_t* wfscfg):
	stream(wfscfg->stream), rand(0),
	dtrat(wfscfg->parms->powfs[wfscfg->ipowfs].dtrat){}
cushwfs_t::cushwfs_t(wfscfg_t* wfscfg)
	:cuwfs_base(wfscfg), nsa(wfscfg->powfs[wfscfg->ipowfs].pts->nsa),
	dxsa(wfscfg->powfs[wfscfg->ipowfs].pts->dx),
	nxsa(wfscfg->powfs[wfscfg->ipowfs].pts->nx){}
cushgeom_t::cushgeom_t(wfscfg_t* wfscfg)
	:nea(0), gradacc(0), cushwfs_t(wfscfg){
	const int ipowfs=wfscfg->ipowfs;
	const int wfsind=wfscfg->wfsind;
	const parms_t* parms=wfscfg->parms;
	const powfs_t* powfs=wfscfg->powfs;
	cp2gpu(&nea, powfs[ipowfs].neasim->p[wfsind]);
	gradacc=curnew(nsa*2, 1);
	if(!parms->powfs[ipowfs].pistatout
		||parms->powfs[ipowfs].pistatstc
		||parms->powfs[ipowfs].dtrat==1){
		gradcalc=curef(gradacc);
	} else{
		gradcalc=curnew(nsa*2, 1);
	}
}
void cushgeom_t::output(){
	error("To implement\n");
	cuzero(gradacc);
}
__global__ void add_geom_noise_do(Real* restrict g, const Real* restrict nea,
	int nsa, curandState* restrict rstat);
void cushgeom_t::addnoise(){
	add_geom_noise_do<<<rand->nb, rand->nt, 0, stream>>>
		(gradacc->p, nea->p, nsa, rand->rstat);
}
void cushgeom_t::acc(curmat* opd){
	if(gradcalc->p!=gradacc->p){
		calcg(opd, 1);
		curadd(&gradacc, 1, gradcalc, 1.f/(Real)dtrat, stream);
	} else{
		calcg(opd, 1.f/(Real)dtrat);
	}
	error("To implement\n");
}
cushg_t::cushg_t(wfscfg_t* wfscfg)
	:cushgeom_t(wfscfg), GS0(0){
	const powfs_t* powfs=wfscfg->powfs;
	const int ipowfs=wfscfg->ipowfs;
	const int wfsind=wfscfg->wfsind;
	GS0=new cusp(powfs[ipowfs].GS0->p[powfs[ipowfs].GS0->nx>1?wfsind:0], 1);
}
void cushg_t::calcg(curmat& opd, Real ratio){
	cuspmul(gradcalc->p, GS0, opd->p, 1, 'n', ratio, stream);
}
cushz_t::cushz_t(wfscfg_t* wfscfg)
	:cushgeom_t(wfscfg), imcc(0){
	const powfs_t* powfs=wfscfg->powfs;
	const int ipowfs=wfscfg->ipowfs;
	const int wfsind=wfscfg->wfsind;
	void* _imcc[nsa];
	for(int isa=0; isa<nsa; isa++){
		_imcc[isa]=NULL;
		cp2gpu((Real**)&(_imcc[isa]), powfs[ipowfs].saimcc->p[powfs[ipowfs].saimcc->nx>1?wfsind:0]->p[isa]);
	}
	cp2gpu(&imcc, _imcc, nsa, 1);
}
void cushz_t::calcg(curmat& opd, Real ratio){
	error("need to implement. Decide where to put pts info.\n");
	/*cuztilt<<<nsa, dim3(16,16), 0, stream>>>
	(gradcalc->p, opd->p, opd->nx, dxsa, nxsa, imcc,
	 cupowfs[ipowfs].pts->p, cuwfs[iwfs].amp,
	 ratio);*/
}
cullt_t::cullt_t(wfscfg_t* wfscfg)
	:pts(0), loc(0), amp(0), ncpa(0), imcc(0){
	const parms_t* parms=wfscfg->parms;
	const powfs_t* powfs=wfscfg->powfs;
	const int ipowfs=wfscfg->ipowfs;
	const int wfsind=wfscfg->wfsind;
	pts=new cupts_t(powfs[ipowfs].llt->pts);
	loc=new culoc_t(powfs[ipowfs].llt->loc);
	if(powfs[ipowfs].llt->ncpa){
		cp2gpu(&ncpa, powfs[ipowfs].llt->ncpa->p[powfs[ipowfs].llt->ncpa->nx>1?wfsind:0]);
	}
	{
		Real* _imcc[1]={0};
		cp2gpu((Real**)&_imcc[0], powfs[ipowfs].llt->imcc->p[0]);
		cp2gpu(&imcc, _imcc, 1, 1);
	}
	cp2gpu(&amp, powfs[ipowfs].llt->amp);
	{
		int nlwvf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
		int nlwvf2[2]={nlwvf, nlwvf};
		if(cufftPlanMany(&plan_wvf, 2, nlwvf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, 1)){
			error("CUFFT plan failed\n");
		}
		cufftSetStream(plan_wvf, wfscfg->stream);
		const int ncompx=powfs[ipowfs].ncompx;
		const int ncompy=powfs[ipowfs].ncompy;
		const int notf=MAX(ncompx, ncompy);
		if(notf==nlwvf){
			plan_otf=plan_wvf;
		} else{
			int notf2[2]={notf, notf};
			if(cufftPlanMany(&plan_otf, 2, notf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, 1)){
				error("CUFFT plan failed\n");
			}
			cufftSetStream(plan_otf, wfscfg->stream);
		}
	}
}
/**Concatenate identical arrays to single array*/
#define def_concat(z,c)					\
    z##mat *concat_##c##mat(c##mat **p, int nc){	\
    int n=p[0]->nx*p[0]->ny;				\
    z##mat *tmp=z##new(n, nc);				\
    for(int ic=0; ic<nc; ic++){				\
	type_convert(tmp->p+n*ic, p[ic]->p, n);		\
    }							\
    return tmp;						\
}
def_concat(z, c);
def_concat(s, d);

cushphy_t::cushphy_t(wfscfg_t* wfscfg)
	:dtf(0), srot(0), bkgrnd2(0), bkgrnd2c(0), mtche(0), i0sum(0), ints(0), pistatout(0),
	cushwfs_t(wfscfg){
	const parms_t* parms=wfscfg->parms;
	const powfs_t* powfs=wfscfg->powfs;
	const int ipowfs=wfscfg->ipowfs;
	const int wfsind=wfscfg->wfsind;
	{ /*FFT plans*/
	/*CUFFTW is row major. */
		int nwvf=powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac;/*size of fft */
		int nwvf2[2]={nwvf, nwvf};
		const int ncompx=powfs[ipowfs].ncompx;
		const int ncompy=powfs[ipowfs].ncompy;
		const int notf=MAX(ncompx, ncompy);
		int ncomp2[2]={ncompy, ncompx};
		int notf2[2]={notf, notf};
		/*limit the number of subapertures in each batch to less than 1024
		  to save memory. The speed is actually a tiny bit faster for NFIRAOS.*/
		msa=nsa>1024?((int)ceil((Real)nsa/(Real)(nsa/800))):nsa;
		if(cufftPlanMany(&plan1, 2, nwvf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, msa)){
			error("CUFFT plan failed\n");
		}
		cufftSetStream(plan1, stream);

		if(notf==nwvf){
			plan2=plan1;
		} else{
			if(cufftPlanMany(&plan2, 2, notf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, msa)){
				error("CUFFT plan failed\n");
			}
			cufftSetStream(plan2, stream);
		}
		if(notf==ncompx&&notf==ncompy){
			plan3=plan2;
		} else{
			if(cufftPlanMany(&plan3, 2, ncomp2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, msa)){
				error("CUFFT plan failed\n");
			}
			cufftSetStream(plan3, stream);
		}
	}
	/*DTF, ETF*/
	{
		int nwvl=parms->powfs[ipowfs].nwvl;
		dtf=new cuccell(nwvl, 1);
		etf=new cuccell(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			if(!powfs[ipowfs].dtf[iwvl].fused){//not fused to ETF
				int ndtf=powfs[ipowfs].dtf[iwvl].nominal->nx;
				int icol=powfs[ipowfs].dtf[iwvl].nominal->ny>1?wfsind:0;
				zmat* dtf0=concat_cmat(powfs[ipowfs].dtf[iwvl].nominal->p+ndtf*icol, ndtf);
				cp2gpu(&dtf->p[iwvl], dtf0);
				zfree(dtf0);
			}
			if(parms->powfs[ipowfs].llt){//Elongation
				cmat* (*petf)[nsa]=NULL;
				if(powfs[ipowfs].etfsim[iwvl].p1){
					petf=(cmat*(*)[nsa])powfs[ipowfs].etfsim[iwvl].p1->p;
					etfis1d=1;
				} else{
					petf=(cmat*(*)[nsa])powfs[ipowfs].etfsim[iwvl].p2->p;
					etfis1d=0;
				}
				cmat** petfi=petf[parms->powfs[ipowfs].llt->n>1?wfsind:0];
				zmat* etf0=concat_cmat(petfi, nsa);
				cp2gpu(&etf->p[iwvl], etf0);
				zfree(etf0);
			}
		}
		if(parms->powfs[ipowfs].llt){
			cp2gpu(&srot, powfs[ipowfs].srot->p[parms->powfs[ipowfs].llt->n>1?wfsind:0]);
		}
	}

	//Gradient operator
	switch(parms->powfs[ipowfs].phytypesim){
	case PTYPE_MF:{//mtche
		int icol=powfs[ipowfs].intstat->mtche->ny>1?wfsind:0;
		X(mat)* mtche0=concat_dmat(powfs[ipowfs].intstat->mtche->p+nsa*icol, nsa);
		cp2gpu(&mtche, mtche0);
		X(free)(mtche0);
		icol=(powfs[ipowfs].intstat->i0sum->ny>1?wfsind:0);
		cp2gpu(&i0sum, powfs[ipowfs].intstat->i0sum->p+nsa*icol, nsa, 1);
	}
		  break;
	default:error("Invalid phytypesim\n");
	}

	int dtrat=parms->powfs[ipowfs].dtrat;
	//Pixel averaged background and correction 
	bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	bkgrndc=parms->powfs[ipowfs].bkgrndc*bkgrnd;
	if(powfs[ipowfs].bkgrnd){
	//Pixel specific background and correction
		int icol=(powfs[ipowfs].bkgrnd->ny==1?wfsind:0);
		X(mat)* tmp=concat_dmat(powfs[ipowfs].bkgrnd->p+icol*nsa, nsa);
		cp2gpu(&bkgrnd2, tmp);X(free)(tmp);
	}
	if(powfs[ipowfs].bkgrndc){
		int icol=(powfs[ipowfs].bkgrndc->ny==1?wfsind:0);
		X(mat)* tmp=concat_dmat(powfs[ipowfs].bkgrndc->p+icol*nsa, nsa);
		cp2gpu(&bkgrnd2c, tmp);X(free)(tmp);
	}
	//runtime data
	ints=curcellnew(nsa, 1, powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
	if(parms->powfs[ipowfs].pistatout){
		const int notfx=powfs[ipowfs].ncompx;/*necessary size to build detector image. */
		const int notfy=powfs[ipowfs].ncompy;
		const int npsf=MAX(notfx, notfy);
		pistatout=curcellnew(nsa, parms->powfs[ipowfs].nwvl, npsf, npsf);
	}
}

cuwfs_t::cuwfs_t(const parms_t* parms, const powfs_t* powfs, int iwfs, int igpu)
	:gradoff(0), opdadd(0), opd(0), geom(0), phy(0), llt(0){
	gpu_set(igpu);
	wfsinfo=new cuwfs_dbg(parms, powfs, iwfs, igpu);
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	wfscfg_t wfscfg(parms, powfs, iwfs, stream);
	fieldstop=new cufieldstop_t(powfs[ipowfs].fieldstop,
		parms->powfs[ipowfs].wvl, parms->powfs[ipowfs].nwvl, stream);
	int has_geom=(!parms->powfs[ipowfs].usephy||parms->save.gradgeom->p[iwfs]||parms->powfs[ipowfs].pistatout);
	if(has_geom){
		switch(parms->powfs[ipowfs].gtype_sim){
		case GTYPE_G://gtilt
			geom=new cushg_t(&wfscfg);
			break;
		case GTYPE_Z://ztilt
			geom=new cushz_t(&wfscfg);
			break;
		default:
			error("Invalid gtype_sim\n");
		}
	}
	if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout){
		phy=new cushphy_t(&wfscfg);
		if(powfs[ipowfs].llt&&parms->powfs[ipowfs].trs){
			llt=new cullt_t(&wfscfg);
		}
	}

	if(powfs[ipowfs].opdadd&&powfs[ipowfs].opdadd->p[wfsind]){
		cp2gpu(&opdadd, powfs[ipowfs].opdadd->p[wfsind]);
	}
	if(powfs[ipowfs].gradoff&&powfs[ipowfs].gradoff->p[wfsind]){
		cp2gpu(&gradoff, powfs[ipowfs].gradoff->p[wfsind]);
	}
}
void cushphy_t::acc(curmat* opd){
	error("To impelemnt\n");
}
void cuwfs_t::initsim(){
	gpu_set(wfsinfo->igpu);
	if(geom) geom->initsim();
	if(phy) phy->initsim();
}
void cuwfs_t::seeding(int seed){
	gpu_set(wfsinfo->igpu);
	geom->seeding(seed);
	phy->seeding(seed);
}
}//namespace
void gpu_wfs_init_new(const parms_t* parms, const powfs_t* powfs){
	cuda_wfs::cuwfs_t** cuwfs;
	cuwfs=new cuda_wfs::cuwfs_t*[parms->nwfs];
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int igpu=gpu_next();
		if(NGPU>4&&igpu==gpu_recon){
			igpu=gpu_next();
		}
		cuwfs[iwfs]=new cuda_wfs::cuwfs_t(parms, powfs, iwfs, igpu);
	}
}
void gpu_wfs_initsim(const parms_t* parms, const powfs_t* powfs){
	error("obtain cuwfs from cudata\n");
	/*cuda_wfs::cuwfs_t **cuwfs;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	cuwfs[iwfs]->initsim();
	}*/
}
void gpu_wfs_seeing(const parms_t* parms, const powfs_t* powfs, rand_t* rstat){
	/*cuda_wfs::cuwfs_t **cuwfs;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	cuwfs[iwfs]->seeding(lrand(rstat));
	}*/
}
