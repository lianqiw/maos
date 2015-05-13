/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
#include "perf.h"

int   *cuperf_t::nembed=0;
int   *cuperf_t::psfsize=0;
Real *cuperf_t::wvls=0;    
cudaStream_t    *cuperf_t::stream=0;
cublasHandle_t  *cuperf_t::handle=0;
cufftHandle     *cuperf_t::plan=0;

curcell *cuperf_t::surf=0;
curcell *cuperf_t::opd=0;
curcell *cuperf_t::psfcl=0;
curcell *cuperf_t::psfcl_ngsr=0;
curcell *cuperf_t::opdcov=0;
curcell *cuperf_t::opdcov_ngsr=0;
curcell *cuperf_t::opdmean=0;
curcell *cuperf_t::opdmean_ngsr=0;
curcell *cuperf_t::cc_cl=0;
curcell *cuperf_t::cc_ol=0;
curcell *cuperf_t::coeff=0;
Real **cuperf_t::ccb_cl=0;
Real **cuperf_t::ccb_ol=0;

/**
   Initialize perfevl
*/
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper){
    if(!parms->gpu.evl){
	return;
    }
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    /*The following lives in CPU. */
    if(parms->evl.psfmean || parms->evl.psfhist){
	cuperf_t::nembed =(int*)  calloc(nwvl, sizeof(int));
	cuperf_t::psfsize=(int*)  calloc(nwvl, sizeof(int));
	cuperf_t::wvls   =(Real*)calloc(nwvl, sizeof(Real));
    
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cuperf_t::nembed[iwvl]=(int)aper->embed->nembed->p[iwvl];
	    cuperf_t::psfsize[iwvl]=parms->evl.psfsize->p[iwvl];
	    cuperf_t::wvls[iwvl]=parms->evl.wvl->p[iwvl];
	}
    }
    /*The following lives in GPU. */
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->perf=new cuperf_t;
	pthread_mutex_init(&cudata->perf->mutex, 0);
	cudata->perf->locs=new culoc_t(aper->locs);
	cp2gpu(&cudata->perf->amp, aper->amp);
	cp2gpu(&cudata->perf->imcc, aper->imcc);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cudata->perf->embed    = (int**) calloc(nwvl, sizeof(int*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cp2gpu(&cudata->perf->embed[iwvl], aper->embed->embed->p[iwvl]->p, aper->locs->nloc, 1);
	    }
	}
    }/*for igpu */
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	if(!cudata->perf->locs_dm){
	    cudata->perf->locs_dm=new culoc_t**[nevl];
	}
	cudata->perf->locs_dm[ievl]=new culoc_t*[parms->ndm];
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_t *loc_dm;
	    if(aper->locs_dm && aper->locs_dm->p[ievl+idm*nevl]){
		loc_dm=aper->locs_dm->p[ievl+idm*nevl];
	    }else{
		loc_dm=aper->locs;
	    }
	    cudata->perf->locs_dm[ievl][idm]=new culoc_t(loc_dm);
	}
    }
    cuperf_t::stream=(cudaStream_t*)calloc(nevl, sizeof(cudaStream_t));
    cuperf_t::handle=(cublasHandle_t*)calloc(nevl, sizeof(cublasHandle_t));
    if(parms->evl.psfmean || parms->evl.psfhist){
	cuperf_t::plan  = (cufftHandle*)calloc(nwvl*nevl, sizeof(cufftHandle));
    }
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	STREAM_NEW(cuperf_t::stream[ievl]);
	cublasCreate(&cuperf_t::handle[ievl]);
	cublasSetStream(cuperf_t::handle[ievl], cuperf_t::stream[ievl]);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(iwvl>0 && cuperf_t::nembed[iwvl]==cuperf_t::nembed[0]){
		    cuperf_t::plan[iwvl+nwvl*ievl]=cuperf_t::plan[0+nwvl*ievl];
		}else{
		    DO(cufftPlan2d(&cuperf_t::plan[iwvl+nwvl*ievl],cuperf_t::nembed[iwvl],
				   cuperf_t::nembed[iwvl],FFT_T_C2C));
		    DO(cufftSetStream(cuperf_t::plan[iwvl+nwvl*ievl], cuperf_t::stream[ievl]));
		}
	    }/*for iwvl */
	}
    }
    cuperf_t::opd=curcellnew(nevl,1);
    cuperf_t::cc_cl=curcellnew(nevl, 1);
    cuperf_t::cc_ol=curcellnew(nevl, 1);
    cuperf_t::coeff=curcellnew(nevl, 1);
    cuperf_t::ccb_ol=(Real**)malloc(sizeof(Real*)*nevl);
    cuperf_t::ccb_cl=(Real**)malloc(sizeof(Real*)*nevl);
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	cuperf_t::ccb_ol[ievl]=(Real*)malloc4async(sizeof(Real)*7);
	cuperf_t::ccb_cl[ievl]=(Real*)malloc4async(sizeof(Real)*7);
	cuperf_t::cc_cl->p[ievl]=curnew(7, 1);
	cuperf_t::cc_ol->p[ievl]=curnew(7, 1);
	cuperf_t::coeff->p[ievl]=curnew(7, 1);
	cuperf_t::opd->p[ievl]=curnew(aper->locs->nloc, 1);
    }
    if(!parms->sim.evlol){
	if(parms->evl.cov && parms->gpu.psf){
	    cuperf_t::opdcov=curcellnew(nevl, 1);
	    cuperf_t::opdmean=curcellnew(nevl, 1);
	    cuperf_t::opdcov_ngsr=curcellnew(nevl, 1);
	    cuperf_t::opdmean_ngsr=curcellnew(nevl, 1);
	}
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cuperf_t::psfcl = curcellnew(nwvl, parms->evl.nevl);
	    cuperf_t::psfcl_ngsr = curcellnew(nwvl, parms->evl.nevl);
	}
    }
    if(aper->opdadd){
	cuperf_t::surf=curcellnew(nevl, 1);
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(cudata_t::evlgpu[ievl]);
	    cp2gpu(&cuperf_t::surf->p[ievl], aper->opdadd->p[ievl]);
	}
    }
    gpu_print_mem("perf init");
}
/*
  Initialize simulation data. Seed dependent. Create for the first seed and zero for the next.
*/
void gpu_perfevl_init_sim(const PARMS_T *parms, APER_T *aper){
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    int nloc=aper->locs->nloc;
    if(!parms->gpu.evl){
	return;
    }
    /*first open loop ones are on every GPU.*/
    if(parms->evl.psfol){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(parms->evl.cov && parms->gpu.psf){ /*do OL opd cov*/
		initzero(&cudata->perf->opdcovol, nloc, nloc);
		initzero(&cudata->perf->opdmeanol, nloc, 1);
	    }
	    if(parms->evl.psfmean || parms->evl.psfhist){
		if(cudata->perf->psfol){
		    curcellzero(cudata->perf->psfol);
		}else{
		    cudata->perf->psfol=curcellnew(nwvl,1);
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cudata->perf->psfol->p[iwvl]=curnew(cudata->perf->psfsize[iwvl], 
							    cudata->perf->psfsize[iwvl]);
		    }
		}
	    }
	}
    }

    if(parms->evl.cov && parms->gpu.psf && !parms->sim.evlol){
	for(int ievl=0; ievl<nevl; ievl++){
	    if(parms->evl.psf->p[ievl]==0){
		continue;
	    }
	    gpu_set(cudata_t::evlgpu[ievl]);
	    if(parms->evl.psfngsr->p[ievl]){
		initzero(&cuperf_t::opdcov_ngsr->p[ievl], nloc,nloc);
		initzero(&cuperf_t::opdmean_ngsr->p[ievl], nloc,1);
	    }
	    if(parms->evl.psfngsr->p[ievl]!=2){
		initzero(&cuperf_t::opdcov->p[ievl],nloc,nloc);
		initzero(&cuperf_t::opdmean->p[ievl],nloc,1);
	    }
	}
    }
	
    if(parms->evl.psfmean || parms->evl.psfhist){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(!cudata->perf->wvf){/*temporary. no need to zero*/
		cudata->perf->wvf=cuccellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    if(!parms->evl.psfhist && iwvl>0 && cuperf_t::nembed[iwvl] == cuperf_t::nembed[iwvl-1]){
			cudata->perf->wvf->p[iwvl]=cucref(cudata->perf->wvf->p[iwvl-1]);
		    }else{
			cudata->perf->wvf->p[iwvl]=cucnew(cuperf_t::nembed[iwvl], cuperf_t::nembed[iwvl]);
		    }
		}
	    }
	    if(parms->evl.psfmean){
		cudata->perf->psfs=cuccellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    if(cuperf_t::psfsize[iwvl]<cuperf_t::nembed[iwvl]){
			cudata->perf->psfs->p[iwvl]=cucnew(cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
		    }else{
			cudata->perf->psfs->p[iwvl]=cucref(cudata->perf->wvf->p[iwvl]);
		    }
		}
	    }
	}
	if(!parms->sim.evlol){
	    for(int ievl=0; ievl<nevl; ievl++){
		if(parms->evl.psf->p[ievl]==0){
		    continue;
		}
		gpu_set(cudata_t::evlgpu[ievl]);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    if(parms->evl.psfngsr->p[ievl]){
			initzero(&cuperf_t::psfcl_ngsr->p[iwvl+nwvl*ievl], 
				 cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
		    }
		    if(parms->evl.psfngsr->p[ievl]!=2){
			initzero(&cuperf_t::psfcl->p[iwvl+nwvl*ievl],
				 cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
		    }
		}	
	    }
	}
    }
    CUDA_SYNC_DEVICE;
}
