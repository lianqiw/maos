/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	cuglobal->perf.nembed =Array<int>(nwvl, 1);
	cuglobal->perf.psfsize=Array<int>(nwvl, 1);
	cuglobal->perf.wvls   =Array<Real>(nwvl, 1);
    
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cuglobal->perf.nembed[iwvl]=(int)aper->embed->nembed->p[iwvl];
	    cuglobal->perf.psfsize[iwvl]=parms->evl.psfsize->p[iwvl];
	    cuglobal->perf.wvls[iwvl]=parms->evl.wvl->p[iwvl];
	}
    }
    /*The following lives in GPU. */
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->perf.locs=culoc_t(aper->locs);
	cp2gpu(cudata->perf.amp, aper->amp);
	cp2gpu(cudata->perf.imcc, aper->imcc);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cudata->perf.embed    = Array<cuimat>(nwvl, 1);//(int**) calloc(nwvl, sizeof(int*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cp2gpu(cudata->perf.embed[iwvl], aper->embed->embed->p[iwvl]->p, aper->locs->nloc, 1);
	    }
	}
    }/*for igpu */
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cuglobal->evlgpu[ievl]);
	if(!cudata->perf.locs_dm){
	    cudata->perf.locs_dm=Array<Array<culoc_t> >(nevl, 1);
	}
	cudata->perf.locs_dm[ievl]=Array<culoc_t>(parms->ndm,1);
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_t *loc_dm;
	    if(aper->locs_dm && aper->locs_dm->p[ievl+idm*nevl]){
		loc_dm=aper->locs_dm->p[ievl+idm*nevl];
	    }else{
		loc_dm=aper->locs;
	    }
	    cudata->perf.locs_dm[ievl][idm]=culoc_t(loc_dm);
	}
    }
    //cuglobal->perf.stream=Array<stream_t>(nevl, 1);
    if(parms->evl.psfmean || parms->evl.psfhist){
	cuglobal->perf.plan  = Array<cufftHandle>(nwvl*nevl,1);
    }
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cuglobal->evlgpu[ievl]);
	//STREAM_NEW(cuglobal->perf.stream[ievl]);
	//Use stream created per GPU in order to share resource within GPU between different evl dir.
	//cuglobal->perf.stream[ievl]=cudata->perf_stream;
	if(parms->evl.psfmean || parms->evl.psfhist){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		DO(cufftPlan2d(&cuglobal->perf.plan[iwvl+nwvl*ievl],cuglobal->perf.nembed[iwvl],
			       cuglobal->perf.nembed[iwvl],FFT_T_C2C));
		DO(cufftSetStream(cuglobal->perf.plan[iwvl+nwvl*ievl], cudata->perf_stream));
	    }/*for iwvl */
	}
    }
    cuglobal->perf.nevl=nevl;
    cuglobal->perf.opd=curcell(nevl,1);
    cuglobal->perf.cc_cl=curcell(nevl, 1);
    cuglobal->perf.cc_ol=curcell(nevl, 1);
    cuglobal->perf.coeff=curcell(nevl, 1);
    cuglobal->perf.ccb_ol.init(nevl,1);
    cuglobal->perf.ccb_cl.init(nevl,1);
    //cuglobal->perf.ccb_ol=(Real**)malloc(sizeof(Real*)*nevl);
    //cuglobal->perf.ccb_cl=(Real**)malloc(sizeof(Real*)*nevl);
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cuglobal->evlgpu[ievl]);
	cuglobal->perf.ccb_ol[ievl].init(7,1);//=(Real*)malloc4async(sizeof(Real)*7);
	cuglobal->perf.ccb_cl[ievl].init(7,1);//=(Real*)malloc4async(sizeof(Real)*7);
	cuglobal->perf.cc_cl[ievl].init(7,1);
	cuglobal->perf.cc_ol[ievl].init(7,1);
	cuglobal->perf.coeff[ievl].init(7,1);
	cuglobal->perf.opd[ievl].init(aper->locs->nloc, 1);
    }
    if(!parms->sim.evlol){
	if(parms->evl.cov && parms->gpu.psf){
	    cuglobal->perf.opdcov=curcell(nevl, 1);
	    cuglobal->perf.opdmean=curcell(nevl, 1);
	    cuglobal->perf.opdcov_ngsr=curcell(nevl, 1);
	    cuglobal->perf.opdmean_ngsr=curcell(nevl, 1);
	}
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cuglobal->perf.psfcl = curcell(nwvl, parms->evl.nevl);
	    cuglobal->perf.psfcl_ngsr = curcell(nwvl, parms->evl.nevl);
	}
    }
    if(aper->opdadd){
	cuglobal->perf.surf=curcell(nevl, 1);
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(cuglobal->evlgpu[ievl]);
	    cp2gpu(cuglobal->perf.surf[ievl], aper->opdadd->p[ievl]);
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
		initzero(cudata->perf.opdcovol, nloc, nloc);
		initzero(cudata->perf.opdmeanol, nloc, 1);
	    }
	    if(parms->evl.psfmean || parms->evl.psfhist){
		if(cudata->perf.psfol){
		    cuzero(cudata->perf.psfol);
		}else{
		    cudata->perf.psfol=curcell(nwvl,1);
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cudata->perf.psfol[iwvl].init(cuglobal->perf.psfsize[iwvl], 
						      cuglobal->perf.psfsize[iwvl]);
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
	    gpu_set(cuglobal->evlgpu[ievl]);
	    if(parms->evl.psfngsr->p[ievl]){
		initzero(cuglobal->perf.opdcov_ngsr[ievl], nloc,nloc);
		initzero(cuglobal->perf.opdmean_ngsr[ievl], nloc,1);
	    }
	    if(parms->evl.psfngsr->p[ievl]!=2){
		initzero(cuglobal->perf.opdcov[ievl],nloc,nloc);
		initzero(cuglobal->perf.opdmean[ievl],nloc,1);
	    }
	}
    }
	
    if((parms->evl.psfmean || parms->evl.psfhist) && !parms->sim.evlol){
	for(int ievl=0; ievl<nevl; ievl++){
	    if(parms->evl.psf->p[ievl]==0){
		continue;
	    }
	    gpu_set(cuglobal->evlgpu[ievl]);
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(parms->evl.psfngsr->p[ievl]){
		    initzero(cuglobal->perf.psfcl_ngsr[iwvl+nwvl*ievl], 
			     cuglobal->perf.psfsize[iwvl], cuglobal->perf.psfsize[iwvl]);
		}
		if(parms->evl.psfngsr->p[ievl]!=2){
		    initzero(cuglobal->perf.psfcl[iwvl+nwvl*ievl],
			     cuglobal->perf.psfsize[iwvl], cuglobal->perf.psfsize[iwvl]);
		}
	    }	
	}
    }
}
