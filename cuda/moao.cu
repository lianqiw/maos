/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "recon.h"
#include "accphi.h"
#include "pcg.h"
#define CALL_ONCE\
    do {static int count=0; count++; if(count>1) warning("This function should only be called once\n");} while(0)
extern int *wfsgpu;
extern int *evlgpu;
void gpu_setup_moao(const PARMS_T *parms, RECON_T *recon){
    gpu_set(gpu_recon);
    if(!curecon){
	curecon=(curecon_t*)calloc(1, sizeof(curecon_t));
    }
    curecon->moao=(cumoao_t*)calloc(parms->nmoao,sizeof(cumoao_t));
    for(int imoao=0; imoao<parms->nmoao; imoao++){
	cumoao_t* cumoao=&curecon->moao[imoao];
	gpu_dcell2cu(&cumoao->fitNW,recon->moao[imoao].NW);
	gpu_spcell2dev(&cumoao->actslave, recon->moao[imoao].actslave);
	if(parms->moao[imoao].cubic){
	    cumoao->cubic_cc=gpu_dmcubic_cc(parms->moao[imoao].iac);
	}
	cumoao->dxa=recon->moao[imoao].amap->dx;
	cumoao->oxa=recon->moao[imoao].amap->ox;
	cumoao->oya=recon->moao[imoao].amap->oy;
	cumoao->nxa=recon->moao[imoao].amap->nx;
	cumoao->nya=recon->moao[imoao].amap->ny;
	
	cumoao->dxf=recon->fmap->dx;
	cumoao->oxf=recon->fmap->ox;
	cumoao->oyf=recon->fmap->oy;
	cumoao->nxf=recon->fmap->nx;
	cumoao->nyf=recon->fmap->ny;
	if(curecon->W01){
	    cumoao->W01=curecon->W01;/*same configuration.*/
	}else{
	    cumoao->W01=gpu_get_W01(recon->W0, recon->W1);
	}
	cumoao->rhs=curcellnew(1,1);
	cumoao->rhs->p[0]=curnew(cumoao->nxa, cumoao->nya);
	cudaMalloc(&cumoao->pis, sizeof(float));
	cumoao->xp=curnew(cumoao->nxf, cumoao->nyf);
	cumoao->xp2=curnew(cumoao->nxf, cumoao->nyf);
	if(cumoao->fitNW){
	    cumoao->tmpNW=curnew(cumoao->fitNW->p[0]->nx, 1);
	}
    }
    if(!curecon->moao_wfs){
	int nwfs=parms->nwfs;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0) continue;
	    cumoao_t *cumoao=&curecon->moao[imoao];
	    if(!curecon->moao_wfs){
		curecon->moao_wfs=(curcell**)calloc(nwfs, sizeof(curcell*));
	    }
	    curecon->moao_wfs[iwfs]=curcellnew(1,1);
	    curecon->moao_wfs[iwfs]->p[0]=new cumap_t(cumoao->nxa, cumoao->nya, 
						      cumoao->oxa, cumoao->oya,
						      cumoao->dxa, 0, 0, 0);
	    if(parms->gpu.wfs){
		gpu_set(wfsgpu[iwfs]);
	    }
	    if(!cudata->moao_wfs){
		cudata->moao_wfs=curcellnew(nwfs, 1);
	    }
	    if(parms->sim.closeloop){
		cudata->moao_wfs->p[iwfs]=new cumap_t(cumoao->nxa, cumoao->nya, 
						      cumoao->oxa, cumoao->oya,
						      cumoao->dxa, 0, 0, 0);
	    }else{
		cudata->moao_wfs->p[iwfs]=curecon->moao_wfs[iwfs]->p[0]->ref(); 
	    }
	    gpu_set(gpu_recon);
	}
    }
    if(!curecon->moao_evl && parms->evl.moao!=-1){
	const int imoao=parms->evl.moao;
	const int nevl=parms->evl.nevl;
	curecon->moao_evl=(curcell**)calloc(nevl, sizeof(curcell*));
	cumoao_t *cumoao=&curecon->moao[imoao];
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(gpu_recon);
	    curecon->moao_evl[ievl]=curcellnew(1,1);
	    curecon->moao_evl[ievl]->p[0]=new cumap_t(cumoao->nxa, cumoao->nya,
						      cumoao->oxa, cumoao->oya,
						      cumoao->dxa, 0, 0, 0);
	    if(parms->gpu.evl){
		gpu_set(evlgpu[ievl]);
	    }
	    if(!cudata->moao_evl){
		cudata->moao_evl=curcellnew(nevl,1);
	    }
	    if(parms->sim.closeloop){
		cudata->moao_evl->p[ievl]=new cumap_t(cumoao->nxa, cumoao->nya,
						      cumoao->oxa, cumoao->oya,
						      cumoao->dxa, 0, 0, 0);
	    }else{
		cudata->moao_evl->p[ievl]=curecon->moao_evl[ievl]->p[0]->ref();
	    }
	}
    }
    gpu_set(gpu_recon);
    STREAM_NEW(curecon->moao_stream);
    HANDLE_NEW(curecon->moao_handle, curecon->moao_stream);
    SPHANDLE_NEW(curecon->moao_sphandle, curecon->moao_stream);
}

#define DO_W								\
    gpu_inn(cumoao->pis, xp->p, cumoao->W01->W1->p, np, stream);	\
    add2_do<<<DIM(np, 256), 0, stream>>>(xp2->p, cumoao->W01->W1->p, cumoao->pis, -1.f, np); \
    cuspmul(xp2->p, cumoao->W01->W0p, xp->p, 1.f, sphandle);		\
    if(cumoao->W01->nW0f){						\
	apply_W_do<<<DIM(np, 256),0,stream>>>(xp2->p, xp->p, cumoao->W01->W0f, cumoao->W01->W0v, \
					      cumoao->nxf, cumoao->W01->nW0f); \
    }

#define DO_HAT /*Apply HAT, from xp2 to xout.*/				\
    if(!*xout) *xout=curcellnew(1,1);					\
    if(!(*xout)->p[0]) (*xout)->p[0]=curnew(cumoao->nxa, cumoao->nya);	\
    if(cumoao->cubic_cc){						\
	gpu_prop_grid_cubic(xp2, cumoao->oxf,cumoao->oyf, cumoao->dxf,	\
			    (*xout)->p[0], cumoao->oxa, cumoao->oya, cumoao->dxa, \
			    0,0, cumoao->cubic_cc, alpha, 't', stream);	\
    }else{								\
	gpu_prop_grid(xp2, cumoao->oxf, cumoao->oyf, cumoao->dxf,	\
		      (*xout)->p[0], cumoao->oxa, cumoao->oya, cumoao->dxa, \
		      0,0, alpha, 't', stream);				\
    }

/*Right hand size vector.*/
void gpu_moao_FitR(curcell **xout, SIM_T *simu, cumoao_t *cumoao, float thetax, float thetay, float hs, const float alpha){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int npsr=recon->npsr;
    const int np=cumoao->nxf*cumoao->nyf;
    cudaStream_t stream=curecon->moao_stream;
    cusparseHandle_t sphandle=curecon->moao_sphandle;
    curmat *xp=cumoao->xp;
    curmat *xp2=cumoao->xp2;
    curzero(xp, stream);
    curzero(xp2, stream);

    /*do HX operation, from curecon->opdr to xp. */
    if(!*xout) *xout=curcellnew(1,1);
    curzero((*xout)->p[0], stream);
    for(int ips=0; ips<npsr; ips++){
	const float ht = (float)recon->ht->p[ips];
	const float scale=1.f-ht/hs;
	const float dispx=thetax*ht;
	const float dispy=thetay*ht;
	gpu_prop_grid(xp, cumoao->oxf*scale, cumoao->oyf*scale, cumoao->dxf*scale, 
		      curecon->opdr->p[ips], recon->xmap[ips]->ox, recon->xmap[ips]->oy,
		      recon->xmap[ips]->dx,
		      dispx, dispy, 1.f,'n', stream);
    }
    /*do HA operation, from curecon->dmfit to xp */
    for(int idm=0; idm<recon->ndm; idm++){
	const float ht = (float)parms->dm[idm].ht;
	const float scale=1.f-ht/hs;
	const float dispx=thetax*ht;
	const float dispy=thetay*ht;
	if(curecon->cubic_cc[idm]){
	    gpu_prop_grid_cubic
		(xp, cumoao->oxf*scale,cumoao->oyf*scale, cumoao->dxf*scale, 
		 curecon->dmfit->p[idm], recon->amap[idm]->ox, recon->amap[idm]->oy, recon->amap[idm]->dx,
		 dispx, dispy, curecon->cubic_cc[idm], -1.f, 'n', stream);
	}else{
	    gpu_prop_grid
		(xp, cumoao->oxf*scale, cumoao->oyf*scale, cumoao->dxf*scale, 
		 curecon->dmfit->p[idm], recon->amap[idm]->ox, recon->amap[idm]->oy, recon->amap[idm]->dx,
		 dispx, dispy, -1.f, 'n', stream);
	}
    }
    CUDA_SYNC_DEVICE;
    /*apply W, from xp to xp2*/
    DO_W;
    /*do HAT operation, from xp2 to xout*/
    DO_HAT;
    CUDA_SYNC_STREAM;
}

void gpu_moao_FitL(curcell **xout, const void *A, const curcell *xin, const float alpha){
    cumoao_t *cumoao=(cumoao_t*)A;
    cudaStream_t stream=curecon->moao_stream;
    cublasHandle_t handle=curecon->moao_handle;
    cusparseHandle_t sphandle=curecon->moao_sphandle;

    const int np=cumoao->nxf*cumoao->nyf;
    curmat *xp=cumoao->xp;
    curmat *xp2=cumoao->xp2;
    curzero(xp, stream);
    curzero(xp2, stream);
    /*Apply HA, from xin to xp.*/
    if(cumoao->cubic_cc){
	gpu_prop_grid_cubic(xp, cumoao->oxf,cumoao->oyf, cumoao->dxf, 
			    xin->p[0], cumoao->oxa, cumoao->oya, cumoao->dxa,
			    0,0, cumoao->cubic_cc, 1.f, 'n', stream);
    }else{
	gpu_prop_grid(xp, cumoao->oxf, cumoao->oyf, cumoao->dxf, 
		      xin->p[0], cumoao->oxa, cumoao->oya, cumoao->dxa,
		      0,0, 1.f, 'n', stream);	
    }
    /*Apply W, from xp to xp2*/
    DO_W;
    /*Apply Hat, from xp2 to xout*/
    DO_HAT;
    /*Additional terms*/
    if(cumoao->fitNW){
	curmat *tmp=cumoao->tmpNW;
	curmv(tmp->p, 0, cumoao->fitNW->p[0], xin->p[0]->p, 't', 1, handle);
	curmv((*xout)->p[0]->p, 1, cumoao->fitNW->p[0], tmp->p, 'n', alpha, handle);
    }
    if(cumoao->actslave){
	cuspmul((*xout)->p[0]->p, cumoao->actslave->p[0], xin->p[0]->p, alpha, sphandle);
    }
    CUDA_SYNC_STREAM;
}
/**
   MOAO reconstruction.

   Do not output directly to cudata->moao since that is being used by wfsgrad
and perfevl. The new result is supposed to be used next time step. The input
based on opdr, dmfit is on gradients from last time step. So two cycle delay is
maintained.  */
void gpu_moao_recon(SIM_T *simu){
    gpu_set(gpu_recon);
    curcell *dmcommon=NULL;
    const PARMS_T *parms=simu->parms;
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    if(parms->gpu.fit){
	dmcommon=curecon->dmfit;
    }else{
	gpu_dcell2cu(&dmcommon, simu->dmfit_hi);
    }
    cudaStream_t stream=curecon->moao_stream;
    if(curecon->moao_wfs){/*There is MOAO DM for WFS */
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0){
		continue;
	    }
	    cumoao_t *cumoao=&curecon->moao[imoao];
	    gpu_moao_FitR(&cumoao->rhs, simu, cumoao,
			  parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, 
			  parms->powfs[ipowfs].hs, 1);
	    if(gpu_pcg(&curecon->moao_wfs[iwfs], gpu_moao_FitL, cumoao, NULL, NULL, cumoao->rhs,
		       simu->parms->recon.warm_restart, parms->fit.maxit, stream)){
		error("PCG failed\n");
	    }
	}
    }
    if(curecon->moao_evl){
	int imoao=parms->evl.moao;
	if(imoao==-1) error("Inconsistent\n");
	cumoao_t *cumoao=&curecon->moao[imoao];
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_moao_FitR(&cumoao->rhs, simu, cumoao,
			  parms->evl.thetax[ievl], parms->evl.thetay[ievl], 
			  parms->evl.hs[ievl], 1);
	    if(gpu_pcg(&curecon->moao_evl[ievl], gpu_moao_FitL, cumoao, NULL, NULL, cumoao->rhs,
		       simu->parms->recon.warm_restart, parms->fit.maxit, stream)){
		error("PCG failed\n");
	    }
	}
    }
    CUDA_SYNC_STREAM;
    if(dmcommon!=curecon->dmfit){
	curcellfree(dmcommon);
    }
}
void gpu_moao_filter(SIM_T *simu){
    gpu_set(gpu_recon);
    const PARMS_T *parms=simu->parms;
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    if(curecon->moao_wfs){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0) continue;
	    double g=parms->moao[imoao].gdm;
	    curmat *temp=NULL;
	    if(parms->gpu.wfs && wfsgpu[iwfs]!=gpu_recon){//copy between GPUs
		gpu_set(wfsgpu[iwfs]);
		curcp(&temp, curecon->moao_wfs[iwfs]->p[0]);
	    }else{
		temp=curecon->moao_wfs[iwfs]->p[0];
	    }
	    /*use 0 as stream because this is gpu specific*/
	    if(parms->sim.closeloop){
		curadd(&cudata->moao_wfs->p[iwfs], 1.-g, temp, g, 0);
	    }
	    if(!parms->gpu.wfs ){
		gpu_cur2d(&simu->moao_wfs->p[iwfs], 0, temp, 1, 0);
	    }
	    cudaStreamSynchronize(0);
	    if(temp!=curecon->moao_wfs[iwfs]->p[0]){
		curfree(temp);
	    }
	}
    }
    if(curecon->moao_evl){
	int imoao=parms->evl.moao;
	double g=parms->moao[imoao].gdm;
	for(int ievl=0; ievl<nevl; ievl++){
	    curmat *temp=NULL;
	    if(parms->gpu.evl && evlgpu[ievl]!=gpu_recon){//copy between GPUs
		gpu_set(evlgpu[ievl]);
		curcp(&temp, curecon->moao_evl[ievl]->p[0]);
	    }else{
		gpu_set(gpu_recon);
		temp=curecon->moao_evl[ievl]->p[0];
	    }
	    /*use 0 as stream because this is gpu specific*/
	    if(parms->sim.closeloop){
		curadd(&cudata->moao_evl->p[ievl], 1.-g, temp, g, 0);
	    }
	    if(!parms->gpu.evl){
		gpu_cur2d(&simu->moao_evl->p[ievl], 0, temp, 1, 0);
	    }
	    cudaStreamSynchronize(0);
	    if(temp!=curecon->moao_evl[ievl]->p[0]){
		curfree(temp);
	    }
	}
    }
    gpu_set(gpu_recon);
}
/**
   Copy MOAO DM commands from CPU to GPU.*/
void gpu_moao_2gpu(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->gpu.moao){
	error("Invalid use\n");
    }
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    if(parms->gpu.wfs && simu->moao_wfs){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0) continue;
	    gpu_set(wfsgpu[iwfs]);
	    if(!cudata->moao_wfs){
		cudata->moao_wfs=curcellnew(nwfs, 1);
	    }
	    if(!cudata->moao_wfs->p[iwfs]){
		double dxa=recon->moao[imoao].amap->dx;
		double oxa=recon->moao[imoao].amap->ox;
		double oya=recon->moao[imoao].amap->oy;
		int nxa=recon->moao[imoao].amap->nx;
		int nya=recon->moao[imoao].amap->ny;
		cudata->moao_wfs->p[iwfs]=new cumap_t(nxa, nya,
						      oxa, oya,
						      dxa, 0, 0, 0);
	    }
	    gpu_dmat2cu(&cudata->moao_wfs->p[iwfs], simu->moao_wfs->p[iwfs]);
	}
    }
    if(parms->gpu.evl && simu->moao_evl){
	int imoao=parms->evl.moao;
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(evlgpu[ievl]);
	    if(!cudata->moao_evl){
		cudata->moao_evl=curcellnew(nevl, 1);
	    }
	    if(!cudata->moao_evl->p[ievl]){
		double dxa=recon->moao[imoao].amap->dx;
		double oxa=recon->moao[imoao].amap->ox;
		double oya=recon->moao[imoao].amap->oy;
		int nxa=recon->moao[imoao].amap->nx;
		int nya=recon->moao[imoao].amap->ny;
		cudata->moao_evl->p[ievl]=new cumap_t(nxa, nya,
						      oxa, oya,
						      dxa, 0, 0, 0);
	    }
	    gpu_dmat2cu(&cudata->moao_evl->p[ievl], simu->moao_evl->p[ievl]);
	}
    }
}
