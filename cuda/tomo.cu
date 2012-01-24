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
#include "wfs.h"
#include "recon.h"
#include "accphi.h"
#include "cucmat.h"
#define SYNC_PS  for(int ips=0; ips<recon->npsr; ips++){ cudaStreamSynchronize(curecon->psstream[ips]); }
#define SYNC_WFS  for(int iwfs=0; iwfs<nwfs; iwfs++){ cudaStreamSynchronize(curecon->wfsstream[iwfs]); }

#define DIM_REDUCE 128 /*dimension to use in reduction. */

#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
#define ctoc(A)
#else
#define ctoc(A) CUDA_SYNC_STREAM; toc2(A);tic
#endif
__global__ static void ptt_proj_do(float *restrict out, float (*restrict PTT)[2], float *restrict grad, int ng){
    __shared__ float gx[DIM_REDUCE];
    __shared__ float gy[DIM_REDUCE];
    gx[threadIdx.x]=0;
    gy[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x;
    for(int ig=blockIdx.x * blockDim.x + threadIdx.x; ig<ng; ig+=step){/*ng is nsa*2. */
	gx[threadIdx.x]+=PTT[ig][0]*grad[ig];
	gy[threadIdx.x]+=PTT[ig][1]*grad[ig];
    }
    __syncthreads();
    for(int step=(DIM_REDUCE>>1); step>0; step>>=1){
	if(threadIdx.x<step){
	    gx[threadIdx.x]+=gx[threadIdx.x+step];
	    gy[threadIdx.x]+=gy[threadIdx.x+step];
	}
	__syncthreads();
    }
    if(threadIdx.x==0){
	atomicAdd(&out[0], -gx[0]);
	atomicAdd(&out[1], -gy[0]);
    }
}

/**
   Multiply nea to gradients inplace.
*/
__global__ static void gpu_tt_nea_do(float *restrict g, const float (*neai)[3], const float *const tt, const int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	float gx=g[isa]+tt[0]; 
	float gy=g[isa+nsa]+tt[1];
	float cx=neai[isa][0];
	float cy=neai[isa][1];
	float cxy=neai[isa][2];
	g[isa]=cx*gx+cxy*gy;
	g[isa+nsa]=cxy*gx+cy*gy;
    }
}

/**
   apply GP';
*/
__global__ static void gpu_gpt_o2_fuse_do(float *restrict map, int nx, const float *restrict g, 
					  const float (*neai)[3], const float *const tt,
					  const int (*restrict saptr)[2], float dsa, 
					  float *px, float *py, int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	float cx=neai[isa][0];
	float cy=neai[isa][1];
	float cxy=neai[isa][2];
	float gx=g[isa]+tt[0];
	float gy=g[isa+nsa]+tt[1];
	float tmp=cxy*gx;
	gx=cx*gx+cxy*gy;
	gy=tmp+cy*gy;
	float *restrict px2=px+isa*9;
	float *restrict py2=py+isa*9;
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	atomicAdd(&map[iy    *nx+ix],   gx*px2[0] + gy*py2[0]);
	atomicAdd(&map[iy    *nx+ix+1], gx*px2[1] + gy*py2[1]);
	atomicAdd(&map[iy    *nx+ix+2], gx*px2[2] + gy*py2[2]);
	atomicAdd(&map[(iy+1)*nx+ix],   gx*px2[3] + gy*py2[3]);
	atomicAdd(&map[(iy+1)*nx+ix+1], gx*px2[4] + gy*py2[4]);
	atomicAdd(&map[(iy+1)*nx+ix+2], gx*px2[5] + gy*py2[5]);
	atomicAdd(&map[(iy+2)*nx+ix],   gx*px2[6] + gy*py2[6]);
	atomicAdd(&map[(iy+2)*nx+ix+1], gx*px2[7] + gy*py2[7]);
	atomicAdd(&map[(iy+2)*nx+ix+2], gx*px2[8] + gy*py2[8]);
    }
}

/**
   apply GP; Has some difference due to edge effects from gradients.
*/
__global__ static void gpu_gp_o2_fuse_do(const float *restrict map, int nx, float *restrict g, 
					 const int (*restrict saptr)[2], float dsa, 
					 float *px, float *py, int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	
	float *restrict px2=px+isa*9;
	float *restrict py2=py+isa*9;
	g[isa]=
	    +map[iy*nx+ix  ]*px2[0]
	    +map[iy*nx+ix+1]*px2[1]
	    +map[iy*nx+ix+2]*px2[2]
	    +map[(iy+1)*nx+ix ] *px2[3]
	    +map[(iy+1)*nx+ix+1]*px2[4]
	    +map[(iy+1)*nx+ix+2]*px2[5]
	    +map[(iy+2)*nx+ix ] *px2[6]
	    +map[(iy+2)*nx+ix+1]*px2[7]
	    +map[(iy+2)*nx+ix+2]*px2[8];
	g[isa+nsa]=
	    +map[iy*nx+ix  ]*py2[0]
	    +map[iy*nx+ix+1]*py2[1]
	    +map[iy*nx+ix+2]*py2[2]
	    +map[(iy+1)*nx+ix ] *py2[3]
	    +map[(iy+1)*nx+ix+1]*py2[4]
	    +map[(iy+1)*nx+ix+2]*py2[5]
	    +map[(iy+2)*nx+ix ] *py2[6]
	    +map[(iy+2)*nx+ix+1]*py2[7]
	    +map[(iy+2)*nx+ix+2]*py2[8];
    }/*for isa */
}



#define DO_HX								\
    const float hs = parms->powfs[ipowfs].hs;				\
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);			\
    for(int ips=0; ips<recon->npsr; ips++){				\
	const float ht=recon->ht->p[ips];				\
	const float scale = 1.f - ht/hs;				\
	const float oxx=recon->xmap[ips]->ox;				\
	const float oyx=recon->xmap[ips]->oy;				\
	float dispx=parms->wfsr[iwfs].thetax*ht;			\
	float dispy=parms->wfsr[iwfs].thetay*ht;			\
	if(parms->tomo.predict){					\
	    int ips0=parms->atmr.indps[ips];				\
	    dispx+=simu->atm[ips0]->vx*simu->dt*2;			\
	    dispy+=simu->atm[ips0]->vy*simu->dt*2;			\
	}								\
	gpu_prop_grid(opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, \
		      xin->p[ips], oxx, oyx,recon->xmap[ips]->dx,	\
		      dispx, dispy,					\
		      1.f, 'n', curecon->wfsstream[iwfs]);		\
    }/*for ips*/

#define DO_HXT								\
    const float ht=recon->ht->p[ips];					\
    const float oxx=recon->xmap[ips]->ox;				\
    const float oyx=recon->xmap[ips]->oy;				\
    for(int iwfs=0; iwfs<nwfs; iwfs++){					\
	const int ipowfs = parms->wfsr[iwfs].powfs;			\
	if(parms->powfs[ipowfs].skip) continue;				\
	const float hs = parms->powfs[ipowfs].hs;			\
	const float scale = 1.f - ht/hs;				\
	float dispx=parms->wfsr[iwfs].thetax*ht;			\
	float dispy=parms->wfsr[iwfs].thetay*ht;			\
	if(parms->tomo.predict){					\
	    int ips0=parms->atmr.indps[ips];				\
	    dispx+=simu->atm[ips0]->vx*simu->dt*2;			\
	    dispy+=simu->atm[ips0]->vy*simu->dt*2;			\
	}								\
	gpu_prop_grid(opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, \
		      opdx->p[ips], oxx, oyx,recon->xmap[ips]->dx,	\
		      dispx, dispy,					\
		      alpha, 't', curecon->psstream[ips]);		\
    }

/*
  With ptt_proj_do, it is much faster than curmm.  Now adding ptt is
  twice as slow as the projection because grad is accessed as write
  after read. Don't split proj and add_ptt to two loops. Depth wise is
  actually faster. Final timing: 0.260 ms per call.
*/

#define DO_GP								\
    if(cupowfs[ipowfs].GP){						\
	curzero(grad->p[iwfs], curecon->wfsstream[iwfs]);		\
	cuspmul(grad->p[iwfs]->p, cupowfs[ipowfs].GP, opdwfs->p[iwfs]->p, 1.f, curecon->wfssphandle[iwfs]); \
	gpu_tt_nea_do<<<DIM(nsa,256),0,curecon->wfsstream[iwfs]>>>	\
	    (grad->p[iwfs]->p, (float(*)[3])curecon->neai->p[iwfs]->p, &ttf->p[iwfs*2], nsa); \
    }else{								\
	gpu_gp_o2_fuse_do<<<DIM(nsa,64),0,curecon->wfsstream[iwfs]>>>	\
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cupowfs[ipowfs].saptr, dsa, \
	     cupowfs[ipowfs].GPpx->p, cupowfs[ipowfs].GPpy->p, nsa);	\
    } 

#define DO_PTT								\
    if(ptt && curecon->PTT && curecon->PTT->p[iwfs+iwfs*nwfs]){		\
	/*Using ptt_proj_do is much faster than using curmm.*/		\
    	ptt_proj_do<<<DIM(nsa*2, DIM_REDUCE), 0, curecon->wfsstream[iwfs]>>> \
	    (&ttf->p[iwfs*2], (float(*)[2])curecon->PTT->p[iwfs+iwfs*nwfs]->p, grad->p[iwfs]->p, nsa*2); \
    }									

#define DO_NEA_GPT								\
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);			\
    if(cupowfs[ipowfs].GP){						\
	cusptmul(opdwfs->p[iwfs]->p, cupowfs[ipowfs].GP, grad->p[iwfs]->p, 1.f, curecon->wfssphandle[iwfs]); \
    }else{								\
	gpu_gpt_o2_fuse_do<<<DIM(nsa,64),0,curecon->wfsstream[iwfs]>>> \
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, \
	     (float(*)[3])curecon->neai->p[iwfs]->p, &ttf->p[iwfs*2], cupowfs[ipowfs].saptr, dsa, \
	     cupowfs[ipowfs].GPpx->p, cupowfs[ipowfs].GPpy->p, nsa);	\
    }									

__global__ static void laplacian_do(float *restrict out, const float *in, int nx, int ny, const float alpha){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int nx1=nx-1;
    int ny1=ny-1;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	int iy1 =iy+1; if(iy1>ny1) iy1-=ny;
	int iy2 =iy+2; if(iy2>ny1) iy2-=ny;
	int iy1p=iy-1; if(iy1p<0) iy1p+=ny;
	int iy2p=iy-2; if(iy2p<0) iy2p+=ny;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    int ix1 =ix+1; if(ix1>nx1) ix1-=nx;
	    int ix2 =ix+2; if(ix2>nx1) ix2-=nx;
	    int ix1p=ix-1; if(ix1p<0) ix1p+=nx;
	    int ix2p=ix-2; if(ix2p<0) ix2p+=nx;
	    out[ix+iy*nx]+=alpha*(20.f*in[ix+iy*nx]
				  -8.f*(in[ix1p+iy*nx]+in[ix1+iy*nx]+in[ix+iy1p*nx]+in[ix+iy1*nx])
				  +2.f*(in[ix1p+iy1p*nx]+in[ix1+iy1p*nx]+in[ix1p+iy1*nx]+in[ix1+iy1*nx])
				  +(in[ix+iy2p*nx]+in[ix2p+iy*nx]+in[ix2+iy*nx]+in[ix+iy2*nx]));
	}
    }
}
__global__ static void zzt_do(float *restrict out, const float *in, int ix, float val){
    if(threadIdx.x==0 && threadIdx.y==0){
	out[ix]+=in[ix]*val;
    }
}
static inline curcell *new_xout(const RECON_T *recon){
    /*curcell *xout=curcellnew(recon->npsr, 1);
    for(int ips=0; ips<recon->npsr; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	xout->p[ips]=curnew(nxo, nyo);
	}*/
    curcell *xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    return xout;
}
void gpu_TomoR(curcell **xout, const void *A, curcell *grad, const float alpha){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(!*xout){
	*xout=new_xout(recon);
    }
    curcell *opdx=*xout;
    int ptt=1;
    const int nwfs=grad->nx;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const int nxp=recon->pmap->nx;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
    curmat *ttf=curnew(2, nwfs);
    cuwloc_t *cupowfs=cudata->powfs;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	const int nsa=cupowfs[ipowfs].nsa;
	if(parms->powfs[ipowfs].skip) continue;
	const float dsa=cupowfs[ipowfs].dsa;
	DO_PTT;
	/*curwrite(grad->p[iwfs], "grad_nea_%d", iwfs); */
	DO_NEA_GPT;
	/*curwrite(opdwfs->p[iwfs], "gpt_%d", iwfs); */
    }
    SYNC_WFS;
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;

    }
    SYNC_PS;
    curfree(ttf);
    toc("TomoR");
}


void gpu_TomoL(curcell **xout, const float beta, const void *A, const curcell *xin, const float alpha){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int nwfs=parms->nwfsr;
    const int nps=recon->npsr;
    curcell *grad=curecon->grad;
    if(!*xout){
	*xout=new_xout(recon);
    }
    curcell *opdx=*xout;
    const int nxp=recon->pmap->nx;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
    int ptt=(!parms->recon.split || parms->dbg.splitlrt); 
    curmat *ttf=curnew(2, nwfs);
    cuwloc_t *cupowfs=cudata->powfs;
#if TIMING==2
    static cudaEvent_t (*wfsevent)[5]=NULL;
    static cudaEvent_t (*psevent)[3]=NULL;
    if(!wfsevent){
	wfsevent=(cudaEvent_t (*)[5])calloc(5*nwfs, sizeof(cudaEvent_t));
	psevent =(cudaEvent_t (*)[3])calloc(3*nps, sizeof(cudaEvent_t));
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    for(int it=0; it<5; it++){
		DO(cudaEventCreate(&wfsevent[iwfs][it]));
	    }
	}
	for(int ips=0; ips<nps; ips++){
	    for(int it=0; it<3; it++){
		DO(cudaEventCreate(&psevent[ips][it]));
	    }
	}
    }
#endif
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const int nsa=cupowfs[ipowfs].nsa;
	const float dsa=cupowfs[ipowfs].dsa;
#if TIMING==2
	cudaEventRecord(wfsevent[iwfs][0], curecon->wfsstream[iwfs]);
	DO_HX;
	cudaEventRecord(wfsevent[iwfs][1], curecon->wfsstream[iwfs]);
	DO_GP;
	cudaEventRecord(wfsevent[iwfs][2], curecon->wfsstream[iwfs]);
	DO_PTT;
	cudaEventRecord(wfsevent[iwfs][3], curecon->wfsstream[iwfs]);
	DO_NEA_GPT;
	cudaEventRecord(wfsevent[iwfs][4], curecon->wfsstream[iwfs]);
#else
	DO_HX;
	DO_GP;
	DO_PTT;
	DO_NEA_GPT;
#endif
    }
    SYNC_WFS;
    /*#endif */
   
    for(int ips=0; ips<nps; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
#if TIMING == 2
	cudaEventRecord(psevent[ips][0], curecon->psstream[ips]);
	curscale((*xout)->p[ips], beta, curecon->psstream[ips]);
	laplacian_do<<<DIM2(nxo, nyo, 16), 0, curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, nxo, nyo, curecon->l2c[ips]*alpha);
	if(parms->tomo.piston_cr){
	    zzt_do<<<1,1,0,curecon->psstream[ips]>>>
		(opdx->p[ips]->p, xin->p[ips]->p, curecon->zzi[ips], curecon->zzv[ips]*alpha);
	}
	cudaEventRecord(psevent[ips][1], curecon->psstream[ips]);
	DO_HXT;
	cudaEventRecord(psevent[ips][2], curecon->psstream[ips]);
#else
	curscale((*xout)->p[ips], beta, curecon->psstream[ips]);
	laplacian_do<<<DIM2(nxo, nyo, 16), 0, curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, nxo, nyo, curecon->l2c[ips]*alpha);
	if(parms->tomo.piston_cr){
	    zzt_do<<<1,1,0,curecon->psstream[ips]>>>
		(opdx->p[ips]->p, xin->p[ips]->p, curecon->zzi[ips], curecon->zzv[ips]*alpha);
	}
	DO_HXT;
#endif
    }
    SYNC_PS;
    curfree(ttf);
#if TIMING==2
    static char *wfstimc[]={"TomoL:HX", "TomoL:GP", "TomoL:PTT", "TomoL:NEA_GPT"};
    static int count=0;
    static float wfstim[4]={0};
    count++;
    float tmp;
    for(int it=0; it<4; it++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    DO(cudaEventElapsedTime(&tmp, wfsevent[iwfs][it], wfsevent[iwfs][it+1]));
	    wfstim[it]+=tmp;
	}
	info2("%-20s takes %.3f ms\n", wfstimc[it], wfstim[it]/count);
    }
    static char *pstimc[]={"TomoL:L2", "TomoL:HXT"};
    static float pstim[2]={0};
    for(int it=0; it<2; it++){
	for(int ips=0; ips<nps; ips++){
	    DO(cudaEventElapsedTime(&tmp, psevent[ips][it], psevent[ips][it+1]));
	    pstim[it]+=tmp;
	}
	info2("%-20s takes %.3f ms\n", pstimc[it], pstim[it]/count);
    }
#elif TIMING ==1
    toc("TomoL:");
#endif
}
/* embed real to complex data.*/
__global__ static void fdpcg_embed(fcomplex *out, float *in, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=make_cuComplex(in[ix], 0);
    }
}
/* extract real from complex data.*/
__global__ static void fdpcg_extract(float *out, fcomplex *in, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=cuCrealf(in[ix]);
    }
}
__global__ static void fdpcg_perm(fcomplex *out, fcomplex *in, int *perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=in[perm[ix]];
    }
}
__global__ static void fdpcg_perm_i(fcomplex *out, fcomplex *in, int *perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[perm[ix]]=in[ix];
    }
}
/* Each thread block is bsxbs*/
/*
__global__ static void fdpcg_mul_block_slow(fcomplex *xout, fcomplex *xin, fcomplex *M, int nx){
    extern __shared__ fcomplex v[];
    int ib=blockIdx.x;
    int bs=blockDim.x;
    int ix=threadIdx.x;
    int iy=threadIdx.y;
    fcomplex *vin=v;
    fcomplex *vout=v+bs;
    xin+=ib*bs;
    xout+=ib*bs;
    M+=ib*bs*bs;
    if(iy==0){
	vin[ix]=xin[ix];
	vout[ix]=make_cuComplex(0,0);
    }
    __syncthreads();
    fcomplex tmp=cuCmulf(M[ix+iy*bs],vin[iy]);
    atomicAdd((float*)(&vout[ix]), cuCrealf(tmp));
    atomicAdd(((float*)(&vout[ix]))+1, cuCimagf(tmp));
    __syncthreads();
    if(iy==0){
	xout[ix]=vout[ix];
    }
    }*/
/* Each thread block is bs x 1. no need to lock.*/
__global__ static void fdpcg_mul_block_sync(fcomplex *xout, fcomplex *xin, fcomplex *M, int nx){
    extern __shared__ fcomplex v[];
    int ib=blockIdx.x;
    int bs=blockDim.x;
    int ix=threadIdx.x;
    fcomplex *vin=v;
    fcomplex *vout=v+bs;
    xin+=ib*bs;
    xout+=ib*bs;
    M+=ib*bs*bs;
    vin[ix]=xin[ix];
    vout[ix]=make_cuComplex(0,0);
    __syncthreads();
    for(int iy=0; iy<bs; iy++){
	vout[ix]=cuCfmaf(M[ix+iy*bs], vin[iy], vout[ix]);
    }
    xout[ix]=vout[ix];
}

/* Each thread block is bs x 1. no need to lock. unroll the loop.*/
template <int bs>
__global__ static void fdpcg_mul_block(fcomplex *xout, fcomplex *xin, fcomplex *M, int nx){
    int ib=blockIdx.x;
    int ix=threadIdx.x;
    __shared__ fcomplex vin[bs];
    __shared__ fcomplex vout[bs];
    xin+=ib*bs;
    xout+=ib*bs;
    M+=ib*bs*bs;
    vin[ix]=xin[ix];
    vout[ix]=make_cuComplex(0,0);
#pragma unroll
    for(int iy=0; iy<bs; iy++){
	vout[ix]=cuCfmaf(M[ix+iy*bs], vin[iy], vout[ix]);
    }
    xout[ix]=vout[ix];
}
__global__ static void fdpcg_scale(fcomplex *x, float alpha, int nx){
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nx; i+=step){
	x[i]=make_cuComplex(cuCrealf(x[i])*alpha, cuCimagf(x[i])*alpha);
    }
}
/**
   Timing for full reconstruction (tomo+fit)
   0.0238 with mul_block.
   0.0224 with mul_block_nosync
   0.0225 with mul_block_sync

   Detailed timing here (in micro-second). os4n means not unrolled. others are unrolled.
   Done in orion with nfiraos with 6 layer.
                 os4n os4  os1   os6   os2
   Copy:          72   72   42    93    52
   FFT+scale:    152  144   89   165   109
   Permutation:   90   86   38   162    54
   Multiply:     964  563  232   568   288
   Inverse Perm: 161  156   47   242    77
   Inverse FFT:  148  140   83   158   107
   Copy back:     44   38   27    43    38
   Total:       1619 1199  558  1431   725
   Total FDPCG3:     6982 5760  7434  5759 (excluding right hand size)
   Try to use the hermitian property when trying to optimize for the speed. 
   Be careful about oversampling:

   The positive (right half) base frequency (os=1) couples to both positive and
   negative frequencies to layers with os=2. Only apply the block matrix to
   positive frequencies maynot be good.

*/
void gpu_Tomo_fdprecond(curcell **xout, const void *A, const curcell *xin, cudaStream_t stream){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    RECON_T *recon=simu->recon;
    if(!xin->m){
	error("xin is not continuous");
    }
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }else if(!(*xout)->m){
	error("xout is not continuous");
    }
    fdpcg_embed<<<DIM(curecon->fd_nxtot, 256),0,stream>>>
	(curecon->fd_xhat1->p[0]->p, xin->p[0]->p, curecon->fd_nxtot);
    ctoc("fdpcg: Copy");
    for(int ic=0; ic<curecon->fd_fftnc; ic++){
	int ips=curecon->fd_fftips[ic];
	CUFFT(curecon->fd_fft[ic], curecon->fd_xhat1->p[ips]->p, CUFFT_FORWARD);
	int nps=curecon->fd_fftips[ic+1]-curecon->fd_fftips[ic];
	int nx=curecon->fd_xhat1->p[ips]->nx;
	int ny=curecon->fd_xhat1->p[ips]->ny;
	fdpcg_scale<<<DIM(nps*nx*ny, 256), 0, stream>>>
	    (curecon->fd_xhat1->p[ips]->p, 1.f/sqrtf((float)(nx*ny)), nps*nx*ny);
    }
    ctoc("fdpcg: FFT+scale");
    fdpcg_perm<<<DIM(curecon->fd_nxtot, 256),0,stream>>>
	(curecon->fd_xhat2->p[0]->p, curecon->fd_xhat1->p[0]->p, curecon->fd_perm, curecon->fd_nxtot);
    ctoc("fdpcg: Permutation");
    int nb=curecon->fd_Mb->nx;
    int bs=curecon->fd_Mb->p[0]->nx;

#define MUL(N)								\
    fdpcg_mul_block<N><<<nb, N, 0, stream>>>				\
	(curecon->fd_xhat1->p[0]->p,curecon->fd_xhat2->p[0]->p,		\
	 curecon->fd_Mb->p[0]->p, curecon->fd_nxtot)
    switch(bs){
    case 1: MUL(1); break;
    case 2: MUL(2); break;
    case 3: MUL(3); break;
    case 4: MUL(4); break;
    case 5: MUL(5); break;
    case 6: MUL(6); break;
    case 7: MUL(7); break;
    case 8: MUL(8); break;
    case 9: MUL(9); break;
    case 10: MUL(10); break;
    case 11: MUL(11); break;
    case 12: MUL(12); break;
    case 13: MUL(13); break;
    case 14: MUL(14); break;
    case 15: MUL(15); break;
    case 16: MUL(16); break;
    case 17: MUL(17); break;
    case 18: MUL(18); break;
    case 19: MUL(19); break;
    case 20: MUL(20); break;
    case 21: MUL(21); break;
    case 22: MUL(22); break;
    case 23: MUL(23); break;
    case 24: MUL(24); break;
    case 25: MUL(25); break;
    case 26: MUL(26); break;
    case 27: MUL(27); break;
    case 28: MUL(28); break;
    case 29: MUL(29); break;
    case 30: MUL(30); break;
    case 31: MUL(31); break;
    default:
	fdpcg_mul_block_sync<<<nb, bs, sizeof(fcomplex)*bs*2, stream>>>
	    (curecon->fd_xhat1->p[0]->p,curecon->fd_xhat2->p[0]->p, curecon->fd_Mb->p[0]->p, curecon->fd_nxtot);
    }
    ctoc("fdpcg: mul block");
    fdpcg_perm_i<<<DIM(curecon->fd_nxtot, 256),0,stream>>>
	(curecon->fd_xhat2->p[0]->p, curecon->fd_xhat1->p[0]->p, curecon->fd_perm, curecon->fd_nxtot);
    ctoc("fdpcg: Inverse Permutation");
    for(int ic=0; ic<curecon->fd_fftnc; ic++){
	int ips=curecon->fd_fftips[ic];
	CUFFT(curecon->fd_fft[ic], curecon->fd_xhat2->p[ips]->p, CUFFT_INVERSE);
	int nps=curecon->fd_fftips[ic+1]-curecon->fd_fftips[ic];
	int nx=curecon->fd_xhat1->p[ips]->nx;
	int ny=curecon->fd_xhat1->p[ips]->ny;
	fdpcg_scale<<<DIM(nps*nx*ny, 256), 0, stream>>>
	    (curecon->fd_xhat2->p[ips]->p, 1.f/sqrtf((float)(nx*ny)), nps*nx*ny);
    }
    ctoc("fdpcg: Inverse FFT + scale");
    fdpcg_extract<<<DIM(curecon->fd_nxtot, 256),0,stream>>>
	((*xout)->p[0]->p, curecon->fd_xhat2->p[0]->p, curecon->fd_nxtot);
    ctoc("fdpcg: Copy back");
}
