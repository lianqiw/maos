extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "accphi.h"
#define SYNC_PS  for(int ips=0; ips<recon->npsr; ips++){ cudaStreamSynchronize(curecon->psstream[ips]); }
#define SYNC_WFS  for(int iwfs=0; iwfs<nwfs; iwfs++){ cudaStreamSynchronize(curecon->wfsstream[iwfs]); }

#define DIM_REDUCE 128 //dimension to use in reduction.

#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
#endif
__global__ static void ptt_proj_do(float *restrict out, float (*restrict PTT)[2], float *restrict grad, int ng){
    __shared__ float gx[DIM_REDUCE];
    __shared__ float gy[DIM_REDUCE];
    gx[threadIdx.x]=0;
    gy[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x;
    for(int ig=blockIdx.x * blockDim.x + threadIdx.x; ig<ng; ig+=step){//ng is nsa*2.
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
    }//for isa
}



#define DO_HX								\
    const float hs = parms->powfs[ipowfs].hs;				\
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);			\
    for(int ips=0; ips<recon->npsr; ips++){				\
	const float ht=recon->ht->p[ips];				\
	const float scale = 1.f - ht/hs;				\
	const float oxx=recon->xmap[ips]->ox;				\
	const float oyx=recon->xmap[ips]->oy;				\
	gpu_prop_grid(opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, \
		      xin->p[ips], oxx, oyx,recon->xmap[ips]->dx,	\
		      parms->wfsr[iwfs].thetax*ht, parms->wfsr[iwfs].thetay*ht, \
		      1.f, 'n', curecon->wfsstream[iwfs]);		\
    }/*for ips*/

#define DO_HXT								\
    const float ht=recon->ht->p[ips];					\
    const float oxx=recon->xmap[ips]->ox;				\
    const float oyx=recon->xmap[ips]->oy;				\
    for(int iwfs=0; iwfs<nwfs; iwfs++){				\
	const int ipowfs = parms->wfsr[iwfs].powfs;			\
	if(parms->powfs[ipowfs].skip) continue;				\
	const float hs = parms->powfs[ipowfs].hs;			\
	const float scale = 1.f - ht/hs;				\
	gpu_prop_grid(opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, \
		      opdx->p[ips], oxx, oyx,recon->xmap[ips]->dx,	\
		      parms->wfsr[iwfs].thetax*ht, parms->wfsr[iwfs].thetay*ht, \
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
	    (grad->p[iwfs]->p, (float(*)[3])curecon->neai->p[iwfs]->p, ttf->p[iwfs]->p, nsa); \
    }else{								\
	gpu_gp_o2_fuse_do<<<DIM(nsa,64),0,curecon->wfsstream[iwfs]>>>	\
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, \
	     cupowfs[ipowfs].GPpx->p, cupowfs[ipowfs].GPpy->p, nsa);	\
    } 

#define DO_PTT								\
    ttf->p[iwfs]=curnew(2,1);						\
    if(ptt && curecon->PTT && curecon->PTT->p[iwfs+iwfs*nwfs]){		\
	/*Using ptt_proj_do is much faster than using curmm.*/		\
    	ptt_proj_do<<<DIM(nsa*2, DIM_REDUCE), 0, curecon->wfsstream[iwfs]>>> \
	    (ttf->p[iwfs]->p, (float(*)[2])curecon->PTT->p[iwfs+iwfs*nwfs]->p, grad->p[iwfs]->p, nsa*2); \
    }									

#define DO_NEA_GPT								\
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);			\
    if(cupowfs[ipowfs].GP){						\
	cusptmul(opdwfs->p[iwfs]->p, cupowfs[ipowfs].GP, grad->p[iwfs]->p, 1.f, curecon->wfssphandle[iwfs]); \
    }else{								\
	gpu_gpt_o2_fuse_do<<<DIM(nsa,64),0,curecon->wfsstream[iwfs]>>> \
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, \
	     (float(*)[3])curecon->neai->p[iwfs]->p, ttf->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, \
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
    atomicAdd(&out[ix], in[ix]*val);
}
static inline curcell *new_xout(const RECON_T *recon){
    curcell *xout=curcellnew(recon->npsr, 1);
    for(int ips=0; ips<recon->npsr; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	xout->p[ips]=curnew(nxo, nyo);
    }
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
    curcell *ttf=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const float dsa=cuwfs[iwfs].powfs->dsa;
	DO_PTT;
	//curwrite(grad->p[iwfs], "grad_nea_%d", iwfs);
	DO_NEA_GPT;
	//curwrite(opdwfs->p[iwfs], "gpt_%d", iwfs);
    }
    SYNC_WFS;
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;
    }
    SYNC_PS;
    curcellfree(ttf);
    toc("TomoR");
}


void gpu_TomoL(curcell **xout, const void *A, const curcell *xin, const float alpha){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int nwfs=parms->nwfsr;
    if(!*xout){
	*xout=new_xout(recon);
    }
    curcell *opdx=*xout;
    curcell *grad=curecon->grad;
  
    for(int ips=0; ips<recon->npsr; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	laplacian_do<<<DIM2(nxo, nyo, 16), 0, curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, nxo, nyo, curecon->l2c[ips]*alpha);
	zzt_do<<<1,1,0,curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, curecon->zzi[ips], curecon->zzv[ips]*alpha);
    }
#if TIMING == 2
    toc("TomoL:L2");tic;
#endif
    const int nxp=recon->pmap->nx;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
    int ptt=(!parms->recon.split || parms->dbg.splitlrt); 
    curcell *ttf=curcellnew(nwfs, 1);
#if TIMING == 2
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	DO_HX;
    }
    SYNC_WFS;
    toc("TomoL:HX");tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const float dsa=cuwfs[iwfs].powfs->dsa;
	DO_GP;
    }
    SYNC_WFS;
    toc("TomoL:GP"); tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	DO_PTT;
    }
    SYNC_WFS;
    toc("TomoL:NEA"); tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const float dsa=cuwfs[iwfs].powfs->dsa;
	DO_NEA_GPT;
    }
    SYNC_WFS;
    toc("TomoL:GPT"); tic;
#else
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const float dsa=cuwfs[iwfs].powfs->dsa;
	DO_HX;
	DO_GP;
	DO_PTT;
	DO_NEA_GPT;
    }
    SYNC_WFS;
#endif
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;
    }
    SYNC_PS;
    curcellfree(ttf);
#if TIMING==2
    toc("TomoL:HXT");
#elif TIMING ==1
    toc("TomoL:");
#endif
}

/*void gpu_fdpcg_precond(curcell **xout, const void *A, const curcell *xin){
    SIM_T *simu=(SIM_T*)A;
    
    }*/
