extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "wfs.h"
#include "recon.h"
#include "pcg.h"
curecon_t *curecon;
#define SCALE 1 //Scale both NEA and L2 to balance the dynamic range.
#define USE_GP 0
#define SYNC_PS  for(int ips=0; ips<recon->npsr; ips++){ cudaStreamSynchronize(curecon->psstream[ips]); }
#define SYNC_WFS  for(int iwfs=0; iwfs<nwfs; iwfs++){ cudaStreamSynchronize(curecon->wfsstream[iwfs]); }

__global__ static void saloc2ptr_do(int (*restrict saptr)[2], float (*restrict saloc)[2], 
			 int nsa, float ox, float oy, float dx){
    const int step=blockDim.x * gridDim.x;
    const float dx1=1./dx;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	saptr[isa][0]=(int)roundf((saloc[isa][0]-ox)*dx1);
	saptr[isa][1]=(int)roundf((saloc[isa][1]-oy)*dx1);
    }
}

void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    curecon=(curecon_t*)calloc(1, sizeof(curecon_t));
    curecon->wfsstream=(cudaStream_t*)calloc(parms->nwfsr, sizeof(cudaStream_t));
    curecon->wfshandle=(cublasHandle_t*)calloc(parms->nwfsr, sizeof(cublasHandle_t));
    curecon->wfssphandle=(cusparseHandle_t*)calloc(parms->nwfsr, sizeof(cusparseHandle_t));
    curecon->psstream=(cudaStream_t*)calloc(recon->npsr, sizeof(cudaStream_t));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	STREAM_NEW(curecon->wfsstream[iwfs]);
	DO(cublasCreate(&curecon->wfshandle[iwfs]));
	DO(cusparseCreate(&curecon->wfssphandle[iwfs]));
	DO(cublasSetStream(curecon->wfshandle[iwfs], curecon->wfsstream[iwfs]));
	DO(cusparseSetKernelStream(curecon->wfssphandle[iwfs], curecon->wfsstream[iwfs]));
    }
    for(int ips=0; ips<recon->npsr; ips++){
	STREAM_NEW(curecon->psstream[ips]);
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].skip) continue;
	int nsa=powfs[ipowfs].pts->nsa;
	cudaMalloc(&cupowfs[ipowfs].saptr, nsa*2*sizeof(int));
	saloc2ptr_do<<<MAX(MIN(nsa/256,16),1), MIN(256, nsa)>>>
	    (cupowfs[ipowfs].saptr, cupowfs[ipowfs].saloc, nsa, 
	     recon->pmap->ox, recon->pmap->oy, recon->pmap->dx);
	gpu_sp2dev(&cupowfs[ipowfs].GP, recon->GP->p[ipowfs]);
    }
    CUDA_SYNC_DEVICE;
    cudaStream_t stream;
    STREAM_NEW(stream);
    cublasHandle_t handle;
    DO(cublasCreate(&handle));
    DO(cublasSetStream(handle, stream));
    curecon->l2c=(float*)calloc(recon->npsr, sizeof(float));
    for(int ips=0; ips<recon->npsr; ips++){
	curecon->l2c[ips]=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
	curecon->l2c[ips]=curecon->l2c[ips]*curecon->l2c[ips]*SCALE;
    }
    if(parms->tomo.piston_cr){
	curecon->zzi=(int*)calloc(recon->npsr, sizeof(int));
	curecon->zzv=(float*)calloc(recon->npsr, sizeof(float));
	for(int ips=0; ips<recon->npsr; ips++){
	    double r0=recon->r0;
	    double dx=recon->xloc[ips]->dx;
	    double wt=recon->wt->p[ips];
	    int icenter=loccenter(recon->xloc[ips]);
	    curecon->zzi[ips]=icenter;
	    curecon->zzv[ips]=pow(laplacian_coef(r0,wt,dx),2);
	}
    }
    curecon->neai=curcellnew(parms->nwfsr, 1);
    //convert recon->saneai to our format.
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	int iwfs0=parms->powfs[ipowfs].wfs[0];//first wfs in this group.
	dsp *nea=recon->saneai->p[iwfs+iwfs*parms->nwfsr];
	spint *pp=nea->p;
	spint *pi=nea->i;
	double *px=nea->x;
	float (*neai)[2]=(float(*)[2])calloc(2*nsa*2, sizeof(float));
	if(nea->n!=2*nsa) error("nea doesn't have 2nsa x 2nsa dimension\n");
	for(int ic=0; ic<nea->n; ic++){
	    for(int ir=pp[ic]; ir<pp[ic+1]; ir++){
		int ix=pi[ir];
		if(ix==ic){//diagonal part.
		    neai[ic][0]=(float)px[ir]*SCALE;
		}else if(ix==ic-nsa || ix==ic+nsa){//cross part
		    neai[ic][1]=(float)px[ir]*SCALE;
		}else{
		    error("saneai has invalid format\n");
		}
	    }
	}
	CUDA_SYNC_DEVICE;
	curecon->neai->p[iwfs]=curnew(2, 2*nsa);
	DO(cudaMemcpy(curecon->neai->p[iwfs]->p, neai, 4*nsa*sizeof(float), cudaMemcpyDefault));
	CUDA_SYNC_DEVICE;
	free(neai);
	if(iwfs!=iwfs0 && nsa>4){//don't check tt. overflows.
	    TIC;tic;
	    float diff=curinn(curecon->neai->p[iwfs], curecon->neai->p[iwfs0], stream);
	    float diff2=curinn(curecon->neai->p[iwfs0], curecon->neai->p[iwfs0],  stream);
	    toc("curinn");
	    if((diff-diff2)<1e-4*diff2){
		warning2("%d and %d has the same neai. reuse the data\n", iwfs, iwfs0);
		curfree(curecon->neai->p[iwfs]);
		curecon->neai->p[iwfs]=curecon->neai->p[iwfs0];
	    }
	}
	CUDA_SYNC_DEVICE;
    }//for iwfs
    CUDA_SYNC_DEVICE;
    if(recon->PTT && !curecon->PTT){
	gpu_dcell2cu(&curecon->PTT, recon->PTT);
    }
    if(recon->PDF && !curecon->PDF){
	gpu_dcell2cu(&curecon->PDF, recon->PDF);
    }
    STREAM_DONE(stream);
    cublasDestroy(handle);
}

__global__ static void ptt_proj_do(float *restrict out, float (*restrict PTT)[2], float *restrict grad, int ng){
    const int step=blockDim.x * gridDim.x;
    __shared__ float g[2];
    if(threadIdx.x<2){
	g[threadIdx.x]=0;
    }
    float gi[2]={0,0};
    for(int ig=blockIdx.x * blockDim.x + threadIdx.x; ig<ng; ig+=step){
	gi[0]+=PTT[ig][0]*grad[ig];
	gi[1]+=PTT[ig][1]*grad[ig];
    }
    atomicAdd(&g[0], gi[0]);
    atomicAdd(&g[1], gi[1]);
    __syncthreads();
    if(threadIdx.x<2){
	atomicAdd(&out[threadIdx.x], g[threadIdx.x]);
    }
}

__global__ static void ptt_add_do(float *restrict grad, const float *const tt, int nsa){
    const int step=blockDim.x * gridDim.x;
    float *restrict grady=grad+nsa;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nsa; i+=step){
	grad[i] +=tt[0];
	grady[i]+=tt[1];
    }
}

/**
   Multiply nea to gradients inplace.
*/
__global__ static void gpu_nea_do(float *restrict g, const float (*neai)[2], const int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	float gx=g[isa]; 
	float gy=g[isa+nsa];
	g[isa]=neai[isa][0]*gx+neai[isa+nsa][1]*gy;
	g[isa+nsa]=neai[isa][1]*gx+neai[isa+nsa][0]*gy;
    }
}

/**
   apply GP';
*/
__global__ static void gpu_gpt_o2_do(float *restrict map, int nx, const float *restrict g, 
				     const int (*restrict saptr)[2], float dsa, int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	float gx=g[isa];
	float gy=g[isa+nsa];
	float a1=gx/dsa*0.5f;
	float a2=a1*0.5f;
	atomicAdd(&map[iy*nx+ix],  -a2);
	atomicAdd(&map[iy*nx+ix+2], a2);
	atomicAdd(&map[(iy+1)*nx+ix], -a1);
	atomicAdd(&map[(iy+1)*nx+ix+2], a1);
	atomicAdd(&map[(iy+2)*nx+ix], -a2);
	atomicAdd(&map[(iy+2)*nx+ix+2], a2);
	
	a1=gy/dsa*0.5f;
	a2=a1*0.5f;

	atomicAdd(&map[iy*nx+ix],   -a2);
	atomicAdd(&map[iy*nx+ix+1], -a1);
	atomicAdd(&map[iy*nx+ix+2], -a2);
	atomicAdd(&map[(iy+2)*nx+ix],   a2);
	atomicAdd(&map[(iy+2)*nx+ix+1], a1);
	atomicAdd(&map[(iy+2)*nx+ix+2], a2);
    }
}
/**
   apply GP; Has some difference due to edge effects from gradients.
*/
__global__ static void gpu_gp_o2_do(const float *restrict map, int nx, float *restrict g, 
				     const int (*restrict saptr)[2], float dsa, int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	float a2=0.25f/dsa;
	g[isa]=(map[iy*nx+ix+2]-map[iy*nx+ix]+map[(iy+2)*nx+ix+2]-map[(iy+2)*nx+ix]
		+2.f*(map[(iy+1)*nx+ix+2]-map[(iy+1)*nx+ix]))*a2;

	g[isa+nsa]=(map[(iy+2)*nx+ix]+map[(iy+2)*nx+ix+2]-map[iy*nx+ix]-map[iy*nx+ix+2]
		    +2.f*(map[(iy+2)*nx+ix+1]-map[iy*nx+ix+1]))*a2;
    }
}
/*
  Ray tracing with matched spacing. Reverse, from out to in. out is xloc, in is ploc.
*/
__global__ static void gpu_hx_do_1(float *restrict out, int nxout,
				   const float *restrict in, int nxin, 
				   float fracx, float fracy,
				   int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    float fracx1=1.f-fracx;
    float fracy1=1.f-fracy;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    out[ix+iy*nxout]+=
		+(in[ix+    iy*nxin]*fracx1+in[ix+1+    iy*nxin]*fracx)*fracy1
		+(in[ix+(iy+1)*nxin]*fracx1+in[ix+1+(iy+1)*nxin]*fracx)*fracy;
	}
    }
}
/*
  Ray tracing with over sampling. Reverse, from out to in. out is xloc, in is
ploc. confirmed to agree with HXW'.  */
__global__ static void gpu_hxt_do_2(const float *restrict out, int nxout,
				    float *restrict in, int nxin, 
				    float fracx, float fracy,
				    int nx, int ny, float alpha){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int ax=fracx<0.5f?0:1;
    int ay=fracy<0.5f?0:1;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
	    //odd and even points are different.
#pragma unroll 
	    for(int by=0; by<2; by++){
		int iy2=iy*2+by;
		int iy3=iy+ay*by;
		float fracy2=fracy+0.5f-ay*by;
		float fracy21=1.f-fracy2;
		for(int bx=0; bx<2; bx++){
		    int ix2=ix*2+bx;
		    int ix3=ix+ax*bx;
		    float fracx2=fracx+0.5f-ax*bx;
		    float fracx21=1.f-fracx2;
		    if(ix2<nx && iy2<ny){
			float a=out[ix2+(iy2)*nxout]*alpha;
			atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
			atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
			atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
			atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
		    }
		}
	    }
	}
    }
}
/*
  Ray tracing with over sampling. Forward, from out to in. out is xloc, in is
ploc. confirmed to agree with HXW'.  */
__global__ static void gpu_hx_do_2(float *restrict out, int nxout,
				   const float *restrict in, int nxin, 
				   float fracx, float fracy,
				   int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int ax=fracx<0.5f?0:1;
    int ay=fracy<0.5f?0:1;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
	    //odd and even points are different.
#pragma unroll 
	    for(int by=0; by<2; by++){
		int iy2=iy*2+by;
		int iy3=iy+ay*by;
		float fracy2=fracy+0.5f-ay*by;
		float fracy21=1.f-fracy2;
		for(int bx=0; bx<2; bx++){
		    int ix2=ix*2+bx;
		    int ix3=ix+ax*bx;
		    float fracx2=fracx+0.5f-ax*bx;
		    float fracx21=1.f-fracx2;
		    if(ix2<nx && iy2<ny){
			out[ix2+(iy2)*nxout]+=
			    +in[ix3+    (iy3)*nxin]*fracx21*fracy21
			    +in[ix3+1+  (iy3)*nxin]*fracx2*fracy21
			    +in[ix3+  (iy3+1)*nxin]*fracx21*fracy2
			    +in[ix3+1+(iy3+1)*nxin]*fracx2*fracy2;
		    }
		}
	    }
	}
    }
}
/* Do the ray tracing when in/out grid matches in sampling. */
static inline void prop_grid_match(curmat *out, float oxo, float oyo,
			      const curmat *in, float oxi, float oyi, float dxi,
			      float dispx, float dispy,
			      float alpha, cudaStream_t stream){
    const float dx1=1.f/dxi;
    dispx=(dispx-oxi+oxo)*dx1;
    dispy=(dispy-oyi+oyo)*dx1;
    const int offx=(int)floorf(dispx); dispx-=offx;
    const int offy=(int)floorf(dispy); dispy-=offy;
    const int nxo=out->nx;
    const int nyo=out->ny;
    const int nxi=in->nx;
    const int nyi=in->ny;
    int offx1, offx2;
    int offy1, offy2;
    int nx, ny;
    if(offx>0){
	offx1=0;
	offx2=offx;
	nx=nxi-offx-1; if(nx>nxo) nx=nxo;
    }else{
	offx1=-offx;
	offx2=0;
	nx=nxo+offx; if(nx>nxi-1) nx=nxi-1;
    }
    if(offy>0){
	offy1=0;
	offy2=offy;
	ny=nyi-offy-1; if(ny>nyo) ny=nyo;
    }else{
	offy1=-offy;
	offy2=0;
	ny=nyo+offy; if(ny>nyi-1) ny=nyi-1;
    }
    gpu_hx_do_1<<<DIM2(nx,ny,16,16),0,stream>>>
	(out->p+offx1+offy1*nxo, nxo, in->p+offx2+offy2*nxi, nxi, dispx, dispy, nx, ny);
}
/* Do the ray tracing when in/out grid does not match in sampling. */
static inline void prop_grid(curmat *out, float oxo, float oyo,
			      const curmat *in, float oxi, float oyi, float dxi,
			      float dispx, float dispy,
			      float alpha, cudaStream_t stream){

}
#define DO_HX								\
    const float hs = parms->powfs[ipowfs].hs;				\
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);			\
    for(int ips=0; ips<recon->npsr; ips++){				\
	const float ht=recon->ht->p[ips];				\
	const float scale = 1.f - ht/hs;				\
	const float dx1=1./recon->xmap[ips]->dx;			\
	const float oxx=recon->xmap[ips]->ox;				\
	const float oyx=recon->xmap[ips]->oy;				\
									\
	const float ratio=dxp*scale*dx1;				\
	if(fabs(ratio-1.f)<1.e-4){					\
	    prop_grid_match(opdwfs->p[iwfs], oxp*scale, oyp*scale,	\
			    xin->p[ips], oxx, oyx,recon->xmap[ips]->dx,	\
			    parms->wfsr[iwfs].thetax*ht, parms->wfsr[iwfs].thetay*ht, \
			    1.f, curecon->wfsstream[iwfs]);		\
	}else if(fabs(ratio-0.5f)<1.e-4){				\
	    const int nxx=recon->xmap[ips]->nx;				\
	    const int nyx=recon->xmap[ips]->ny;				\
	    float dispx=(parms->wfsr[iwfs].thetax*ht-oxx+oxp*scale)*dx1; \
	    float dispy=(parms->wfsr[iwfs].thetay*ht-oyx+oyp*scale)*dx1; \
	    int offx=(int)floorf(dispx); dispx-=offx;			\
	    int offy=(int)floorf(dispy); dispy-=offy;			\
	    int nx=(nxx-offx-1)*2-((dispx<0.5f)?0:1); if(nx>nxp) nx=nxp; \
	    int ny=(nyx-offy-1)*2-((dispy<0.5f)?0:1); if(ny>nyp) ny=nyp; \
	    gpu_hx_do_2<<<DIM2((nx+1)/2,(ny+1)/2,16,16),0,curecon->wfsstream[iwfs]>>> \
		(opdwfs->p[iwfs]->p, nxp, xin->p[ips]->p+offx+offy*nxx, nxx, dispx, dispy, nx, ny); \
	}else{								\
	    TO_IMPLEMENT;						\
	}								\
    }/*for ips*/

#define DO_HXT								\
    float ht=recon->ht->p[ips];						\
    float dx1=1./recon->xmap[ips]->dx;					\
    const int nxo=recon->xmap[ips]->nx;					\
    const int nyo=recon->xmap[ips]->ny;					\
    const float oxo=recon->xmap[ips]->ox;				\
    const float oyo=recon->xmap[ips]->oy;				\
    for(int iwfs=0; iwfs<nwfs; iwfs++){					\
	const int ipowfs = parms->wfsr[iwfs].powfs;			\
	if(parms->powfs[ipowfs].skip) continue;				\
	const float hs = parms->powfs[ipowfs].hs;			\
	const float scale = 1.f - ht/hs;				\
	const float ratio=dxp*scale*dx1;				\
	if(fabs(ratio-1.f)<1.e-4){/*matched*/				\
	    prop_grid_match(opdx->p[ips], oxo, oyo,			\
			    opdwfs->p[iwfs], oxp*scale, oyp*scale, recon->pmap->dx*scale, \
			    -parms->wfsr[iwfs].thetax*ht, -parms->wfsr[iwfs].thetay*ht, \
			    1.f, curecon->psstream[ips]);		\
	}else if(fabs(ratio-0.5f)<1.e-4){/*double sampled.*/		\
	    float dispx=(parms->wfsr[iwfs].thetax*ht-oxo+oxp*scale)*dx1; \
	    float dispy=(parms->wfsr[iwfs].thetay*ht-oyo+oyp*scale)*dx1; \
	    int offx=(int)floorf(dispx); dispx-=offx;			\
	    int offy=(int)floorf(dispy); dispy-=offy;			\
	    int nx=(nxo-offx-1)*2-((dispx<0.5f)?0:1); if(nx>nxp) nx=nxp; \
	    int ny=(nyo-offy-1)*2-((dispy<0.5f)?0:1); if(ny>nyp) ny=nyp; \
	    gpu_hxt_do_2<<<DIM2((nx+1)/2,(ny+1)/2,16,8),0,curecon->psstream[ips]>>> \
		(opdwfs->p[iwfs]->p, nxp, opdx->p[ips]->p+offx+offy*nxo, nxo, dispx, dispy, nx, ny, alpha); \
	}else{								\
	    TO_IMPLEMENT;						\
	}								\
    }

/*
  With ptt_proj_do, it is much faster than curmm.  Now adding ptt is
  twice as slow as the projection because grad is accessed as write
  after read. Don't split proj and add_ptt to two loops. Depth wise is
  actually faster. Final timing: 0.260 ms per call.
*/

#define DO_NEA								\
    if(ptt && curecon->PTT && curecon->PTT->p[iwfs+iwfs*nwfs]){		\
	/*Using ptt_proj_do is much faster than using curmm.*/		\
	ttf->p[iwfs]=curnew(2,1);					\
	ptt_proj_do<<<DIM(nsa*2, 256, 32), 0, curecon->wfsstream[iwfs]>>> \
	    (ttf->p[iwfs]->p, (float(*)[2])curecon->PTT->p[iwfs+iwfs*nwfs]->p, grad->p[iwfs]->p, nsa*2); \
	curscale(ttf->p[iwfs], -1.f, curecon->wfsstream[iwfs]);		\
	ptt_add_do<<<DIM(nsa,256,16), 0, curecon->wfsstream[iwfs]>>>	\
	    (grad->p[iwfs]->p, ttf->p[iwfs]->p, nsa);			\
    }									\
    gpu_nea_do<<<DIM(nsa,256,16),0,curecon->wfsstream[iwfs]>>>		\
	(grad->p[iwfs]->p, (float(*)[2])curecon->neai->p[iwfs]->p, nsa); \
    curzero(opdwfs->p[iwfs], curecon->wfsstream[iwfs]);

#if USE_GP
#define DO_GP								\
    cuspmul(grad->p[iwfs]->p, cupowfs[ipowfs].GP, opdwfs->p[iwfs]->p, 1, curecon->wfssphandle[iwfs]);
#define DO_GPT								\
    cusptmul(opdwfs->p[iwfs]->p, cupowfs[ipowfs].GP, grad->p[iwfs]->p, 1, curecon->wfssphandle[iwfs]);
#else
#define DO_GP								\
    gpu_gp_o2_do<<<DIM(nsa,256,32),0,curecon->wfsstream[iwfs]>>>	\
	(opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, nsa);
#define DO_GPT								\
    gpu_gpt_o2_do<<<DIM(nsa,256,16),0,curecon->wfsstream[iwfs]>>>	\
	(opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, nsa);
#endif
/**
   Apply opdx=GX'*W*grad (GX=GP*HX). We use the stencil specified in NFIRAOS RTC, instead
   of using sparse matrix.  Specified for oversampling of 2 in ploc.*/
static void gpu_iprop(const PARMS_T *parms, const RECON_T *recon, curcell *opdx,
		      curcell *grad, float alpha, int ptt){
    TIC;tic;
    const int nwfs=grad->nx;
    int nxp=recon->pmap->nx;
    int nyp=recon->pmap->ny;
    float oxp=recon->pmap->ox;
    float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
   
    /*
      With ptt_proj_do, it is much faster than curmm.  Now adding ptt is
      twice as slow as the projection because grad is accessed as write
      after read. Don't split proj and add_ptt to two loops. Depth wise is
      actually faster. Final timing: 0.260 ms per call.
    */
    curcell *ttf=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
#if !USE_GP 
	const float dsa=cuwfs[iwfs].powfs->dsa;
#endif
	DO_NEA;
	DO_GPT;
    }
    SYNC_WFS;
    curcellfree(ttf);
    toc("gpu_iprop:ptt+nea+gp");tic;
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;
    }
    SYNC_PS;
    toc("gpu_iprop:hxt");
}

void gpu_TomoR(curcell **xout, const void *A, curcell *grad, const float alpha){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    curcell *opdx;
    if(!*xout){
	opdx=curcellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    const int nxo=recon->xmap[ips]->nx;
	    const int nyo=recon->xmap[ips]->ny;
	    opdx->p[ips]=curnew(nxo, nyo);
	}
	*xout=opdx;
    }else{
	opdx=*xout;
    }
    //gpu_iprop(parms, recon, opdx, grad, alpha, 1);
    int ptt=1;
    const int nwfs=grad->nx;
    int nxp=recon->pmap->nx;
    int nyp=recon->pmap->ny;
    float oxp=recon->pmap->ox;
    float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
    curcell *ttf=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
#if !USE_GP 
	const float dsa=cuwfs[iwfs].powfs->dsa;
#endif
	DO_NEA;
	DO_GPT;
    }
    SYNC_WFS;
    curcellfree(ttf);
    toc("gpu_iprop:ptt+nea+gp");tic;
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;
    }
    SYNC_PS;
    toc("gpu_iprop:hxt");
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
	    atomicAdd(&out[ix+iy*nx],alpha*
		      (20.f*in[ix+iy*nx]
		       -8.f*(in[ix1p+iy*nx]+in[ix1+iy*nx]+in[ix+iy1p*nx]+in[ix+iy1*nx])
		       +2.f*(in[ix1p+iy1p*nx]+in[ix1+iy1p*nx]+in[ix1p+iy1*nx]+in[ix1+iy1*nx])
		       +(in[ix+iy2p*nx]+in[ix2p+iy*nx]+in[ix2+iy*nx]+in[ix+iy2*nx])));
	}
    }
}
__global__ static void zzt_do(float *restrict out, const float *in, int ix, float val){
    atomicAdd(&out[ix], in[ix]*val);
}
void gpu_TomoL(curcell **xout, const void *A, const curcell *xin, const float alpha){
    TIC;tic;
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int nwfs=parms->nwfsr;
    curcell *opdx;
    if(!*xout){
	opdx=curcellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    const int nxo=recon->xmap[ips]->nx;
	    const int nyo=recon->xmap[ips]->ny;
	    opdx->p[ips]=curnew(nxo, nyo);
	}
	*xout=opdx;
    }else{
	opdx=*xout;
    }
    /*no need to initialize grad to zero. override the value.*/
    curcell *grad=curecon->grad;
  
    for(int ips=0; ips<recon->npsr; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	laplacian_do<<<DIM2(nxo, nyo, 16, 16), 0, curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, nxo, nyo, curecon->l2c[ips]);
	zzt_do<<<1,1,0,curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, curecon->zzi[ips], curecon->zzv[ips]);
    }
    const int nxp=recon->pmap->nx;
    const int nyp=recon->pmap->ny;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curecon->opdwfs;
    int ptt=(!parms->tomo.split || parms->dbg.splitlrt); 
    curcell *ttf=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
#if !USE_GP 
	const float dsa=cuwfs[iwfs].powfs->dsa;
#endif
	DO_HX;
	DO_GP;
	DO_NEA;
	DO_GPT;
    }
    SYNC_WFS;
    curcellfree(ttf);
    
    for(int ips=0; ips<recon->npsr; ips++){
	DO_HXT;
    }
    SYNC_PS;
}
void gpu_FitR(curcell **xout, const void *A, const curcell *xin, const float alpha){
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int nfit=parms->fit.nfit;
    const int npsr=recon->npsr;
    curcell *opdfit=curcellnew(nfit, 1);
    for(int ifit=0; ifit<nfit; ifit++){
	double hs=parms->fit.ht[ifit];
	float thetax=parms->fit.thetax[ifit];
	float thetay=parms->fit.thetay[ifit];
	opdfit->p[ifit]=curnew(recon->ploc->nloc,1);
	for(int ips=0; ips<npsr; ips++){
	    const double ht = recon->ht->p[ips];
	    float scale=1.f-ht/hs;
	    float dispx=(ht*thetax-recon->xmap[ips]->ox)/recon->xmap[ips]->dx;
	    float dispy=(ht*thetay-recon->xmap[ips]->oy)/recon->xmap[ips]->dx;
	    //prop_linear<<<>>>(opdfit->p[ifit], recn->xloc[ips]->nx, recon->xloc[ips]->ny,
			      
	}
    }
}
void gpu_tomofit(SIM_T *simu){
    cudaStream_t stream;
    STREAM_NEW(stream);
    cublasHandle_t handle;
    DO(cublasCreate(&handle));
    DO(cublasSetStream(handle, stream));
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    //if(!parms->tomo.split){
    //TO_IMPLEMENT;
    //}
    if(parms->tomo.pos!=2){
	TO_IMPLEMENT;
    }
    if(curecon->PDF){
	TO_IMPLEMENT;
    }
    {
	const int nwfs=parms->nwfsr;
	int nxp=recon->pmap->nx;
	int nyp=recon->pmap->ny;
	//Create temporary memory
	curecon->opdwfs=curcellnew(nwfs, 1);
	curecon->grad=curcellnew(nwfs, 1);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    curecon->opdwfs->p[iwfs]=curnew(nxp, nyp);
	    const int nsa=cuwfs[iwfs].powfs->nsa;
	    curecon->grad->p[iwfs]=curnew(nsa*2,1);
	}

    }
    //first send gradients to GPU. can be skipped if keep grad in gpu. fast though.
    TIC;tic;
    gpu_dcell2cu(&curecon->gradin, simu->gradlastol);
    toc("grad");
    curcell *rhs=NULL;
    gpu_TomoR(&rhs, simu, curecon->gradin, 1);
    toc("TomoR");
    curcell *opdx=NULL;
    gpu_pcg(&opdx, gpu_TomoL, simu, NULL, NULL, rhs, simu->recon->warm_restart, parms->tomo.maxit, stream, handle);
    toc("pcg");
    curcellfree(rhs);
    cublasDestroy(handle);
    STREAM_DONE(stream);
    {
	curcellfree(curecon->opdwfs);
	curcellfree(curecon->grad);
    }
}
