extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "wfs.h"
#include "recon.h"
#define DIM(nsa,nb,ng) MAX(1,MIN(nsa/nb,ng)),MIN(nsa,nb)
#define DIM2(nx,ny,nb,ng) dim3(MAX(1,MIN((nx)/(nb),ng)),MAX(1,MIN((nx)/(nb),ng))),dim3(nb,nb)
curecon_t *curecon;
#define SCALE 1 //Scale both NEA and L2 to balance the dynamic range.
#define USE_GP 0
typedef void (*G_CGFUN)(curcell**, const void*, const curcell*, float);
typedef void (*G_PREFUN)(curcell**, const void*, const curcell*);
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
	info("ips=%d l2c=%g\n", ips, curecon->l2c[ips]);
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
	    float diff=0, diff2=0;
	    curinn(&diff, curecon->neai->p[iwfs], curecon->neai->p[iwfs0], handle);
	    curinn(&diff2, curecon->neai->p[iwfs0], curecon->neai->p[iwfs0],  handle);
	    CUDA_SYNC_STREAM;
	    toc("curdot");
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

__global__ static void add_gtt_do(float *restrict grad, const float *const tt, int nsa){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nsa; i+=step){
	grad[i]+=tt[0];
	grad[i+nsa]+=tt[1];
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
		+2*(map[(iy+1)*nx+ix+2]-map[(iy+1)*nx+ix]))*a2;

	g[isa+nsa]=(map[(iy+2)*nx+ix]+map[(iy+2)*nx+ix+2]-map[iy*nx+ix]-map[iy*nx+ix+2]
		    +2*(map[(iy+2)*nx+ix+1]-map[iy*nx+ix+1]))*a2;
    }
}
/*
  Ray tracing with matched spacing. Reverse, from out to in. out is xloc, in is ploc.
*/
__global__ static void gpu_hxt_do_1(const float *restrict out, int nxout,
				    float *restrict in, int nxin, 
				    float fracx, float fracy,
				    int nx, int ny, float alpha){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    float fracx1=1.-fracx;
    float fracy1=1.-fracy;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    float a=out[ix+iy*nxout]*alpha;
	    atomicAdd(&in[ix+iy*nxin], a*fracx1*fracy1);
	    atomicAdd(&in[ix+1+iy*nxin], a*fracx*fracy1);
	    atomicAdd(&in[ix+(iy+1)*nxin],a*fracx1*fracy);
	    atomicAdd(&in[ix+1+(iy+1)*nxin],a*fracx*fracy);
	}
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
		+(in[ix+      iy*nxin]*fracx1+in[ix+1+    iy*nxin]*fracx)*fracy1
		+(in[ix+  (iy+1)*nxin]*fracx1+in[ix+1+(iy+1)*nxin]*fracx)*fracy;
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
			atomicAdd(&out[ix2+(iy2)*nxout], 
				  +in[ix3+    (iy3)*nxin]*fracx21*fracy21
				  +in[ix3+1+  (iy3)*nxin]*fracx2*fracy21
				  +in[ix3+  (iy3+1)*nxin]*fracx21*fracy2
				  +in[ix3+1+(iy3+1)*nxin]*fracx2*fracy2);
		    }
		}
	    }
	}
    }
}
/** Forward Ray Tracing. */
__global__ static void gpu_h_do(){
    

}
/**
   Apply GX'*W (GX=GP*HX). We use the stencil specified in NFIRAOS RTC, instead
   of using sparse matrix.  Specified for oversampling of 2 in ploc.*/
static void gpu_iprop(const PARMS_T *parms, const RECON_T *recon, curcell *opdx,
		      curcell *grad, float alpha, int ptt,
		      cudaStream_t* psstream, cudaStream_t *wfsstream, cublasHandle_t *wfshandle){
    TIC;tic;
    const int nwfs=grad->nx;
    int nxp=recon->pmap->nx;
    int nyp=recon->pmap->ny;
    float oxp=recon->pmap->ox;
    float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	opdwfs->p[iwfs]=curnew(nxp, nyp);
    }
    curcell *ttf=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const float dsa=cuwfs[iwfs].powfs->dsa;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	if(ptt && curecon->PTT && curecon->PTT->p[iwfs+iwfs*nwfs]){
	    curmm(&ttf->p[iwfs], 1, curecon->PTT->p[iwfs+iwfs*nwfs], 
		  grad->p[iwfs], "nn", 1, wfshandle[iwfs]);
	    curscale(ttf->p[iwfs], -1.f, wfsstream[iwfs]);
	    add_gtt_do<<<DIM(nsa,256,16), 0, wfsstream[iwfs]>>>
		(grad->p[iwfs]->p, ttf->p[iwfs]->p, nsa);
	}
	gpu_nea_do<<<DIM(nsa,256,16),0,wfsstream[iwfs]>>>
	    (grad->p[iwfs]->p, (float(*)[2])curecon->neai->p[iwfs]->p, nsa);
	curzero(opdwfs->p[iwfs], wfsstream[iwfs]);
#if USE_GP
	cusptmul(opdwfs->p[iwfs]->p, cupowfs[ipowfs].GP, grad->p[iwfs]->p, 1, curecon->wfssphandle[iwfs]);
#else
	gpu_gpt_o2_do<<<DIM(nsa,256,16),0,wfsstream[iwfs]>>>
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, nsa);
#endif
    }
    info("Synchronize wfs...");
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	cudaStreamSynchronize(wfsstream[iwfs]);
    }
    toc("wfs");
    info2("done\n");
    curcellfree(ttf);
    for(int ips=0; ips<recon->npsr; ips++){
	float ht=recon->ht->p[ips];
	float dx1=1./recon->xmap[ips]->dx;
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	const float oxo=recon->xmap[ips]->ox;
	const float oyo=recon->xmap[ips]->oy;
	
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    const float hs = parms->powfs[ipowfs].hs;
	    const float scale = 1.f - ht/hs;
	    float dispx=(parms->wfsr[iwfs].thetax*ht-oxo+oxp*scale)*dx1;
	    float dispy=(parms->wfsr[iwfs].thetay*ht-oyo+oyp*scale)*dx1;
	    int offx=(int)floorf(dispx); dispx-=offx;
	    int offy=(int)floorf(dispy); dispy-=offy;
	    const float ratio=dxp*scale*dx1;
	    if(fabs(ratio-1.f)<1.e-4){//matched
		int nx=nxo-offx-1; if(nx>nxp) nx=nxp;
		int ny=nyo-offy-1; if(ny>nyp) ny=nyp;
		gpu_hxt_do_1<<<DIM2(nx,ny,16,8),0,psstream[ips]>>>
		    (opdwfs->p[iwfs]->p, nxp, opdx->p[ips]->p+offx+offy*nxo, nxo, dispx, dispy, nx, ny, alpha);
	    }else if(fabs(ratio-0.5f)<1.e-4){//double sampled.
		int nx=(nxo-offx-1)*2-((dispx<0.5f)?0:1); if(nx>nxp) nx=nxp; 
		int ny=(nyo-offy-1)*2-((dispy<0.5f)?0:1); if(ny>nyp) ny=nyp;
		gpu_hxt_do_2<<<DIM2((nx+1)/2,(ny+1)/2,16,8),0,psstream[ips]>>>
		    (opdwfs->p[iwfs]->p, nxp, opdx->p[ips]->p+offx+offy*nxo, nxo, dispx, dispy, nx, ny, alpha);
	    }else{
		TO_IMPLEMENT;
	    }
	}
    }
    info("Synchronize ps ...");
    for(int ips=0; ips<recon->npsr; ips++){
	cudaStreamSynchronize(psstream[ips]);
    }
    toc("hxt");
    info2("done\n");
    curcellfree(opdwfs);
}
/**
   Apply GX. GX=GP*HX.
 */
static void gpu_prop(const PARMS_T *parms, const RECON_T *recon, const curcell *opdx, 
		     curcell *restrict grad, cudaStream_t *wfsstream){
    const int nwfs=grad->nx;
    const int nxp=recon->pmap->nx;
    const int nyp=recon->pmap->ny;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    curcell *opdwfs=curcellnew(nwfs, 1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const float hs = parms->powfs[ipowfs].hs;
	opdwfs->p[iwfs]=curnew(nxp, nyp);
	cudaMemsetAsync(opdwfs->p[iwfs]->p, 0, nxp*nyp*sizeof(float), wfsstream[iwfs]);
	for(int ips=0; ips<recon->npsr; ips++){
	    const float ht=recon->ht->p[ips];
	    const float scale = 1.f - ht/hs;
	    const float dx1=1./recon->xmap[ips]->dx;
	    const int nxo=recon->xmap[ips]->nx;
	    const int nyo=recon->xmap[ips]->ny;
	    const float oxo=recon->xmap[ips]->ox;
	    const float oyo=recon->xmap[ips]->oy;
	    float dispx=(parms->wfsr[iwfs].thetax*ht-oxo+oxp*scale)*dx1;
	    float dispy=(parms->wfsr[iwfs].thetay*ht-oyo+oyp*scale)*dx1;
	    int offx=(int)floorf(dispx); dispx-=offx;
	    int offy=(int)floorf(dispy); dispy-=offy;
	    const float ratio=dxp*scale*dx1;
	    if(fabs(ratio-1.f)<1.e-4){//matched
		int nx=nxo-offx-1; if(nx>nxp) nx=nxp;
		int ny=nyo-offy-1; if(ny>nyp) ny=nyp;
		gpu_hx_do_1<<<DIM2(nx,ny,16,16),0,wfsstream[iwfs]>>>
		    (opdwfs->p[iwfs]->p, nxp, opdx->p[ips]->p+offx+offy*nxo, nxo, dispx, dispy, nx, ny);
	    }else if(fabs(ratio-0.5f)<1.e-4){//double sampled.
		int nx=(nxo-offx-1)*2-((dispx<0.5f)?0:1); if(nx>nxp) nx=nxp; 
		int ny=(nyo-offy-1)*2-((dispy<0.5f)?0:1); if(ny>nyp) ny=nyp;
		gpu_hx_do_2<<<DIM2((nx+1)/2,(ny+1)/2,16,16),0,wfsstream[iwfs]>>>
		    (opdwfs->p[iwfs]->p, nxp, opdx->p[ips]->p+offx+offy*nxo, nxo, dispx, dispy, nx, ny);
	    }else{
		TO_IMPLEMENT;
	    }
	}//for ips
	const float dsa=cuwfs[iwfs].powfs->dsa;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	//no need to initialize grad to zero. override the value.
#if USE_GP
	cuspmul(grad->p[iwfs]->p, cupowfs[ipowfs].GP, opdwfs->p[iwfs]->p, 1, curecon->wfssphandle[iwfs]);
#else
	gpu_gp_o2_do<<<DIM(nsa,256,16),0,wfsstream[iwfs]>>>
	    (opdwfs->p[iwfs]->p, nxp, grad->p[iwfs]->p, cuwfs[iwfs].powfs->saptr, dsa, nsa);
#endif
    }
    info("Synchronize wfs...");
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	cudaStreamSynchronize(wfsstream[iwfs]);
    }
    info2("done\n");
    curcellfree(opdwfs);
}
#define PRINT_RES 0
/** The PCG algorithm. Copy of lib/pcg.c, but replacing dcell with curcell. */
void gpu_pcg(curcell **px, 
	     G_CGFUN Amul, const void *A, 
	     G_PREFUN Mmul, const void *M, 
	     const curcell *b, int warm, int maxiter,
	     cudaStream_t stream,
	     cublasHandle_t handle){
    curcell *r0=NULL;
    curcell *x0=NULL;//The initial vector.
    curcell *z0=NULL;//Is reference or preconditioned value.
    //computes r0=b-A*x0
    curcellcp(&r0, b, stream);
    if(!*px || !warm){//start from zero guess.
	x0=curcellnew2(b);
	if(!*px) curcellcp(px, x0, stream);//initialize the output;
    }else{
	curcellcp(&x0, *px, stream);
	CUDA_SYNC_STREAM;
	Amul(&r0, A, x0, -1);//r0=r0+(-1)*A*x0
    }
    double r0z1,r0z2,r0zmin;
    curcell *p0=NULL;
    if(Mmul){
	Mmul(&z0,M,r0);
    }else{
	z0=r0;
    }
    curcellcp(&p0, z0, stream);
    r0z1=0;
    r0z1=curcellinn(r0, z0, handle);
    curcell *Ap=NULL;
    double ak,bk;
    r0zmin=r0z1;
#if PRINT_RES == 1
    double r0z0=0;
    r0z0=curcellinn(b, b, handle);//|b|
    double res[maxiter+1];
    if(Mmul){
	double tmp=curcellinn(r0, r0, handle);
	res[0]=sqrt(tmp/r0z0);
    }else{
	res[0]=sqrt(r0z1/r0z0);
    }
    int kres=0;
    info("Step %d, res=%g\n", kres, res[kres]);
#endif
    for(int k=0; k<maxiter; k++){
	if(Ap) curcellzero(Ap, stream);
	Amul(&Ap, A, p0, 1);
	ak=r0z1/curcellinn(p0,Ap,handle);
	curcelladd(&x0, 1, p0, ak, handle);//x0=x0+ak*p0
	curcelladd(&r0, 1, Ap, -ak, handle);//r0=r0-ak*Ap
	CUDA_SYNC_STREAM;
	if(Mmul){
	    Mmul(&z0,M,r0);
	}else{
	    z0=r0;
	}
	r0z2=curcellinn(r0,z0,handle);//r0z2=r0'*z0
	bk=r0z2/r0z1;
	curcelladd(&p0, bk, z0, 1., handle);//p0=bk*p0+z0
	if(r0z2<r0zmin){
	    curcellcp(px, x0, stream);//record new result
	    r0zmin=r0z2;
#if PRINT_RES == 1
	    kres=k;
#endif
	}
	r0z1=r0z2;
#if PRINT_RES == 1
	if(Mmul){
	    res[k+1]=sqrt(curcellinn(r0, r0, handle)/r0z0);
	}else{
	    res[k+1]=sqrt(r0z2/r0z0);
	}
	info("Step %d, res=%g\n", k+1, res[k+1]);
#endif
    }
#if PRINT_RES == 1
    info("Solution found at step %d with residual %g\n", kres+1, res[kres+1]);
#endif
    curcellfree(r0); 
    if(Mmul){
	curcellfree(z0);
    }
    curcellfree(x0);
    curcellfree(Ap);
    curcellfree(p0);
}
void gpu_TomoR(curcell **xout, const void *A, const curcell *gin, const float alpha){
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
    gpu_iprop(parms, recon, opdx, curecon->grad, alpha, 1, 
	      curecon->psstream, curecon->wfsstream, curecon->wfshandle);
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
    curcell *grad=curcellnew(parms->nwfsr, 1);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const int nsa=cuwfs[iwfs].powfs->nsa;
	grad->p[iwfs]=curnew(nsa*2,1);
    }
    toc("new");
    for(int ips=0; ips<recon->npsr; ips++){
	const int nxo=recon->xmap[ips]->nx;
	const int nyo=recon->xmap[ips]->ny;
	laplacian_do<<<DIM2(nxo, nyo, 16, 16), 0, curecon->psstream[ips]>>>
	    (opdx->p[ips]->p, xin->p[ips]->p, nxo, nyo, curecon->l2c[ips]);
	zzt_do<<<1,1,0,curecon->psstream[ips]>>>(opdx->p[ips]->p, xin->p[ips]->p,
	curecon->zzi[ips], curecon->zzv[ips]);
    }
    toc("L2");
    gpu_prop(parms, recon, xin, grad, curecon->wfsstream);
    toc("prop");
    int ttr=(!parms->tomo.split || parms->dbg.splitlrt);
 
    gpu_iprop(parms, recon, opdx, grad, alpha, ttr, 
	      curecon->psstream, curecon->wfsstream, curecon->wfshandle);
    toc("iprop");
    curcellfree(grad);
    toc("TomoL");
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

    //if(!parms->tomo.split){
    //TO_IMPLEMENT;
    //}
    if(parms->tomo.pos!=2){
	TO_IMPLEMENT;
    }
    if(curecon->PDF){
	TO_IMPLEMENT;
    }
    //first send gradients to GPU. can be skipped if keep grad in gpu. fast though.
    TIC;tic;
    gpu_dcell2cu(&curecon->grad, simu->gradlastol);
    toc("grad");
    curcell *rhs=NULL;
    gpu_TomoR(&rhs, simu, curecon->grad, 1);
    toc("TomoR");
    curcell *opdx=NULL;
    gpu_pcg(&opdx, gpu_TomoL, simu, NULL, NULL, rhs, simu->recon->warm_restart, parms->tomo.maxit, stream, handle);
    toc("pcg");
    curcellfree(rhs);
    cublasDestroy(handle);
    STREAM_DONE(stream);
}
