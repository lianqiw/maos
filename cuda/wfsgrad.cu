extern "C"
{
#include <cuda.h>
#include "gpu.h"
#include "utils.h"
#include "accphi.h"
}
#include "curand_kernel.h"
#include "cusparse.h"
typedef curandState_t curandStat;
static int      cunpowfs;
static culoc_t *cuwfsloc=NULL;
static float  **cuwfsamp=NULL;
static cusp_t **cuGS0=NULL;
static float  **cugradacc=NULL;
static curandStat **custat=NULL;
static int     *custatb=NULL;//allocated block
static int     *custatt=NULL;//allocated thread
static float  **cunea=NULL;//the noise equivalent angles for each grad.
static cusparseHandle_t *cusphandle=NULL;
static cusparseMatDescr_t cuspdesc;
#define RAND_BLOCK 16
#define RAND_THREAD 32
typedef struct{
    float (*orig)[2];
    float dx;
    int nsa;
    int nx;
    float (***imcc)[3];
}cupts_t;
cupts_t *cupts;
/*
  Save POWFS loc to GPU.
*/
static void gpu_wfsloc2gpu(int npowfs, int ipowfs, loc_t *loc){
    DO(cudaThreadSynchronize());
    if(!cuwfsloc){
	cuwfsloc=(culoc_t*)calloc(npowfs, sizeof(culoc_t));
	cunpowfs=npowfs;
    }else if(cunpowfs!=npowfs){
	error("npowfs mismatch\n");
    }
    if(cuwfsloc[ipowfs].loc){
	error("Already initialized\n");
    }
    
    cuwfsloc[ipowfs].nloc=loc->nloc;
    cuwfsloc[ipowfs].dx=loc->dx;
    gpu_loc2dev(&cuwfsloc[ipowfs].loc, loc);
}

/**
   save GS0 to gpu. Multiple the GS0 directly with opd is extremely slow (not
   progressing) using cuspmul or cusparseScsrmv. 
   After fixing -DLONG
   cusparseScsrmv takes 0.008 seconds for a 60x60 WFS.
   cusptmul takes 0.038 seconds.
*/
static void gpu_GS02gpu(int npowfs, int ipowfs, spcell *GS0){
    if(!cuGS0){
	cuGS0=(cusp_t**)calloc(npowfs, sizeof(cusp_t*));
    }
    if(cuGS0[ipowfs]) error("GS0 is already copied for powfs %d\n", ipowfs);
    if(!GS0 || !GS0->p) return;
    cudaMallocHost(&cuGS0[ipowfs], (GS0->nx*GS0->ny)*sizeof(cusp_t));
    for(int i=0; i<GS0->nx*GS0->ny; i++){
	dsp *t=sptrans(GS0->p[i]);
	sp2gpu(&cuGS0[ipowfs][i], t);
	spfree(t);
    }
}
/**
   Initialize other arrays
*/
void gpu_wfsgrad_init(const PARMS_T *parms, POWFS_T *powfs){
    if(cusphandle) error("Already initialized");
    cudaCallocHost(cusphandle, parms->nwfs * sizeof(cusparseHandle_t));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	DO(cusparseCreate(&cusphandle[iwfs]));
    }
    DO(cusparseCreateMatDescr(&cuspdesc));
    cusparseSetMatType(cuspdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(cuspdesc, CUSPARSE_INDEX_BASE_ZERO);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	gpu_wfsloc2gpu(parms->npowfs, ipowfs, powfs[ipowfs].loc);
	gpu_GS02gpu(parms->npowfs, ipowfs, powfs[ipowfs].GS0);
    }
    if(!cupts){
	cupts=(cupts_t*)calloc(parms->npowfs, sizeof(cupts_t));
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	pts_t *pts=powfs[ipowfs].pts;
	cupts[ipowfs].nx=pts->nx;
	cupts[ipowfs].nsa=pts->nsa;
	cupts[ipowfs].dx=pts->dx;
	gpu_loc2dev(&cupts[ipowfs].orig, (loc_t*)pts);
	int nsa=pts->nsa;
	int mwfs=powfs[ipowfs].nimcc>1 ? parms->powfs[ipowfs].nwfs : 1;
	if(parms->powfs[ipowfs].gtype_recon==1 ||parms->powfs[ipowfs].gtype_sim==1){
	    cupts[ipowfs].imcc=(float(***)[3])malloc(mwfs*sizeof(void*));
	    for(int jwfs=0; jwfs<mwfs; jwfs++){
		cudaMallocHost(&cupts[ipowfs].imcc[jwfs], nsa*sizeof(void*));
		for(int isa=0; isa<nsa; isa++){
		    double *src=powfs[ipowfs].saimcc[jwfs]->p[isa]->p;
		    cupts[ipowfs].imcc[jwfs][isa]=NULL;
		    gpu_dbl2dev((float**)&(cupts[ipowfs].imcc[jwfs][isa]),src, 3*3);
		    cudaDeviceSynchronize();
		}
	    }
	}
    }
    //Do not use cudaMallocHost here
    cugradacc=(float**)calloc(parms->nwfs, sizeof(float*));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	cudaCalloc(cugradacc[iwfs], powfs[ipowfs].pts->nsa*2*sizeof(float));
    }
    //Do not use cudaMallocHost here
    if(!cunea){
	cunea=(float**)calloc(parms->nwfs, sizeof(float*));
    }

    //NEA
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	if(!powfs[ipowfs].neasim) continue;
	dmat *nea=powfs[ipowfs].neasim->p[wfsind];
	if(!nea) continue;
	gpu_dbl2dev(&cunea[iwfs], nea->p, nea->nx*nea->ny);
    }
    if(!cuwfsamp){
	cudaMallocHost(&cuwfsamp, parms->nwfs*sizeof(float*));
    }
    //Amp
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	if(powfs[ipowfs].locm || wfsind==0){
	    cuwfsamp[iwfs]=NULL;
	    gpu_dbl2dev(&cuwfsamp[iwfs], powfs[ipowfs].realamp[wfsind], powfs[ipowfs].loc->nloc);
	}else{
	    cuwfsamp[iwfs]=cuwfsamp[parms->powfs[ipowfs].wfs[0]];
	    if(cuwfsamp[iwfs]==0){
		error("Why?");
	    }
	}
    }

}
__global__ static void setup_rand(curandStat *rstat, int seed){
    int id=threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, id, 0, &rstat[id]);
}
/**
   Seed the random number genrator
 */
void gpu_wfsgrad_seeding(const PARMS_T *parms, POWFS_T *powfs, rand_t *rstat){
    custat=(curandStat**)calloc(parms->nwfs, sizeof(curandStat*));
    custatb=(int*)calloc(parms->nwfs, sizeof(int));
    custatt=(int*)calloc(parms->nwfs, sizeof(int));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int seed=lrand(rstat);//don't put this after continue.
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].noisy) continue;
	int nsa=powfs[ipowfs].pts->nsa*2;
	if(nsa<RAND_THREAD){
	    custatt[iwfs]=nsa;
	    custatb[iwfs]=1;
	}else if(nsa<RAND_THREAD*RAND_BLOCK){
	    custatt[iwfs]=RAND_THREAD;
	    custatb[iwfs]=nsa/RAND_THREAD+(nsa%RAND_THREAD)?1:0;
	}else{
	    custatt[iwfs]=RAND_THREAD;
	    custatb[iwfs]=RAND_BLOCK;
	}
	cudaMalloc(&custat[iwfs], (custatt[iwfs]*custatb[iwfs])*sizeof(curandStat));
	setup_rand<<<custatb[iwfs], custatt[iwfs]>>>(custat[iwfs], seed);
    }
    cudaDeviceSynchronize();
}

__global__ static void add_geom_noise(float *restrict g, const float *restrict nea, int ng, curandStat *restrict rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    const int nstep=blockDim.x * gridDim.x;
    curandStat lstat=rstat[id];
    for(int i=id; i<ng; i+=nstep){
	g[i]+=curand_normal(&lstat)*nea[i];
    }
    rstat[id]=lstat;
}

/**
   Compute ztilt
*/
__global__ static void cuztilt(float *restrict g, float *restrict opd, 
			       const int nsa, const float dx, const int nx, float (**imcc)[3],
			       const float (*orig)[2], const float*restrict amp, float alpha){
    __shared__ float a0,a1,a2;
    if(threadIdx.x==0 && threadIdx.y==0){
	a0=0.f;
	a1=0.f;
	a2=0.f;
    }
    const int isa=blockIdx.x;
    float b0=0.f;
    float b1=0.f;
    float b2=0.f;
    const int skip=isa*nx*nx;
    const float ox=orig[isa][0];
    const float oy=orig[isa][1];
    for(int iy=threadIdx.y; iy<nx; iy+=blockDim.y){
	const int skip2=skip+iy*nx;
	const float y=iy*dx+oy;
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    const int ind=skip2+ix;
	    const float tmp=amp[ind]*opd[ind];
	    b0+=tmp;
	    b1+=tmp*(dx*ix+ox);
	    b2+=tmp*y;
	}
    }
    atomicAdd(&a0,b0);
    atomicAdd(&a1,b1);
    atomicAdd(&a2,b2);
    __syncthreads();//Wait until all threads in this block is done.
    if(threadIdx.x==0 && threadIdx.y==0){
	float (*restrict A)[3]=imcc[isa];
	g[isa]    +=alpha*(A[0][1]*a0+A[1][1]*a1+A[2][1]*a2);
	g[isa+nsa]+=alpha*(A[0][2]*a0+A[1][2]*a1+A[2][2]*a2);
    }
}

/**
   Ray tracing and gradient computation for WFS. \todo Expand to do gradients in GPU without transfering
   data back to CPU.
*/
void gpu_wfsgrad(thread_t *info){
    TIC;tic;
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
    const int iwfs=info->start;
    assert(info->end==info->start+1);//only 1 WFS.
    assert(iwfs<parms->nwfs);
    const POWFS_T *powfs=simu->powfs;
    //output
    const int CL=parms->sim.closeloop;
    const int isim=simu->isim;
    //The following are truly constants for this powfs
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int pixpsa=powfs[ipowfs].pts->nx*powfs[ipowfs].pts->nx;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    const double hs=parms->powfs[ipowfs].hs;
    const int npix=pixpsa*nsa;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom[iwfs];
    const int save_grad=parms->save.grad[iwfs];
    const int save_opd =parms->save.wfsopd[iwfs];
    const int save_ints=parms->save.ints[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    const int nthread=powfs[ipowfs].nthread;
    //The following depends on isim
    const int dtrat_output=((isim+1)%dtrat==0);
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_geom=!do_phy || save_gradgeom;
    const double *realamp=powfs[ipowfs].realamp[wfsind];
    const float thetax=parms->wfs[iwfs].thetax;
    const float thetay=parms->wfs[iwfs].thetay;
    const float mispx=powfs[ipowfs].misreg[wfsind][0];
    const float mispy=powfs[ipowfs].misreg[wfsind][1];
    const float dtisim=parms->sim.dt*isim;
    const float (*loc)[2]=cuwfsloc[ipowfs].loc;
    const int nloc=cuwfsloc[ipowfs].nloc;
    const int jwfs=powfs[ipowfs].nlocm>1?wfsind:0;
    //Out to host for now. \todo : keep grad in device when do reconstruction on device.
    dmat *gradout=simu->gradcl->p[iwfs];
    float *phiout;
    float *gradacc=cugradacc[iwfs];
    cudaCalloc(phiout, nloc*sizeof(float));
    cudaStream_t stream;
    STREAM_NEW(stream);

    if(simu->surfwfs && simu->surfwfs->p[iwfs] || powfs[ipowfs].ncpa){
    //copy to GPU and copy to phiout.
	TO_IMPLEMENT;
    }
    if(parms->sim.idealwfs){
	TO_IMPLEMENT;
    }else{
	gpu_atm2loc(phiout, loc, nloc, hs, thetax, thetay, mispx, mispy, dtisim, 1, stream);
	if(parms->sim.wfsalias){
	    TO_IMPLEMENT;
	}
    }
    if(CL){
	gpu_dm2loc(phiout, loc, nloc, hs, thetax, thetay, mispx, mispy, -1, stream);
    }
    CUDA_SYNC_STREAM;
    //toc2("accphi");tic;
    if(imoao>-1){
	TO_IMPLEMENT;
    }
    if(simu->telws){
	TO_IMPLEMENT;
    }
   
    if(powfs[ipowfs].focus){
	TO_IMPLEMENT;
    }
    if(parms->powfs[ipowfs].llt && simu->focusint && simu->focusint->p[iwfs]){
	TO_IMPLEMENT;
    }
    if(parms->wfs[iwfs].psfmean && isim>=parms->evl.psfisim){
	TO_IMPLEMENT;
    }

    if(do_geom){
	if(parms->powfs[ipowfs].gtype_sim==1){
	    dim3 nth(16,16);
	    cuztilt<<<nsa, nth, 0, stream>>>(gradacc, phiout, cupts[ipowfs].nsa, cupts[ipowfs].dx, 
					   cupts[ipowfs].nx, cupts[ipowfs].imcc[jwfs],
					   cupts[ipowfs].orig, cuwfsamp[iwfs], 1.f/(float)dtrat);
	}else{
	    const int iGS0=powfs[ipowfs].GS0->nx==1?0:wfsind;
#if MYSPARSE==1
	    cusptmul(gradacc, &cuGS0[ipowfs][iGS0], phiout, 1.f/(float)dtrat, stream);
#else
	    cusparseSetKernelStream(cusphandle[iwfs], stream);
	    int status=cusparseScsrmv(cusphandle[iwfs], 
			      CUSPARSE_OPERATION_NON_TRANSPOSE, 
			      cuGS0[ipowfs][iGS0].ny,
			      cuGS0[ipowfs][iGS0].nx,
			      1.f,
			      cuspdesc,
			      cuGS0[ipowfs][iGS0].x,
			      cuGS0[ipowfs][iGS0].p,
			      cuGS0[ipowfs][iGS0].i,
			      phiout, 
			      1.f/(float)dtrat,
			      gradacc);
	    if(status!=0){
		error("cusparseScsrmv failed with status %d\n", status);
	    }
#endif
	}
    }
    CUDA_SYNC_STREAM;
    //toc2("grad");tic;
    if(dtrat_output){
	if(do_phy){
	    TO_IMPLEMENT;
	}else{
	    if(noisy){
		add_geom_noise<<<custatb[iwfs], custatt[iwfs], 0, stream>>>(gradacc, cunea[iwfs], nsa*2, custat[iwfs]);
		CUDA_SYNC_STREAM;
		//toc2("noise");
	    }
	    gpu_dev2dbl(&gradout->p, gradacc, nsa*2);
	    //toc2("dev2dbl");
	    cudaMemset(gradacc, 0, nsa*2*sizeof(float));
	    //toc2("zero");
	}
	if(save_grad || save_gradgeom){
	    TO_IMPLEMENT;
	}
    }
    //toc2("done");
    STREAM_DONE(stream);
    cudaFree(phiout);
    //toc2("final");
}
