extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "curand_kernel.h"
#include "cusparse.h"
#include "cufft.h"
#include "wfs.h"

#define ctoc(A) //CUDA_SYNC_STREAM; toc(A)

cusparseMatDescr_t cuspdesc;
cuwloc_t *cupowfs=NULL;
cuwfs_t *cuwfs=NULL;

/*
  Notice that both blocks and threads are partitioning isa
 */
__global__ static void add_geom_noise_do(float *restrict g, const float *restrict nea, 
				      int ng, curandStat *restrict rstat){
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
__global__ void cuztilt(float *restrict g, float *restrict opd, 
			const int nsa, const float dx, const int nx, float (**imcc)[3],
			const float (*orig)[2], const float*restrict amp, float alpha){
    __shared__ float a0,a1,a2;
    __syncthreads();
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
   Apply matched filter. \todo this implementation relies on shared variable. It is probably causing competition.
*/
__global__ static void mtche_do(float *restrict grad, float (*restrict *restrict mtches)[2], const float *restrict ints, int pixpsa, int nsa){
    __shared__ float g[2];//shared by threads in the same block (with the same isa).
    int isa=blockIdx.x;
    ints+=isa*pixpsa;
    const float (*const restrict mtche)[2]=mtches[isa];
    __syncthreads();
    if(threadIdx.x==0){
	g[0]=g[1]=0.f;
    }
    float gp[2]={0.f,0.f};
    for (int ipix=threadIdx.x; ipix<pixpsa; ipix+=blockDim.x){
	gp[0]+=mtche[ipix][0]*ints[ipix];
	gp[1]+=mtche[ipix][1]*ints[ipix];
    }
    if(fabsf(gp[0])>1e-5 || fabsf(gp[1])>1e-5){
	printf("gp=%g %g. g=%g %g\n", gp[0], gp[1], g[0], g[1]);
    }
    atomicAdd(&g[0], gp[0]);
    atomicAdd(&g[1], gp[1]);
    __syncthreads();
    if(fabsf(g[0])>1e-5 || fabsf(g[1])>1e-5){
	if(threadIdx.x==0){
	    printf("isa=%d id=%d g=%g %g\n", isa, threadIdx.x, g[0], g[1]);
	}
	printf("id=%d gp=%g %g.\n", threadIdx.x, gp[0], gp[1]);
    }
    __syncthreads();
    if(threadIdx.x==0){
	grad[isa]=g[0];
	grad[isa+nsa]=g[1];
    }
}

/**
   Poisson random generator.
*/
__device__ static float curandp(curandStat *rstat, float xm){
    float g, t, xmu;
    int x=0, xu;
    if(xm>200){
	x=(int)round(xm+curand_normal(rstat)*sqrt(xm));
    }else{
	while(xm>0){
	    xmu = xm > 12.f ? 12.f : xm;
	    xm -= xmu;
	    g   = __expf(-xmu);
	    xu  = -1;
	    t   = 1.f;
	    while(t>g){
		xu++;
		t *= curand_uniform(rstat);
	    }
	    x += xu;
	}
    }
    return x;
}
/**
   Add noise to pix images.
*/
__global__ static void addnoise_do(float *restrict ints0, int nsa, int pixpsa, float bkgrnd, float pcalib, 
				   float *const restrict *restrict bkgrnd2s, float *const restrict *restrict bkgrnd2cs,
				   float rne, curandStat *rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    const int nstep=blockDim.x * gridDim.x;
    curandStat lstat=rstat[id];
    for(int isa=id; isa<nsa; isa+=nstep){
	float *restrict ints=ints0+isa*pixpsa;
	const float *restrict bkgrnd2=bkgrnd2s?bkgrnd2s[isa]:NULL;
	const float *restrict bkgrnd2c=bkgrnd2cs?bkgrnd2cs[isa]:NULL;
	for(int ipix=0; ipix<pixpsa; ipix++){
	    if(bkgrnd2){
		ints[ipix]=curandp(&lstat, ints[ipix]+bkgrnd+bkgrnd2[ipix])+rne*curand_normal(&lstat)-bkgrnd*pcalib;
	    }else{
		ints[ipix]=curandp(&lstat, ints[ipix]+bkgrnd)+rne*curand_normal(&lstat)-bkgrnd*pcalib;
	    }
	    if(bkgrnd2c){
		ints[ipix]-=bkgrnd2c[ipix];
	    }
	}
    }
    rstat[id]=lstat;
}
/**
   Collect the noise we actually applied.
 */
__global__ static void collect_noise_do(float *restrict neareal, const float *restrict gnf, const float *restrict gny, int nsa){
    for(int isa=threadIdx.x+blockIdx.x*blockDim.x; isa<nsa; isa+=blockDim.x*gridDim.x){
        float dx=gny[isa]-gnf[isa];
	float dy=gny[isa+nsa]-gnf[isa+nsa];
	float *restrict nea=neareal+isa*4;
	nea[0]+=dx*dx;
	nea[1]+=dx*dy;
	nea[2]+=dx*dy;
	nea[3]+=dy*dy;
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
    const RECON_T *recon=simu->recon;
    //output
    const int CL=parms->sim.closeloop;
    const int isim=simu->isim;
    //The following are truly constants for this powfs
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    const float hs=parms->powfs[ipowfs].hs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom[iwfs];
    const int save_grad=parms->save.grad[iwfs];
    const int save_ints=parms->save.ints[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    //The following depends on isim
    const int dtrat_output=((isim+1)%dtrat==0);
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_geom=!do_phy || save_gradgeom;
    const float thetax=parms->wfs[iwfs].thetax;
    const float thetay=parms->wfs[iwfs].thetay;
    const float mispx=powfs[ipowfs].misreg[wfsind][0];
    const float mispy=powfs[ipowfs].misreg[wfsind][1];
    const float dtisim=parms->sim.dt*isim;
    float (*loc)[2]=cupowfs[ipowfs].loc;
    const int nloc=cupowfs[ipowfs].nloc;
    //Out to host for now. \todo : keep grad in device when do reconstruction on device.
    cudaStream_t stream=cuwfs[iwfs].stream;
    dmat *gradout=simu->gradcl->p[iwfs];
    curmat *phiout=curnew(nloc, 1);
    curzero(phiout, stream);
    float *gradacc=cuwfs[iwfs].gradacc;
    if(cuwfs[iwfs].opdadd){ //copy to phiout.
	curcp(&phiout, cuwfs[iwfs].opdadd, stream);
    }
    if(parms->sim.idealwfs){
	TO_IMPLEMENT;
    }else{
	gpu_atm2loc(phiout->p, loc, nloc, hs, thetax, thetay, mispx, mispy, dtisim, 1, stream);
	if(parms->sim.wfsalias){
	    TO_IMPLEMENT;
	}
    }
    if(CL){
	gpu_dm2loc(phiout->p, loc, nloc, hs, thetax, thetay, mispx, mispy, -1, stream);
    }
    //CUDA_SYNC_STREAM;

    if(imoao>-1){
	TO_IMPLEMENT;
    }
    if(simu->telws){
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(phiout, loc, tt*cosf(angle), tt*sinf(angle), stream);
    }
   
    if(powfs[ipowfs].focus){
	TO_IMPLEMENT;
    }
    if(parms->powfs[ipowfs].llt && simu->focusint && simu->focusint->p[iwfs]){
	TO_IMPLEMENT;
    }
    if(do_geom){
	if(parms->powfs[ipowfs].gtype_sim==1){
	    cuztilt<<<nsa, dim3(16,16), 0, stream>>>(gradacc, phiout->p, cupowfs[ipowfs].nsa, 
					     cupowfs[ipowfs].dx, 
					     cupowfs[ipowfs].nxsa, cuwfs[iwfs].imcc,
					     cupowfs[ipowfs].pts, cuwfs[iwfs].amp, 
					     1.f/(float)dtrat);
	}else{
	    cusp *GS0=cuwfs[iwfs].GS0t;
	    cusptmul(gradacc, GS0, phiout->p, 1.f/(float)dtrat, cuwfs[iwfs].sphandle);
	}
    }   
    //CUDA_SYNC_STREAM;
    if(do_phy || parms->powfs[ipowfs].psfout || (parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart)){//physical optics
	wfsints(simu, phiout->p, iwfs, isim, stream);
    }//do phy
    //CUDA_SYNC_STREAM;
    ctoc("grad");
    if(dtrat_output){
	if(do_phy){
	    //signal level was already multiplied in ints.
	    float *restrict const ints=cuwfs[iwfs].ints;
	    float *gradnf, *gradny=NULL;
	    const int pixpsa=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
	    cudaMalloc(&gradnf, 2*nsa*sizeof(float));
	    switch(parms->powfs[ipowfs].phytypesim){
	    case 1:
		//use 32 instead of pixpsa here. using pixpsa causes random error in g. is this due to lack of ECC?
		mtche_do<<<nsa, 32,0,stream>>>(gradnf, cuwfs[iwfs].mtche, ints, pixpsa, nsa);
		break;
	    default:
		TO_IMPLEMENT;
	    }
	    //CUDA_SYNC_STREAM;
	    ctoc("mtche");tic;
	    if(noisy){
		float rne=parms->powfs[ipowfs].rne;
		float bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
		addnoise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (ints, nsa, pixpsa, bkgrnd, parms->powfs[ipowfs].bkgrndc,
		     cuwfs[iwfs].bkgrnd2, cuwfs[iwfs].bkgrnd2c, 
		     rne, cuwfs[iwfs].custat);
		ctoc("noise");tic;
		DO(cudaMalloc(&gradny, 2*nsa*sizeof(float)));
		switch(parms->powfs[ipowfs].phytypesim){
		case 1:
		    mtche_do<<<nsa, 16, 0, stream>>>(gradny, cuwfs[iwfs].mtche, ints, pixpsa, nsa);
		    break;
		default:
		    TO_IMPLEMENT;
		}
		collect_noise_do<<<MAX(nsa/256,1), MIN(256, nsa), 0, stream>>>(cuwfs[iwfs].neareal, gradnf, gradny, nsa);
		//CUDA_SYNC_STREAM;
	    }
	    //send grad to CPU.
	    gpu_dev2dbl(&gradout->p, gradny?gradny:gradnf, nsa*2, stream);
	    
	    if(save_ints){
		dmat *intstmp=dnew(pixpsa, nsa);
		gpu_dev2dbl(&intstmp->p, ints, pixpsa*nsa, stream);
		CUDA_SYNC_STREAM;
		cellarr_dmat(simu->save->intsnf[iwfs], intstmp);
		dfree(intstmp);
	    }
	    if(save_grad && noisy){
		dmat *gradnftmp=dnew(nsa*2,1);
		gpu_dev2dbl(&(gradnftmp->p), gradnf, nsa*2, stream);
		CUDA_SYNC_STREAM;
		cellarr_dmat(simu->save->gradnf[iwfs], gradnftmp);
		dfree(gradnftmp);
	    }
	    //CUDA_SYNC_STREAM;
	    if(gradny) cudaFree(gradny);
	    cudaFree(gradnf);
	    cudaMemsetAsync(ints, 0, nsa*pixpsa*sizeof(float), stream);
	    ctoc("mtche");
	    if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		if(!recon->PTT){
		    error("powfs %d has llt, but recon->PTT is NULL",ipowfs);
		}
		dmat *PTT=NULL;
		if(parms->sim.glao){
		    PTT=recon->PTT->p[ipowfs+ipowfs*parms->npowfs];
		}else{
		    PTT=recon->PTT->p[iwfs+iwfs*parms->nwfs];
		}
		if(!PTT){
		    error("powfs %d has llt, but TT removal is empty\n", ipowfs);
		}
		/* Compute LGS Uplink error. */
		dzero(simu->upterr->p[iwfs]);
		dmm(&simu->upterr->p[iwfs], PTT, gradout, "nn", 1);
		/* copy upterr to output. */
		PDMAT(simu->upterrs->p[iwfs], pupterrs);
		pupterrs[isim][0]=simu->upterr->p[iwfs]->p[0];
		pupterrs[isim][1]=simu->upterr->p[iwfs]->p[1];
	    }
	}else{
	    if(noisy){
		add_geom_noise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (gradacc, cuwfs[iwfs].neasim, nsa*2,cuwfs[iwfs].custat);
		CUDA_SYNC_STREAM;
		ctoc("noise");
	    }
	    gpu_dev2dbl(&gradout->p, gradacc, nsa*2, stream);
	    ctoc("dev2dbl");
	    cudaMemsetAsync(gradacc, 0, nsa*2*sizeof(float), stream);
	    ctoc("zero");
	}
	CUDA_SYNC_STREAM;
	if(powfs[ipowfs].ncpa_grad){
	    warning("Applying ncpa_grad to gradout\n");
	    dadd(&gradout, 1., powfs[ipowfs].ncpa_grad->p[wfsind], -1.);
	}
	if(save_grad){
	    cellarr_dmat(simu->save->gradcl[iwfs], gradout);
	}
	if(save_gradgeom){
	    dmat *gradtmp=dnew(nsa*2,1);
	    gpu_dev2dbl(&gradtmp->p, gradacc, nsa*2, stream);
	    CUDA_SYNC_STREAM;
	    if(dtrat!=1) dscale(gradtmp, 1./dtrat);
	    cellarr_dmat(simu->save->gradgeom[iwfs], gradtmp);//noise free.
	    dfree(gradtmp);
	}
    }
    ctoc("done");
    curfree(phiout);
}
