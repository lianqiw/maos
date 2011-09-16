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

#undef TIMING
#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
#endif
cusparseMatDescr_t cuspdesc;
cuwloc_t *cupowfs=NULL;
cuwfs_t *cuwfs=NULL;

/*
  Notice that both blocks and threads are partitioning isa
 */
__global__ static void add_geom_noise_do(float *restrict g, const float *restrict nea, 
				      int nsa, curandStat *restrict rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    curandStat lstat=rstat[id];
    const int nstep=blockDim.x * gridDim.x;
    for(int i=id; i<nsa; i+=nstep){
	float n1=curand_normal(&lstat);
	float n2=curand_normal(&lstat);
	g[i]+=n1*nea[i];
	g[i+nsa]+=n2*nea[i+nsa]+n1*nea[i+nsa*2];//cross term.
    }
    rstat[id]=lstat;
}

/**
   Compute ztilt
*/
__global__ void cuztilt(float *restrict g, float *restrict opd, 
			const int nsa, const float dx, const int nx, float (**imcc)[3],
			const float (*orig)[2], const float*restrict amp, float alpha){
    __shared__ float a[3];
    if(threadIdx.x<3 && threadIdx.y==0){
	a[threadIdx.x]=0.f;
    }
    __syncthreads();
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
    atomicAdd(&a[0],b0);
    atomicAdd(&a[1],b1);
    atomicAdd(&a[2],b2);
    __syncthreads();//Wait until all threads in this block is done.
    if(threadIdx.x<3 && threadIdx.y==0){
	float (*restrict A)[3]=imcc[isa];
	atomicAdd(&g[isa],     alpha*(A[threadIdx.x][1]*a[threadIdx.x]));
	atomicAdd(&g[isa+nsa], alpha*(A[threadIdx.x][2]*a[threadIdx.x]));
    }
    /*
      if(threadIdx.x==0 && threadIdx.y==0){
      g[isa]    +=alpha*(A[0][1]*a[0]+A[1][1]*a[1]+A[2][1]*a[2]);
      g[isa+nsa]+=alpha*(A[0][2]*a[0]+A[1][2]*a[1]+A[2][2]*a[2]);
      }*/
}
/**
   Apply matched filter. \todo this implementation relies on shared variable. It is probably causing competition.
*/
__global__ static void mtche_do(float *restrict grad, float (*restrict *restrict mtches)[2], 
				const float *restrict ints, int pixpsa, int nsa){
    __shared__ float g[2];//shared by threads in the same block (with the same isa).
    if(threadIdx.x<2){
	g[threadIdx.x]=0.f;
    }
    __syncthreads();
    int isa=blockIdx.x;
    ints+=isa*pixpsa;
    const float (*const restrict mtche)[2]=mtches[isa];
 
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
    if(threadIdx.x<2){
	grad[isa+nsa*threadIdx.x]=g[threadIdx.x];
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
    curmat *gradacc=cuwfs[iwfs].gradacc;
    if(cuwfs[iwfs].opdadd){ //copy to phiout.
	curcp(&phiout, cuwfs[iwfs].opdadd, stream);
    }
    if(parms->sim.idealwfs){
	gpu_dm2loc(phiout->p, loc, nloc, cudmproj, hs, thetax, thetay, mispx, mispy, 1, stream);
    }else{
	gpu_atm2loc(phiout->p, loc, nloc, hs, thetax, thetay, mispx, mispy, dtisim, 1, stream);
	if(parms->sim.wfsalias){
	    gpu_dm2loc(phiout->p, loc, nloc, cudmproj, hs, thetax, thetay, mispx, mispy, -1, stream);
	}
    }
    if(CL){
	gpu_dm2loc(phiout->p, loc, nloc, cudmreal, hs, thetax, thetay, mispx, mispy, -1, stream);
    }
    //CUDA_SYNC_STREAM;
    
    if(imoao>-1){
	TO_IMPLEMENT;
    }
    if(simu->telws){
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(phiout, loc, 0, tt*cosf(angle), tt*sinf(angle), stream);
    }
   
    if(powfs[ipowfs].focus){
	TO_IMPLEMENT;
    }
    if(parms->powfs[ipowfs].llt && simu->focusint && simu->focusint->p[iwfs]){
	TO_IMPLEMENT;
    }
    if(do_geom){
	if(parms->powfs[ipowfs].gtype_sim==1){
	    cuztilt<<<nsa, dim3(16,16), 0, stream>>>(gradacc->p, phiout->p, cupowfs[ipowfs].nsa, 
					     cupowfs[ipowfs].dx, 
					     cupowfs[ipowfs].nxsa, cuwfs[iwfs].imcc,
					     cupowfs[ipowfs].pts, cuwfs[iwfs].amp, 
					     1.f/(float)dtrat);
	}else{
	    cusp *GS0=cuwfs[iwfs].GS0t;
	    cusptmul(gradacc->p, GS0, phiout->p, 1.f/(float)dtrat, cuwfs[iwfs].sphandle);
	}
    }   
    //CUDA_SYNC_STREAM;
    if(do_phy || parms->powfs[ipowfs].psfout 
       || (parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart)){//physical optics
	wfsints(simu, phiout->p, iwfs, isim, stream);
    }//do phy
    //CUDA_SYNC_STREAM;
    toc("grad");
    if(dtrat_output){
	if(do_phy){
	    //signal level was already multiplied in ints.
	    curmat *ints=cuwfs[iwfs].ints;
	    curmat *gradnf=curnew(nsa,2);
	    curmat *gradny=NULL;
	    const int pixpsa=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
	    switch(parms->powfs[ipowfs].phytypesim){
	    case 1:
		/*use 32 instead of pixpsa here. using pixpsa causes random
		  error in g. is this due to lack of ECC?*/
		mtche_do<<<nsa, 32,0,stream>>>(gradnf->p, cuwfs[iwfs].mtche, ints->p, pixpsa, nsa);
		break;
	    default:
		TO_IMPLEMENT;
	    }
	    if(save_ints){
		cellarr_cur(simu->save->intsnf[iwfs], ints, stream);
	    }
	    //CUDA_SYNC_STREAM;
	    toc("mtche");tic;
	    if(noisy){
		float rne=parms->powfs[ipowfs].rne;
		float bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
		addnoise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (ints->p, nsa, pixpsa, bkgrnd, parms->powfs[ipowfs].bkgrndc,
		     cuwfs[iwfs].bkgrnd2, cuwfs[iwfs].bkgrnd2c, 
		     rne, cuwfs[iwfs].custat);
		toc("noise");tic;
		gradny=curnew(nsa,2);
		switch(parms->powfs[ipowfs].phytypesim){
		case 1:
		    mtche_do<<<nsa, 16, 0, stream>>>(gradny->p, cuwfs[iwfs].mtche, ints->p, pixpsa, nsa);
		    break;
		default:
		    TO_IMPLEMENT;
		}
		collect_noise_do<<<DIM(nsa,256), 0, stream>>>
		    (cuwfs[iwfs].neareal, gradnf->p, gradny->p, nsa);
		if(save_ints){
		    cellarr_cur(simu->save->intsny[iwfs], ints, stream);
		}
		if(save_grad){
		    cellarr_cur(simu->save->gradnf[iwfs], gradnf, stream);
		}
	    }
	    //send grad to CPU.
	    gpu_dev2dbl(&gradout->p, gradny?gradny->p:gradnf->p, nsa*2, stream);
	    curzero(ints, stream);
	    CUDA_SYNC_STREAM;//necessary.
	    curfree(gradny);
	    curfree(gradnf);
	    toc("mtche");
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
		if(save_grad){
		    cellarr_cur(simu->save->gradnf[iwfs], gradacc, stream);
		}
		add_geom_noise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (gradacc->p, cuwfs[iwfs].neasim, nsa,cuwfs[iwfs].custat);
		toc("noise");
	    }
	    gpu_cur2d(&gradout, gradacc, stream);
	    //gpu_dev2dbl(&gradout->p, gradacc->p, nsa*2, stream);
	    toc("dev2dbl");
	    curzero(gradacc, stream);
	    //cudaMemsetAsync(gradacc, 0, nsa*2*sizeof(float), stream);
	    toc("zero");
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
	    if(dtrat!=1) curscale(gradacc, 1./dtrat, stream);
	    cellarr_cur(simu->save->gradgeom[iwfs], gradacc, stream);
	    curzero(gradacc, stream);
	}
    }//dtrat_output
    toc("done");
    CUDA_SYNC_STREAM;
    curfree(phiout);
}
