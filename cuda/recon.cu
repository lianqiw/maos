extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "pcg.h"
curecon_t *curecon;
#define SCALE 1e-12 //Scale both NEA and L2 to balance the dynamic range. Does not work yet (strange).
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
    CUDA_SYNC_DEVICE;
    cudaStream_t stream;
    STREAM_NEW(stream);
    cublasHandle_t handle;
    DO(cublasCreate(&handle));
    DO(cublasSetStream(handle, stream));

    curecon=(curecon_t*)calloc(1, sizeof(curecon_t));
    curecon->wfsstream=(cudaStream_t*)calloc(parms->nwfsr, sizeof(cudaStream_t));
    curecon->wfshandle=(cublasHandle_t*)calloc(parms->nwfsr, sizeof(cublasHandle_t));
    curecon->wfssphandle=(cusparseHandle_t*)calloc(parms->nwfsr, sizeof(cusparseHandle_t));
    curecon->wfsevent=(cudaEvent_t*)calloc(parms->nwfsr, sizeof(cudaEvent_t));

    curecon->fitstream=(cudaStream_t*)calloc(parms->fit.nfit, sizeof(cudaStream_t));
    curecon->fitsphandle=(cusparseHandle_t*)calloc(parms->fit.nfit, sizeof(cusparseHandle_t));
    curecon->fithandle=(cublasHandle_t*)calloc(parms->fit.nfit, sizeof(cublasHandle_t));

    curecon->psstream=(cudaStream_t*)calloc(recon->npsr, sizeof(cudaStream_t));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	STREAM_NEW(curecon->wfsstream[iwfs]);
	DO(cublasCreate(&curecon->wfshandle[iwfs]));
	DO(cublasSetStream(curecon->wfshandle[iwfs], curecon->wfsstream[iwfs]));
	DO(cusparseCreate(&curecon->wfssphandle[iwfs]));
	DO(cusparseSetKernelStream(curecon->wfssphandle[iwfs], curecon->wfsstream[iwfs]));
	DO(cudaEventCreateWithFlags(&curecon->wfsevent[iwfs], cudaEventDisableTiming));
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	STREAM_NEW(curecon->fitstream[ifit]);
	DO(cusparseCreate(&curecon->fitsphandle[ifit]));
	DO(cusparseSetKernelStream(curecon->fitsphandle[ifit], curecon->fitstream[ifit]));
	DO(cublasCreate(&curecon->fithandle[ifit]));
	DO(cublasSetStream(curecon->fithandle[ifit], curecon->fitstream[ifit]));
    }
    for(int ips=0; ips<recon->npsr; ips++){
	STREAM_NEW(curecon->psstream[ips]);
    }
    STREAM_NEW(curecon->cgstream);
    DO(cublasCreate(&curecon->cghandle));
    DO(cublasSetStream(curecon->cghandle, curecon->cgstream));
    if(parms->gpu.tomo){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].skip) continue;
	    int nsa=powfs[ipowfs].pts->nsa;
	    cudaMalloc(&cupowfs[ipowfs].saptr, nsa*2*sizeof(int));
	    saloc2ptr_do<<<MAX(MIN(nsa/256,16),1), MIN(256, nsa)>>>
		(cupowfs[ipowfs].saptr, cupowfs[ipowfs].saloc, nsa, 
		 recon->pmap->ox, recon->pmap->oy, recon->pmap->dx);
	    if(recon->GP->p[ipowfs]){
		const int use_mat=parms->tomo.pos==2;
		if(use_mat){
		    dsp *GP=sptrans(recon->GP->p[ipowfs]);
		    spint *pp=GP->p;
		    spint *pi=GP->i;
		    double *px=GP->x;
		    dmat *partx=NULL;
		    dmat *party=NULL;
		    partx=dnew(9, nsa);
		    party=dnew(9, nsa);

		    int nsa=powfs[ipowfs].pts->nsa;
		    double dx1=1./recon->ploc->dx;
		
		    for(int ic=0; ic<GP->n; ic++){
			int isa=(ic<nsa)?ic:(ic-nsa);
			for(int ir=pp[ic]; ir<pp[ic+1]; ir++){
			    int ix=pi[ir];
			    double lx=recon->ploc->locx[ix];
			    double ly=recon->ploc->locy[ix];
			    double sx=powfs[ipowfs].saloc->locx[isa];
			    double sy=powfs[ipowfs].saloc->locy[isa];
			    int zx=(int)round((lx-sx)*dx1);
			    int zy=(int)round((ly-sy)*dx1);
			    /**
			       Squeeze the weights to closed points in 3x3 in the subaperture.
			       Does not work well. simply drop these points doesn't work either.
			    */
			    if(zx<0 || zx>2 || zy<0 || zy>2){
				warning("isa=%d, zxy=%d %d\n", isa, zx, zy);
			    }
			    if(zx<0) zx=0;
			    if(zx>2) zx=2;
			    if(zy<0) zy=0;
			    if(zy>2) zy=2;
			    if(ic<nsa){//x
				partx->p[9*isa+zx+zy*3]+=px[ir];
			    }else{//y
				party->p[9*isa+zx+zy*3]+=px[ir];
			    }
			}
		    }
		    gpu_dmat2cu(&cupowfs[ipowfs].GPpx, partx);
		    gpu_dmat2cu(&cupowfs[ipowfs].GPpy, party);
		    dfree(partx);
		    dfree(party);
		    spfree(GP);
		}else{//use sparse
		    gpu_sp2dev(&cupowfs[ipowfs].GP, recon->GP->p[ipowfs]);
		}
	    }else{
		error("GP is required\n");
	    }
	}
 
	curecon->l2c=(float*)calloc(recon->npsr, sizeof(float));
	for(int ips=0; ips<recon->npsr; ips++){
	    float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
	    curecon->l2c[ips]=tmp*tmp*SCALE;
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
		curecon->zzv[ips]=pow(laplacian_coef(r0,wt,dx),2)*SCALE;
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
	    float (*neai)[3]=(float(*)[3])calloc(3*nsa, sizeof(float));
	    if(nea->n!=2*nsa) error("nea doesn't have 2nsa x 2nsa dimension\n");
	    for(int ic=0; ic<nea->n; ic++){
		for(int ir=pp[ic]; ir<pp[ic+1]; ir++){
		    int ix=pi[ir];
		    int isa=ic<nsa?ic:ic-nsa;
		    if(ix==ic){//diagonal part.
			if(ic==isa){//x
			    neai[isa][0]=(float)px[ir]*SCALE;
			}else{//y
			    neai[isa][1]=(float)px[ir]*SCALE;
			}
		    }else if(ix==ic-nsa || ix==ic+nsa){//cross part. symmetric.
			neai[isa][2]=(float)px[ir]*SCALE;
		    }else{
			error("saneai has invalid format\n");
		    }
		}
	    }
	    curecon->neai->p[iwfs]=curnew(3, nsa);
	    DO(cudaMemcpy(curecon->neai->p[iwfs]->p, neai, 3*nsa*sizeof(float), cudaMemcpyDefault));
	    free(neai);
	    if(iwfs!=iwfs0 && nsa>4){//don't check tt. overflows.
		float diff=curinn(curecon->neai->p[iwfs], curecon->neai->p[iwfs0], stream);
		float diff2=curinn(curecon->neai->p[iwfs0], curecon->neai->p[iwfs0],  stream);
		if((diff-diff2)<1e-4*diff2){
		    curfree(curecon->neai->p[iwfs]);
		    curecon->neai->p[iwfs]=curecon->neai->p[iwfs0];
		}
	    }
	}//for iwfs
	CUDA_SYNC_DEVICE;
	if(recon->PTT && !curecon->PTT){
	    gpu_dcell2cu(&curecon->PTT, recon->PTT);
	}
	if(recon->PDF && !curecon->PDF){
	    gpu_dcell2cu(&curecon->PDF, recon->PDF);
	}
    }
    if(parms->gpu.fit){
	//For fitting
	/*
	  gpu_dmat2cu(&curecon->W1, recon->W1);
	  if(recon->W0){
	  spint *pp=recon->W0->p;
	  spint *pi=recon->W0->i;
	  double *px=recon->W0->x;
	  dsp *W0new=spnew(recon->W0->m, recon->W0->n, recon->W0->nzmax);
	  spint *pp2=W0new->p;
	  spint *pi2=W0new->i;
	  double *px2=W0new->x;
	  int *full;
	  cudaMallocHost(&full, recon->W0->n*sizeof(int));
	  double W1max=dmax(recon->W1);
	  double thres=W1max*(1.f-1e-6);
	  curecon->W0v=(float)(W1max*4./9.);//max of W0 is 4/9 of max of W1.
	  info("W0v=%g\n", curecon->W0v);
	  int count=0;
	  int count2=0;
	  for(int ic=0; ic<recon->W0->n; ic++){
	  pp2[ic]=count;
	  if(recon->W1->p[ic]>thres){//full
	  full[count2]=ic;
	  count2++;
	  }else{
	  int nv=pp[ic+1]-pp[ic];
	  memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
	  memcpy(px2+count, px+pp[ic], sizeof(double)*nv);
	  count+=nv;
	  }
	  }
	  pp2[recon->W0->n]=count;
	  W0new->nzmax=count;
	  dsp *W0new2=sptrans(W0new);
	  gpu_sp2dev(&curecon->W0p, W0new2);
	  gpu_int2dev(&curecon->W0f, full, count2);
	  curecon->nW0f=count2;
	  spfree(W0new);
	  spfree(W0new2);
	  cudaFreeHost(full);
	  }else{
	  error("W0, W1 is required\n");
	  }*/
	gpu_muv2dev(&curecon->FR, &recon->FR);
	gpu_muv2dev(&curecon->FL, &recon->FL);
    }
    STREAM_DONE(stream);
    cublasDestroy(handle);
}

void cumuv(curcell **out, cumuv_t *A, const curcell *in, float alpha){
    if(!A->Mt) error("A->M Can not be empty\n");
    if(A->U->ny>1 || A->V->ny>1) error("Not handled yet\n");
    if(!*out){
	*out=curcellnew(A->Mt->ny, 1);
	for(int ix=0; ix<A->Mt->ny; ix++){
	    (*out)->p[ix]=curnew(A->Mt->p[ix*A->Mt->nx]->ny, 1);
	}
    }
    cusparseHandle_t sphandle=curecon->fitsphandle[0];
    cublasHandle_t handle=curecon->fithandle[0];
    for(int iy=0; iy<A->Mt->nx; iy++){
	for(int ix=0; ix<A->Mt->ny; ix++){
	    cusptmul((*out)->p[ix]->p, A->Mt->p[iy+ix*A->Mt->nx], in->p[iy]->p, alpha, sphandle);
	}
    }
    curmat *tmp=curnew(A->V->p[0]->ny, 1);
    for(int iy=0; iy<A->V->nx; iy++){
	curmv(&tmp, 1.f, A->V->p[iy], in->p[iy], 't', 1, handle);
    }
    for(int ix=0; ix<A->U->nx; ix++){
	curmv(&(*out)->p[ix], 1.f, A->U->p[ix], tmp, 'n', -alpha, handle);
    }
    curfree(tmp);
    cudaStreamSynchronize(curecon->fitstream[0]);
}
__global__ void nothing(void){

}
void gpu_tomofit(SIM_T *simu){
    TIC;tic;
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
    if(0){
	//Debugging.
	dcell *rhsc=NULL;
	dcell *lc=NULL;
	muv(&rhsc, &recon->RR, simu->gradlastol, 1);
	muv(&lc, &recon->RL, rhsc, -1);
	dcellwrite(rhsc, "CPU_TomoR");
	dcellwrite(lc, "CPU_TomoL");

	gpu_dcell2cu(&curecon->gradin, simu->gradlastol);
	curcell *rhsg=NULL;
	curcell *lg=NULL;
	gpu_TomoR(&rhsg, simu, curecon->gradin, 1);
	curcellwrite(rhsg, "GPU_TomoR");
	gpu_TomoL(&lg, simu, rhsg, -1);
	curcellwrite(lg, "GPU_TomoL");
	CUDA_SYNC_DEVICE;
	exit(0);
    }
 
    //first send gradients to GPU. can be skipped if keep grad in gpu. fast though.
    curcell *rhs=NULL;
    if(parms->gpu.tomo){
	toc("Before gradin");
	gpu_dcell2cu(&curecon->gradin, simu->gradlastol);
	toc("Gradin");
	gpu_TomoR(&rhs, simu, curecon->gradin, 1);
	toc("TomoR");
	gpu_pcg(&curecon->opdr, gpu_TomoL, simu, NULL, NULL, rhs, 
		simu->parms->recon.warm_restart, parms->tomo.maxit, curecon->cgstream);
	toc("TomoL CG");
	curcellfree(rhs); rhs=NULL;
	curcellfree(curecon->opdwfs); curecon->opdwfs=NULL;
	curcellfree(curecon->grad);   curecon->grad=NULL;
	curcellfree(curecon->gradin); curecon->gradin=NULL;
	if(parms->save.opdr){
	    gpu_curcell2d(&simu->opdr, curecon->opdr, curecon->cgstream);
	}
    }else{
	//Use CPU tomography to compare fit with cpu.
	muv_solve(&simu->opdr, &recon->RL, &recon->RR, simu->gradlastol);
	gpu_dcell2cu(&curecon->opdr, simu->opdr);
	for(int i=0; i<curecon->opdr->nx; i++){
	    curecon->opdr->p[i]->nx=recon->xmap[i]->nx;
	    curecon->opdr->p[i]->ny=recon->xmap[i]->ny;
	}
	toc("CPU Tomo");
    }

    if(parms->gpu.fit){
	int nfit=parms->fit.nfit;
	curecon->opdfit=curcellnew(nfit, 1);
	curecon->opdfit2=curcellnew(nfit,1);
	for(int ifit=0; ifit<nfit; ifit++){
	    curecon->opdfit->p[ifit]=curnew(nxp, nyp);
	    curecon->opdfit2->p[ifit]=curnew(nxp, nyp);
	}
	toc("Before FitR");
	cumuv(&rhs, &curecon->FR, curecon->opdr, 1);
	toc("FitR");
	gpu_pcg(&curecon->dmfit, (G_CGFUN)cumuv, &curecon->FL, NULL, NULL, rhs,
		simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream);
	toc("FitL CG");
	gpu_curcell2d(&simu->dmfit_hi, curecon->dmfit, curecon->cgstream);
	curcellfree(curecon->opdfit);
	curcellfree(curecon->opdfit2);
    }else{
	gpu_curcell2d(&simu->opdr, curecon->opdr, curecon->cgstream);
	for(int i=0; i<simu->opdr->nx; i++){
	    simu->opdr->p[i]->nx=simu->opdr->p[i]->nx*simu->opdr->p[i]->ny;
	    simu->opdr->p[i]->ny=1;
	}
	muv_solve(&simu->dmfit_hi, &recon->FL, &recon->FR, simu->opdr);
	toc("CPU Fit");
    }
    //Don't free opdr.
    curcellfree(rhs); rhs=NULL;
    toc("Final");
}
