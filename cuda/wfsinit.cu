/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "accphi.h"
#include <curand_kernel.h>
#include <cusparse.h>
#include <cufft.h>
#include "wfs.h"
static void gpu_pts2cuwloc(cuwloc_t *wloc, pts_t *pts, loc_t *loc){
    wloc->nxsa=pts->nx;
    wloc->nsa=pts->nsa;
    wloc->dx=pts->dx;
    wloc->nloc=loc->nloc;
    gpu_loc2dev(&wloc->pts, (loc_t*)pts);
    gpu_loc2dev(&wloc->loc, loc);
}
int *wfsgpu=NULL;/*assign GPU to wfs statically. */
/**
   Initialize other arrays
*/
void gpu_wfsgrad_init(const PARMS_T *parms, const POWFS_T *powfs){
    wfsgpu=(int*)calloc(parms->nwfs, sizeof(int));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	wfsgpu[iwfs]=gpu_next();
    }
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->powfs=(cuwloc_t*)calloc(parms->npowfs, sizeof(cuwloc_t));
	cudata->wfs=(cuwfs_t*)calloc(parms->nwfs, sizeof(cuwfs_t));
	cuwloc_t *cupowfs=cudata->powfs;
	
	/*DO(cudaSetDeviceFlags(cudaDeviceBlockingSync)); */
	DO(cusparseCreateMatDescr(&cudata->wfsspdesc));
	cusparseSetMatType(cudata->wfsspdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(cudata->wfsspdesc, CUSPARSE_INDEX_BASE_ZERO);
	
	/* Setup information that are same for wfs in each powfs*/
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs==0) continue;
	    pts_t *pts=powfs[ipowfs].pts;
	    loc_t *loc=powfs[ipowfs].loc;
	    gpu_pts2cuwloc(&cupowfs[ipowfs], pts, loc);
	    gpu_loc2dev(&cupowfs[ipowfs].saloc, powfs[ipowfs].saloc);
	    cupowfs[ipowfs].dsa=pts->dsa;
	    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		pts=powfs[ipowfs].llt->pts;
		loc=powfs[ipowfs].llt->loc;
		cupowfs[ipowfs].llt=(cuwloc_t*)calloc(1, sizeof(cuwloc_t));
		gpu_pts2cuwloc(cupowfs[ipowfs].llt, pts, loc);
	    }
	    /*cupowfs[ipowfs].skip=parms->powfs[ipowfs].skip; */
	}
    }

    /* setup information that maybe different for wfs in same powfs due to
       misregistration or NCPA.*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwloc_t *cupowfs=cudata->powfs;
	cuwfs_t *cuwfs=cudata->wfs;
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	int iwfs0=parms->powfs[ipowfs].wfs[0];/*first wfs in this group. */
	/*imcc for ztilt. */
	STREAM_NEW(cuwfs[iwfs].stream);
	DO(cusparseCreate(&cuwfs[iwfs].sphandle));
	DO(cusparseSetKernelStream(cuwfs[iwfs].sphandle, cuwfs[iwfs].stream));
	DO(cublasCreate(&cuwfs[iwfs].handle));
	DO(cublasSetStream(cuwfs[iwfs].handle, cuwfs[iwfs].stream));
	if(powfs[ipowfs].saimcc){
	    if(powfs[ipowfs].nsaimcc>1 || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		cudaMallocHost(&cuwfs[iwfs].imcc, nsa*sizeof(void*));
		for(int isa=0; isa<nsa; isa++){
		    cuwfs[iwfs].imcc[isa]=NULL;
		    gpu_dmat2dev((float**)&(cuwfs[iwfs].imcc[isa]),
				 powfs[ipowfs].saimcc[powfs[ipowfs].nsaimcc>1?wfsind:0]->p[isa]);
		}
	    }else{
		cuwfs[iwfs].imcc=cuwfs[iwfs0].imcc;
	    }
	}
	cuwfs[iwfs].powfs=cupowfs+ipowfs;
	cudaDeviceSynchronize();
	/*GS0 for gtilt. */
	if(powfs[ipowfs].GS0){
	    if(powfs[ipowfs].GS0->nx>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		dsp *t=sptrans(powfs[ipowfs].GS0->p[powfs[ipowfs].GS0->nx>1?wfsind:0]);
		gpu_sp2dev(&cuwfs[iwfs].GS0t, t);
		spfree(t);
	    }else{
		cuwfs[iwfs].GS0t=cuwfs[iwfs0].GS0t;
	    }
	}
	/*wfs amplitude map on loc */
	if(powfs[ipowfs].nlocm>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
	    gpu_dbl2dev(&cuwfs[iwfs].amp, powfs[ipowfs].realamp[powfs[ipowfs].nlocm>1?wfsind:0], powfs[ipowfs].loc->nloc);
	}else{
	    cuwfs[iwfs].amp=cuwfs[iwfs0].amp;
	}
	

	dmat *nea=powfs[ipowfs].neasim->p[wfsind];
	if(nea){
	    gpu_dmat2dev(&cuwfs[iwfs].neasim, nea);
	}

	/* * Now start physical optics setup * */

	if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout){
	    /*If there is llt. */
	    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		if(powfs[ipowfs].llt->ncpa){
		    if(powfs[ipowfs].llt->ncpa->nx>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			gpu_dmat2dev(&cuwfs[iwfs].lltncpa, powfs[ipowfs].llt->ncpa->p[powfs[ipowfs].llt->ncpa->nx>1?wfsind:0]);
		    }else{
			cuwfs[iwfs].lltncpa=cuwfs[iwfs0].lltncpa;
		    }
		}
		if(wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    cudaMallocHost(&cuwfs[iwfs].lltimcc, 1*sizeof(void*));
		    cuwfs[iwfs].lltimcc[0]=NULL;
		    gpu_dmat2dev((float**)&cuwfs[iwfs].lltimcc[0], powfs[ipowfs].llt->imcc->p[0]);
		    gpu_dmat2dev((float**)&cuwfs[iwfs].lltamp, powfs[ipowfs].llt->amp);
		}else{
		    cuwfs[iwfs].lltimcc=cuwfs[iwfs0].lltimcc;
		    cuwfs[iwfs].lltamp=cuwfs[iwfs0].lltamp;
		}
	    }
	    cudaDeviceSynchronize();

	    /*CUFFTW is row major. */
	    int npsf=powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac;/*size of fft */
	    int npsf2[2]={npsf, npsf};
	    const int ncompx=powfs[ipowfs].ncompx;
	    const int ncompy=powfs[ipowfs].ncompy;
	    const int ncompm=MAX(ncompx, ncompy);
	    int ncomp[2]={ncompx, ncompy};
	    int ncompm2[2]={ncompm, ncompm};
	    /*
	      int inembed[2]; inembed[0]=nx; inembed[1]=nx;
	      int istride=1;
	      int idist=nx*nx;
	      DO(cufftPlanMany(&cuwfs[iwfs].plan, 2, nr, 
	      inembed, istride, idist, 
	      inembed, istride, idist, 
	      CUFFT_C2C, nsa));
	    */
	    /*limit the number of subapertures in each batch to less than 1024
	      to save memory. The speed is actually a tiny bit faster for NFIRAOS.*/
	    cuwfs[iwfs].msa=nsa>1024?((int)ceil((float)nsa/(float)(nsa/800))):nsa;
	    if(cufftPlanMany(&cuwfs[iwfs].plan1, 2, npsf2, NULL, 1, 0, NULL, 1, 0, 
			     CUFFT_C2C, cuwfs[iwfs].msa)){
		error("CUFFT plan failed\n");
	    }
	    cufftSetStream(cuwfs[iwfs].plan1, cuwfs[iwfs].stream);
	    if(ncompx==npsf && ncompy==npsf){
		cuwfs[iwfs].plan2=cuwfs[iwfs].plan1;
	    }else{
		if(cufftPlanMany(&cuwfs[iwfs].plan2, 2, ncomp, NULL, 1, 0, NULL, 1, 0, 
				 CUFFT_C2C, cuwfs[iwfs].msa)){
		    error("CUFFT plan failed\n");
		}
		cufftSetStream(cuwfs[iwfs].plan2, cuwfs[iwfs].stream);
	    }
	    if(ncompm==ncompx && ncompm==ncompy){
		cuwfs[iwfs].plan3=cuwfs[iwfs].plan2;
	    }else{
		if(cufftPlanMany(&cuwfs[iwfs].plan3, 2, ncompm2, NULL, 1, 0, NULL, 1, 0, 
				 CUFFT_C2C, cuwfs[iwfs].msa)){
		    error("CUFFT plan failed\n");
		}
		cufftSetStream(cuwfs[iwfs].plan3, cuwfs[iwfs].stream);
	    }

	    if(parms->powfs[ipowfs].llt){
		int nlpsf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
		int nlpsf2[2]={nlpsf, nlpsf};
		if(cufftPlanMany(&cuwfs[iwfs].lltplan_wvf, 2, nlpsf2, NULL, 1,0, NULL, 1, 0, 
				 CUFFT_C2C, 1)){
		    error("CUFFT plan failed\n");
		    cufftSetStream(cuwfs[iwfs].lltplan_wvf, cuwfs[iwfs].stream);
		}
		if(npsf==nlpsf){
		    cuwfs[iwfs].lltplan_otf=cuwfs[iwfs].lltplan_wvf;
		}else{
		    if(cufftPlanMany(&cuwfs[iwfs].lltplan_otf, 2, npsf2, NULL, 1, 0, NULL, 1, 0, 
				     CUFFT_C2C, 1)){
			error("CUFFT plan failed\n");
		    }
		    cufftSetStream(cuwfs[iwfs].lltplan_otf, cuwfs[iwfs].stream);
		}
	    }
	    /*DTF. */
	    if(parms->powfs[ipowfs].usephy){
		if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].llt->n>1 
		   || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    /*Need one per wfs in this powfs, or the first wfs. */
		    int nwvl=parms->powfs[ipowfs].nwvl;
		    cuwfs[iwfs].dtf=(cudtf_t*)calloc(nwvl, sizeof(cudtf_t));
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			int notfused=!powfs[ipowfs].dtf[iwvl].fused;
			if(notfused){
			    cudaCallocHostBlock(cuwfs[iwfs].dtf[iwvl].nominal, nsa*sizeof(void*));
			}
			/*cudaCallocHostBlock(cuwfs[iwfs].dtf[iwvl].si, nsa*sizeof(void*)); */
			int multi_nominal=(powfs[ipowfs].dtf[iwvl].si->nx==nsa);
			for(int isa=0; isa<nsa; isa++){
			    if(multi_nominal || isa==0){
				if(notfused){
				    gpu_cmat2dev(&cuwfs[iwfs].dtf[iwvl].nominal[isa], 
						 powfs[ipowfs].dtf[iwvl].nominal->p[isa+nsa*(powfs[ipowfs].dtf[iwvl].nominal->ny>1?wfsind:0)]);
				}
			    }else{
				cuwfs[iwfs].dtf[iwvl].nominal[isa]=cuwfs[iwfs].dtf[iwvl].nominal[0];
			    }
			}
		    
			if(parms->powfs[ipowfs].llt){
			    cudaCallocHostBlock(cuwfs[iwfs].dtf[iwvl].etf, nsa*sizeof(void*));
			    cmat *(*petf)[nsa]=NULL;
			    if(powfs[ipowfs].etfsim[iwvl].p1){
				petf=(cmat *(*)[nsa])powfs[ipowfs].etfsim[iwvl].p1->p;
				cuwfs[iwfs].dtf[iwvl].etfis1d=1;
			    }else{
				petf=(cmat *(*)[nsa])powfs[ipowfs].etfsim[iwvl].p2->p;
				cuwfs[iwfs].dtf[iwvl].etfis1d=0;
			    }
			    for(int isa=0; isa<nsa; isa++){
				gpu_cmat2dev(&cuwfs[iwfs].dtf[iwvl].etf[isa], petf[parms->powfs[ipowfs].llt->n>1?wfsind:0][isa]);
			    }
			}
		    }/*for iwvl. */
		    if(parms->powfs[ipowfs].llt){
			gpu_dmat2dev(&cuwfs[iwfs].srot, powfs[ipowfs].srot->p[parms->powfs[ipowfs].llt->n>1?wfsind:0]);
		    }
		}else{
		    cuwfs[iwfs].dtf  = cuwfs[iwfs0].dtf;
		    cuwfs[iwfs].srot = cuwfs[iwfs0].srot;
		}
		/*Matched filter */
		if(parms->powfs[ipowfs].phytypesim==1){
		    if(powfs[ipowfs].intstat->mtche->ny>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cudaCallocHostBlock(cuwfs[iwfs].mtche, nsa*sizeof(void*));
			dmat **mtche=powfs[ipowfs].intstat->mtche->p+nsa*(powfs[ipowfs].intstat->mtche->ny>1?wfsind:0);
			for(int isa=0; isa<nsa; isa++){
			    gpu_dmat2dev((float**)&cuwfs[iwfs].mtche[isa], mtche[isa]);
			}
		    }else{
			cuwfs[iwfs].mtche=cuwfs[iwfs0].mtche;
		    }
		}
		if(powfs[ipowfs].bkgrnd){
		    if(powfs[ipowfs].bkgrnd->ny==1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cudaCallocHostBlock(cuwfs[iwfs].bkgrnd2, nsa*sizeof(void*));
			dmat **bkgrnd=powfs[ipowfs].bkgrnd->p+nsa*(powfs[ipowfs].bkgrnd->ny==1?wfsind:0);
			for(int isa=0; isa<nsa; isa++){
			    gpu_dmat2dev((float**)&cuwfs[iwfs].bkgrnd2[isa], bkgrnd[isa]);
			}
		    }else{
			cuwfs[iwfs].bkgrnd2=cuwfs[iwfs0].bkgrnd2;
		    }
		}
		if(powfs[ipowfs].bkgrndc){
		    if(powfs[ipowfs].bkgrndc->ny==1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cudaCallocHostBlock(cuwfs[iwfs].bkgrnd2c, nsa*sizeof(void*));
			dmat **bkgrndc=powfs[ipowfs].bkgrndc->p+nsa*(powfs[ipowfs].bkgrndc->ny==1?wfsind:0);
			for(int isa=0; isa<nsa; isa++){
			    gpu_dmat2dev((float**)&cuwfs[iwfs].bkgrnd2c[isa], bkgrndc[isa]);
			}
		    }else{
			cuwfs[iwfs].bkgrnd2c=cuwfs[iwfs0].bkgrnd2c;
		    }	
		}
	    }
	}/*if phy */
	CUDA_SYNC_DEVICE;
	gpu_print_mem("wfs init");
    }/*for iwfs */
}
void gpu_wfs_init_sim(const PARMS_T *parms, POWFS_T *powfs){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwfs_t *cuwfs=cudata->wfs;
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	if(parms->powfs[ipowfs].phystep!=0
	   ||parms->save.gradgeom[iwfs]
	   ||parms->powfs[ipowfs].pistatout){
	    /*gradacc for geom wfs accumulation */
	    curfree(cuwfs[iwfs].gradacc);
	    cuwfs[iwfs].gradacc=curnew(nsa*2,1);
	}
	if(parms->powfs[ipowfs].usephy){
	    curcellfree(cuwfs[iwfs].ints);
	    cuwfs[iwfs].ints=curcellnew(nsa,1,powfs[ipowfs].pixpsax,powfs[ipowfs].pixpsay);
	    if(parms->powfs[ipowfs].noisy){
		curfree(cuwfs[iwfs].neareal);
		cuwfs[iwfs].neareal=curnew(nsa*4,1);
	    }
	}
	if(parms->powfs[ipowfs].pistatout){
	    if(parms->powfs[ipowfs].pistatstc){
		error("pistatstc is not supported yet.\n");
	    }
	    curcellfree(cuwfs[iwfs].pistatout);
	    const int notfx=powfs[ipowfs].ncompx;/*necessary size to build detector image. */
	    const int notfy=powfs[ipowfs].ncompy;
	    const int npsf=MAX(notfx,notfy);
	    cuwfs[iwfs].pistatout=curcellnew(nsa, parms->powfs[ipowfs].nwvl, npsf, npsf);
	}
	CUDA_SYNC_DEVICE;
    }
}
void gpu_wfssurf2gpu(const PARMS_T *parms, POWFS_T *powfs){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(GPUS[wfsgpu[iwfs]]);
	cuwfs_t *cuwfs=cudata->wfs;
	int ipowfs=parms->wfs[iwfs].powfs;
	int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	if(powfs[ipowfs].ncpa){
	    gpu_dmat2cu(&cuwfs[iwfs].opdadd, powfs[ipowfs].ncpa->p[wfsind]);
	}else{
	    curzero(cuwfs[iwfs].opdadd, cuwfs[iwfs].stream);
	}
	if(powfs[ipowfs].opdadd && powfs[ipowfs].opdadd->p[wfsind]){
	    curmat *temp=NULL;
	    gpu_dmat2cu(&temp, powfs[ipowfs].opdadd->p[wfsind]);
	    curadd(&cuwfs[iwfs].opdadd, 1, temp, 1, cuwfs[iwfs].stream);
	    curfree(temp);
	    dfree(powfs[ipowfs].opdadd->p[wfsind]);/*no longer need it in CPU memory. */
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
void gpu_wfsgrad_seeding(const PARMS_T *parms, const POWFS_T *powfs, rand_t *rstat){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(GPUS[wfsgpu[iwfs]]);
	cuwfs_t *cuwfs=cudata->wfs;
	int seed=lrand(rstat);/*don't put this after continue. */
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].noisy) continue;
	int nsa=powfs[ipowfs].pts->nsa*2;
	if(nsa<RAND_THREAD){
	    cuwfs[iwfs].custatt=nsa;
	    cuwfs[iwfs].custatb=1;
	}else if(nsa<RAND_THREAD*RAND_BLOCK){
	    cuwfs[iwfs].custatt=RAND_THREAD;
	    cuwfs[iwfs].custatb=nsa/RAND_THREAD+(nsa%RAND_THREAD)?1:0;
	}else{
	    cuwfs[iwfs].custatt=RAND_THREAD;
	    cuwfs[iwfs].custatb=RAND_BLOCK;
	}
	cudaMalloc(&cuwfs[iwfs].custat, (cuwfs[iwfs].custatt*cuwfs[iwfs].custatb)*sizeof(curandStat));
	setup_rand<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt>>>(cuwfs[iwfs].custat, seed);
    }
    CUDA_SYNC_DEVICE;
    gpu_print_mem("wfs seeding");
}
