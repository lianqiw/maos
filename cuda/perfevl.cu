extern "C"
{
#include <cuda.h>
#include "gpu.h"
#include "utils.h"
#include "accphi.h"
}
static float (*cuplocs)[2]=NULL;
static float *cupamp=NULL;

/**
  save aper_locs, aper_amp to GPU.
*/
static void gpu_plocs2gpu(loc_t *loc, dmat *amp){
    if(cuplocs) error("plocs is already Copied.\n");
    int nloc=loc->nloc;
    gpu_loc2dev(&cuplocs, loc);
    gpu_dbl2dev(&cupamp, amp->p, nloc);
}
/**
   Initialize perfevl
*/
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper){
    (void)parms;
    gpu_plocs2gpu(aper->locs, aper->amp);
}
/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
 */
void gpu_perfevl(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int ievl=info->start;
    assert(info->end==info->start+1);//only one evl.
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int do_psf=(parms->evl.psfmean || parms->evl.psfhist) && isim>=parms->evl.psfisim;
    const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
    const int nloc=aper->locs->nloc;
    const double thetax=parms->evl.thetax[ievl];
    const double thetay=parms->evl.thetay[ievl];
    //Setup pointers for easy usage
    PDMAT(simu->olmp->p[ievl],polmp);//OL mode for each dir
    PDMAT(simu->olep->p[ievl],polep);//OL error for each dir
    PDMAT(simu->clmp->p[ievl],pclmp);
    PDMAT(simu->clep->p[ievl],pclep);

    float *iopdevl; 
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    cudaCalloc(iopdevl, aper->locs->nloc*sizeof(float));
    cudaStream_t stream;
    STREAM_NEW(stream);

    if(parms->sim.idealevl){
	error("Please finished by: \n"
	      "1) make aloc square, \n"
	      "2) send dmproj to this file by calling gpu_dm2gpu\n");
    }else if(simu->atm && !parms->sim.wfsalias){
	gpu_atm2loc(iopdevl, cuplocs, nloc, parms->evl.hs[ievl], thetax, thetay, 
		parms->evl.misreg[0], parms->evl.misreg[1], isim*dt, 1, stream);
    }
    CUDA_SYNC_STREAM;
    if(simu->telws){//Wind shake
	TO_IMPLEMENT;
    }
    if(simu->surfevl && simu->surfevl->p[ievl]){
	TO_IMPLEMENT;
    }
    if(save_evlopd){
	TO_IMPLEMENT;
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(nmod==3){
	gpu_calc_ptt(polep[isim], polmp[isim], aper->ipcc, aper->imcc,
		 cuplocs, nloc, iopdevl, cupamp, stream);
    }else{
	TO_IMPLEMENT;
    }
  
    if(parms->evl.psfmean &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
			     ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	TO_IMPLEMENT;
    }
    
    if(parms->sim.evlol) goto end;
    
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	gpu_dm2loc(iopdevl, cuplocs, nloc, parms->evl.hs[ievl], thetax, thetay,
	       parms->evl.misreg[0], parms->evl.misreg[1], -1, stream);
	CUDA_SYNC_STREAM;
	if(imoao>-1){
	    TO_IMPLEMENT;
	}
    }
    if(save_evlopd){
	TO_IMPLEMENT;
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(parms->tomo.split){
	if(parms->ndm<=2){
	    PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
	    if(nmod==3){
		gpu_calc_ngsmod(pclep[isim], pclmp[isim], pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl, cupamp, stream);
	    }else{
		gpu_calc_ngsmod(NULL, NULL, pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl, cupamp, stream);	
		TO_IMPLEMENT;
	    }
	}
    }else{
	if(nmod==3){
	    gpu_calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,
			 cuplocs, nloc, iopdevl, cupamp, stream);
	}else{
	    TO_IMPLEMENT;
	}
    }
    if(parms->evl.psf[ievl] && isim>=parms->evl.psfisim && do_psf){
	TO_IMPLEMENT;
    }
 end:
    STREAM_DONE(stream);
    cudaFree(iopdevl);
}
