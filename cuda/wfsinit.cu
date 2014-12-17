/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include "cucmat.h"
#include "wfs.h"
#include "cudata.h"

/**
   Reshape each cell into column vector and concatenate each cell into a column
*/
static cmat *concat_ccell_as_vector(ccell *input){
    if(input->ny!=1){
	error("Invalid format\n");
    }
    int npix=input->p[0]->nx*input->p[0]->ny;
    int nsa=input->nx;
    cmat *output=cnew(npix, nsa);
    for(long isa=0; isa<nsa; isa++){
	memcpy(output->p+isa*npix, input->p[isa]->p, sizeof(dcomplex)*npix);
    }
    return output;
}
/**
   Reshape each cell into column vector and concatenate each cell into a column
*/
static dmat *concat_dcell_as_vector(dcell *input){
    if(input->ny!=1){
	error("Invalid format\n");
    }
    int npix=input->p[0]->nx*input->p[0]->ny;
    int nsa=input->nx;
    dmat *output=dnew(npix, nsa);
    for(long isa=0; isa<nsa; isa++){
	memcpy(output->p+isa*npix, input->p[isa]->p, sizeof(double)*npix);
    }
    return output;
}
static void etf2gpu(cucmat **cuetf, ETF_T *etf, int icol, int *etfis1d){
    ccell *etfc=0;
    if(etf->p1){
	etfc=ccellsub(etf->p1, 0, 0, icol, 1);
	*etfis1d=1;
    }else{
	etfc=ccellsub(etf->p2, 0, 0, icol, 1);
	*etfis1d=0;
    }
    cmat *etfm=concat_ccell_as_vector(etfc);
    cp2gpu(cuetf, etfm);
    cfree(etfm);
    ccellfree(etfc);
}
/**
   Initialize or update etf.
*/
void gpu_wfsgrad_update_etf(const PARMS_T *parms, const POWFS_T *powfs){
    const int *wfsgpu=cudata_t::wfsgpu;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwfs_t *cuwfs=cudata_t::wfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int nwvl=parms->powfs[ipowfs].nwvl;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
	if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout){
	    if(parms->powfs[ipowfs].usephy){
		if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].llt->n>1 || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    if(parms->powfs[ipowfs].llt){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
			    int icol=parms->powfs[ipowfs].llt->n>1?wfsind:0;
			    if(powfs[ipowfs].etfsim){
				etf2gpu(&cuwfs[iwfs].dtf[iwvl].etf, &powfs[ipowfs].etfsim[iwvl], icol, &cuwfs[iwfs].dtf[iwvl].etfis1d);
			    }
			    if(powfs[ipowfs].etfsim2){
				etf2gpu(&cuwfs[iwfs].dtf[iwvl].etf2, &powfs[ipowfs].etfsim2[iwvl], icol, &cuwfs[iwfs].dtf[iwvl].etfis1d);
			    }
			}
		    }
		}
	    }
	}
    }
}
void gpu_wfsgrad_update_mtche(const PARMS_T *parms, const POWFS_T *powfs){
    const int *wfsgpu=cudata_t::wfsgpu;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwfs_t *cuwfs=cudata_t::wfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
	const int nsa=powfs[ipowfs].pts->nsa;
	if(parms->powfs[ipowfs].usephy){
	    if(parms->powfs[ipowfs].phytypesim==1){
		if(powfs[ipowfs].intstat->mtche->ny>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    int icol=powfs[ipowfs].intstat->mtche->ny>1?wfsind:0;
		    dcell *mtchec=dcellsub(powfs[ipowfs].intstat->mtche, 0, 0, icol, 1);
		    dmat *mtche=concat_dcell_as_vector(mtchec);
		    dcellfree(mtchec);
		    if(iwfs!=iwfs0 && cuwfs[iwfs].mtche==cuwfs[iwfs0].mtche){
			info("Reset mtche to 0\n");
			cuwfs[iwfs].mtche=0;
			cuwfs[iwfs].i0sum=0;
		    }
		    cp2gpu(&cuwfs[iwfs].mtche, mtche);
		    dfree(mtche);
		    cp2gpu(&cuwfs[iwfs].i0sum, powfs[ipowfs].intstat->i0sum->p+nsa*icol, nsa, 1);
		}else{
		    cuwfs[iwfs].mtche=cuwfs[iwfs0].mtche;
		    cuwfs[iwfs].i0sum=cuwfs[iwfs0].i0sum;
		}
	    }
	}
    }
}
/**
   Initialize or update mtched filter
*/

/**
   Initialize other arrays
*/
void gpu_wfsgrad_init(const PARMS_T *parms, const POWFS_T *powfs){
    const int *wfsgpu=cudata_t::wfsgpu;
    cudata_t::wfs=(cuwfs_t*)calloc(parms->nwfs, sizeof(cuwfs_t));
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->powfs=(cuwloc_t*)calloc(parms->npowfs, sizeof(cuwloc_t));
	cuwloc_t *cupowfs=cudata->powfs;
	/* Setup information that are same for wfs in each powfs*/
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs==0) continue;
	    pts_t *pts=powfs[ipowfs].pts;
	    loc_t *loc=powfs[ipowfs].loc;
	    cupowfs[ipowfs].pts=new cupts_t(pts);
	    cupowfs[ipowfs].loc=new culoc_t(loc);
	    cupowfs[ipowfs].saloc=new culoc_t(powfs[ipowfs].saloc);
	
	    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		pts=powfs[ipowfs].llt->pts;
		loc=powfs[ipowfs].llt->loc;
		cupowfs[ipowfs].llt=new cullt_t;
		cupowfs[ipowfs].llt->pts=new cupts_t(pts);
		cupowfs[ipowfs].llt->loc=new culoc_t(loc);
	    }
	    /*cupowfs[ipowfs].skip=parms->powfs[ipowfs].skip; */
	    if(parms->powfs[ipowfs].fieldstop){
		cp2gpu(&cupowfs[ipowfs].embed, powfs[ipowfs].fieldstop->embed->p[0]->p, powfs[ipowfs].loc->nloc, 1);
		cupowfs[ipowfs].nembed=powfs[ipowfs].fieldstop->nembed->p[0];
		cp2gpu(&cupowfs[ipowfs].fieldstop, powfs[ipowfs].fieldstop->fieldmask->p[0]);
	    }
	}
    }

    /* setup information that maybe different for wfs in same powfs due to
       misregistration or NCPA.*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwloc_t *cupowfs=cudata->powfs;
	cuwfs_t *cuwfs=cudata_t::wfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int nsa=powfs[ipowfs].pts->nsa;
	const int nwvl=parms->powfs[ipowfs].nwvl;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
	const int nwfsp=parms->powfs[ipowfs].nwfs;
	const int ndm=parms->ndm;
	/*imcc for ztilt. */
	cuwfs[iwfs].stream=new stream_t;
	cuwfs[iwfs].loc_dm=new culoc_t*[ndm];
	for(int idm=0; idm<ndm; idm++){
	    if(powfs[ipowfs].loc_dm){
		cuwfs[iwfs].loc_dm[idm]=new culoc_t(powfs[ipowfs].loc_dm->p[wfsind+idm*nwfsp]);
	    }else{
		cuwfs[iwfs].loc_dm[idm]=new culoc_t(powfs[ipowfs].loc);
	    }
	}
	if(powfs[ipowfs].loc_tel){
	    cuwfs[iwfs].loc_tel=new culoc_t(powfs[ipowfs].loc_tel->p[wfsind]);
	}else{
	    cuwfs[iwfs].loc_tel=new culoc_t(powfs[ipowfs].loc);
	}
	cuwfs[iwfs].phiout=curnew(powfs[ipowfs].loc->nloc, 1);

	if(parms->powfs[ipowfs].fieldstop){
	    DO(cufftPlan2d(&cuwfs[iwfs].plan_fs, cupowfs[ipowfs].nembed, cupowfs[ipowfs].nembed, FFT_T_C2C));
	    cufftSetStream(cuwfs[iwfs].plan_fs, *cuwfs[iwfs].stream);
	}
	if(powfs[ipowfs].saimcc){
	    if(powfs[ipowfs].saimcc->nx>1 || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		Real *imcc[nsa];
		for(int isa=0; isa<nsa; isa++){
		    imcc[isa]=NULL;
		    cp2gpu(&(imcc[isa]),
			   powfs[ipowfs].saimcc->p[powfs[ipowfs].saimcc->nx>1?wfsind:0]->p[isa]);
		}
		DO(cudaMalloc(&cuwfs[iwfs].imcc, nsa*sizeof(void*)));
		DO(cudaMemcpy(cuwfs[iwfs].imcc, imcc, nsa*sizeof(void*), cudaMemcpyHostToDevice));
	    }else{
		cuwfs[iwfs].imcc=cuwfs[iwfs0].imcc;
	    }
	}
	cuwfs[iwfs].powfs=cupowfs+ipowfs;
	cudaDeviceSynchronize();
	/*GS0 for gtilt. */
	if(powfs[ipowfs].GS0){
	    dsp *t=powfs[ipowfs].GS0->p[powfs[ipowfs].GS0->nx>1?wfsind:0];
	    cuwfs[iwfs].GS0=new cusp(t, 1);
	}
	/*wfs amplitude map on loc */
	cp2gpu(&cuwfs[iwfs].amp, powfs[ipowfs].realamp->p[wfsind]);
	dmat *nea=powfs[ipowfs].neasim->p[wfsind];
	if(nea){
	    cp2gpu(&cuwfs[iwfs].neasim, nea);
	}

	/* * Now start physical optics setup * */

	if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout){
	    /*If there is llt. */
	    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		if(powfs[ipowfs].llt->ncpa){
		    if(powfs[ipowfs].llt->ncpa->nx>1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cp2gpu(&cuwfs[iwfs].lltncpa, 
			       powfs[ipowfs].llt->ncpa->p[powfs[ipowfs].llt->ncpa->nx>1?wfsind:0]);
		    }else{
			cuwfs[iwfs].lltncpa=cuwfs[iwfs0].lltncpa;
		    }
		}
		if(wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    DO(cudaMallocHost(&cuwfs[iwfs].lltimcc, 1*sizeof(void*)));
		    cuwfs[iwfs].lltimcc[0]=NULL;
		    cp2gpu((Real**)&cuwfs[iwfs].lltimcc[0], powfs[ipowfs].llt->imcc->p[0]);
		    cp2gpu((Real**)&cuwfs[iwfs].lltamp, powfs[ipowfs].llt->amp);
		}else{
		    cuwfs[iwfs].lltimcc=cuwfs[iwfs0].lltimcc;
		    cuwfs[iwfs].lltamp=cuwfs[iwfs0].lltamp;
		}
	    }
	    cudaDeviceSynchronize();

	    /*CUFFTW is row major. */
	    int nwvf=powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac;/*size of fft */
	    int nwvf2[2]={nwvf, nwvf};
	    const int ncompx=powfs[ipowfs].ncompx;
	    const int ncompy=powfs[ipowfs].ncompy;
	    const int notf=MAX(ncompx, ncompy);
	    int ncomp2[2]={ncompx, ncompy};
	    int notf2[2]={notf, notf};
	    /*limit the number of subapertures in each batch to less than 1024
	      to save memory. The speed is actually a tiny bit faster for NFIRAOS.*/
	    cuwfs[iwfs].msa=nsa>1024?((int)ceil((Real)nsa/(Real)(nsa/800))):nsa;
	    if(cufftPlanMany(&cuwfs[iwfs].plan1, 2, nwvf2, NULL, 1, 0, NULL, 1, 0, 
			     FFT_T_C2C, cuwfs[iwfs].msa)){
		error("CUFFT plan failed\n");
	    }
	    cufftSetStream(cuwfs[iwfs].plan1, *cuwfs[iwfs].stream);

	    if(notf==nwvf){
		cuwfs[iwfs].plan2=cuwfs[iwfs].plan1;
	    }else{
		if(cufftPlanMany(&cuwfs[iwfs].plan2, 2, notf2, NULL, 1, 0, NULL, 1, 0, 
				 FFT_T_C2C, cuwfs[iwfs].msa)){
		    error("CUFFT plan failed\n");
		}
		cufftSetStream(cuwfs[iwfs].plan2, *cuwfs[iwfs].stream);
	    }
	    if(notf==ncompx && notf==ncompy){
		cuwfs[iwfs].plan3=cuwfs[iwfs].plan2;
	    }else{
		if(cufftPlanMany(&cuwfs[iwfs].plan3, 2, ncomp2, NULL, 1, 0, NULL, 1, 0, 
				 FFT_T_C2C, cuwfs[iwfs].msa)){
		    error("CUFFT plan failed\n");
		}
		cufftSetStream(cuwfs[iwfs].plan3, *cuwfs[iwfs].stream);
	    }

	    if(parms->powfs[ipowfs].llt){
		int nlwvf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
		int nlwvf2[2]={nlwvf, nlwvf};
		if(cufftPlanMany(&cuwfs[iwfs].lltplan_wvf, 2, nlwvf2, NULL, 1,0, NULL, 1, 0, 
				 FFT_T_C2C, 1)){
		    error("CUFFT plan failed\n");
		}
		cufftSetStream(cuwfs[iwfs].lltplan_wvf, *cuwfs[iwfs].stream);
		if(notf==nlwvf){
		    cuwfs[iwfs].lltplan_otf=cuwfs[iwfs].lltplan_wvf;
		}else{
		    if(cufftPlanMany(&cuwfs[iwfs].lltplan_otf, 2, notf2, NULL, 1, 0, NULL, 1, 0, 
				     FFT_T_C2C, 1)){
			error("CUFFT plan failed\n");
		    }
		    cufftSetStream(cuwfs[iwfs].lltplan_otf, *cuwfs[iwfs].stream);
		}
	    }
	    /*DTF. */
	    if(parms->powfs[ipowfs].usephy){
		if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].llt->n>1 
		   || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
		    /*Need one per wfs in this powfs, or the first wfs. */
		    cuwfs[iwfs].dtf=(cudtf_t*)calloc(nwvl, sizeof(cudtf_t));
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			int notfused=!powfs[ipowfs].dtf[iwvl].fused;
			if(notfused){
			    int icol=powfs[ipowfs].dtf[iwvl].nominal->ny>1?wfsind:0;
			    ccell *nominalc=ccellsub(powfs[ipowfs].dtf[iwvl].nominal, 0, 0, icol, 1);
			    cmat *nominal=concat_ccell_as_vector(nominalc);
			    ccellfree(nominalc);
			    cp2gpu(&cuwfs[iwfs].dtf[iwvl].nominal, nominal);
			    cfree(nominal);
			}
			//ETF moved to gpu_wfsgrad_update_etf();
		    }/*for iwvl. */
		    if(parms->powfs[ipowfs].llt){
			cp2gpu(&cuwfs[iwfs].srot, powfs[ipowfs].srot->p[parms->powfs[ipowfs].llt->n>1?wfsind:0]);
		    }
		}else{
		    cuwfs[iwfs].dtf  = cuwfs[iwfs0].dtf;
		    cuwfs[iwfs].srot = cuwfs[iwfs0].srot;
		}
		/*Matched filter */
		if(parms->powfs[ipowfs].phytypesim==1){
		    //Separated with gpu_wfsgrad_upate_mtche();
		}else if(parms->powfs[ipowfs].phytypesim==2){/*cog*/
		    if(powfs[ipowfs].intstat->cogcoeff->nx>1 || wfsind==0 || wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cp2gpu(&cuwfs[iwfs].cogcoeff, 
			       powfs[ipowfs].intstat->cogcoeff->p[powfs[ipowfs].intstat->cogcoeff->nx>1?wfsind:0]->p, nsa*2, 1);
		    }else{
			cuwfs[iwfs].cogcoeff=cuwfs[iwfs0].cogcoeff;
		    }
		}
		if(powfs[ipowfs].bkgrnd){
		    if(powfs[ipowfs].bkgrnd->ny==1 || wfsind==0|| wfsgpu[iwfs]!=wfsgpu[iwfs0]){
			cudaCallocHostBlock(cuwfs[iwfs].bkgrnd2, nsa*sizeof(void*));
			dmat **bkgrnd=powfs[ipowfs].bkgrnd->p+nsa*(powfs[ipowfs].bkgrnd->ny==1?wfsind:0);
			for(int isa=0; isa<nsa; isa++){
			    cp2gpu((Real**)&cuwfs[iwfs].bkgrnd2[isa], bkgrnd[isa]);
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
			    cp2gpu((Real**)&cuwfs[iwfs].bkgrnd2c[isa], bkgrndc[isa]);
			}
		    }else{
			cuwfs[iwfs].bkgrnd2c=cuwfs[iwfs0].bkgrnd2c;
		    }	
		}
		if(parms->powfs[ipowfs].dither){
		    cuwfs[iwfs].dither=new dither_t(nsa,powfs[ipowfs].pixpsax,powfs[ipowfs].pixpsay);
		}
	    }
	    const int msa=cuwfs[iwfs].msa;
	    cuwfs[iwfs].wvf=cucnew(nwvf*nwvf,msa);
	    if(nwvf!=notf){
		cuwfs[iwfs].psf=cucnew(notf*notf,msa);
	    }
	    if(parms->powfs[ipowfs].radrot || ncompx!=notf || ncompy!=notf){
		cuwfs[iwfs].otf=cucnew(ncompx*ncompy,msa);
	    }
	    if(parms->powfs[ipowfs].psfout){
		const int wvf_n=notf/2+2;
		cuwfs[iwfs].wvfout=cuccellnew(nsa, nwvl, wvf_n, wvf_n);
		cuwfs[iwfs].psfout=cucnew(notf*notf, msa);
	    }
	    if(parms->powfs[ipowfs].pistatout){
		cuwfs[iwfs].psfstat=cucnew(notf*notf, msa);
	    }
	    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		int nlx=powfs[ipowfs].llt->pts->nx;
		int nlwvf=nlx*parms->powfs[ipowfs].embfac;
		cuwfs[iwfs].lltopd=curnew(nlx, nlx);
		if(parms->powfs[ipowfs].pistatout || parms->sim.uptideal){
		    DO(cudaMallocHost(&cuwfs[iwfs].lltg, 2*sizeof(Real)));
		}
		cuwfs[iwfs].lltwvf=cucnew(nlwvf, nlwvf);
		if(nlwvf!=notf){
		    cuwfs[iwfs].lltotfc=cucnew(notf, notf);
		}
	    }
	
	
	}/*if phy */
	CUDA_SYNC_DEVICE;
	gpu_print_mem("wfs init");
    }/*for iwfs */
    gpu_wfsgrad_update_etf(parms, powfs);
    gpu_wfsgrad_update_mtche(parms, powfs);
}
void gpu_wfs_init_sim(const PARMS_T *parms, POWFS_T *powfs){
    int *wfsgpu=cudata_t::wfsgpu;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
	cuwfs_t *cuwfs=cudata_t::wfs;
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	//gradacc is used for accumulation in geom mode and for output in phy mode
	initzero(&cuwfs[iwfs].gradacc, nsa*2, 1);
	initzero(&cuwfs[iwfs].gradcalc,nsa*2,1);
	if(parms->powfs[ipowfs].usephy || parms->powfs[ipowfs].dither){
	    if(!cuwfs[iwfs].ints){
		cuwfs[iwfs].ints=curcellnew(nsa,1,powfs[ipowfs].pixpsax,powfs[ipowfs].pixpsay);
	    }else{
		curcellzero(cuwfs[iwfs].ints);
	    }
	}
	if(parms->powfs[ipowfs].pistatout){
	    if(parms->powfs[ipowfs].pistatstc){
		error("pistatstc is not supported yet.\n");
	    }
	    if(!cuwfs[iwfs].pistatout){
		const int notfx=powfs[ipowfs].ncompx;/*necessary size to build detector image. */
		const int notfy=powfs[ipowfs].ncompy;
		const int npsf=MAX(notfx,notfy);
		cuwfs[iwfs].pistatout=curcellnew(nsa, parms->powfs[ipowfs].nwvl, npsf, npsf);
	    }else{
		curcellzero(cuwfs[iwfs].pistatout);
	    }
	}
	CUDA_SYNC_DEVICE;
    }
}
void gpu_wfssurf2gpu(const PARMS_T *parms, POWFS_T *powfs){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(cudata_t::wfsgpu[iwfs]);
	cuwfs_t *cuwfs=cudata_t::wfs;
	int ipowfs=parms->wfs[iwfs].powfs;
	int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	if(powfs[ipowfs].opdadd && powfs[ipowfs].opdadd->p[wfsind]){
	    cp2gpu(&cuwfs[iwfs].opdadd, powfs[ipowfs].opdadd->p[wfsind]);
	    dfree(powfs[ipowfs].opdadd->p[wfsind]);/*no longer need it in CPU memory. */
	}
    }
}
__global__ static void setup_rand(curandState *rstat, int seed){
    int id=threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, id, 0, &rstat[id]);
}
/**
   Seed the random number genrator
*/
void gpu_wfsgrad_seeding(const PARMS_T *parms, const POWFS_T *powfs, rand_t *rstat){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	gpu_set(cudata_t::wfsgpu[iwfs]);
	cuwfs_t *cuwfs=cudata_t::wfs;
	int seed=lrand(rstat);/*don't put this after continue. */
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].noisy) continue;
	int nsa=powfs[ipowfs].pts->nsa*2;
	if(nsa<RAND_THREAD){
	    cuwfs[iwfs].custatt=nsa;//number of threads
	    cuwfs[iwfs].custatb=1;//number of blocks
	}else if(nsa<RAND_THREAD*RAND_BLOCK){
	    cuwfs[iwfs].custatt=RAND_THREAD;
	    cuwfs[iwfs].custatb=nsa/RAND_THREAD+(nsa%RAND_THREAD)?1:0;
	}else{
	    cuwfs[iwfs].custatt=RAND_THREAD;
	    cuwfs[iwfs].custatb=RAND_BLOCK;
	}
	DO(cudaMalloc(&cuwfs[iwfs].custat, (cuwfs[iwfs].custatt*cuwfs[iwfs].custatb)*sizeof(curandState)));
	setup_rand<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt>>>(cuwfs[iwfs].custat, seed);
    }
    CUDA_SYNC_DEVICE;
    gpu_print_mem("wfs seeding");
}
