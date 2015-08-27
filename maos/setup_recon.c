/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "common.h"
#include "setup_recon.h"
#include "recon.h"
#include "fdpcg.h"
#include "ahst.h"
#include "recon_utils.h"
#include "moao.h"
#include "setup_powfs.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

/**
   \file setup_recon.c Contains routines that setup the wavefront reconstructor
   and DM fitting.  Use parms->wfsr instead of parms->wfs for wfs information,
   which hands GLAO mode correctly.x 

   TOMOSCALE is used for RLM, RRM, in MVST for M, and in CUDA due to
   limited dynamic range of single precision floating point numbers.

   2014-04-01: Scale saneai, cxx, zzt by TOMOSCALE.

   All routines in this file depends on saneai and may be called repeatedly during simulation.
*/


/**
   Setup the matrix of the inverse of gradient measurement noise equivalent
   angle covariance matrix. For physical optics wfs, the NEA is computed using
   matched filter output. For geometric optics, the NEA is from input.
*/
static void
setup_recon_saneai(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    const int nwfs=parms->nwfsr;
    dspcellfree(recon->sanea);
    dspcellfree(recon->saneai);
    dspcellfree(recon->saneal);
    dspcell *sanea=recon->sanea=cellnew(nwfs,nwfs);//The subaperture NEA
    dspcell *saneal=recon->saneal=cellnew(nwfs,nwfs);//The LL' decomposition of sanea
    dspcell *saneai=recon->saneai=cellnew(nwfs,nwfs);//The inverse of sanea
    info2("Recon NEA:\n");
    
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip==3) continue;
	int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs->p[0];
	int nsa=powfs[ipowfs].saloc->nloc;
	int do_ref=0;
	lmat *samask=0;
	if(parms->wfs[iwfs].sabad){
	    samask=loc_coord2ind(powfs[ipowfs].saloc, parms->wfs[iwfs].sabad);
	}
	if(parms->powfs[ipowfs].neareconfile){/*load sanea from file */
	    dmat *nea=dread("%s_wfs%d",parms->powfs[ipowfs].neareconfile,iwfs);/*rad */
	    if(nea && nea->p[0]<1e-11) {
		error("Supplied NEA is too small\n");
	    }
	    if(samask){
		for(int isa=0; isa<nsa; isa++){
		    if(samask->p[isa]){
			warning("wfs %d sa %d is masked\n", iwfs, isa);
			nea->p[isa]=INFINITY;
			nea->p[isa+nsa]=INFINITY;
		    }
		}
	    }
	    /*rad */
	    saneal->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea, 2);/*rad^2 */
	    sanea->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea,-1);/*rad^-2 */
	    saneai->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
	    dfree(nea);
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && !parms->powfs[ipowfs].phyusenea){
	    /*Physical optics use nea from intstat*/
	    if(!powfs[ipowfs].saneaxy){
		error("saneaxy cannot be null\n");
	    }
	    dcell *saneaxy=powfs[ipowfs].saneaxy;
	    if(parms->recon.glao && saneaxy->ny!=1){
		error("Please average powfs[ipowfs].saneaxy for GLAO mode.\n");
	    }
	    if(saneaxy->ny>1 || iwfs==iwfs0 || parms->wfs[iwfs].sabad){
		int ind=saneaxy->ny>1?parms->powfs[ipowfs].wfsind->p[iwfs]:0;
		dcell *saneaxyi=cellnew(nsa, 1);
		dcell *saneaxyl=cellnew(nsa, 1);
		dcell *saneaixy=cellnew(nsa, 1);
		for(int isa=0; isa<nsa; isa++){
		    dcp(&saneaxyi->p[isa], IND(saneaxy, isa, ind));
		    saneaxyl->p[isa]=dchol(saneaxyi->p[isa]);
		    saneaixy->p[isa]=dinvspd(saneaxyi->p[isa]);
		    if(samask && samask->p[isa]){
			warning("wfs %d sa %d is masked\n", iwfs, isa);
			dset(saneaxyi->p[isa], INFINITY);
			dset(saneaxyl->p[isa], INFINITY);
			dset(saneaixy->p[isa], 0);
		    }
		}
		sanea->p[iwfs+iwfs*nwfs]=nea2sp(saneaxyi->p, nsa);
		saneal->p[iwfs+iwfs*nwfs]=nea2sp(saneaxyl->p,nsa);
		saneai->p[iwfs+iwfs*nwfs]=nea2sp(saneaixy->p,nsa);
		dcellfree(saneaxyi);
		dcellfree(saneaxyl);
		dcellfree(saneaixy);
	    }else{
		do_ref=1;
	    }
	}else{
	    /*compute nea from nearecon, scaled by area and dtrat. nea scales as sqrt(1/dtrat) */
	    if(iwfs==iwfs0  || parms->wfs[iwfs].sabad){
		const double neasq=pow(parms->powfs[ipowfs].nearecon/206265000.,2)/parms->powfs[ipowfs].dtrat;
		if(neasq<1.e-30) error("nea is too small\n");
		dmat *nea=dnew(nsa,2);
		/*scale neasq by area^-1 if seeing limited */
		/*scale neasq by area^-2 if diffraction limited */
		/*only implementing seeing limited here. */
		PDMAT(nea, neap);
		double *area=powfs[ipowfs].saa->p;
		for(int isa=0; isa<nsa; isa++){
		    double tmp=neasq/area[isa];
		    if(samask && samask->p[isa]){
			warning("wfs %d sa %d is masked\n", iwfs, isa);
			tmp=INFINITY;
		    }
		    neap[0][isa]=neap[1][isa]=tmp;
		}
		sanea->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
		dcwpow(nea, -1);
		saneai->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
		dcwpow(nea, -0.5);
		saneal->p[iwfs+iwfs*nwfs]=dspnewdiag(nsa*2,nea->p,1.);
		dfree(nea);
	    }else{
		do_ref=1;
	    }
	}
	

	if(do_ref){
	    sanea->p[iwfs+iwfs*nwfs] =dspref( sanea->p[iwfs0+iwfs0*nwfs]);
	    saneal->p[iwfs+iwfs*nwfs]=dspref(saneal->p[iwfs0+iwfs0*nwfs]);
	    saneai->p[iwfs+iwfs*nwfs]=dspref(saneai->p[iwfs0+iwfs0*nwfs]);
	}else{
	    dspscale(recon->saneai->p[iwfs+iwfs*nwfs], TOMOSCALE);
	}
	lfree(samask);
    }/*iwfs */
    
    /*Compute the averaged SANEA for each WFS */
    dfree(recon->neam);
    recon->neam=dnew(parms->nwfsr, 1);
    double neamhi=0; 
    int counthi=0;
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip==3) continue;
	const int nsa=powfs[ipowfs].saloc->nloc;
	dmat *sanea_iwfs=dspdiag(recon->sanea->p[iwfs+iwfs*parms->nwfsr]);
	double area_thres;
	if(nsa>4){
	    area_thres=0.9*parms->powfs[ipowfs].safill2d;
	}else{
	    area_thres=0;
	}
	double nea2_sum=0;
	int count=0;
	for(int isa=0; isa<nsa; isa++){
	    if(powfs[ipowfs].saa->p[isa]>area_thres){
		nea2_sum+=(sanea_iwfs->p[isa])+(sanea_iwfs->p[isa+nsa]);
		count++;
	    }
	}
	dfree(sanea_iwfs);
	recon->neam->p[iwfs]=sqrt(nea2_sum/count/2);/*average sanea in radian */
	double pixtheta=parms->powfs[ipowfs].pixtheta;
	if(recon->neam->p[iwfs]>pixtheta*4 
	   && parms->powfs[ipowfs].usephy
	   && parms->powfs[ipowfs].order==1
	    ){
	    warning("WFS %d has too much measurement error. Ignore it\n", iwfs);
	    //Neglecting WFS whos NEA is greater than twice pixel size in
	    //physical optics mode.
	    dspfree(recon->saneai->p[iwfs+iwfs*parms->nwfsr]);
	    recon->saneai->p[iwfs+iwfs*parms->nwfsr]=dspnewdiag(nsa*2,NULL,0);
	    dspfree(recon->saneal->p[iwfs+iwfs*parms->nwfsr]);
	    recon->saneal->p[iwfs+iwfs*parms->nwfsr]=dspnewdiag(nsa*2,NULL,0);
	    dspfree(recon->sanea->p[iwfs+iwfs*parms->nwfsr]);
	    recon->sanea->p[iwfs+iwfs*parms->nwfsr]=dspnewdiag(nsa*2,NULL, pixtheta*1e4);
	    recon->neam->p[iwfs]=INFINITY;
	}
	char *neatype;
	if(parms->powfs[ipowfs].neareconfile){
	    neatype="FILE";
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    neatype="mtch";
	}else{
	    neatype="geom";
	}
	info2("%s(%.2f) ", neatype, recon->neam->p[iwfs]*206265000);
    }
    dscale(recon->neam, 1./sqrt(TOMOSCALE));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip==3) continue;
	if(!parms->powfs[ipowfs].lo){
	    neamhi+=pow(recon->neam->p[iwfs],2);
	    counthi++;
	}
    }
    recon->neamhi=sqrt(neamhi/counthi);
    info2("\n");
    if(parms->save.setup){
	writebin(recon->sanea, "sanea");
	writebin(recon->saneal,"saneal");
	writebin(recon->saneai,"saneai");
    }
}

/**
   wrapps setup_recon_TTR() and setup_recon_DFR() to removal global tip/tilt and
   differential focus.
*/
static void
setup_recon_TTFR(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    dcellfree(recon->PTT);
    dcellfree(recon->PDF);
    cellfree(recon->PTTF);

    recon->PTT=dcellpinv(recon->TT,recon->saneai);
    recon->PDF=dcellpinv(recon->DF,recon->saneai);
    recon->PTTF=dcellpinv(recon->TTF,recon->saneai);
    if(parms->save.setup){
	writebin(recon->TTF, "TTF");	
	writebin(recon->PTT, "PTT");
	writebin(recon->PTTF, "PTTF");
    }
    /*dcellfree(recon->DF);//don't free DF to use in PDF. */
    /*Keep TT, PTT, used in fsm pointing or dithering. */
}
/**
   Frees recon->invpsd or recon->fractal
*/
static void free_cxx(RECON_T *recon){
    if(recon->invpsd){
	dcellfree(recon->invpsd->invpsd);
	ccellfree(recon->invpsd->fftxopd);
	free(recon->invpsd);
	recon->invpsd=NULL;
    }
    if(recon->fractal){
	dcellfree(recon->fractal->xopd);
	free(recon->fractal);
	recon->fractal=NULL;
    }
}
/**
   Prepares for tomography. ALlow it to be called multiple times for Cn2 update.
*/
void
setup_recon_tomo_prep(RECON_T *recon, const PARMS_T *parms){
    /*Free existing struct if already exist.  */
    free_cxx(recon);
    if(parms->tomo.assemble){
	/*We need the old copy of L2 when we update the turbulence profile. */
	dspcellfree(recon->L2save);
	recon->L2save=recon->L2;
    }else{
	dspcellfree(recon->L2);
    }
    recon->L2=NULL;
    /*When layers get a weight less than 1%, we put it at 1% to avoid
      regularization unstability issues.*/
    dclip(recon->wt, 0.01, 1);
    /*normalize the weights to sum to 1. */
    normalize_sum(recon->wt->p, recon->npsr, 1);
    const int npsr=recon->npsr;
    recon->cxx=parms->tomo.cxx;
    /*test_cxx(recon, parms); */
    if(parms->tomo.cxx==0){
	if(parms->load.cxx){
	    recon->L2=dspcellread("%s",parms->load.cxx);
	    if(recon->L2->nx!=npsr || recon->L2->ny!=npsr){
		error("Wrong format of loaded L2\n");
	    }
	}else{
	    recon->L2=cellnew(npsr,npsr);
	    for(int ips=0; ips<npsr; ips++){
		if(parms->tomo.square){/*periodic bc */
		    recon->L2->p[ips+npsr*ips]=mklaplacian_map
			(recon->xmap->p[ips]->nx, recon->xmap->p[ips]->nx,
			 recon->xloc->p[ips]->dx, recon->r0,
			 recon->wt->p[ips]);
		}else{/*reflecive bc */
		    recon->L2->p[ips+npsr*ips]=mklaplacian_loc
			(recon->xloc->p[ips], recon->r0, 
			 recon->wt->p[ips]);
		}
	    }
	}
	if(parms->save.setup){
	    writebin(recon->L2, "L2");
	}
	dspcellscale(recon->L2, sqrt(parms->tomo.cxxscale*TOMOSCALE));
    }
    if(parms->tomo.cxx==1 || (parms->tomo.cxx==2 && parms->tomo.precond==1)){
	recon->invpsd=calloc(1, sizeof(INVPSD_T));
	if(parms->load.cxx){
	    recon->invpsd->invpsd=dcellread("%s",parms->load.cxx);
	    if(recon->invpsd->invpsd->nx!=npsr || recon->invpsd->invpsd->ny!=1){
		error("Wrong format of loaded invpsd\n");
	    }
	}else{
	    dcell* invpsd=recon->invpsd->invpsd=cellnew(npsr,1);
	    for(int ips=0; ips<npsr; ips++){
		long nx=recon->xmap->p[ips]->nx;
		long ny=recon->xmap->p[ips]->ny;
		double r0i=recon->r0*pow(recon->wt->p[ips],-3./5.);
		invpsd->p[ips]=turbpsd(nx, ny, recon->xloc->p[ips]->dx, r0i,
				       recon->L0, 0, -1);
		dscale(invpsd->p[ips], pow((double)(nx*ny),-2));
	    }
	}
	if(parms->save.setup){
	    writebin(recon->invpsd->invpsd, "invpsd");
	}
	dcellscale(recon->invpsd->invpsd, sqrt(parms->tomo.cxxscale*TOMOSCALE));
	
	ccell* fftxopd=recon->invpsd->fftxopd=cellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    fftxopd->p[ips]=cnew(recon->xmap->p[ips]->nx, recon->xmap->p[ips]->ny);
	    //cfft2plan(fftxopd->p[ips],-1);
	    //cfft2plan(fftxopd->p[ips],1);
	}
	recon->invpsd->xloc = recon->xloc;
	recon->invpsd->square = parms->tomo.square;
    }
    if(parms->tomo.cxx==2){
	recon->fractal=calloc(1, sizeof(FRACTAL_T));
	recon->fractal->xloc=recon->xloc;
	recon->fractal->r0=parms->atmr.r0;
	recon->fractal->L0=parms->atmr.L0;
	recon->fractal->wt=parms->atmr.wt->p;
	recon->fractal->scale=sqrt(parms->tomo.cxxscale*TOMOSCALE);
	recon->fractal->ninit=parms->tomo.ninit;
	dcell *xopd=recon->fractal->xopd=cellnew(npsr, 1);
	for(int ips=0; ips<npsr; ips++){
	    int nn=nextfftsize(MAX(recon->xmap->p[ips]->nx, recon->xmap->p[ips]->ny))+1;
	    xopd->p[ips]=dnew(nn,nn);
	}
    }
   
    if(parms->tomo.piston_cr){
	/*when add constraint, make sure the order of
	  magnitude are at the same range.*/
	dspcellfree(recon->ZZT);
	recon->ZZT=cellnew(npsr,npsr);
	for(int ips=0; ips<npsr; ips++){
	    double r0=recon->r0;
	    double dx=recon->xloc->p[ips]->dx;
	    double wt=recon->wt->p[ips];
	    double val=pow(laplacian_coef(r0,wt,dx),2)*1e-6;
	    /*info("Scaling of ZZT is %g\n",val); */
	    /*piston mode eq 47 in Brent 2002 paper */
	    int icenter=loccenter(recon->xloc->p[ips]);
	    int nloc=recon->xloc->p[ips]->nloc;
	    dsp *ZZT=recon->ZZT->p[ips+npsr*ips]
		=dspnew(nloc,nloc,1);
	    int icol;
	    int count=0;
	    for(icol=0; icol<nloc; icol++){
		ZZT->p[icol]=count;
		if(icol==icenter){
		    ZZT->i[count]=icenter;
		    ZZT->x[count]=val;
		    count++;
		}
	    }
	    ZZT->p[nloc]=count;
	}
	if(parms->save.setup){
	    writebin(recon->ZZT, "ZZT");
	}
	dspcellscale(recon->ZZT, parms->tomo.cxxscale*TOMOSCALE);
    }
}
/**
   assemble tomography matrix. In CG mode, this function is not executed if
   tomo.assemble=0, Instead, the algorithm is contained in recon.c. When you
   modify anything, make sure you also do it there.  

   For integrated tomograhy: 
   
   \f$\hat{x}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1}+G_{ngs}^{T}C_{ngs}^{-1}
   G_{ngs})^{-1}(G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}+G_{ngs}^{T}C_{ngs}^{-1}s_{ngs}){\equiv}RL^{-1}RR s.\f$

   For split tomography, the terms regarding NGS are dropped:

   \f$\hat{x}_{lgs}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1})^{-1}
   G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}\equiv R_L^{-1} R_R s.\f$

   The left hand side of linear equation (inside the inverse) is stored in RL.
   The right hand side of the linear equation (outside of the inverse) is stored
   in RR. The terms regarding the NGS are handled using low rank terms. The
   gradients from LGS and NGS and concatenated to \f$s\f$.

   In the LGS part, there is a global tip/tilt removal operator because LGS is
   insensitive to global tip/tilts.

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803

*/
void setup_recon_tomo_matrix(RECON_T *recon, const PARMS_T *parms){
    /*if not cg or forced, build explicitly the tomography matrix. */
    int npsr=recon->npsr;
    int nwfs=parms->nwfsr;
    /*Free OLD matrices if any. */
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    info2("Before assembling tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);

    if(parms->load.tomo){
	/*right hand side. */
	warning("Loading saved recon->RR\n");
	recon->RR.M=readbin("RRM");
	if(recon->has_ttr){
	    recon->RR.U=dcellread("RRU");
	    recon->RR.V=dcellread("RRV");
	}
	/*Left hand side */
	warning("Loading saved recon->RL\n");
	recon->RL.M=readbin("RLM");
	recon->RL.U=dcellread("RLU");
	recon->RL.V=dcellread("RLV");
	if(parms->tomo.alg==0 && zfexist("RLC")){
	    recon->RL.C=chol_read("RLC");
	}
	if(parms->tomo.alg==2 && zfexist("RLMI")){
	    recon->RL.MI=dread("RLMI");
	}
    }else{
	info2("Building recon->RR\n");
	PDSPCELL(recon->GX,GX);
	const dspcell *saneai=recon->saneai;
	/*
	  Reconstruction Right hand side matrix. In split tomography mode, low
	  order NGS are skipped. recon->GXtomo contains GXs that only
	  participate in tomography.
	*/
	dspcell *GXtomoT=dspcelltrans(recon->GXtomo);
	recon->RR.M=dcellmm2(GXtomoT, saneai, "nn");
	PDSPCELL(recon->RR.M, RRM);
	/*
	  Tip/tilt and diff focus removal low rand terms for LGS WFS.
	*/
	if(recon->TTF){
	    dcellmm(&recon->RR.U, recon->RR.M, recon->TTF, "nn", 1);
	    recon->RR.V=dcelltrans(recon->PTTF);
	}
 
	info2("Building recon->RL\n"); /*left hand side matrix */
	recon->RL.M=dcellmm2(recon->RR.M,recon->GXtomo, "nn");
	PDSPCELL(recon->RL.M,RLM);
	if(parms->tomo.piston_cr){ 
	    /*single point piston constraint. no need tikholnov.*/
	    info2("Adding ZZT to RLM\n");
	    for(int ips=0; ips<npsr; ips++){
		dspadd(&RLM[ips][ips], 1, recon->ZZT->p[ips+ips*npsr], 1);
	    }
	    dspcellfree(recon->ZZT);
	}
	/*Apply tikholnov regularization.*/
	if(fabs(parms->tomo.tikcr)>1.e-15){
	    /*Estimated from the Formula */
	    double maxeig=pow(recon->neamhi * recon->xloc->p[0]->dx, -2);
	    double tikcr=parms->tomo.tikcr;
	    info2("Adding tikhonov constraint of %g to RLM\n",tikcr);
	    info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	    dcelladdI(recon->RL.M, tikcr*maxeig);
	}
	/*add L2 and ZZT */
	switch(parms->tomo.cxx){
	case 0:/*Add L2'*L2 to RL.M */
	    for(int ips=0; ips<npsr; ips++){
		dsp* tmp=dspmulsp(recon->L2->p[ips+npsr*ips], recon->L2->p[ips+npsr*ips],"tn");
		if(!tmp){
		    error("L2 is empty!!\n");
		}
		dspadd(&RLM[ips][ips], 1, tmp, 1);
		dspfree(tmp);
	    }
	    break;
	case 1:/*Need to apply invpsd separately */
	    recon->RL.extra = recon->invpsd;
	    recon->RL.exfun = apply_invpsd;
	    break;
	case 2:/*Need to apply fractal separately */
	    recon->RL.extra = recon->fractal;
	    recon->RL.exfun = apply_fractal;
	}

	/*Symmetricize, remove values below 1e-15*max and sort RLM (optional). */
	/*dspcellsym(recon->RL.M); */

	/*Low rank terms for low order wfs. Only in Integrated tomography. */
	dcell *ULo=cellnew(npsr,nwfs);
	PDCELL(ULo, pULo);
	dcell *VLo=cellnew(npsr,nwfs);
	PDCELL(VLo, pVLo);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip){
		continue;
	    }
	    if(parms->powfs[ipowfs].lo){
		for(int ips=0; ips<npsr; ips++){
		    dspfull(&pULo[iwfs][ips], RRM[iwfs][ips],'n',-1);
		    dspfull(&pVLo[iwfs][ips], GX[ips][iwfs],'t',1);
		}
	    }
	}
	if(!parms->recon.split || parms->tomo.splitlrt || parms->recon.split==2){
	    recon->RL.U=dcellcat(recon->RR.U, ULo, 2);
	    dcell *GPTTDF=NULL;
	    dcellmm(&GPTTDF, recon->GX, recon->RR.V, "tn", 1);
	    recon->RL.V=dcellcat(GPTTDF, VLo, 2);
	    dcellfree(GPTTDF);
	}else{
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	}
	dcellfree(ULo);
	dcellfree(VLo);
	/*Remove empty cells. */
	//dcelldropempty(&recon->RR.U,2);
	//dcelldropempty(&recon->RR.V,2);
	dcelldropempty(&recon->RL.U,2);
	dcelldropempty(&recon->RL.V,2);

	if(recon->RL.U){
	    /* balance UV. may not be necessary. Just to compare well against
	       laos. */
	    double r0=recon->r0;
	    double dx=recon->xloc->p[0]->dx;
	    double val=laplacian_coef(r0,1,dx);/*needs to be a constant */
	    dcellscale(recon->RL.U, 1./val);
	    dcellscale(recon->RL.V, val);
	}
	/*collect statistics.*/
	long nll=0,nlr=0;
	if(recon->RR.U){
	    int nx=recon->RR.U->nx;
	    for(int i=0; i<recon->RR.U->ny;i++){
		if(recon->RR.U->p[i*nx]){
		    nlr+=recon->RR.U->p[i*nx]->ny;
		}
	    }
	}
	if(recon->RL.U){
	    int nx=recon->RL.U->nx;
	    for(int i=0; i<recon->RL.U->ny;i++){
		if(recon->RL.U->p[i*nx]){
		    nll+=recon->RL.U->p[i*nx]->ny;
		}
	    }
	}
	info2("Tomography number of Low rank terms: %ld in RHS, %ld in LHS\n", nlr,nll);
	if(parms->save.recon){
	    writebin(recon->RR.M,"RRM");
	    writebin(recon->RR.U,"RRU");
	    writebin(recon->RR.V,"RRV");

	    writebin(recon->RL.M,"RLM.bin");/*disable compression */
	    writebin(recon->RL.U,"RLU");
	    writebin(recon->RL.V,"RLV"); 
	}
	dspcellfree(GXtomoT);
    }
    if((parms->tomo.alg==0 || parms->tomo.alg==2) && parms->tomo.bgs ){
	/* We need cholesky decomposition in CBS or MVST method. */
	muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
    }
    if(((parms->tomo.alg==0 || parms->tomo.alg==2) && !parms->tomo.bgs) || parms->sim.ecnn){
	if(parms->load.tomo){
	    if(parms->tomo.alg==0 && zfexist("RLC")){
		recon->RL.C=chol_read("RLC");
	    }
	    if(parms->tomo.alg==2 && zfexist("RLMI")){
		recon->RL.MI=dread("RLMI");
	    }
	}
	if(!recon->RL.C && !recon->RL.MI){
	    muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }

    if(parms->save.recon){
       	if(recon->RL.C){
	    chol_convert(recon->RL.C, 1);
	    chol_save(recon->RL.C,"RLC.bin");
	}
	if(recon->RL.MI){
	    writebin(recon->RL.MI,"RLMI");
	}
	if(recon->RL.Up){
	    writebin(recon->RL.Up,"RLUp"); 
	    writebin(recon->RL.Vp,"RLVp"); 
	}
	if(recon->RL.CB){
	    for(int ib=0; ib<recon->RL.nb; ib++){
		chol_save(recon->RL.CB[ib],"RLCB_%d.bin", ib);
	    }
	}
	if(recon->RL.MIB){
	    writebin(recon->RL.MIB,"RLMIB");
	}
    }
    /*Don't free PTT. Used in forming LGS uplink err */
    info2("After assemble tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}

static dcell * setup_recon_ecnn(RECON_T *recon, const PARMS_T *parms, loc_t *locs, lmat *mask){
    /**
       We compute the wavefront estimation error covariance in science focal
       plane due to wavefront measurement noise. Basically we compute
       Hx*E*Cnn*E'*Hx' where E is the tomography operator, and Hx is ray
       tracing from tomography grid xloc to science focal plane ploc. Since
       Cnn is symmetrical and sparse, we can decompose it easily into
       Cnn=Cnl*Cnl'; We first compute L=Hx*E*Cnl, and the result is simply
       LL'; This is much faster than computing left and right separately,
       because 1) the number of points in xloc is larger than in Cnn, so
       after the tomography right hand side vector is applied, the number of
       rows is larger than number of columns, this causes the right hand
       side solver to be much slower. 2) Simply double the computation.

       For HX opeation, build the sparse matrix and do multiply is way
       slower than doing ray tracing directly. 

       For ad hoc split tomography, we need to remove the five NGS modes
       from here, as well as in time averaging of estimated turbulence.

       recon->saneal contains Cnl.
    */
    TIC;tic;
    read_self_cpu();
    dmat *t1=NULL;
    if(recon->MVM){//MVM
	dspcell *sanealhi=cellnew(parms->nwfsr, parms->nwfsr);
	for(int iwfsr=0; iwfsr<parms->nwfsr; iwfsr++){
	    int ipowfs=parms->wfsr[iwfsr].powfs;
	    if(!parms->powfs[ipowfs].skip){
		IND(sanealhi, iwfsr, iwfsr)=dspref(IND(recon->saneal, iwfsr, iwfsr));
	    }
	}
	//dsp *tmp=dspcell2sp(sanealhi);  dspcellfree(sanealhi);
	dmat *tmp=dcell2m(sanealhi); dspcellfree(sanealhi);
	dcellmm(&t1, recon->MVM, tmp, "nn", 1);
	cellfree(tmp);
	toc2("MVM ");tic;
    }else if(parms->recon.alg==0){//MV
	dcell *tmp2=NULL;
	dcell *tmp=NULL;
	dspcellfull(&tmp2, recon->saneal, 'n', 1); 
	muv(&tmp, &recon->RR, tmp2, 1); 
	toc2("RR ");tic;
	cellfree(tmp2);
	muv_direct_solve(&tmp2, &recon->RL, tmp); dcellfree(tmp);
	toc2("RL ");tic;
	//2015-06-29: Put in ommited DM fitting operation
	muv(&tmp, &recon->FR, tmp2, 1); dcellfree(tmp2);
	toc2("FR ");tic;
	muv_direct_solve(&tmp2, &recon->FL, tmp); dcellfree(tmp);
	toc2("FL ");tic;
	t1=dcell2m(tmp2); dcellfree(tmp2);
    }else{//LSR
	dcell *tmp2=NULL;
	dcell *tmp=NULL;
	dspcellfull(&tmp2, recon->saneal, 'n', 1); 
	muv(&tmp, &recon->LR, tmp2, 1); cellfree(tmp2);
	toc2("LR ");tic;
	muv_direct_solve(&tmp2, &recon->LL, tmp); dcellfree(tmp);
	toc2("LL ");tic;
	t1=dcell2m(tmp2); dcellfree(tmp2);
    }
    dcell *ecnn=cellnew(parms->evl.nevl, 1);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(mask && !mask->p[ievl]) continue;
	tic;
	/*Build HX for science directions that need ecov.*/
	dmat *x1=dnew(locs->nloc, t1->ny);
	double hs=parms->evl.hs->p[ievl];
	int offset=0;
	for(int idm=0; idm<parms->ndm; idm++){
	    const double ht=parms->dm[idm].ht;
	    const double scale=1.-ht/hs;
	    const double dispx=parms->evl.thetax->p[ievl]*ht;
	    const double dispy=parms->evl.thetay->p[ievl]*ht;
	    for(int icol=0; icol<t1->ny; icol++){
		prop_nongrid(recon->aloc->p[idm], PCOL(t1, icol)+offset,
			     locs, PCOL(x1, icol), 1, dispx, dispy, scale, 0, 0);
	    }
	    offset+=recon->aloc->p[idm]->nloc;
	}
	toc2("Prop ");tic;
	dmm(&ecnn->p[ievl], 0, x1, x1, "nt", 1);
	dfree(x1);
	toc2("MM ");
    }
    dcellfree(t1);
    return ecnn;
}
/**
   Update assembled tomography matrix with new L2. Called from cn2est when new
   profiles are available.
*/
void setup_recon_tomo_update(RECON_T *recon, const PARMS_T *parms){
    setup_recon_tomo_prep(recon, parms); /*redo L2, invpsd */
#if USE_CUDA
    if(parms->gpu.tomo){
	gpu_update_recon_cn2(parms, recon);
    }
#endif
    if(parms->tomo.alg==1&&!parms->tomo.assemble){/*no need to do anything */
	return;
    }
    if(parms->tomo.cxx==0 && recon->L2save){
	/*Need to adjust RLM with the new L2. */
	PDSPCELL(recon->RL.M,RLM);
	const int npsr=recon->npsr;
	for(int ips=0; ips<npsr; ips++){
	    dsp* LL=dspmulsp(recon->L2->p[ips+npsr*ips], 
			     recon->L2->p[ips+npsr*ips],"tn");
	    dsp *LLold=dspmulsp(recon->L2save->p[ips+npsr*ips], 
				recon->L2save->p[ips+npsr*ips],"tn");
	    if(!LL){
		error("L2 is empty!!\n");
	    }
	    dsp *LLdiff=dspadd2(LL,1,LLold,-1);/*adjustment to RLM */
	    dspadd(&RLM[ips][ips], 1, LLdiff, 1);
	    dspfree(LLdiff);
	    dspfree(LL);
	    dspfree(LLold);
	}
    }

    if(parms->tomo.alg==0 || parms->tomo.alg==2){
	/*We need cholesky decomposition in CBS or MVST method. */
	if(!parms->tomo.bgs){/*Full Matrix */
	    muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}else{/*BGS */
	    muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }

    if(parms->recon.split==2){
	setup_recon_mvst(recon,parms);
    }
}


/**
   Create the reconstructor to reconstruct the residual focus error due to LGS
   sodium tracking error. Need to average out the focus error caused by
   atmosphere when applying (a low pass filter is applied to the output).  */
static void
setup_recon_focus(RECON_T *recon, POWFS_T *powfs, const PARMS_T *parms){
    if(!parms->nlgspowfs){
	return;
    }
    
    if(parms->recon.split==2 && parms->sim.mffocus){//For MVST.
	dmat *GMGngs=NULL;
	dcell *GMngs=cellnew(1, parms->nwfsr);
	/*Compute focus reconstructor from NGS Grads. fuse grads
	  together to construct a single focus measurement*/
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].trs==0 && parms->powfs[ipowfs].order>1 && parms->powfs[ipowfs].skip!=2){
		info2("wfs %d will be used to track focus\n", iwfs);
	    }else{
		continue;
	    }
	    dspmm(&GMngs->p[iwfs], recon->saneai->p[iwfs+parms->nwfsr*iwfs], 
		  recon->GFall->p[ipowfs],"nn",1);
	    dmm(&GMGngs,1,recon->GFall->p[ipowfs], GMngs->p[iwfs], "tn",1);
	}
	dinvspd_inplace(GMGngs);
	/*A focus reconstructor from all NGS measurements.*/
	dcell *RFngsg=recon->RFngsg=cellnew(1, parms->nwfsr);
  
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(!recon->GFall->p[ipowfs]) continue;
	    //NGS gradient to Focus mode reconstructor.
	    dmm(&RFngsg->p[iwfs], 0, GMGngs, GMngs->p[iwfs],"nt",1);
	}
	dfree(GMGngs);
	dcellfree(GMngs);
    }
    /*
      Compute focus constructor from LGS grads. A constructor for each LGS
      because each LGS may have different range error. Applyes to parms->wfs, not parms->wfsr.
    */
    cellfree(recon->RFlgsg);
    recon->RFlgsg=cellnew(parms->nwfs, parms->nwfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].llt) continue;
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    int iwfs0=parms->powfs[ipowfs].wfs->p[0];
	    if(iwfs==iwfs0 || !parms->recon.glao){
		IND(recon->RFlgsg, iwfs, iwfs)=dpinv(recon->GFall->p[ipowfs], IND(recon->saneai, iwfs, iwfs));
	    }else{
		IND(recon->RFlgsg, iwfs, iwfs)=dref(IND(recon->RFlgsg, iwfs0, iwfs0));
	    }
	}
    }
     
    if(parms->save.setup){
	writebin(recon->RFngsg,"RFngsg");
	writebin(recon->RFlgsg,"RFlgsg");
    }
}

/**
   Setup reconstructor for TWFS
*/
static void
setup_recon_twfs(RECON_T *recon, POWFS_T *powfs, const PARMS_T *parms){
    cellfree(recon->RRtwfs);
    dcell *GRtwfs=cellnew(parms->nwfsr, 1);
    dspcell *neai=cellnew(parms->nwfsr, parms->nwfsr);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip==2){//twfs
	    GRtwfs->p[iwfs]=dref(recon->GRall->p[ipowfs]);
	    neai->p[iwfs+iwfs*parms->nwfsr]=dspref(recon->saneai->p[iwfs+parms->nwfsr*iwfs]);
	}
    }
    recon->RRtwfs=dcellpinv(GRtwfs, neai);
    if(parms->save.setup){
	writebin(recon->GRall, "GRall");
	writebin(recon->RRtwfs, "RRtwfs");
    }
    /*if(parms->recon.fnsphpsd){
	recon->eptwfs=twfs_gain_optim(parms, recon, powfs);
	warning("eptwfs is reset to %g\n", parms->sim.eptwfs);
	}else{*/
    recon->eptwfs=parms->sim.eptwfs;
//}
    cellfree(GRtwfs);
    cellfree(neai);
}

/**
   compute the MVST split tomography NGS mode reconstructor.

   2010-03-16:
   New Implementation.

   Definition:

   - \f$G_{lgs}\f$ is LGS gradient operator from xloc, contained in recon->GXtomo as dspcell
   - \f$G_{ngs}\f$ is NGS gradient operator from xloc, contained in recon->GXL as dense matrix
   - \f$C_{lgs}, C_{ngs}\f$ is the LGS, and NGS measurement noise covariance matrix.
   - \f$\hat{x}_{lgs}\f$ is LGS tomography output
   - \f$\hat{x}_{ngs}\f$ is NGS minimum variance split tomography output
   - \f$a_{ngs}\f$ is NGS minimum variance DM output.
   - \f$F\f$ is the fitting operator
   - \f$H_A\f$ is the ray tracing operator from aloc to ploc, contained in recon->HA.
   - \f$MVModes\f$ is the MVST NGS modes
   - \f$MVRngs\f$ is the MVST NGS reconstructor

   We have
 
   \f{eqnarray*}{
   \hat{x}_{lgs}&=&A^{-1}G^{T}_{lgs}C_{lgs}^{-1}s_{lgs}\\
   Uw&=&A^{-1}G_{ngs}^TC_{ngs}^{-1}\\
   \hat{x}_{ngs}&=&Uw(1+G_{ngs}Uw)^{-1}(s_{ngs}-G_{ngs}\hat{x}_{lgs});\\
   a_{NGS}&=&F\hat{x}_{ngs}\\
   MVModes&=&F\cdot Uw\\
   MVRngs&=&(1+G_{ngs}Uw)^{-1}
   \f}

   If we want to orthnormalize the NGS modes. Propagate the NGS modes \f$F\cdot Uw\f$ to
   fitting directions and compute the cross-coupling matrix with weighting on \f$W_0, W_1\f$, we have

   \f{eqnarray*}{
   Q&=&H_A\cdot F \cdot Uw\\
   M_{CC}&=&Q^T(W_0-W_1 W_1^T)Q\\
   \f}

   Do SVD on \f$MCC\f$ we have \f$M_{CC}=u\Sigma v^T \f$ where \f$u{\equiv}v\f$
   because \f$MCC\f$ is symmetric. Redefine the NGS modes and reconstructor as

   \f{eqnarray*}{
   MVModes&=&F\cdot Uw\cdot u \Sigma^{-1/2}\\
   MVRngs&=&\Sigma^{1/2}u^T (1+G_{ngs}A^{-1}G_{ngs}^T C_{ngs}^{-1})^{-1}
   \f}
*/

void
setup_recon_mvst(RECON_T *recon, const PARMS_T *parms){
    TIC;tic;
    /*
      Notice that: Solve Fitting on Uw and using FUw to form Rngs gives
      slightly different answer than solve fitting after assemble the
      reconstructor. 10^-6 relative difference.
      
      2010-03-10: Bug found and fixed: The MVST with CBS-CBS method gives worst
      performance than integrated tomography. The probelm is in PUm
      computing. I mistakenly called chol_solve, while I should have
      called muv_direct_solve. The former doesn ot apply the low rank
      terms.
    */
    if(parms->recon.split!=2){
	return;
    }
    cellfree(recon->MVRngs);
    cellfree(recon->MVModes);
    cellfree(recon->MVGM);
    cellfree(recon->MVFM);

    dcellfree(recon->GXL);
    dcelladd(&recon->GXL, 1, recon->GXlo, 1);
    //NEA of low order WFS.
    dcell *neailo=cellnew(parms->nwfsr, parms->nwfsr);
    dcell *nealo=cellnew(parms->nwfsr, parms->nwfsr);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].lo){
	    dspfull(&neailo->p[iwfs*(1+parms->nwfsr)],
		    recon->saneai->p[iwfs*(1+parms->nwfsr)], 'n', 1);
	    dspfull(&nealo->p[iwfs*(1+parms->nwfsr)],
		    recon->sanea->p[iwfs*(1+parms->nwfsr)], 'n', 1);
	}
    }
    /* 2012-03-21: Remove focus mode from GL and NEA so no focus is
       measured/estimated by mvst. It is then estimate separately*/
    if(parms->sim.mffocus){
	dcell *focus=NULL;
	dcellmm(&focus, recon->RFngsg, recon->GXL, "nn", 1);
	dcellmm(&recon->GXL, recon->GFngs, focus, "nn", -1);
	//dcelldropzero(recon->GXL, 1e-12);
	dcellfree(focus);
	dcellmm(&focus, recon->RFngsg, neailo, "nn", 1);
	dcellmm(&neailo, recon->GFngs, focus, "nn", -1);
	//dcelldropzero(neailo, 1e-8);
	dcellfree(focus);
    }

    dcell *U=NULL; 
    dcell *FU=NULL;
    if(parms->load.mvst){
	U=dcellread("mvst_U");
	FU=dcellread("mvst_FU");
    }else{
	/*Prepare CBS if not already done */
	if(!recon->RL.C && !recon->RL.MI){
	    muv_direct_prep(&(recon->RL), 0);
	}
	if(!recon->FL.C && !recon->FL.MI){
	    muv_direct_prep(&(recon->FL), 0);
	}
	toc2("MVST: svd prep");
	dcell *GXLT=dcelltrans(recon->GXL);
	muv_direct_solve(&U, &recon->RL, GXLT);
	dcellfree(GXLT);
	dcell *rhs=NULL;
	muv(&rhs, &recon->FR, U, 1);
	muv_direct_solve(&FU, &recon->FL, rhs);
	dcellfree(rhs);
	toc2("MVST: U, FU");
    
	if(parms->save.mvst || parms->save.setup){
	    writebin(U, "mvst_U");
	    writebin(FU, "mvst_FU");
	}
    }
    dcell *Uw=NULL;
    dcell *FUw=NULL;

    dcellmm(&Uw, U, neailo, "nn", 1);
    dcellmm(&FUw, FU, neailo, "nn", 1);

    dcell *M=NULL;
    dcellmm(&M, recon->GXL, Uw, "nn", 1);
    dcelladdI(M, 1);
    dcell *Minv=dcellinv(M);
    dcellfree(M);
    if(parms->sim.mffocus){
	//Embed a focus removal. Necessary!
	dcell *focus=NULL;
	dcellmm(&focus, Minv, recon->GFngs, "nn", 1);
	dcellmm(&Minv, focus, recon->RFngsg, "nn", -1);
    }
    if(parms->save.setup){
	writebin(Minv, "mvst_Rngs_0");
	writebin(FUw,  "mvst_Modes_0");
    }
    /*Compute the cross coupling matrix of the Modes:
      FUw'*Ha'*W*Fuw*Ha. Re-verified on 2013-03-24.*/
    dcell *QwQc=NULL;
    {
	dcell *Q=NULL;/*the NGS modes in ploc. */
	dcellmm(&Q, recon->HA, FUw, "nn", 1);
	QwQc=calcWmcc(Q,Q,recon->W0,recon->W1,recon->fitwt);
	dcellfree(Q);
    }
    /*Compute the wavefront error due to measurement noise. Verified on
      2013-03-24. The gain optimization is yet a temporary hack because the way
      PSDs are input.*/
    if(0){
	dcell *RC=NULL;
	dcellmm(&RC, Minv, nealo, "nn", 1);
	dcell *RCRt=NULL;
	dcellmm(&RCRt, RC, Minv, "nt", 1);
	dcell *RCRtQwQ=NULL;
	dcellmm(&RCRtQwQ, RCRt, QwQc, "nn", 1);
	dmat *tmp=dcell2m(RCRtQwQ);
	PDMAT(tmp, ptmp);
	double rss=0;
	for(int i=0; i<tmp->nx; i++){
	    rss+=ptmp[i][i];
	}
	dfree(tmp);
	dcellfree(RCRtQwQ);
	dcellfree(RCRt);
	dcellfree(RC);
	recon->sigmanlo=rss;
	info("rms=%g nm\n", sqrt(rss)*1e9);

	if(zfexist("../../psd_ngs.bin")){
	    warning("Temporary solution for testing\n");
	    dmat *psd_ngs=dread("../../psd_ngs.bin");
	    if(parms->sim.wspsd){//windshake
		//need to convert from rad to m2.
		dmat *psd_ws=dread("%s", parms->sim.wspsd);
		dmat *psd_ws_m=ddup(psd_ws); 
		dfree(psd_ws);
		dmat *psd_ws_y=dnew_ref(psd_ws_m->nx,1,psd_ws_m->p+psd_ws_m->nx);
		dscale(psd_ws_y, 4./parms->aper.d); dfree(psd_ws_y);
		add_psd2(&psd_ngs, psd_ws_m); dfree(psd_ws_m);
	    }
	    writebin(psd_ngs, "psd_ngs_servo");
	    dmat *rss2=dnew(1,1); rss2->p[0]=rss;
	    int dtrat=parms->powfs[parms->lopowfs->p[0]].dtrat;
	    dcell *res=servo_optim(psd_ngs, parms->sim.dt, 
				   dtrat, M_PI/4, rss2, 2); 
	    dfree(rss2);
	    info("dtrat=%d\n", dtrat);
	    info("g,a,T was %g,%g,%g\n", parms->sim.eplo->p[0], parms->sim.eplo->p[1], parms->sim.eplo->p[2]);
	    memcpy(parms->sim.eplo->p, res->p[0]->p, 3*sizeof(double));
	    info("g,a,T=%g,%g,%g\n", parms->sim.eplo->p[0], parms->sim.eplo->p[1], parms->sim.eplo->p[2]);
	    info("res=%g, resn=%g nm\n", sqrt(res->p[0]->p[3])*1e9, sqrt(res->p[0]->p[4])*1e9);
	    dcellfree(res); 
	    dfree(psd_ngs);
	}
    }
    if(1){/*Orthnormalize the Modes.*/
	/*
	  Change FUw*Minv -> FUw*(U*sigma^-1/2) * (U*sigma^1/2)'*Minv
	  columes of FUw*(U*sigma^-1/2) are the eigen vectors.

	  U, sigma is the eigen value decomposition of <FUw' HA' W HA FUw>
	*/

	dmat *QSdiag=NULL, *QU=NULL, *QVt=NULL;
	{
	    dmat *QwQ=dcell2m(QwQc);
	    dsvd(&QU, &QSdiag, &QVt, QwQ);
	    dfree(QwQ);
	}
	if(parms->save.setup) {
	    writebin(QSdiag,"mvst_QSdiag");
	}
	dcwpow_thres(QSdiag, -1./2., 1e-14);
	dmuldiag(QU,QSdiag);/*U*sigma^-1/2 */
	d2cell(&QwQc,QU,NULL);
	dcell *FUw_keep=FUw;FUw=NULL; 
	dcellmm(&FUw, FUw_keep, QwQc, "nn", 1);
	dcellfree(FUw_keep);
	dcwpow_thres(QSdiag,-2, 1e-14);
	dmuldiag(QU,QSdiag);/*U*sigma^1/2 (From U*sigma^(-1/2)*sigma) */
	d2cell(&QwQc,QU,NULL);
	dcell *Minv_keep=Minv; Minv=NULL;
	dcellmm(&Minv,QwQc,Minv_keep,"tn", 1);
	dcellfree(Minv_keep);
	dfree(QSdiag);
	dfree(QVt);
	dfree(QU);
    }
    dcellfree(QwQc);
  
    recon->MVRngs=dcellreduce(Minv,1);/*1xnwfs cell */
    recon->MVModes=dcellreduce(FUw,2);/*ndmx1 cell */
    dcellmm(&recon->MVGM, recon->GAlo, recon->MVModes, "nn", 1);
    dcellmm(&recon->MVFM, recon->RFngsg, recon->MVGM, "nn", 1);
    dcellfree(neailo);
    dcellfree(nealo);
    dcellfree(Minv);
    dcellfree(U);
    dcellfree(FU);
    dcellfree(FUw);
    dcellfree(Uw);
    if(parms->save.setup){
	dcell *Qn=NULL;
	dcellmm(&Qn, recon->HA, recon->MVModes, "nn", 1);
	dcell *Qntt=cellnew(Qn->nx,Qn->ny);
	dmat *TTploc=loc2mat(recon->floc,1);/*TT mode. need piston mode too! */
	dmat *PTTploc=dpinv(TTploc,recon->W0);/*TT projector. no need w1 since we have piston. */
	dfree(TTploc);
	for(int ix=0; ix<Qn->nx*Qn->ny; ix++){
	    if(!Qn->p[ix]) continue;
	    dmm(&Qntt->p[ix], 0, PTTploc, Qn->p[ix],"nn",1);
	}
	writebin(Qntt,"mvst_modptt");
	dcellfree(Qn);
	dcellfree(Qntt);
	dfree(PTTploc);
	writebin(recon->MVRngs, "mvst_Rngs");
	writebin(recon->MVModes,"mvst_Modes");
    }
    if(parms->dbg.mvstlimit>0){/*limit number of modes used. */
	warning("MVST: Correction is limited to %d modes\n", parms->dbg.mvstlimit);
	dmat *tmp;
	for(int iy=0; iy<recon->MVRngs->ny; iy++){
	    tmp=recon->MVRngs->p[iy];
	    if(tmp){
		recon->MVRngs->p[iy]=dsub(tmp,0,parms->dbg.mvstlimit,0,tmp->ny);
		dfree(tmp);
	    }
	}
	for(int ix=0; ix<recon->MVModes->nx; ix++){
	    tmp=recon->MVModes->p[ix];
	    if(tmp){
		recon->MVModes->p[ix]=dsub(tmp,0,tmp->nx,0,parms->dbg.mvstlimit);
		dfree(tmp);
	    }
	}
	if(parms->save.setup){
	    writebin(recon->MVRngs, "mvst_Rngs_limit");
	    writebin(recon->MVModes,"mvst_Modes_limit");
	}
    }
    /*
    if(parms->save.setup){
	dcell *QQ=NULL;
	dcellmm(&QQ, recon->HA, recon->MVModes,"nn", 1);
	dcell *MCC=calcWmcc(QQ,QQ,recon->W0,recon->W1,recon->fitwt);
	writebin(MCC,"mvst_MCC");
    
	dcellfree(MCC);
	dcellfree(QQ);
	}*/
    toc2("MVST");
}

/**
   Sets up the tomogrpahy turbulence reconstruction structs including wavefront
   reconstructor and DM fitting operator \callgraph
   
   Calls setup_recon_tomo_matrix() 
   and setup_recon_fit_matrix() to setup the tomography and DM
   fitting matrix. 
 
   AHST is handled in setup_ngsmod().

   MVST is handled in setup_recon_mvst().
   
   MOAO is handled in setup_recon_moao().
   
*/
void setup_recon_tomo(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, const APER_T *aper){
    TIC;tic;
    /*setup inverse noise covariance matrix. */
    /*prepare for tomography setup */
    setup_recon_tomo_prep(recon,parms);
    if(parms->tomo.assemble || parms->recon.split==2 || parms->sim.psfr || parms->recon.psd){
	/*assemble the matrix only if not using CG CG apply the
	  individual matrices on fly to speed up and save memory. */
	setup_recon_tomo_matrix(recon,parms);
    }
	
    if(parms->tomo.precond==1){
	fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
	recon->fdpcg=fdpcg_prepare(parms, recon, powfs, NULL);
    }
	
    /*Fall back function method if .M is NULL */
    recon->RL.Mfun=TomoL;
    recon->RL.Mdata=recon;
    recon->RR.Mfun=TomoR;
    recon->RR.Mtfun=TomoRt;
    recon->RR.Mdata=recon;
    if(parms->tomo.alg==1){/*CG */
	switch(parms->tomo.precond){
	case 0:/*no preconditioner */
	    recon->RL.pfun=NULL;
	    recon->RL.pdata=NULL;
	    break;
	case 1:
	    recon->RL.pfun=fdpcg_precond;
	    recon->RL.pdata=(void*)recon;
	    break;
	default:
	    error("Invalid tomo.precond");
	}
    }
    recon->RL.alg = parms->tomo.alg;
    recon->RL.bgs = parms->tomo.bgs;
    recon->RL.warm  = parms->recon.warm_restart;
    recon->RL.maxit = parms->tomo.maxit;
    
    toc2("setup_recon_tomo");
}
/**
   Dither using command path (DM) aberration
 */
void setup_recon_dither_dm(RECON_T *recon, POWFS_T *powfs, const PARMS_T *parms){
    int any=0;
    int dither_mode=0;
    double dither_amp=0;
    int dither_npoint=0;
    int dither_dtrat=0;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].dither>1){//common path dithering
	    if(any){//already found, check consistency
		if(dither_mode!=-parms->powfs[ipowfs].dither ||
		   fabs(dither_amp-parms->powfs[ipowfs].dither_amp)>dither_amp*1e-5
		    || dither_npoint!=parms->powfs[ipowfs].dither_npoint
		    || dither_dtrat!=parms->powfs[ipowfs].dtrat){
		    error("Multiple dither with different configuration is not supported\n");
		}
	    }
	    dither_mode=-parms->powfs[ipowfs].dither;
	    dither_amp=parms->powfs[ipowfs].dither_amp;
	    dither_npoint=parms->powfs[ipowfs].dither_npoint;
	    dither_dtrat=parms->powfs[ipowfs].dtrat;
	    any=1;
	}
    }
    if(any){
	const int idm=parms->idmground;
	recon->dither_npoint=dither_npoint;
	recon->dither_dtrat=dither_dtrat;
	recon->dither_m=dcellnew(parms->ndm, 1);
	recon->dither_m->p[idm]=zernike(recon->aloc->p[idm], parms->aper.d, 0, 0, dither_mode);
	dscale(recon->dither_m->p[idm], dither_amp);
	recon->dither_ra=dcellnew(parms->ndm, parms->ndm);
	IND(recon->dither_ra, idm, idm)=dpinv(recon->dither_m->p[idm], 0);
	recon->dither_rg=dcellnew(parms->nwfsr, parms->nwfsr);
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].dither>1){
		dmat *opd=dnew(powfs[ipowfs].loc->nloc, 1);
		double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
		double scale=1.-ht/parms->powfs[ipowfs].hs;
		double dispx=ht*parms->wfsr[iwfs].thetax;
		double dispy=ht*parms->wfsr[iwfs].thetay;
		prop_nongrid(recon->aloc->p[idm], recon->dither_m->p[idm]->p, 
			     powfs[ipowfs].loc, opd->p, 
			     -1, dispx, dispy, scale, 0, 0);
		dmat *ints=0;
		dmat *grad=0;
		pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
		pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
		IND(recon->dither_rg, iwfs, iwfs)=dpinv(grad, IND(recon->saneai, iwfs, iwfs));
		if(0){//test linearity
		    dscale(opd, 1./4.);
		    dmat *tmp=0;
		    dmat *res=dnew(10,1);
		    for(int i=0; i<10; i++){
			dscale(opd, 2);
			dzero(ints);
			pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
			pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
			dmm(&tmp, 0, IND(recon->dither_rg, iwfs, iwfs), grad, "nn", 1);
			res->p[i]=tmp->p[0];
		    }
		    writebin(res, "linearity");
		    dfree(tmp);
		    dfree(res);
		    exit(0);
		}
		dfree(ints);
		dfree(opd);
		dfree(grad);
	    }
	}
	if(parms->save.setup){
	    writebin(recon->dither_m, "dither_m");
	    writebin(recon->dither_ra, "dither_ra");
	    writebin(recon->dither_rg, "dither_rg");
	}
    }
}

/**
   Setup either the minimum variance reconstructor by calling setup_recon_mvr()
   or least square reconstructor by calling setup_recon_lsr() */
void setup_recon(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, const APER_T *aper){
    TIC;tic;
    /*assemble noise equiva angle inverse from powfs information */
    setup_recon_saneai(recon,parms,powfs);
    /*setup LGS tip/tilt/diff focus removal */
    setup_recon_TTFR(recon,parms,powfs);
    if(parms->nlgspowfs){
	//if(parms->ilgspowfs!=-1 && (parms->sim.mffocus || parms->sim.ahstfocus || parms->dither)){
	/*mvst uses information here*/
	setup_recon_focus(recon, powfs, parms);
    }
    if(parms->itpowfs!=-1){ /*setup Truth wfs*/
	setup_recon_twfs(recon,powfs,parms);
    }
    if(!parms->sim.idealfit){
	if(parms->recon.mvm && parms->load.mvm){
	    recon->MVM=dread("%s", parms->load.mvm);
	}else{
	    switch(parms->recon.alg){
	    case 0:
		setup_recon_tomo(recon, parms, powfs, aper);
		break;
	    case 1:
		setup_recon_lsr(recon, parms, powfs);
		break;
	    default:
		error("recon.alg=%d is not recognized\n", parms->recon.alg);
	    }
	}
    }
    cellfree(recon->gloc);
    cellfree(recon->gamp);
    if(parms->recon.split){
	/*split tomography */
	setup_ngsmod_recon(parms,recon,aper,powfs);
	if(!parms->sim.idealfit && parms->recon.split==2 && parms->recon.alg==0){/*Need to be after fit */
	    setup_recon_mvst(recon, parms);
	}
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(!parms->powfs[ipowfs].needGS0 && powfs[ipowfs].GS0){
	    dspcellfree(powfs[ipowfs].GS0);
	    powfs[ipowfs].GS0=NULL;
	}
    }
    setup_recon_dither_dm(recon, powfs, parms);
    toc2("setup_recon");
}

/**
   PSD computation for gian update
 */
void setup_recon_psd(RECON_T *recon, const PARMS_T *parms){
    if(!parms->recon.psd) return;
    recon->Herr=cellnew(parms->evl.nevl, parms->ndm);
    double d1=parms->aper.d-parms->dm[0].dx;
    double d2=parms->aper.din+parms->dm[0].dx;
    double dx=(d1-d2)*0.1;
    if(dx<parms->dm[0].dx){
	dx=parms->dm[0].dx;
    }
    loc_t *eloc=mkannloc(d1, d2, dx, 0.8);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	for(int idm=0; idm<parms->ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dispx=parms->evl.thetax->p[ievl]*ht;
	    double dispy=parms->evl.thetay->p[ievl]*ht;
	    double scale=1-ht/parms->evl.hs->p[ievl];
	    IND(recon->Herr, ievl, idm)=mkh(recon->aloc->p[idm], eloc, dispx, dispy, scale);
	}
    }
    if(parms->recon.psd==2){//don't use signanhi by default
	dcell *ecnn=setup_recon_ecnn(recon, parms, eloc, 0);
	double sigma2e=0;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    double sigma2i=0;
	    for(int iloc=0; iloc<eloc->nloc; iloc++){
		sigma2i+=IND(ecnn->p[ievl], iloc, iloc);
	    }
	    sigma2e+=parms->evl.wt->p[ievl]*(sigma2i/eloc->nloc);
	}
	recon->sigmanhi=sigma2e;
	info2("High order WFS mean noise propagation is %g nm\n", sqrt(sigma2e)*1e9);
	if(parms->save.setup){
	    writebin(ecnn, "psd_ecnn");
	}
	dcellfree(ecnn);
    }
    if(parms->save.setup){
	writebin(recon->Herr, "Herr");
	writebin(eloc, "eloc.bin");
    }
    locfree(eloc);
}
	    
/**
   A few further operations that needs MVM.
 */
void setup_recon_post(RECON_T *recon, const PARMS_T *parms, const APER_T *aper){
    TIC;tic;
    if(parms->sim.ecnn){
	recon->ecnn=setup_recon_ecnn(recon, parms, aper->locs, parms->evl.psfr);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    char strht[24];
	    if(!isinf(parms->evl.hs->p[ievl])){
		snprintf(strht, 24, "_%g", parms->evl.hs->p[ievl]);
	    }else{
		strht[0]='\0';
	    }
	    writebin(recon->ecnn->p[ievl],  "ecnn_x%g_y%g%s.bin", 
		     parms->evl.thetax->p[ievl]*206265,
		     parms->evl.thetay->p[ievl]*206265, strht);
	}
    }
    if(parms->recon.psd){
	setup_recon_psd(recon, parms);
    }
    toc2("setup_recon_post");
}

/**
   Free unused object in recon struct after preparation is done.
 */
void free_recon_unused(const PARMS_T *parms, RECON_T *recon){
    if(!recon) return;
    /* Free arrays that will no longer be used after reconstruction setup is done. */
    dspcellfree(recon->sanea); 
    dspcellfree(recon->saneal);
    if(!(parms->tomo.assemble && parms->tomo.alg==1) && !parms->cn2.tomo && !parms->tomo.bgs){
	/*We no longer need RL.M,U,V */
	cellfree(recon->RL.M);
	dcellfree(recon->RL.U);
	dcellfree(recon->RL.V);
	cellfree(recon->RR.M);
	dcellfree(recon->RR.U);
	dcellfree(recon->RR.V);
    }else{
	dcellfree(recon->TTF);
	dcellfree(recon->PTTF);
    }
    if(parms->fit.alg!=1 && !parms->fit.bgs){
	cellfree(recon->FL.M);
	dcellfree(recon->FL.U);
	dcellfree(recon->FL.V);
    }
    if(parms->tomo.alg==1){
	muv_direct_free(&recon->RL);
    }
    if(parms->fit.alg==1){
	muv_direct_free(&recon->FL);
    }

    if(recon->RR.M){
	dspcellfree(recon->GP);
	dspcellfree(recon->GP2);
    }

    /*The following have been used in fit matrix. */
    if(parms->fit.assemble || parms->gpu.fit){
	dcellfree(recon->fitNW);
	dspcellfree(recon->actslave);
    }
    /* when sim.dmproj=1, we need these matrices to use in FR.Mfun*/
    if(recon->FR.M && !parms->sim.dmproj && parms->fit.assemble && !parms->gpu.moao && !parms->sim.ncpa_calib){
	dspfree(recon->W0); 
	dfree(recon->W1); 
	dspcellfree(recon->HA); 
	dspcellfree(recon->HXF); 
    }
    /*
      The following arrys are not used after preparation is done.
    */
    dspcellfree(recon->GX);
    dspcellfree(recon->GXtomo);/*we use HXWtomo instead. faster */
    if(!(parms->cn2.tomo && parms->recon.split==2)){/*mvst needs GXlo when updating. */
	dspcellfree(recon->GXlo);
    }
    if(parms->tomo.square && !parms->dbg.tomo_hxw && recon->RR.M){
	dspcellfree(recon->HXWtomo);
    }
    if(parms->recon.mvm){
	muv_free(&recon->RR);
	muv_free(&recon->RL);
	muv_free(&recon->FR);
	muv_free(&recon->FL);
	muv_free(&recon->LR);
	muv_free(&recon->LL);
	fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
	if(parms->gpu.tomo && parms->gpu.fit){
	    dfree(recon->MVM);//keep GPU copy.
	}
    }
}
/**
   Free the recon struct.
*/
void free_recon(const PARMS_T *parms, RECON_T *recon){
    if(!recon) return;
    ngsmod_free(recon->ngsmod); recon->ngsmod=0;
    free_recon_unused(parms, recon);
    free_recon_moao(recon, parms);
    dfree(recon->ht);
    dfree(recon->os);
    dfree(recon->wt);
    dfree(recon->dx);
    dcellfree(recon->MVRngs);
    dcellfree(recon->MVGM);
    dcellfree(recon->MVFM);
    dcellfree(recon->MVModes);
    dspcellfree(recon->GX);
    dspcellfree(recon->GXlo);
    dspcellfree(recon->GXtomo);
    dspcellfree(recon->GP);
    dspcellfree(recon->GP2);
    dspcellfree(recon->GA); 
    dspcellfree(recon->GAlo);
    dspcellfree(recon->GAhi);
    dcellfree(recon->GXL);
    dspcellfree(recon->L2);
    dspcellfree(recon->L2save);
    free_cxx(recon);
    dcellfree(recon->TT);
    dcellfree(recon->PTT);
    dcellfree(recon->DF);
    dcellfree(recon->PDF);
    dcellfree(recon->TTF);
    dcellfree(recon->PTTF);
    dcellfree(recon->GFngs);
    dcellfree(recon->GFall);
    dcellfree(recon->RFlgsg);
    dcellfree(recon->RFlgsa);
    dcellfree(recon->RFlgsx);
    dcellfree(recon->RFngsg);
    dcellfree(recon->RFngsa);
    dcellfree(recon->RFngsx);
    dcellfree(recon->GRall);
    dcellfree(recon->RRtwfs);
    dspcellfree(recon->ZZT);
    dspcellfree(recon->HXF); 
    dspcellfree(recon->HXW);
    dspcellfree(recon->HXWtomo);
    dspcellfree(recon->HA); 
    dspfree(recon->W0); 
    dfree(recon->W1); 
    dcellfree(recon->fitNW);
    dfree(recon->fitwt);

    cellfree(recon->xloc);
    cellfree(recon->xmap);
    cellfree(recon->xcmap);
    lfree(recon->xnx);
    lfree(recon->xny);
    lfree(recon->xnloc);
    dcellfree(recon->xmcc);

    lfree(recon->anx);
    lfree(recon->any);
    lfree(recon->anloc);
    lfree(recon->ngrad);
    locfree(recon->floc); 
    locfree(recon->ploc);
    cellfree(recon->ploc_tel);
    free(recon->amap->p);free(recon->amap);//data is referenced
    cellfree(recon->acmap);
    cellfree(recon->aloc);
    cellfree(recon->actstuck);
    cellfree(recon->actfloat);
    cellfree(recon->actinterp);
    cellfree(recon->actslave);
    cellfree(recon->actcpl);
    cellfree(recon->aimcc);/*used in filter.c */
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    muv_free(&recon->FR);
    muv_free(&recon->FL);
    muv_free(&recon->LR);
    muv_free(&recon->LL);
    dfree(recon->MVM);
    dspcellfree(recon->sanea);
    dspcellfree(recon->saneal);
    dspcellfree(recon->saneai);
    dfree(recon->neam); 
    fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
    cn2est_free(recon->cn2est);
    cellfree(recon->gloc);
    cellfree(recon->gamp);
    cellfree(recon->DMTT);
    cellfree(recon->DMPTT);
    free(recon);
}


/*
  some tips.
  1) UV balance need to be done carefully. The scaling must be a single number
  2) When add regularizations, the number has to be in the same order as the matrix.
  This is supper important!!!
  3) single point piston constraint in Tomography is not implemented correctly.

  2009-11-22
  Found a difference: Ha in LAOS computes the weighting even if not enough points.
  HA in MAOS only computes the weighting only if 4 coupled points all exist.
  AOS used laos geometry that are not enough to cover all points. Will modify mkH/accphi.
  H modified. Performance now agree with laos with AOS W0/W1 and HA

*/
