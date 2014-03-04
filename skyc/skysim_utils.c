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
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "skysim_utils.h"
/**
   \file skyc/skysim_utils.c
   Utilities for skysim.c
*/

/**
   Compute Open loop NGS mode wavefront error from mode vectors.  */
double calc_rms(const dmat *mod, const dmat *mcc, int istep0){
    double rms=0;
    PDMAT(mod, pmod);
    for(long istep=istep0; istep<mod->ny; istep++){
	rms+=dwdot(pmod[istep], mcc, pmod[istep]);
    }
    return rms/(mod->ny-istep0);
}

/**
   add photon and read out noise.  pcaclib part of bkgrnd is calibrated out.
   set to 1 usually.  */
static void addnoise(dmat *A, rand_t* rstat, 
		     const double bkgrnd, const double pcalib, 
		     const double rne){
    for(int ix=0; ix<A->nx*A->ny; ix++){
	A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)
	    -bkgrnd*pcalib+rne*randn(rstat);
    }
}

/**
   convert mod vector to ngs WFS cloc and add to it's opd or complex pupil function.
   Notice the the coordinate for each subaperture is different for TTF.
*/
void ngsmod2wvf(cmat *wvf,            /**<[in/out] complex pupil function*/
		double wvl,           /**<[in] the wavelength*/
		const dmat *modm,     /**<[in] the NGS mode vector*/
		const POWFS_S *powfs,       /**<[in] the powfs configuration*/
		int isa,              /**<[in] index of subaperture*/
		double thetax,        /**<[in] direction of WFS*/
		double thetay,        /**<[in] direction of WFS*/
		const PARMS_S *parms  /**<[in] the parms*/
    ){
    const double *mod=modm->p;
    const dcomplex ik=2*M_PI/wvl*I;
    double dx=powfs->dxwvf;
    double ox=powfs->saloc->locx[isa]+dx*0.5;
    double oy=powfs->saloc->locy[isa]+dx*0.5;
    int nx=powfs->nxwvf;
    assert(wvf->nx==wvf->ny);
    dcomplex *p=wvf->p+(wvf->nx-nx)/2*(1+wvf->nx);
    if(modm->nx==2){
	for(int iy=0; iy<nx; iy++){
	    double ym=(oy+iy*dx)*mod[1];
	    for(int ix=0; ix<nx; ix++){
		double x=ox+ix*dx;
		double tmp=x*mod[0]+ym;
		p[ix+wvf->nx*iy]*=cexp(ik*tmp);
	    }
	}
    }else{
	const double hc=parms->maos.hc;
	const double hs=parms->maos.hs;
	const double scale=pow(1.-hc/hs, -2);
	const double scale1=1.-scale;
	double focus;
	if(modm->nx>5){
	    focus=mod[5];
	    if(!parms->maos.ahstfocus){
		focus+=mod[2]*scale1;
	    }
	}else{
	    focus=mod[2]*scale1;
	}
	for(int iy=0; iy<nx; iy++){
	    double y=oy+iy*dx;
	    for(int ix=0; ix<nx; ix++){
		double x=ox+ix*dx;
		double xy=x*y;
		double x2=x*x;
		double y2=y*y;
		double tmp= 
		    +x*mod[0]
		    +y*mod[1]
		    +focus*(x2+y2)
		    +mod[2]*(-2*scale*hc*(thetax*x+thetay*y))
		    +mod[3]*((x2-y2)*scale1 - 2*scale*hc*(thetax*x-thetay*y))
		    +mod[4]*(xy*scale1-scale*hc*(thetay*x+thetax*y));
		p[ix+wvf->nx*iy]*=cexp(ik*tmp);
	    }
	}
    }
}


/**
   Time domain physical simulation.
   
   noisy: 
   - 0: no noise at all; 
   - 1: poisson and read out noise. 
   - 2: only poisson noise.   
*/
dmat *skysim_phy(dmat **mresout, const dmat *mideal, const dmat *mideal_oa, double ngsol, 
		 ASTER_S *aster, const POWFS_S *powfs, 
		 const PARMS_S *parms, int idtratc, int noisy, int phystart){
    int dtratc=0;
    if(!parms->skyc.multirate){
	dtratc=parms->skyc.dtrats[idtratc];
    }
    int hasphy;
    if(phystart>-1 && phystart<aster->nstep){
	hasphy=1;
    }else{
	hasphy=0;
    }
    const int nmod=mideal->nx;
    PDMAT(mideal,pmideal);
    PDMAT(mideal_oa, pmideal_oa);
    dmat *res=dnew(6,1);/*Results. 1-2: NGS and TT modes., 
			  3-4:On axis NGS and TT modes,
			  4-6: On axis NGS and TT wihtout considering un-orthogonality.*/
    dmat *mreal=NULL;/*modal correction at this step. */
    dmat *merr=dnew(nmod,1);/*modal error */
    dcell *merrm=dcellnew(1,1);dcell *pmerrm=NULL;
    const int nstep=aster->nstep?aster->nstep:parms->maos.nstep;
    dmat *mres=dnew(nmod,nstep);
    PDMAT(mres,pmres);
    PDMAT(parms->skyc.rnefs,rnefs);
    dcell *zgradc=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
    dcell *gradout=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
    dmat *gradsave=0;
    if(parms->skyc.dbg){
	gradsave=dnew(aster->tsa*2,nstep);
    }
   
    
    SERVO_T *st2t=0;
    kalman_t *kalman=0;
    dcell *mpsol=0;
    dmat *pgm=0;
    dmat *dtrats=0;
    int multirate=parms->skyc.multirate;
    if(multirate){
	kalman=aster->kalman[0];
	dtrats=aster->dtrats;
    }else{
	if(parms->skyc.servo>0){
	    const double dtngs=parms->maos.dt*dtratc;
	    st2t=servo_new(merrm, NULL, 0, dtngs, aster->gain->p[idtratc]);
	    pgm=aster->pgm->p[idtratc];
	}else{
	    kalman=aster->kalman[idtratc];
	}
    }
    if(kalman){
	kalman_init(kalman);
	mpsol=dcellnew(aster->nwfs, 1); //for psol grad.
    }
    const long nwvl=parms->maos.nwvl;
    dcell **psf=0, **mtche=0, **ints=0;
    ccell *wvf=0,*wvfc=0, *otf=0;
    if(hasphy){
	psf=calloc(aster->nwfs, sizeof(dcell*));
	wvf=ccellnew(aster->nwfs,1);
	wvfc=ccellnew(aster->nwfs,1);
	mtche=calloc(aster->nwfs, sizeof(dcell*));
	ints=calloc(aster->nwfs, sizeof(dcell*));
	otf=ccellnew(aster->nwfs,1);
    
	for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
	    const int ipowfs=aster->wfs[iwfs].ipowfs;
	    const long ncomp=parms->maos.ncomp[ipowfs];
	    const long nsa=parms->maos.nsa[ipowfs];
	    wvf->p[iwfs]=cnew(ncomp,ncomp);
	    wvfc->p[iwfs]=NULL;
	    psf[iwfs]=dcellnew(nsa,nwvl);
	    cfft2plan(wvf->p[iwfs], -1);
	    if(parms->skyc.multirate){
		mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[(int)aster->idtrats->p[iwfs]];
	    }else{
		mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[idtratc];
	    }
	    otf->p[iwfs]=cnew(ncomp,ncomp);
	    cfft2plan(otf->p[iwfs],-1);
	    cfft2plan(otf->p[iwfs],1);
	    ints[iwfs]=dcellnew(nsa,1);
	    int pixpsa=parms->skyc.pixpsa[ipowfs];
	    for(long isa=0; isa<nsa; isa++){
		ints[iwfs]->p[isa]=dnew(pixpsa,pixpsa);
	    }
	}
    }
    for(int irep=0; irep<parms->skyc.navg; irep++){
	if(kalman){
	    kalman_init(kalman);
	}else{
	    servo_reset(st2t);
	}
	dcellzero(zgradc);
	dcellzero(gradout);
	if(ints){
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dcellzero(ints[iwfs]);
	    }
	}
	for(int istep=0; istep<nstep; istep++){
	    memcpy(merr->p, pmideal[istep], nmod*sizeof(double));
	    dadd(&merr, 1, mreal, -1);/*form NGS mode error; */
	    memcpy(pmres[istep],merr->p,sizeof(double)*nmod);
	    if(mpsol){//collect averaged modes for PSOL.
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    dadd(&mpsol->p[iwfs], 1, mreal, 1);
		}
	    }
	    pmerrm=0;
	    if(istep>=parms->skyc.evlstart){/*performance evaluation*/
		double res_ngs=dwdot(merr->p,parms->maos.mcc,merr->p);
		if(res_ngs>ngsol*100){
		    dfree(res); res=NULL;
		    break;
		}
		{
		    res->p[0]+=res_ngs;
		    res->p[1]+=dwdot2(merr->p,parms->maos.mcc_tt,merr->p);
		    double dot_oa=dwdot(merr->p, parms->maos.mcc_oa, merr->p);
		    double dot_res_ideal=dwdot(merr->p, parms->maos.mcc_oa, pmideal[istep]);
		    double dot_res_oa=0;
		    for(int imod=0; imod<nmod; imod++){
			dot_res_oa+=merr->p[imod]*pmideal_oa[istep][imod];
		    }
		    res->p[2]+=dot_oa-2*dot_res_ideal+2*dot_res_oa;
		    res->p[4]+=dot_oa;
		}
		{
		    double dot_oa_tt=dwdot2(merr->p, parms->maos.mcc_oa_tt, merr->p);
		    /*Notice that mcc_oa_tt2 is 2x5 marix. */
		    double dot_res_ideal_tt=dwdot(merr->p, parms->maos.mcc_oa_tt2, pmideal[istep]);
		    double dot_res_oa_tt=0;
		    for(int imod=0; imod<2; imod++){
			dot_res_oa_tt+=merr->p[imod]*pmideal_oa[istep][imod];
		    }
		    res->p[3]+=dot_oa_tt-2*dot_res_ideal_tt+2*dot_res_oa_tt;
		    res->p[5]+=dot_oa_tt;
		}
	    }//if evl

	    if(istep<phystart || phystart<0){
		/*Ztilt, noise free simulation for acquisition. */
		dmm(&zgradc->m, 1, aster->gm, merr, "nn", 1);/*grad due to residual NGS mode. */
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    const int ipowfs=aster->wfs[iwfs].ipowfs;
		    const long ng=parms->maos.nsa[ipowfs]*2;
		    for(long ig=0; ig<ng; ig++){
			zgradc->p[iwfs]->p[ig]+=aster->wfs[iwfs].ztiltout->p[istep*ng+ig];
		    }
		}
	
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=(multirate?(int)dtrats->p[iwfs]:dtratc);
		    if((istep+1) % dtrati==0){
			dadd(&gradout->p[iwfs], 0, zgradc->p[iwfs], 1./dtrati);
			dzero(zgradc->p[iwfs]);
			if(noisy){
			    int idtrati=(multirate?(int)aster->idtrats->p[iwfs]:idtratc);
			    dmat *nea=aster->wfs[iwfs].pistat->sanea->p[idtrati];
			    for(int i=0; i<nea->nx; i++){
				gradout->p[iwfs]->p[i]+=nea->p[i]*randn(&aster->rand);
			    }
			}
			pmerrm=merrm;//record output.
		    }
		}
	    }else{
		/*Accumulate PSF intensities*/
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    const double thetax=aster->wfs[iwfs].thetax;
		    const double thetay=aster->wfs[iwfs].thetay;
		    const int ipowfs=aster->wfs[iwfs].ipowfs;
		    const long nsa=parms->maos.nsa[ipowfs];
		    PCCELL(aster->wfs[iwfs].wvfout[istep],wvfout);
		    for(long iwvl=0; iwvl<nwvl; iwvl++){
			double wvl=parms->maos.wvl[iwvl];
			for(long isa=0; isa<nsa; isa++){
			    ccp(&wvfc->p[iwfs], wvfout[iwvl][isa]);
			    /*Apply NGS mode error to PSF. */
			    ngsmod2wvf(wvfc->p[iwfs], wvl, merr, powfs+ipowfs, isa,
				       thetax, thetay, parms);
			    cembedc(wvf->p[iwfs],wvfc->p[iwfs],0,C_FULL);
			    cfft2(wvf->p[iwfs],-1);
			    /*peak in corner. */
			    cabs22d(&psf[iwfs]->p[isa+nsa*iwvl], 1., wvf->p[iwfs], 1.);
			}/*isa */
		    }/*iwvl */
		}/*iwfs */
	
		/*Form detector image from accumulated PSFs*/
		double igrad[2];
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=dtratc, idtrat=idtratc;
		    if(multirate){//multirate
			idtrat=aster->idtrats->p[iwfs];
			dtrati=dtrats->p[iwfs];
		    }
		    if((istep+1) % dtrati == 0){/*has output */
			dcellzero(ints[iwfs]);
			const int ipowfs=aster->wfs[iwfs].ipowfs;
			const long nsa=parms->maos.nsa[ipowfs];
			for(long isa=0; isa<nsa; isa++){
			    for(long iwvl=0; iwvl<nwvl; iwvl++){
				double siglev=aster->wfs[iwfs].siglev->p[iwvl];
				ccpd(&otf->p[iwfs],psf[iwfs]->p[isa+nsa*iwvl]);
				cfft2i(otf->p[iwfs], 1); /*turn to OTF, peak in corner */
				ccwm(otf->p[iwfs], powfs[ipowfs].dtf[iwvl].nominal);
				cfft2(otf->p[iwfs], -1);
				spmulcreal(ints[iwfs]->p[isa]->p, powfs[ipowfs].dtf[iwvl].si, 
					   otf->p[iwfs]->p, siglev);
			    }
		
			    /*Add noise and apply matched filter. */
#if _OPENMP >= 200805 
#pragma omp critical 
#endif
			    switch(noisy){
			    case 0:/*no noise at all. */
				break;
			    case 1:/*both poisson and read out noise. */
				addnoise(ints[iwfs]->p[isa], &aster->rand, aster->wfs[iwfs].bkgrnd*dtrati,
					 1, rnefs[ipowfs][idtrat]);
				break;
			    case 2:/*there is still poisson noise. */
				addnoise(ints[iwfs]->p[isa], &aster->rand, 0, 1, 0);
				break;
			    default:
				error("Invalid noisy\n");
			    }
		
			    igrad[0]=0;
			    igrad[1]=0;
			    double pixtheta=parms->skyc.pixtheta[ipowfs];
			    if(parms->skyc.mtch){
				dmulvec(igrad, mtche[iwfs]->p[isa], ints[iwfs]->p[isa]->p, 1);
			    }
			    if(!parms->skyc.mtch || fabs(igrad[0])>pixtheta || fabs(igrad[1])>pixtheta){
				if(!parms->skyc.mtch){
				    warning2("fall back to cog\n");
				}else{
				    warning_once("mtch is out of range\n");
				}
				dcog(igrad, ints[iwfs]->p[isa], 0, 0, 0, 3*rnefs[ipowfs][idtrat]); 
				igrad[0]*=pixtheta;
				igrad[1]*=pixtheta;
			    }
			    gradout->p[iwfs]->p[isa]=igrad[0];
			    gradout->p[iwfs]->p[isa+nsa]=igrad[1];
			}/*isa */
			pmerrm=merrm;
			dcellzero(psf[iwfs]);/*reset accumulation.*/
		    }/*if iwfs has output*/
		}/*for wfs*/
	    }/*if phystart */
	    //output to mreal after using it to ensure two cycle delay.
	    if(st2t){//Type I or II control.
		if(st2t->mint[0]){//has output.
		    dcp(&mreal, st2t->mint[0]->p[0]);
		}
	    }else{//LQG control
		kalman_output(kalman, &mreal, 0, 1);
	    }
	    if(kalman){//LQG control
		int indk=0;
		//Form PSOL grads and obtain index to LQG M
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		    int dtrati=(multirate?(int)dtrats->p[iwfs]:dtratc);
		    if((istep+1) % dtrati==0){
			indk|=1<<iwfs;
			dmm(&gradout->p[iwfs], 1, aster->g->p[iwfs], mpsol->p[iwfs], "nn", 1./dtrati);
			dzero(mpsol->p[iwfs]);
		    }
		}
		if(indk){
		    kalman_update(kalman, gradout->m, indk-1);
		}
	    }else if(st2t){
		if(pmerrm){
		    dmm(&merrm->p[0], 0, pgm, gradout->m, "nn", 1);	
		}
		servo_filter(st2t, pmerrm);//do even if merrm is zero. to simulate additional latency
	    }
	    if(parms->skyc.dbg){
		memcpy(gradsave->p+istep*gradsave->nx, gradout->m->p, sizeof(double)*gradsave->nx);
	    }
	}/*istep; */
    }
    if(parms->skyc.dbg){
	int dtrati=(multirate?(int)dtrats->p[0]:dtratc);
	dwrite(gradsave,"%s/skysim_grads_aster%d_dtrat%d",dirsetup, aster->iaster,dtrati);
	dwrite(mres,"%s/skysim_phy_mres_aster%d_dtrat%d",dirsetup,aster->iaster,dtrati);
    }
  
    dfree(mreal);
    dcellfree(mpsol);
    dfree(merr);
    dcellfree(merrm);
    dcellfree(zgradc);
    dcellfree(gradout);
    dfree(gradsave);
    if(hasphy){
	dcellfreearr(psf, aster->nwfs);
	dcellfreearr(ints, aster->nwfs);
        ccellfree(wvf);
	ccellfree(wvfc);
	ccellfree(otf);
	free(mtche);
    }
    servo_free(st2t);
    /*dfree(mres); */
    if(mresout) {
	*mresout=mres;
    }else{
	dfree(mres);
    }
    dscale(res, 1./((nstep-parms->skyc.evlstart)*parms->skyc.navg));
    return res;
}

/**
   Save NGS WFS and other information for later use in MAOS simulations.*/
void skysim_save(const SIM_S *simu, const ASTER_S *aster, const double *ipres, int selaster, int seldtrat, int isky){
    const PARMS_S* parms=simu->parms;
    const int nwvl=parms->maos.nwvl;
    char path[PATH_MAX];
    snprintf(path,PATH_MAX,"Res%d_%d_maos/sky%d",simu->seed_maos,parms->skyc.seed,isky);
    mymkdir("%s",path);
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	dcell *sepsf=dcelldup(aster[selaster].wfs[iwfs].pistat->psf);
	for(int ic=0; ic<sepsf->nx*sepsf->ny; ic++){
	    dfftshift(sepsf->p[ic]);/*put peak in center. required by MAOS. */
	}
	dcellwrite(sepsf, "%s/pistat_wfs%d",path,iwfs+6);
	dcellfree(sepsf);
	dwrite(aster->wfs[iwfs].pistat->sanea->p[seldtrat], "%s/nea_tot_wfs%d",path,iwfs+6);
	dwrite(aster[selaster].wfs[iwfs].pistat->sanea->p[seldtrat], 
	       "%s/nea_wfs%d",path,iwfs+6);
	dcellwrite(aster[selaster].wfs[iwfs].pistat->sanea, 
		   "%s/neafull_wfs%d",path,iwfs+6);
    }
    dwrite(aster[selaster].gain->p[seldtrat], "%s/gain",path);
    dwrite(simu->mres->p[isky], "%s/mres",path);
    dwrite(simu->psd_tt,"%s/psd_tt",path);
    dwrite(simu->psd_ps,"%s/psd_ps",path);
    char fnconf[PATH_MAX];
    snprintf(fnconf,PATH_MAX,"%s/base.conf",path);
    FILE *fp=fopen(fnconf,"w");

    fprintf(fp,"sim.seeds=[%d]\n",simu->seed_maos);
    fprintf(fp,"sim.end=%d\n", parms->maos.nstep);
    fprintf(fp,"sim.dt=%g\n", parms->maos.dt);
    fprintf(fp,"sim.zadeg=%g\n", parms->maos.zadeg);
    fprintf(fp,"sim.mffocus=%d\n", parms->maos.mffocus);
    fprintf(fp,"sim.ahstfocus=%d\n", parms->maos.ahstfocus);
    fprintf(fp,"tomo.ahst_wt=3\n");
    fprintf(fp,"sim.servotype_lo=2\n");/*type II */
    fprintf(fp,"sim.eplo='gain.bin'\n");
    fprintf(fp,"powfs0_llt.fnrange='%s'\n", parms->maos.fnrange);
    fprintf(fp,"atm.r0z=%.4f\n", parms->maos.r0z);
    fprintf(fp,"atm.size=[128 128]\n");
    if(parms->maos.wddeg){
	fprintf(fp, "atm.wddeg=[");
	for(int ips=0; ips<parms->maos.nwddeg; ips++){
	    fprintf(fp, "%.2f ", parms->maos.wddeg[ips]);
	}
	fprintf(fp, "]\n");
    }
    fprintf(fp,"wfs.thetax=[0 0  -33.287 -20.5725  20.5725 33.287");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetax*206265);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.thetay=[0 35 10.8156 -28.3156 -28.3156 10.8156");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetay*206265);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"wfs.siglev=[900 900 900 900 900 900");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp, " %.2f", aster[selaster].wfs[iwfs].siglevtot);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.wvlwts=[1 1 1 1 1 1");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    fprintf(fp," %.2f ", aster[selaster].wfs[iwfs].siglev->p[iwvl]
		    /aster[selaster].wfs[iwfs].siglevtot);
	}
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.powfs=[0 0 0 0 0 0");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp, " %d", aster[selaster].wfs[iwfs].ipowfs+1);
    }
    fprintf(fp,"]\n");

    PDMAT(parms->skyc.rnefs,rnefs);
    double rne=rnefs[0][seldtrat];
    double bkgrnd=aster[selaster].wfs[0].bkgrnd;

    if(parms->maos.npowfs==1){
	fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" ]\n");
	fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\"]\n");
	fprintf(fp, "powfs.phyusenea=[0 1]\n");
	fprintf(fp, "powfs.dtrat=[1 %d]\n", parms->skyc.dtrats[seldtrat]);
	fprintf(fp, "powfs.bkgrnd=[0 %.2f]\n", bkgrnd);
	fprintf(fp, "powfs.rne=[3 %.2f]\n", rne);
	fprintf(fp, "powfs.phystep=[0 %ld]\n", (long)50+parms->skyc.dtrats[seldtrat]*20);
	fprintf(fp, "powfs.noisy=[1 1 ]\n");
	fprintf(fp, "powfs.pixtheta=[0.5/206265 %g/206265000]\n", parms->skyc.pixtheta[1]*206265000);
	fprintf(fp, "powfs.pixpsa=[6 %d]\n", parms->skyc.pixpsa[0]);
	fprintf(fp, "powfs.ncomp=[64 %d]\n", parms->maos.ncomp[0]);
	fprintf(fp, "powfs.nwvl=[1 %d]\n",nwvl);
	fprintf(fp, "powfs.wvl=[0.589e-6");
	for(int ip=0; ip<1; ip++){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
	    }
	}
	fprintf(fp,"]\n");
    }else if(parms->maos.npowfs==2){
	fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" \"pistat\"]\n");
	fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\" \"nea_tot\"]\n");
	fprintf(fp, "powfs.phyusenea=[0 1 1]\n");
	fprintf(fp, "powfs.dtrat=[1 %d %d]\n", parms->skyc.dtrats[seldtrat],
		parms->skyc.dtrats[seldtrat]);
	fprintf(fp, "powfs.bkgrnd=[0 %.2f %.2f]\n", bkgrnd, bkgrnd);
	fprintf(fp, "powfs.rne=[3 %.2f %.2f]\n", rne,rne);
	fprintf(fp, "powfs.phystep=[0 %ld %ld]\n", 
		(long)50+parms->skyc.dtrats[seldtrat]*20, 
		(long)50+parms->skyc.dtrats[seldtrat]*20);
	fprintf(fp, "powfs.noisy=[1 1 1]\n");
	fprintf(fp, "powfs.pixtheta=[0.5/206265 %g/206265000 %g/206265000]\n",
		parms->skyc.pixtheta[0]*206265000,
		parms->skyc.pixtheta[1]*206265000);
	fprintf(fp, "powfs.pixpsa=[6 %d %d]\n",
		parms->skyc.pixpsa[0], 
		parms->skyc.pixpsa[1]);
	fprintf(fp, "powfs.ncomp=[64 %d %d]\n", 
		parms->maos.ncomp[0], parms->maos.ncomp[1]);
	fprintf(fp, "powfs.nwvl=[1 %d %d]\n",nwvl,nwvl);
	fprintf(fp, "powfs.wvl=[0.589e-6");
	for(int ip=0; ip<2; ip++){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
	    }
	}
	fprintf(fp,"]\n");
	fprintf(fp, "powfs.wvlwts=[]\n");
    }else{
	error("Fill this out please\n");
    }
 
    fclose(fp);
    snprintf(fnconf,PATH_MAX,"%s/skyres.txt",path);
    fp=fopen(fnconf,"w");
    fprintf(fp, "TotAll\tNGS\tTT\n");
    fprintf(fp, "%g\t%g\t%g\n",
	    sqrt(ipres[0])*1e9, sqrt(ipres[1])*1e9, sqrt(ipres[2])*1e9);
    fclose(fp);
}
