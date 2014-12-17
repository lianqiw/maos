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
#include "common.h"
#include "sim.h"
#include "sim_utils.h"
#include "ahst.h"
#include "mtch.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file wfsgrad.c
   contains functions that computes WFS gradients in geometric or physical optics mode.
*/
#define TIMING 0
#if TIMING == 1
#define TIM(A) double tk##A=myclockd()
#else
#define TIM(A)
#endif

static void wfs_ideal_atm(SIM_T *simu, dmat *opd, int iwfs, double alpha){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const double hs=parms->wfs[iwfs].hs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    for(int idm=0; idm<parms->ndm; idm++){
	loc_t *loc=powfs[ipowfs].loc_dm?powfs[ipowfs].loc_dm->p[wfsind+idm*parms->nwfs]:powfs[ipowfs].loc;
	double *amp=powfs[ipowfs].realamp->p[wfsind]->p;
	const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	double dispx=ht*parms->wfs[iwfs].thetax;
	double dispy=ht*parms->wfs[iwfs].thetay;
	double scale=1.-ht/hs;
	if(parms->dm[idm].cubic){
	    prop_grid_cubic(simu->dmprojsq->p[idm],
			    loc, amp, opd->p, 
			    alpha, dispx, dispy, scale, parms->dm[idm].iac, 
			    0, 0);
	}else{
	    prop_grid(simu->dmprojsq->p[idm],
		      loc, amp, opd->p, 
		      alpha, dispx, dispy, scale, 0,
		      0, 0);
	}
    }
}

/**
   computes close loop and pseudo open loop gradidents for both gometric and
   physical optics WFS. Calls wfsints() to accumulate WFS subapertures images in
   physical optics mode.  */

void wfsgrad_iwfs(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int isim=simu->isim;
    const int iwfs=info->start;
    const PARMS_T *parms=simu->parms;
    const int ipowfs=parms->wfs[iwfs].powfs;
    if(isim<parms->powfs[ipowfs].step) return;
    assert(iwfs<parms->nwfs);
    /*
      simu->gradcl is CL grad output (also for warm-restart of maxapriori
      simu->gradacc is internal, to accumulate geometric grads.
      do not accumulate opd. accumate ints for phy, g for GS
    */
    /*input */
    
    mapcell *atm=simu->atm;
    const RECON_T *recon=simu->recon;
    const POWFS_T *powfs=simu->powfs;
    /*output */
    const int CL=parms->sim.closeloop;
    const int nps=parms->atm.nps;
    const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
    const double dt=simu->dt;
    TIM(0);
    /*The following are truly constants for this powfs */
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const double hs=parms->wfs[iwfs].hs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom->p[iwfs];
    const int save_opd =parms->save.wfsopd->p[iwfs];
    const int save_ints=parms->save.ints->p[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    /*The following depends on isim */
    /*const int dtrat_reset=(isim%dtrat==0); */
    const int dtrat_output=(isim+1-parms->powfs[ipowfs].step)%dtrat==0;
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_pistatout=parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart;
    const int do_geom=!do_phy || save_gradgeom || do_pistatout;
    const double *realamp=powfs[ipowfs].realamp->p[wfsind]->p;
    double *srot=(do_phy && parms->powfs[ipowfs].radpix)?
	powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p:NULL;
    dmat *gradcalc=NULL;
    dmat **gradacc=&simu->gradacc->p[iwfs];
    dmat **gradout=&simu->gradcl->p[iwfs];
    dcell *ints=simu->ints->p[iwfs];
    dmat  *opd=simu->wfsopd->p[iwfs];
    dzero(opd);
    if((isim-parms->powfs[ipowfs].step)%dtrat==0){
	dcellzero(ints);
	dzero(*gradacc);
    }
    if(simu->telws){/*Wind shake */
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir?simu->winddir->p[0]:0;
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
    }

    /* Add surface error*/
    if(powfs[ipowfs].opdadd && powfs[ipowfs].opdadd->p[wfsind]){
	dadd(&opd,1, powfs[ipowfs].opdadd->p[wfsind],1);
    }

    /* Now begin ray tracing. */
    if(parms->sim.idealwfs && !parms->powfs[ipowfs].lo){
	wfs_ideal_atm(simu, opd, iwfs, 1);
    }else if(atm){
	for(int ips=0; ips<nps; ips++){
	    thread_t *wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
	    wfs_propdata->phiout=opd->p;
	    wfs_propdata->displacex1=-atm->p[ips]->vx*dt*isim;
	    wfs_propdata->displacey1=-atm->p[ips]->vy*dt*isim;
	    wfs_propdata->alpha=atmscale;
	    /* have to wait to finish before another phase screen. */
	    CALL_THREAD(wfs_prop, 0);
	}
	/* most expensive 0.10 per LGS for*/
	if(parms->sim.wfsalias){
	    /* Remove subspace of atm projected onto range of DM.*/
	    wfs_ideal_atm(simu, opd, iwfs, -1);
	}
    }
    if(save_opd){
	cellarr_dmat(simu->save->wfsopdol[iwfs], isim, opd);
    }
 
    if(CL){
	for(int idm=0; idm<parms->ndm; idm++){
	    thread_t *wfs_prop=simu->wfs_prop_dm[iwfs+parms->nwfs*idm];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
	    wfs_propdata->phiout=opd->p;
	    CALL_THREAD(wfs_prop, 0);
	}/*idm */
	if(simu->ttmreal){
	    double ptt[3]={0, -simu->ttmreal->p[0], -simu->ttmreal->p[1]};
	    loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
	}
    }
    if(parms->powfs[ipowfs].skip && parms->tomo.ahst_idealngs){
	//apply ideal NGS modes to NGS WFS
	ngsmod2science(opd, powfs[ipowfs].loc, recon->ngsmod, 
		       parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay,
		       simu->cleNGSm->p+isim*recon->ngsmod->nmod, -1);
    }
    if(imoao>-1){
	dmat **dmwfs=simu->dm_wfs->p;
	if(dmwfs[iwfs]){
	    /* No need to do mis registration here since the MOAO DM is attached
	       to close to the WFS.*/
	    if(parms->moao[imoao].cubic){
		prop_nongrid_pts_cubic(recon->moao[imoao].aloc->p[0], dmwfs[iwfs]->p,
				       powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				       parms->moao[imoao].iac, 0, 0);
	    }else{
		prop_nongrid_pts(recon->moao[imoao].aloc->p[0], dmwfs[iwfs]->p,
				 powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				 0, 0);
	    }
	}
    }
    /* Add defocus to OPD if needed. */
    if(parms->powfs[ipowfs].llt){
	double focus=wfsfocusadj(simu, iwfs);
	if(fabs(focus)>1e-20){
	    loc_add_focus(opd->p, powfs[ipowfs].loc, focus);
	}
    }
    if(parms->powfs[ipowfs].fieldstop>0){
	locfft_fieldstop(powfs[ipowfs].fieldstop, opd, parms->powfs[ipowfs].wvlwts);
    }

    if(save_opd){
	cellarr_dmat(simu->save->wfsopd[iwfs], isim, opd);
    }
    if(parms->plot.run){
	drawopdamp("wfsopd",powfs[ipowfs].loc,opd->p,realamp,NULL,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
    }
    if(do_geom){
	/* Now Geometric Optics gradient calculations. if dtrat==1, we compute
	   gradients directly to gradacc, which is the same as gradcalc. If
	   dtrat>1, we compute gradients to gradcalc, and accumulate to
	   gradacc. gradcalc is used to shift pistat. We DONOT include gradoff
	   adjustment to gradref, but only do it on gradcl. This will make the
	   pistat always peak in center no matter what NCPA is present.
	*/
	if(!do_pistatout || parms->powfs[ipowfs].pistatstc || dtrat==1){
	    gradcalc=dref(*gradacc);
	}//else: calculate first to gradcalc then add to gradacc
	if(parms->powfs[ipowfs].gtype_sim==1){ /*compute ztilt. */
	    pts_ztilt(&gradcalc,powfs[ipowfs].pts,
		      powfs[ipowfs].saimcc->p[powfs[ipowfs].nsaimcc>1?wfsind:0], 
		      realamp, opd->p);
	}else{/*G tilt */
	    dspmm(&gradcalc,adpind(powfs[ipowfs].GS0,wfsind),opd,'n',1);
	}
	if(gradcalc->p!=(*gradacc)->p){
	    dadd(gradacc, 1, gradcalc, 1);
	}
    }

    ccell *psfout=NULL;
    cellarr *psfoutcellarr=NULL;
    cellarr *ztiltoutcellarr=NULL;
    if(parms->powfs[ipowfs].psfout){
	psfout=simu->wfspsfout->p[iwfs];
	psfoutcellarr=simu->save->wfspsfout[iwfs];
	ztiltoutcellarr=simu->save->ztiltout[iwfs];
    }
    TIM(1);
    /* Now begin Physical Optics Intensity calculations */
    if(do_phy || psfout || do_pistatout || parms->powfs[ipowfs].dither){
	dmat *lltopd=NULL;
	if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
	    if(powfs[ipowfs].llt->ncpa){
		int iotf=powfs[ipowfs].llt->ncpa->nx==1?0:wfsind;
		lltopd=ddup(powfs[ipowfs].llt->ncpa->p[iotf]);
	    }else{
		lltopd=dnew(powfs[ipowfs].llt->pts->nx,
			    powfs[ipowfs].llt->pts->nx);
	    }
	    const long illt=parms->powfs[ipowfs].llt->i->p[wfsind];
	    if(atm){/*LLT OPD */
		for(int ips=0; ips<nps; ips++){
		    const double hl=atm->p[ips]->h;
		    const double scale=1.-hl/hs;
		    const double thetax=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox->p[illt]/hs;
		    const double thetay=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy->p[illt]/hs;
		    const double displacex=-atm->p[ips]->vx*isim*dt+thetax*hl+parms->powfs[ipowfs].llt->misreg->p[0];
		    const double displacey=-atm->p[ips]->vy*isim*dt+thetay*hl+parms->powfs[ipowfs].llt->misreg->p[1];
		    prop_grid_pts(atm->p[ips],powfs[ipowfs].llt->pts,NULL,
				  lltopd->p,atmscale,displacex,displacey,
				  scale, 1., 0, 0);
		}
	    }
	    double ttx=0, tty=0;//FSM + wind shake induced jitter
	    if((simu->uptreal && simu->uptreal->p[iwfs]) ||do_pistatout||parms->sim.uptideal){
		if(do_pistatout||parms->sim.uptideal){
		    /* remove tip/tilt completely */
		    dmat *lltg=dnew(2,1);
		    pts_ztilt(&lltg,powfs[ipowfs].llt->pts,
			      powfs[ipowfs].llt->imcc,
			      powfs[ipowfs].llt->amp->p,
			      lltopd->p);
		    ttx=-lltg->p[0];
		    tty=-lltg->p[1];
		    dfree(lltg);
		}else{
		    ttx=simu->uptreal->p[iwfs]->p[0];
		    tty=simu->uptreal->p[iwfs]->p[1];
		}
		/* copy uptreal to output */
		PDMAT(simu->uptcmds->p[iwfs], puptcmds);
		puptcmds[isim][0]=ttx;
		puptcmds[isim][1]=tty;
	    }
	    if(simu->telws){
		double tmp=simu->telws->p[isim]*parms->powfs[ipowfs].llt->ttrat;
		double angle=simu->winddir?simu->winddir->p[0]:0;
		ttx+=tmp*cos(angle);
		tty+=tmp*sin(angle);
	    }
	    if(ttx !=0 || tty != 0){ /* add tip/tilt to llt opd */
		double ptt[3]={0, ttx, tty};
		loc_add_ptt(lltopd->p, ptt, powfs[ipowfs].llt->loc);
	    }
	    if(save_opd){
		cellarr_dmat(simu->save->wfslltopd[iwfs], isim, lltopd);
	    }
	}
	WFSINTS_T *intsdata=simu->wfs_intsdata+iwfs;
	intsdata->ints=ints;
	intsdata->psfout=psfout;
	intsdata->pistatout=simu->pistatout->p[iwfs];
	if(parms->powfs[ipowfs].pistatout==1){
	    intsdata->gradref=gradcalc;
	}
	intsdata->opd=opd;
	intsdata->lltopd=lltopd;
	intsdata->isim=isim;
	CALL_THREAD(simu->wfs_ints[iwfs], 0);
	dfree(lltopd);
	intsdata->opd=0;
	intsdata->lltopd=0;
	if(psfout){
	    cellarr_ccell(psfoutcellarr, isim, psfout);
	    cellarr_dmat(ztiltoutcellarr, isim, *gradacc);
	}
    }
    TIM(2);
    if(dtrat_output){
	const double rne=parms->powfs[ipowfs].rne;
	const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	if(do_phy || parms->powfs[ipowfs].dither){
	    /* In Physical optics mode, do integration and compute
	       gradients. The matched filter are in x/y coordinate even if
	       radpix=1. */
	    if(save_ints){
		cellarr_dcell(simu->save->intsnf[iwfs], isim, ints);
	    }
	    if(noisy){/*add noise */
		const double bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
		dmat **bkgrnd2=NULL;
		dmat **bkgrnd2c=NULL;
		if(powfs[ipowfs].bkgrnd){
		    if(powfs[ipowfs].bkgrnd->ny==1){
			bkgrnd2=powfs[ipowfs].bkgrnd->p;
		    }else{
			bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*wfsind;
		    }
		}
		if(powfs[ipowfs].bkgrndc){
		    if(powfs[ipowfs].bkgrndc->ny==1){
			bkgrnd2c=powfs[ipowfs].bkgrndc->p;
		    }else{
			bkgrnd2c=powfs[ipowfs].bkgrndc->p+nsa*wfsind;
		    }
		}
		for(int isa=0; isa<nsa; isa++){
		    dmat *bkgrnd2i=(bkgrnd2)?bkgrnd2[isa]:NULL;
		    dmat *bkgrnd2ic=(bkgrnd2c)?bkgrnd2c[isa]:NULL;
		    addnoise(ints->p[isa], &simu->wfs_rand[iwfs],
			     bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, rne);
		}
		if(save_ints){
		    cellarr_dcell(simu->save->intsny[iwfs], isim, ints);
		}
	    }
	    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_nskip){
		/*Collect statistics with dithering*/
		DITHER_T *pd=simu->dither[iwfs];
		double cs, ss;
		dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
		dcelladd(&pd->imb, 1, ints, 1.);
		dcelladd(&pd->imx, 1, ints, cs);
		dcelladd(&pd->imy, 1, ints, ss);
	    }
	}
	if(do_phy){
	    dmat **mtche=NULL;
	    double *i0sum=NULL;
	    if(parms->powfs[ipowfs].phytypesim==1){
		if(powfs[ipowfs].intstat->mtche->ny==1){
		    mtche=powfs[ipowfs].intstat->mtche->p;
		    i0sum=powfs[ipowfs].intstat->i0sum->p;
		}else{
		    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
		    i0sum=powfs[ipowfs].intstat->i0sum->p+nsa*wfsind;
		}
	    }
	    double pixthetax=parms->powfs[ipowfs].radpixtheta;
	    double pixthetay=parms->powfs[ipowfs].pixtheta;
	    /*output directly to simu->gradcl. replace */
	    double *pgradx=(*gradout)->p;
	    double *pgrady=pgradx+nsa;
	    for(int isa=0; isa<nsa; isa++){
		double geach[3]={0,0,1};
		switch(parms->powfs[ipowfs].phytypesim){
		case 1:{
		    dmulvec(geach, mtche[isa],ints->p[isa]->p,1.);
		    if(parms->powfs[ipowfs].mtchscl){
			double scale=i0sum[isa]/dsum(ints->p[isa]);
			geach[0]*=scale;
			geach[1]*=scale;
		    }
		}
		    break;
		case 2:{
		    dcog(geach,ints->p[isa],0.,0.,
			 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2],
			 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2+1]);
		    geach[0]*=pixthetax;
		    geach[1]*=pixthetay;
		    if(srot){
			double theta=srot[isa];
			double cx=cos(theta);
			double sx=sin(theta);
			double tmp=geach[0]*cx-geach[1]*sx;
			geach[1]=geach[0]*sx+geach[1]*cx;
			geach[0]=tmp;
		    }
		}
		    break;
		case 3:{
		    geach[0]=pgradx[isa];//warm restart
		    geach[1]=pgrady[isa];
		    maxapriori(geach, ints->p[isa], parms, powfs, iwfs, isa, 1, bkgrnd, rne);
		}
		    break;
		default:
		    error("Invalid");
		}
	
		pgradx[isa]=geach[0];
		pgrady[isa]=geach[1];
	    };/*isa */
	}else{
	    /* geomtric optics accumulation mode. scale and copy results to output. */
	    dcp(gradout,*gradacc);
	    if(dtrat!=1){
		dscale(*gradout,1./dtrat);/*average */
	    }
	    if(noisy && !parms->powfs[ipowfs].usephy){
		const dmat *nea=powfs[ipowfs].neasim->p[wfsind];
		const double *neax=nea->p;
		const double *neay=nea->p+nsa;
		const double *neaxy=nea->p+nsa*2;
		double *restrict ggx=(*gradout)->p;
		double *restrict ggy=(*gradout)->p+nsa;
		for(int isa=0; isa<nsa; isa++){
		    /*Preserve the random sequence. */
		    double n1=randn(&simu->wfs_rand[iwfs]);
		    double n2=randn(&simu->wfs_rand[iwfs]);
		    double errx=neax[isa]*n1;
		    double erry=neay[isa]*n2+neaxy[isa]*n1;/*cross term. */
		    ggx[isa]+=errx;
		    ggy[isa]+=erry;
		}
	    }
	}
	if(save_gradgeom){
	    dmat *gradtmp=NULL;
	    dadd(&gradtmp, 1, *gradacc, 1./dtrat);
	    cellarr_dmat(simu->save->gradgeom[iwfs], isim, gradtmp);/*noise free. */
	    dfree(gradtmp);
	}
    }//dtrat_out
    dfree(gradcalc);
    TIM(3);
#if TIMING==1
    info("wfs %d grad timing: ray %.2f ints %.2f grad %.2f\n",iwfs,tk1-tk0,tk2-tk1,tk3-tk2);
#endif
}
/**
   Save telemetry. TODO: copy from GPU to CPU.
*/
static void wfsgrad_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    if(isim+1==parms->sim.end){
	for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	    const int ipowfs=simu->parms->wfs[iwfs].powfs;
	    const int dtrat=parms->powfs[ipowfs].dtrat;
	    double scale;
	    if(parms->powfs[ipowfs].usephy){
		scale=(simu->isim+1-simu->parms->powfs[ipowfs].phystep)/dtrat;
	    }else{
		scale=(simu->isim+1)/dtrat;
	    }
	    if(scale<=0) continue;	    
	    if(simu->pistatout && simu->pistatout->p[iwfs]){
		int nstep=isim+1-parms->powfs[ipowfs].pistatstart;
		scale=1./(double)nstep;
		dcell *pp=simu->pistatout->p[iwfs];
		dcellscale(pp,scale);
		if(parms->sim.skysim){/*need peak in corner */
		    for(long ic=0; ic<pp->nx*pp->ny; ic++){
			dfftshift(pp->p[ic]);
		    }
		    writebin(pp,"%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
			       dirskysim,simu->seed,
			       parms->powfs[ipowfs].order,
			       parms->wfs[iwfs].thetax*206265,
			       parms->wfs[iwfs].thetay*206265);
		    for(long ic=0; ic<pp->nx*pp->ny; ic++){
			dfftshift(pp->p[ic]);
		    }
		}else{/*need peak in center */
		    writebin(pp,"pistat_seed%d_wfs%d.bin", simu->seed,iwfs);
		}
		dcellscale(pp,1./scale);
	    }
	}
    }
}

/*Fast steering mirror for each WFS*/
void wfsgrad_fsm(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    int isim=simu->isim;
    /*Uplink FSM*/
   
    if(!recon->PTT){
	error("powfs %d has llt, but recon->PTT is NULL",ipowfs);
    }
    dmat *PTT=recon->PTT->p[parms->recon.glao
			    ?(ipowfs+ipowfs*parms->npowfs)
			    :(iwfs+iwfs*parms->nwfs)];
    if(!PTT){
	error("powfs %d has FSM, but PTT is empty\n", ipowfs);
    }
    /* Compute LGS Uplink error. */
    simu->upterr=simu->upterr_store;
    dmm(&simu->upterr->p[iwfs], 0, PTT, simu->gradcl->p[iwfs], "nn", 1);
    /* PLL loop.*/
    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_pllskip){
	DITHER_T *pd=simu->dither[iwfs];
	{
	    double err, cd, sd;
	    dither_position(&cd, &sd, parms, ipowfs, isim, pd->deltam);
	    //Use angle of expdected averaged position during integration, correlate with error signal
	    err=(-sd*(simu->upterr->p[iwfs]->p[0])
		 +cd*(simu->upterr->p[iwfs]->p[1]))/(parms->powfs[ipowfs].dither_amp);
	    pd->delta+=parms->powfs[ipowfs].dither_gpll*err;
	}
	/*2014-10-31: 
	  To estimate the actual dithering amplitude.

	  Optionally to use LPF instead of averaging. The dithering amplitude
	  does not change over the course. However, with 240 steps of * averaging,
	  USE_SUM=1 is a touch better USE_SUM=0.

	  The RTC ADD uses FSM command (uptreal) instead of FSM error signal
	  (upterr) because the amplitude of FSM error signal is affected by the
	  gain of the current gradient estimation algorithm and is therefore
	  should not be used.
	*/
	{
	    //Don't use cd, sd above, which doesnot have unit amplitude.
	    const int adjust=parms->sim.alupt+1-parms->powfs[ipowfs].dtrat;
	    const double angle=M_PI*0.5*((isim-adjust)/parms->powfs[ipowfs].dtrat);
	    double cs=cos(angle);
	    double ss=sin(angle);
#define USE_SUM 1
	    double *fsmpos=simu->uptreal->p[iwfs]->p;
	    double ipv=(fsmpos[0]*cs+fsmpos[1]*ss);
	    double qdv=(fsmpos[0]*ss-fsmpos[1]*cs);
#if USE_SUM
	    pd->ipv+=ipv;
	    pd->qdv+=qdv;
#else
	    double gpll=parms->powfs[ipowfs].dither_gpll;
	    pd->ipv=pd->ipv*(1-gpll)+ipv*gpll;
	    pd->qdv=pd->qdv*(1-gpll)+qdv*gpll;
#endif
	}
	/*Update DLL loop measurement. The delay is about 0.2 of a
	 * cycle, according to closed loop transfer function*/
	const int npll=parms->powfs[ipowfs].dither_npll;
	int npllacc=(simu->isim-parms->powfs[ipowfs].dither_pllskip+1)/parms->powfs[ipowfs].dtrat;
	if(npllacc>0 && npllacc%npll==0){
	    pd->deltam=pd->delta;
#if USE_SUM
	    pd->a2m=sqrt(pd->ipv*pd->ipv+pd->qdv*pd->qdv)/npll;
	    if(iwfs==parms->powfs[ipowfs].wfs->p[0]){
		info2("PLL step%d, wfs%d: deltam=%.2f frame, a2m=%.1f mas\n",
		      isim, iwfs, pd->deltam/(0.5*M_PI), pd->a2m*206265000);
	    }
	    pd->ipv=pd->qdv=0;
#else
	    pd->a2m=sqrt(pd->ipv*pd->ipv+pd->qdv*pd->qdv);
#endif
	    if(simu->resdither){
		int ic=(npllacc-1)/(npll);
		simu->resdither->p[iwfs]->p[ic*2+0]=pd->deltam;
		simu->resdither->p[iwfs]->p[ic*2+1]=pd->a2m;
	    }
	}
	/* Update drift mode computation*/
	int ndrift=parms->powfs[ipowfs].dither_ndrift;
	int ndriftacc=(simu->isim-parms->powfs[ipowfs].dither_nskip+1)/parms->powfs[ipowfs].dtrat;
	if(ndriftacc>0 && ndriftacc % ndrift==0){
	    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
	    const double pixthetay=parms->powfs[ipowfs].pixtheta;
	    const int nsa=powfs[ipowfs].pts->nsa;
	    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	    double *srot=NULL;
	    if(parms->powfs[ipowfs].radpix){
		srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p;
	    }
	    dmat *ibgrad=dnew(nsa*2, 1);
	    double scale1=1./ndrift;
	    dcellscale(pd->imb, scale1);
	    double *gx=ibgrad->p;
	    double *gy=ibgrad->p+nsa;
	    double geach[3]={0,0,1};
	    for(int isa=0; isa<nsa; isa++){
		dcog(geach,pd->imb->p[isa],0.,0.,
		     powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2],
		     powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2+1]);
		geach[0]*=pixthetax;
		geach[1]*=pixthetay;
		if(srot){
		    double theta=srot[isa];
		    double cx=cos(theta);
		    double sx=sin(theta);
		    double tmp=geach[0]*cx-geach[1]*sx;
		    geach[1]=geach[0]*sx+geach[1]*cx;
		    geach[0]=tmp;
		}
		gx[isa]=geach[0];
		gy[isa]=geach[1];
	    }
	    {
		dmat *tt=dnew(2,1);
		dmm(&tt, 0, PTT, ibgrad, "nn", 1);
		simu->upterr->p[iwfs]->p[0]+=tt->p[0];
		simu->upterr->p[iwfs]->p[1]+=tt->p[1];
		dfree(tt);
	    }
	    //Smooth trombone movement by provide continuous err.
	    if(parms->powfs[ipowfs].llt){
		dmat *focus=dnew(1,1);
		dmat *RFlgsg=recon->RFlgsg->p[parms->recon.glao
					      ?(ipowfs+ipowfs*parms->npowfs)
					      :(iwfs+iwfs*parms->nwfs)];
		dmm(&focus, 0, RFlgsg, ibgrad, "nn", 1);
		simu->zoomerr->p[iwfs]=focus->p[0]/ndrift;
		dfree(focus);
	    }
	    dfree(ibgrad);
	    double scale2=scale1*2./(pd->a2m);
	    dcelladd(&pd->i0, 1, pd->imb, 1);//imb was already scaled
	    dcelladd(&pd->gx, 1, pd->imx, scale2);
	    dcelladd(&pd->gy, 1, pd->imy, scale2);
	    dcellzero(pd->imb);
	    dcellzero(pd->imx);
	    dcellzero(pd->imy);
	    if(powfs[ipowfs].gradoff){
		dadd(&pd->goff, 1, powfs[ipowfs].gradoff->p[wfsind], 1);
	    }
	}//Drift mode computation
    }//PLL loop
    /* copy upterr to output. */
    PDMAT(simu->upterrs->p[iwfs], pupterrs);
    pupterrs[isim][0]=simu->upterr->p[iwfs]->p[0];
    pupterrs[isim][1]=simu->upterr->p[iwfs]->p[1];
}

/**
   Accomplish Two tasks:
   
   1) High pass filter lgs focus to remove sodium range variation effact.
   2) Average LGS focus measurement to drive the trombone.
   
   We trust the focus measurement of the LGS WFS at high temporal frequency
   where NGS cannot provide due to low frame rate. After the HPF on lgs
   gradients, our system is NO LONGER affected by sodium layer variation.

   if sim.mffocus==1: The HPF is applied to each LGS WFS indepently. This largely
   removes the effect of differential focus. powfs.dfrs is no longer needed. (preferred).

   if sim.mffocus==2: We apply a LPF on the average focus from six LGS WFS, and
   then remove this value from all LGS WFS. The differential focus is still
   present and powfs.dfrs need to be set to 1 to handle it in tomography. This
   is the original focus tracking method.
*/
void wfsgrad_mffocus(SIM_T* simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    
    /*New plate mode focus offset for LGS WFS. Not needed*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].llt && parms->sim.ahstfocus==2 
	   && simu->Mint_lo && simu->Mint_lo->mint->p[1]
	   && (simu->isim+1)%parms->powfs[ipowfs].dtrat==0){
	    /*In new ahst mode, the first plate scale mode contains focus for
	      lgs. But it turns out to be not necessary to remove it because the
	      HPF in the LGS path removed the influence of this focus mode. set
	      sim.ahstfocus=2 to enable adjust gradients.*/
	    double scale=simu->recon->ngsmod->scale;
	    double focus=-simu->Mint_lo->mint->p[1]->p[0]->p[2]*(scale-1);
	    dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], focus);
	}
    }

    int hi_output=(!parms->sim.closeloop || (simu->isim+1)%parms->sim.dtrat_hi==0);
    if(hi_output){
	dcell *LGSfocus=NULL;/*residual focus along ngs estimated from LGS measurement.*/
	dcellmm(&LGSfocus, recon->RFlgsg, simu->gradcl,"nn",1);
	long nwfsllt=0; 
	double lpfocusm=0;
	double lgsfocusm=0;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(!LGSfocus->p[iwfs]) continue;
	    int ipowfs=parms->wfs[iwfs].powfs;
	    
	    //LPF can be put after using the value to put it off critical path. Put here to simulate case where lpfocus=1
	    double lpfocus=parms->sim.lpfocus;
	    simu->lgsfocuslpf->p[iwfs]=simu->lgsfocuslpf->p[iwfs]*(1-lpfocus)+LGSfocus->p[iwfs]->p[0]*lpfocus;

	    if(parms->sim.mffocus==1){//remove LPF focus from each lgs
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], -simu->lgsfocuslpf->p[iwfs]);
	    }
	    lgsfocusm+=LGSfocus->p[iwfs]->p[0];
	    //Averaged LPF focus
	    lpfocusm+=simu->lgsfocuslpf->p[iwfs]; 
	    nwfsllt++;
	    
	}
	lgsfocusm/=nwfsllt;
	if(parms->sim.mffocus==2){//remove LPF GLOBAL focus from each lgs
	    lpfocusm/=nwfsllt;
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		if(!LGSfocus->p[iwfs]) continue;
		int ipowfs=parms->wfs[iwfs].powfs;
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[ipowfs], -lpfocusm);
	    }
	}
    
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(LGSfocus->p[iwfs]){
		if(parms->sim.zoomshare){
		    //all lgs share the same trombone so take the average value.
		    simu->zoomavg->p[iwfs]+=lgsfocusm;
		}else{
		    simu->zoomavg->p[iwfs]+=LGSfocus->p[iwfs]->p[0];
		}
	    }
	}
	dcellfree(LGSfocus);
	/*zoom error is zero order hold even if no output from averager*/
	if((simu->reconisim+1)%parms->sim.zoomdtrat==0){
	    int dtrat=parms->sim.zoomdtrat;
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].llt &&
		   (!parms->powfs[ipowfs].dither || parms->powfs[ipowfs].phytype!=1)){
		    //For those WFS that dither and run mtch, use focus error from ib instead
		    simu->zoomerr->p[iwfs]=simu->zoomavg->p[iwfs]/(dtrat*dtrat);
		}
	    }
	    dzero(simu->zoomavg);
	}
    }
}

/**
   Every operation here should be in the Simulator not the Controller 
*/
void wfsgrad_post(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
#if USE_CUDA
    if(parms->gpu.wfs){
	gpu_wfsgrad_sync(info);
    }
#endif
    //Postprocessing gradients
    const POWFS_T *powfs=simu->powfs;
    const int isim=simu->isim;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const int dtrat_output=((isim+1-parms->powfs[ipowfs].step)%dtrat==0);
	const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
	dmat **gradout=&simu->gradcl->p[iwfs];
	if(dtrat_output){
	    if(do_phy){
		if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		    wfsgrad_fsm(simu, iwfs);
		}
		if(parms->powfs[ipowfs].phytypesim==2 && powfs[ipowfs].gradoffcog){
		    /*Gradient offset due to CoG*/
		    dadd(gradout, 1, powfs[ipowfs].gradoffcog->p[wfsind], -1);
		}
	    }
	    //Gradient offset due to mainly NCPA calibration
	    if(simu->powfs[ipowfs].gradoff){
		dadd(gradout, 1, simu->powfs[ipowfs].gradoff->p[wfsind], -1);
	    }
	    if(parms->save.grad->p[iwfs]){
		cellarr_dmat(simu->save->gradcl[iwfs], isim, simu->gradcl->p[iwfs]);
	    }
	    if(parms->plot.run){
		drawopd("Gclx",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p, NULL,
			"WFS Closeloop Gradients (x)","x (m)", "y (m)",
			"x %d",  iwfs);
		drawopd("Gcly",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p+powfs[ipowfs].pts->nsa, NULL,
			"WFS Closeloop Gradients (y)","x (m)", "y (m)",
			"y %d",  iwfs);
	    }
	}
    }//for iwfs
}
/**
   Dither update: zoom corrector, matched filter, TWFS
*/
static void dither_update(SIM_T *simu){
    POWFS_T *powfs=simu->powfs;
    const PARMS_T *parms=simu->parms;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].dither) continue;
	if((simu->isim+1)%parms->powfs[ipowfs].dtrat!=0) continue;
	const int nacc=(simu->isim-parms->powfs[ipowfs].dither_nskip+1)/parms->powfs[ipowfs].dtrat;
	const int nwfs=parms->powfs[ipowfs].nwfs;
	const int ndrift=parms->powfs[ipowfs].dither_ndrift;
	
	if(nacc>0 && nacc % ndrift==0){//There is drift mode computation
	    if(parms->sim.zoomshare){//average focus error
		double sum=0;
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    sum+=simu->zoomerr->p[iwfs];
		}
		sum/=(double)nwfs;
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    simu->zoomerr->p[iwfs]=sum;
		}
	    }
	}
	int nmtch=parms->powfs[ipowfs].dither_nmtch;
	if(nacc>0 && nacc % nmtch==0){//This is matched filter update
	    info("Step %d: Update matched filter for powfs %d\n", simu->isim, ipowfs);
	    if(!powfs[ipowfs].intstat){
		powfs[ipowfs].intstat=calloc(1, sizeof(INTSTAT_T));
	    }
	    const int nsa=powfs[ipowfs].pts->nsa;
	    parms->powfs[ipowfs].radgx=0;
	    if(!powfs[ipowfs].intstat->i0 || powfs[ipowfs].intstat->i0->ny!=nwfs){
		dcellfree(powfs[ipowfs].intstat->i0);
		dcellfree(powfs[ipowfs].intstat->gx);
		dcellfree(powfs[ipowfs].intstat->gy);

		powfs[ipowfs].intstat->i0=dcellnew(nsa, nwfs);
		powfs[ipowfs].intstat->gx=dcellnew(nsa, nwfs);
		powfs[ipowfs].intstat->gy=dcellnew(nsa, nwfs);
	    
	    }
	    double scale1=(double)parms->powfs[ipowfs].dither_ndrift/(double)parms->powfs[ipowfs].dither_nmtch;
	    for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		//End of accumulation
		DITHER_T *pd=simu->dither[iwfs];
		//Scale the output due to accumulation
		for(int isa=0; isa<nsa; isa++){
		    dadd(powfs[ipowfs].intstat->i0->p+isa+jwfs*nsa, 0, pd->i0->p[isa], scale1);
		    dadd(powfs[ipowfs].intstat->gx->p+isa+jwfs*nsa, 0, pd->gx->p[isa], scale1);
		    dadd(powfs[ipowfs].intstat->gy->p+isa+jwfs*nsa, 0, pd->gy->p[isa], scale1);
		}
		dcellzero(pd->i0);
		dcellzero(pd->gx);
		dcellzero(pd->gy);
		if(powfs[ipowfs].gradoff){
		    //goff is what matched filter remembered during last cycle.
		    dadd(&powfs[ipowfs].gradoff->p[jwfs], 1, simu->dither[iwfs]->goff, -scale1);
		}
	    }
	    genmtch(parms, powfs, ipowfs);
	    if(parms->powfs[ipowfs].phytypesim!=1){
		warning("powfs%d: switch to matched filter\n", ipowfs);
		parms->powfs[ipowfs].phytypesim=1;//set to matched filter
	    }
	    if(parms->powfs[ipowfs].phystep>simu->isim+1){
		warning("powfs%d: switch to physical optics wfs\n", ipowfs);
		parms->powfs[ipowfs].phystep=simu->isim+1;
	    }
	    if(parms->sim.epdm->p[0]<=0){
		parms->sim.epdm->p[0]=0.5;
		warning("set ephi to 0.5\n");
	    }
#if USE_CUDA
	    if(parms->gpu.wfs){
		gpu_wfsgrad_update_mtche(parms, powfs);
	    }
#endif
	    if(parms->save.dither){
		writebin(powfs[ipowfs].intstat->i0, "powfs%d_i0_%d", ipowfs, simu->isim);
		writebin(powfs[ipowfs].intstat->gx, "powfs%d_gx_%d", ipowfs, simu->isim);
		writebin(powfs[ipowfs].intstat->gy, "powfs%d_gy_%d", ipowfs, simu->isim);
		writebin(powfs[ipowfs].intstat->mtche, "powfs%d_mtche_%d", ipowfs, simu->isim);
		writebin(powfs[ipowfs].intstat->i0sum, "powfs%d_i0sum_%d", ipowfs, simu->isim);
	    }
	}
	int itpowfs=parms->itpowfs;
	if(itpowfs!=-1){
	    int ntacc=(simu->isim-parms->powfs[itpowfs].step+1);
	    if(ntacc>0 && ntacc%parms->powfs[itpowfs].dtrat==0){
		info2("Step %d: TWFS has output\n", simu->isim);
		dcell *Rmod=0;
		//Build radial mode error using closed loop TWFS measurements from this time step.
		dcellmm(&Rmod, simu->recon->RRtwfs, simu->gradcl, "nn", 1);
		writebin(simu->gradcl->p[parms->nwfsr-1], "twfs_gcl_%d", simu->isim);
		writebin(Rmod, "twfs_rmod_%d", simu->isim);
		if(!powfs[ipowfs].gradoff){
		    powfs[ipowfs].gradoff=dcellnew(nwfs, 1);
		}
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    dmm(&powfs[ipowfs].gradoff->p[jwfs], 1, simu->recon->GRall->p[ipowfs], Rmod->p[0], "nn", -0.5);//temp: set gain to 0.5
		    if(parms->plot.run){
			const int nsa=powfs[ipowfs].pts->nsa;
			drawopd("Goffx",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p,NULL,
				"WFS Offset (x)","x (m)", "y (m)", "x %d",  iwfs);
			drawopd("Goffy",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p+nsa, NULL,
				"WFS Offset (y)","x (m)", "y (m)", "y %d",  iwfs);
		    }
		}
		writebin(powfs[ipowfs].gradoff, "powfs%d_gradoff_%d", ipowfs, simu->isim);
		dcellfree(Rmod);
	    }
	}
    }
}
/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
   It also includes operations on Gradients before tomography.
*/
void wfsgrad(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(parms->sim.idealfit || parms->sim.evlol) return;
    // call the task in parallel and wait for them to finish. It may be done in CPU or GPU.
    extern int PARALLEL;
    if(!PARALLEL || parms->tomo.ahst_idealngs || !parms->gpu.wfs){
	CALL_THREAD(simu->wfs_grad_pre, 0);
    }
    CALL_THREAD(simu->wfs_grad_post, 0);
    dither_update(simu);
    if(parms->sim.mffocus){
	//high pass filter lgs focus to remove sodium range variation effect
	wfsgrad_mffocus(simu);
    }
#if USE_CUDA
    if(parms->gpu.wfs){
	gpu_wfsgrad_save(simu);
    }else
#endif
    wfsgrad_save(simu);
    simu->tk_wfs=myclockd()-tk_start;
}
