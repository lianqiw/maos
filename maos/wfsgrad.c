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
#include "maos.h"
#include "sim.h"
#include "sim_utils.h"
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

static void wfs_ideal_correction(SIM_T *simu, dmat *opd, int iwfs, double alpha){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    RECON_T *recon=simu->recon;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const double hs=parms->powfs[ipowfs].hs;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
   
    for(int idm=0; idm<parms->ndm; idm++){
	const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	double dispx=ht*parms->wfs[iwfs].thetax+powfs[ipowfs].misreg[wfsind][0];
	double dispy=ht*parms->wfs[iwfs].thetay+powfs[ipowfs].misreg[wfsind][1];
	double scale=1.-ht/hs;
	if(parms->dm[idm].cubic){
	    prop_nongrid_cubic(recon->aloc[idm], simu->dmproj->p[idm]->p,
			       powfs[ipowfs].loc, powfs[ipowfs].amp->p, opd->p, 
			       alpha, dispx, dispy, scale, parms->dm[idm].iac, 
			       0, 0);
	}else{
	    prop_nongrid(recon->aloc[idm], simu->dmproj->p[idm]->p,
			 powfs[ipowfs].loc, powfs[ipowfs].amp->p, opd->p, 
			 alpha, dispx, dispy, scale, 
			 0, 0);
	}
    }
}

/**
   computes close loop and pseudo open loop gradidents for both gometric and
   physical optics WFS. Calls wfsints() to accumulate WFS subapertures images in
   physical optics mode.  */

void wfsgrad_iwfs(thread_t *info){
    SIM_T *simu=info->data;
    const PARMS_T *parms=simu->parms;
    const int iwfs=info->start;
    assert(info->end==info->start+1);//only 1 WFS.
    assert(iwfs<parms->nwfs);
    /*
      simu->gradcl is CL grad output
      simu->gradacc is internal, to accumulate geometric grads.
      do not accumulate opd. accumate ints for phy, g for GS
    */
    //input
    
    map_t **atm=simu->atm;
    const RECON_T *recon=simu->recon;
    const POWFS_T *powfs=simu->powfs;
    //output
    const int CL=parms->sim.closeloop;
    const int isim=simu->isim;
    const int nwfs=parms->nwfs;
    const int nps=parms->atm.nps;
    const double dt=simu->dt;
    TIM(0);
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
    //const int dtrat_reset=(isim%dtrat==0);
    const int dtrat_output=((isim+1)%dtrat==0);
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_geom=!do_phy || save_gradgeom;
    const double *realamp=powfs[ipowfs].realamp[wfsind];
    dmat **gradacc=&simu->gradacc->p[iwfs];
    dmat **gradout=&simu->gradcl->p[iwfs];
    dcell *ints=simu->ints[iwfs];
    dmat  *opd=dnew(npix,1);

    if(simu->telws){//Wind shake
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir?simu->winddir->p[0]:0;
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
    }

    /* Add surface error*/
    if(simu->surfwfs && simu->surfwfs->p[iwfs]){
	dadd(&opd,1, simu->surfwfs->p[iwfs],1);
    }

    /* Add NCPA to WFS as needed. Todo: merge surfwfs with ncpa. be careful
       about ncpa calibration. */
    if(powfs[ipowfs].ncpa){
	dadd(&opd, 1, powfs[ipowfs].ncpa->p[wfsind], 1);
    }

    /* Now begin ray tracing. */
    if(parms->sim.idealwfs){
	wfs_ideal_correction(simu, opd, iwfs, 1);
    }else if(atm){
	for(int ips=0; ips<nps; ips++){
	    thread_t *wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
	    wfs_propdata->phiout=opd->p;
	    wfs_propdata->displacex1=-atm[ips]->vx*dt*isim;
	    wfs_propdata->displacey1=-atm[ips]->vy*dt*isim;
	    /* have to wait to finish before another phase screen. */
	    CALL_THREAD(wfs_prop, nthread, 0);
	}//ips
	/* most expensive 0.10 per LGS for*/
	if(parms->sim.wfsalias){
	    /* Remove subspace of atm projected onto range of DM.*/
	    wfs_ideal_correction(simu, opd, iwfs,-1);
	}
    }
    if(save_opd){
	cellarr_dmat(simu->save->wfsopdol[iwfs], opd);
    }
 
    if(CL){
	for(int idm=0; idm<parms->ndm; idm++){
	    thread_t *wfs_prop=simu->wfs_prop_dm[iwfs+parms->nwfs*idm];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
	    wfs_propdata->phiout=opd->p;
	    CALL_THREAD(wfs_prop, nthread, 0);
	}//idm
    }
 
    if(imoao>-1){
	dmat **dmwfs=simu->moao_wfs->p;
	if(dmwfs[iwfs]){
	    info("iwfs %d: Adding MOAO correction\n", iwfs);
	    /* No need to do mis registration here since the MOAO DM is attached
	       to close to the WFS.*/
	    if(parms->moao[imoao].cubic){
		prop_nongrid_pts_cubic(recon->moao[imoao].aloc, dmwfs[iwfs]->p,
				       powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				       parms->moao[imoao].iac, 0, 0);
	    }else{
		prop_nongrid_pts(recon->moao[imoao].aloc, dmwfs[iwfs]->p,
				 powfs[ipowfs].pts, realamp, opd->p, -1, 0, 0, 1, 
				 0, 0);
	    }
	}
    }
    /* Add defocus to OPD if needed. */
    double focus=0;
    if(powfs[ipowfs].focus){
	int iy=0;
	int nx=powfs[ipowfs].focus->nx;
	int ny=powfs[ipowfs].focus->ny;
	if(ny==1){
	    iy=0;
	}else if(ny==parms->powfs[ipowfs].nwfs){
	    iy=wfsind;
	}else{
	    error("powfs[%d].focus wrong format\n",ipowfs);
	}
	int ix=isim%nx;
	double focusadd=powfs[ipowfs].focus->p[ix+nx*iy];
	if(fabs(focusadd)>1.e-200){
	    info("WFS %d: adding %g focus from input.\n", iwfs, focusadd);
	    focus+=focusadd;
	}
    }
    if(parms->powfs[ipowfs].llt && simu->focusint && simu->focusint->p[iwfs]){
	info("WFS %d: Adding focus adjust to %g\n", 
	     iwfs, simu->focusint->p[iwfs]->p[0]);
	focus+=simu->focusint->p[iwfs]->p[0];
    }
    if(fabs(focus)>1.e-200){
	loc_add_focus(opd->p, powfs[ipowfs].loc, focus);
    }
    if(save_opd){
	cellarr_dmat(simu->save->wfsopd[iwfs], opd);
    }

    if(do_geom){
	/* Now Geometric Optics gradient calculations */
	if(parms->powfs[ipowfs].gtype_sim==1){
	    //compute ztilt.
	    pts_ztilt(gradacc,powfs[ipowfs].pts,
		      powfs[ipowfs].saimcc[powfs[ipowfs].nsaimcc>1?wfsind:0], 
		      realamp, opd->p);
	}else{//G tilt
	    spmulmat(gradacc,adpind(powfs[ipowfs].GS0,wfsind),opd,1);
	}
    }

    ccell *psfout=NULL;
    cellarr *psfoutcellarr=NULL;
    cellarr *ztiltoutcellarr=NULL;
    if(parms->powfs[ipowfs].psfout){
	psfout=simu->wfspsfout[iwfs];
	psfoutcellarr=simu->save->wfspsfout[iwfs];
	ztiltoutcellarr=simu->save->ztiltout[iwfs];
    }
    dcell *pistatout=NULL;
    if(parms->powfs[ipowfs].pistatout
       &&isim>=parms->powfs[ipowfs].pistatstart){
	if(!simu->pistatout[iwfs]){
	    simu->pistatout[iwfs]
		=dcellnew(nsa,parms->powfs[ipowfs].nwvl);
	}
	pistatout=simu->pistatout[iwfs];
    }
    TIM(1);
    /* Now begin Physical Optics Intensity calculations */
    if(do_phy || psfout || pistatout){
	dmat *lltopd=NULL;
	if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
	    if(powfs[ipowfs].llt->ncpa){
		int iotf=powfs[ipowfs].llt->ncpa->nx==1?0:wfsind;
		lltopd=ddup(powfs[ipowfs].llt->ncpa->p[iotf]);
	    }else{
		lltopd=dnew(powfs[ipowfs].llt->pts->nx,
			    powfs[ipowfs].llt->pts->nx);
	    }
	    const int illt=parms->powfs[ipowfs].llt->i[wfsind];
	    if(atm){//LLT OPD
		for(int ips=0; ips<nps; ips++){
		    const double hl=atm[ips]->h;
		    const double scale=1.-hl/hs;
		    /*
		      Bug fixed: 2009-01-05: multiply ox,oy to hl/hs.
		      Bug fixed: 2011-02-14: -ox/hs instead of ox/hs.
		    */
		    const double thetax=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox[illt]/hs;
		    const double thetay=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy[illt]/hs;
		    const double displacex=-atm[ips]->vx*isim*dt+thetax*hl+parms->powfs[ipowfs].llt->misreg[0];
		    const double displacey=-atm[ips]->vy*isim*dt+thetay*hl+parms->powfs[ipowfs].llt->misreg[1];
		    prop_grid_pts(atm[ips],powfs[ipowfs].llt->pts,NULL,
				  lltopd->p,1,displacex,displacey,
				  scale, 1., 0, 0);
		}
	    }
	    if((simu->uptreal && simu->uptreal->p[iwfs]) ||pistatout||parms->sim.uptideal){
		const double dx=powfs[ipowfs].llt->pts->dx;
		const double ox=powfs[ipowfs].llt->pts->origx[0];
		const double oy=powfs[ipowfs].llt->pts->origy[0];
		double ttx;
		double tty;
		if(pistatout||parms->sim.uptideal){
		    warning("Remove tip/tilt in uplink ideally\n");
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
		/* add tip/tilt to opd */
		PDMAT(lltopd, pp);
		for(int iy=0; iy<lltopd->ny; iy++){
		    double vty=(oy+iy*dx)*tty;
		    for(int ix=0; ix<lltopd->nx; ix++){
			pp[iy][ix]+=vty+(ox+ix*dx)*ttx;
		    }
		}
	    }
	    if(save_opd){
		cellarr_dmat(simu->save->wfslltopd[iwfs],lltopd);
	    }
	}
	dmat *gradref=NULL;
	if(pistatout && !parms->powfs[ipowfs].pistatstc && do_geom){
	    gradref=*gradacc;
	}
	WFSINTS_T *intsdata=simu->wfs_intsdata+iwfs;
	intsdata->ints=ints;
	intsdata->psfout=psfout;
	intsdata->pistatout=pistatout;
	intsdata->gradref=gradref;
	intsdata->opd=opd;
	intsdata->lltopd=lltopd;
	CALL_THREAD(simu->wfs_ints[iwfs], nthread, 0);
	dfree(lltopd);
	if(psfout){
	    cellarr_ccell(psfoutcellarr,psfout);
	    cellarr_dmat(ztiltoutcellarr, *gradacc);
	}
    }
    TIM(2);
 
    if(parms->plot.run){
	drawopdamp("wfsopd",powfs[ipowfs].loc,opd->p,realamp,NULL,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
    }
    dfree(opd);

    if(dtrat_output){
	if(do_phy){
	    /* In Physical optics mode, do integration and compute
	       gradients. The matched filter are for x/y directions even if
	       radpix=1. */
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
	    double pixtheta=parms->powfs[ipowfs].pixtheta;
	    
	    //Rayleigh scattering (bkgrnd)
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
	    
	    //output directly to simu->gradcl. replace
	    const double rne=parms->powfs[ipowfs].rne;
	    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	    const double siglev=parms->wfs[iwfs].siglevsim;//don't multiply to dtrat.
	    double *pgradx=(*gradout)->p;
	    double *pgrady=pgradx+nsa;
	    dmat *gradnf=NULL;
	    if(save_grad){
		gradnf=dnew(nsa*2,1);//save noise free gradients.
	    }
	    dcellscale(ints, siglev);
	    if(save_ints){
		cellarr_dcell(simu->save->intsnf[iwfs], ints);
	    }
	    for(int isa=0; isa<nsa; isa++){
		/* TODO: Do something to remove negative pixels. shift image or
		   mask out. This is important when bkgrndfnc is greater than
		   1. */
		double gnf[2]={0,0};
		double gny[2]={0,0};
		switch(parms->powfs[ipowfs].phytypesim){
		case 1:
		    dmulvec(gnf, mtche[isa], ints->p[isa]->p,1.);
		    break;
		case 2:{
		    double pmax=dmax(ints->p[isa]);
		    dcog(gnf,ints->p[isa],0.,0.,0.1*pmax,0.1*pmax);
		    gnf[0]*=pixtheta;
		    gnf[1]*=pixtheta;
		}
		    break;
		default:
		    error("Invalid");
		}
		if(noisy){//add noise
		    double *bkgrnd2i=(bkgrnd2 && bkgrnd2[isa])?bkgrnd2[isa]->p:NULL;
		    double *bkgrnd2ic=(bkgrnd2c && bkgrnd2c[isa])?bkgrnd2c[isa]->p:NULL;
		    addnoise(ints->p[isa], &simu->wfs_rand[iwfs],
			     bkgrnd,parms->powfs[ipowfs].bkgrndc,
			     bkgrnd2i, bkgrnd2ic, rne);
		    switch(parms->powfs[ipowfs].phytypesim){
		    case 1:
			dmulvec(gny, mtche[isa],ints->p[isa]->p,1.);
			if(parms->powfs[ipowfs].mtchscl){
			    double scale=i0sum[isa]/dsum(ints->p[isa]);
			    gny[0]*=scale;
			    gny[1]*=scale;
			}
			break;
		    case 2:{
			double pmax=dmax(ints->p[isa]);
			dcog(gny,ints->p[isa],0.,0.,0.1*pmax,0.1*pmax);
			gny[0]*=pixtheta;
			gny[1]*=pixtheta;
		    }
			break;
		    default:
			error("Invalid");
		    }
		    double errx=gny[0]-gnf[0];
		    double erry=gny[1]-gnf[1];
		    simu->sanea_sim[iwfs]->p[isa]->p[0]+=errx*errx;
		    simu->sanea_sim[iwfs]->p[isa]->p[1]+=errx*erry;
		    simu->sanea_sim[iwfs]->p[isa]->p[2]+=errx*erry;
		    simu->sanea_sim[iwfs]->p[isa]->p[3]+=erry*erry;
		}else{
		    gny[0]=gnf[0];
		    gny[1]=gnf[1];
		}
		if(save_grad){
		    gradnf->p[isa]=gnf[0];
		    gradnf->p[isa+nsa]=gnf[1];
		}
		pgradx[isa]=gny[0];
		pgrady[isa]=gny[1];
	    };//isa

	    if(save_ints){
		cellarr_dcell(simu->save->intsny[iwfs], ints);
	    }
	    if(save_grad && noisy){
		cellarr_dmat(simu->save->gradnf[iwfs], gradnf);
		dfree(gradnf);
	    }
	    if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
		if(!recon->PTT){
		    error("powfs %d has llt, but recon->PTT is NULL",ipowfs);
		}
		dmat *PTT=NULL;
		if(parms->sim.glao){
		    PTT=recon->PTT->p[ipowfs+ipowfs*parms->npowfs];
		}else{
		    PTT=recon->PTT->p[iwfs+iwfs*nwfs];
		}
		if(!PTT){
		    error("powfs %d has llt, but TT removal is empty\n", ipowfs);
		}
		/* Compute LGS Uplink error. */
		dzero(simu->upterr->p[iwfs]);
		dmm(&simu->upterr->p[iwfs], PTT, *gradout, "nn", 1);
		/* copy upterr to output. */
		PDMAT(simu->upterrs->p[iwfs], pupterrs);
		pupterrs[isim][0]=simu->upterr->p[iwfs]->p[0];
		pupterrs[isim][1]=simu->upterr->p[iwfs]->p[1];
	    }
	}else{
	    /* geomtric optics accumulation mode. scale and copy results to
	       output. */
	    dcp(gradout,*gradacc);
	    if(dtrat!=1)
		dscale(*gradout,1./dtrat);//average
	    if(noisy){
		if(save_grad){//save noise free gradient.
		    cellarr_dmat(simu->save->gradnf[iwfs], *gradout);
		}
		if(!parms->powfs[ipowfs].usephy){
		    const dmat *nea=powfs[ipowfs].neasim->p[wfsind];
		    const double *neax=nea->p;
		    const double *neay=nea->p+nsa;
		    double *ggx=(*gradout)->p;
		    double *ggy=(*gradout)->p+nsa;
		    for(int isa=0; isa<nsa; isa++){
			//Preserve the random sequence.
			double errx=neax[isa]*randn(&simu->wfs_rand[iwfs]);
			double erry=neay[isa]*randn(&simu->wfs_rand[iwfs]);
			ggx[isa]+=errx;
			ggy[isa]+=erry;
			simu->sanea_sim[iwfs]->p[isa]->p[0]+=errx*errx;
			simu->sanea_sim[iwfs]->p[isa]->p[1]+=errx*erry;
			simu->sanea_sim[iwfs]->p[isa]->p[2]+=errx*erry;
			simu->sanea_sim[iwfs]->p[isa]->p[3]+=erry*erry;
		    }
		}else if(isim==0){
		    info2("Will not add noise at acquisition proccess for physical optics\n");
		}
	    }
	}
  
	if(powfs[ipowfs].ncpa_grad){
	    warning("Applying ncpa_grad to gradout\n");
	    dadd(gradout, 1, powfs[ipowfs].ncpa_grad->p[wfsind], -1);
	}
	if(save_grad){
	    cellarr_dmat(simu->save->gradcl[iwfs], simu->gradcl->p[iwfs]);
	}
	if(save_gradgeom){
	    dmat *gradtmp=NULL;
	    dadd(&gradtmp, 1, *gradacc, 1./dtrat);
	    cellarr_dmat(simu->save->gradgeom[iwfs], gradtmp);//noise free.
	    dfree(gradtmp);
	}
	if(parms->plot.run){
	    drawopd("Gclx",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p, NULL,
		    "WFS Closeloop Gradients (x)","x (m)", "y (m)",
		    "x %d",  iwfs);
	    drawopd("Gcly",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p+nsa, NULL,
		    "WFS Closeloop Gradients (y)","x (m)", "y (m)",
		    "y %d",  iwfs);
	}
	if(do_geom){
	    dzero(*gradacc);
	}
	if(do_phy){
	    dcellzero(simu->ints[iwfs]);
	}
    }
    
    TIM(3);
#if TIMING==1
    info("wfs %d grad timing: ray %.2f ints %.2f grad %.2f\n",iwfs,tk1-tk0,tk2-tk1,tk3-tk2);
#endif
}
/**
   Save telemetry
 */
static void wfsgrad_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    const int seed=simu->seed;
    if((isim % 50 ==0) || isim+1==parms->sim.end){
	if(simu->wfspsfmean && simu->isim>=parms->evl.psfisim){
	    double scalewfs=1./(double)(simu->isim+1-parms->evl.psfisim);
	    dcellswrite(simu->wfspsfmean, scalewfs, "wfspsfmean_%d.bin", seed);
	}
	for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	    if(!simu->sanea_sim[iwfs]) continue;
	    dcell *sanea=NULL;
	    dcellcp(&sanea, simu->sanea_sim[iwfs]);
	    const int ipowfs=simu->parms->wfs[iwfs].powfs;
	    const int dtrat=parms->powfs[ipowfs].dtrat;
	    if(sanea && simu->isim >=simu->parms->powfs[ipowfs].phystep){
		int nstep=(simu->isim+1-simu->parms->powfs[ipowfs].phystep)/dtrat;
		dcellscale(sanea,1./nstep);
	    }
	    dcellwrite(sanea,"sanea_sim_wfs%d_%d.bin",iwfs,seed);
	    dcellfree(sanea);
	}
	if(simu->pistatout){
	    for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
		const int ipowfs=simu->parms->wfs[iwfs].powfs;
		if(simu->pistatout[iwfs]){
		    int nstep=isim+1-parms->powfs[ipowfs].pistatstart;
		    dcell* tmp=NULL;
		    dcelladd(&tmp,0,simu->pistatout[iwfs],1./(double)nstep);
		    if(parms->sim.skysim){//need peak in corner
			for(long ic=0; ic<tmp->nx*tmp->ny; ic++){
			    dfftshift(tmp->p[ic]);
			}
			dcellwrite(tmp,"%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
				   dirskysim,simu->seed,
				   parms->powfs[ipowfs].order,
				   parms->wfs[iwfs].thetax*206265,
				   parms->wfs[iwfs].thetay*206265);
		    }else{//need peak in center
			dcellwrite(tmp,"pistat_seed%d_wfs%d.bin", simu->seed,iwfs);
		    }
		    dcellfree(tmp);
		}
	    }
	}
    }
}
/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
*/
void wfsgrad(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.idealfit || parms->sim.evlol) return;
    double tk_start=myclockd();
  
    /* call the task in parallel and wait for them to finish. */
    CALL_THREAD(simu->wfs_grad, parms->nwfs, 0);
    /* Uplink pointing servo. Moved to here from filter.c because of
       synchronization issue. dcellcp before integrator changes because wfsgrad
       updates upterr with current gradient. */
    dcellcp(&simu->uptreal, simu->uptint[0]);
    if(simu->upterr){
	/* uplink tip/tilt mirror. use Integrator/Derivative control
	   update command for next step.*/
	shift_inte(parms->sim.napupt, parms->sim.apupt, simu->uptint);
	double gain1=parms->sim.epupt+parms->sim.dpupt;
	double gain2=-parms->sim.dpupt;
	dcelladd(&simu->uptint[0], 1., simu->upterr,gain1);
	if(fabs(gain2)>EPS){
	    dcelladd(&simu->uptint[0], 1., simu->upterrlast,gain2);
	    //save to use in next servo time step.
	    dcellcp(&simu->upterrlast,simu->upterr);
	}
    }
    wfsgrad_save(simu);
    simu->tk_wfs=myclockd()-tk_start;
}
