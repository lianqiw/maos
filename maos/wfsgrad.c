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
#include "cn2est.h"
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
      simu->gradol is pseudo OL grad output.
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
    const dsp *(*GA)[nwfs] = (const dsp *(*)[nwfs]) recon->GA->p;
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
    const int dtrat_reset=(!CL||dtrat==1||isim%dtrat==0);
    const int dtrat_output=(!CL||dtrat==1||(isim+1)%dtrat==0);
    const int do_phy=(parms->powfs[ipowfs].usephy&&isim>=parms->powfs[ipowfs].phystep);
    const int do_geom=!do_phy || save_gradgeom;

    dmat **gradacc=NULL;
    dmat **gradout=&simu->gradcl->p[iwfs];
    dcell *ints=NULL;

    if(do_geom){ //output to gradacc
	gradacc=&simu->gradacc->p[iwfs];
	if(!*gradacc)
	    *gradacc=dnew(nsa*2,1);
	else if(dtrat_reset)
	    dzero(*gradacc);
    }
    dmat *opd=dnew(npix,1);
    /*
      Add surface
    */
    if(simu->surfwfs && simu->surfwfs->p[iwfs]){
	dadd(&opd,1, simu->surfwfs->p[iwfs],1);
    }
    /*
      Add NCPA to WFS as needed. Todo: merge surfwfs
      with ncpa. be careful about ncpa calibration.
    */
    if(powfs[ipowfs].ncpa){
	//info2("Adding NCPA to wfs %d\n",iwfs);
	dadd(&opd, 1, powfs[ipowfs].ncpa->p[wfsind], 1);
    }
    double *realamp=powfs[ipowfs].realamp[wfsind];
    /*
      Now begin ray tracing.
    */
    if(atm){
	for(int ips=0; ips<nps; ips++){
	    thread_t *wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
	    wfs_propdata->phiout=opd->p;
	    wfs_propdata->displacex1=-atm[ips]->vx*dt*isim;//frozen flow.
	    wfs_propdata->displacey1=-atm[ips]->vy*dt*isim;
	    CALL_THREAD(wfs_prop, nthread, 0);//have to wait to finish before another phase screen.
	}//ips
	//most expensive 0.10 per LGS for
    }
    if(simu->telws){//Wind shake
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir->p[0];
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
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
	dmat **dmwfs=simu->moao_r_wfs->p;
	if(dmwfs[iwfs]){
	    info("iwfs %d: Adding MOAO correction\n", iwfs);
	    /*No need to do mis registration here since the MOAO DM is
	      attached to close to the WFS.*/
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
    /*
      Add defocus to OPD if needed.
    */
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
    if(parms->powfs[ipowfs].hasllt && simu->focusint && simu->focusint->p[iwfs]){
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
    /*
      calculate and save full aperture PSF is requested. The wavelength are
      evl.psfwvl, but powfs.wvl since our purpose is to compare with science PSF.
     */
    if(parms->wfs[iwfs].psfmean && isim>=parms->evl.psfisim){
	dmat *opdcopy=NULL;
	if(parms->powfs[ipowfs].trs){
	    opdcopy=ddup(opd);
	    double ptt[3]={0,0,0};
	    int indmcc=powfs[ipowfs].nlocm>1?wfsind:0;
	    loc_calc_ptt(NULL, ptt, powfs[ipowfs].loc, powfs[ipowfs].ipcc->p[indmcc],
			 powfs[ipowfs].imcc->p[indmcc], powfs[ipowfs].amp->p, opdcopy->p);
	    loc_remove_ptt(opdcopy->p, ptt, powfs[ipowfs].loc);
	}else{
	    opdcopy=dref(opd);
	}
	ccell *psf2s=psfcomp(opdcopy, powfs[ipowfs].realamp[wfsind], 
			     powfs[ipowfs].embed, powfs[ipowfs].nembed, 
			     parms->evl.psfsize, parms->evl.nwvl, parms->evl.psfwvl);
	dfree(opdcopy);
	int nwvl=parms->evl.nwvl;
	PDCELL(simu->wfspsfmean, wfspsfmean);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&wfspsfmean[iwfs][iwvl], 1,
		    psf2s->p[iwvl], pow(powfs[ipowfs].sumamp->p[wfsind], -2));
	}
	ccellfree(psf2s);
    }
    /*
      Now begin Physical Optics Intensity calculations
    */
    if(do_phy){ //initialize ints is not already initlialized
	if(dtrat_reset){
	    dcellzero(simu->ints[iwfs]);
	}
	ints=simu->ints[iwfs];
    }
    if(do_geom){
	if(parms->powfs[ipowfs].gtype_sim==1){
	    //compute ztilt.
	    pts_ztilt((*gradacc)->p,powfs[ipowfs].pts,
		      powfs[ipowfs].saimcc[powfs[ipowfs].nimcc>1?wfsind:0], 
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
	psfoutcellarr=simu->wfspsfoutcellarr[iwfs];
	ztiltoutcellarr=simu->ztiltoutcellarr[iwfs];
    }
    dcell *pistatout=NULL;
    if(parms->powfs[ipowfs].pistatout
       &&isim>=parms->powfs[ipowfs].pistatstart){
	//assert(parms->powfs[ipowfs].dtrat==1);
	if(!simu->pistatout[iwfs]){
	    simu->pistatout[iwfs]
		=dcellnew(nsa,parms->powfs[ipowfs].nwvl);
	}
	pistatout=simu->pistatout[iwfs];
    }
    TIM(1);
    if(ints || psfout || pistatout){
	dmat *lltopd=NULL;
	if(powfs[ipowfs].llt){
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
		    prop_grid_pts(atm[ips],powfs[ipowfs].llt->pts,
				  lltopd->p,1,displacex,displacey,
				  scale, 1., 0, 0);
		}
	    }
	    if((simu->uptreal && simu->uptreal->p[iwfs]) ||pistatout){
		const int nx=lltopd->nx;
		const double dx=powfs[ipowfs].llt->pts->dx;
		const double ox=powfs[ipowfs].llt->pts->origx[0];
		const double oy=powfs[ipowfs].llt->pts->origy[0];
		double ttx;
		double tty;
		if(pistatout||parms->sim.uptideal){
		    warning("Remove tip/tilt in uplink ideally\n");
		    //remove tip/tilt completely
		    dmat *lltg=dnew(2,1);
		    pts_ztilt(lltg->p,powfs[ipowfs].llt->pts,
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
		//copy uptreal to output
		PDMAT(simu->uptcmds->p[iwfs], puptcmds);
		puptcmds[isim][0]=ttx;
		puptcmds[isim][1]=tty;
		double vty=0;
		PDMAT(lltopd, pp);
		for(int iy=0; iy<lltopd->ny; iy++){
		    vty=(oy+iy*dx)*tty;
		    for(int ix=0; ix<nx; ix++){
			pp[iy][ix]+=vty+(ox+ix*dx)*ttx;
		    }
		}
	    }
	    if(save_opd){
		cellarr_dmat(simu->save->wfslltopd[iwfs],lltopd);
	    }
	}
	dmat *gradref=NULL;
	if(pistatout){
	    if(!parms->powfs[ipowfs].pistatstc && do_geom){
		gradref=*gradacc;
		//info("Using ztilt to shift pistat\n");
	    }else{
		//info("Using fft to shift pistat\n");
	    }
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
    if(dtrat_output && do_phy){
	/*
	  In Physical optics mode, finish integration and
	  output gradients to grad.
	*/
	//matched filter
	//index into matched filter.
	//most useful for multillt
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
	if(powfs[ipowfs].bkgrnd){
	    if(powfs[ipowfs].bkgrnd->ny==1){
		bkgrnd2=powfs[ipowfs].bkgrnd->p;
	    }else{
		bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*wfsind;
	    }
	}
	    
	double *srot=NULL;
	if(powfs[ipowfs].srot){
	    const int illt=parms->powfs[ipowfs].llt->i[wfsind];
	    //this is in r/a coordinate, get angles
	    srot=powfs[ipowfs].srot->p[illt]->p;
	}
	//output directly to simu->gradcl. replace
	double gnf[2]={0,0};
	double gny[2]={0,0};
	const double rne=parms->powfs[ipowfs].rne;
	const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	const double siglev=parms->wfs[iwfs].siglevsim;//don't multiply to dtrat.
	dmat *gradnf=dnew(nsa*2,1);//save noise free gradients.
	double *pgradx=(*gradout)->p;
	double *pgrady=pgradx+nsa;
	double *pgradnfx=gradnf->p;
	double *pgradnfy=pgradnfx+nsa;
	dcellscale(ints, siglev);
	if(save_ints){
	    cellarr_dcell(simu->save->intsnf[iwfs], ints);
	}
	for(int isa=0; isa<nsa; isa++){
	    /*
	      TODO: Do something to remove negative pixels. shift image
	      or mask out. This is important when bkgrndfnc is greater
	      than 1.
	    */
	    double pmax;
	    gnf[0]=0; gnf[1]=0;
	    switch(parms->powfs[ipowfs].phytypesim){
	    case 1:
		dmulvec(gnf, mtche[isa], ints->p[isa]->p,1.);
		break;
	    case 2:
		pmax=dmax(ints->p[isa]);
		dcog(gnf,ints->p[isa],0.,0.,0.1*pmax,0.1*pmax);
		gnf[0]*=pixtheta;
		gnf[1]*=pixtheta;
		break;
	    default:
		error("Invalid");
	    }
	    if(noisy){//add noise
		double *bkgrnd2i=(bkgrnd2 && bkgrnd2[isa])?bkgrnd2[isa]->p:NULL;
		
		addnoise(ints->p[isa], &simu->wfs_rand[iwfs],
			 bkgrnd,parms->powfs[ipowfs].bkgrndc,
			 bkgrnd2i, parms->powfs[ipowfs].bkgrndfnc, rne);
		gny[0]=0; gny[1]=0;
		switch(parms->powfs[ipowfs].phytypesim){
		case 1:
		    dmulvec(gny, mtche[isa],ints->p[isa]->p,1.);
		    if(parms->powfs[ipowfs].mtchscl){
			double scale=i0sum[isa]/dsum(ints->p[isa]);
			//info("scale wfs %d, isa %d by %g\n",iwfs,isa,scale);
			gny[0]*=scale;
			gny[1]*=scale;
		    }
		    break;
		case 2:
		    pmax=dmax(ints->p[isa]);
		    dcog(gny,ints->p[isa],0.,0.,0.1*pmax,0.1*pmax);
		    gny[0]*=pixtheta;
		    gny[1]*=pixtheta;
		    break;
		default:
		    error("Invalid");
		}
		simu->sanea_sim->p[iwfs]->p[isa]+=pow(gny[0]-gnf[0],2);
		simu->sanea_sim->p[iwfs]->p[isa+nsa]+=pow(gny[1]-gnf[1],2);
	    }else{
		gny[0]=gnf[0];
		gny[1]=gnf[1];
	    }
	    if(parms->powfs[ipowfs].radpix){
		/*
		  rotate gradients from ra to xy 
		  coordinate rotates cw by theta from ra to xy.
		  therefore vector rotates ccw by theta.
		  use standard R(theta)=[cos -sin; sin cos]
		  |xnew|                 |x|
		  |ynew|   =  R(theta) * |y|
		*/
		double theta;
		theta=srot[isa]; 
		    
		const double sth=sin(theta);
		const double cth=cos(theta);
		pgradnfx[isa]=gnf[0]*cth-gnf[1]*sth;
		pgradnfy[isa]=gnf[0]*sth+gnf[1]*cth;
		pgradx[isa]=gny[0]*cth-gny[1]*sth;
		pgrady[isa]=gny[0]*sth+gny[1]*cth;
	    }else{
		//already xy
		pgradnfx[isa]=gnf[0];
		pgradnfy[isa]=gnf[1];
		pgradx[isa]=gny[0];
		pgrady[isa]=gny[1];
	    }
	};//isa

	if(save_ints && ints){
	    cellarr_dcell(simu->save->intsny[iwfs], ints);
	}
	if(save_grad && noisy){
	    cellarr_dmat(simu->save->gradnf[iwfs], gradnf);
	}
	dfree(gradnf);
	if(parms->powfs[ipowfs].llt){
	    if(!recon->PTT){
		error("powfs %d has llt, but recon->PTT is NULL",ipowfs);
	    }
	    dmat *PTT=NULL;
	    if(parms->sim.recon==2){
		PTT=recon->PTT->p[ipowfs+ipowfs*parms->npowfs];
	    }else{
		PTT=recon->PTT->p[iwfs+iwfs*nwfs];
	    }
	    if(!PTT){
		error("powfs %d has llt, but TT removal is empty\n", ipowfs);
	    }
	    /*Compute LGS Uplink error*/
	    dzero(simu->upterr->p[iwfs]);
	    dmm(&simu->upterr->p[iwfs], PTT, simu->gradcl->p[iwfs], "nn", 1);
	    //copy upterr to output.
	    PDMAT(simu->upterrs->p[iwfs], pupterrs);
	    pupterrs[isim][0]=simu->upterr->p[iwfs]->p[0];
	    pupterrs[isim][1]=simu->upterr->p[iwfs]->p[1];
	
	}
    }else if (dtrat_output && !do_phy){
	//geomtric optics accumulation mode. scale and copy results to output.
	dcp(gradout,*gradacc);
	if(dtrat!=1)
	    dscale(*gradout,1./dtrat);//average
	if(noisy){
	    if(save_grad){//save noise free gradient.
		cellarr_dmat(simu->save->gradnf[iwfs], *gradout);
	    }
	    if(parms->powfs[ipowfs].usephy){
		info2("Will not add noise at acquisition proccess for physical optics\n");
	    }else{
		const dmat *nea=powfs[ipowfs].neasim->p[wfsind];
		const double *neax=nea->p;
		const double *neay=nea->p+nsa;
		double *ggx=(*gradout)->p;
		double *ggy=(*gradout)->p+nsa;
		for(int isa=0; isa<nsa; isa++){
		    //Preserve the random sequence.
		    double noisex=neax[isa]*randn(&simu->wfs_rand[iwfs]);
		    double noisey=neay[isa]*randn(&simu->wfs_rand[iwfs]);
		    ggx[isa]+=noisex;
		    ggy[isa]+=noisey;
		    simu->sanea_sim->p[iwfs]->p[isa]+=noisex*noisex;
		    simu->sanea_sim->p[iwfs]->p[isa+nsa]+=noisey*noisey;
		}
	    }
	}
    }
    if(dtrat_output && save_gradgeom){
	dmat *gradtmp=NULL;
	dcp(&gradtmp,*gradacc);
	if(dtrat!=1){
	    dscale(gradtmp,1./dtrat);
	}
	cellarr_dmat(simu->save->gradgeom[iwfs], gradtmp);//noise free.
    }

    if(parms->plot.run){
	drawopdamp("wfsopd",powfs[ipowfs].loc,opd->p,realamp,NULL,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
	drawopd("Gclx",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p, NULL,
		"WFS Closeloop Gradients (x)","x (m)", "y (m)",
		"x %d",  iwfs);
	drawopd("Gcly",(loc_t*)powfs[ipowfs].pts, simu->gradcl->p[iwfs]->p+nsa, NULL,
		"WFS Closeloop Gradients (y)","x (m)", "y (m)",
		"y %d",  iwfs);
    }
    dfree(opd);
    if(dtrat_output && powfs[ipowfs].ncpa_grad){
	warning("Applying ncpa_grad to gradout\n");
	dadd(gradout, 1, powfs[ipowfs].ncpa_grad->p[wfsind], -1);
    }
    //create pseudo open loop gradients. in split mode 1, only do for high order wfs.
    if(dtrat_output && parms->powfs[ipowfs].psol){
	dcp(&simu->gradol->p[iwfs], *gradout);
	if(CL){
	    if(simu->dmpsol[ipowfs]){
		for(int idm=0; idm<parms->ndm; idm++){
		    spmulmat(&simu->gradol->p[iwfs], GA[idm][iwfs],
			     simu->dmpsol[ipowfs]->p[idm], 1.);
		}
	    }
	}
	if(parms->plot.run){
	    drawopd("Gpolx",(loc_t*)powfs[ipowfs].pts, simu->gradol->p[iwfs]->p,NULL,
		    "WFS Pseudo Openloop Gradients (x)","x (m)", "y (m)",
		    "x %d",  iwfs);
	    drawopd("Gpoly",(loc_t*)powfs[ipowfs].pts, simu->gradol->p[iwfs]->p+nsa, NULL,
		    "WFS Pseudo Openloop Gradients (y)","x (m)", "y (m)",
		    "y %d",  iwfs);
	}
    }
    if(save_grad && dtrat_output){
	cellarr_dmat(simu->save->gradcl[iwfs], simu->gradcl->p[iwfs]);
	if(simu->gradol->p[iwfs]){
	    cellarr_dmat(simu->save->gradol[iwfs], simu->gradol->p[iwfs]);
	}
    }
    if(parms->cn2.pair && recon->cn2est->wfscov[iwfs]){
	cn2est_embed(recon->cn2est, simu->gradol->p[iwfs], iwfs);
    }
    TIM(3);
#if TIMING==1
    info("wfs %d grad timing: ray %.2f ints %.2f grad %.2f\n",iwfs,tk1-tk0,tk2-tk1,tk3-tk2);
#endif
}

/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
*/
void wfsgrad(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    if(parms->dbg.fitonly) return;
    //Updating dmpsol averaging.
    if(parms->sim.closeloop){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->tomo.split==2 || !parms->powfs[ipowfs].skip){
		const int dtrat=parms->powfs[ipowfs].dtrat;
		if(dtrat==1 || simu->isim%dtrat==0){
		    dcellfree(simu->dmpsol[ipowfs]);
		    //info("Freeing dmpsol for powfs %d\n", ipowfs);
		}
		dcell *dmpsol;
		if(parms->dbg.psol){
		    if(parms->sim.fuseint){
			dmpsol=simu->dmint[0];
		    }else{
			dmpsol=simu->dmint_hi[0];
		    }
		    //warning("Using dmreal next step for psol\n");
		}else{//add DM command for the same time step.
		    dmpsol=simu->dmreal;
		}
		//info("Accumulating dmpsol for powfs %d\n", ipowfs);
		dcelladd(&simu->dmpsol[ipowfs], 1, dmpsol, 1./parms->powfs[ipowfs].dtrat);
	    }
	}
    }
    //call the task in parallel and wait for them to finish.
    CALL_THREAD(simu->wfs_grad, parms->nwfs, 0);
    /*
      Uplink pointing servo. Moved to here from filter.c because of
      synchronization issue. dcellcp before integrator changes because wfsgrad
      updates upterr with current gradient.
    */
    dcellcp(&simu->uptreal, simu->uptint[0]);
    if(simu->upterr){
	//uplink tip/tilt mirror. use Integrator/Derivative control
	//update command for next step.
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

    if(parms->save.ngcov>0){
	//Outputing psol gradient covariance.
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    int iwfs1=parms->save.gcov[igcov*2];
	    int iwfs2=parms->save.gcov[igcov*2+1];
	    info("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
	    dmm(&simu->gcov->p[igcov], simu->gradol->p[iwfs1], simu->gradol->p[iwfs2],"nt",1);
	}
    }
    if(parms->cn2.pair){
	CN2EST_T *cn2est=recon->cn2est;
	cn2est_cov(cn2est);//convert gradients to cross covariance.
	if((simu->isim+1-parms->sim.start)%parms->cn2.step == 0){
	    dcellswrite(cn2est->cc, 1./cn2est->nstep, "cc_%d",simu->isim+1);
	    cn2est_est(cn2est, parms);//do the CN2 estimation
	    if(parms->cn2.moveht){
		cn2est_moveht(recon);
	    }
	    if(parms->cn2.tomo){
		cn2est_updatetomo(recon,parms);
	    }
	}
    }//if cn2est
    simu->tk_wfs=myclockd()-tk_start;
}
