/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "setup_surf.h"
#include "setup_powfs.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file maos/sim_utils.c
   Contains a few support functions for simulation.
*/
/*static real opdzlim[2]={-3e-5,3e-5}; */
static real *opdzlim=NULL;
extern int disable_save;
static mapcell *genatm_do(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const ATM_CFG_T *atm=&parms->atm;
    TIC;
    mapcell *screens;
    if(!parms->dbg.atm){
	GENATM_T *gs=simu->atmcfg;
	if(!gs){
	    simu->atmcfg=mycalloc(1,GENATM_T);/*the data for generating screens. */
	    gs=simu->atmcfg;
	    gs->rstat  = simu->atm_rand;
	    gs->wt     = atm->wt->p;
	    gs->r0     = atm->r0;
	    gs->L0     = atm->L0->p;
	    gs->dx     = atm->dx;
	    gs->nx     = atm->nx;
	    gs->ny     = atm->ny;
	    gs->nlayer = atm->nps;
	    gs->ninit  = atm->ninit;
	    gs->share  = atm->share;
	    if((parms->atm.r0evolve & 1) == 1){
		dbg("Scaling turbulence screen spatially\n");
		gs->r0logpsds=atm->r0logpsds;
	    }
	}
	info("Generating Atmospheric Screen...\n");
	tic;
	screens = parms->atm.fun(gs);
	toc("Atmosphere ");
    }else{
	/*
	  create screens on two layers that produce pure
	  tip/tilt for LGS to debug split tomography test
	  pass in open loop mode for both Z and G tilt for
	  NGS.
	      
	  The residual error is pure fit error if
	  layer 5 is at dm layer. otherwise the error is a
	  little larger.
	*/
	int nx=atm->nx;
	int ny=atm->ny;
	screens=mapcellnew(atm->nps, 1);
	real dx=atm->dx;
	int iratio=0;
	if(parms->dbg.atm->nx>=atm->nps){
	    iratio=1;
	}
	loc_t *psloc=0;
	const real strength=sqrt(1.0299*pow(parms->aper.d/atm->r0, 5./3.))*(0.5e-6/(2*M_PI));//PR WFE.
	for(int ips=0; ips<atm->nps; ips++){
	    screens->p[ips]=mapnew(nx, ny, dx, dx);
	    screens->p[ips]->h=atm->ht->p[ips];
	    real dbgatm=parms->dbg.atm->p[ips*iratio];
	    if(dbgatm>0){//zernike mode
		if(!psloc){
		    psloc=mksqloc_auto(nx, ny, atm->dx, atm->dx);
		}
		info("Generating Testing Atmosphere Screen with zernike %g, RMS~=%g nm\n", -dbgatm, strength*1e9);
		dmat *opd=zernike(psloc, nx*atm->dx, 0, 0, -dbgatm);
		dmat *opd2=dref_reshape(opd, nx, ny);
		dadd((dmat**)&screens->p[ips], 0, opd2, atm->wt->p[ips]*strength);
		dfree(opd);
		dfree(opd2);
	    }else if(dbgatm<0){//Fourier mode;
		info("Generating Testing Atmosphere Screen with Fourier mode %g, RMS~=%g nm\n", -dbgatm, strength*1e9);

		real kk=2*M_PI/dbgatm*dx;
		long nn=MAX(nx, ny);
		dmat *sind=dnew(nn, 1);
		for(int ii=0; ii<nn; ii++){
		    sind->p[ii]=sin((ii-nn/2)*kk);
		}
		const real alpha=2*strength*atm->wt->p[ips];
		for(int iy=0; iy<ny; iy++){
		    for(int ix=0; ix<nx; ix++){
			P(screens->p[ips], ix, iy)=alpha*sind->p[ix]*sind->p[iy];
		    }
		}
		dfree(sind);
	    }else{
		//empty screen.
	    }
	}
	locfree(psloc);
    }
 
    return screens;
}

/**
   overlay atm2 with atm2 according to wind direction angle and required
   overlapping region of at least overx*overy.
*/
static void blend_screen_side(map_t *atm1, map_t *atm2, long overx, long overy){
    const long nx=atm1->nx;
    const long ny=atm1->ny;
    int ca=0;
    if(atm1->vx>EPS){
	ca=-1;/*reverse sign of vx */
    }else if(atm1->vx<-EPS){
	ca=1;
    }
    int sa=0;
    if(atm1->vy>EPS){
	sa=-1;/*reverse sign of vy */
    }else if(atm1->vy<-EPS){
	sa=1;
    }
    long rr;
    long offx=nx-overx;
    long offy=(ny-overy)*nx;
    if(ca==0){/*along y. */
	rr=ny-overy;/*distance between the origins. */
	atm2->oy=atm1->oy + rr*sa*atm1->dx;
	atm2->ox=atm1->ox;
	real wty=sa<0?1:0;
	real *p1=atm1->p+(1-(long)wty)*offy;
	real *p2=atm2->p+(long)wty*offy;
	real (*pp1)[nx]=(real(*)[nx])p1;
	real (*pp2)[nx]=(real(*)[nx])p2;
	real overyd=(real)overy;
	for(long iy=0; iy<overy; iy++){
	    real wt1=fabs(wty-(real)(iy+1)/overyd);
	    for(long ix=0; ix<nx; ix++){
		pp1[iy][ix]=(1-wt1)*pp1[iy][ix]+wt1*pp2[iy][ix];
		pp2[iy][ix]=pp1[iy][ix];
	    }
	}
    }else if(sa==0){
	rr=nx-overx;/*distance between the origins. */
	atm2->ox=atm1->ox + rr*ca*atm1->dx;
	atm2->oy=atm1->oy;
	real wtx=ca<0?1:0;
	real *p1=atm1->p+(1-(long)wtx)*offx;
	real *p2=atm2->p+(long)wtx*offx;
	real (*pp1)[nx]=(real(*)[nx])p1;
	real (*pp2)[nx]=(real(*)[nx])p2;
	real wts[overx];
	real overxd=(real)overx;
	for(long ix=0; ix<overx; ix++){
	    wts[ix]=fabs(wtx-(real)(ix+1)/overxd);
	}
	for(long iy=0; iy<ny; iy++){
	    for(long ix=0; ix<overx; ix++){
		pp1[iy][ix]=(1-wts[ix])*pp1[iy][ix]+wts[ix]*pp2[iy][ix];
		pp2[iy][ix]=pp1[iy][ix];
	    }
	}
    }else{
	error("We do not support this wind direction: ca=%d, sa=%d\n", ca, sa);
    }
}
/**
   wrap of the generic vonkarman_genatm to generate turbulence screens. Wind
   velocities are set for each screen.  \callgraph */
void genatm(SIM_T *simu){ 
    const PARMS_T *parms=simu->parms;
    const ATM_CFG_T *atm=&(simu->parms->atm);
    if(simu->atm){
	cellfree(simu->atm);
	dfree(simu->winddir);
    }
    if(simu->parms->sim.noatm){
	warning("sim.noatm flag is on. will not generate atmoshere\n");
	return;
    }
    info("Wind dir:");/*initialize wind direction one time only for each seed in frozen flow mode. */
    simu->winddir=dnew(atm->nps,1);
    int wdnz=0;
    for(int i=0; i<atm->nps; i++){
	real angle;
	if(atm->wdrand){
	    if(fabs(atm->wddeg->p[i])>EPS){
		wdnz=1;
	    }
	    angle=randu(simu->atmwd_rand)*M_PI*2;
	}else{
	    angle=atm->wddeg->p[i]*M_PI/180;
	}
	simu->winddir->p[i]=angle;
	info(" %5.1f", angle*180/M_PI);
    }
    info(" deg\n");
    if(wdnz){
	error("wdrand is specified, but wddeg are not all zero. \n"
	      "possible confliction of intension!\n");
    }
    if(simu->parms->load.atm){
	const char *fn=simu->parms->load.atm;
	info("loading atm from %s\n",fn);
	simu->atm = mapcellread("%s",fn);
	if(simu->atm->nx!=atm->nps) error("ATM Mismatch\n");
    }else{
	simu->atm=genatm_do(simu);
    }
    if(!parms->dbg.atm){
	for(int i=0; i<atm->nps; i++){
	    real angle=simu->winddir->p[i];
	    simu->atm->p[i]->h=parms->atm.ht->p[i];
	    simu->atm->p[i]->vx=cos(angle)*parms->atm.ws->p[i];
	    simu->atm->p[i]->vy=sin(angle)*parms->atm.ws->p[i];
	}
    }
    if(simu->parms->save.atm){
	writebin(simu->atm,"atm_%d.bin",simu->seed);
    }
    
    if(parms->plot.atm && simu->atm){
	for(int ips=0; ips<atm->nps; ips++){
	    drawmap("atm", simu->atm->p[ips],opdzlim,
		    "Atmosphere OPD","x (m)","y (m)","layer%d",ips);
	}
    }
  
    info("After genatm:\t%.2f MiB\n",get_job_mem()/1024.);
    if((parms->atm.r0evolve&2)==2){
	dbg("Scaling OPD temporarily\n");
	simu->atmscale=psd2time(parms->atm.r0logpsdt, simu->atm_rand, parms->sim.dt, parms->sim.end);
	const real r02wt=(-5./3.);//layer weight is prop to r0^(-5/3)
	for(long i=0; i<simu->atmscale->nx; i++){
	    simu->atmscale->p[i]=exp((simu->atmscale->p[i])*r02wt);//convert to cn2dh
	}
	writebin(simu->atmscale, "atmscale_%d", simu->seed);
    }
    if(simu->wfs_prop_atm){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    for(int ips=0; ips<parms->atm.nps; ips++){
		PROPDATA_T *data=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
		data->mapin=simu->atm->p[ips];
	    }
	}
    }
    if(simu->evl_prop_atm){
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    for(int ips=0; ips<parms->atm.nps; ips++){
		PROPDATA_T *data=&simu->evl_propdata_atm[ievl+parms->evl.nevl*ips];
		data->mapin=simu->atm->p[ips];
	    }
	}
    }
}
/**
   Setup ray tracing operator from xloc to ploc, with predictive offsetting
*/
void setup_recon_HXW_predict(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    info("Generating Predictive HXW\n");
    dspcell* HXWtomo=recon->HXWtomo/*PDSPCELL*/;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){/*for tomography */
	    const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
	    const real hs=parms->wfs[iwfs].hs;
	    const real hc=parms->wfs[iwfs].hc;
	    for(int ips=0; ips<npsr; ips++){
		dspfree(P(HXWtomo,iwfs,ips));
		real  ht = recon->ht->p[ips];
		real  scale=1. - (ht-hc)/hs;
		real  displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht;
		displace[1]=parms->wfsr[iwfs].thetay*ht;
		if(parms->tomo.predict){
		    int ips0=parms->atmr.indps->p[ips];
		    displace[0]+=simu->atm->p[ips0]->vx*delay;
		    displace[1]+=simu->atm->p[ips0]->vy*delay;
		}
		P(HXWtomo,iwfs,ips)=mkh(recon->xloc->p[ips], ploc, 
				       displace[0],displace[1],scale);
	    }
	}
    }
}
/**
   Propagate the atmosphere to closest xloc. skip wavefront sensing and
   reconstruction. 

   2011-04-26: opdx was incorrectly computed when atm.ht and atmr.ht does not
   match in number. Fixed. Do not do scaling even if fit.ht is less.

*/
void atm2xloc(dcell **opdx, const SIM_T *simu){
    const RECON_T *recon=simu->recon;
    const PARMS_T *parms=simu->parms;
    if(parms->recon.glao){
	return;
    }
    /*in close loop mode, opdr is from last time step. */
    int isim=simu->reconisim;
    if(!*opdx){
	*opdx=dcellnew(recon->npsr,1);
    }
    for(int ipsr=0; ipsr<recon->npsr; ipsr++){
	if(!(*opdx)->p[ipsr]){
	    (*opdx)->p[ipsr]=dnew(recon->xloc->p[ipsr]->nloc,1);
	}else{
	    dzero((*opdx)->p[ipsr]);
	}
    }
    if(simu->atm){
	for(int ips=0; ips<parms->atm.nps; ips++){
	    real disx=-simu->atm->p[ips]->vx*isim*parms->sim.dt;
	    real disy=-simu->atm->p[ips]->vy*isim*parms->sim.dt;
	    int ipsr=parms->atm.ipsr->p[ips];
	    prop_grid(simu->atm->p[ips],recon->xloc->p[ipsr],(*opdx)->p[ipsr]->p,
		      1,disx,disy,1,1,0,0);
	}
    }
}

/**
   Evolving the Sodium layer by updating the elongation transfer function.
*/
void sim_update_etf(SIM_T *simu){
    int isim=simu->wfsisim;
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	/* Update ETF if necessary. */
	if((parms->powfs[ipowfs].usephy
	    ||parms->powfs[ipowfs].psfout
	    ||parms->powfs[ipowfs].pistatout) 
	   && parms->powfs[ipowfs].llt
	   && parms->powfs[ipowfs].llt->coldtrat>0
	   && isim %parms->powfs[ipowfs].llt->coldtrat == 0){
	    info("Step %d: powfs %d: Updating ETF\n",isim, ipowfs);
	    int dtrat=parms->powfs[ipowfs].llt->coldtrat;
	    int colsim=parms->powfs[ipowfs].llt->colsim;
	    setup_powfs_etf(powfs,parms,ipowfs,1, colsim+isim/dtrat);
	    setup_powfs_etf(powfs,parms,ipowfs,2, colsim+isim/dtrat+1);
#if USE_CUDA
	    if(parms->gpu.wfs){
		gpu_wfsgrad_update_etf(parms, powfs);
	    }
#endif
	}
    }
}

/**
   use random number dirived from input seed to seed other stream.  necessary to
   have independant streams for different wfs in threading routines to avoid
   race condition and have consitent result */
void seeding(SIM_T *simu){
    info("Running seed %d\n",simu->seed);
    simu->init_rand=mycalloc(1,rand_t);
    simu->atm_rand=mycalloc(1,rand_t);
    simu->atmwd_rand=mycalloc(1,rand_t);
    simu->telws_rand=mycalloc(1,rand_t);
    simu->misc_rand=mycalloc(1,rand_t);
    seed_rand(simu->init_rand,simu->seed);
    seed_rand(simu->atm_rand,   lrand(simu->init_rand));
    /*2011-02-02: changed to wdrand-1 so that when wdrand=1, we reproduce old directions. */
    seed_rand(simu->atmwd_rand, lrand(simu->init_rand)+(simu->parms->atm.wdrand-1));
    const PARMS_T *parms=simu->parms;
    simu->wfs_rand=mycalloc(parms->nwfs,rand_t);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	seed_rand(&simu->wfs_rand[iwfs],lrand(simu->init_rand));
    }
    seed_rand(simu->telws_rand, lrand(simu->init_rand)*parms->sim.wsseq);
#if USE_CUDA
    if(parms->gpu.wfs && !parms->sim.evlol){
	gpu_wfsgrad_seeding(parms,simu->powfs, simu->init_rand);
    }
#endif
    seed_rand(simu->misc_rand, simu->seed);
}

static void init_simu_evl(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    RECON_T *recon=simu->recon;
    const int nsim=parms->sim.end;
    const int nevl=parms->evl.nevl;
    const int nmod=parms->evl.nmod;
    const int seed=simu->seed;
    SIM_SAVE_T *save=simu->save;
    simu->evlopd=dcellnew(nevl, 1);
    simu->perfevl_iground=parms->atm.iground;
    if(!disable_save && parms->save.extra){
	simu->timing=dnew_mmap(5, nsim, NULL, "Timing_%d.bin", seed);
    }
    {/*MMAP the main result file */
	int ahst_nmod=0;
	if(parms->recon.split){
	    ahst_nmod=3;
	    if(simu->recon->ngsmod->indfocus) ahst_nmod++;
	}
	long nnx[4]={nmod,0,nmod,ahst_nmod};
	long nny[4]={nsim,0,nsim,nsim};

	simu->res = dcellnew_mmap(4,1,nnx,nny,NULL,NULL,"Res_%d.bin",seed);
	//Do not reference. Just assign. Don't free.
	simu->ole     = simu->res->p[0];
	simu->cle     = simu->res->p[2];
	simu->clem    = simu->res->p[3];
	dcellset(simu->res, NAN);
    }

    {/*USE MMAP for data that need to save at every time step */
	long nnx[nevl];
	long nny[nevl];
	for(int ievl=0; ievl<nevl; ievl++){
	    nnx[ievl]=nmod;
	    nny[ievl]=nsim;
	}
	if(parms->save.extra){
	    simu->olep=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resolep_%d.bin",seed);
	    simu->olmp=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resolmp_%d.bin",seed);
	    simu->clep=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resclep_%d.bin",seed);
	    simu->clmp=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resclmp_%d.bin",seed);
	}else{
	    simu->olep=dcellnew3(nevl,1,nnx,nny);
	    simu->olmp=dcellnew3(nevl,1,nnx,nny);
	    simu->clep=dcellnew3(nevl,1,nnx,nny);
	    simu->clmp=dcellnew3(nevl,1,nnx,nny);
	}
	if(parms->recon.split){
	    long nnx_split[nevl];
	    long nnx_3[nevl];
	    for(int ievl=0; ievl<nevl; ievl++){
		nnx_3[ievl]=3;
		nnx_split[ievl]=recon->ngsmod->nmod;
	    }
	    if(parms->save.extra){
		simu->oleNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,NULL,"ResoleNGSm_%d.bin",seed);
		simu->oleNGSmp=dcellnew_mmap(nevl,1,nnx_split,nny,NULL,NULL,"ResoleNGSmp_%d.bin",seed);
		if(!parms->sim.evlol){
		    simu->clemp=dcellnew_mmap(nevl,1, nnx_3, nny, NULL,NULL,"Resclemp_%d.bin",seed);
		    simu->cleNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,NULL,"RescleNGSm_%d.bin",seed);
		    simu->cleNGSmp=dcellnew_mmap(nevl,1,nnx_split,nny,NULL,NULL,"RescleNGSmp_%d.bin",seed);
		    if(parms->recon.split==1 && !parms->sim.fuseint){
			simu->corrNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,NULL,"RescorrNGSm_%d.bin",seed);
		    }
		}
	    }else{
		simu->oleNGSm=dnew(recon->ngsmod->nmod,nsim);
		simu->oleNGSmp=dcellnew3(nevl,1,nnx_split,nny);
		if(!parms->sim.evlol){
		    simu->clemp=dcellnew3(nevl,1, nnx_3, nny);
		    simu->cleNGSm=dnew(recon->ngsmod->nmod,nsim);
		    simu->cleNGSmp=dcellnew3(nevl,1,nnx_split,nny);
		    if(parms->recon.split==1 && !parms->sim.fuseint){
			simu->corrNGSm=dnew(recon->ngsmod->nmod,nsim);
		    }
		}
	    }
	}
    }
    if(parms->sim.skysim){
	char fnold[PATH_MAX];
	char fnnew[PATH_MAX];
	snprintf(fnnew, PATH_MAX, "%s/RescleNGSm_%d.bin",dirskysim,seed);
	snprintf(fnold, PATH_MAX, "RescleNGSm_%d.bin", seed);
	if(exist(fnnew)) remove(fnnew);
	if(link(fnold, fnnew)){
	    warning("Error link\n");
	}
	snprintf(fnnew, PATH_MAX, "%s/RescleNGSmp_%d.bin",dirskysim,seed);
	snprintf(fnold, PATH_MAX, "RescleNGSmp_%d.bin", seed);
	if(exist(fnnew)) remove(fnnew);
	if(link(fnold, fnnew)){
	    warning("Error link\n");
	}
    }

    if(parms->evl.psfmean || parms->evl.psfhist || parms->evl.cov){
	char header[800];
	header[0]='\0';
	if(parms->evl.psfmean || parms->evl.psfhist ){
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		char headeri[00];
		int nembed=aper->embed->nembed->p[iwvl];
		snprintf(headeri, sizeof(headeri), "wvl=%g\nPSF sampling is %.15g radian\nPSF will sum to %.15g\n",
			 parms->evl.wvl->p[iwvl],
			 parms->evl.wvl->p[iwvl]/(nembed*parms->evl.dx),
			 aper->sumamp2*nembed*nembed); 
		strcat(header, headeri);
	    }
	}
	header[sizeof(header)-1]='\0';
	//The saved PSF and COVs are padded by empty cells.
	long nframepsf=parms->sim.end;
	long nframecov=parms->sim.end;
	char strht[24];
	if(parms->evl.psfmean && !parms->sim.evlol){
	    simu->evlpsfmean=dcellnew(parms->evl.nwvl,nevl);
	    simu->evlpsfmean->header=strdup(header);
	    save->evlpsfmean=mycalloc(nevl,zfarr*);
	    simu->evlpsfmean_ngsr=dcellnew(parms->evl.nwvl,nevl);
	    simu->evlpsfmean_ngsr->header=strdup(header);
	    save->evlpsfmean_ngsr=mycalloc(nevl,zfarr*);
	}
	if(parms->evl.psfhist){
	    save->evlpsfhist=mycalloc(nevl,zfarr*);
	    save->evlpsfhist_ngsr=mycalloc(nevl,zfarr*);
	}
	if(parms->evl.cov){
	    simu->evlopdcov=dcellnew(nevl,1);
	    simu->evlopdcov_ngsr=dcellnew(nevl,1);
	    save->evlopdcov=mycalloc(nevl,zfarr*);
	    save->evlopdcov_ngsr=mycalloc(nevl,zfarr*);
	    simu->evlopdmean=dcellnew(nevl,1);
	    simu->evlopdmean_ngsr=dcellnew(nevl,1);
	    save->evlopdmean=mycalloc(nevl,zfarr*);
	    save->evlopdmean_ngsr=mycalloc(nevl,zfarr*);
	}
	for(int ievl=0; ievl<nevl; ievl++){
	    if(!parms->evl.psf->p[ievl]) continue;
	    if(isfinite(parms->evl.hs->p[ievl])){
		snprintf(strht, 24, "_%g", parms->evl.hs->p[ievl]);
	    }else{
		strht[0]='\0';
	    }
	    if(parms->evl.psfngsr->p[ievl]!=2 && !parms->sim.evlol){
		if(parms->evl.psfmean){
		    save->evlpsfmean[ievl]=zfarr_init(parms->evl.nwvl,nframepsf, 
							"evlpsfcl_%d_x%g_y%g%s.fits", seed,
							parms->evl.thetax->p[ievl]*206265,
							parms->evl.thetay->p[ievl]*206265, strht);
		}
		if(parms->evl.psfhist){
		    save->evlpsfhist[ievl]=zfarr_init(parms->sim.end, 1,
							"evlpsfhist_%d_x%g_y%g%s.bin", seed,
							parms->evl.thetax->p[ievl]*206265,
							parms->evl.thetay->p[ievl]*206265,strht);
		}
		if(parms->evl.cov){
		    save->evlopdcov[ievl]=zfarr_init(nframecov, 1,
						       "evlopdcov_%d_x%g_y%g%s.bin", seed,
						       parms->evl.thetax->p[ievl]*206265,
						       parms->evl.thetay->p[ievl]*206265, strht);
		    save->evlopdmean[ievl]=zfarr_init(nframecov, 1,
							"evlopdmean_%d_x%g_y%g%s.bin", seed,
							parms->evl.thetax->p[ievl]*206265,
							parms->evl.thetay->p[ievl]*206265, strht);
		}
	    }
	    if(parms->evl.psfngsr->p[ievl]!=0 && !parms->sim.evlol){
		if(parms->evl.psfmean){
		    save->evlpsfmean_ngsr[ievl]=zfarr_init(parms->evl.nwvl,nframepsf, 
							     "evlpsfcl_ngsr_%d_x%g_y%g%s.fits", seed,
							     parms->evl.thetax->p[ievl]*206265,
							     parms->evl.thetay->p[ievl]*206265, strht);
		}
		if(parms->evl.psfhist){
		    save->evlpsfhist_ngsr[ievl]=zfarr_init(parms->sim.end, 1,
							     "evlpsfhist_ngsr_%d_x%g_y%g%s.bin", seed,
							     parms->evl.thetax->p[ievl]*206265,
							     parms->evl.thetay->p[ievl]*206265,strht);
		}
		if(parms->evl.cov){
		    save->evlopdcov_ngsr[ievl]=zfarr_init(nframecov, 1,
							    "evlopdcov_ngsr_%d_x%g_y%g%s.bin", seed,
							    parms->evl.thetax->p[ievl]*206265,
							    parms->evl.thetay->p[ievl]*206265, strht);
		    save->evlopdmean_ngsr[ievl]=zfarr_init(nframecov, 1,
							     "evlopdmean_ngsr_%d_x%g_y%g%s.bin", seed,
							     parms->evl.thetax->p[ievl]*206265,
							     parms->evl.thetay->p[ievl]*206265, strht);
		}
	    }
	}/*for ievl */
	
	if(parms->evl.psfol){
	    if(parms->evl.psfmean){
		simu->evlpsfolmean=dcellnew(parms->evl.nwvl,1);
		simu->evlpsfolmean->header=strdup(header);
		save->evlpsfolmean=zfarr_init(parms->evl.nwvl, nframepsf, "evlpsfol_%d.fits", seed);
	    }
	    if(parms->evl.cov){
		save->evlopdcovol=zfarr_init(nframecov, 1, "evlopdcovol_%d", seed);
		save->evlopdmeanol=zfarr_init(nframecov, 1, "evlopdmean_%d", seed);
	    }
	}
    }

    if(parms->save.ecov){
	if(!parms->dbg.useopdr || parms->sim.idealfit){
	    simu->ecov=dcellnew(parms->ndm,parms->ndm);/*not really need. */
	}else{/*deprecated */
	    simu->ecov=dcellnew(nevl, 1);
	    if(parms->dbg.ecovxx){/*temporary. */
		warning("Saving history of ecov_xx\n");
		save->ecovxx=mycalloc(nevl,zfarr *);
		char strht[24];
		for(int ievl=0; ievl<nevl; ievl++){
		    if(!parms->evl.psf->p[ievl]) continue;
		    if(isfinite(parms->evl.hs->p[ievl])){
			snprintf(strht, 24, "_%g", parms->evl.hs->p[ievl]);
		    }else{
			strht[0]='\0';
		    }
	    
		    save->ecovxx[ievl]=zfarr_init(parms->sim.end, 1,
						    "ecovxx_%d_x%g_y%g%s.bin", seed,
						    parms->evl.thetax->p[ievl]*206265,
						    parms->evl.thetay->p[ievl]*206265,strht);
		}
	    }
	}
    }

    if(parms->save.evlopd){
	int nstep=parms->sim.end;
	save->evlopdol=mycalloc(nevl,zfarr*);
	save->evlopdcl=mycalloc(nevl,zfarr*);

	for(int ievl=0; ievl<nevl; ievl++){
	    save->evlopdol[ievl]=zfarr_init(nstep,1, 
					      "evl%d_opdol_%d.bin",ievl,seed);
	    save->evlopdcl[ievl]=zfarr_init(nstep,1,
					      "evl%d_opdcl_%d.bin",ievl,seed);
	}
    }

    /* For threading */
    simu->evl_prop_atm=mycalloc(nevl*parms->atm.nps,thread_t*);
    simu->evl_propdata_atm=mycalloc(nevl*parms->atm.nps,PROPDATA_T);
    simu->evl_prop_dm=mycalloc(nevl*parms->ndm,thread_t*);
    simu->evl_propdata_dm=mycalloc(nevl*parms->ndm,PROPDATA_T);
    for(int ievl=0; ievl<nevl; ievl++){
	const int nthread=1;
	int tot;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    const int ind=ievl+nevl*ips;
	    PROPDATA_T *data=&simu->evl_propdata_atm[ind];
	    const real ht=parms->atm.ht->p[ips];
	    data->displacex0=ht*parms->evl.thetax->p[ievl];
	    data->displacey0=ht*parms->evl.thetay->p[ievl];
	    data->scale=1-ht/parms->evl.hs->p[ievl];
	    data->alpha=1;
	    data->wrap=1;
	    data->ostat=aper->locs->stat;
	    tot=aper->locs->stat->ncol;
	    simu->evl_prop_atm[ind]=mycalloc(nthread,thread_t);
	    thread_prep(simu->evl_prop_atm[ind], 0, tot, nthread, prop, data);
	}
	for(int idm=0; idm<parms->ndm && !parms->sim.evlol; idm++){
	    const int ind=ievl+nevl*idm;
	    PROPDATA_T *data=&simu->evl_propdata_dm[ind];
	    const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
	    data->displacex0=ht*parms->evl.thetax->p[ievl];
	    data->displacey0=ht*parms->evl.thetay->p[ievl];
	    data->scale=1-ht/parms->evl.hs->p[ievl];
	    data->alpha=-1;
	    data->wrap=0;
	    if(parms->sim.cachedm){
		data->mapin=simu->cachedm->p[idm];
	    }else{
		if(simu->dmrealsq){
		    data->mapin=simu->dmrealsq->p[idm];
		}else{
		    data->locin=recon->aloc->p[idm];
		    data->phiin=simu->dmreal->p[idm]->p;
		}
	    }
	    if(aper->locs_dm){
		data->locout=aper->locs_dm->p[ind];
		tot=data->locout->nloc;
	    }else{
		data->ostat=aper->locs->stat;
		tot=aper->locs->stat->ncol;
	    }

	    data->phiout=(real*)1;/*replace later in simulation. */
	    simu->evl_prop_dm[ind]=mycalloc(nthread,thread_t);
	    thread_prep(simu->evl_prop_dm[ind], 0, tot, nthread, prop, data);
	}
    }
}

static void init_simu_wfs(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.idealfit || parms->sim.idealtomo) return;
    POWFS_T *powfs=simu->powfs;
    RECON_T *recon=simu->recon;
    SIM_SAVE_T *save=simu->save;
    const int nwfs=parms->nwfs;
    const int nsim=parms->sim.end;
    const int seed=simu->seed;
    simu->eptwfs=recon->eptwfs;
    simu->ints=dccellnew(nwfs, 1);
    simu->intsout=dccellnew(nwfs, 1);
    simu->wfspsfout=cccellnew(nwfs, 1);
    simu->pistatout=dccellnew(nwfs, 1);
    save->wfspsfout=mycalloc(nwfs,zfarr*);
    save->ztiltout =mycalloc(nwfs,zfarr*);
    simu->gradcl=dcellnew(nwfs,1);
    simu->wfsopd=dcellnew(nwfs,1);
    /*Do not initialize gradlastcl. Do not initialize gradlastol in open
      loop. They are used for testing*/
    if(parms->sim.closeloop){
	long nnx[parms->nwfsr];
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    nnx[iwfs]=powfs[ipowfs].saloc->nloc*2;
	}
	simu->gradlastol=dcellnew3(parms->nwfsr, 1, nnx, NULL);
    }
    simu->gradacc=dcellnew(nwfs,1);/*wfsgrad internal */
    simu->gradoff=dcellnew(nwfs,1);
    simu->gradscale=dcellnew(nwfs,1);//leave empty first
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].saloc->nloc;
	simu->gradcl->p[iwfs]=dnew(nsa*2,1);
	simu->wfsopd->p[iwfs]=dnew(powfs[ipowfs].loc->nloc, 1);
	if(parms->powfs[ipowfs].usephy){
	    if(parms->powfs[ipowfs].type==0){
		simu->ints->p[iwfs]=dcellnew_same(nsa,1,powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
	    }else{
		simu->ints->p[iwfs]=dcellnew_same(1,1,powfs[ipowfs].saloc->nloc,powfs[ipowfs].pywfs->nside);
	    }
	}
	if(parms->powfs[ipowfs].phystep!=0 || parms->save.gradgeom->p[iwfs] || parms->powfs[ipowfs].pistatout){
	    simu->gradacc->p[iwfs]=dnew(nsa*2,1);
	}
	if(parms->powfs[ipowfs].pistatout){
	    simu->pistatout->p[iwfs]=dcellnew(nsa,parms->powfs[ipowfs].nwvl);
	}
	if(powfs[ipowfs].gradncpa && !(parms->powfs[ipowfs].phytype_sim==1 && parms->powfs[ipowfs].ncpa_method==2)){
	    //CMF has gradncpa with in matched filter
	    int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	    dadd(&simu->gradoff->p[iwfs], 1, PR(powfs[ipowfs].gradncpa, wfsind,1), 1);
	}

    }
    if(parms->sim.mffocus){
	if(fabs(parms->sim.lpfocushi)<1.e-15){
	    warning("sim.mffocus is enabled but sim.lpfocus is zero.\n");
	}
	simu->lgsfocuslpf=dnew(parms->nwfs, 1);
	simu->ngsfocuslpf=0;
    }
    if(parms->nphypowfs){
	simu->fsmerr_store=dcellnew(nwfs,1);
	simu->fsmreal=dcellnew(nwfs,1);
	simu->fsmsho=mycalloc(nwfs*2,SHO_T*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].llt || parms->powfs[ipowfs].dither==1){
		simu->fsmerr_store->p[iwfs]=dnew(2,1);
		simu->fsmreal->p[iwfs]=dnew(2,1);
		if(parms->sim.f0fsm>0){
		    simu->fsmsho[iwfs]=sho_new(parms->sim.f0fsm, parms->sim.zetafsm);//x
		    simu->fsmsho[iwfs+nwfs]=sho_new(parms->sim.f0fsm, parms->sim.zetafsm);//y
		}
	    }
	}
	simu->fsmint=servo_new(simu->fsmreal, parms->sim.apfsm, parms->sim.alfsm,
			       parms->sim.dthi, parms->sim.epfsm);


    }
  
    {/*MMAP the LGS fsmlink error/command output */
	long nnx[nwfs];
	long nny[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    nnx[iwfs]=0;
	    if(parms->powfs[ipowfs].llt || parms->powfs[ipowfs].dither==1){
		nnx[iwfs]+=2;//tip/tilt of FSM
	    }
	    if(parms->powfs[ipowfs].dither>1){
		nnx[iwfs]+=1;//signal in common path dithering
	    }
	    if(nnx[iwfs]){
		nny[iwfs]=nsim;
	    }else{
		nny[iwfs]=0;
	    }
	}
	if(parms->save.extra){
	    simu->fsmerrs = dcellnew_mmap(nwfs, 1, nnx, nny, NULL,NULL,"Resfsmerr_%d.bin", seed);
	    simu->fsmcmds = dcellnew_mmap(nwfs, 1, nnx, nny, NULL,NULL,"Resfsmcmd_%d.bin", seed);
	}else{
	    simu->fsmerrs = dcellnew3(nwfs, 1, nnx, nny);
	    simu->fsmcmds = dcellnew3(nwfs, 1, nnx, nny);
	}
    }
   
    /* For sky coverage telemetry output */
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	long ncompx=powfs[ipowfs].ncompx;
	long ncompy=powfs[ipowfs].ncompy;
	long notf=MAX(ncompx, ncompy);
	if(parms->powfs[ipowfs].psfout){
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    /*The PSFs here are PSFs of each subaperture. */
	    simu->wfspsfout->p[iwfs]=ccellnew_same(nsa,parms->powfs[ipowfs].nwvl,notf/2+2,notf/2+2);
	    mymkdir("%s/wvfout/", dirskysim);
	    mymkdir("%s/ztiltout/", dirskysim);
	    save->wfspsfout[iwfs]=zfarr_init
		(parms->sim.end-parms->sim.start,1,
		 "%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g.bin",
		 dirskysim, seed,
		 parms->powfs[ipowfs].order,
		 parms->wfs[iwfs].thetax*206265,
		 parms->wfs[iwfs].thetay*206265);
	
	    save->ztiltout[iwfs]=zfarr_init
		(parms->sim.end-parms->sim.start,1,
		 "%s/ztiltout/ztiltout_seed%d_sa%d_x%g_y%g.bin",
		 dirskysim, seed,
		 parms->powfs[ipowfs].order,
		 parms->wfs[iwfs].thetax*206265,
		 parms->wfs[iwfs].thetay*206265);

	}
	if(parms->sim.skysim && parms->powfs[ipowfs].pistatout){
	    mymkdir("%s/pistat/",dirskysim);
	}
    }
 
    if(parms->save.ngcov>0){
	simu->gcov=dcellnew(parms->save.ngcov,1);
    }
    int nstep=parms->sim.end;
    if(parms->save.wfsopd){
	save->wfsopd=mycalloc(nwfs,zfarr*);
	save->wfsopdol=mycalloc(nwfs,zfarr*);
	save->wfslltopd=mycalloc(nwfs,zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(!parms->save.wfsopd->p[iwfs]){
		continue;
	    }
	    save->wfsopd[iwfs]=zfarr_init(nstep, 1,"wfs%d_opd_%d.bin", iwfs, seed);
	    save->wfsopdol[iwfs]=zfarr_init(nstep, 1,"wfs%d_opdol_%d.bin", iwfs, seed);
	    if(powfs[ipowfs].llt){
		save->wfslltopd[iwfs]=zfarr_init(nstep,1, "wfs%d_lltopd_%d.bin", iwfs, seed);
	    }
	}
    }
    if(parms->save.ints){
	save->intsny=mycalloc(nwfs,zfarr*);
	save->intsnf=mycalloc(nwfs,zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->save.ints->p[iwfs] && parms->powfs[ipowfs].usephy){
		save->intsnf[iwfs]=zfarr_init(nstep,1, "wfs%d_intsnf_%d.bin", iwfs, seed);
		if(parms->powfs[ipowfs].noisy){
		    save->intsny[iwfs]=zfarr_init(nstep,1, "wfs%d_intsny_%d.bin", iwfs, seed);
		}
	    }
	}
    }

    save->gradcl=mycalloc(nwfs,zfarr*);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(parms->save.grad->p[iwfs]){
	    save->gradcl[iwfs]=zfarr_init(nstep,1, "wfs%d_gradcl_%d.bin", iwfs, seed);
	}
    }

    save->gradol=mycalloc(nwfs,zfarr*);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].psol && parms->save.gradpsol->p[iwfs]){
	    save->gradol[iwfs]=zfarr_init(nstep,1, 
					  "wfs%d_gradpsol_%d.bin", iwfs, seed);
	}
    }

    save->gradnf=mycalloc(nwfs,zfarr*);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->save.gradnf->p[iwfs] && parms->powfs[ipowfs].noisy){
	    save->gradnf[iwfs]=zfarr_init(nstep,1, "wfs%d_gradnf_%d.bin", iwfs, seed);
	}
    }

    save->gradgeom=mycalloc(nwfs,zfarr*);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->save.gradgeom->p[iwfs] && parms->powfs[ipowfs].phystep>-1){
	    save->gradgeom[iwfs]=zfarr_init(nstep,1, "wfs%d_gradgeom_%d.bin", iwfs, seed);
	}
    }
    
   
    /* threading */
    simu->wfs_prop_atm=mycalloc(nwfs*parms->atm.nps,thread_t*);
    simu->wfs_propdata_atm=mycalloc(nwfs*parms->atm.nps,PROPDATA_T);
    simu->wfs_prop_dm=mycalloc(nwfs*parms->ndm,thread_t*);
    simu->wfs_propdata_dm=mycalloc(nwfs*parms->ndm,PROPDATA_T);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int nwfsp=parms->powfs[ipowfs].nwfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const real hs=parms->wfs[iwfs].hs;
	const real hc=parms->wfs[iwfs].hc;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    const real ht=parms->atm.ht->p[ips];
	    PROPDATA_T *data=&simu->wfs_propdata_atm[iwfs+nwfs*ips];
	    data->displacex0=ht*parms->wfs[iwfs].thetax;
	    data->displacey0=ht*parms->wfs[iwfs].thetay;
	    data->scale=1.-(ht-hc)/hs;
	    data->alpha=1;
	    data->wrap=1;
	    data->mapin=(map_t*)1;/*need to update this in genatm. */
	    data->phiout=(real*)1;/*replace later in simulation. */
	    int tot=0;
	    if(powfs[ipowfs].loc_tel){/*misregistration. */
		data->locout=powfs[ipowfs].loc_tel->p[wfsind];
		tot=data->locout->nloc;
	    }else if(parms->powfs[ipowfs].type==1 || powfs[ipowfs].saloc->nloc<NTHREAD){
		data->locout=powfs[ipowfs].loc;
		tot=data->locout->nloc;
	    }else{
		data->ptsout=powfs[ipowfs].pts;
		tot=data->ptsout->nsa;
	    }
	    int nthread=1;//NTHREAD;
	    simu->wfs_prop_atm[iwfs+nwfs*ips]=mycalloc(nthread,thread_t);
	    thread_prep(simu->wfs_prop_atm[iwfs+nwfs*ips],0,tot,nthread,prop,data);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    const real ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	    PROPDATA_T *data=&simu->wfs_propdata_dm[iwfs+nwfs*idm];
	    int tot;
	    data->displacex0=ht*parms->wfs[iwfs].thetax;
	    data->displacey0=ht*parms->wfs[iwfs].thetay;
	    data->scale=1.-(ht-hc)/hs;
	    data->alpha=-1;/*remove dm contribution. */
	    data->wrap=0;
	    if(parms->sim.cachedm){
		data->mapin=simu->cachedm->p[idm];
	    }else{
		if(simu->dmrealsq){
		    data->mapin=simu->dmrealsq->p[idm];
		}else{
		    data->locin=recon->aloc->p[idm];
		    data->phiin=simu->dmreal->p[idm]->p;
		}
	    }
	    data->phiout=(real*)1;/*replace later in simulation */
	    if(powfs[ipowfs].loc_dm){/*misregistration. */
		data->locout=powfs[ipowfs].loc_dm->p[wfsind+idm*nwfsp];
		tot=data->locout->nloc;
	    }else if(parms->powfs[ipowfs].type==1 || powfs[ipowfs].saloc->nloc<NTHREAD){
		data->locout=powfs[ipowfs].loc;
		tot=data->locout->nloc;
	    }else{
		data->ptsout=powfs[ipowfs].pts;
		tot=data->ptsout->nsa;
	    }
	    int nthread=1;//NTHREAD
	    simu->wfs_prop_dm[iwfs+nwfs*idm]=mycalloc(nthread,thread_t);
	    thread_prep(simu->wfs_prop_dm[iwfs+nwfs*idm], 0, tot, nthread, prop,data);
	}/*idm */
    }/*iwfs */
    simu->wfs_intsdata=mycalloc(nwfs,WFSINTS_T);
    simu->wfs_ints=mycalloc(nwfs,thread_t*);

    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int tot=powfs[ipowfs].saloc->nloc;
	WFSINTS_T *data=simu->wfs_intsdata+iwfs;
	data->iwfs=iwfs;
	int nthread=NTHREAD;
	simu->wfs_ints[iwfs]=mycalloc(nthread,thread_t);
	thread_prep(simu->wfs_ints[iwfs], 0, tot, nthread, wfsints,data);
    }
    if(parms->nlgspowfs){
	simu->llt_tt=dcellnew(parms->nwfs, 1);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].llt && parms->powfs[ipowfs].llt->ttpsd){
		dmat *psdin=dread("%s", parms->powfs[ipowfs].llt->ttpsd);
		simu->llt_tt->p[iwfs]=psd2time(psdin, simu->misc_rand, parms->sim.dt, parms->sim.end);
	    }
	}
	//writebin(simu->llt_tt, "llt_tt_%d", seed);
    }
    if(parms->nlgspowfs){
	simu->LGSfocus=dcellnew(parms->nwfs,1);
	simu->zoomerr=dnew(parms->nwfs,1);
	simu->zoomint=dnew(parms->nwfs,1);
	simu->zoomavg=dnew(parms->nwfs, 1);
	if(!disable_save){
	    long nnx[parms->nwfs];
	    long nny[parms->nwfs];
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].llt){
		    nnx[iwfs]=parms->sim.end;
		    nny[iwfs]=1;
		}else{
		    nnx[iwfs]=0;
		    nny[iwfs]=0;
		}
	    }
	    if(parms->save.extra){
		simu->zoompos=dcellnew_mmap(parms->nwfs, 1, nnx, nny, NULL, NULL, "Reszoompos_%d.bin", seed);
	    }else{
		simu->zoompos=dcellnew3(parms->nwfs, 1, nnx, nny);
	    }
	}
	
    }
    if(parms->dither){
	simu->dither=mycalloc(nwfs,DITHER_T*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].dither){
		simu->dither[iwfs]=mycalloc(1,DITHER_T);
		if(parms->powfs[ipowfs].dither!=1){
		    simu->dither[iwfs]->mr=dcellnewsame_mmap(2, 1, 1, nsim, NULL, "Resdithermr_%d", seed);
		}
	    }
	}
	if(parms->save.extra){
	    long nnx[nwfs];
	    long nny[nwfs];
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].dither){
		    nnx[iwfs]=4;
		    nny[iwfs]=(nsim-parms->powfs[ipowfs].dither_pllskip)/(parms->powfs[ipowfs].dtrat*parms->powfs[ipowfs].dither_pllrat);
		}else{
		    nnx[iwfs]=0;
		    nny[iwfs]=0;
		}
	    }
	    simu->resdither = dcellnew_mmap(nwfs, 1, nnx, nny, NULL,NULL,"Resdither_%d.bin", seed);
	}
    }
    if(recon->cn2est){
	simu->cn2res=dcellnew(2,1);
	int ncn2=(parms->sim.end-1)/parms->cn2.step;
	long nnx[2]={ncn2, recon->cn2est->htrecon->nx};
	long nny[2]={1, ncn2};
	simu->cn2res=dcellnew_mmap(2, 1, nnx, nny, NULL, NULL, "Rescn2_%d.bin", seed);
    }
    if(parms->itpowfs!=-1){
	const int ipowfs=parms->itpowfs;
	const int iwfs=parms->powfs[ipowfs].wfs->p[0];
	const int nmod=recon->RRtwfs->p[iwfs]->nx;
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const int nacc=parms->sim.end-parms->powfs[ipowfs].step;
	if(nacc>0){
	    simu->restwfs=dnew_mmap(nmod, nacc/dtrat, NULL, "Restwfs_%d.bin", seed);
	}
    }
}

static void init_simu_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    SIM_SAVE_T *save=simu->save;
    const int seed=simu->seed;
    /*Setup hysterisis */
    for(int idm=0; idm<parms->ndm; idm++){
	if(parms->dm[idm].hyst){
	    if(!simu->hyst){
		simu->hyst = mycalloc(parms->ndm,HYST_T*);
	    }
	    dmat *coeff=dread("%s",parms->dm[idm].hyst);
	    simu->hyst[idm]=hyst_new(coeff, recon->aloc->p[idm]->nloc);
	    hyst_calib(simu->hyst[idm], idm);
	    dfree(coeff);
	}
    }
    /*we initialize dmreal, so that wfs_prop_dm can reference dmreal. */
    if(1){
	simu->dmerr_store=dcellnew3(parms->ndm,1, parms->recon.modal?recon->anmod->p:recon->anloc->p, NULL);
    }else{
	simu->dmerr_store=dcellnew_mmap(parms->ndm,1, parms->recon.modal?recon->anmod->p:recon->anloc->p, NULL, NULL, NULL, "/dmerr");
    }
    simu->dmcmd=dcellnew(parms->ndm,1);
    simu->dmreal=dcellnew(parms->ndm,1);
    if(parms->fit.cgwarm){
	simu->dmfit=dcellnew3(parms->ndm,1, parms->recon.modal?recon->anmod->p:recon->anloc->p, NULL);
    }else{
	simu->dmfit=simu->dmerr_store;
    }
    if(parms->sim.lpttm>EPS){
	simu->ttmreal=dnew(2,1);
    }
    simu->dmrealsq=mapcellnew(parms->ndm, 1);

    if(parms->sim.dmproj){
	simu->dmproj=dcellnew3(parms->ndm,1, recon->anloc->p, NULL);
	simu->dmprojsq=mapcellnew(parms->ndm, 1);
    }
    for(int idm=0; idm<parms->ndm; idm++){
	simu->dmcmd->p[idm]=dnew(recon->anloc->p[idm],1);
	//Do not reference for the new synchronization scheme.
	simu->dmreal->p[idm]=dnew(recon->anloc->p[idm],1);
	if(simu->dmrealsq){
	    simu->dmrealsq->p[idm]=mapnew2(recon->amap->p[idm]);
	    dset((dmat*)simu->dmrealsq->p[idm], NAN);
	}
	if(simu->dmprojsq){
	    simu->dmprojsq->p[idm]=mapnew2(recon->amap->p[idm]);
	    dset((dmat*)simu->dmprojsq->p[idm], NAN);
	}
	if(parms->fit.square){/*dmreal is already square.*/
//	    free(simu->dmrealsq->p[idm]->p);
	    mem_unref(simu->dmrealsq->p[idm]->mem);
	    simu->dmrealsq->p[idm]->p=simu->dmreal->p[idm]->p;
	    mem_ref(simu->dmreal->p[idm]->mem);
//	    free(simu->dmrealsq->p[idm]->nref);
//	    simu->dmrealsq->p[idm]->nref=NULL;
	    if(simu->dmprojsq){
		mem_unref(simu->dmprojsq->p[idm]->mem);
		simu->dmprojsq->p[idm]->p=simu->dmproj->p[idm]->p;
		mem_ref(simu->dmprojsq->p[idm]->mem);
//		simu->dmprojsq->p[idm]->nref=NULL;
	    }
	}
    }
    if(parms->sim.dmadd){
	simu->dmadd=dcellread("%s", parms->sim.dmadd);
	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmadd->p[idm] && simu->dmadd->p[idm]->nx!=simu->dmreal->p[idm]->nx){
		if(!parms->fit.square || simu->dmadd->p[idm]->nx!=recon->amap->p[idm]->nx*recon->amap->p[idm]->ny){
		    error("DM[%d]: dmadd does not match dm geometry\n", idm);
		}
	    }
	}
    }
    if(parms->sim.cachedm){
	prep_cachedm(simu);
    }
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	gpu_dmreal2gpu(simu->dmrealsq);
	if(simu->dmprojsq){
	    gpu_dmproj2gpu(simu->dmprojsq);
	}
    }
#endif
    simu->wfspsol=dccellnew(parms->npowfs, 1);
    simu->dmint=servo_new(simu->dmerr_store, parms->sim.aphi, parms->sim.alhi, 
			  parms->sim.dthi, parms->sim.ephi);
    if(parms->dbg.ncpa_preload && recon->dm_ncpa){//set the integrator
	warning_once("Preload integrator with NCPA\n");
	dcelladd(&simu->dmint->mint->p[0], 1, recon->dm_ncpa, 1);
    }
    if(parms->recon.split){
	simu->Merr_lo_store=dcellnew_same(1,1,recon->ngsmod->nmod,1);
	simu->Merr_lo2=dcellnew_same(1,1,recon->ngsmod->nmod,1);
	simu->Mint_lo=servo_new(simu->Merr_lo_store, 
				parms->sim.aplo, parms->sim.allo,
				parms->sim.dtlo, parms->sim.eplo);
    }
    {/* History */
	int dm_hist=0;
	for(int idm=0; idm<parms->ndm; idm++){
	    if(parms->dm[idm].hist){
		dm_hist=1; 
		break;
	    }
	}
	if(dm_hist){
	    long nnx[parms->ndm];
	    long nny[parms->ndm];
	    for(int idm=0; idm<parms->ndm; idm++){
		if(parms->dm[idm].hist){
		    nnx[idm]=parms->dm[idm].histn;
		    nny[idm]=recon->aloc->p[idm]->nloc;
		}else{
		    nnx[idm]=0;
		    nny[idm]=0;
		}
	    }
	    simu->dmhist=dcellnew_mmap(parms->ndm, 1, nnx, nny, 
				       NULL, NULL,"dmhist_%d.bin",simu->seed);
	}
    }
    if(parms->recon.psd){
	simu->dmerrts=dcellnew_same(parms->evl.nevl, 1, recon->Herr->p[0]->nx,
				   parms->recon.psddtrat_hi);
	if(parms->recon.split){
	    simu->Merrts=dnew(recon->ngsmod->nmod, parms->recon.psddtrat_lo);
	}
    }
    int nstep=parms->sim.end;
    int nrstep=nstep-(parms->sim.closeloop?1:0);
    if(parms->save.dm){
	save->dmerr=zfarr_init(nrstep, 1,"dmerr_%d.bin", seed);
	save->dmfit=zfarr_init(nrstep, 1, "dmfit_%d.bin", seed);
	if(parms->recon.split){
	    save->Merr_lo=zfarr_init(nstep, 1, "Merr_lo_%d.bin", seed);
	    if(!parms->sim.fuseint){
		save->Mint_lo=zfarr_init(nstep, 1, "Mint_lo_%d.bin", seed);
	    }
	}
	if(parms->sim.lpttm>EPS){
	    save->ttmreal= dnew_mmap(2, nstep, NULL, "ttmreal_%d.bin", seed);
	}
	save->dmint=zfarr_init(nstep, 1,"dmint_%d.bin", seed);
	save->dmreal = zfarr_init(nstep, 1, "dmreal_%d.bin", seed);
	save->dmcmd  = zfarr_init(nstep, 1, "dmcmd_%d.bin", seed);
	if(parms->sim.dmproj){
	    save->dmproj = zfarr_init(nstep, 1, "dmproj_%d.bin", seed);
	}
    }
}

/** MOAO **/
static void init_simu_moao(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    SIM_SAVE_T *save=simu->save;
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    const int seed=simu->seed;
    int nstep=parms->sim.end;
    int ny=parms->sim.closeloop && !parms->gpu.moao ? 2 : 1;
    if(!parms->gpu.wfs || ! parms->gpu.moao || parms->plot.run){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao>-1){
		if(!simu->dm_wfs){
		    simu->dm_wfs=dcellnew(nwfs, ny);
		}
		for(int iy=0; iy<ny; iy++){
		    simu->dm_wfs->p[iwfs+iy*nwfs]=dnew(recon->moao[imoao].aloc->p[0]->nloc, 1);
		}
	    }
	}
    }
    if(!parms->gpu.evl || ! parms->gpu.moao || parms->plot.run){ 
	int imoao=parms->evl.moao;
	if(imoao!=-1){
	    simu->dm_evl=dcellnew(nevl, ny);
	    for(int ievl=0; ievl<nevl; ievl++){
		for(int iy=0; iy<ny; iy++){
		    simu->dm_evl->p[ievl+iy*nevl]=dnew(recon->moao[imoao].aloc->p[0]->nloc, 1);
		}
	    }
	}
    }
/*#if USE_CUDA
    if(parms->gpu.moao && recon->moao){
	gpu_moao_2gpu(simu);//initilization.
    }
    #endif*/
    if(parms->save.dm){
	if(simu->dm_wfs){
	    save->dm_wfs=mycalloc(nwfs,zfarr*);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int imoao=parms->powfs[ipowfs].moao;
		if(imoao>-1){
		    save->dm_wfs[iwfs]=zfarr_init(nstep,1,"wfs%d_moaofit_%d.bin",iwfs,seed);
		}
	    }
	}
	if(simu->dm_evl){
	    save->dm_evl=mycalloc(nwfs,zfarr*);
	    for(int ievl=0; ievl<nevl; ievl++){
		save->dm_evl[ievl]=zfarr_init(nstep,1, "evl%d_moaofit_%d.bin",ievl,seed);
	    }
	}
    }
}

/**
   Initialize simu (of type SIM_T) and various simulation data structs. Called
   for every seed.
*/
SIM_T* init_simu(const PARMS_T *parms,POWFS_T *powfs, 
		 APER_T *aper, RECON_T *recon, int iseed){
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    SIM_T *simu=mycalloc(1,SIM_T);
    global->simu=simu;
    simu->pause=parms->sim.pause;
    pthread_cond_init(&simu->dmreal_condr,0);
    pthread_cond_init(&simu->dmreal_condw,0);
    pthread_mutex_init(&simu->dmreal_mutex,0);
    simu->dmreal_isim=-1;

    pthread_cond_init(&simu->wfsgrad_condr,0);
    pthread_cond_init(&simu->wfsgrad_condw,0);
    pthread_mutex_init(&simu->wfsgrad_mutex,0);
    simu->wfsgrad_isim=-1;

    SIM_SAVE_T *save=simu->save=mycalloc(1,SIM_SAVE_T);
    simu->wfsisim=-1;
    simu->perfisim=-1;
    simu->reconisim=-2;
    simu->parms=parms;
    simu->powfs=powfs;
    simu->recon=recon;
    simu->aper=aper;    
    simu->iseed=iseed;
    simu->seed=parms->sim.seeds->p[iseed];
    if(simu->seed==0){
	simu->seed=myclocki();
	warning("Seed is 0, using time as seed: %d\n",simu->seed);
    }
    int seed=simu->seed;
    seeding(simu);

    simu->status=mycalloc(1,STATUS_T);
    simu->status->iseed=iseed;
    simu->status->nseed=parms->sim.nseed;
    simu->status->simstart=parms->sim.start;
    simu->status->simend=parms->sim.end;
    simu->status->nthread=NTHREAD;
    simu->status->timstart=myclocki();
    simu->status->info=S_RUNNING;

    if(parms->sim.wspsd){
	/* Telescope wind shake added to TT input. */
	dmat *psdin=dread("%s", parms->sim.wspsd);
	info("Loading windshake PSD from file %s\n", parms->sim.wspsd);
	simu->telws = psd2time(psdin, simu->telws_rand, parms->sim.dt, parms->sim.end);
	dfree(psdin);
	writebin(simu->telws, "telws_%d", seed);
    }

    /* Select GPU or CPU for the tasks.*/
    simu->wfsgrad_pre=mycalloc(nwfs,thread_t);
    simu->wfsgrad_post=mycalloc(nwfs,thread_t);
    simu->perfevl_pre=mycalloc(nevl,thread_t);
#if USE_CUDA
    if(parms->gpu.evl){
	thread_prep(simu->perfevl_pre, 0, nevl, nevl, gpu_perfevl_queue, simu);
	simu->perfevl_post=mycalloc(nevl,thread_t);
	thread_prep(simu->perfevl_post, 0, nevl, nevl, gpu_perfevl_sync, simu);
    }else
#endif
    {
	thread_prep(simu->perfevl_pre, 0, nevl, nevl, perfevl_ievl, simu);
    }
#if USE_CUDA
    if(parms->gpu.wfs){
	thread_prep(simu->wfsgrad_pre, 0, nwfs, nwfs, gpu_wfsgrad_queue, simu);
    }else
#endif
    {
	thread_prep(simu->wfsgrad_pre, 0, nwfs, nwfs, wfsgrad_iwfs, simu);
    }
    thread_prep(simu->wfsgrad_post, 0, nwfs, nwfs, wfsgrad_post, simu);
    
    if(!parms->sim.evlol){
	init_simu_dm(simu);
	init_simu_moao(simu);
	if(!parms->sim.idealfit && !parms->sim.idealtomo){
	    init_simu_wfs(simu);
	}
	if(parms->recon.alg==0){
	    int nstep=parms->sim.end;
	    if(parms->save.opdr){
		save->opdr=zfarr_init(nstep-1,1, "opdr_%d.bin", seed);
	    }
	    if(parms->save.opdx){
		save->opdx=zfarr_init(nstep, 1,"opdx_%d.bin", seed);
	    }
	}
	{
	    long nnx[2], nny[2];
	    nnx[1]=nnx[0]=parms->sim.end;
	    nny[1]=nny[0]=1;
	    const char *header[2]={"Tomography", "DM Fit"};
	    if(parms->save.extra){
		simu->cgres=dcellnew_mmap(2, 1, nnx, nny, "CG Residual", header,
					  "ResCG_%d.bin", seed);
	    }else{
		simu->cgres=dcellnew3(2, 1, nnx, nny);
	    }
	}
    }
    init_simu_evl(simu);
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	if(parms->gpu.evl){
	    gpu_perfevl_init_sim(parms, aper);
	}
	if(parms->gpu.wfs && !parms->sim.evlol){
	    gpu_wfs_init_sim(parms, powfs);
	}
    }
    if(parms->gpu.tomo || parms->gpu.fit){
	gpu_recon_reset(parms);
    }
#endif
    OMPTASK_SINGLE
	filter_dm(simu);//2014-03-31. //so that dm_ncpa is effective at first cycle. replaced by copy dm_ncpa to dmreal.
    return simu;
}
/**
   Release memory of simu (of type SIM_T) and close files.
*/
void free_simu(SIM_T *simu){
    if(!simu) return;
    global->simu=0;
    const PARMS_T *parms=simu->parms;
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    free(simu->init_rand);
    free(simu->atm_rand);
    free(simu->atmwd_rand);
    free(simu->wfs_rand);
    free(simu->telws_rand);
    free(simu->misc_rand);
    dfree(simu->telws);
    if(simu->atmcfg){
	free(simu->atmcfg);
    }
    cellfree(simu->atm);
    dfree(simu->atmscale);
    if(parms->sim.cachedm && !parms->sim.evlol){
	cellfree(simu->cachedm);
	for(int i=0; i<parms->ndm; i++){
	    free(simu->cachedm_prop[i]);
	}
	free(simu->cachedm_prop);
	free(simu->cachedm_propdata);
	
    }
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	for(int ips=0; ips<parms->atm.nps; ips++){
	    free(simu->wfs_prop_atm[iwfs+nwfs*ips]);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    free(simu->wfs_prop_dm[iwfs+nwfs*idm]);
	}
	free(simu->wfs_ints[iwfs]);
    }
    for(int ievl=0; ievl<nevl; ievl++){
	for(int ips=0; ips<parms->atm.nps; ips++){
	    int ind=ievl+nevl*ips;
	    free(simu->evl_prop_atm[ind]);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    int ind=ievl+nevl*idm;
	    free(simu->evl_prop_dm[ind]);
	}
    }
    if(simu->hyst){
	for(int idm=0; idm<parms->ndm; idm++){
	    hyst_free(simu->hyst[idm]);
	}
	free(simu->hyst);
    }
    free(simu->wfs_prop_atm);
    free(simu->wfs_prop_dm);
    free(simu->wfs_ints);
    free(simu->wfs_propdata_atm);
    free(simu->wfs_propdata_dm);
    free(simu->wfs_intsdata);
    free(simu->evl_prop_atm);
    free(simu->evl_propdata_atm);
    free(simu->evl_prop_dm);
    free(simu->evl_propdata_dm);
    free(simu->wfsgrad_pre);
    free(simu->wfsgrad_post);
    free(simu->perfevl_pre);
    free(simu->perfevl_post);
    free(simu->status);
    dfree(simu->opdevlground);
    dcellfree(simu->wfsopd);
    dcellfree(simu->gradcl);
    dcellfree(simu->gradacc);
    dcellfree(simu->gradlastcl);
    dcellfree(simu->gradlastol);
    dcellfree(simu->gradoff);
    dcellfree(simu->gradscale);
    dcellfree(simu->opdr);
    dcellfree(simu->gngsmvst);
    dcellfree(simu->cn2res);
    dcellfree(simu->dmreal);
    dcellfree(simu->dmpsol);
    dfree(simu->ttmreal);
    cellfree(simu->dmrealsq);
    dcellfree(simu->dmproj);
    cellfree(simu->wfspsol);
    dcellfree(simu->dmcmd);
    dcellfree(simu->dmtmp);
    dcellfree(simu->dmadd);
    servo_free(simu->dmint);
    servo_free(simu->Mint_lo);
    dcellfree(simu->Mngs);
    dcellfree(simu->dmerrts);
    dcellfree(simu->gcov);
    dcellfree(simu->ecov);
    if(simu->dmfit!=simu->dmerr_store){
	dcellfree(simu->dmfit);
    }
    dcellfree(simu->dmerr_store);
    dcellfree(simu->dmhist);
    dcellfree(simu->Merr_lo_store);
    dcellfree(simu->Merr_lo2);
    dcellfree(simu->fsmerr_store);
    dcellfree(simu->fsmreal);
    servo_free(simu->fsmint);
    if(simu->fsmsho){
	for(int i=0; i<parms->nwfs*2; i++){
	    free(simu->fsmsho[i]);
	}
	free(simu->fsmsho);
    }
    dcellfree(simu->cgres);
    dcellfree(simu->dm_wfs);
    dcellfree(simu->dm_evl);
    dcellfree(simu->res);
    dfree(simu->timing);
    dcellfree(simu->olep);
    dcellfree(simu->olmp);
    dcellfree(simu->clep);
    dcellfree(simu->clmp);
    dcellfree(simu->fsmerrs);
    dcellfree(simu->fsmcmds);
    dcellfree(simu->LGSfocus);
    dcellfree(simu->telfocusint);
    dcellfree(simu->telfocusreal);
    dfree(simu->zoomerr);
    dfree(simu->zoomreal);
    dfree(simu->zoomavg);
    dfree(simu->zoomint);
    if(parms->recon.split){
	dcellfree(simu->clemp);
	dfree(simu->cleNGSm);
	dfree(simu->oleNGSm);
	dfree(simu->corrNGSm);
	dcellfree(simu->cleNGSmp);
	dcellfree(simu->oleNGSmp);
    }
    SIM_SAVE_T *save=simu->save;
    zfarr_close_n(save->evlpsfhist, nevl);
    zfarr_close_n(save->evlpsfhist_ngsr, nevl);
    dcellfree(simu->evlpsfmean);
    dcellfree(simu->evlpsfolmean);
    dcellfree(simu->evlopdcov);
    dcellfree(simu->evlopdmean);
    dfree(simu->evlopdcovol);
    dfree(simu->evlopdmeanol);
    dcellfree(simu->evlpsfmean_ngsr);
    dcellfree(simu->evlopdcov_ngsr);
    dcellfree(simu->evlopdmean_ngsr);
    zfarr_close_n(save->evlpsfmean, nevl);
    zfarr_close_n(save->evlopdcov, nevl);
    zfarr_close_n(save->evlopdmean, nevl);
    zfarr_close(save->evlopdcovol);
    zfarr_close(save->evlopdmeanol);
    zfarr_close_n(save->evlpsfmean_ngsr, nevl);
    zfarr_close_n(save->evlopdcov_ngsr, nevl);
    zfarr_close_n(save->evlopdmean_ngsr,nevl);
    zfarr_close(save->evlpsfolmean);
    dcellfree(simu->evlopd);    
    dfree(simu->lgsfocuslpf);
    cellfree(simu->ints);
    cellfree(simu->intsout);
    cellfree(simu->wfspsfout);
    cellfree(simu->pistatout);
    if(simu->dither){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(simu->dither[iwfs]){
		cellfree(simu->dither[iwfs]->imx);
		cellfree(simu->dither[iwfs]->imy);
		cellfree(simu->dither[iwfs]->imb);
		cellfree(simu->dither[iwfs]->i0);
		cellfree(simu->dither[iwfs]->gx);
		cellfree(simu->dither[iwfs]->gy);
		cellfree(simu->dither[iwfs]->ggm);
		cellfree(simu->dither[iwfs]->gg0);
		cellfree(simu->dither[iwfs]->mr);
		
		free(simu->dither[iwfs]);
	    }
	}
	free(simu->dither);
    }
    cellfree(simu->resdither);
    cellfree(simu->restwfs);
    cellfree(simu->zoompos);
    cellfree(simu->llt_tt);
    /*Close all files */
    
    zfarr_close_n(save->wfspsfout, nwfs);
    zfarr_close_n(save->ztiltout, nwfs);
    zfarr_close(save->dmerr);
    zfarr_close(save->dmint);
    zfarr_close(save->dmfit);
    zfarr_close(save->dmreal);    
    zfarr_close(save->dmcmd);
    dfree(save->ttmreal);
    zfarr_close(save->dmproj);
    zfarr_close(save->Merr_lo);
    zfarr_close(save->Mint_lo);
    zfarr_close(save->opdr);
    zfarr_close(save->opdx);
    zfarr_close_n(save->evlopdcl, nevl);
    zfarr_close_n(save->evlopdol, nevl);
    zfarr_close_n(save->ecovxx, nevl);
    zfarr_close_n(save->wfsopd, nwfs);
    zfarr_close_n(save->wfsopdol, nwfs);
    zfarr_close_n(save->wfslltopd, nwfs);
    zfarr_close_n(save->gradcl, nwfs);
    zfarr_close_n(save->gradnf, nwfs);
    zfarr_close_n(save->gradgeom, nwfs);
    zfarr_close_n(save->gradol, nwfs);
    zfarr_close_n(save->intsny, nwfs);
    zfarr_close_n(save->intsnf, nwfs);
    zfarr_close_n(save->dm_evl, nevl);
    zfarr_close_n(save->dm_wfs, nwfs);
 
    dfree(simu->winddir);
    if(parms->fdlock){
	/*release the lock and close the file. */
	char fn[80];
	char fnnew[80];
	snprintf(fn, 80, "Res_%d.lock",simu->seed);
	snprintf(fnnew, 80, "Res_%d.done",simu->seed);
	(void)rename(fn, fnnew);
	close(parms->fdlock->p[simu->iseed]);
    }
    free(simu->save);
    free(simu);
}
/**
   Print out wavefront error information and timing at each time step.
*/
void print_progress(SIM_T *simu){
    const int isim=simu->perfisim;
    const PARMS_T *parms=simu->parms;
    simu->status->wfs  =simu->tk_wfs;
    simu->status->recon=simu->tk_recon;
    simu->status->other=simu->tk_cache;
    simu->status->eval =simu->tk_eval;
    simu->status->scale=1;
    if(simu->timing){
	simu->timing->p[isim*simu->timing->nx]=get_job_mem();
	simu->timing->p[isim*simu->timing->nx+1]=simu->status->tot;
	simu->timing->p[isim*simu->timing->nx+2]=simu->status->wfs;
	simu->timing->p[isim*simu->timing->nx+3]=simu->status->recon;
	simu->timing->p[isim*simu->timing->nx+4]=simu->status->eval;
    }

    real this_time=myclockd();
    if(this_time>simu->last_report_time+1 || isim+1==parms->sim.end){
	/*we don't print out or report too frequently. */
	simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
	scheduler_report(simu->status);
#endif
	


	const STATUS_T *status=simu->status;
	const real tkmean=status->scale;
	const long rest=status->rest;
	const long laps=status->laps;
	const long resth=rest/3600;
	const long restm=(rest-resth*3600)/60;
	const long lapsh=laps/3600;
	const long lapsm=(laps-lapsh*3600)/60;
	const int nmod=parms->evl.nmod;
    
	if(parms->sim.evlol){
	    info("%sStep %5d: OL: %6.1f %6.1f %6.1f nm%s\n", GREEN,
		 isim,
		 mysqrt(simu->ole->p[isim*nmod])*1e9,
		 mysqrt(simu->ole->p[1+isim*nmod])*1e9,
		 mysqrt(simu->ole->p[2+isim*nmod])*1e9, BLACK);
	
	    info("Timing: Tot:%5.2f Mean:%5.2f Used %ld:%02ld Left %ld:%02ld\n",
		 status->tot*tkmean, status->mean*tkmean, lapsh,lapsm,resth,restm);
	}else{    
	    info("%sStep %5d: OL: %6.1f %6.1f %6.1f nm CL %6.1f %6.1f %6.1f nm",
		 GREEN, isim,
		 mysqrt(P(simu->ole, 0, isim))*1e9,
		 mysqrt(P(simu->ole, 1, isim))*1e9,
		 mysqrt(P(simu->ole, 2, isim))*1e9,
		 mysqrt(P(simu->cle, 0, isim))*1e9,
		 mysqrt(P(simu->cle, 1, isim))*1e9,
		 mysqrt(P(simu->cle, 2, isim))*1e9);
	    if(parms->recon.split){
		info(" Split %6.1f %6.1f %6.1f nm",
		     mysqrt(P(simu->clem, 0, isim))*1e9,
		     mysqrt(P(simu->clem, 1, isim))*1e9,
		     mysqrt(P(simu->clem, 2, isim))*1e9);
	    }
	    info("%s\n", BLACK);
    
	    info("Timing: WFS %.3f Recon %.3f EVAL %.3f Other %.3f Tot %.3f Mean %.3f."
		 " Used %ld:%02ld, Left %ld:%02ld\n",
		 status->wfs*tkmean, status->recon*tkmean, 
		 status->eval*tkmean, status->other*tkmean,
		 status->tot*tkmean, status->mean*tkmean,
		 lapsh,lapsm,resth,restm);
	    extern int NO_EVL;
	    if(!NO_EVL && !isfinite(simu->cle->p[isim*nmod])){
		error("Step %d: NaN/inf found: cle is %g\n", isim, simu->cle->p[isim*nmod]);
	    }
	}
    
	if(parms->plot.run && draw_current("Res", "Close loop")){
	    dmat *tmp=parms->recon.split?simu->clem:simu->cle;
	    const char *legs[4]={0,0,0,0};
	    int nline=2;
	    legs[0]="High Order";
	    legs[1]="Tip/Tilt";
	    if(parms->recon.split){
		if(simu->recon->ngsmod->indps){
		    legs[nline]="Plate Scale";
		    nline++;
		}else if(simu->recon->ngsmod->indastig){
		    legs[nline]="Astig";
		    nline++;
		}
		if(tmp->nx==4){
		    legs[nline]="Focus";
		    nline++;
		}
	    }
	    dcell *res=dcellnew_same(nline,1,simu->perfisim+1,1);
	    if(parms->recon.split){
		for(int i=0; i<=simu->perfisim; i++){
		    P(res->p[0], i)=P(tmp, 0, i);//LGS
		    P(res->p[1], i)=P(tmp, 1, i);//TT
		    if(nline>2){
			P(res->p[2], i)=P(tmp, 2, i)-P(tmp, 1, i);//PS
		    }
		    if(nline>3){
			P(res->p[3], i)=P(tmp, 3, i);//Focus
		    }
		}
	    }else{
		for(int i=0; i<simu->perfisim; i++){
		    P(res->p[0], i)=P(tmp, 2, i);//PTTR
		    P(res->p[1], i)=P(tmp, 0, i)-P(tmp,2,i);//TT
		}
	    }
	    dcellcwpow(res, 0.5); dcellscale(res, 1e9);
	    plot_points("Res", res->nx, NULL, res, NULL, NULL, "nn", NULL, legs,
			"Wavefront Error", "Time Step", "Wavefront Error (nm)", "Close loop");
	    dcellfree(res);
	}
    }
}
/**
   Output parameters necessary to run postproc using skyc/skyc.c
*/
void save_skyc(POWFS_T *powfs, RECON_T *recon, const PARMS_T *parms){
    char fn[PATH_MAX];
    real zadeg=parms->sim.za*180/M_PI;
    snprintf(fn,PATH_MAX,"%s/maos.conf",dirskysim);
    FILE *fp=fopen(fn,"w");
    fprintf(fp,"maos.r0z=%g\n",parms->atm.r0z);
    fprintf(fp,"maos.dt=%g\n",parms->sim.dt);
    fprintf(fp,"maos.zadeg=%g\n",zadeg);
    fprintf(fp,"maos.hc=%g\n",parms->dm[parms->ndm-1].ht);
    fprintf(fp,"maos.hs=%g\n",recon->ngsmod->hs);
    fprintf(fp,"maos.nmod=%d\n",recon->ngsmod->nmod);
    fprintf(fp,"maos.D=%g\n",parms->aper.d);
    fprintf(fp,"maos.wvl=[");
    int nwvl=0;
    int npowfs_ngs=0;
    int powfs_ngs[parms->npowfs];
    real ngsgrid=0;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].lo){
	    powfs_ngs[npowfs_ngs++]=ipowfs;
	    if(nwvl==0){
		nwvl=parms->powfs[ipowfs].nwvl;
	    }else{
		if(nwvl!=parms->powfs[ipowfs].nwvl){
		    error("Different low order WFS has different wvl\n");
		}
	    }
	    real sepmean=0;
	    for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs-1; indwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[indwfs];
		int iwfs2=parms->powfs[ipowfs].wfs->p[indwfs+1];
		if(fabs(parms->wfs[iwfs].thetax-parms->wfs[iwfs2].thetax)<1.e-10){
		    real sep=parms->wfs[iwfs2].thetay-parms->wfs[iwfs].thetay;
		    if(fabs(sepmean)<1e-20){
			sepmean=sep;
		    }else if(fabs(sep-sepmean)>1.e-10){
			sepmean=0;
			warning("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
		    }
		}
		if(fabs(parms->wfs[iwfs].thetay-parms->wfs[iwfs2].thetay)<1.e-10){
		    real sep=parms->wfs[iwfs2].thetax-parms->wfs[iwfs].thetax;
		    if(fabs(sepmean)<1e-20){
			sepmean=sep;
		    }else if(fabs(sep-sepmean)>1.e-10){
			sepmean=0;
			warning("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
		    }
		}
	    }
	    if(fabs(sepmean)>1.e-20){
		if(fabs(ngsgrid)<1.e-20){
		    ngsgrid=sepmean;
		}else{
		    if(fabs(sepmean-ngsgrid)>1.e-10){
			error("Different NGS POWFS have different spacing\n");
		    }
		}
	    }
	}
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	fprintf(fp,"%g ",parms->powfs[powfs_ngs[0]].wvl->p[iwvl]);
    }
    fprintf(fp,"]\n");  
    fprintf(fp,"maos.ngsgrid=%g\n",ngsgrid>0?ngsgrid*206265:1);
    fprintf(fp,"maos.npowfs=%d\n",npowfs_ngs);
    fprintf(fp,"maos.msa=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%d ",parms->powfs[powfs_ngs[ipowfs]].order);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.nsa=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%ld ",powfs[powfs_ngs[ipowfs]].saloc->nloc);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.ncomp=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%d ",powfs[powfs_ngs[ipowfs]].ncompx);
	if(powfs[powfs_ngs[ipowfs]].ncompx!=powfs[powfs_ngs[ipowfs]].ncompy){
	    error("Invalid\n");
	}
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.embfac=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%d ",parms->powfs[powfs_ngs[ipowfs]].embfac);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.dxsa=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%g ",powfs[powfs_ngs[ipowfs]].saloc->dx);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnwfsloc=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_loc\" " ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnwfsamp=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_amp\" " ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnsaloc=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_saloc\" " ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnmideal=\"RescleNGSm\"\n");
    fprintf(fp,"maos.fnmidealp=\"RescleNGSmp\"\n");
    fprintf(fp,"maos.evlindoa=%d\n",parms->evl.indoa);
    fprintf(fp,"maos.fnmcc=\"MCC_za%g\"\n",zadeg);
    fprintf(fp,"maos.fnmcc_oa=\"MCC_OA_za%g\"\n",zadeg);
    
    fprintf(fp,"maos.seeds=[");
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	fprintf(fp,"%ld " ,parms->sim.seeds->p[iseed]);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"include=\"skyc.conf\"\n");
    fprintf(fp,"include=\"skyc_za%g.conf\"\n",zadeg);
    fprintf(fp,"maos.wddeg=[");
    for(int ips=0; ips<parms->atm.nps; ips++){
	fprintf(fp, "%.2f ", parms->atm.wddeg->p[ips]);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.nstep=%d\n",parms->sim.end);
    fprintf(fp,"maos.ahstfocus=%d\n", parms->tomo.ahst_focus);
    fprintf(fp,"maos.mffocus=%d\n", parms->sim.mffocus);
    fprintf(fp,"maos.fnrange=%s\n", parms->powfs[parms->hipowfs->p[0]].llt->fnrange);
    fprintf(fp,"maos.indps=%d\n", recon->ngsmod->indps);
    fprintf(fp,"maos.indastig=%d\n", recon->ngsmod->indastig);
    fprintf(fp,"maos.indfocus=%d\n", recon->ngsmod->indfocus);
    fclose(fp);
    for(int jpowfs=0; jpowfs<npowfs_ngs; jpowfs++){
	int ipowfs=powfs_ngs[jpowfs];
	locwrite(powfs[ipowfs].loc,"%s/powfs%d_loc",dirskysim,ipowfs);
	writebin(powfs[ipowfs].amp, "%s/powfs%d_amp",dirskysim,ipowfs);
	locwrite(powfs[ipowfs].saloc,"%s/powfs%d_saloc",dirskysim,ipowfs);
	if(powfs[ipowfs].gradncpa){
	    int nsa=parms->powfs[ipowfs].order;
	    mymkdir("%s/gradoff/", dirskysim);
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		writebin(powfs[ipowfs].gradncpa->p[jwfs], 
		       "%s/gradoff/gradoff_sa%d_x%.0f_y%.0f.bin", 
		       dirskysim, nsa, 
		       parms->wfs[iwfs].thetax*206265, 
		       parms->wfs[iwfs].thetay*206265);
	    }
	}
    }
    writebin(recon->ngsmod->MCC,"%s/MCC_za%g", dirskysim,zadeg);
    writebin(recon->ngsmod->MCCP->p[parms->evl.indoa],"%s/MCC_OA_za%g", dirskysim,zadeg);
}
