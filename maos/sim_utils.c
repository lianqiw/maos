/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <sys/stat.h>
#include <sys/types.h>
#include "maos.h"
#include "sim.h"
#include "sim_utils.h"
#include "setup_surf.h"
#include "setup_powfs.h"
/**
   \file maos/sim_utils.c
   Contains a few support functions for simulation.
*/
/**
   Just do open loop error evalution. Usually not used.
   \callgraph
*/

void sim_evlol(const PARMS_T *parms,  POWFS_T *powfs, 
	       APER_T *aper,  RECON_T *recon){
    int simend=parms->sim.end;
    int simstart=parms->sim.start;
    if(simstart>=simend){
	maos_done(0);
	return;
    }
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
	if(!simu) continue;//skip
	genscreen(simu);
	simu->dt=parms->sim.dt;
	double tk_1=myclockd();
	for(int isim=simstart; isim<simend; isim++){
	    double ck_0=myclockd();
	    simu->isim=isim;
	    simu->status->isim=isim;
	    perfevl(simu);
	    double ck_1=myclockd();
	    
	    simu->status->rest=((ck_1-tk_1)*(double)(simend-isim-1)
				/(double)(isim-simstart+1));
	    simu->status->laps=(ck_1-tk_1);
	    
	    simu->status->eval=(double)(ck_1-ck_0);
	    simu->status->tot=(double)(ck_1-ck_0);
	    simu->status->scale=1;
	    simu->status->info=S_RUNNING;
#if defined(__linux__)
	    scheduler_report(simu->status);
#endif
	    print_progress(simu);
	}
	free_simu(simu);
    }
    maos_done(0);
}
/**
   wrap of the generic vonkarman_genscreen to generate turbulence screens. Wind
   velocities are set for each screen.  \callgraph */
void genscreen(SIM_T *simu){ 
    rand_t *rstat=simu->atm_rand;
    const PARMS_T *parms=simu->parms;
    const ATM_CFG_T *atm=&(simu->parms->atm);
    if(simu->atm){
	sqmaparrfree(simu->atm, parms->atm.nps);
    }
    if(simu->parms->dbg.noatm){
	warning("dbg.noatm flag is on. will not generate atmoshere\n");
	return;
    }
  
    map_t **screens=NULL;
    if(simu->parms->load.atm){
	const char *fn=simu->parms->load.atm;
	info2("loading atm from %s\n",fn);
	int nlayer;
	screens = sqmaparrread(&nlayer,"%s",fn);
	if(nlayer!=atm->nps)
	    error("Mismatch\n");
    }else{
	TIC;
	if(parms->dbg.atm==0){
	    info2("Generating Atmospheric Screen...");
	    tic;
	    screens = vonkarman_screen(rstat,atm->nx,atm->ny,atm->dx,atm->r0,
				       atm->l0,atm->wt,atm->nps,simu->nthread);
	    toc2("done");
	}else if(parms->dbg.atm==-1){
	    info2("Generating Biharmonic Atmospheric Screen...");
	    tic;
	    screens = biharmonic_screen(rstat,atm->nx,atm->ny,atm->dx,atm->r0,
					atm->l0,atm->wt,atm->nps,simu->nthread);
	    toc2("done");
	}else{
	    info2("Generating Testing Atmosphere Screen\n");
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
	    screens=calloc(atm->nps,sizeof(map_t*));
	    double hs=90000;
	    double dx=atm->dx;
	    for(int is=0; is<atm->nps; is++){
		screens[is]=calloc(1, sizeof(map_t));
		screens[is]->p=calloc(nx*ny,sizeof(double));
		screens[is]->nx=nx;
		screens[is]->ny=ny;
		screens[is]->dx=dx;
		screens[is]->ox=-nx/2*dx;
		screens[is]->oy=-ny/2*dx;
		screens[is]->h=atm->ht[is];
	    }
	    double scale=-pow(1.-screens[5]->h/hs,-2);
	    double strength=1./20626500.;
	    info("strength=%g, scale=%g\n",strength,scale);
	    switch(parms->dbg.atm){
	    case 1:{
		for(int iy=0; iy<ny; iy++){
		    double *p0=screens[0]->p+iy*nx;
		    for(int ix=0; ix<nx; ix++){
			double x=(ix-nx/2)*dx;
			p0[ix]=x*strength;
		    }
		}
	    }
		break;
	    case 2:{
		for(int iy=0; iy<ny; iy++){
		    double *p0=screens[0]->p+iy*nx;
		    double *p1=screens[5]->p+iy*nx;
		    double y=(iy-ny/2)*dx;
		    double yy=y*y;
		    for(int ix=0; ix<nx; ix++){
			double x=(ix-nx/2)*dx;
			double xx=x*x;
			//double xy=x*y;
			//p0[ix]=(iy-nx/2)*dx*strength;
			//p0[ix]=(x*0.2+y*0.1+xx*0.3-yy*0.7+xy*0.3)*strength;
			p0[ix]=(xx+yy)*strength;
			p1[ix]=scale*(p0[ix]);
		    }
		}
	    }
		break;
	    default:
		error("Invalid\n");
	    }
	}
	if(simu->parms->save.atm){
	    sqmaparrwrite(screens,atm->nps,"atm_%d.bin",simu->seed);
	}
    }
    int i;
    info2("Wind dir:");
    simu->winddir=dnew(atm->nps,1);
    int wdnz=0;
    for(i=0; i<atm->nps; i++){
	screens[i]->h=atm->ht[i];
	double angle;
	if(atm->wdrand){
	    if(fabs(atm->wddeg[i])>EPS){
		wdnz=1;
	    }
	    angle=randu(simu->atmwd_rand)*M_PI*2;
	}else{
	    angle=atm->wddeg[i]*M_PI/180;
	}
	simu->winddir->p[i]=angle;
	info2(" %5.1f", angle*180/M_PI);
	screens[i]->vx=cos(angle)*atm->ws[i];
	screens[i]->vy=sin(angle)*atm->ws[i];
    }
    info2(" deg\n");
    if(wdnz){
	error("wdrand is specified, but wddeg are not all zero. \n"
	      "possible confliction of intension!\n");
    }
    simu->atm=screens;
    info2("After genscreen:\t%.2f MiB\n",get_job_mem()/1024.);
    if(parms->plot.atm && simu->atm){
	for(int ips=0; ips<atm->nps; ips++){
	    drawmap("atm", simu->atm[ips],
		    "Atmosphere OPD","x (m)","y (m)","layer%d",ips);
	}
    }
    if(!parms->sim.frozenflow && parms->sim.closeloop){
	warning("Creating new screen in CL mode will not work\n");
	warning("Creating new screen in CL mode will not work\n");
	warning("Creating new screen in CL mode will not work\n");
    }
}

/**
   Propagate the atmosphere to closest xloc. skip wavefront sensing and
   reconstruction.
*/
dcell *atm2xloc(const SIM_T *simu){
    const RECON_T *recon=simu->recon;
    const PARMS_T *parms=simu->parms;
    if(!simu->atm)
	return NULL;
    dcell *opdx=dcellnew(recon->npsr,1);
    int isim=simu->isim;
    for(int ips=0; ips<parms->atm.nps; ips++){
	double disx=-simu->atm[ips]->vx*isim*simu->dt;
	double disy=-simu->atm[ips]->vy*isim*simu->dt;
	int ipsr=parms->atm.ipsr[ips];
	opdx->p[ipsr]=dnew(recon->xloc[ipsr]->nloc,1);
	prop_grid(simu->atm[ips],recon->xloc[ipsr],opdx->p[ipsr]->p,
		  1,disx,disy,1,1,0,0);
    }
    return opdx;
}

/**
   Evolving the Sodium layer by updating the elongation transfer function.
*/
void sim_update_etf(SIM_T *simu){
    int isim=simu->isim;
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    if(isim>0){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    // Update ETF if necessary.
	    if((parms->powfs[ipowfs].usephy
		||parms->powfs[ipowfs].psfout
		||parms->powfs[ipowfs].pistatout) 
	       && parms->powfs[ipowfs].hasllt
	       && parms->powfs[ipowfs].llt->colsimdtrat>0
	       && isim %parms->powfs[ipowfs].llt->colsimdtrat == 0){
		warning("powfs %d: Updating ETF\n",ipowfs);
		setup_powfs_etf(powfs,parms,ipowfs,1,
				isim/parms->powfs[ipowfs].llt->colsimdtrat);
	    }
	}
    }
}


/**
   Scale a dcell array and save to file.
*/
void dcell_mean_and_save(dcell *A, double scale, const char *format, ...){
    format2fn;
    dcell *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    dcelladd(&tmp, 0, A, scale);
    dcellwrite(tmp,"%s",fn);
    dcellfree(tmp);
}

/**
   Scale a dcell array and save to file.
*/
void dmat_mean_and_save(dmat *A, double scale, const char *format, ...){
    format2fn;
    dmat *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    dadd(&tmp, 0, A, scale);
    dwrite(tmp,"%s",fn);
    dfree(tmp);
}
/**
   use random number dirived from input seed to seed other stream.  necessary to
   have independant streams for different wfs in threading routines to avoid
   race condition and have consitent result */
void seeding(SIM_T *simu){
    info2("Running seed %d\n",simu->seed);
    simu->init=calloc(1, sizeof(rand_t));
    simu->atm_rand=calloc(1, sizeof(rand_t));
    simu->atmwd_rand=calloc(1, sizeof(rand_t));
    seed_rand(simu->init,simu->seed);
    seed_rand(simu->atm_rand,   lrand(simu->init));
    seed_rand(simu->atmwd_rand, lrand(simu->init));
    const PARMS_T *parms=simu->parms;
    simu->wfs_rand=calloc(parms->nwfs, sizeof(rand_t));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	seed_rand(&simu->wfs_rand[iwfs],lrand(simu->init));
    }
}
/**
   Initialize simu (of type SIM_T) and various simulation data structs. Called
   for every seed.
 */
SIM_T* init_simu(const PARMS_T *parms,POWFS_T *powfs, 
		 APER_T *aper,RECON_T *recon, int iseed){
    if(parms->fdlock[iseed]<0){
	warning("Another MAOS is already running. Skip seed %d\n", parms->sim.seeds[iseed]);
	return NULL;
    }
   
    SIM_T *simu=calloc(1, sizeof(SIM_T));
    simu->save=calloc(1, sizeof(SIM_SAVE_T));
    PINIT(simu->mutex_plot);
    PINIT(simu->mutex_wfsgrad);
    PINIT(simu->mutex_perfevl);
    PINIT(simu->mutex_cachedm);
    simu->parms=parms;
    simu->powfs=powfs;
    simu->recon=recon;
    simu->aper=aper;    
    simu->iseed=iseed;
    simu->seed=parms->sim.seeds[iseed];
    if(simu->seed==0){
	simu->seed=myclocki();
	warning("Seed is 0, using time as seed: %d\n",simu->seed);
    }
    int seed=simu->seed;
    simu->ints=calloc(parms->nwfs,sizeof(dcell*));
    simu->wfspsfout=calloc(parms->nwfs,sizeof(dcell*));
    simu->wfspsfoutcellarr=calloc(parms->nwfs,sizeof(cellarr*));
    simu->ztiltoutcellarr=calloc(parms->nwfs,sizeof(cellarr*));
    simu->pistatout=calloc(parms->nwfs,sizeof(dcell*));
    simu->sanea_sim=dcellnew(parms->nwfs,1);
    simu->gradcl=dcellnew(parms->nwfs,1);//output
    simu->gradpsol=dcellnew(parms->nwfs,1);//output
    simu->gradacc=dcellnew(parms->nwfs,1);//wfsgrad internal
    simu->nthread=parms->sim.nthread;
    if(parms->sim.servotype_lo==2){
	simu->gtypeII_lo=dread("%s",parms->sim.gtypeII_lo);
	simu->MtypeII_lo=calloc(1, sizeof(TYPEII_T));
	if(simu->gtypeII_lo->nx!=3){
	    error("%s has wrong format. should have 3 rows\n", parms->sim.gtypeII_lo);
	}
    }
    if(parms->sim.fuseint){
	simu->dmint=calloc(parms->sim.napdm, sizeof(dcell*));
    }else{
	simu->dmint_hi=calloc(parms->sim.napdm, sizeof(dcell*));
	simu->Mint_lo=calloc(parms->sim.napdm, sizeof(dcell*));
    }
    simu->uptint=calloc(parms->sim.napupt, sizeof(dcell*));
    if(parms->evl.psfmean){
	simu->evlpsfmean=dcellnew(parms->evl.nwvl,parms->evl.nevl);
	simu->evlpsfolmean=dcellnew(parms->evl.nwvl,1);
	if(parms->evl.tomo){
	    simu->evlpsftomomean=dcellnew(parms->evl.nwvl,parms->evl.nevl);
	}
    }
    simu->opdevl=dcellnew(parms->evl.nevl,1);
    if(parms->evl.psfhist){
	const int nevl=parms->evl.nevl;
	simu->evlpsfhist=calloc(nevl, sizeof(cellarr*));
	for(int ievl=0; ievl<nevl; ievl++){
	    if(!parms->evl.psf[ievl]) continue;
	    if(parms->evl.tomo!=2){//only evaluate tomography result.
		simu->evlpsfhist[ievl]=cellarr_init(parms->sim.end-parms->evl.psfisim, 
						    "evlpsfhist_%d_ievl%d.bin",
						    seed,ievl);
	    }
	    if(parms->evl.tomo){
		simu->evlpsftomohist[ievl]=cellarr_init(parms->sim.end-parms->evl.psfisim, 
							"evlpsftomohist_%d_ievl%d.bin",
							seed,ievl);
	    }
	}
    }
    if(parms->evl.psfmean || parms->evl.psfhist){
	//compute diffraction limited PSF.
	dmat *iopdevl=dnew(simu->aper->locs->nloc,1);
	ccell *psf2s=psfcomp(iopdevl, aper->amp, aper->embed, aper->nembed,
			     parms->evl.psfsize, parms->evl.nwvl, parms->evl.psfwvl);
	dfree(iopdevl);
	int nwvl=parms->evl.nwvl;
	dcell *evlpsfdl=dcellnew(nwvl,1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&evlpsfdl->p[iwvl], 1, psf2s->p[iwvl], 1);
	}
	ccellfree(psf2s);
	dcellwrite(evlpsfdl, "evlpsfdl_%d.bin",seed);
	dcellfree(evlpsfdl);
    }

    simu->has_upt=0;//flag for uplink tip/tilt control
    if(parms->sim.epfocus>1.e-15){
	if(parms->sim.lpfocus<1.e-15){
	    error("When epfocus is nonzero, lpfocus need to be zero\n");
	}
	warning("Creating simu->focus\n");
	simu->focusint=dcellnew(parms->nwfs,1);
	simu->focuslpf=dcellnew(parms->nwfs,1);
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].usephy && parms->powfs[ipowfs].llt){
	    //has LGS pointing loop.
	    simu->has_upt=1;
	    if(parms->powfs[ipowfs].trs && !recon->PTT) 
		//needed to form upterr from LGS gradients.
		error("PTT is not specified for an LGS!");
	    if(!simu->upterr){
		simu->upterr=dcellnew(parms->nwfs,1);
		simu->uptreal=dcellnew(parms->nwfs,1);
		break;
	    }
	}
    }
    {
	int moao_wfs=0;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao>-1){
		moao_wfs=1;
		break;
	    }
	}
	if(moao_wfs){
	    if(parms->sim.closeloop){
		simu->moao_wfs=dcellnew(3,parms->nwfs);
	    }else{
		simu->moao_wfs=dcellnew(1,parms->nwfs);
	    }
	}
	if(parms->evl.moao>-1){
	    if(parms->sim.closeloop){
		//we only need 2 here because perfevl is ahead of moao_recon
		simu->moao_evl=dcellnew(2,parms->evl.nevl);
	    }else{
		simu->moao_evl=dcellnew(1,parms->evl.nevl);
	    }
	}
    }
    simu->perfevl_iground=parms->atm.iground;
    if(parms->sim.cachedm)
	prep_cachedm(simu);

    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	simu->gradcl->p[iwfs]=dnew(nsa*2,1);
	if(parms->powfs[ipowfs].usephy){
	    simu->sanea_sim->p[iwfs]=dnew(nsa*2,1);
	}
    }
    simu->dtrat_hi=1;
    simu->dtrat_lo=1;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].dtrat>1){
	    if(parms->powfs[ipowfs].lo){
		if(simu->dtrat_lo==1){
		    simu->dtrat_lo=parms->powfs[ipowfs].dtrat;
		}else if(simu->dtrat_lo!=parms->powfs[ipowfs].dtrat){
		    error("We don't handle multiple framerate of the LO WFS yet\n");
		}
	    }else{
		warning("This is high order wfs, but dtrat=%d", parms->powfs[ipowfs].dtrat);
		if(simu->dtrat_hi==1){
		    simu->dtrat_hi=parms->powfs[ipowfs].dtrat;
		}else if(simu->dtrat_hi!=parms->powfs[ipowfs].dtrat){
		    error("We don't handle multiple framerate of the LO WFS yet\n");
		}
	    }
	}
    }

    //prepare data for ray tracing in wfsgrad.c
    simu->wfs_prop_dm=calloc(parms->nwfs*parms->ndm, sizeof(thread_t));
    simu->wfs_propdata_dm=calloc(parms->nwfs*parms->ndm, sizeof(PROPDATA_T));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	const int indwfs=parms->powfs[ipowfs].indwfs[iwfs];
	const double hs=parms->powfs[ipowfs].hs;
	int ilocm=-1;
	if(simu->powfs[ipowfs].locm){
	    ilocm=(simu->powfs[ipowfs].nlocm>1)?indwfs:0;
	    info("ilocm=%d\n",ilocm);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    const double ht = parms->dm[idm].ht;
	    PROPDATA_T *data=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
	    data->displacex=ht*parms->wfs[iwfs].thetax;
	    data->displacey=ht*parms->wfs[iwfs].thetay;
	    data->scale=1.-ht/hs;
	    data->alpha=-1;//remove dm contribution.
	    data->wrap=1;
	    if(simu->cachedm){
		int isc=parms->powfs[ipowfs].scalegroup[idm];
		data->mapin=&simu->cachedm[idm][isc];
		data->cubic=0;//already accounted for in cachedm.
		data->cubic_iac=0;//not needed
	    }else{
		data->locin=recon->aloc[idm];
		data->phiin=NULL;//replace later in simulation
		data->cubic=parms->dm[idm].cubic;
		data->cubic_iac=parms->dm[idm].iac;
	    }
	    data->phiout=NULL;//replace later in simulation
	    if(powfs[ipowfs].locm){//misregistration.
		data->locout=powfs[ipowfs].locm[ilocm];
	    }else{
		data->ptsout=powfs[ipowfs].pts;
	    }
	    thread_prep(&simu->wfs_prop_dm[iwfs+parms->nwfs*idm], 0, 0, 0, 1, data);
	}//idm
    }//iwfs

    if(parms->sim.frozenflow){
	simu->dt=parms->sim.dt;
    }else{
	simu->dt=0;
    }
    simu->dtlo=simu->dtrat_lo*simu->dt;
    simu->dthi=simu->dtrat_hi*simu->dt;
    //evaluation.
    int nsim=parms->sim.end;
    const int nevl=parms->evl.nevl;
    const int nmod=parms->evl.nmod;
    
    {//USE MMAP for data that need to save at every time step
	long nnx[nevl];
	long nny[nevl];
	for(int ievl=0; ievl<nevl; ievl++){
	    nnx[ievl]=nmod;
	    nny[ievl]=nsim;
	}
	simu->olep=dcellnew_mmap(nevl,1,nnx,nny,"Resolep_%d.bin",seed);
	simu->clep=dcellnew_mmap(nevl,1,nnx,nny,"Resclep_%d.bin",seed);
	simu->olmp=dcellnew_mmap(nevl,1,nnx,nny,"Resolmp_%d.bin",seed);
	simu->clmp=dcellnew_mmap(nevl,1,nnx,nny,"Resclmp_%d.bin",seed);
	if(parms->evl.tomo){
	    simu->cleptomo=dcellnew_mmap(nevl,1,nnx,nny,"Rescleptomo_%d.bin",seed); 
	    simu->clmptomo=dcellnew_mmap(nevl,1,nnx,nny,"Resclmptomo_%d.bin",seed);  
	}
	if(parms->tomo.split && parms->ndm<=2){
	    long nnx_split[nevl];
	    long nnx_3[nevl];
	    for(int ievl=0; ievl<nevl; ievl++){
		nnx_3[ievl]=3;
		nnx_split[ievl]=recon->ngsmod->nmod;
	    }
	    simu->clemp=dcellnew_mmap(nevl,1, nnx_3, nny, "Resclemp_%d.bin",seed);
	    simu->cleNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,"RescleNGSm_%d.bin",seed);
	    simu->cleNGSmp=dcellnew_mmap(nevl,1,nnx_split,nny,"RescleNGSmp_%d.bin",seed);
	    if(parms->tomo.split==1 && !parms->sim.fuseint){
		simu->corrNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,"RescorrNGSm_%d.bin",seed);
	    }
	    if(parms->sim.skysim){
		char fnold[PATH_MAX];
		char fnnew[PATH_MAX];
		snprintf(fnnew, PATH_MAX, "%s/RescleNGSm_%d.bin",dirskysim,seed);
		snprintf(fnold, PATH_MAX, "RescleNGSm_%d.bin", seed);
		if(link(fnold, fnnew)){
		    error("Error link\n");
		}
		snprintf(fnnew, PATH_MAX, "%s/RescleNGSmp_%d.bin",dirskysim,seed);
		snprintf(fnold, PATH_MAX, "RescleNGSmp_%d.bin", seed);
		if(link(fnold, fnnew)){
		    error("Error link\n");
		}
	    }
	}
    }
    {//MMAP the main result file
	long nnx[4]={nmod,nmod,nmod,nmod};
	long nny[4]={nsim,nsim,nsim,nsim};
	if(!parms->evl.tomo){
	    nnx[1]=0;
	    nny[1]=0;
	}
	if(!parms->tomo.split){
	    nnx[3]=0; 
	    nny[3]=0;
	}
	simu->res     = dcellnew_mmap(4,1,nnx,nny,"Res_%d.bin",seed);
	simu->ole     = dref(simu->res->p[0]);
	simu->cletomo = dref(simu->res->p[1]);
	simu->cle     = dref(simu->res->p[2]);
	simu->clem    = dref(simu->res->p[3]);
    }
    {//MMAP the LGS uptlink error/command output
	long nnx[parms->nwfs];
	long nny[parms->nwfs];
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].llt){
		nnx[iwfs]=2;
		nny[iwfs]=nsim;
	    }else{
		nnx[iwfs]=0;
		nny[iwfs]=0;
	    }
	}
	simu->upterrs = dcellnew_mmap(parms->nwfs, 1, nnx, nny, "Resupterr_%d.bin", seed);
	simu->uptcmds = dcellnew_mmap(parms->nwfs, 1, nnx, nny, "Resuptcmd_%d.bin", seed);

    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	long ncompx=powfs[ipowfs].ncompx;
	long ncompy=powfs[ipowfs].ncompy;
	if(parms->powfs[ipowfs].psfout){
	    const int nsa=powfs[ipowfs].pts->nsa;
	    //assert(parms->powfs[ipowfs].dtrat==1);
	    simu->wfspsfout[iwfs]=ccellnew(nsa,parms->powfs[ipowfs].nwvl);
	    for(long ipsf=0; ipsf<simu->wfspsfout[iwfs]->nx*simu->wfspsfout[iwfs]->ny; ipsf++){
		simu->wfspsfout[iwfs]->p[ipsf]=cnew(ncompx/2,ncompy/2);
	    }
	    mymkdir("%s/wvfout/", dirskysim);
	    mymkdir("%s/ztiltout/", dirskysim);
	    simu->wfspsfoutcellarr[iwfs]=cellarr_init
		(parms->sim.end-parms->sim.start,
		 "%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g.bin",
		 dirskysim, seed,
		 parms->powfs[ipowfs].order,
		 parms->wfs[iwfs].thetax*206265,
		 parms->wfs[iwfs].thetay*206265);
	
	    simu->ztiltoutcellarr[iwfs]=cellarr_init
		(parms->sim.end-parms->sim.start,
		 "%s/ztiltout/ztiltout_seed%d_sa%d_x%g_y%g.bin",
		 dirskysim, seed,
		 parms->powfs[ipowfs].order,
		 parms->wfs[iwfs].thetax*206265,
		 parms->wfs[iwfs].thetay*206265);

	}
	if(parms->powfs[ipowfs].pistatout){
	    mymkdir("%s/pistat/",dirskysim);
	}
    }
 
    if(parms->save.ngcov>0){
	simu->gcov=dcellnew(parms->save.ngcov,1);
    }
    SIM_SAVE_T *save=simu->save;
    int nstep=parms->sim.end-parms->sim.start;
    if(parms->save.dm){
	save->dmerr_hi=cellarr_init(nstep, "dmerr_hi_%d.bin", seed);
	save->dmfit_hi=cellarr_init(nstep, "dmfit_hi_%d.bin", seed);
	save->dmreal=cellarr_init(nstep, "dmreal_%d.bin", seed);
	if(parms->sim.fuseint){
	    save->dmint =cellarr_init(nstep, "dmint_%d.bin", seed);
	}else{
	    save->dmint_hi=cellarr_init(nstep, "dmint_hi_%d.bin", seed);
	}
	if(parms->tomo.split){
	    save->Merr_lo=cellarr_init(nstep, "Merr_lo_%d.bin", seed);
	    if(!parms->sim.fuseint){
		save->Mint_lo=cellarr_init(nstep, "Mint_lo_%d.bin", seed);
	    }
	}
	if(simu->moao_wfs){
	    save->moao_wfs=calloc(parms->nwfs, sizeof(cellarr*));
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int imoao=parms->powfs[ipowfs].moao;
		if(imoao>-1){
		    save->moao_wfs[iwfs]=cellarr_init(nstep,"wfs%d_moaofit_%d.bin",iwfs,seed);
		}
	    }
	}
	if(simu->moao_evl){
	    save->moao_evl=calloc(parms->nwfs, sizeof(cellarr*));
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		save->moao_evl[ievl]=cellarr_init(nstep, "evl%d_moaofit_%d.bin",ievl,seed);
	    }
	}
    }
    if(parms->save.dmpttr){
	save->dmpttr=cellarr_init(nstep, "dmpttr_%d.bin", seed);
    }

    if(parms->save.opdr){
	save->opdr=cellarr_init(nstep, "opdr_%d.bin", seed);
    }
    if(parms->save.opdx){
	save->opdx=cellarr_init(nstep, "opdx_%d.bin", seed);
    }
    if(parms->save.wfsopdhi || parms->save.wfsopdlo){
	save->wfsopd=calloc(parms->nwfs, sizeof(cellarr*));
	save->wfslltopd=calloc(parms->nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(!parms->save.powfs_opd[ipowfs]){
		continue;
	    }
	    save->wfsopd[iwfs]=cellarr_init(nstep, "wfs%d_opd_%d.bin", iwfs, seed);
	    if(powfs[ipowfs].lotf){
		save->wfslltopd[iwfs]=cellarr_init(nstep, "wfs%d_lltopd_%d.bin", iwfs, seed);
	    }
	}
    }
    if(parms->save.ints){
	save->intsny=calloc(parms->nwfs, sizeof(cellarr*));
	save->intsnf=calloc(parms->nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    int phystep=parms->powfs[ipowfs].phystep;
	    int nstep2=(parms->sim.end-phystep)/dtrat;
	    if(parms->save.powfs_ints[ipowfs]){
		save->intsny[iwfs]=cellarr_init(nstep2, "wfs%d_intsny_%d.bin", iwfs, seed);
		save->intsnf[iwfs]=cellarr_init(nstep2, "wfs%d_intsnf_%d.bin", iwfs, seed);
	    }
	}
    }
    if(parms->save.grad){
	save->gradcl=calloc(parms->nwfs, sizeof(cellarr*));
	save->gradnf=calloc(parms->nwfs, sizeof(cellarr*));
	save->gradpsol=calloc(parms->nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    if(parms->save.powfs_grad[ipowfs]){
		save->gradcl[iwfs]=cellarr_init(nstep/dtrat, "wfs%d_gradcl_%d.bin", iwfs, seed);
		save->gradnf[iwfs]=cellarr_init(nstep/dtrat, "wfs%d_gradnf_%d.bin", iwfs, seed);
		if(parms->tomo.split==2 || !parms->wfs[iwfs].skip){
		    save->gradpsol[iwfs]=cellarr_init(nstep/dtrat, "wfs%d_gradpsol_%d.bin", iwfs, seed);
		}
	    }
	}
    }
    if(parms->save.gradgeom){
	save->gradgeom=calloc(parms->nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    if(parms->save.powfs_gradgeom[ipowfs]){
		save->gradgeom[iwfs]=cellarr_init(nstep/dtrat, "wfs%d_gradgeom_%d.bin", iwfs, seed);
	    }
	}
    }
   
    if(parms->save.evlopd){
	save->evlopdol=calloc(parms->evl.nevl, sizeof(cellarr*));
	save->evlopdcl=calloc(parms->evl.nevl, sizeof(cellarr*));
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    save->evlopdol[ievl]=cellarr_init(nstep, "evl%d_opdol_%d.bin",ievl,seed);
	    save->evlopdcl[ievl]=cellarr_init(nstep, "evl%d_opdcl_%d.bin",ievl,seed);
	}
    }
    simu->dmpsol=calloc(parms->npowfs, sizeof(dcell*));
    simu->status=calloc(1, sizeof(STATUS_T));
    simu->status->iseed=iseed;
    simu->status->nseed=parms->sim.nseed;
    simu->status->simstart=parms->sim.start;
    simu->status->simend=parms->sim.end;
    simu->status->nthread=parms->sim.nthread;
    simu->status->timstart=myclocki();
    simu->status->info=S_RUNNING;
    seeding(simu);
    setup_tsurf(simu);//setting up M3 tilted surf.
    setup_surf(simu);//setting up M1/M2/M3 surface OPD.
    return simu;
}
/**
   Release memory of simu (of type SIM_T) and close files.
 */
void free_simu(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    free(simu->init);
    free(simu->atm_rand);
    free(simu->atmwd_rand);
    free(simu->wfs_rand);
    sqmaparrfree(simu->atm, parms->atm.nps);
    PDEINIT(simu->mutex_plot);
    PDEINIT(simu->mutex_wfsgrad);
    PDEINIT(simu->mutex_perfevl);
    PDEINIT(simu->mutex_cachedm);
    if(parms->sim.cachedm){
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int iscale=0; 
		iscale<parms->dm[idm].ncache; 
		iscale++){
		free(simu->cachedm[idm][iscale].p);
	    }
	    free(simu->cachedm[idm]);
	}
	free(simu->cachedm);
	free(simu->pcachedm);
	free(simu->cachedm_prop);
	free(simu->cachedm_propdata);
	
    }
    free(simu->wfs_prop_dm);
    free(simu->wfs_propdata_dm);
    free(simu->status);
    dcellfree(simu->gradcl);
    dcellfree(simu->gradacc);
    dcellfree(simu->gradpsol);
    dcellfree(simu->opdr);
    dcellfree(simu->opdrmvst);
    dcellfree(simu->opdevl);
    dcellfree(simu->dmreal);
    dcellfreearr(simu->dmpsol, parms->npowfs);
    dfree(simu->gtypeII_lo);
    if(parms->sim.fuseint){
	dcellfreearr(simu->dmint, parms->sim.napdm);
    }else{
	dcellfreearr(simu->dmint_hi, parms->sim.napdm);
	dcellfreearr(simu->Mint_lo, parms->sim.napngs);
    }
    dcellfree(simu->gcov);
    dcellfree(simu->dmerr_hi);
    dcellfree(simu->dmfit_hi);
    dcellfree(simu->dmhist);
    dcellfree(simu->Merr_lo);
    
    if(simu->MtypeII_lo){
	dcellfree(simu->MtypeII_lo->lead);
	dcellfree(simu->MtypeII_lo->firstint);
	dcellfree(simu->MtypeII_lo->errlast);
    }
    dcellfree(simu->upterr);
    dcellfree(simu->upterrlast);
    dcellfree(simu->uptreal);
    dcellfreearr(simu->uptint, parms->sim.napupt);

    dcellfree(simu->moao_wfs);
    dcellfree(simu->moao_evl);

    dcellfree(simu->res);
    dfree(simu->ole);
    dfree(simu->cle);
    dfree(simu->cletomo);
    dfree(simu->clem);
    dcellfree(simu->olep);
    dcellfree(simu->olmp);
    dcellfree(simu->clep);
    dcellfree(simu->clmp);
    dcellfree(simu->sanea_sim);
    dcellfree(simu->upterrs);
    dcellfree(simu->uptcmds);
    if(parms->tomo.split){
	dcellfree(simu->clemp);
	dfree(simu->cleNGSm);
	dfree(simu->corrNGSm);
	dcellfree(simu->cleNGSmp);
    }
    if(parms->evl.tomo){
	dcellfree(simu->cleptomo);
	dcellfree(simu->clmptomo);
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	if(simu->ints[iwfs])
	    dcellfree(simu->ints[iwfs]);
	if(simu->wfspsfout[iwfs])
	    ccellfree(simu->wfspsfout[iwfs]);
	if(simu->pistatout[iwfs])
	    dcellfree(simu->pistatout[iwfs]);
    }
    if(simu->evlpsfhist){
	const int nevl=parms->evl.nevl;
	for(int ievl=0; ievl<nevl; ievl++){
	    if(simu->evlpsfhist[ievl]){
		cellarr_close(simu->evlpsfhist[ievl]);
	    }
	    if(simu->evlpsftomohist[ievl]){
		cellarr_close(simu->evlpsftomohist[ievl]);
	    }
	}
	free(simu->evlpsfhist);
    }
    if(simu->evlpsfmean){
	dcellfree(simu->evlpsfmean);
	dcellfree(simu->evlpsfolmean);
	dcellfree(simu->evlpsftomomean);
    }
    free(simu->ints);
    free(simu->wfspsfout);
    free(simu->pistatout);
    //Close all files
    
    cellarr_close_n(simu->wfspsfoutcellarr, parms->nwfs);
    cellarr_close_n(simu->ztiltoutcellarr, parms->nwfs);
    SIM_SAVE_T *save=simu->save;
    cellarr_close(save->dmerr_hi);
    cellarr_close(save->dmint_hi);
    cellarr_close(save->dmfit_hi);
    cellarr_close(save->dmint);
    cellarr_close(save->dmpttr);
    cellarr_close(save->dmreal);
    cellarr_close(save->Merr_lo);
    cellarr_close(save->Mint_lo);
    cellarr_close(save->opdr);
    cellarr_close(save->opdx);
    cellarr_close_n(save->evlopdcl, parms->evl.nevl);
    cellarr_close_n(save->evlopdol, parms->evl.nevl);
    cellarr_close_n(save->wfsopd, parms->nwfs);
    cellarr_close_n(save->wfslltopd, parms->nwfs);
    cellarr_close_n(save->gradcl, parms->nwfs);
    cellarr_close_n(save->gradgeom, parms->nwfs);
    cellarr_close_n(save->gradnf, parms->nwfs);
    cellarr_close_n(save->gradpsol, parms->nwfs);
    cellarr_close_n(save->intsny, parms->nwfs);
    cellarr_close_n(save->intsnf, parms->nwfs);
    cellarr_close_n(save->moao_evl, parms->evl.nevl);
    cellarr_close_n(save->moao_wfs, parms->nwfs);

    dcellfree(simu->surfevl);
    dcellfree(simu->surfwfs);
    dfree(simu->winddir);
    dfree(simu->windest);
    spcellfree(simu->windshift);
    {
	//release the lock and close the file.
	close(parms->fdlock[simu->iseed]);
	char fn[80];
	snprintf(fn, 80, "Res_%d.lock",simu->seed);
	(void)remove(fn);
    }
    free(simu->save);
    free(simu);
}
/**
  Save telemetry data during simulation every 50 time steps.
*/
void save_simu(const SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    const int seed=simu->seed;
    if((isim % 50 ==0) || isim+1==parms->sim.end){
	double scale;
	if(parms->evl.psfmean && simu->isim>=parms->evl.psfisim){
	    if(parms->evl.tomo!=2){
		scale=1./(double)(simu->isim-parms->evl.psfisim+1);
		if(simu->evlpsfmean){
		    dcell_mean_and_save(simu->evlpsfmean, scale,
					"evlpsfcl_%d.bin",seed);
		}
		scale=1./(double)(simu->isim-parms->sim.start+1);
		if(parms->evl.psfol==2){
		    scale=scale/parms->evl.npsf;
		}
		if(parms->evl.psfol && simu->evlpsfolmean){
		    dcell_mean_and_save(simu->evlpsfolmean,scale,
					"evlpsfol_%d.bin",seed);
		}
	    }
	    if(parms->evl.tomo && simu->evlpsfmean){
		scale=1./(double)(simu->isim-parms->evl.psfisim+1);
		dcell_mean_and_save(simu->evlpsftomomean, scale,
				    "evlpsftomo_%d.bin",seed);
	    }
	}
	
	dcell *sanea=NULL;
	dcellcp(&sanea, simu->sanea_sim);
	for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	    const int ipowfs=simu->parms->wfs[iwfs].powfs;
	    const int dtrat=parms->powfs[ipowfs].dtrat;
	    if(sanea->p[iwfs] && simu->isim >=simu->parms->powfs[ipowfs].phystep){
		int nstep=(simu->isim+1-simu->parms->powfs[ipowfs].phystep)/dtrat;
		dscale(sanea->p[iwfs],1./nstep);
	    }
	}
	dcellwrite(sanea,"sanea_sim_%d.bin",seed);
	dcellfree(sanea);
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
    if(parms->save.ngcov>0 && ((simu->isim+1-parms->sim.start) % parms->save.gcovp) == 0){
	double scale=1./(double)(simu->isim-parms->sim.start+1);
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    dmat_mean_and_save(simu->gcov->p[igcov], scale, "gcov_wfs%d_%d_%d_%d.bin",
			       parms->save.gcov[igcov*2], parms->save.gcov[igcov*2+1],
			       simu->isim+1, seed);
	}
    }
}
/**
   Print out wavefront error information and timing at each time step.
 */
void print_progress(const SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const STATUS_T *status=simu->status;
    const double tkmean=status->scale;
    const long rest=status->rest;
    const long laps=status->laps;
    const long resth=rest/3600;
    const long restm=(rest-resth*3600)/60;
    const long lapsh=laps/3600;
    const long lapsm=(laps-lapsh*3600)/60;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
  
    
    
    fprintf(stderr,"\033[00;32mStep %3d: OL: %6.1f %6.1f %6.1f nm; CL %6.1f %6.1f %6.1f nm;",
	    isim,
	    mysqrt(simu->ole->p[isim*nmod])*1e9,
	    mysqrt(simu->ole->p[1+isim*nmod])*1e9,
	    mysqrt(simu->ole->p[2+isim*nmod])*1e9,
	    mysqrt(simu->cle->p[isim*nmod])*1e9,
	    mysqrt(simu->cle->p[1+isim*nmod])*1e9,
	    mysqrt(simu->cle->p[2+isim*nmod])*1e9);
    if(parms->tomo.split){
	fprintf(stderr," Split %6.1f %6.1f %6.1f nm;",
		mysqrt(simu->clem->p[isim*3])*1e9,
		mysqrt(simu->clem->p[1+isim*3])*1e9,
		mysqrt(simu->clem->p[2+isim*3])*1e9);
    }
    if(parms->evl.tomo){
	fprintf(stderr," TOMO %6.1f nm %6.1f nm;",
		mysqrt(simu->cletomo->p[isim*nmod])*1e9,
		mysqrt(simu->cletomo->p[2+isim*nmod])*1e9);
    }
    fprintf(stderr,"\033[00;00m\n");

    fprintf(stderr,"Timing: WFS:%5.2f Recon:%5.2f CACHE:%5.2f EVAL:%5.2f Tot:%5.2f Mean:%5.2f."
	    " Used %ld:%02ld, Left %ld:%02ld\n",
	    status->wfs*tkmean, status->recon*tkmean, 
	    status->cache*tkmean, status->eval*tkmean, 
	    status->tot*tkmean, status->mean*tkmean,
	    lapsh,lapsm,resth,restm);
}
/**
   Output parameters necessary to run postproc using skyc/skyc.c
 */
void save_skyc(POWFS_T *powfs, RECON_T *recon, const PARMS_T *parms){
    char fn[PATH_MAX];
    int zadeg=(int)round(parms->sim.za*180/M_PI);
    snprintf(fn,PATH_MAX,"%s/maos.conf",dirskysim);
    FILE *fp=fopen(fn,"w");
    fprintf(fp,"maos.r0z=%g\n",parms->atm.r0z);
    fprintf(fp,"maos.dt=%g\n",parms->sim.dt);
    fprintf(fp,"maos.zadeg=%d\n",zadeg);
    if(parms->ndm==2){
	fprintf(fp,"maos.hc=%g\n",parms->dm[1].ht);
    }else{
	error("Invalid");
    }
    fprintf(fp,"maos.hs=%g\n",recon->ngsmod->hs);
    fprintf(fp,"maos.nmod=%d\n",recon->ngsmod->nmod);
    fprintf(fp,"maos.D=%g\n",parms->aper.d);
    fprintf(fp,"maos.wvl=[");
    int nwvl=0;
    int npowfs_ngs=0;
    int powfs_ngs[parms->npowfs];
    double ngsgrid=0;
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
	    double sepmean=0;
	    for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs-1; indwfs++){
		int iwfs=parms->powfs[ipowfs].indwfs[indwfs];
		int iwfs2=parms->powfs[ipowfs].indwfs[indwfs+1];
		if(fabs(parms->wfs[iwfs].thetax-parms->wfs[iwfs2].thetax)<1.e-10){
		    double sep=parms->wfs[iwfs2].thetay-parms->wfs[iwfs].thetay;
		    if(fabs(sepmean)<1e-20){
			sepmean=sep;
		    }else if(fabs(sep-sepmean)>1.e-10){
			error("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
		    }
		}
		if(fabs(parms->wfs[iwfs].thetay-parms->wfs[iwfs2].thetay)<1.e-10){
		    double sep=parms->wfs[iwfs2].thetax-parms->wfs[iwfs].thetax;
		    if(fabs(sepmean)<1e-20){
			sepmean=sep;
		    }else if(fabs(sep-sepmean)>1.e-10){
			error("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
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
	fprintf(fp,"%g ",parms->powfs[powfs_ngs[0]].wvl[iwvl]);
    }
    fprintf(fp,"]\n");  
    fprintf(fp,"maos.ngsgrid=%g\n",ngsgrid*206265);
    fprintf(fp,"maos.npowfs=%d\n",npowfs_ngs);
    fprintf(fp,"maos.msa=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%d ",parms->powfs[powfs_ngs[ipowfs]].order);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.nsa=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"%ld ",powfs[powfs_ngs[ipowfs]].pts->nsa);
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
	fprintf(fp,"%g ",powfs[powfs_ngs[ipowfs]].pts->dsa);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnwfsloc=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_loc.bin.gz\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnwfsamp=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_amp.bin.gz\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnsaloc=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_saloc.bin.gz\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnmideal=\"RescleNGSm\"\n");
    fprintf(fp,"maos.fnmidealp=\"RescleNGSmp\"\n");
    fprintf(fp,"maos.evlindoa=%d\n",parms->evl.indoa);
    fprintf(fp,"maos.fnmcc=\"MCC_za%d.bin.gz\"\n",zadeg);
    fprintf(fp,"maos.fnmcc_oa=\"MCC_OA_za%d.bin.gz\"\n",zadeg);
    
    fprintf(fp,"maos.seeds=[");
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	fprintf(fp,"%d " ,parms->sim.seeds[iseed]);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"include=\"skyc.conf\"\n");
    fprintf(fp,"include=\"skyc_za%d.conf\"\n",zadeg);
    fprintf(fp,"maos.wddeg=[");
    for(int ips=0; ips<parms->atm.nps; ips++){
	fprintf(fp, "%.2f ", parms->atm.wddeg[ips]);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.nstep=%d\n",parms->sim.end);
    fclose(fp);
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	int jpowfs=powfs_ngs[ipowfs];
	locwrite(powfs[jpowfs].loc,"%s/powfs%d_loc.bin.gz",dirskysim,jpowfs);
	writedbl(powfs[jpowfs].amp,powfs[jpowfs].loc->nloc,1,
		 "%s/powfs%d_amp.bin.gz",dirskysim,jpowfs);
	locwrite(powfs[jpowfs].saloc,"%s/powfs%d_saloc.bin.gz",dirskysim,jpowfs);
    }
    dwrite(recon->ngsmod->MCC,"%s/MCC_za%d.bin.gz", dirskysim,zadeg);
    dwrite(recon->ngsmod->MCC_OA,"%s/MCC_OA_za%d.bin.gz", dirskysim,zadeg);
}
