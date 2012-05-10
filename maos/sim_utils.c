/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file maos/sim_utils.c
   Contains a few support functions for simulation.
*/
/*static double opdzlim[2]={-3e-5,3e-5}; */
static double *opdzlim=NULL;

static map_t **genscreen_do(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const ATM_CFG_T *atm=&parms->atm;
    const int nthread=simu->nthread;
    TIC;
    map_t **screens;
    if(parms->dbg.atm == 0 || parms->dbg.atm == -1){
	GENSCREEN_T *gs=simu->genscreen;
	if(!gs){
	    simu->genscreen=calloc(1, sizeof(GENSCREEN_T));/*the data for generating screens. */
	    gs=simu->genscreen;
	    gs->rstat  = simu->atm_rand;
	    gs->wt     = atm->wt;
	    gs->r0     = atm->r0;
	    gs->l0     = atm->l0;
	    gs->dx     = atm->dx;
	    gs->nx     = atm->nx;
	    gs->ny     = atm->ny;
	    gs->nlayer = atm->nps;
	    gs->ninit  = atm->ninit;
	    gs->share  = atm->share;
	    gs->nthread= nthread;
	}
	info2("Generating Atmospheric Screen...\n");
	tic;
	screens = parms->atm.fun(gs);
	toc2("Atmosphere ");
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
	    screens[is]=mapnew(nx, ny, dx, NULL);
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
		    /*double xy=x*y; */
		    /*p0[ix]=(iy-nx/2)*dx*strength; */
		    /*p0[ix]=(x*0.2+y*0.1+xx*0.3-yy*0.7+xy*0.3)*strength; */
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
    for(int i=0; i<atm->nps; i++){
	screens[i]->h=atm->ht[i];
	double angle=simu->winddir->p[i];
	screens[i]->vx=cos(angle)*atm->ws[i];
	screens[i]->vy=sin(angle)*atm->ws[i];
	/*misregistration */
	screens[i]->ox+=parms->aper.misreg[0];
	screens[i]->oy+=parms->aper.misreg[1];
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
	double wty=sa<0?1:0;
	double *p1=atm1->p+(1-(long)wty)*offy;
	double *p2=atm2->p+(long)wty*offy;
	double (*pp1)[nx]=(void*)p1;
	double (*pp2)[nx]=(void*)p2;
	double overyd=(double)overy;
	for(long iy=0; iy<overy; iy++){
	    double wt1=fabs(wty-(double)(iy+1)/overyd);
	    for(long ix=0; ix<nx; ix++){
		pp1[iy][ix]=(1-wt1)*pp1[iy][ix]+wt1*pp2[iy][ix];
		pp2[iy][ix]=pp1[iy][ix];
	    }
	}
    }else if(sa==0){
	rr=nx-overx;/*distance between the origins. */
	atm2->ox=atm1->ox + rr*ca*atm1->dx;
	atm2->oy=atm1->oy;
	double wtx=ca<0?1:0;
	double *p1=atm1->p+(1-(long)wtx)*offx;
	double *p2=atm2->p+(long)wtx*offx;
	double (*pp1)[nx]=(void*)p1;
	double (*pp2)[nx]=(void*)p2;
	double wts[overx];
	double overxd=(double)overx;
	for(long ix=0; ix<overx; ix++){
	    wts[ix]=fabs(wtx-(double)(ix+1)/overxd);
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
static void map_crop(map_t *atm, long overx, long overy){
    /*offset the atm to start from the corner and crop the screen to
      keep only the part that actually participates in ray tracing.*/
    
    double dx=atm->dx;
    long nx=atm->nx;
    long ny=atm->ny;
    long nxnew, nynew;
    double ox, oy;
    if(fabs(atm->vx)<EPS){/*along y */
	atm->vx=0;
	nxnew=overx;
	nynew=atm->ny;
	ox=-overx/2*dx;
	if(atm->vy<0){
	    oy=-overy/2*dx;
	}else{
	    oy=(overy/2-ny)*dx;
	}
    }else{
	atm->vy=0;
	nxnew=atm->nx;
	nynew=overy;
	oy=-overy/2*dx;
	if(atm->vx<0){
	    ox=-overx/2*dx;
	}else{
	    ox=(overx/2-nx)*dx;
	}
    }
    dmat *tmp=dsub((dmat*)atm, 0, nxnew, 0, nynew);
    free(atm->p);
    atm->p=tmp->p;
    atm->nx=tmp->nx;
    atm->ny=tmp->ny;
    atm->ox=ox;
    atm->oy=oy;
    dfree_keepdata(tmp);
}
/**
   wrap of the generic vonkarman_genscreen to generate turbulence screens. Wind
   velocities are set for each screen.  \callgraph */
void genscreen(SIM_T *simu){ 
    const PARMS_T *parms=simu->parms;
    const ATM_CFG_T *atm=&(simu->parms->atm);
    if(simu->atm){
	maparrfree(simu->atm, parms->atm.nps); simu->atm=NULL;
	dfree(simu->winddir);
        if(simu->atm2){
	    maparrfree(simu->atm2, parms->atm.nps);
	}
    }
    if(simu->parms->sim.noatm){
	warning("sim.noatm flag is on. will not generate atmoshere\n");
	return;
    }
    info2("Wind dir:");/*initialize wind direction one time only for each seed in frozen flow mode. */
    simu->winddir=dnew(atm->nps,1);
    int wdnz=0;
    for(int i=0; i<atm->nps; i++){
	double angle;
	if(atm->wdrand){
	    if(fabs(atm->wddeg[i])>EPS){
		wdnz=1;
	    }
	    angle=randu(simu->atmwd_rand)*M_PI*2;
	}else{
	    angle=atm->wddeg[i]*M_PI/180;
	}
	if(atm->evolve){
	    angle=round(angle*2/M_PI)*(M_PI/2);
	}
	simu->winddir->p[i]=angle;
	info2(" %5.1f", angle*180/M_PI);
    }
    info2(" deg\n");
    if(atm->evolve){
	warning("evolving screen requries direction to align along x/y.\n");
    }
    if(wdnz){
	error("wdrand is specified, but wddeg are not all zero. \n"
	      "possible confliction of intension!\n");
    }
    if(simu->parms->load.atm){
	const char *fn=simu->parms->load.atm;
	info2("loading atm from %s\n",fn);
	int nlayer;
	simu->atm = maparrread(&nlayer,"%s",fn);
	if(nlayer!=atm->nps)
	    error("Mismatch\n");
    }else{
	simu->atm=genscreen_do(simu);
    }
    if(simu->parms->save.atm){
	maparrwrite(simu->atm,atm->nps,"atm_%d.bin",simu->seed);
    }
    
    if(parms->plot.atm && simu->atm){
	for(int ips=0; ips<atm->nps; ips++){
	    drawmap("atm", simu->atm[ips],opdzlim,
		    "Atmosphere OPD","x (m)","y (m)","layer%d",ips);
	}
    }
    if(parms->atm.evolve){
	simu->atm2=genscreen_do(simu);
	for(int ips=0; ips<atm->nps; ips++){
	    long overx=parms->atm.overx[ips];
	    long overy=parms->atm.overy[ips];
	    map_crop(simu->atm[ips], overx, overy);
	    map_crop(simu->atm2[ips], overx, overy);
	    /*blend with the new screen. */
	    blend_screen_side(simu->atm[ips], simu->atm2[ips],
			      parms->atm.overx[ips],
			      parms->atm.overy[ips]);
	}

        if(parms->plot.atm && simu->atm){
	    for(int ips=0; ips<atm->nps; ips++){
		drawmap("atm1", simu->atm[ips],opdzlim,
			"Atmosphere OPD 1","x (m)","y (m)","layer%d",ips);
	    }
	    for(int ips=0; ips<atm->nps; ips++){
		drawmap("atm2", simu->atm2[ips],opdzlim,
			"Atmosphere OPD 2","x (m)","y (m)","layer%d",ips);
	    }
	}
    }
    info2("After genscreen:\t%.2f MiB\n",get_job_mem()/1024.);

    if(!parms->atm.frozenflow && parms->sim.closeloop){
	warning("Creating new screen in CL mode will not work\n");
	warning("Creating new screen in CL mode will not work\n");
	warning("Creating new screen in CL mode will not work\n");
    }
    if(simu->wfs_prop_atm){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    for(int ips=0; ips<parms->atm.nps; ips++){
		PROPDATA_T *data=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
		data->mapin=simu->atm[ips];
	    }
	}
    }
    if(simu->evl_prop_atm){
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    for(int ips=0; ips<parms->atm.nps; ips++){
		PROPDATA_T *data=&simu->evl_propdata_atm[ievl+parms->evl.nevl*ips];
		data->mapin=simu->atm[ips];
	    }
	}
    }
}
/**
   Evolve the turbulence screen when the ray goes near the edge by generating an
   new screen and blend into the old screen.  */
void evolve_screen(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const long isim=simu->isim;
    const double dt=parms->sim.dt;
    const int nps=parms->atm.nps;
    for(int ips=0; ips<nps; ips++){
	/*center of beam */
	double dx=simu->atm[ips]->dx;
	int do_evolve=0;
	if(fabs(simu->atm[ips]->vx)<EPS){/*along y */
	    long ny=simu->atm[ips]->ny;
	    long ry=(parms->atm.overy[ips]>>1);
	    double ay=-simu->atm[ips]->vy*isim*dt;
	    if(simu->atm[ips]->vy<0){
		if(ay>simu->atm[ips]->oy+(ny-ry)*dx){
		    do_evolve=1;
		}
	    }else{
		if(ay<simu->atm[ips]->oy+ry*dx){
		    do_evolve=1;
		}
	    }
	}else{/*along x. */
	    long nx=simu->atm[ips]->nx;
	    long rx=parms->atm.overx[ips]>>1;
	    double ax=-simu->atm[ips]->vx*isim*dt;
	    if(simu->atm[ips]->vx<0){
		if(ax>simu->atm[ips]->ox+(nx-rx)*dx){
		    do_evolve=1;
		}
	    }else{
		if(ax<simu->atm[ips]->ox+rx*dx){
		    do_evolve=1;
		}
	    }
	}
	if(do_evolve){
	    /*need to evolve screen. */
	    long overx=parms->atm.overx[ips];
	    long overy=parms->atm.overy[ips];
	    info("Evolving screen %d\n", ips);
	    mapfree(simu->atm[ips]);
	    simu->atm[ips]=simu->atm2[ips];
	    map_t **screen=parms->atm.fun(simu->genscreen);	
	    simu->atm2[ips]=screen[0]; 
	    simu->atm2[ips]->vx=simu->atm[ips]->vx;
	    simu->atm2[ips]->vy=simu->atm[ips]->vy;
	    simu->atm2[ips]->h=simu->atm[ips]->h;
	    free(screen);
	    map_crop(simu->atm2[ips], overx, overy);
	    blend_screen_side(simu->atm[ips], simu->atm2[ips], 
			      parms->atm.overx[ips],
			      parms->atm.overy[ips]);
	    if(simu->wfs_prop_atm){
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		    PROPDATA_T *data=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
		    data->mapin=simu->atm[ips];
		}
	    }
	    if(simu->evl_prop_atm){
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    PROPDATA_T *data=&simu->evl_propdata_atm[ievl+parms->evl.nevl*ips];
		    data->mapin=simu->atm[ips];
		}
	    }
	    if(parms->plot.atm){
		drawmap("atm1", simu->atm[ips],opdzlim,
			"Atmosphere OPD 1","x (m)","y (m)","layer%d_%d",ips,simu->isim);
		drawmap("atm2", simu->atm2[ips],opdzlim,
			"Atmosphere OPD 2","x (m)","y (m)","layer%d_%d",ips,simu->isim);
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
    warning("Generating Predictive HXW\n");
    PDSPCELL(recon->HXWtomo,HXWtomo);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){/*for tomography */
	    double  hs = parms->powfs[ipowfs].hs;
	    for(int ips=0; ips<npsr; ips++){
		spfree(HXWtomo[ips][iwfs]);
		double  ht = recon->ht->p[ips];
		double  scale=1. - ht/hs;
		double  displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht;
		displace[1]=parms->wfsr[iwfs].thetay*ht;
		if(parms->tomo.predict){
		    int ips0=parms->atmr.indps[ips];
		    displace[0]+=simu->atm[ips0]->vx*simu->dt*2;
		    displace[1]+=simu->atm[ips0]->vy*simu->dt*2;
		}
		HXWtomo[ips][iwfs]=mkh(recon->xloc[ips], ploc, NULL, 
				       displace[0],displace[1],scale,
				       parms->tomo.cubic, parms->tomo.iac);
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
    /*in close loop mode, opdr is from last time step. */
    int isim=parms->sim.closeloop?simu->isim-1:simu->isim;
    if(!*opdx){
	*opdx=dcellnew(recon->npsr,1);
    }
    for(int ipsr=0; ipsr<recon->npsr; ipsr++){
	if(!(*opdx)->p[ipsr]){
	    (*opdx)->p[ipsr]=dnew(recon->xloc[ipsr]->nloc,1);
	}else{
	    dzero((*opdx)->p[ipsr]);
	}
    }
    if(simu->atm){
	for(int ips=0; ips<parms->atm.nps; ips++){
	    double disx=-simu->atm[ips]->vx*isim*simu->dt;
	    double disy=-simu->atm[ips]->vy*isim*simu->dt;
	    int ipsr=parms->atm.ipsr[ips];
	    prop_grid(simu->atm[ips],recon->xloc[ipsr],NULL,(*opdx)->p[ipsr]->p,
		      1,disx,disy,1,1,0,0);
	}
    }
    if(recon->opdxadd){
	dcelladd(opdx, 1, recon->opdxadd, 1);
    }
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
	    /* Update ETF if necessary. */
	    if((parms->powfs[ipowfs].usephy
		||parms->powfs[ipowfs].psfout
		||parms->powfs[ipowfs].pistatout) 
	       && parms->powfs[ipowfs].llt
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
   use random number dirived from input seed to seed other stream.  necessary to
   have independant streams for different wfs in threading routines to avoid
   race condition and have consitent result */
void seeding(SIM_T *simu){
    info2("Running seed %d\n",simu->seed);
    simu->init=calloc(1, sizeof(rand_t));
    simu->atm_rand=calloc(1, sizeof(rand_t));
    simu->atmwd_rand=calloc(1, sizeof(rand_t));
    simu->telws_rand=calloc(1, sizeof(rand_t));
    seed_rand(simu->init,simu->seed);
    seed_rand(simu->atm_rand,   lrand(simu->init));
    /*2011-02-02: changed to wdrand-1 so that when wdrand=1, we reproduce old directions. */
    seed_rand(simu->atmwd_rand, lrand(simu->init)+(simu->parms->atm.wdrand-1));
    const PARMS_T *parms=simu->parms;
    simu->wfs_rand=calloc(parms->nwfs, sizeof(rand_t));
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	seed_rand(&simu->wfs_rand[iwfs],lrand(simu->init));
    }
    seed_rand(simu->telws_rand, lrand(simu->init)*parms->sim.wsseq);
#if USE_CUDA
    if(parms->gpu.wfs && !parms->sim.evlol){
	gpu_wfsgrad_seeding(parms,simu->powfs, simu->init);
    }
#endif
}

static void init_simu_evl(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    RECON_T *recon=simu->recon;
    const int nsim=parms->sim.end;
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    const int nmod=parms->evl.nmod;
    const int seed=simu->seed;
    SIM_SAVE_T *save=simu->save;
    simu->evlopd=dcellnew(nevl, 1);
    simu->perfevl_iground=parms->atm.iground;
    
    {/*MMAP the main result file */
	long nnx[4]={nmod,nmod,nmod,nmod};
	long nny[4]={nsim,nsim,nsim,nsim};
	nnx[1]=0;
	nny[1]=0;
	if(!parms->recon.split || parms->ndm>2){
	    nnx[3]=0; 
	    nny[3]=0;
	}
	simu->res     = dcellnew_mmap(4,1,nnx,nny,NULL,NULL,"Res_%d.bin",seed);
	simu->ole     = dref(simu->res->p[0]);
	simu->cle     = dref(simu->res->p[2]);
	simu->clem    = dref(simu->res->p[3]);
    }

    {/*USE MMAP for data that need to save at every time step */
	long nnx[nevl];
	long nny[nevl];
	for(int ievl=0; ievl<nevl; ievl++){
	    nnx[ievl]=nmod;
	    nny[ievl]=nsim;
	}
	simu->olep=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resolep_%d.bin",seed);
	simu->olmp=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resolmp_%d.bin",seed);
	simu->clep=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resclep_%d.bin",seed);
	simu->clmp=dcellnew_mmap(nevl,1,nnx,nny,NULL,NULL,"Resclmp_%d.bin",seed);
	if(!parms->sim.evlol){
	    if(parms->recon.split && parms->ndm<=2){
		long nnx_split[nevl];
		long nnx_3[nevl];
		for(int ievl=0; ievl<nevl; ievl++){
		    nnx_3[ievl]=3;
		    nnx_split[ievl]=recon->ngsmod->nmod;
		}
		simu->clemp=dcellnew_mmap(nevl,1, nnx_3, nny, NULL,NULL,"Resclemp_%d.bin",seed);
		simu->cleNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,NULL,"RescleNGSm_%d.bin",seed);
		simu->cleNGSmp=dcellnew_mmap(nevl,1,nnx_split,nny,NULL,NULL,"RescleNGSmp_%d.bin",seed);
		if(parms->recon.split==1 && !parms->sim.fuseint){
		    simu->corrNGSm=dnew_mmap(recon->ngsmod->nmod,nsim,NULL,"RescorrNGSm_%d.bin",seed);
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
	    }
	}
    }
        

    if(parms->evl.psfmean || parms->evl.psfhist || parms->evl.opdcov){
	char header[800];
	header[0]='\0';
	if(parms->evl.psfmean || parms->evl.psfhist ){
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		char headeri[80];
		snprintf(headeri, 80, "wvl=%g\nPSF sampling is %.15g radian\nPSF will sum to %.15g\n",
			 parms->evl.wvl[iwvl],
			 parms->evl.wvl[iwvl]/(aper->nembed[iwvl]*parms->evl.dx),
			 aper->sumamp2*aper->nembed[iwvl]*aper->nembed[iwvl]); 
		strncat(header, headeri, 800-strlen(header)-2);
	    }
	}
	long nframepsf=1;
	long nframecov=1;
	if(parms->evl.psfmean>1){
	    long nstep=(parms->sim.end-parms->evl.psfisim);
	    nframepsf=nstep/parms->evl.psfmean;
	    if(nstep > nframepsf*parms->evl.psfmean) nframepsf++;
	}
	if(parms->evl.opdcov>1){
	    long nstep=(parms->sim.end-parms->evl.psfisim);
	    nframecov=nstep/parms->evl.opdcov;
	    if(nstep > nframecov*parms->evl.opdcov) nframecov++;
	}
	char strht[24];
	if(parms->evl.psfmean){
	    simu->evlpsfmean=dcellnew(parms->evl.nwvl,nevl);
	    simu->evlpsfmean->header=strdup(header);
	    save->evlpsfmean=calloc(nevl, sizeof(cellarr*));
	    simu->evlpsfmean_ngsr=dcellnew(parms->evl.nwvl,nevl);
	    simu->evlpsfmean_ngsr->header=strdup(header);
	    save->evlpsfmean_ngsr=calloc(nevl, sizeof(cellarr*));
	}
	if(parms->evl.psfhist){
	    save->evlpsfhist=calloc(nevl, sizeof(cellarr*));
	    save->evlpsfhist_ngsr=calloc(nevl, sizeof(cellarr*));
	}
	if(parms->evl.opdcov){
	    simu->evlopdcov=dcellnew(nevl,1);
	    simu->evlopdcov_ngsr=dcellnew(nevl,1);
	    save->evlopdcov=calloc(nevl, sizeof(cellarr*));
	    save->evlopdcov_ngsr=calloc(nevl, sizeof(cellarr*));
	    simu->evlopdmean=dcellnew(nevl,1);
	    simu->evlopdmean_ngsr=dcellnew(nevl,1);
	    save->evlopdmean=calloc(nevl, sizeof(cellarr*));
	    save->evlopdmean_ngsr=calloc(nevl, sizeof(cellarr*));
	}
	for(int ievl=0; ievl<nevl; ievl++){
	    if(!parms->evl.psf[ievl]) continue;
	    if(isfinite(parms->evl.hs[ievl])){
		snprintf(strht, 24, "_%g", parms->evl.hs[ievl]);
	    }else{
		strht[0]='\0';
	    }
	    if(parms->evl.psfngsr[ievl]!=2){
		if(parms->evl.psfmean){
		    save->evlpsfmean[ievl]=cellarr_init(parms->evl.nwvl,nframepsf, 
							"evlpsfcl_%d_x%g_y%g%s.fits", seed,
							parms->evl.thetax[ievl]*206265,
							parms->evl.thetay[ievl]*206265, strht);
		}
		if(parms->evl.psfhist){
		    save->evlpsfhist[ievl]=cellarr_init(parms->sim.end-parms->evl.psfisim, 1,
							"evlpsfhist_%d_x%g_y%g%s.bin", seed,
							parms->evl.thetax[ievl]*206265,
							parms->evl.thetay[ievl]*206265,strht);
		}
		if(parms->evl.opdcov){
		    save->evlopdcov[ievl]=cellarr_init(nframecov, 1,
						       "evlopdcov_%d_x%g_y%g%s.bin", seed,
						       parms->evl.thetax[ievl]*206265,
						       parms->evl.thetay[ievl]*206265, strht);
		    save->evlopdmean[ievl]=cellarr_init(nframecov, 1,
							"evlopdmean_%d_x%g_y%g%s.bin", seed,
							parms->evl.thetax[ievl]*206265,
							parms->evl.thetay[ievl]*206265, strht);
		}
	    }
	    if(parms->evl.psfngsr[ievl]!=0){
		if(parms->evl.psfmean){
		    save->evlpsfmean_ngsr[ievl]=cellarr_init(parms->evl.nwvl,nframepsf, 
							     "evlpsfcl_ngsr_%d_x%g_y%g%s.fits", seed,
							     parms->evl.thetax[ievl]*206265,
							     parms->evl.thetay[ievl]*206265, strht);
		}
		if(parms->evl.psfhist){
		    save->evlpsfhist_ngsr[ievl]=cellarr_init(parms->sim.end-parms->evl.psfisim, 1,
							     "evlpsfhist_ngsr_%d_x%g_y%g%s.bin", seed,
							     parms->evl.thetax[ievl]*206265,
							     parms->evl.thetay[ievl]*206265,strht);
		}
		if(parms->evl.opdcov){
		    save->evlopdcov_ngsr[ievl]=cellarr_init(nframecov, 1,
							    "evlopdcov_ngsr_%d_x%g_y%g%s.bin", seed,
							    parms->evl.thetax[ievl]*206265,
							    parms->evl.thetay[ievl]*206265, strht);
		    save->evlopdmean_ngsr[ievl]=cellarr_init(nframecov, 1,
							     "evlopdmean_ngsr_%d_x%g_y%g%s.bin", seed,
							     parms->evl.thetax[ievl]*206265,
							     parms->evl.thetay[ievl]*206265, strht);
		}
	    }
	}/*for ievl */
	
	if(parms->evl.psfol){
	    if(parms->evl.psfmean){
		simu->evlpsfolmean=dcellnew(parms->evl.nwvl,1);
		simu->evlpsfolmean->header=strdup(header);
		save->evlpsfolmean=cellarr_init(parms->evl.nwvl, nframepsf, "evlpsfol_%d.fits", seed);
	    }
	    if(parms->evl.opdcov){
		save->evlopdcovol=cellarr_init(nframecov, 1, "evlopdcovol_%d", seed);
		save->evlopdmeanol=cellarr_init(nframecov, 1, "evlopdmean_%d", seed);
	    }
	}
    }
  
    if(parms->evl.psfmean || parms->evl.psfhist){
	/*compute diffraction limited PSF. */
	dmat *iopdevl=dnew(aper->locs->nloc,1);
	ccell *psf2s=psfcomp(iopdevl, aper->amp->p, aper->embed, aper->nembed,
			     parms->evl.psfsize, parms->evl.nwvl, parms->evl.wvl);
	dfree(iopdevl);
	dcell *evlpsfdl=dcellnew(nwvl,1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&evlpsfdl->p[iwvl], 1, psf2s->p[iwvl], 1);
	    evlpsfdl->p[iwvl]->header=evl_header(parms, aper, -1, iwvl);
	}
	ccellfree(psf2s);
	dcellwrite(evlpsfdl, "evlpsfdl_%d.fits",seed);
	dcellfree(evlpsfdl);
    }

    if(parms->sim.psfr){
	if(!parms->dbg.useopdr || parms->sim.idealfit){
	    simu->ecov=dcellnew(parms->ndm,parms->ndm);/*not really need. */
	}else{/*deprecated */
	    simu->ecov=dcellnew(nevl, 1);
	    if(parms->dbg.ecovxx){/*temporary. */
		warning("Saving history of ecov_xx\n");
		save->ecovxx=calloc(nevl, sizeof(cellarr *));
		char strht[24];
		for(int ievl=0; ievl<nevl; ievl++){
		    if(!parms->evl.psf[ievl]) continue;
		    if(isfinite(parms->evl.hs[ievl])){
			snprintf(strht, 24, "_%g", parms->evl.hs[ievl]);
		    }else{
			strht[0]='\0';
		    }
	    
		    save->ecovxx[ievl]=cellarr_init(parms->sim.end-parms->evl.psfisim-(parms->sim.closeloop?1:0), 1,
						    "ecovxx_%d_x%g_y%g%s.bin", seed,
						    parms->evl.thetax[ievl]*206265,
						    parms->evl.thetay[ievl]*206265,strht);
		}
	    }
	}
    }

    if(parms->save.evlopd){
	int nstep=parms->sim.end-parms->sim.start;
	save->evlopdol=calloc(nevl, sizeof(cellarr*));
	save->evlopdcl=calloc(nevl, sizeof(cellarr*));

	for(int ievl=0; ievl<nevl; ievl++){
	    save->evlopdol[ievl]=cellarr_init(nstep/parms->save.evlopd,1, 
					      "evl%d_opdol_%d.bin",ievl,seed);
	    save->evlopdcl[ievl]=cellarr_init(nstep/parms->save.evlopd,1,
					      "evl%d_opdcl_%d.bin",ievl,seed);
	}
    }

    /* For threading */
    simu->evl_prop_atm=calloc(nevl*parms->atm.nps, sizeof(thread_t*));
    simu->evl_propdata_atm=calloc(nevl*parms->atm.nps, sizeof(PROPDATA_T));
    simu->evl_prop_dm=calloc(nevl*parms->ndm, sizeof(thread_t*));
    simu->evl_propdata_dm=calloc(nevl*parms->ndm, sizeof(PROPDATA_T));
    for(int ievl=0; ievl<nevl; ievl++){
	const int nthread=parms->evl.nthread;
	int tot;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    const int ind=ievl+nevl*ips;
	    PROPDATA_T *data=&simu->evl_propdata_atm[ind];
	    const double ht=parms->atm.ht[ips];
	    data->displacex0=ht*parms->evl.thetax[ievl]+parms->evl.misreg[0];
	    data->displacey0=ht*parms->evl.thetay[ievl]+parms->evl.misreg[1];
	    data->scale=1-ht/parms->evl.hs[ievl];
	    data->alpha=1;
	    data->wrap=1;
	    data->mapin=(void*)1;/*need to update this in genscreen. */
	    data->phiout=(void*)1;/*replace later in simulation. */
	    data->ostat=aper->locs->stat;
	    prop_index(data);
	    tot=aper->locs->stat->ncol;
	    simu->evl_prop_atm[ind]=calloc(nthread, sizeof(thread_t));
	    thread_prep(simu->evl_prop_atm[ind], 0, tot, nthread, prop, data);
	}
	for(int idm=0; idm<parms->ndm && !parms->sim.evlol; idm++){
	    const int ind=ievl+nevl*idm;
	    PROPDATA_T *data=&simu->evl_propdata_dm[ind];
	    const double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
	    data->displacex0=ht*parms->evl.thetax[ievl]+parms->evl.misreg[0];
	    data->displacey0=ht*parms->evl.thetay[ievl]+parms->evl.misreg[1];
	    data->scale=1-ht/parms->evl.hs[ievl];
	    data->alpha=-1;
	    data->wrap=0;
	    if(simu->cachedm && parms->evl.scalegroup){
		int isc=parms->evl.scalegroup[idm+ievl*parms->ndm];
		data->mapin=simu->cachedm[idm][isc];
		data->cubic=0;/*already accounted for in cachedm. */
		data->cubic_iac=0;/*not needed */
		data->ostat=aper->locs->stat;
		tot=aper->locs->stat->ncol;
	    }else{
		if(simu->dmrealsq){
		    data->mapin=simu->dmrealsq[idm];
		}else{
		    data->locin=recon->alocm[idm];
		    data->phiin=simu->dmreal->p[idm]->p;
		}
		data->cubic=parms->dm[idm].cubic;
		data->cubic_iac=parms->dm[idm].iac;
		data->locout=aper->locs;/*propagate to locs if no cachedm. */
		tot=aper->locs->nloc;
	    }
	    data->phiout=(void*)1;/*replace later in simulation. */
	    prop_index(data);
	    simu->evl_prop_dm[ind]=calloc(nthread, sizeof(thread_t));
	    thread_prep(simu->evl_prop_dm[ind], 0, tot, nthread, prop, data);
	}
    }
}

static void init_simu_wfs(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    RECON_T *recon=simu->recon;
    SIM_SAVE_T *save=simu->save;
    const int nwfs=parms->nwfs;
    const int nsim=parms->sim.end;
    const int seed=simu->seed;
    simu->ints=calloc(nwfs,sizeof(dcell*));
    simu->wfspsfout=calloc(nwfs,sizeof(dcell*));
    save->wfspsfout=calloc(nwfs,sizeof(cellarr*));
    save->ztiltout =calloc(nwfs,sizeof(cellarr*));
    simu->pistatout=calloc(nwfs,sizeof(dcell*));
    simu->sanea_sim=dcellnew(nwfs, 1);
    simu->gradcl=dcellnew(nwfs,1);
    simu->gradnf=dcellnew(nwfs,1);
    /*Do not initialize gradlastcl. Do not initialize gradlastol in open
      loop. They are used for testing*/
    if(parms->sim.closeloop){
	simu->gradlastol=dcellnew(parms->nwfsr, 1);
    }
    simu->gradacc=dcellnew(nwfs,1);/*wfsgrad internal */

    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	simu->gradcl->p[iwfs]=dnew(nsa*2,1);
	if(parms->powfs[ipowfs].phytypesim==3 || parms->save.grad[iwfs]){
	    simu->gradnf->p[iwfs]=dnew(nsa*2,1);
	}
	if(parms->powfs[ipowfs].noisy){
	    simu->sanea_sim->p[iwfs]=dnew(4, nsa);
	}
	if(parms->powfs[ipowfs].usephy){
	    simu->ints[iwfs]=dcellnew(nsa,1);
	    for(int isa=0; isa<nsa; isa++){
		simu->ints[iwfs]->p[isa]=dnew(powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
	    }
	}
	if(parms->powfs[ipowfs].phystep!=0 || parms->save.gradgeom[iwfs] || parms->powfs[ipowfs].pistatout){
	    simu->gradacc->p[iwfs]=dnew(nsa*2,1);
	}
	if(parms->powfs[ipowfs].pistatout){
	    simu->pistatout[iwfs]=dcellnew(nsa,parms->powfs[ipowfs].nwvl);
	}
    }

    if(parms->sim.epfocus>1.e-15){
	if(parms->sim.lpfocus<1.e-15){
	    error("When epfocus is nonzero, lpfocus need to be zero\n");
	}
	warning("Creating simu->focus\n");
	simu->focusint=dcellnew(nwfs,1);
	simu->focuslpf=dcellnew(nwfs,1);
    }

    simu->has_upt=0;/*flag for uplink tip/tilt control */
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].usephy && parms->powfs[ipowfs].llt){
	    /*has LGS pointing loop. */
	    simu->has_upt=1;
	    if(parms->powfs[ipowfs].trs && !recon->PTT) 
		/*needed to form upterr from LGS gradients. */
		error("PTT is not specified for an LGS!");
	    if(!simu->upterr){
		simu->upterr=dcellnew(nwfs,1);
		simu->uptreal=dcellnew(nwfs,1);
		break;
	    }
	}
    }

    {/*MMAP the LGS uptlink error/command output */
	long nnx[nwfs];
	long nny[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].llt){
		nnx[iwfs]=2;
		nny[iwfs]=nsim;
	    }else{
		nnx[iwfs]=0;
		nny[iwfs]=0;
	    }
	}
	simu->upterrs = dcellnew_mmap(nwfs, 1, nnx, nny, NULL,NULL,"Resupterr_%d.bin", seed);
	simu->uptcmds = dcellnew_mmap(nwfs, 1, nnx, nny, NULL,NULL,"Resuptcmd_%d.bin", seed);
    }
   
    /* For sky coverage telemetry output */
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	long ncompx=powfs[ipowfs].ncompx;
	long ncompy=powfs[ipowfs].ncompy;
	long notf=MAX(ncompx, ncompy);
	if(parms->powfs[ipowfs].psfout){
	    const int nsa=powfs[ipowfs].pts->nsa;
	    /*The PSFs here are PSFs of each subaperture. */
	    simu->wfspsfout[iwfs]=ccellnew(nsa,parms->powfs[ipowfs].nwvl);
	    for(long ipsf=0; ipsf<simu->wfspsfout[iwfs]->nx*simu->wfspsfout[iwfs]->ny; ipsf++){
		simu->wfspsfout[iwfs]->p[ipsf]=cnew(notf/2,notf/2);
	    }
	    mymkdir("%s/wvfout/", dirskysim);
	    mymkdir("%s/ztiltout/", dirskysim);
	    save->wfspsfout[iwfs]=cellarr_init
		(parms->sim.end-parms->sim.start,1,
		 "%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g.bin",
		 dirskysim, seed,
		 parms->powfs[ipowfs].order,
		 parms->wfs[iwfs].thetax*206265,
		 parms->wfs[iwfs].thetay*206265);
	
	    save->ztiltout[iwfs]=cellarr_init
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
    int nstep=parms->sim.end-parms->sim.start;
    if(parms->save.wfsopd){
	save->wfsopd=calloc(nwfs, sizeof(cellarr*));
	save->wfsopdol=calloc(nwfs, sizeof(cellarr*));
	save->wfslltopd=calloc(nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(!parms->save.wfsopd[iwfs]){
		continue;
	    }
	    save->wfsopd[iwfs]=cellarr_init(nstep, 1,"wfs%d_opd_%d.bin", iwfs, seed);
	    save->wfsopdol[iwfs]=cellarr_init(nstep, 1,"wfs%d_opdol_%d.bin", iwfs, seed);
	    if(powfs[ipowfs].llt){
		save->wfslltopd[iwfs]=cellarr_init(nstep,1, "wfs%d_lltopd_%d.bin", iwfs, seed);
	    }
	}
    }
    if(parms->save.ints){
	save->intsny=calloc(nwfs, sizeof(cellarr*));
	save->intsnf=calloc(nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    int phystep=parms->powfs[ipowfs].phystep;
	    int nstep2=(parms->sim.end-phystep)/dtrat;
	    if(parms->save.ints[iwfs] && parms->powfs[ipowfs].usephy){
		save->intsny[iwfs]=cellarr_init(nstep2,1, "wfs%d_intsny_%d.bin", iwfs, seed);
		save->intsnf[iwfs]=cellarr_init(nstep2,1, "wfs%d_intsnf_%d.bin", iwfs, seed);
	    }
	}
    }
    if(parms->save.grad && !parms->sim.idealfit){
	save->gradcl=calloc(nwfs, sizeof(cellarr*));
	save->gradnf=calloc(nwfs, sizeof(cellarr*));
	save->gradol=calloc(nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    if(parms->save.grad[iwfs]){
		save->gradcl[iwfs]=cellarr_init(nstep/dtrat,1, "wfs%d_gradcl_%d.bin", iwfs, seed);
		if(parms->powfs[ipowfs].noisy)
		    save->gradnf[iwfs]=cellarr_init(nstep/dtrat,1, "wfs%d_gradnf_%d.bin", iwfs, seed);
		if(parms->recon.alg==0 &&(parms->recon.split==2 || !parms->powfs[ipowfs].skip)){
		    save->gradol[iwfs]=cellarr_init((nstep-(parms->sim.closeloop?1:0))/dtrat,1, 
						    "wfs%d_gradol_%d.bin", iwfs, seed);
		}
	    }
	}
    }
    if(parms->save.gradgeom && !parms->sim.idealfit){
	save->gradgeom=calloc(nwfs, sizeof(cellarr*));
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int dtrat=parms->powfs[ipowfs].dtrat;
	    if(parms->save.gradgeom[iwfs]){
		save->gradgeom[iwfs]=cellarr_init(nstep/dtrat,1, "wfs%d_gradgeom_%d.bin", iwfs, seed);
	    }
	}
    }
   
    /* threading */
    simu->wfs_prop_atm=calloc(nwfs*parms->atm.nps, sizeof(thread_t*));
    simu->wfs_propdata_atm=calloc(nwfs*parms->atm.nps, sizeof(PROPDATA_T));
    simu->wfs_prop_dm=calloc(nwfs*parms->ndm, sizeof(thread_t*));
    simu->wfs_propdata_dm=calloc(nwfs*parms->ndm, sizeof(PROPDATA_T));
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int nthread=powfs[ipowfs].nthread;
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	const double hs=parms->powfs[ipowfs].hs;
	const int ilocm=(simu->powfs[ipowfs].nlocm>1)?wfsind:0;

	for(int ips=0; ips<parms->atm.nps; ips++){
	    const double ht=parms->atm.ht[ips];
	    if(ht>hs){
		error("Layer is above guide star\n");
	    }
	    PROPDATA_T *data=&simu->wfs_propdata_atm[iwfs+nwfs*ips];
	    data->displacex0=ht*parms->wfs[iwfs].thetax+powfs[ipowfs].misreg[wfsind][0];
	    data->displacey0=ht*parms->wfs[iwfs].thetay+powfs[ipowfs].misreg[wfsind][1];
	    data->scale=1.-ht/hs;
	    data->alpha=1;
	    data->wrap=1;
	    data->mapin=(void*)1;/*need to update this in genscreen. */
	    data->phiout=(void*)1;/*replace later in simulation. */
	    int tot=0;
	    if(powfs[ipowfs].locm && powfs[ipowfs].locm[ilocm]){/*misregistration. */
		data->locout=powfs[ipowfs].locm[ilocm];
		tot=data->locout->nloc;
	    }else{
		data->ptsout=powfs[ipowfs].pts;
		tot=data->ptsout->nsa;
	    }
	    prop_index(data);
	    simu->wfs_prop_atm[iwfs+nwfs*ips]=calloc(nthread, sizeof(thread_t));
	    thread_prep(simu->wfs_prop_atm[iwfs+nwfs*ips],0,tot,nthread,prop,data);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	    PROPDATA_T *data=&simu->wfs_propdata_dm[iwfs+nwfs*idm];
	    int tot;
	    data->displacex0=ht*parms->wfs[iwfs].thetax+powfs[ipowfs].misreg[wfsind][0];
	    data->displacey0=ht*parms->wfs[iwfs].thetay+powfs[ipowfs].misreg[wfsind][1];
	    data->scale=1.-ht/hs;
	    data->alpha=-1;/*remove dm contribution. */
	    data->wrap=0;
	    if(simu->cachedm && parms->powfs[ipowfs].scalegroup[idm]!=-1){
		int isc=parms->powfs[ipowfs].scalegroup[idm];
		data->mapin=simu->cachedm[idm][isc];
		data->cubic=0;/*already accounted for in cachedm. */
		data->cubic_iac=0;/*not needed */
	    }else{
		if(simu->dmrealsq){
		    data->mapin=simu->dmrealsq[idm];
		}else{
		    data->locin=recon->alocm[idm];
		    data->phiin=simu->dmreal->p[idm]->p;
		}
		data->cubic=parms->dm[idm].cubic;
		data->cubic_iac=parms->dm[idm].iac;
	    }
	    data->phiout=(void*)1;/*replace later in simulation */
	    if(powfs[ipowfs].locm && powfs[ipowfs].locm[ilocm]){/*misregistration. */
		data->locout=powfs[ipowfs].locm[ilocm];
		tot=data->locout->nloc;
	    }else{
		data->ptsout=powfs[ipowfs].pts;
		tot=data->ptsout->nsa;
	    }
	    prop_index(data);
	    simu->wfs_prop_dm[iwfs+nwfs*idm]=calloc(nthread, sizeof(thread_t));
	    thread_prep(simu->wfs_prop_dm[iwfs+nwfs*idm], 0, tot, nthread, prop,data);
	}/*idm */
    }/*iwfs */
    simu->wfs_intsdata=calloc(nwfs, sizeof(WFSINTS_T));
    simu->wfs_ints=calloc(nwfs, sizeof(thread_t*));

    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int nthread=powfs[ipowfs].nthread;
	int tot=powfs[ipowfs].pts->nsa;
	WFSINTS_T *data=simu->wfs_intsdata+iwfs;
	data->iwfs=iwfs;
	data->parms=parms;
	data->powfs=powfs;
	simu->wfs_ints[iwfs]=calloc(parms->sim.nthread, sizeof(thread_t));
	thread_prep(simu->wfs_ints[iwfs], 0, tot, nthread, wfsints,data);
    }
}

static void init_simu_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    SIM_SAVE_T *save=simu->save;
    const int seed=simu->seed;
    /*we initialize dmreal, so that wfs_prop_dm can reference dmreal. */
    simu->dmcmd=dcellnew(parms->ndm,1);
    simu->dmreal=dcellnew(parms->ndm,1);
    simu->dmrealsq=calloc(parms->ndm,sizeof(map_t*));
    if(parms->sim.dmproj){
	simu->dmproj=dcellnew(parms->ndm,1);
	simu->dmprojsq=calloc(parms->ndm,sizeof(map_t*));
    }
    for(int idm=0; idm<parms->ndm; idm++){
	simu->dmcmd->p[idm]=dnew(recon->aloc[idm]->nloc,1);
	if(simu->hyst){
	    simu->dmreal->p[idm]=dnew(recon->aloc[idm]->nloc,1);
	}else{
	    simu->dmreal->p[idm]=dref(simu->dmcmd->p[idm]);
	}
	simu->dmrealsq[idm]=mapnew2(recon->amap[idm]);
	if(simu->dmprojsq){
	    simu->dmproj->p[idm]=dnew(recon->aloc[idm]->nloc,1);
	    simu->dmprojsq[idm]=mapnew2(recon->amap[idm]);
	}
	if(parms->fit.square){/*dmreal is already square.*/
	    free(simu->dmrealsq[idm]->p);
	    simu->dmrealsq[idm]->p=simu->dmreal->p[idm]->p;
	    simu->dmrealsq[idm]->nref=simu->dmreal->p[idm]->nref;
	    simu->dmrealsq[idm]->nref[0]++;
	    if(simu->dmprojsq){
		free(simu->dmprojsq[idm]->p);
		simu->dmprojsq[idm]->p=simu->dmproj->p[idm]->p;
		simu->dmprojsq[idm]->nref=simu->dmproj->p[idm]->nref;
		simu->dmprojsq[idm]->nref[0]++;
	    }
	}
    }
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	gpu_dmreal2gpu(simu->dmrealsq, parms->ndm, parms->dm);
	if(simu->dmprojsq){
	    gpu_dmproj2gpu(simu->dmprojsq, parms->ndm, parms->dm);
	}
    }
#endif
    simu->dmpsol=calloc(parms->npowfs, sizeof(dcell*));
    simu->dmint=servo_new(simu->dmreal, parms->sim.epdm);
    simu->Mint_lo=servo_new(NULL, parms->sim.eplo);
    simu->uptint=servo_new(NULL, parms->sim.epupt);

    /*Setup hysterisis */
    int anyhyst=0;
    simu->hyst = calloc(parms->ndm, sizeof(HYST_T*));
    for(int idm=0; idm<parms->ndm; idm++){
	if(parms->dm[idm].hyst){
	    simu->hyst[idm]=calloc(1, sizeof(HYST_T));
	    simu->hyst[idm]->coeff=dread("%s",parms->dm[idm].hyst);
	    int nhmod=simu->hyst[idm]->coeff->ny;
	    if(simu->hyst[idm]->coeff->nx!=3 || nhmod<1){
		error("DM hystereis file %s has wrong format. Expect 3 rows\n",parms->dm[idm].hyst);
	    }
	    int naloc=recon->aloc[idm]->nloc;
	    simu->hyst[idm]->xlast=dnew(naloc,1);
	    simu->hyst[idm]->ylast=dnew(nhmod,naloc);
	    simu->hyst[idm]->dxlast=dnew(naloc,1);
	    simu->hyst[idm]->x0=dnew(naloc,1);
	    simu->hyst[idm]->y0=dnew(nhmod,naloc);
	    anyhyst=1;
	}
    }
    if(!anyhyst){
	free(simu->hyst); simu->hyst=NULL;
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
		    nny[idm]=recon->aloc[idm]->nloc;
		}else{
		    nnx[idm]=0;
		    nny[idm]=0;
		}
	    }
	    simu->dmhist=dcellnew_mmap(parms->ndm, 1, nnx, nny, 
				       NULL, NULL,"dmhist_%d.bin",simu->seed);
	}
    }
    int nstep=parms->sim.end-parms->sim.start;
    if(parms->save.dm){
	int nrstep=nstep-(parms->sim.closeloop?1:0);
	save->dmerr=cellarr_init(nrstep, 1,"dmerr_%d.bin", seed);
	if(parms->recon.alg==0){
	    save->dmfit=cellarr_init(nrstep, 1, "dmfit_%d.bin", seed);
	}
	if(parms->recon.split){
	    save->Merr_lo=cellarr_init(nrstep, 1, "Merr_lo_%d.bin", seed);
	}	
	save->dmreal = cellarr_init(nstep, 1, "dmreal_%d.bin", seed);
	save->dmcmd  = cellarr_init(nstep, 1, "dmcmd_%d.bin", seed);
	if(parms->sim.dmproj){
	    save->dmproj = cellarr_init(nstep, 1, "dmproj_%d.bin", seed);
	}

    }
  
    if(parms->save.dmpttr){
	save->dmpttr=cellarr_init(nstep, 1,"dmpttr_%d.bin", seed);
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
    int nstep=parms->sim.end-parms->sim.start;
    int ny=parms->sim.closeloop && !parms->gpu.moao ? 2 : 1;
    if(!parms->gpu.wfs || ! parms->gpu.moao){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao>-1){
		if(!simu->dm_wfs){
		    simu->dm_wfs=dcellnew(nwfs, ny);
		}
		for(int iy=0; iy<ny; iy++){
		    simu->dm_wfs->p[iwfs+iy*nwfs]=dnew(recon->moao[imoao].aloc->nloc, 1);
		}
	    }
	}
    }
    if(!parms->gpu.evl || ! parms->gpu.moao){ 
	int imoao=parms->evl.moao;
	if(imoao!=-1){
	    simu->dm_evl=dcellnew(nevl, ny);
	    for(int ievl=0; ievl<nevl; ievl++){
		for(int iy=0; iy<ny; iy++){
		    simu->dm_evl->p[ievl+iy*nevl]=dnew(recon->moao[imoao].aloc->nloc, 1);
		}
	    }
	}
    }
#if USE_CUDA
    if(!parms->gpu.moao){
	gpu_moao_2gpu(simu);//initilization.
    }
#endif
    if(parms->save.dm){
	if(simu->dm_wfs){
	    save->dm_wfs=calloc(nwfs, sizeof(cellarr*));
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int imoao=parms->powfs[ipowfs].moao;
		if(imoao>-1){
		    save->dm_wfs[iwfs]=cellarr_init(nstep,1,"wfs%d_moaofit_%d.bin",iwfs,seed);
		}
	    }
	}
	if(simu->dm_evl){
	    save->dm_evl=calloc(nwfs, sizeof(cellarr*));
	    for(int ievl=0; ievl<nevl; ievl++){
		save->dm_evl[ievl]=cellarr_init(nstep,1, "evl%d_moaofit_%d.bin",ievl,seed);
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
    if(parms->fdlock[iseed]<0){
	warning("Another MAOS is already running. Skip seed %d\n", parms->sim.seeds[iseed]);
	return NULL;
    }
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    SIM_T *simu=calloc(1, sizeof(SIM_T));
    SIM_SAVE_T *save=simu->save=calloc(1, sizeof(SIM_SAVE_T));
  
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
    seeding(simu);

    simu->status=calloc(1, sizeof(STATUS_T));
    simu->status->iseed=iseed;
    simu->status->nseed=parms->sim.nseed;
    simu->status->simstart=parms->sim.start;
    simu->status->simend=parms->sim.end;
    simu->status->nthread=parms->sim.nthread;
    simu->status->timstart=myclocki();
    simu->status->info=S_RUNNING;

    {   /* dtrat */
	if(parms->atm.frozenflow){
	    simu->dt=parms->sim.dt;
	}else{
	    simu->dt=0;
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
		    if(simu->dtrat_hi==1){
			simu->dtrat_hi=parms->powfs[ipowfs].dtrat;
		    }else if(simu->dtrat_hi!=parms->powfs[ipowfs].dtrat){
			error("We don't handle multiple framerate of the LO WFS yet\n");
		    }
		}
	    }
	}
	simu->dtlo=simu->dtrat_lo*simu->dt;
	simu->dthi=simu->dtrat_hi*simu->dt;
    }

    if(parms->sim.wspsd){
	/* Telescope wind shake added to TT input. */
	dmat *psdin=dread("%s", parms->sim.wspsd);
	info2("Loading windshake PSD from file %s\n", parms->sim.wspsd);
	simu->telws = psd2time(psdin, simu->telws_rand, parms->sim.dt, parms->sim.end);
	dfree(psdin);
	dwrite(simu->telws, "telws_%d", seed);
    }

    /* Select GPU or CPU for the tasks.*/
    simu->wfs_grad=calloc(nwfs, sizeof(thread_t));
#if USE_CUDA
    if(parms->gpu.wfs){
	thread_prep(simu->wfs_grad, 0, nwfs, nwfs, gpu_wfsgrad, simu);
    }else{
#endif
	thread_prep(simu->wfs_grad, 0, nwfs, nwfs, wfsgrad_iwfs, simu);
#if USE_CUDA
    }
#endif

    simu->perf_evl=calloc(nevl, sizeof(thread_t));
#if USE_CUDA
    if(parms->gpu.evl){
	thread_prep(simu->perf_evl, 0, nevl, nevl, gpu_perfevl, simu);
    }else{
#endif
	thread_prep(simu->perf_evl, 0, nevl, nevl, perfevl_ievl, simu);
#if USE_CUDA
    }
#endif
    simu->nthread=parms->sim.nthread;

    if(!parms->sim.evlol){
	init_simu_dm(simu);
	init_simu_moao(simu);
	init_simu_wfs(simu);
	if(parms->sim.cachedm){
	    prep_cachedm(simu);
	}
	if(parms->recon.alg==0){
	    int nstep=parms->sim.end-parms->sim.start;
	    if(parms->save.opdr){
		save->opdr=cellarr_init(nstep,1, "opdr_%d.bin", seed);
	    }
	    if(parms->save.opdx){
		save->opdx=cellarr_init(nstep, 1,"opdx_%d.bin", seed);
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
    return simu;
}
/**
   Release memory of simu (of type SIM_T) and close files.
*/
void free_simu(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    free(simu->init);
    free(simu->atm_rand);
    free(simu->atmwd_rand);
    free(simu->wfs_rand);
    free(simu->telws_rand);
    dfree(simu->telws);
    if(simu->genscreen){
	dfree(simu->genscreen->spect);
	free(simu->genscreen);
    }
    maparrfree(simu->atm, parms->atm.nps);
    maparrfree(simu->atm2, parms->atm.nps);

    if(parms->sim.cachedm && !parms->sim.evlol){
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int iscale=0; 
		iscale<parms->dm[idm].ncache; 
		iscale++){
		mapfree(simu->cachedm[idm][iscale]);
	    }
	    free(simu->cachedm[idm]);
	}
	free(simu->cachedm);
	free(simu->pcachedm);
	for(int i=0; i<simu->cachedm_n; i++){
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
	    if(simu->hyst[idm]){
		dfree(simu->hyst[idm]->coeff);
		dfree(simu->hyst[idm]->xlast);
		dfree(simu->hyst[idm]->ylast);
		dfree(simu->hyst[idm]->dxlast);
		dfree(simu->hyst[idm]->x0);
		dfree(simu->hyst[idm]->y0);
	    }
	    free(simu->hyst[idm]);
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
    free(simu->wfs_grad);
    free(simu->perf_evl);
    free(simu->status);
    dcellfree(simu->gradcl);
    dcellfree(simu->gradnf);
    dcellfree(simu->gradacc);
    dcellfree(simu->gradlastcl);
    dcellfree(simu->gradlastol);
    dcellfree(simu->opdr);
    dcellfree(simu->opdrmvst);
    dcellfree(simu->dmreal);
    maparrfree(simu->dmrealsq, parms->ndm);
    dcellfree(simu->dmproj);
    dcellfreearr(simu->dmpsol, parms->npowfs);
    dcellfree(simu->dmcmd);
    dcellfree(simu->dmcmdlast);
    servo_free(simu->dmint);
    servo_free(simu->Mint_lo);
    dcellfree(simu->gcov);
    dcellfree(simu->ecov);
    dcellfree(simu->dmerr);
    dcellfree(simu->dmfit);
    dcellfree(simu->dmhist);
    dcellfree(simu->Merr_lo);
    dcellfree(simu->Merr_lo_keep);
    dcellfree(simu->upterr);
    dcellfree(simu->uptreal);
    servo_free(simu->uptint);

    dcellfree(simu->dm_wfs);
    dcellfree(simu->dm_evl);

    dcellfree(simu->res);
    dfree(simu->ole);
    dfree(simu->cle);
    dfree(simu->clem);
    dcellfree(simu->olep);
    dcellfree(simu->olmp);
    dcellfree(simu->clep);
    dcellfree(simu->clmp);
    dcellfree(simu->sanea_sim);
    dcellfree(simu->upterrs);
    dcellfree(simu->uptcmds);
    if(parms->recon.split){
	dcellfree(simu->clemp);
	dfree(simu->cleNGSm);
	dfree(simu->corrNGSm);
	dcellfree(simu->cleNGSmp);
    }
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(simu->ints[iwfs])
	    dcellfree(simu->ints[iwfs]);
	if(simu->wfspsfout[iwfs])
	    ccellfree(simu->wfspsfout[iwfs]);
	if(simu->pistatout[iwfs])
	    dcellfree(simu->pistatout[iwfs]);
    }
    SIM_SAVE_T *save=simu->save;
    cellarr_close_n(save->evlpsfhist, nevl);
    cellarr_close_n(save->evlpsfhist_ngsr, nevl);
    dcellfree(simu->evlpsfmean);
    dcellfree(simu->evlpsfolmean);
    dcellfree(simu->evlopdcov);
    dcellfree(simu->evlopdmean);
    dfree(simu->evlopdcovol);
    dfree(simu->evlopdmeanol);
    dcellfree(simu->evlpsfmean_ngsr);
    dcellfree(simu->evlopdcov_ngsr);
    dcellfree(simu->evlopdmean_ngsr);
    cellarr_close_n(save->evlpsfmean, nevl);
    cellarr_close_n(save->evlopdcov, nevl);
    cellarr_close_n(save->evlopdmean, nevl);
    cellarr_close(save->evlopdcovol);
    cellarr_close(save->evlopdmeanol);
    cellarr_close_n(save->evlpsfmean_ngsr, nevl);
    cellarr_close_n(save->evlopdcov_ngsr, nevl);
    cellarr_close_n(save->evlopdmean_ngsr,nevl);
    cellarr_close(save->evlpsfolmean);
    dcellfree(simu->evlopd);    
    free(simu->ints);
    free(simu->wfspsfout);
    free(simu->pistatout);
    /*Close all files */
    
    cellarr_close_n(save->wfspsfout, nwfs);
    cellarr_close_n(save->ztiltout, nwfs);

    cellarr_close(save->dmerr);
    cellarr_close(save->dmfit);
    cellarr_close(save->dmpttr);
    cellarr_close(save->dmreal);
    cellarr_close(save->dmproj);
    cellarr_close(save->Merr_lo);
    cellarr_close(save->opdr);
    cellarr_close(save->opdx);
    cellarr_close_n(save->evlopdcl, nevl);
    cellarr_close_n(save->evlopdol, nevl);
    cellarr_close_n(save->ecovxx, nevl);
    cellarr_close_n(save->wfsopd, nwfs);
    cellarr_close_n(save->wfsopdol, nwfs);
    cellarr_close_n(save->wfslltopd, nwfs);
    cellarr_close_n(save->gradcl, nwfs);
    cellarr_close_n(save->gradgeom, nwfs);
    cellarr_close_n(save->gradnf, nwfs);
    cellarr_close_n(save->gradol, nwfs);
    cellarr_close_n(save->intsny, nwfs);
    cellarr_close_n(save->intsnf, nwfs);
    cellarr_close_n(save->dm_evl, nevl);
    cellarr_close_n(save->dm_wfs, nwfs);
 
    dfree(simu->winddir);
    {
	/*release the lock and close the file. */
	close(parms->fdlock[simu->iseed]);
	char fn[80];
	char fnnew[80];
	snprintf(fn, 80, "Res_%d.lock",simu->seed);
	snprintf(fnnew, 80, "Res_%d.done",simu->seed);
	(void)rename(fn, fnnew);
    }
    free(simu->save);
    free(simu);
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
    
    if(parms->sim.evlol){
	info2("\033[00;32mStep %5d: OL: %6.1f %6.1f %6.1f nm;\033[00;00m\n",
	      isim,
	      mysqrt(simu->ole->p[isim*nmod])*1e9,
	      mysqrt(simu->ole->p[1+isim*nmod])*1e9,
	      mysqrt(simu->ole->p[2+isim*nmod])*1e9);
	
	info2("Timing: Tot:%5.2f Mean:%5.2f. Used %ld:%02ld, Left %ld:%02ld\n",
	      status->tot*tkmean, status->mean*tkmean, lapsh,lapsm,resth,restm);
    }else{    
	info2("\033[00;32mStep %5d: OL: %6.1f %6.1f %6.1f nm; CL %6.1f %6.1f %6.1f nm;",
	      isim,
	      mysqrt(simu->ole->p[isim*nmod])*1e9,
	      mysqrt(simu->ole->p[1+isim*nmod])*1e9,
	      mysqrt(simu->ole->p[2+isim*nmod])*1e9,
	      mysqrt(simu->cle->p[isim*nmod])*1e9,
	      mysqrt(simu->cle->p[1+isim*nmod])*1e9,
	      mysqrt(simu->cle->p[2+isim*nmod])*1e9);
	if(parms->recon.split && parms->ndm<=2){
	    info2(" Split %6.1f %6.1f %6.1f nm;",
		  mysqrt(simu->clem->p[isim*3])*1e9,
		  mysqrt(simu->clem->p[1+isim*3])*1e9,
		  mysqrt(simu->clem->p[2+isim*3])*1e9);
	}
	info2("\033[00;00m\n");
    
	info2("Timing: WFS:%5.2f Recon:%6.4f CACHE:%5.2f EVAL:%5.2f Tot:%5.2f Mean:%5.2f."
	      " Used %ld:%02ld, Left %ld:%02ld\n",
	      status->wfs*tkmean, status->recon*tkmean, 
	      status->cache*tkmean, status->eval*tkmean, 
	      status->tot*tkmean, status->mean*tkmean,
	      lapsh,lapsm,resth,restm);
	if(isnan(mysqrt(simu->cle->p[isim*nmod])*1e9)){
	    error("NaN found\n");
	}
    }
}
/**
   Output parameters necessary to run postproc using skyc/skyc.c
*/
void save_skyc(POWFS_T *powfs, RECON_T *recon, const PARMS_T *parms){
    char fn[PATH_MAX];
    double zadeg=parms->sim.za*180/M_PI;
    snprintf(fn,PATH_MAX,"%s/maos.conf",dirskysim);
    FILE *fp=fopen(fn,"w");
    fprintf(fp,"maos.r0z=%g\n",parms->atm.r0z);
    fprintf(fp,"maos.dt=%g\n",parms->sim.dt);
    fprintf(fp,"maos.zadeg=%g\n",zadeg);
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
		int iwfs=parms->powfs[ipowfs].wfs[indwfs];
		int iwfs2=parms->powfs[ipowfs].wfs[indwfs+1];
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
	fprintf(fp,"\"powfs%d_loc\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnwfsamp=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_amp\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnsaloc=[");
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	fprintf(fp,"\"powfs%d_saloc\"," ,powfs_ngs[ipowfs]);
    }
    fprintf(fp,"]\n");
    
    fprintf(fp,"maos.fnmideal=\"RescleNGSm\"\n");
    fprintf(fp,"maos.fnmidealp=\"RescleNGSmp\"\n");
    fprintf(fp,"maos.evlindoa=%d\n",parms->evl.indoa);
    fprintf(fp,"maos.fnmcc=\"MCC_za%g\"\n",zadeg);
    fprintf(fp,"maos.fnmcc_oa=\"MCC_OA_za%g\"\n",zadeg);
    
    fprintf(fp,"maos.seeds=[");
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	fprintf(fp,"%d " ,parms->sim.seeds[iseed]);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"include=\"skyc.conf\"\n");
    fprintf(fp,"include=\"skyc_za%g.conf\"\n",zadeg);
    fprintf(fp,"maos.wddeg=[");
    for(int ips=0; ips<parms->atm.nps; ips++){
	fprintf(fp, "%.2f ", parms->atm.wddeg[ips]);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"maos.nstep=%d\n",parms->sim.end);
    fclose(fp);
    for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
	int jpowfs=powfs_ngs[ipowfs];
	locwrite(powfs[jpowfs].loc,"%s/powfs%d_loc",dirskysim,jpowfs);
	dwrite(powfs[jpowfs].amp, "%s/powfs%d_amp",dirskysim,jpowfs);
	locwrite(powfs[jpowfs].saloc,"%s/powfs%d_saloc",dirskysim,jpowfs);
    }
    dwrite(recon->ngsmod->MCC,"%s/MCC_za%g", dirskysim,zadeg);
    dwrite(recon->ngsmod->MCCP->p[parms->evl.indoa],"%s/MCC_OA_za%g", dirskysim,zadeg);
}
