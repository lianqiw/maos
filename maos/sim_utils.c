/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <unistd.h>

#include "common.h"
#include "sim.h"
#include "sim_utils.h"
#include "surf.h"
#include "powfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

/**
   \file sim_utils.h

   Contains a few support functions for simulation.
*/
extern int disable_save;
real tk_setup;   /**<Start time of setup */
static mapcell* genatm_do(sim_t* simu){
	const parms_t* parms=simu->parms;
	const atm_cfg_t* atm=&parms->atm;
	TIC;
	mapcell* screens;
	if(!parms->dbg.atm){
		genatm_t* gs=simu->atmcfg;
		if(!gs){
			simu->atmcfg=mycalloc(1, genatm_t);/*the data for generating screens. */
			gs=simu->atmcfg;
			gs->rstat=simu->atm_rand;
			gs->wt=P(atm->wt);
			gs->r0=atm->r0;
			gs->L0=P(atm->L0);
			gs->dx=atm->dx;
			gs->nx=NX(atm);
			gs->ny=NY(atm);
			gs->nlayer=atm->nps;
			gs->ninit=atm->ninit;
			gs->share=atm->share;
			if((parms->atm.r0evolve&1)==1){
				dbg("Scaling turbulence screen spatially\n");
				gs->r0logpsds=atm->r0logpsds;
			}
			switch(parms->atm.method){
			case 0:
				gs->slope=-11./3.;
				break;
			case 1:
				gs->slope=-1.;
				break;
			case 2:
				gs->slope=-4.;
				break;
			default:
				gs->slope=parms->atm.method;
			}
		}
		tic;
		screens=genscreen(gs);
		toc2("Atmosphere ");
	} else{
	/*
	  create screens on two layers that produce pure
	  tip/tilt for LGS to debug split tomography test
	  pass in open loop mode for both Z and G tilt for
	  NGS.

	  The residual error is pure fit error if
	  layer 5 is at dm layer. otherwise the error is a
	  little larger.
	*/
		int nx=NX(atm);
		int ny=NY(atm);
		screens=mapcellnew(atm->nps, 1);
		real dx=atm->dx;
		loc_t* psloc=0;
		const real strength=sqrt(1.0299*pow(parms->aper.d/atm->r0, 5./3.))*(0.5e-6/(2*M_PI));//PR WFE.
		for(int ips=0; ips<atm->nps; ips++){
			P(screens, ips)=mapnew(nx, ny, dx, dx);
			P(screens, ips)->h=P(atm->ht, ips);
			if(!P(atm->wt, ips)) continue;
			real dbgatm=PR(parms->dbg.atm, ips);
			if(dbgatm>0){//zernike mode
				if(!psloc){
					psloc=mksqloc_auto(nx, ny, atm->dx, atm->dx);
				}
				info2("Generating Testing Atmosphere Screen with zernike %g, RMS~=%g nm\n", dbgatm, P(atm->wt, ips)*strength*1e9);
				dmat* opd=zernike(psloc, nx*atm->dx, 0, 0, -dbgatm);
				dmat* opd2=dref_reshape(opd, nx, ny);
				dadd((dmat**)&P(screens, ips), 0, opd2, P(atm->wt, ips)*strength);
				dfree(opd);
				dfree(opd2);
			} else if(dbgatm<0){//Fourier mode;
				info2("Generating Testing Atmosphere Screen with Fourier mode %g, RMS~=%g nm\n", -dbgatm, P(atm->wt, ips)*strength*1e9);

				real kk=2*M_PI/dbgatm*dx;
				long nn=MAX(nx, ny);
				dmat* sind=dnew(nn, 1);
				for(int ii=0; ii<nn; ii++){
					P(sind, ii)=sin((ii-nn/2)*kk);
				}
				const real alpha=2*strength*P(atm->wt, ips);
				for(int iy=0; iy<ny; iy++){
					for(int ix=0; ix<nx; ix++){
						P(P(screens, ips), ix, iy)=alpha*P(sind, ix)*P(sind, iy);
					}
				}
				dfree(sind);
			} else{
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

void blend_screen_side(map_t* atm1, map_t* atm2, long overx, long overy){
	const long nx=NX(atm1);
	const long ny=NY(atm1);
	int ca=0;
	if(atm1->vx>EPS){
		ca=-1;/*reverse sign of vx */
	} else if(atm1->vx<-EPS){
		ca=1;
	}
	int sa=0;
	if(atm1->vy>EPS){
		sa=-1;/*reverse sign of vy */
	} else if(atm1->vy<-EPS){
		sa=1;
	}
	long rr;
	long offx=nx-overx;
	long offy=(ny-overy)*nx;
	if(ca==0){/*along y. */
		rr=ny-overy;/*distance between the origins. */
		atm2->oy=atm1->oy+rr*sa*atm1->dx;
		atm2->ox=atm1->ox;
		real wty=sa<0?1:0;
		real* p1=P(atm1)+(1-(long)wty)*offy;
		real* p2=P(atm2)+(long)wty*offy;
		real(*pp1)[nx]=(real(*)[nx])p1;
		real(*pp2)[nx]=(real(*)[nx])p2;
		real overyd=(real)overy;
		for(long iy=0; iy<overy; iy++){
			real wt1=fabs(wty-(real)(iy+1)/overyd);
			for(long ix=0; ix<nx; ix++){
				pp1[iy][ix]=(1-wt1)*pp1[iy][ix]+wt1*pp2[iy][ix];
				pp2[iy][ix]=pp1[iy][ix];
			}
		}
	} else if(sa==0){
		rr=nx-overx;/*distance between the origins. */
		atm2->ox=atm1->ox+rr*ca*atm1->dx;
		atm2->oy=atm1->oy;
		real wtx=ca<0?1:0;
		real* p1=P(atm1)+(1-(long)wtx)*offx;
		real* p2=P(atm2)+(long)wtx*offx;
		real(*pp1)[nx]=(real(*)[nx])p1;
		real(*pp2)[nx]=(real(*)[nx])p2;
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
	} else{
		error("We do not support this wind direction: ca=%d, sa=%d\n", ca, sa);
	}
}
/**
   wrap of the generic vonkarman_genatm to generate turbulence screens. Wind
   velocities are set for each screen.  \callgraph */
void genatm(sim_t* simu){
	const parms_t* parms=simu->parms;
	const atm_cfg_t* atm=&(simu->parms->atm);
	if(simu->atm){
		cellfree(simu->atm);
		dfree(simu->winddir);
	}
	if(simu->parms->sim.noatm){
		warning("sim.noatm flag is on. will not generate atmoshere\n");
		return;
	}
	info("Wind dir:");/*initialize wind direction one time only for each seed in frozen flow mode. */
	simu->winddir=dnew(atm->nps, 1);
	for(int i=0; i<atm->nps; i++){
		real angle;
		if(atm->wdrand){
			angle=randu(simu->atmwd_rand)*M_PI*2;
		} else{
			angle=P(atm->wddeg, i)*M_PI/180;
		}
		P(simu->winddir, i)=angle;
		info(" %5.1f", angle*180/M_PI);
	}
	info(" deg\n");
	int atm_movie=0;//whether loaded atm is movie to be playback.
	if(simu->parms->load.atm){
		const char* fn=simu->parms->load.atm;
		info("loading atm from %s\n", fn);
		simu->atm=mapcellread("%s", fn);
		if(parms->aper.rot){
			dmaprot(simu->atm, parms->aper.rot);
		}
		if(NX(simu->atm)!=atm->nps){
			if(parms->atm.dtrat && NX(simu->atm)>1){
				real h=P(simu->atm, 0)->h;
				int sameh=1;
				for(int ips=1; ips<NX(simu->atm); ips++){
					if(h!=P(simu->atm ,ips)->h){
						sameh=0;
						break;
					}
				}
				P(((parms_t *)parms)->atm.ht, 0)=h;
				P(((parms_t *)parms)->atm.ht, 1)=h;
				if(sameh){
					atm_movie=1;
				}
			}
			if(!atm_movie){
				error("Loaded turbulence does not match atm specification.\n");
			}
		}
	} else{
		simu->atm=genatm_do(simu);
	}
	if(parms->atm.dtrat&&!atm_movie){
		error("atm.dtrat is set but atm is not movie.\n");
	}
	if(!parms->dbg.atm&&atm->nps==parms->atm.nps&&!atm_movie){
		for(int i=0; i<atm->nps; i++){
			real angle=P(simu->winddir, i);
			P(simu->atm, i)->h=P(parms->atm.ht, i);
			P(simu->atm, i)->vx=cos(angle)*P(parms->atm.ws, i);
			P(simu->atm, i)->vy=sin(angle)*P(parms->atm.ws, i);
		}
	}else{
		dbg("Turbulence frozen flow velocity are set to 0.\n");
	}
	if(simu->parms->save.atm){
		writebin(simu->atm, "atm_%d.bin", simu->seed);
	}

	if(parms->plot.atm&&simu->atm){
		for(int ips=0; ips<atm->nps; ips++){
			drawmap("Atm", P(simu->atm, ips), 0,
				"Atmosphere OPD", "x (m)", "y (m)", "layer %d", ips);
		}
	}

	print_mem("After genatm");
	if((parms->atm.r0evolve&2)==2){
		dbg("Scaling OPD temporarily\n");
		simu->atmscale=psd2ts(parms->atm.r0logpsdt, simu->atm_rand, parms->sim.dt, parms->sim.end);
		const real r02wt=(-5./3.);//layer weight is prop to r0^(-5/3)
		for(long i=0; i<NX(simu->atmscale); i++){
			P(simu->atmscale, i)=exp((P(simu->atmscale, i))*r02wt);//convert to cn2dh
		}
		writebin(simu->atmscale, "atmscale_%d", simu->seed);
	}
	if(simu->wfs_prop_atm){
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			for(int ips=0; ips<parms->atm.nps; ips++){
				propdata_t* data=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
				data->mapin=P(simu->atm, ips);
			}
		}
	}
	if(simu->evl_prop_atm){
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			for(int ips=0; ips<parms->atm.nps; ips++){
				propdata_t* data=&simu->evl_propdata_atm[ievl+parms->evl.nevl*ips];
				data->mapin=P(simu->atm, ips);
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
void atm2xloc(dcell** opdx, const sim_t* simu){
	const recon_t* recon=simu->recon;
	const parms_t* parms=simu->parms;
	if(parms->recon.glao){
		return;
	}
	/*in close loop mode, opdr is from last time step. */
	int isim=simu->reconisim;
	if(!*opdx){
		*opdx=dcellnew(recon->npsr, 1);
	}
	for(int ipsr=0; ipsr<recon->npsr; ipsr++){
		if(!P(*opdx, ipsr)){
			P(*opdx, ipsr)=dnew(P(recon->xloc, ipsr)->nloc, 1);
		} else{
			dzero(P(*opdx, ipsr));
		}
	}
	if(simu->atm){
		for(int ips=0; ips<parms->atm.nps; ips++){
			real disx=-P(simu->atm, ips)->vx*isim*parms->sim.dt;
			real disy=-P(simu->atm, ips)->vy*isim*parms->sim.dt;
			int ipsr=P(parms->atm.ipsr, ips);
			prop_grid(P(simu->atm, ips), P(recon->xloc, ipsr), P(P(*opdx, ipsr)),
				1, disx, disy, 1, 1, 0, 0);
		}
	}
}
/**
   Evolving the Sodium layer by updating the elongation transfer function.
*/
void sim_update_etf(sim_t* simu){
	int isim=simu->wfsisim;
	const parms_t* parms=simu->parms;
	powfs_t* powfs=simu->powfs;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		/* Update ETF if necessary. */
		if(!parms->powfs[ipowfs].llt) continue;
		//Needs ETF for imaging
		const int has_phy=parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout;
		//Need to initialize trombone position
		const int zoomset=parms->powfs[ipowfs].zoomset
			&&simu->wfsisim==parms->sim.start+parms->powfs[ipowfs].zoomset
			&&parms->powfs[ipowfs].phytype_sim!=PTYPE_MF;
		//Time for sodium profile update
		const int na_update=parms->powfs[ipowfs].llt->coldtrat>0&&isim%parms->powfs[ipowfs].llt->coldtrat==0;
		if(has_phy &&(na_update||zoomset)){
			int icol=0, icol2=0;
			if(parms->powfs[ipowfs].llt->coldtrat>0){
				int dtrat=parms->powfs[ipowfs].llt->coldtrat;
				int colsim=parms->powfs[ipowfs].llt->colsim;
				icol=colsim+isim/dtrat;
				icol2=icol+1;
			}
			static real deltah_last=0;
			real deltah=0;
			//factor converts our definition of focus mode (x^2*y^2)*alpha to LGS height error (delta_h*D^2)/(8*hs^2)
			//for larger distance, the convertion should use 1/(2*h1)-1/(2*h2)=alpha
			const real factor=-2*pow(parms->powfs[ipowfs].hs, 2);
			if(simu->zoomint){
				if(!parms->powfs[ipowfs].zoomshare){
					error("Please implement\n");
				}
				int iwfs0=P(parms->powfs[ipowfs].wfs, 0);
				real zoomerr1=0, zoomerr2=0;
				if(P(simu->zoomavg_count, iwfs0)){//gradient based
					average_powfs(simu->zoomavg, parms->powfs[ipowfs].wfs, 1);
					zoomerr1=P(simu->zoomavg, iwfs0)/P(simu->zoomavg_count, iwfs0);
					dzero(simu->zoomavg);
					lzero(simu->zoomavg_count);
				}
				if(P(simu->zoomdrift_count, iwfs0)){//i0 based
					average_powfs(simu->zoomdrift, parms->powfs[ipowfs].wfs, 1);
					zoomerr2=P(simu->zoomdrift, iwfs0)/P(simu->zoomdrift_count, iwfs0);
					dzero(simu->zoomdrift);
					lzero(simu->zoomdrift_count);
				}
				if(zoomset){
					P(simu->zoomint, iwfs0)=(zoomerr1+zoomerr2);
					//more accurate method:
					real hold=parms->powfs[ipowfs].hs;
					real hnew=1./(1./hold+2*P(simu->zoomint, iwfs0));
					deltah=hnew-hold;
				}else{
					P(simu->zoomint, iwfs0)+=parms->powfs[ipowfs].zoomgain*zoomerr1
						+parms->powfs[ipowfs].zoomgain_drift*zoomerr2;
					deltah=P(simu->zoomint, iwfs0)*factor;//convert focus to height
				}
				if(simu->zoompos&&simu->zoompos_icol<PN(simu->zoompos, iwfs0)){
					P(P(simu->zoompos, iwfs0), simu->zoompos_icol)=P(simu->zoomint, iwfs0);
					simu->zoompos_icol++;
				}
				dbg("Step %d: powfs %d: trombone error is %g %g zoompos is %g zoomset=%d.\n", isim, ipowfs, zoomerr1*factor, zoomerr2*factor, deltah, zoomset);
				if(fabs(deltah-deltah_last)<parms->powfs[ipowfs].llt->na_thres){//limit trombone travel minimum step size.
					deltah=deltah_last;
				}
			}

			info("Step %d: powfs %d: Updating ETF.\n", isim, ipowfs);
			TIC;tic;
			setup_shwfs_etf(powfs, parms, ipowfs, 1, icol, deltah, zoomset?0:100);
			if(icol2!=icol){
				setup_shwfs_etf(powfs, parms, ipowfs, 2, icol2, deltah, 0);
			}
			toc2("ETF");
#if USE_CUDA
			if(parms->gpu.wfs){
				gpu_wfsgrad_update_etf(parms, powfs, ipowfs);
			}
#endif
		}
	}
}
/**
 * Update flags
 * */
void update_wfsflags(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(!simu->wfsflags){
		simu->wfsflags=mycalloc(parms->npowfs, wfsflags_t);
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		wfsflags_t* wfsflags=simu->wfsflags+ipowfs;
		const int dtrat=parms->powfs[ipowfs].dtrat;
		const int wfsframe=simu->wfsisim+1-parms->powfs[ipowfs].step;
		wfsflags->gradout=(wfsframe>0&&wfsframe%dtrat==0)?(wfsframe/dtrat):0;
		wfsflags->do_phy=parms->powfs[ipowfs].usephy&&simu->wfsisim>=parms->powfs[ipowfs].phystep;
		wfsflags->do_pistat=parms->powfs[ipowfs].pistatout&&simu->wfsisim>=parms->powfs[ipowfs].pistatstart;
		if(parms->powfs[ipowfs].dither){
			const int pllrat=parms->powfs[ipowfs].dither_pllrat;
			const int pllcount=(simu->wfsisim-parms->powfs[ipowfs].dither_pllskip+1)/dtrat;
			wfsflags->pllout=((pllcount>0)&&(pllcount%pllrat)==0)?(pllcount/pllrat):0;
			const int ograt=parms->powfs[ipowfs].dither_ograt;//multiple of pllrat
			const int ogcount=(simu->wfsisim-parms->powfs[ipowfs].dither_ogskip+1)/dtrat;
			wfsflags->ogacc=(ogcount>0&&(ogcount%pllrat)==0)?(ogcount/pllrat):0;
			wfsflags->ogout=(ogcount>0&&(ogcount%ograt)==0)?(ogcount/ograt):0;
			if(simu->gradoffisim0<=0
				&&simu->wfsisim>=parms->powfs[ipowfs].dither_ogskip
				&&parms->dbg.gradoff_reset==2
				&&parms->powfs[ipowfs].phytype_sim2==PTYPE_MF){
				simu->gradoffisim0=simu->wfsisim;
				simu->gradoffisim=simu->wfsisim;
			}
		}
	}
}

/**
   Shift gradient when new gradients are ready (in the end of parallel section
   in sim in CL or wfsgrad in OL). Do not execute in parallel with other
   routines. In GLAO mode, also averaged gradients from the same type of powfs.
*/
void shift_grad(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(parms->nwfs==0||parms->sim.evlol||parms->sim.idealtomo) return;
	if(PARALLEL==2){
		pthread_mutex_lock(&simu->wfsgrad_mutex);
		if(simu->wfsisim>0){
			while(simu->wfsgrad_count<1){//not being consumed yet
				//dbg("waiting: wfsgrad_count is %d, need %d\n", simu->wfsgrad_count, 1);
				struct timespec ts;
				clock_gettime(CLOCK_REALTIME, &ts);
				ts.tv_nsec+=1e6;
				pthread_cond_timedwait(&simu->wfsgrad_condw, &simu->wfsgrad_mutex, &ts);
			}
			//dbg("ready: wfsgrad_count is ready: %d\n", simu->wfsgrad_count);
		}
		pthread_mutex_unlock(&simu->wfsgrad_mutex);
	}
	if(parms->recon.glao){
		/* Average the gradients in GLAO mode. */
		if(simu->gradlastcl){
			dcellzero(simu->gradlastcl);
		} else{
			long nnx[parms->nwfsr];
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				int ipowfs=parms->wfsr[iwfs].powfs;
				nnx[iwfs]=simu->powfs[ipowfs].saloc->nloc*2;
			}
			simu->gradlastcl=dcellnew3(parms->nwfsr, 1, nnx, NULL);
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			const real scale=1./parms->powfs[ipowfs].nwfs;
			for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, indwfs);
				dadd(&P(simu->gradlastcl, ipowfs), 1., P(simu->gradcl, iwfs), scale);
			}
		}
	} else{
		dcellcp(&simu->gradlastcl, simu->gradcl);
	}
	if(PARALLEL==2){
		pthread_mutex_lock(&simu->wfsgrad_mutex);
		//Signal recon wfsgrad is ready/
		simu->wfsgrad_isim=simu->wfsisim;
		simu->wfsgrad_count=0;//reset the counter
		pthread_cond_broadcast(&simu->wfsgrad_condr);
		//dbg("wfsgrad_isim is set to %d\n", simu->wfsgrad_isim);
		pthread_mutex_unlock(&simu->wfsgrad_mutex);
	}
}

/**
   use random number dirived from input seed to seed other stream.  necessary to
   have independant streams for different wfs in threading routines to avoid
   race condition and have consitent result */
void seeding(sim_t* simu){
	info2("Running seed %d\n", simu->seed);
	simu->init_rand=mycalloc(1, rand_t);
	simu->atm_rand=mycalloc(1, rand_t);
	simu->atmwd_rand=mycalloc(1, rand_t);
	simu->telws_rand=mycalloc(1, rand_t);
	simu->misc_rand=mycalloc(1, rand_t);
	seed_rand(simu->init_rand, simu->seed);
	seed_rand(simu->atm_rand, lrand(simu->init_rand));
	/*2011-02-02: changed to wdrand-1 so that when wdrand=1, we reproduce old directions. */
	seed_rand(simu->atmwd_rand, lrand(simu->init_rand)+(simu->parms->atm.wdrand-1));
	const parms_t* parms=simu->parms;
	simu->wfs_rand=mycalloc(parms->nwfs, rand_t);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		seed_rand(&simu->wfs_rand[iwfs], lrand(simu->init_rand));
	}
	seed_rand(simu->telws_rand, lrand(simu->init_rand)*parms->sim.wsseq);
#if USE_CUDA
	if(parms->gpu.wfs&&!parms->sim.evlol){
		gpu_wfsgrad_seeding(parms, simu->powfs, simu->init_rand);
	}
#endif
	seed_rand(simu->misc_rand, simu->seed);
}
const char *fnextra="-";//replaced when save.extra is set
static void init_simu_evl(sim_t* simu){
	const parms_t* parms=simu->parms;
	const aper_t* aper=simu->aper;
	recon_t* recon=simu->recon;
	const int nsim=parms->sim.end;
	const int nevl=parms->evl.nevl;
	const int nmod=parms->evl.nmod;
	const int seed=simu->seed;
	sim_save_t* save=simu->save;
	simu->evlopd=dcellnew(nevl, 1);
	simu->perfevl_iground=parms->atm.iground;
	/*if(!disable_save&&parms->save.extra>1){
		simu->timing=dnew_file(5, nsim, "Simulation timing per time step", "%s/Timing_%d.bin", fnextra, seed);
	}*/
	//2021-07-2: do not use mmap to write file. It causes a lot of disk activity.
	/*MMAP the main result file */
	{
		int ahst_nmod=0;
		if(parms->evl.split){
			ahst_nmod=3;
			if(parms->nlgspowfs) ahst_nmod++;//focus
		}
		long nnx[4]={nmod,0,nmod,ahst_nmod};
		long nny[4]={nsim,0,nsim,nsim};

		simu->res=dcellnew_file(4, 1, nnx, nny, "Field averaged open and closed loop WFE.", "Res_%d.bin", seed);
		//Do not reference. Just assign. Don't free.
		simu->ole=P(simu->res, 0);
		//P(simu->res, 1) is no longer used.
		simu->cle=P(simu->res, 2);
		simu->clem=P(simu->res, 3);
		dcellset(simu->res, NAN);
	}

	{/*USE async write for data that need to save at every time step */
		char keywords[1024];
		{
			char *pk=keywords;
			char *pkend=keywords+sizeof(keywords);
			for(int j=0; j<2; j++){
				char elem[16];
				dmat *theta=j==0?parms->evl.thetax:parms->evl.thetay;
				for(int i=-1; i<nevl+1; i++){
					if(i==-1){
						snprintf(elem, sizeof(elem), "evl.theta%s=[", j==0?"x":"y");
					}else if(i==nevl){
						snprintf(elem, sizeof(elem), "]; ");
					}else{
						snprintf(elem, sizeof(elem), " %.2f", P(theta, i)*RAD2AS);
					}
					int ne=strnlen(elem, sizeof(elem));
					if(pk+ne<pkend){
						memcpy(pk, elem, ne+1);
						pk+=ne;
					}
				}
			}
		}
		simu->resp=dcellnewsame_file(nevl, 4, nmod, nsim, keywords, "Resp_%d.bin", seed);

		simu->olmp=dcellsub(simu->resp, 0, nevl, 0, 1);
		simu->clmp=dcellsub(simu->resp, 0, nevl, 1, 1);
		simu->olep=dcellsub(simu->resp, 0, nevl, 2, 1);
		simu->clep=dcellsub(simu->resp, 0, nevl, 3, 1);

		if(parms->evl.split){
			simu->oleNGSm=dnew_file(recon->ngsmod->nmod, nsim, keywords, "%s/ResoleNGSm_%d.bin", fnextra, seed);
			simu->oleNGSmp=dcellnewsame_file(nevl, 1, recon->ngsmod->nmod, nsim, keywords, "%s/ResoleNGSmp_%d.bin", fnextra, seed);
			if(!parms->sim.evlol){
				simu->clemp=dcellnewsame_file(nevl, 1, 3, nsim, keywords, "%s/Resclemp_%d.bin", fnextra, seed);
				simu->cleNGSm=dnew_file(recon->ngsmod->nmod, nsim, keywords, "%s/RescleNGSm_%d.bin", fnextra, seed);
				simu->cleNGSmp=dcellnewsame_file(nevl, 1, recon->ngsmod->nmod, nsim, keywords, "%s/RescleNGSmp_%d.bin", fnextra, seed);
			}
		}
		if(parms->recon.split==1&&!parms->sim.fuseint&&!parms->sim.evlol){
			simu->corrNGSm=dnew_file(recon->ngsmod->nmod, nsim, keywords, "%s/RescorrNGSm_%d.bin", fnextra, seed);
		}
	}
	if(parms->sim.skysim){
		char fnold[PATH_MAX];
		char fnnew[PATH_MAX];
		snprintf(fnnew, PATH_MAX, "%s/RescleNGSm_%d.bin", dirskysim, seed);
		snprintf(fnold, PATH_MAX, "%s/RescleNGSm_%d.bin", fnextra, seed);
		if(exist(fnnew)) remove(fnnew);
		if(link(fnold, fnnew)){
			warning("Error link\n");
		}
		snprintf(fnnew, PATH_MAX, "%s/RescleNGSmp_%d.bin", dirskysim, seed);
		snprintf(fnold, PATH_MAX, "%s/RescleNGSmp_%d.bin", fnextra, seed);
		if(exist(fnnew)) remove(fnnew);
		if(link(fnold, fnnew)){
			warning("Error link\n");
		}
	}

	if(parms->evl.psfmean||parms->evl.psfhist||parms->evl.cov||parms->evl.opdmean){
		//The saved PSF and COVs are padded by empty cells.
		long nframepsf=parms->sim.end;
		char strht[24];
		if(parms->evl.psfmean&&!parms->sim.evlol){
			simu->evlpsfmean=dcellnew(parms->evl.nwvl, nevl);
			save->evlpsfmean=mycalloc(nevl, zfarr*);
			simu->evlpsfmean_ngsr=dcellnew(parms->evl.nwvl, nevl);
			save->evlpsfmean_ngsr=mycalloc(nevl, zfarr*);
		}
		if(parms->evl.psfhist){
			save->evlpsfhist=mycalloc(nevl, zfarr*);
			save->evlpsfhist_ngsr=mycalloc(nevl, zfarr*);
		}
		if(parms->evl.cov||parms->evl.opdmean){//need cell array of both
			simu->evlopdcov=dcellnew(nevl, 1);
			simu->evlopdcov_ngsr=dcellnew(nevl, 1);
			save->evlopdcov=mycalloc(nevl, zfarr*);
			save->evlopdcov_ngsr=mycalloc(nevl, zfarr*);
			simu->evlopdmean=dcellnew(nevl, 1);
			simu->evlopdmean_ngsr=dcellnew(nevl, 1);
			save->evlopdmean=mycalloc(nevl, zfarr*);
			save->evlopdmean_ngsr=mycalloc(nevl, zfarr*);
		}
		for(int ievl=0; ievl<nevl; ievl++){
			if(!P(parms->evl.psf, ievl)||parms->sim.evlol) continue;
			if(!isinf(P(parms->evl.hs, ievl))){
				snprintf(strht, 24, "_%g", P(parms->evl.hs, ievl));
			} else{
				strht[0]='\0';
			}
#define DIR_SUFFIX "_%d_x%g_y%g%s.fits", seed,P(parms->evl.thetax, ievl)*RAD2AS, P(parms->evl.thetay, ievl)*RAD2AS, strht
			if((P(parms->evl.psf, ievl)&1)){
				if(parms->evl.psfmean){
					save->evlpsfmean[ievl]=zfarr_init(parms->evl.nwvl, nframepsf, "evlpsfcl" DIR_SUFFIX);
				}
				if(parms->evl.psfhist){
					save->evlpsfhist[ievl]=zfarr_init(parms->sim.end, 1, "evlpsfhist" DIR_SUFFIX);
				}
				if(parms->evl.cov){
					save->evlopdcov[ievl]=zfarr_init(0, 0, "evlopdcov" DIR_SUFFIX);
				}
				if(parms->evl.cov||parms->evl.opdmean){
					save->evlopdmean[ievl]=zfarr_init(0, 0, "evlopdmean" DIR_SUFFIX);
				}
			}
			if((P(parms->evl.psf, ievl)&2)){
				if(parms->evl.psfmean){
					save->evlpsfmean_ngsr[ievl]=zfarr_init(parms->evl.nwvl, nframepsf, "evlpsfcl_ngsr" DIR_SUFFIX);
				}
				if(parms->evl.psfhist){
					save->evlpsfhist_ngsr[ievl]=zfarr_init(parms->sim.end, 1, "evlpsfhist_ngsr" DIR_SUFFIX);
				}
				if(parms->evl.cov){
					save->evlopdcov_ngsr[ievl]=zfarr_init(0, 0, "evlopdcov_ngsr" DIR_SUFFIX);
				}
				if(parms->evl.cov||parms->evl.opdmean){
					save->evlopdmean_ngsr[ievl]=zfarr_init(0, 0, "evlopdmean_ngsr" DIR_SUFFIX);
				}
			}
#undef DIR_SUFFIX
		}/*for ievl */
#define DIR_SUFFIX_OL "_%d.fits", seed
		if(parms->evl.psfol){
			if(parms->evl.psfmean){
				simu->evlpsfolmean=dcellnew(parms->evl.nwvl, 1);
				save->evlpsfolmean=zfarr_init(parms->evl.nwvl, nframepsf, "evlpsfol" DIR_SUFFIX_OL);
			}
			if(parms->evl.cov){
				save->evlopdcovol=zfarr_init(0, 0, "evlopdcovol" DIR_SUFFIX_OL);
			}
			if(parms->evl.cov||parms->evl.opdmean){
				save->evlopdmeanol=zfarr_init(0, 0, "evlopdmeanol" DIR_SUFFIX_OL);
			}
		}
#undef DIR_SUFFIX_OL
	}

	if(parms->save.ecov){
		if(!parms->dbg.useopdr||parms->sim.idealtomo){
			simu->ecov=dcellnew(parms->ndm, parms->ndm);/*not really need. */
		} else{/*deprecated */
			simu->ecov=dcellnew(nevl, 1);
			if(parms->dbg.ecovxx){/*temporary. */
				warning("Saving history of ecov_xx\n");
				save->ecovxx=mycalloc(nevl, zfarr*);
				char strht[24];
				for(int ievl=0; ievl<nevl; ievl++){
					if(!P(parms->evl.psf, ievl)) continue;
					if(!isinf(P(parms->evl.hs, ievl))){
						snprintf(strht, 24, "_%g", P(parms->evl.hs, ievl));
					} else{
						strht[0]='\0';
					}

					save->ecovxx[ievl]=zfarr_init(parms->sim.end, 1,
						"ecovxx_%d_x%g_y%g%s.bin", seed,
						P(parms->evl.thetax, ievl)*RAD2AS,
						P(parms->evl.thetay, ievl)*RAD2AS, strht);
				}
			}
		}
	}

	if(parms->save.evlopd){
		int nstep=parms->sim.end;
		save->evlopdol=mycalloc(nevl, zfarr*);
		save->evlopdcl=mycalloc(nevl, zfarr*);

		for(int ievl=0; ievl<nevl; ievl++){
			save->evlopdol[ievl]=zfarr_init(nstep, 1, "evl%d_opdol_%d.fits", ievl, seed);
			save->evlopdcl[ievl]=zfarr_init(nstep, 1, "evl%d_opdcl_%d.fits", ievl, seed);
		}
	}

	/* For threading */
	simu->evl_prop_atm=mycalloc(nevl*parms->atm.nps, thread_t*);
	simu->evl_propdata_atm=mycalloc(nevl*parms->atm.nps, propdata_t);
	simu->evl_prop_dm=mycalloc(nevl*parms->ndm, thread_t*);
	simu->evl_propdata_dm=mycalloc(nevl*parms->ndm, propdata_t);
	for(int ievl=0; ievl<nevl; ievl++){
		int nthread=2;
		int tot;
		for(int ips=0; ips<parms->atm.nps; ips++){
			const int ind=ievl+nevl*ips;
			propdata_t* data=&simu->evl_propdata_atm[ind];
			const real ht=P(parms->atm.ht, ips);
			data->displacex0=ht*P(parms->evl.thetax, ievl);
			data->displacey0=ht*P(parms->evl.thetay, ievl);
			data->scale=1-ht/P(parms->evl.hs, ievl);
			data->alpha=1;
			data->wrap=1;
			data->ostat=aper->locs->stat;
			tot=aper->locs->stat->ncol;
			simu->evl_prop_atm[ind]=thread_prep(0, tot, nthread, prop, data);
		}
		for(int idm=0; idm<parms->ndm&&!parms->sim.evlol; idm++){
			const int ind=ievl+nevl*idm;
			propdata_t* data=&simu->evl_propdata_dm[ind];
			const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
			data->displacex0=ht*P(parms->evl.thetax, ievl);
			data->displacey0=ht*P(parms->evl.thetay, ievl);
			data->scale=1-ht/P(parms->evl.hs, ievl);
			const real theta=RSS(P(parms->evl.thetax, ievl), P(parms->evl.thetay, ievl));
			data->alpha=-cos(theta*parms->dm[idm].dratio);
			data->wrap=0;
			if(parms->sim.cachedm){
				data->mapin=P(simu->cachedm, idm);
			} else{
				if(simu->dmrealsq){
					data->mapin=P(simu->dmrealsq, idm);
				} else{
					data->locin=P(recon->aloc, idm);
					data->phiin=P(simu->dmreal, idm);
				}
			}
			if(aper->locs_dm){
				data->locout=P(aper->locs_dm, ind);
				tot=data->locout->nloc;
			} else{
				data->ostat=aper->locs->stat;
				tot=aper->locs->stat->ncol;
			}
			data->phiout=(dmat*)1;/*replace later in simulation. */
			simu->evl_prop_dm[ind]=thread_prep(0, tot, nthread, prop, data);
		}
	}
}

static void init_simu_wfs(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(parms->sim.idealtomo) return;
	powfs_t* powfs=simu->powfs;
	recon_t* recon=simu->recon;
	sim_save_t* save=simu->save;
	const int nwfs=parms->nwfs;
	const int nsim=parms->sim.end;
	const int seed=simu->seed;
	simu->eptwfs=parms->sim.eptwfs;
	simu->ints=dccellnew(nwfs, 1);
	simu->intsout=dccellnew(nwfs, 1);
	simu->wfspsfout=cccellnew(nwfs, 1);
	simu->pistatout=dccellnew(nwfs, 1);
	save->wfspsfout=mycalloc(nwfs, zfarr*);
	save->ztiltout=mycalloc(nwfs, zfarr*);
	simu->gradcl=dcellnew(nwfs, 1);
	simu->wfsopd=dcellnew(nwfs, 1);

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
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		parms->powfs[ipowfs].phytype_sim=parms->powfs[ipowfs].phytype_sim1;//restore the value
	}
	simu->gradacc=dcellnew(nwfs, 1);/*wfsgrad internal */
	simu->gradoff=dcellnew(nwfs, 1);
	simu->gradscale=dcellnew(nwfs, 1);//leave empty first
	simu->gradscale2=dcellnew(nwfs, 1);//leave empty first.
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		const int nsa=powfs[ipowfs].saloc->nloc;
		const int ng=parms->powfs[ipowfs].ng;
		P(simu->gradcl, iwfs)=dnew(nsa*ng, 1);
		P(simu->wfsopd, iwfs)=dnew(powfs[ipowfs].loc->nloc, 1);
		if(parms->powfs[ipowfs].usephy){
			if(parms->powfs[ipowfs].type==WFS_SH){
				P(simu->ints, iwfs)=dcellnew_same(nsa, 1, powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
			} else{
				P(simu->ints, iwfs)=dcellnew_same(1, 1, powfs[ipowfs].saloc->nloc, powfs[ipowfs].pywfs->cfg->nside);
			}
		}
		if(parms->powfs[ipowfs].phystep!=0||P(parms->save.gradgeom, iwfs)||parms->powfs[ipowfs].pistatout){
			P(simu->gradacc, iwfs)=dnew(nsa*ng, 1);
		}
		if(parms->powfs[ipowfs].pistatout){
			P(simu->pistatout, iwfs)=dcellnew(nsa, parms->powfs[ipowfs].nwvl);
		}
		if(powfs[ipowfs].gradncpa&&!(parms->powfs[ipowfs].phytype_sim==1&&parms->powfs[ipowfs].ncpa_method==NCPA_I0)){
			//CMF has gradncpa with in matched filter
			int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
			dbg3("wfs %d: copying gradncpa to gradoff\n", iwfs);
			dadd(&P(simu->gradoff, iwfs), 1, PR(powfs[ipowfs].gradncpa, wfsind, 1), 1);
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
		//TODO: split implementation for each POWFS.
		long nnx[nwfs];
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			nnx[iwfs]=0;
			if(parms->powfs[ipowfs].llt||parms->powfs[ipowfs].dither==1){
				nnx[iwfs]=2;
			}
		}
		simu->fsmerr_store=dcellnew3(nwfs, 1, nnx, NULL);
		simu->fsmerr_drift=dcellnew3(nwfs, 1, nnx, NULL);
		simu->fsmcmd=dcellnew3(nwfs, 1, nnx, NULL);
		simu->fsmreal=dcellnew3(nwfs, 1, nnx, NULL);
		simu->fsmint=mycalloc(parms->nwfs, servo_t *);
		simu->fsmsho=mycalloc(parms->nwfs, sho_t *);
		if(parms->sim.closeloop){
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				if(parms->powfs[ipowfs].llt||parms->powfs[ipowfs].dither==1){
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
						simu->fsmint[iwfs]=servo_new_scalar(P(simu->fsmreal, iwfs), parms->powfs[ipowfs].apfsm, parms->powfs[ipowfs].alfsm, parms->sim.dt, parms->powfs[ipowfs].epfsm);
						simu->fsmsho[iwfs]=sho_new(parms->powfs[ipowfs].f0fsm, parms->powfs[ipowfs].zetafsm);
					}
				}
			}
		}
	}

	{/*MMAP the LGS fsmlink error/command output */
		long nnx[nwfs];
		long nny[nwfs];
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			nnx[iwfs]=0;
			if(parms->powfs[ipowfs].llt||parms->powfs[ipowfs].dither==1){
				nnx[iwfs]+=2;//tip/tilt of FSM
			}
			if(parms->powfs[ipowfs].dither>1){
				nnx[iwfs]+=1;//signal in common path dithering
			}
			if(nnx[iwfs]){
				nny[iwfs]=nsim;
			} else{
				nny[iwfs]=0;
			}
		}

		simu->save->fsmerrs=dcellnew_file(nwfs, 1, nnx, nny, "Uplink FSM error time history", "%s/Resfsmerr_%d.bin", fnextra, seed);
		simu->save->fsmcmds=dcellnew_file(nwfs, 1, nnx, nny, "Uplink FSM demand time history", "%s/Resfsmcmd_%d.bin", fnextra, seed);
	}

	/* For sky coverage telemetry output */
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		long notfx=powfs[ipowfs].notfx;
		long notfy=powfs[ipowfs].notfy;
		long notf=MAX(notfx, notfy);
		if(parms->powfs[ipowfs].psfout){
			const int nsa=powfs[ipowfs].saloc->nloc;
			/*The PSFs here are PSFs of each subaperture. */
			P(simu->wfspsfout, iwfs)=ccellnew_same(nsa, parms->powfs[ipowfs].nwvl, notf/2+2, notf/2+2);
			mymkdir("%s/wvfout/", dirskysim);
			mymkdir("%s/ztiltout/", dirskysim);
			save->wfspsfout[iwfs]=zfarr_init(parms->sim.end-parms->sim.start, 1,
				"%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g.bin",
				dirskysim, seed,
				parms->powfs[ipowfs].order,
				parms->wfs[iwfs].thetax*RAD2AS,
				parms->wfs[iwfs].thetay*RAD2AS);

			save->ztiltout[iwfs]=zfarr_init(parms->sim.end-parms->sim.start, 1,
				"%s/ztiltout/ztiltout_seed%d_sa%d_x%g_y%g.bin",
				dirskysim, seed,
				parms->powfs[ipowfs].order,
				parms->wfs[iwfs].thetax*RAD2AS,
				parms->wfs[iwfs].thetay*RAD2AS);

		}
		if(parms->sim.skysim&&parms->powfs[ipowfs].pistatout){
			mymkdir("%s/pistat/", dirskysim);
		}
	}

	if(parms->save.ngcov>0){
		simu->gcov=dcellnew(parms->save.ngcov, 1);
	}
	int nstep=parms->sim.end;
	if(parms->save.wfsopd){
		save->wfsopd=mycalloc(nwfs, zfarr*);
		save->wfsopdol=mycalloc(nwfs, zfarr*);
		save->wfslltopd=mycalloc(nwfs, zfarr*);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			if(!P(parms->save.wfsopd, iwfs)){
				continue;
			}
			save->wfsopd[iwfs]=zfarr_init(nstep, 1, "wfs%d_opd_%d.bin", iwfs, seed);
			save->wfsopdol[iwfs]=zfarr_init(nstep, 1, "wfs%d_opdol_%d.bin", iwfs, seed);
			if(powfs[ipowfs].llt){
				save->wfslltopd[iwfs]=zfarr_init(nstep, 1, "wfs%d_lltopd_%d.bin", iwfs, seed);
			}
		}
	}
	if(parms->save.ints){
		save->intsny=mycalloc(nwfs, zfarr*);
		save->intsnf=mycalloc(nwfs, zfarr*);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			if(P(parms->save.ints, iwfs)&&parms->powfs[ipowfs].usephy){
				save->intsnf[iwfs]=zfarr_init(nstep, 1, "wfs%d_intsnf_%d.bin", iwfs, seed);
				if(parms->powfs[ipowfs].noisy){
					save->intsny[iwfs]=zfarr_init(nstep, 1, "wfs%d_intsny_%d.bin", iwfs, seed);
				}
			}
		}
	}

	save->gradcl=mycalloc(nwfs, zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(P(parms->save.grad, iwfs)){
			save->gradcl[iwfs]=zfarr_init(nstep, 1, "wfs%d_gradcl_%d.bin", iwfs, seed);
		}
	}

	save->gradol=mycalloc(nwfs, zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].psol&&P(parms->save.gradpsol, iwfs)){
			save->gradol[iwfs]=zfarr_init(nstep, 1, "wfs%d_gradpsol_%d.bin", iwfs, seed);
		}
	}

	save->gradnf=mycalloc(nwfs, zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(P(parms->save.gradnf, iwfs)&&parms->powfs[ipowfs].noisy){
			save->gradnf[iwfs]=zfarr_init(nstep, 1, "wfs%d_gradnf_%d.bin", iwfs, seed);
		}
	}

	save->gradgeom=mycalloc(nwfs, zfarr*);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(P(parms->save.gradgeom, iwfs)&&parms->powfs[ipowfs].phystep>-1){
			save->gradgeom[iwfs]=zfarr_init(nstep, 1, "wfs%d_gradgeom_%d.bin", iwfs, seed);
		}
	}
	/* threading */
	simu->wfs_prop_atm=mycalloc(nwfs*parms->atm.nps, thread_t*);
	simu->wfs_propdata_atm=mycalloc(nwfs*parms->atm.nps, propdata_t);
	simu->wfs_prop_dm=mycalloc(nwfs*parms->ndm, thread_t*);
	simu->wfs_propdata_dm=mycalloc(nwfs*parms->ndm, propdata_t);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
		const real hs=parms->wfs[iwfs].hs;

		for(int ips=0; ips<parms->atm.nps; ips++){
			const real ht=P(parms->atm.ht, ips);
			propdata_t* data=&simu->wfs_propdata_atm[iwfs+nwfs*ips];
			data->scale=1.-ht/hs;
			data->displacex0=ht*parms->wfs[iwfs].thetax+(parms->powfs[ipowfs].type==WFS_SH?data->scale*parms->wfs[iwfs].misregx:0);
			data->displacey0=ht*parms->wfs[iwfs].thetay+(parms->powfs[ipowfs].type==WFS_SH?data->scale*parms->wfs[iwfs].misregy:0);
			data->alpha=1;
			data->wrap=1;
			data->mapin=(map_t*)1;/*need to update this in genatm. */
			data->phiout=(dmat*)1;/*replace later in simulation. */
			int tot=0;
			if(powfs[ipowfs].loc_tel){/*misregistration. */
				data->locout=P(powfs[ipowfs].loc_tel, wfsind);
				tot=data->locout->nloc;
			} else if(parms->powfs[ipowfs].type==WFS_PY||powfs[ipowfs].saloc->nloc<NTHREAD){
				data->locout=powfs[ipowfs].loc;
				tot=data->locout->nloc;
			} else{
				data->ptsout=powfs[ipowfs].pts;
				tot=data->ptsout->nsa;
			}
			int nthread=data->scale==1?2:8;//use more threads for LGS as ray tracing is slow due to cone coordinate.
			simu->wfs_prop_atm[iwfs+nwfs*ips]=thread_prep(0, tot, nthread, prop, data);
		}
		for(int idm=0; idm<parms->ndm; idm++){
			const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
			propdata_t* data=&simu->wfs_propdata_dm[iwfs+nwfs*idm];
			data->scale=1.-ht/hs;
			data->displacex0=ht*parms->wfs[iwfs].thetax+(parms->powfs[ipowfs].type==WFS_SH?data->scale*parms->wfs[iwfs].misregx:0);
			data->displacey0=ht*parms->wfs[iwfs].thetay+(parms->powfs[ipowfs].type==WFS_SH?data->scale*parms->wfs[iwfs].misregy:0);
			const real theta=RSS(parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay);
			data->alpha=-cos(theta*parms->dm[idm].dratio);/*remove dm contribution. */
			data->wrap=0;
			if(parms->sim.cachedm){
				data->mapin=P(simu->cachedm, idm);
			} else{
				if(simu->dmrealsq){
					data->mapin=P(simu->dmrealsq, idm);
				} else{
					data->locin=P(recon->aloc, idm);
					data->phiin=P(simu->dmreal, idm);
				}
			}
			data->phiout=(dmat*)1;/*replace later in simulation */
			int tot;
			if(powfs[ipowfs].loc_dm || powfs[ipowfs].loc_tel){/*distortion or rotation. */
				data->locout=powfs[ipowfs].loc_dm?P(powfs[ipowfs].loc_dm, wfsind, idm):(powfs[ipowfs].loc_tel?P(powfs[ipowfs].loc_tel, wfsind):powfs[ipowfs].loc);
				tot=data->locout->nloc;
			} else if(parms->powfs[ipowfs].type==WFS_PY||powfs[ipowfs].saloc->nloc<NTHREAD){
				data->locout=powfs[ipowfs].loc;
				tot=data->locout->nloc;
			} else{
				data->ptsout=powfs[ipowfs].pts;
				tot=data->ptsout->nsa;
			}
			int nthread=data->scale==1?2:4;//use more threads for LGS as ray tracing is slow due to cone coordinate.
			simu->wfs_prop_dm[iwfs+nwfs*idm]=thread_prep(0, tot, nthread, prop, data);
		}/*idm */
	}/*iwfs */
	simu->wfs_intsdata=mycalloc(nwfs, wfsints_t);
	simu->wfs_ints=mycalloc(nwfs, thread_t*);

	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int tot=powfs[ipowfs].saloc->nloc;
		wfsints_t* data=simu->wfs_intsdata+iwfs;
		data->iwfs=iwfs;
		int nthread=NTHREAD;
		simu->wfs_ints[iwfs]=thread_prep(0, tot, nthread, wfsints, data);
	}
	if(parms->nlgspowfs){
		simu->llt_ws=dcellnew(parms->npowfs, 1);//windshake t/t per LLT
		simu->ltpm_sho=mycalloc(parms->npowfs, sho_t*); //LLT common FSM SHO filter
		simu->ltpm_real=dcellnew(parms->npowfs, 1);//LLT common FSM state
		simu->ltpm_cmd=dcellnew(parms->npowfs, 1);//LLT common FSM command
		simu->ltpm_lpf=dcellnew(parms->npowfs, 1);//LLT common FSM LPF
		long nnx[parms->npowfs];
		long nny[parms->npowfs];
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].llt){
				if(parms->powfs[ipowfs].llt->ttpsd){
					dmat *psdin=dread("%s", parms->powfs[ipowfs].llt->ttpsd);
					P(simu->llt_ws, ipowfs)=dnew(parms->sim.end, parms->powfs[ipowfs].llt->nllt);
					for(int i=0; i<parms->powfs[ipowfs].llt->nllt; i++){//per LLT.
						dmat *tmp=psd2ts(psdin, simu->misc_rand, parms->sim.dt, parms->sim.end);
						daddcol(P(simu->llt_ws, ipowfs), i, 0, tmp, 1);
						dfree(tmp);
					}
					dfree(psdin);
				}
				//P(simu->ltpm_cmd, ipowfs)=dnew(2, parms->powfs[ipowfs].llt->nllt);
				//P(simu->ltpm_real, ipowfs)=dnew(2, parms->powfs[ipowfs].llt->nllt);
				if(parms->powfs[ipowfs].llt->fcfsm){
					simu->ltpm_sho[ipowfs]=sho_new(parms->powfs[ipowfs].llt->fcfsm, 1);
				}
				nnx[ipowfs]=2;
				nny[ipowfs]=parms->sim.end;
			}else{
				nnx[ipowfs]=0;
				nny[ipowfs]=0;
			}
		}
		save->ltpm_real=dcellnew_file(parms->npowfs, 1, nnx, nny, "LLT center launch common path pointing mirror commands", "%s/Resltpm_%d", fnextra, seed);
	}
	if(parms->nlgspowfs){
		simu->LGSfocus=dcellnew(parms->nwfs, 1);
		simu->LGSfocus_drift=dcellnew(parms->nwfs, 1);
		simu->zoomdrift=dnew(parms->nwfs, 1);
		simu->zoomdrift_count=lnew(parms->nwfs, 1);
		simu->zoomint=dnew(parms->nwfs, 1);
		simu->zoomavg=dnew(parms->nwfs, 1);
		simu->zoomavg_count=lnew(parms->nwfs, 1);
		//To use writebin_async, the number of columns must be related to timestep
		long nnx[parms->nwfs];
		long nny[parms->nwfs];
		long nny2[parms->nwfs];
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].llt){
				nnx[iwfs]=1;
				nny[iwfs]=0;
				if(parms->powfs[ipowfs].llt->coldtrat){
					nny[iwfs]=parms->sim.end/parms->powfs[ipowfs].llt->coldtrat+1;
				}
				nny2[iwfs]=nsim;
			} else{
				nnx[iwfs]=0;
				nny[iwfs]=0;
				nny2[iwfs]=0;
			}
		}
		simu->zoompos=dcellnew_file(parms->nwfs, 1, nnx, nny, "LGS Trombone position", "%s/Reszoompos_%d.bin", fnextra, seed);
		simu->LGSfocusts=dcellnew_file(parms->nwfs, 1, nnx, nny2, "LGS focus time history", "%s/Resfocuserrs_%d.bin", fnextra, seed);
	}
	if(parms->dither){
		simu->dither=mycalloc(nwfs, dither_t*);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].dither){
				simu->dither[iwfs]=mycalloc(1, dither_t);
				if(parms->powfs[ipowfs].dither>1){
					//int nm=parms->powfs[ipowfs].dither_mode2?2:1;
					int nm=NX(P(recon->dither_rg, iwfs, iwfs));
					simu->dither[iwfs]->mr=dcellnewsame_file(2, 1, nm, nsim, "Common path dithering time history", "%s/Resdithererr_wfs%d_%d", fnextra, iwfs, seed);
				}
			}
		}
		if(parms->save.extra||parms->save.dither){
			long nnx[nwfs];
			long nny[nwfs];
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				if(parms->powfs[ipowfs].dither){
					if(parms->powfs[ipowfs].dither==1){
						nnx[iwfs]=4;
					}else{
						int nm=NX(P(recon->dither_rg, iwfs, iwfs));
						nnx[iwfs]=4+2*nm;
					}
					nny[iwfs]=(nsim-parms->powfs[ipowfs].dither_pllskip)/(parms->powfs[ipowfs].dtrat*parms->powfs[ipowfs].dither_pllrat);
				} else{
					nnx[iwfs]=0;
					nny[iwfs]=0;
				}
			}
			simu->resdither=dcellnew_file(nwfs, 1, nnx, nny, "Dithering gain output", "%s/Resditheramp_%d.bin", fnextra, seed);
		}
	}
	if(recon->cn2est){
		int ncn2=(parms->sim.end-1)/parms->cn2.step;
		long nnx[2]={ncn2, NX(recon->cn2est->htrecon)};
		long nny[2]={1, ncn2};
		simu->cn2res=dcellnew_file(2, 1, nnx, nny, "Slodar output", "%s/Rescn2_%d.bin", fnextra, seed);
	}
	if(parms->itpowfs!=-1){
		const int ipowfs=parms->itpowfs;
		const int nacc=parms->sim.end-parms->powfs[ipowfs].step;
		if(nacc>0){
			save->restwfs=zfarr_init(0, 0, "%s/Restwfs_%d.bin", fnextra, seed);
		}
	}
}

static void init_simu_dm(sim_t* simu){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	sim_save_t* save=simu->save;
	const int seed=simu->seed;
	/*Setup hysterisis */
	for(int idm=0; idm<parms->ndm; idm++){
		if(parms->dm[idm].hyst){
			if(!simu->hyst){
				simu->hyst=mycalloc(parms->ndm, hyst_t*);
			}
			simu->hyst[idm]=hyst_new(parms->dm[idm].hyst, parms->dm[idm].hyst_alpha, parms->dm[idm].hyst_stroke, P(recon->aloc, idm)->nloc);
		}
	}
	/*we initialize dmreal, so that wfs_prop_dm can reference dmreal. */
	simu->dmerr_store=dcellnew3(parms->ndm, 1, parms->recon.modal?P(recon->anmod):P(recon->anloc), NULL);
	simu->dmcmd=dcellnew(parms->ndm, 1);
	simu->dmreal=dcellnew(parms->ndm, 1);
	simu->dmrealsq=mapcellnew(parms->ndm, 1);
	if(parms->fit.cgwarm){
		simu->dmrecon=dcellnew3(parms->ndm, 1, parms->recon.modal?P(recon->anmod):P(recon->anloc), NULL);
	} else{
		simu->dmrecon=simu->dmerr_store;
	}
	if(parms->sim.lpttm>EPS){
		simu->ttmreal=dnew(2, 1);
	}

	if(parms->sim.dmproj){
		simu->dmproj=dcellnew(parms->ndm, 1);
		simu->dmprojsq=mapcellnew(parms->ndm, 1);
	}
	for(int idm=0; idm<parms->ndm; idm++){
		P(simu->dmcmd, idm)=dnew(P(recon->anloc, idm), 1);
		//Do not reference for the new synchronization scheme.
		P(simu->dmrealsq, idm)=mapnew2(P(recon->amap, idm));
		dset(P(simu->dmrealsq, idm)->dmat, NAN);
		P(simu->dmrealsq, idm)->dratio=parms->dm[idm].dratio;
		if(parms->fit.square){/*dmreal is also square.*/
			P(simu->dmreal, idm)=dref_reshape(P(simu->dmrealsq, idm)->dmat, P(recon->anloc, idm), 1);
		} else{
			P(simu->dmreal, idm)=dnew(P(recon->anloc, idm), 1);
		}

		if(parms->sim.dmproj){
			P(simu->dmprojsq, idm)=mapnew2(P(recon->amap, idm));
			dset(P(simu->dmprojsq, idm)->dmat, NAN);
			if(parms->fit.square){/*dmreal is also square.*/
				P(simu->dmproj, idm)=dref_reshape(P(simu->dmprojsq, idm)->dmat, P(recon->anloc, idm), 1);
			}else{
				P(simu->dmproj, idm)=dnew(P(recon->anloc, idm), 1);
			}
		}

	}
	if(parms->sim.dmadd){
		simu->dmadd=dcellread("%s", parms->sim.dmadd);
		for(int idm=0; idm<parms->ndm; idm++){
			if(P(simu->dmadd, idm)&&P(simu->dmadd, idm)->nx!=P(simu->dmreal, idm)->nx){
				if(!parms->fit.square||P(simu->dmadd, idm)->nx!=P(recon->amap, idm)->nx*P(recon->amap, idm)->ny){
					error("DM[%d]: dmadd does not match dm geometry\n", idm);
				}
			}
		}
	}
	if(parms->ndm&&(parms->sim.cachedm||parms->plot.run)){
		prep_cachedm(simu);
	}
#if USE_CUDA
	if(parms->gpu.evl||parms->gpu.wfs){
		gpu_dmreal2gpu(simu->dmrealsq);
		if(simu->dmprojsq){
			gpu_dmproj2gpu(simu->dmprojsq);
		}
	}
#endif
	simu->wfspsol=dccellnew(parms->npowfs, 1);
	if(parms->sim.closeloop){
		simu->dmint=servo_new_sho(NULL, parms->sim.aphi, parms->sim.alhi,
			parms->sim.dt, parms->sim.ephi, parms->sim.f0dm, parms->sim.zetadm);
	}
	if(parms->ncpa.preload&&recon->dm_ncpa){//set the integrator
		warning_once("Preload integrator with NCPA\n");
		dcelladd(&P(simu->dmint->mint, 0), 1, recon->dm_ncpa, 1);
	}
	if(parms->recon.split||parms->evl.split>1){
		simu->Merr_lo_store=dcellnew(1, 1);//, recon->ngsmod->nmod, 1);
		simu->Merr_lo2=dcellnew(1, 1);//, recon->ngsmod->nmod, 1);
		if(parms->sim.closeloop){
			simu->Mint_lo=servo_new(NULL,
				parms->sim.aplo, parms->sim.allo, parms->sim.dt, parms->sim.eplo);
		}
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
					nny[idm]=P(recon->aloc, idm)->nloc;
				} else{
					nnx[idm]=0;
					nny[idm]=0;
				}
			}
			simu->dmhist=dcellnew_file(parms->ndm, 1, nnx, nny, "DM position demand time history", "%s/dmhist_%d.bin", fnextra, simu->seed);
		}
	}
	if(parms->recon.psd){
		if(parms->recon.psddtrat_hi){
			simu->dmerrts=dcellnew_same(parms->evl.nevl, 1, P(recon->Herr, 0)->nx,
				parms->recon.psddtrat_hi);
		}
		if((parms->recon.split||parms->evl.split)&&parms->recon.psddtrat_lo){
			simu->Merrts=dnew(recon->ngsmod->nmod, parms->recon.psddtrat_lo);
		}
	}
	int nstep=parms->sim.end;
	int nrstep=nstep-(parms->sim.closeloop?1:0);
	if(parms->save.dm){
		save->dmerr=zfarr_init(nrstep, 1, "dmerr_%d.bin", seed);
		save->dmrecon=zfarr_init(nrstep, 1, "dmrecon_%d.bin", seed);
		if(parms->recon.split||parms->evl.split>1){
			save->Merr_lo=zfarr_init(nstep, 1, "Merr_lo_%d.bin", seed);
			if(!parms->sim.fuseint){
				save->Mint_lo=zfarr_init(nstep, 1, "Mint_lo_%d.bin", seed);
			}
		}
		if(parms->sim.lpttm>EPS){
			save->ttmreal=dnew_file(2, nstep, "Tip/tilt mirror demand time history", "ttmreal_%d.bin", seed);
		}
		save->dmint=zfarr_init(nstep, 1, "dmint_%d.bin", seed);
		save->dmreal=zfarr_init(nstep, 1, "dmreal_%d.bin", seed);
		save->dmcmd=zfarr_init(nstep, 1, "dmcmd_%d.bin", seed);
		if(parms->sim.dmproj){
			save->dmproj=zfarr_init(nstep, 1, "dmproj_%d.bin", seed);
		}
	}
}

/** MOAO **/
static void init_simu_moao(sim_t* simu){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	sim_save_t* save=simu->save;
	const int nwfs=parms->nwfs;
	const int nevl=parms->evl.nevl;
	const int seed=simu->seed;
	int nstep=parms->sim.end;
	int ny=parms->sim.closeloop&&!parms->gpu.moao?2:1;
	if(!parms->gpu.wfs||!parms->gpu.moao||parms->plot.run){
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			int imoao=parms->powfs[ipowfs].moao;
			if(imoao>-1){
				if(!simu->dm_wfs){
					simu->dm_wfs=dcellnew(nwfs, ny);
				}
				for(int iy=0; iy<ny; iy++){
					P(simu->dm_wfs, iwfs, iy)=dnew(P(recon->moao[imoao].aloc, 0)->nloc, 1);
				}
			}
		}
	}
	if(!parms->gpu.evl||!parms->gpu.moao||parms->plot.run){
		int imoao=parms->evl.moao;
		if(imoao!=-1){
			simu->dm_evl=dcellnew(nevl, ny);
			for(int ievl=0; ievl<nevl; ievl++){
				for(int iy=0; iy<ny; iy++){
					P(simu->dm_evl, ievl, iy)=dnew(P(recon->moao[imoao].aloc, 0)->nloc, 1);
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
			save->dm_wfs=mycalloc(nwfs, zfarr*);
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				int imoao=parms->powfs[ipowfs].moao;
				if(imoao>-1){
					save->dm_wfs[iwfs]=zfarr_init(nstep, 1, "wfs%d_moaofit_%d.bin", iwfs, seed);
				}
			}
		}
		if(simu->dm_evl){
			save->dm_evl=mycalloc(nwfs, zfarr*);
			for(int ievl=0; ievl<nevl; ievl++){
				save->dm_evl[ievl]=zfarr_init(nstep, 1, "evl%d_moaofit_%d.bin", ievl, seed);
			}
		}
	}
}

/**
   Initialize simu (of type sim_t) and various simulation data structs. Called
   for every seed.
*/
sim_t* init_simu(const parms_t* parms, powfs_t* powfs,
	aper_t* aper, recon_t* recon, int iseed){
	const int nevl=parms->evl.nevl;
	const int nwfs=parms->nwfs;
	sim_t* simu=mycalloc(1, sim_t);
	global->simu=simu;
	simu->pause=parms->sim.pause;
	pthread_cond_init(&simu->dmreal_condr, 0);
	pthread_cond_init(&simu->dmreal_condw, 0);
	pthread_mutex_init(&simu->dmreal_mutex, 0);
	simu->dmreal_isim=-1;

	pthread_cond_init(&simu->wfsgrad_condr, 0);
	pthread_cond_init(&simu->wfsgrad_condw, 0);
	pthread_mutex_init(&simu->wfsgrad_mutex, 0);
	simu->wfsgrad_isim=-1;

	sim_save_t* save=simu->save=mycalloc(1, sim_save_t);
	simu->wfsisim=-1;
	simu->perfisim=-1;
	simu->reconisim=-2;
	simu->parms=parms;
	simu->powfs=powfs;
	simu->recon=recon;
	simu->aper=aper;
	simu->iseed=iseed;
	simu->seed=P(parms->sim.seeds, iseed);
	if(simu->seed==0){
		simu->seed=myclocki();
		warning("Seed is 0, using time as seed: %d\n", simu->seed);
	}
	int seed=simu->seed;
	seeding(simu);

	simu->status=mycalloc(1, status_t);
	simu->status->iseed=iseed;
	simu->status->nseed=parms->sim.nseed;
	simu->status->simstart=parms->sim.start;
	simu->status->simend=parms->sim.end;
	simu->status->nthread=NTHREAD;
	simu->status->info=S_RUNNING;
	if(parms->save.extra||parms->save.dither){
		mymkdir("extra");
	}
	fnextra=parms->save.extra?"extra":"-";
	if(parms->sim.wspsd){
		if(parms->sim.idealtomo){
			warning("sim.idealtomo is not yet implemented for sim.wspsd. Ignored\n");
		}else{
			/* Telescope wind shake added to TT input. */
			info("Converting windshake PSD to time series.\n");
			simu->telws=psd2ts(parms->sim.wspsd, simu->telws_rand, parms->sim.dt, parms->sim.end);
			if(parms->save.extra) writebin(simu->telws, "%s/telws_%d", fnextra, seed);
		}
	}

	/* Select GPU or CPU for the tasks.*/
#if USE_CUDA
	if(parms->gpu.evl){
		simu->perfevl_pre=thread_prep(0, nevl, nevl, gpu_perfevl_queue, simu);
		simu->perfevl_post=thread_prep(0, nevl, nevl, gpu_perfevl_sync, simu);
	} else
#endif
	{
		simu->perfevl_pre=thread_prep(0, nevl, nevl, perfevl_ievl, simu);
	}

#if USE_CUDA
	if(parms->gpu.wfs){
		simu->wfsgrad_pre=thread_prep(0, nwfs, nwfs, gpu_wfsgrad_queue, simu);
	} else
#endif
	{
		simu->wfsgrad_pre=thread_prep(0, nwfs, nwfs, wfsgrad_iwfs, simu);
	}
	{
		int nthread=2;
		simu->wfsgrad_post=thread_prep(0, nwfs, nthread, wfsgrad_post, simu);
	}
	if(!parms->sim.evlol){
		init_simu_dm(simu);
		init_simu_moao(simu);
		if(!parms->sim.idealtomo){
			init_simu_wfs(simu);
		}
		if(parms->recon.alg==RECON_MVR){
			int nstep=parms->sim.end;
			if(parms->save.opdr){
				save->opdr=zfarr_init(nstep-1, 1, "opdr_%d.bin", seed);
			}
			if(parms->save.opdx){
				save->opdx=zfarr_init(nstep, 1, "opdx_%d.bin", seed);
			}
		}
		{
			const char* keywords="CG residual for Tomography; DM Fit";
			simu->cgres=dcellnewsame_file(2, 1, parms->sim.end, 1, keywords, "%s/ResCG_%d.bin", fnextra, seed);
		}
		if(parms->recon.psd&&parms->save.extra){
			if(parms->recon.psddtrat_hi){
				save->psdcl=zfarr_init(0, 0, "%s/Respsdcl_hi_%d.bin", fnextra, seed);
				save->psdol=zfarr_init(0, 0, "%s/Respsdol_hi_%d.bin", fnextra, seed);
			}
			if(parms->recon.psddtrat_lo){
				save->psdcl_lo=zfarr_init(0, 0, "%s/Respsdcl_lo_%d.bin", fnextra, seed);
				save->psdol_lo=zfarr_init(0, 0, "%s/Respsdol_lo_%d.bin", fnextra, seed);
			}
		}
	}
	init_simu_evl(simu);
#if USE_CUDA
	if(parms->gpu.evl||parms->gpu.wfs){
		if(parms->gpu.evl){
			gpu_perfevl_init_sim(parms, aper);
		}
		if(parms->gpu.wfs&&!parms->sim.evlol){
			gpu_wfs_init_sim(parms, powfs);
		}
	}
	if(parms->gpu.tomo||parms->gpu.fit){
		gpu_recon_reset(parms);
	}
#endif

	filter_dm(simu);//2014-03-31. //so that dm_ncpa is effective at first cycle. replaced by copy dm_ncpa to dmreal.
	return simu;
}
/**
   Release memory of simu (of type sim_t) and close files.
*/
void free_simu(sim_t* simu){
	if(!simu) return;
	global->simu=0;
	const parms_t* parms=simu->parms;
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
	if(simu->cachedm){
		cellfree(simu->cachedm);
		for(int i=0; i<parms->ndm; i++){
			free(simu->cachedm_prop[i]);
		}
		free(simu->cachedm_prop);
		free(simu->cachedm_propdata);

	}
	if(parms->nlgspowfs){
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			sho_free(simu->ltpm_sho[ipowfs]);
		}
		free(simu->ltpm_sho);
	}
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		for(int ips=0; ips<parms->atm.nps; ips++){
			free(simu->wfs_prop_atm[iwfs+nwfs*ips]);
		}
		for(int idm=0; idm<parms->ndm; idm++){
			free(simu->wfs_prop_dm[iwfs+nwfs*idm]);
		}
		free(simu->wfs_ints[iwfs]);
		if(simu->fsmint) servo_free(simu->fsmint[iwfs]);
		if(simu->fsmsho) sho_free(simu->fsmsho[iwfs]);
	}
	free(simu->fsmint);
	free(simu->fsmsho);
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
	dfree(simu->evlopdground);
	dcellfree(simu->wfsopd);
	dcellfree(simu->gradcl);
	dcellfree(simu->gradacc);
	dcellfree(simu->gradlastcl);
	dcellfree(simu->gradlastol);
	dcellfree(simu->gradoff);
	dcellfree(simu->gradoffacc);
	dcellfree(simu->gradoffdrift);
	dcellfree(simu->gradscale);
	dcellfree(simu->gradscale2);
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
	cellfree(simu->Merrts);
	dcellfree(simu->gcov);
	dcellfree(simu->ecov);
	if(simu->dmrecon!=simu->dmerr_store){
		dcellfree(simu->dmrecon);
	}
	dcellfree(simu->dmerr_store);
	dcellfree(simu->dmhist);
	dcellfree(simu->Merr_lo_store);
	dcellfree(simu->Merr_lo2);
	dcellfree(simu->ngsmodlpf);
	dcellfree(simu->fsmerr_store);
	dcellfree(simu->fsmerr_drift);
	dcellfree(simu->fsmreal);
	dcellfree(simu->fsmcmd);
	dcellfree(simu->cgres);
	dcellfree(simu->dm_wfs);
	dcellfree(simu->dm_evl);
	dcellfree(simu->res);
	//dfree(simu->timing);
	dcellfree(simu->resp);
	dcellfree(simu->olep);
	dcellfree(simu->olmp);
	dcellfree(simu->clep);
	dcellfree(simu->clmp);

	dcellfree(simu->LGSfocus);
	dcellfree(simu->LGSfocus_drift);
	dcellfree(simu->LGSfocusts);
	dcellfree(simu->telfocusint);
	dcellfree(simu->telfocusreal);
	dfree(simu->zoomdrift);
	lfree(simu->zoomdrift_count);
	dfree(simu->zoomavg);
	lfree(simu->zoomavg_count);
	dfree(simu->zoomint);
	if(parms->evl.split){
		dcellfree(simu->clemp);
		dfree(simu->cleNGSm);
		dfree(simu->oleNGSm);
		dfree(simu->corrNGSm);
		dcellfree(simu->cleNGSmp);
		dcellfree(simu->oleNGSmp);
	}
	cellfree(simu->petal_i0);
	cellfree(simu->petal_m);
	sim_save_t* save=simu->save;
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
	zfarr_close_n(save->evlopdmean_ngsr, nevl);
	zfarr_close(save->evlpsfolmean);
	zfarr_close(save->psdcl);
	zfarr_close(save->psdol);
	zfarr_close(save->psdcl_lo);
	zfarr_close(save->psdol_lo);
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

	cellfree(simu->zoompos);
	cellfree(simu->llt_ws);
	free(simu->plot_legs);
	cellfree(simu->plot_res);

	/*Close all files */
	zfarr_close(save->restwfs);
	zfarr_close_n(save->wfspsfout, nwfs);
	zfarr_close_n(save->ztiltout, nwfs);
	zfarr_close(save->dmerr);
	zfarr_close(save->dmint);
	zfarr_close(save->dmrecon);
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
	dcellfree(save->fsmerrs);
	dcellfree(save->fsmcmds);
	dcellfree(simu->ltpm_cmd);
	dcellfree(simu->ltpm_lpf);
	dcellfree(simu->ltpm_real);
	dcellfree(save->ltpm_real);
	free(simu->wfsflags);
	dfree(simu->winddir);
	free(simu->save);
	free(simu);
}

int compare_dbl2_ascend(const void *a, const void *b){
	return (int)(((double *)a)[0]-((double *)b)[0]);
}
/**
   Print out wavefront error information and timing at each time step.
*/
void print_progress(sim_t* simu){
	const int isim=simu->perfisim;
	const int iseed=simu->iseed;
	const parms_t* parms=simu->parms;
	const int simstart=parms->sim.start;
	const int simend=parms->sim.end;
	long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
	long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
	simu->status->tot=simu->tk_iend-simu->tk_istart;//total step time
	if(isim==simstart){//first step, rough estimate.
		simu->status->mean=simu->status->tot;
		simu->status->rest=simu->status->mean*parms->sim.nseed*(simend-simstart);
	} else{
		simu->status->rest=(long)((simu->tk_iend-simu->tk_s0-(simu->tk_i1-simu->tk_si)*(iseed+1))/steps_done*steps_rest
			+(simu->tk_i1-simu->tk_si)*(parms->sim.nseed-iseed-1));
		simu->status->mean=(simu->tk_iend-simu->tk_i1)/(real)(isim-simstart);
	}
	simu->status->laps=(long)(simu->tk_iend-tk_setup);//total elapsed time

	simu->status->wfs=simu->tk_wfs;
	simu->status->recon=simu->tk_recon;
	simu->status->other=simu->tk_cache;
	simu->status->eval=simu->tk_eval;
	simu->status->scale=1;
/*if(simu->timing){
		P(simu->timing, 0, isim)=get_job_mem();
		P(simu->timing, 1, isim)=simu->status->tot;
		P(simu->timing, 2, isim)=simu->status->wfs;
		P(simu->timing, 3, isim)=simu->status->recon;
		P(simu->timing, 4, isim)=simu->status->eval;
	}
*/
	real this_time=myclockd();
	if(simu->res&&simu->res->fp){//save res periodically for plotting.
		static real last_save_time=0;
		const real gap=isim<10000?10:60;
		if(this_time>last_save_time+gap){
			writebin_async(simu->res, isim+1);
			if(parms->save.extra){
				writebin_async(simu->resp, isim+1);
				//writebin_async(simu->restwfs, isim+1);//column is different
				//writebin_async(simu->resdither, isim+1);//column is different
				writebin_async(simu->save->fsmerrs, simu->wfsisim+1);
				writebin_async(simu->save->fsmcmds, simu->wfsisim+1);
				writebin_async(simu->save->ltpm_real, simu->wfsisim+1);
				if(parms->nlgspowfs){
					writebin_async(simu->LGSfocusts, simu->wfsisim+1);
					if(simu->zoompos_icol){
						writebin_async(simu->zoompos, simu->zoompos_icol);
					}
				}
			}
			last_save_time=this_time;
		}
	}
	static int nstep=1;
	if(isim%nstep==0||isim+1==parms->sim.end){
		if(isim>2){
			int nstep2=ceil(1./simu->status->tot);
			if(nstep2>5){
				nstep2=round(nstep2*0.1)*10;
			}
			if(nstep2>nstep) nstep=nstep2;
		}
	/*we don't print out or report too frequently. */
		simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
		const status_t* status=simu->status;
		const real tkmean=status->scale;
		const long rest=status->rest;
		const long laps=status->laps;
		const long resth=rest/3600;
		const long restm=(rest-resth*3600)/60;
		const long lapsh=laps/3600;
		const long lapsm=(laps-lapsh*3600)/60;
		if(isim==parms->sim.start){
			const char *hol="Step     #  Open Loop  PR TT  ";
			const char *hoa="On Axis PR    TT Field  PR     TT ";
			const char *hsp="Split High     TT    Low    ";
			const char *htm=" Timing ";
			info2("%s%s%s%s\n", hol, !parms->sim.evlol?hoa:"", (!parms->sim.evlol&&parms->evl.split)?hsp:"", LOG_LEVEL<1?htm:"");
		}
		info2("Step %5d: OL%7.1f %6.1f",
			isim,
			mysqrt(P(simu->ole, 0, isim))*1e9,
			mysqrt(P(simu->ole, 1, isim))*1e9);
		if(!parms->sim.evlol){
			info2("  OA%7.1f %6.1f CL%7.1f %6.1f",
				mysqrt(P(simu->clep, 0, 0, 0, isim))*1e9,
				mysqrt(P(simu->clep, 0, 0, 1, isim))*1e9,
				mysqrt(P(simu->cle, 0, isim))*1e9,
				mysqrt(P(simu->cle, 1, isim))*1e9
			);
			if(parms->evl.split){
				info2(" ST %7.1f %6.1f %6.1f",
					mysqrt(P(simu->clem, 0, isim))*1e9,
					mysqrt(P(simu->clem, 1, isim))*1e9,
					mysqrt(P(simu->clem, 2, isim))*1e9);
			}
			extern int NO_EVL;
			if(!NO_EVL&&isinf(P(simu->cle, 0, isim))){
				error("\nStep %d: NaN/inf found: cle is %g\n", isim, P(simu->cle, 0, isim));
			}
		}
		if(LOG_LEVEL<1){//dbg is inactive
			info2(" nm %6.3fs\n", status->mean*tkmean);
		}else{
			info2(" nm\n");
		}

		dbg("Timing: WFS %.3f Recon %.3f EVAL %.3f Other %.3f Tot %.3f Mean %.3f."
					" Used %ld:%02ld, Left %ld:%02ld\n",
					status->wfs*tkmean, status->recon*tkmean,
					status->eval*tkmean, status->other*tkmean,
					status->tot*tkmean, status->mean*tkmean,
					lapsh, lapsm, resth, restm);

		if(parms->plot.run){
			if(!simu->plot_legs){
				simu->plot_legs=mycalloc(5, const char*);
				const char **legs=simu->plot_legs;
				int nline=0;
				legs[nline++]="Total";
				legs[nline++]="High Order";
				legs[nline++]="Tip/Tilt";
				if(parms->evl.split){
					if(simu->recon->ngsmod->indps){
						legs[nline++]="Plate Scale";
					} else if(simu->recon->ngsmod->indastig){
						legs[nline++]="Astig";
					}
					if(simu->recon->ngsmod->indfocus){
						legs[nline++]="Focus";
					}
				}
				simu->plot_res=dccellnew(5,1);
				P(simu->plot_res, 0)=dcellnew_same(nline, 1, parms->sim.end, 1);//CL vs time step for different modes
				P(simu->plot_res, 1)=dcellnew_same(3, 1, parms->sim.end, 1);//OL vs time step
				P(simu->plot_res, 2)=dcellnew_same(1, 1, parms->evl.nevl, 2);//CL vs dir WEF
				P(simu->plot_res, 3)=dcellnew_same(1, 1, parms->evl.nevl, 1);//CL vs dir indexing
				P(simu->plot_res, 4)=dcellnew_same(parms->evl.nwvl, 1, parms->evl.nevl, 2);//CL vs dir Strehl
				if(parms->evl.nevl>1){	// sort evaluation directions
					//\todo modularize this.
					dmat *tmp=dnew(2, parms->evl.nevl);
					for(int ievl=0; ievl<parms->evl.nevl; ievl++){
						real r=RSS(P(parms->evl.thetax, ievl), P(parms->evl.thetay, ievl))*RAD2AS;
						/*if(P(parms->evl.thetax, ievl)<0 || P(parms->evl.thetay, ievl)<0){
							r=-r;
						}*/
						P(tmp, 0, ievl)=r;
						P(tmp, 1, ievl)=ievl;//mark direction
					}
					//dshow(tmp,"tmp");
					qsort(P(tmp), NY(tmp), 2*sizeof(real), compare_dbl2_ascend);
					//dshow(tmp,"tmp sorted");
					for(int ievl=0; ievl<parms->evl.nevl; ievl++){
						P(P(P(simu->plot_res, 2), 0), ievl, 0)=P(tmp, 0, ievl);
						P(P(P(simu->plot_res, 3), 0), ievl, 0)=P(tmp, 1, ievl);//index
						for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
							P(P(P(simu->plot_res, 4), iwvl), ievl, 0)=P(tmp, 0, ievl);
						}
					}
					dfree(tmp);
				}
				//dset(simu->plot_res->m, NAN);
			}
			if(isim>20&&parms->evl.nevl>1&&(draw_current("Res", "FoV WFE")||draw_current("Res", "FoV Strehl"))){
				dcell *res=P(simu->plot_res, 2);
				int istart=MAX(isim-1000, 20);
				for(int ievl=0; ievl<parms->evl.nevl; ievl++){
					int jevl=(int)P(P(P(simu->plot_res,3), 0), ievl);
					P(P(res,0), ievl, 1)=0;
					for(int i=istart; i<isim; i++){
						P(P(res,0), ievl, 1)+=P(P(simu->clep, jevl),0,i);
					}
					P(P(res, 0), ievl, 1)=sqrt(P(P(res, 0), ievl, 1)/(isim-istart))*1e9;
				}
				//check directions of the same radius and average them.
				for(int ievl=0; ievl<parms->evl.nevl-1; ievl++){
					int jevl=ievl+1;
					for(; jevl<parms->evl.nevl; jevl++){
						if(fabs(P(P(res, 0), ievl, 0)-P(P(res, 0), jevl, 0))>1){
							break;
						}
					}
					if(jevl>ievl+1){//directions of the same radius
						real sum2=0;
						for(int i=ievl; i<jevl; i++){
							sum2+=pow(P(P(res, 0), i, 1),2);
						}
						sum2=sqrt(sum2/(jevl-ievl));
						for(int i=ievl; i<jevl; i++){
							P(P(res, 0), i, 1)=sum2;
						}
					}
				}
				for(int ievl=0; ievl<parms->evl.nevl; ievl++){
					real wfe=P(P(P(simu->plot_res, 2), 0), ievl, 1)*1e-9;
					for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
						P(P(P(simu->plot_res, 4), iwvl), ievl, 1)=exp(-pow(TWOPI*wfe/P(parms->evl.wvl,iwvl),2));
					}
				}
				draw("Res", (plot_opts){ .dc=res},
						"Wavefront Error", "Field Angle (as)", "Wavefront Error (nm)", "FoV WFE");
				draw("Res", (plot_opts){ .dc=P(simu->plot_res, 4), .legend=parms->evl.wvlname},
										"Strehl Ratio", "Field Angle (as)", "Strehl Ratio", "FoV Strehl");
			}
			if(draw_current("Res", "RMS WFE") || draw_current("Res", "RMS WFE OL") ){
				for(;simu->plot_isim<=isim;simu->plot_isim++){
					int i=simu->plot_isim;
					for(int ic=0; ic<2; ic++){//0: CL. 1: OL
						dmat * cle=ic==0?simu->cle:simu->ole;
						dcell* res=P(simu->plot_res,ic);
					
						P(P(res, 0), i)=sqrt(P(cle, 0, i))*1e9;//PR
						if(parms->evl.split && ic==0){
							dmat *tmp=simu->clem;
							P(P(res, 1), i)=sqrt(P(tmp, 0, i))*1e9;//LGS
							P(P(res, 2), i)=sqrt(P(tmp, 1, i))*1e9;//TT
							if(NX(res)>4){
								P(P(res, 4), i)=sqrt(P(tmp, 3, i))*1e9;//Focus
							}
							if(NX(res)>3){
								P(P(res, 3), i)=sqrt(P(tmp, 2, i)-(P(tmp, 3, i)+P(tmp, 1, i)))*1e9;//PS
							}
						} else{
							P(P(res, 1), i)=sqrt(P(cle, 2, i))*1e9;//PTTR
							P(P(res, 2), i)=sqrt(P(cle, 0, i)-P(cle, 2, i))*1e9;//TT
						}
					}
				}
				for(int ic=0; ic<2; ic++){
					dcell *res=P(simu->plot_res, ic);
					draw("Res", (plot_opts){.ngroup=NX(res), .maxlen=simu->plot_isim+1, .dc=res, .xylog="nn", .legend=simu->plot_legs},
						"Wavefront Error", "Time Step", "Wavefront Error (nm)", "RMS WFE%s", ic==0?"":" OL");
				}
			}
		}
	}
}
/**
   Output parameters necessary to run postproc using skyc/skyc.c
*/
void save_skyc(powfs_t* powfs, recon_t* recon, const parms_t* parms){
	char fn[PATH_MAX];
	real zadeg=parms->sim.za*180/M_PI;
	snprintf(fn, PATH_MAX, "%s/maos.conf", dirskysim);
	FILE* fp=fopen(fn, "w");
	fprintf(fp, "maos.r0z=%g\n", parms->atm.r0z);
	fprintf(fp, "maos.dt=%g\n", parms->sim.dt);
	fprintf(fp, "maos.zadeg=%g\n", zadeg);
	fprintf(fp, "maos.hc=%g\n", parms->dm[parms->ndm-1].ht);
	fprintf(fp, "maos.hs=%g\n", recon->ngsmod->hs);
	fprintf(fp, "maos.nmod=%d\n", recon->ngsmod->nmod);
	fprintf(fp, "maos.D=%g\n", parms->aper.d);
	fprintf(fp, "maos.wvl=[");
	int nwvl=0;
	int npowfs_ngs=0;
	int powfs_ngs[parms->npowfs];
	real ngsgrid=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].lo){
			powfs_ngs[npowfs_ngs++]=ipowfs;
			if(nwvl==0){
				nwvl=parms->powfs[ipowfs].nwvl;
			} else{
				if(nwvl!=parms->powfs[ipowfs].nwvl){
					error("Different low order WFS has different wvl\n");
				}
			}
			real sepmean=0;
			for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs-1; indwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, indwfs);
				int iwfs2=P(parms->powfs[ipowfs].wfs, indwfs+1);
				if(fabs(parms->wfs[iwfs].thetax-parms->wfs[iwfs2].thetax)<1.e-10){
					real sep=parms->wfs[iwfs2].thetay-parms->wfs[iwfs].thetay;
					if(fabs(sepmean)<1e-20){
						sepmean=sep;
					} else if(fabs(sep-sepmean)>1.e-10){
						sepmean=0;
						warning("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
					}
				}
				if(fabs(parms->wfs[iwfs].thetay-parms->wfs[iwfs2].thetay)<1.e-10){
					real sep=parms->wfs[iwfs2].thetax-parms->wfs[iwfs].thetax;
					if(fabs(sepmean)<1e-20){
						sepmean=sep;
					} else if(fabs(sep-sepmean)>1.e-10){
						sepmean=0;
						warning("NGS WFS are not evenly spaced. Unable to determine ngs spacing.\n");
					}
				}
			}
			if(fabs(sepmean)>1.e-20){
				if(fabs(ngsgrid)<1.e-20){
					ngsgrid=sepmean;
				} else{
					if(fabs(sepmean-ngsgrid)>1.e-10){
						error("Different NGS POWFS have different spacing\n");
					}
				}
			}
		}
	}
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, "%g ", P(parms->powfs[powfs_ngs[0]].wvl, iwvl));
	}
	fprintf(fp, "]\n");
	fprintf(fp, "maos.ngsgrid=%g\n", ngsgrid>0?ngsgrid*RAD2AS:1);
	fprintf(fp, "maos.npowfs=%d\n", npowfs_ngs);
	fprintf(fp, "maos.msa=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "%d ", parms->powfs[powfs_ngs[ipowfs]].order);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "maos.nsa=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "%ld ", powfs[powfs_ngs[ipowfs]].saloc->nloc);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "maos.ncomp=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "%d ", powfs[powfs_ngs[ipowfs]].notfx);
		if(powfs[powfs_ngs[ipowfs]].notfx!=powfs[powfs_ngs[ipowfs]].notfy){
			error("Invalid\n");
		}
	}
	fprintf(fp, "]\n");
	fprintf(fp, "maos.embfac=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "%d ", parms->powfs[powfs_ngs[ipowfs]].embfac);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "maos.dxsa=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "%g ", powfs[powfs_ngs[ipowfs]].saloc->dx);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "maos.fnwfsloc=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "\"powfs%d_loc\" ", powfs_ngs[ipowfs]);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "maos.fnwfsamp=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "\"powfs%d_amp\" ", powfs_ngs[ipowfs]);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "maos.fnsaloc=[");
	for(int ipowfs=0; ipowfs<npowfs_ngs; ipowfs++){
		fprintf(fp, "\"powfs%d_saloc\" ", powfs_ngs[ipowfs]);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "maos.fnmideal=\"RescleNGSm\"\n");
	fprintf(fp, "maos.fnmidealp=\"RescleNGSmp\"\n");
	fprintf(fp, "maos.evlindoa=%d\n", parms->evl.indoa);
	fprintf(fp, "maos.fnmcc=\"MCC_za%g\"\n", zadeg);
	fprintf(fp, "maos.fnmcc_oa=\"MCC_OA_za%g\"\n", zadeg);

	fprintf(fp, "maos.seeds=[");
	for(int iseed=0; iseed<parms->sim.nseed; iseed++){
		fprintf(fp, "%ld ", P(parms->sim.seeds, iseed));
	}
	fprintf(fp, "]\n");

	fprintf(fp, "include=\"skyc.conf\"\n");
	fprintf(fp, "include=\"skyc_za%g.conf\"\n", zadeg);
	fprintf(fp, "maos.wddeg=[");
	for(int ips=0; ips<parms->atm.nps; ips++){
		fprintf(fp, "%.2f ", parms->atm.wddeg?P(parms->atm.wddeg, ips):0);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "maos.nstep=%d\n", parms->sim.end);
	fprintf(fp, "maos.ahstfocus=%d\n", parms->tomo.ahst_focus);
	fprintf(fp, "maos.mffocus=%d\n", parms->sim.mffocus);
	fprintf(fp, "maos.fnrange=%s\n", parms->powfs[P(parms->hipowfs, 0)].llt->fnrange);
	fprintf(fp, "maos.indps=%d\n", recon->ngsmod->indps);
	fprintf(fp, "maos.indastig=%d\n", recon->ngsmod->indastig);
	fprintf(fp, "maos.indfocus=%d\n", recon->ngsmod->indfocus);
	fclose(fp);
	for(int jpowfs=0; jpowfs<npowfs_ngs; jpowfs++){
		int ipowfs=powfs_ngs[jpowfs];
		locwrite(powfs[ipowfs].loc, "%s/powfs%d_loc", dirskysim, ipowfs);
		writebin(P(powfs[ipowfs].amp,0), "%s/powfs%d_amp", dirskysim, ipowfs);
		locwrite(powfs[ipowfs].saloc, "%s/powfs%d_saloc", dirskysim, ipowfs);
		if(powfs[ipowfs].gradncpa){
			int nsa=parms->powfs[ipowfs].order;
			mymkdir("%s/gradoff/", dirskysim);
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				writebin(P(powfs[ipowfs].gradncpa, jwfs),
					"%s/gradoff/gradoff_sa%d_x%.0f_y%.0f.bin",
					dirskysim, nsa,
					parms->wfs[iwfs].thetax*RAD2AS,
					parms->wfs[iwfs].thetay*RAD2AS);
			}
		}
	}
	writebin(recon->ngsmod->MCC, "%s/MCC_za%g", dirskysim, zadeg);
	writebin(P(recon->ngsmod->MCCP, parms->evl.indoa), "%s/MCC_OA_za%g", dirskysim, zadeg);
}
