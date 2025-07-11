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
#include "common.h"
#include "surf.h"
#include "recon_utils.h"
#include "powfs.h"

/**
   \file surf.h
   Setup NCPA surfaces for WFS and Performance evaluation.
*/

/**
   Setup tilted surface (M3) by ray tracing from the tilted surface to WFS and
   Science grid.

   \todo Merge setup_tsurf, setup_surf, and setup_powfs_ncpa. Be careful about
   which NCPA is compensated by WFS offsets. setup_tsurf and setup_surf act on
   surfaces intersect with most WFS and science, while setup_powfs_ncpa act on
   individual surfaces for each WFS. Maybe it is good idea to keep them
   separate.


   2011-09-07:
   Relocate simu->surfwfs to powfs.opdadd In order to perserve surface for different seeds

   Treated as common path aberration
*/

static void
setup_surf_tilt(const parms_t* parms, aper_t* aper, powfs_t* powfs, recon_t* recon){
	info2("Setting up tilt surface (M3)\n");
	rmapcell* tsurf=rmapcellnew( parms->ncpa.ntsurf, 1);
	for(int itsurf=0; itsurf< parms->ncpa.ntsurf; itsurf++){
		char* fn= parms->ncpa.tsurf[itsurf];
		info2("Loading tilt surface from %s\n", fn);
		P(tsurf, itsurf)=rmapread("%s", fn);
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		for(int itsurf=0; itsurf< parms->ncpa.ntsurf; itsurf++){
			m3proj(P(tsurf, itsurf), P(aper->opdadd, ievl), aper->locs,
				P(parms->evl.thetax, ievl), P(parms->evl.thetay, ievl), P(parms->evl.hs, ievl));
		}
	}
	if( parms->ncpa.calib){
		for(int ifit=0; ifit< parms->ncpa.ndir; ifit++){
			for(int itsurf=0; itsurf< parms->ncpa.ntsurf; itsurf++){
				m3proj(P(tsurf, itsurf), P(aper->opdbias, ifit), recon->floc,
					P(parms->ncpa.thetax, ifit), P(parms->ncpa.thetay, ifit), P(parms->ncpa.hs, ifit));
			}
		}
	}

	for(int iwfs=0; iwfs<parms->nwfs&&powfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
		loc_t* locwfs;
		if(powfs[ipowfs].loc_tel){
			warning("We don't handle this case yet. Think carefully when to apply shift.\n");
			locwfs=P(powfs[ipowfs].loc_tel, wfsind);
		} else{
			locwfs=powfs[ipowfs].loc;
		}
		for(int itsurf=0; itsurf< parms->ncpa.ntsurf; itsurf++){
			m3proj(P(tsurf, itsurf), P(powfs[ipowfs].opdadd, wfsind), locwfs,
				parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, parms->wfs[iwfs].hs);
		}
	}
	cellfree(tsurf);
}

typedef struct{
	const parms_t* parms;
	aper_t* aper;
	powfs_t* powfs;
	recon_t* recon;
	map_t* surf;
	int isurf;
	int nevl;
	int nwfs;
	int nncpa;
	int isncpa;
	int* evlcover;
	int* wfscover;
	int* ncpacover;
	int opdxcover;
}SURF_DATA;

static void* prop_surf_evl(thread_t* info){
	SURF_DATA* data=(SURF_DATA*)info->data;
	const parms_t* parms=data->parms;
	const aper_t* aper=data->aper;
	const map_t* surf=data->surf;
	const real hl=surf->h;
	const int* evlcover=data->evlcover;
	for(int ievl=info->start; ievl<info->end; ievl++){
		if(!evlcover[ievl]){
			continue;
		}
		const real displacex=P(parms->evl.thetax, ievl)*hl;
		const real displacey=P(parms->evl.thetay, ievl)*hl;
		const real scale=1-hl/P(parms->evl.hs, ievl);
		prop_grid_stat(surf, aper->locs->stat, P(P(aper->opdadd, ievl)),
			1, displacex, displacey, scale, 0, 0, 0);
	}
	return NULL;
}

static void* prop_surf_ncpa(thread_t* info){
	SURF_DATA* data=(SURF_DATA*)info->data;
	const parms_t* parms=data->parms;
	const aper_t* aper=data->aper;
	const recon_t* recon=data->recon;
	const map_t* surf=data->surf;
	const real hl=surf->h;
	const int* ncpacover=data->ncpacover;
	if(!data->isncpa){
		dbg2("Ignore CPA surface %d for aper->opdbias\n", data->isurf);
		return NULL;
	}
	for(int idir=info->start; idir<info->end; idir++){
		if(!ncpacover[idir]) continue;
		if(!data->isncpa){
			dbg("Ignore CPA surface %d for aper->opdbias\n", data->isurf);
			continue;
		}
		const real displacex=P(parms->ncpa.thetax, idir)*hl;
		const real displacey=P(parms->ncpa.thetay, idir)*hl;
		const real scale=1.-hl/P(parms->ncpa.hs, idir);
		prop_grid(surf, recon->floc, P(P(aper->opdbias, idir)),
			1, displacex, displacey, scale, 0, 0, 0);
	}
	return NULL;
}

static void* prop_surf_wfs(thread_t* info){
	SURF_DATA* data=(SURF_DATA*)info->data;
	const parms_t* parms=data->parms;
	const powfs_t* powfs=data->powfs;
	const map_t* surf=data->surf;
	const real ht=surf->h;
	const int* wfscover=data->wfscover;
	for(int iwfs=info->start; iwfs<info->end; iwfs++){
		if(!wfscover[iwfs]){
			continue;
		}
		const int ipowfs=parms->wfs[iwfs].powfs;
		const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
		const real hs=parms->wfs[iwfs].hs;
		const real scale=1.-ht/hs;
		const real displacex=parms->wfs[iwfs].thetax*ht;
		const real displacey=parms->wfs[iwfs].thetay*ht;

		loc_t* locwfs;
		if(powfs[ipowfs].loc_tel){
			locwfs=P(powfs[ipowfs].loc_tel, wfsind);
		} else{
			locwfs=powfs[ipowfs].loc;
		}
		dcell *opdout=data->isncpa?powfs[ipowfs].opdbias:powfs[ipowfs].opdadd;
		prop_grid(surf, locwfs, P(P(opdout, wfsind)),
			1, displacex, displacey, scale, 1., 0, 0);
	}
	return NULL;
}

/**
   Setup surface perpendicular to the beam by ray tracing from the surface to
   WFS and Science grid
*/
static void
setup_surf_perp(const parms_t* parms, aper_t* aper, powfs_t* powfs, recon_t* recon){
	info2("Setting up surface OPD\n");

	const int nevl=parms->evl.nevl;
	const int nwfs=parms->nwfs;
	const int nncpa= parms->ncpa.ndir;
	int* evlcover=mymalloc(nevl, int);
	int* wfscover=mymalloc(nwfs, int);
	int* ncpacover=mymalloc(nncpa, int);
	int opdxcover=1;
	SURF_DATA sdata={parms, aper, powfs, recon, NULL, 0, nevl, nwfs, nncpa, 0, evlcover, wfscover, ncpacover, 0};
	const int nthread=NTHREAD;
	thread_t* tdata_evl=thread_prep(0, nevl, nthread, prop_surf_evl, &sdata);
	thread_t* tdata_wfs=thread_prep(0, nwfs, nthread, prop_surf_wfs, &sdata);
	thread_t* tdata_ncpa=nncpa?thread_prep(0, nncpa, nthread, prop_surf_ncpa, &sdata):NULL;

	for(int isurf=0; isurf< parms->ncpa.nsurf; isurf++){
		char* fn= parms->ncpa.surf[isurf];
		if(!fn) continue;
		info("\nSurface %d:\n", isurf);
		mapcell* surfs=genscreen_str(fn);
		if(!surfs) continue;
		for(int isurf2=0; isurf2<PN(surfs); isurf2++){
			map_t* surf=P(surfs, isurf2);
			if(parms->plot.setup>1){
				draw("Surf", (plot_opts){.image=surf->dmat}, "Surface", "x", "y", "surf %d %d", isurf, isurf2);
			}
			const char* strname=search_keyword(surf->keywords, "SURFNAME");
			const char* strevl=search_keyword(surf->keywords, "SURFEVL");
			const char* strwfs=search_keyword(surf->keywords, "SURFWFS");
			const char* stropdx=search_keyword(surf->keywords, "SURFOPDX");
			//pupil rotation. rotate the surface directly
			if(strname&&(!strcmp(strname, "M1")||!strcmp(strname, "M2"))){
				if(parms->aper.rot){
					dmaprot(surf, parms->aper.rot);
				}
			}
			//pupil misregistration
			if(strname&&(!strcmp(strname, "M1"))){
				surf->ox+=P(parms->aper.misreg,0);
				surf->oy+=P(parms->aper.misreg,1);
			}
			int evlct=0;
			int ncover=0;//NCPA if 1, CPA if 2
			if(!strevl){
				info("Does not contain SURFEVL. Assume it covers all evaluation directions.\n");
				for(int ievl=0; ievl<nevl; ievl++){
					evlcover[ievl]=1;
				}
				evlct=nevl;
			} else{
				readstr_numarr((void**)&evlcover, NULL, NULL, nevl, 2, M_INT, "SURFEVL",strevl);
				for(int ievl=0; ievl<nevl; ievl++){
					evlct+=evlcover[ievl]==1?1:0;
				}
				if(evlct>0){
					info2("Covers evl");
					for(int ievl=0; ievl<nevl; ievl++){
						if(evlcover[ievl]){
							info2(" %d", ievl);
						}
					}
					info2("\n");
				}
			}
			if(evlct!=nevl){
				for(int idir=0; idir<nncpa; idir++){
					ncpacover[idir]=0;
				}
			} else if(evlct==nevl){
				for(int idir=0; idir<nncpa; idir++){
					ncpacover[idir]=1;
				}
				ncover++;
			}
			if(!parms->sim.idealtomo){//in idealtomo, everything is treated as NCPA
				if(!strwfs){
					info("Does not contain SURFWFS, Assume it covers all WFS.\n");
					for(int iwfs=0;iwfs<nwfs; iwfs++){
						wfscover[iwfs]=1;
					}
					ncover++;
				} else{
					readstr_numarr((void**)&wfscover, NULL, NULL, nwfs, 2, M_INT, "SURFWFS",strwfs);
					int nwfscover=0;
					for(int i=0; i<nwfs; i++){
						nwfscover+=wfscover[i]?1:0;
					}
					if(nwfscover){
						info2("Covers WFS");
						for(int i=0; i<nwfs; i++){
							if(wfscover[i]){
								info2(" %d", i);
							}
						}
						info2("\n");
					}
					if(nwfscover==nwfs){
						ncover++;
					}
				}
			}
			if(ncover==2){
				sdata.isncpa=0;
				dbg("Surface %d is CPA\n", isurf);
			}else{
				sdata.isncpa=1;
			}
			if(stropdx){
				real val=readstr_num("SURFOPDX", stropdx, NULL);
				if(isnan(val)){
					warning("SURFOPDX=%s: failed to read int\n", stropdx);
				} else{
					opdxcover=(int)val;
				}
			}
			sdata.surf=surf;
			sdata.opdxcover=opdxcover;
			sdata.isurf=isurf;
			CALL_THREAD(tdata_evl, 0);
			if(powfs){
				CALL_THREAD(tdata_wfs, 0);
			}
			if( parms->ncpa.calib){
				CALL_THREAD(tdata_ncpa, 0);
			}
			//writebin(P(powfs[0].opdadd,0), "wfs0_%d_%s", isurf,  fn);
		}
		cellfree(surfs);
	}
	free(evlcover);
	free(wfscover);
	free(ncpacover);
	free(tdata_evl);
	free(tdata_wfs);
	free(tdata_ncpa);
}

/** We trace rays from Science focal plan OPD to ploc along evaluation
	directions (type=1) or on axis only (type=2).*/
static void FitR_NCPA(dcell** xout, recon_t* recon, aper_t* aper){
	const parms_t* parms=global->parms;
	dcell* xp=NULL;
	if(aper->opdbias){
		xp=dcelldup(aper->opdbias);
	}
	applyW(xp, recon->W0, recon->W1, P(parms->ncpa.wt));
	dcellmm(xout, recon->HA_ncpa, xp, "tn", 1);
	dcellfree(xp);
}
void FitL_NCPA(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const recon_t* recon=(const recon_t*)A;
	const parms_t* parms=global->parms;
	dcell* xp=NULL;
	dcellmm(&xp, recon->HA_ncpa, xin, "nn", 1.);
	applyW(xp, recon->W0, recon->W1, P(parms->ncpa.wt));
	dcellmm(xout, recon->HA_ncpa, xp, "tn", alpha);
	dcellfree(xp);xp=NULL;
	/*dcellmm(&xp,recon->fitNW, xin, "tn", 1);
	dcellmm(xout,recon->fitNW, xp, "nn", alpha);
	dcellfree(xp);
	if(recon->actslave){
	dcellmm(xout, recon->actslave, xin, "nn", 1);
	}*/
}
static void setup_recon_HAncpa(recon_t* recon, const parms_t* parms){
	const int nevl= parms->ncpa.ndir;
	const int ndm=parms->ndm;
	dspcell* HA=recon->HA_ncpa=dspcellnew(nevl, ndm);
	TIC;tic;
	for(int ievl=0; ievl<nevl; ievl++){
		real hs=P(parms->ncpa.hs, ievl);
		for(int idm=0; idm<ndm; idm++){
			if( parms->ncpa.calib==2&&idm>0){
				continue;
			}
			const real ht=parms->dm[idm].ht;
			const real scale=1.-ht/hs;
			real displace[2];
			displace[0]=P(parms->ncpa.thetax, ievl)*ht;
			displace[1]=P(parms->ncpa.thetay, ievl)*ht;
			P(HA, ievl, idm)=mkh(P(recon->aloc, idm), recon->floc,
				displace[0], displace[1], scale, 0);
			const real theta=RSS(P(parms->ncpa.thetax, ievl), P(parms->ncpa.thetay, ievl));
			dspscale(P(HA, ievl, idm), cos(theta*parms->dm[idm].dratio));
		}
	}
	//We don't handle float ot stuck actuators in NCPA calibration.
	/*if(recon->actfloat){//avoid commanding float actuators
		act_float(recon->aloc, &recon->HA_ncpa, NULL, recon->actfloat);
	}*/
	toc2("HA_ncpa");
	if(parms->save.setup){
		writebin(recon->HA_ncpa, "HA_ncpa");
	}
}
static void lenslet_saspherical(const parms_t* parms, powfs_t* powfs, int ipowfs){
//Spherical aberration
	if(parms->powfs[ipowfs].saspherical!=0){
		if(fabs(parms->powfs[ipowfs].saspherical)<1){
			error("powfs%d: saspherical=%g should be in nm.\n",
				ipowfs, parms->powfs[ipowfs].saspherical);
		}
		real err=parms->powfs[ipowfs].saspherical*1e-9;

		/*
			The OPD is
			Phi=f(r)-a*r^2-b;
			which r^2=x^2+y^2.
			f(r) is r^4 for spherical aberration, or custom formula based on measurement (2013-05-16)

			To simulate focusing the lenslet, we optimize a and b to minimum
			\Int phi^2. For circular lenslet, the end result will be the same as Z11.

			let
			R2=\Int r^2
			R4=\Int r^4
			Rf=\Int f(r)
			Rf2=\Int f(r)*r^2.

			We have
			#b=(R4Rf-R2R6)/(R4-R2R2);
			b=(R2*Rf2-R4*Rf)/(R2*R2-R4*N)
			a=(R4-b*N)/R2;

		*/
		const int nxsa=powfs[ipowfs].pts->nxsa;
		//Can be upgraded to actual amplitude for fill factor, but need to use a full subaperture.
		dmat* ampw=dnew(nxsa, nxsa);
		for(int ix=0; ix<nxsa*nxsa; ix++){
			P(ampw, ix)=1;
		}
		dnormalize_sumabs(ampw, 1);

		real nx2=(nxsa-1)*0.5;
		real fill1d=sqrt(parms->powfs[ipowfs].safill2d);
		//normalize x,y from -1 to 1 in clear aperture
		real Rx2=pow(fill1d*nxsa/2, -2);
		dmat* opdi=dnew(nxsa, nxsa);
		real R2=0, R4=0, RfR2=0, Rf=0, RN=0;
		for(int iy=0; iy<nxsa; iy++){
			for(int ix=0; ix<nxsa; ix++){
#define MEASURED_LENSLET 0
				real rr=((iy-nx2)*(iy-nx2)+(ix-nx2)*(ix-nx2))*Rx2;
#if MEASURED_LENSLET
				real xx=(iy-nx2)*(iy-nx2)*Rx2;
				real yy=(ix-nx2)*(ix-nx2)*Rx2;
				P(opdi, ix, iy)=-50.8+49.2*yy+316.6*xx+77.0*yy*yy-309*yy*xx-260.4*xx*xx;
#else
				P(opdi, ix, iy)=rr*rr;
#endif
				real amp=P(ampw, ix, iy);
				R2+=rr*amp;
				R4+=rr*rr*amp;
				Rf+=P(opdi, ix, iy)*amp;
				RfR2+=P(opdi, ix, iy)*rr*amp;
				RN+=amp;
			}
		}
		//real b=(R2*Rf2-R4*R4)/(R2*R2-R4);
		//real a=(R4-b)/R2;
		real b=(R2*RfR2-R4*Rf)/(R2*R2-R4*RN);//fixed in 8/27/2019
		real a=(Rf-b*RN)/R2;
		real var=0;
		for(int iy=0; iy<nxsa; iy++){
			for(int ix=0; ix<nxsa; ix++){
				real rr=((iy-nx2)*(iy-nx2)+(ix-nx2)*(ix-nx2))*Rx2;
				P(opdi, ix, iy)-=a*rr+b;
				var+=P(opdi, ix, iy)*P(opdi, ix, iy)*P(ampw, ix, iy);
			}
		}
		dfree(ampw);
		dscale(opdi, err/sqrt(var));
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
				for(int i=0; i<nxsa*nxsa; i++){
					P(P(powfs[ipowfs].opdbias, jwfs), nxsa*nxsa*isa+i)+=P(opdi, i);
				}
			}
		}
		dfree(opdi);
	}
}
static void lenslet_safocuspv(const parms_t* parms, powfs_t* powfs, int ipowfs){
	//defocus specified as P/V
	if(parms->powfs[ipowfs].safocuspv!=0){
		if(fabs(parms->powfs[ipowfs].safocuspv)<1){
			error("powfs%d: safocuspv should be in nm.\n", ipowfs);
		}
		real pv=parms->powfs[ipowfs].safocuspv*1e-9;
		info2("powfs %d: Put in focus p/v value of %g to subaperture\n", ipowfs, pv*1e9);
		const int nxsa=powfs[ipowfs].pts->nxsa;
		const int npsa=nxsa*nxsa;
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			real nx2=(nxsa-1)*0.5;
			real Rx2=pow(nx2, -2);//do not use fill2d here as defocus is caused by misregistration.
			real rmax2=2*nx2*nx2*Rx2;
			real scale=pv/rmax2;
			for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
				for(int iy=0; iy<nxsa; iy++){
					for(int ix=0; ix<nxsa; ix++){
						real rr=((iy-nx2)*(iy-nx2)+(ix-nx2)*(ix-nx2))*Rx2;
						P(P(powfs[ipowfs].opdbias, jwfs), isa*npsa+ix+iy*nxsa)+=rr*scale;
					}
				}
			}
		}
	}
}
static void powfs_ncpa(const parms_t *parms, powfs_t *powfs, int ipowfs){
	if(parms->powfs[ipowfs].ncpa){
		dmat *ncpa=parms->powfs[ipowfs].ncpa;
		dbg("ncpa shape is %ldx%ld\n", NX(ncpa), NY(ncpa));
		if(NX(ncpa)%2==0&&NX(ncpa)>2&&NY(ncpa)==1){//reform nx1 factor to 2x(n/2)
			reshape(ncpa, NX(ncpa)/2, 2);
		}
		if(NX(ncpa)!=2){
			error("NCPA is in wrong format: %ldx%ld\n", NX(ncpa), NY(ncpa));
		}
		rand_t rstat;
		seed_rand(&rstat, 1);
		for(long im=0; im<NY(ncpa); im++){
			const real rms=P(ncpa, 0, im);
			const int mod=(int)P(ncpa, 1, im);
			if(mod>0){
				dmat *zer=zernike(powfs[ipowfs].loc, parms->aper.d, 0, 0, -mod);
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					real rms2=rms<0?-rms:rms*(2*randu(&rstat)-1);
					dadd(&P(powfs[ipowfs].opdbias, jwfs, 0), 1, zer, rms2);
				}
				dfree(zer);
			} else{
				error("power law to be implemented\n");
			}
		}
	}
}
/**
   Setup common and non common OPDs from surf and tsurf
 */
void setup_surf(const parms_t* parms, aper_t* aper, powfs_t* powfs, recon_t* recon){
	if(parms->load.ncpa){
		if( parms->ncpa.nsurf|| parms->ncpa.ntsurf){
			error("Please disable surf and tsurf when load.ncpa is set\n");
		}
		if(zfexist("%s/evl_opdadd.bin", parms->load.ncpa)){
			aper->opdadd=dcellread("%s/evl_opdadd.bin", parms->load.ncpa);
			if(NX(aper->opdadd)!=parms->evl.nevl||P(aper->opdadd, 0)->nx!=aper->locs->nloc){
				error("evl_opdadd is in wrong format\n");
			}
		} else{
			warning("%s/evl_opdadd.bin does not exist.\n", parms->load.ncpa);
		}
		if(zfexist("%s/evl_opdbias.bin", parms->load.ncpa)){
			aper->opdbias=dcellread("%s/evl_opdbias.bin", parms->load.ncpa);
			if(NX(aper->opdbias)!= parms->ncpa.ndir||P(aper->opdbias, 0)->nx!=recon->floc->nloc){
				error("opdbias is in wrong foramt\n");
			}
		} else{
			warning("%s/evl_opdbias.bin does not exist.\n", parms->load.ncpa);
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(zfexist("%s/powfs%d_opdadd.bin", parms->load.ncpa, ipowfs)){
				powfs[ipowfs].opdadd=dcellread("%s/powfs%d_opdadd.bin", parms->load.ncpa, ipowfs);
				powfs[ipowfs].opdbias=dcellread("%s/powfs%d_opdbias.bin", parms->load.ncpa, ipowfs);
				if(NX(powfs[ipowfs].opdadd)!=parms->powfs[ipowfs].nwfs){
					error("powfs%d_opdadd is in wrong format, expect %d, got %ld\n", ipowfs,
						parms->powfs[ipowfs].nwfs, NX(powfs[ipowfs].opdadd));
				} else{
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						if(P(powfs[ipowfs].opdadd, jwfs)->nx!=powfs[ipowfs].loc->nloc){
							error("powfs%d_opdadd cell %d is in wrong format\n", ipowfs, jwfs);
						}
					}
				}
			} else{
				warning("%s/powfs%d_opdadd.bin does not exist.\n", parms->load.ncpa, ipowfs);
			}
		}
	} else{
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->ncpa.nsurf||parms->ncpa.ntsurf
				||parms->powfs[ipowfs].safocuspv!=0||parms->powfs[ipowfs].saspherical!=0
				||parms->powfs[ipowfs].ncpa){
				if(!powfs[ipowfs].opdadd){
					powfs[ipowfs].opdadd=dcellnew_same(parms->powfs[ipowfs].nwfs, 1, powfs[ipowfs].loc->nloc, 1);
				}
				if(!powfs[ipowfs].opdbias){
					powfs[ipowfs].opdbias=dcellnew_same(parms->powfs[ipowfs].nwfs, 1, powfs[ipowfs].loc->nloc, 1);
				}
			}
			powfs_ncpa(parms, powfs, ipowfs);
			//Setup lenslet aberration
			lenslet_safocuspv(parms, powfs, ipowfs);
			lenslet_saspherical(parms, powfs, ipowfs);
		}

		if(parms->ncpa.nsurf||parms->ncpa.ntsurf){
			if(!aper->opdadd){
				aper->opdadd=dcellnew_same(parms->evl.nevl, 1, aper->locs->nloc, 1);
			}
			if(!aper->opdbias&& parms->ncpa.calib){
				aper->opdbias=dcellnew_same(parms->ncpa.ndir, 1, recon->floc->nloc, 1);
			}

			TIC;tic;
			if( parms->ncpa.ntsurf>0){
				setup_surf_tilt(parms, aper, powfs, recon);
			}
			if( parms->ncpa.nsurf>0){
				setup_surf_perp(parms, aper, powfs, recon);
			}
			toc2("surf prop");
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(powfs[ipowfs].opdbias){//add NCPA surfaces to opdadd for ray tracing. Those from dm_ncpa is not added.
				dcelladd(&powfs[ipowfs].opdadd, 1, powfs[ipowfs].opdbias, 1);
			}
		}
		/**
		   note about idealtomo: No need to pass surfaces to DM fitting
		   routine. These are already contained in dm_ncpa which got to add to the
		   DM command. Nothing needs to be done.
		*/
	}

	if( parms->ncpa.calib){//calibrate NCPA
		int any_ncpa=0;
		if(aper->opdbias){
			for(int i=0; i<parms->ncpa.ndir; i++){
				if(dmaxabs(P(aper->opdbias, i))>1e-15){
					any_ncpa=1;
					break;
				}
			}
		}

		if(any_ncpa){
			info2("calibrating NCPA\n");
			setup_recon_HAncpa(recon, parms);
			dcell* rhs=NULL;
			FitR_NCPA(&rhs, recon, aper);
			int maxit=40;
			pcg(&recon->dm_ncpa, FitL_NCPA, recon, NULL, NULL, rhs, 1, maxit);
			//don't extrapolate dm_ncpa here.
			dcellfree(rhs);
			if(parms->save.setup){
				writebin(recon->dm_ncpa, "dm_ncpa");
			}
			dspcellfree(recon->HA_ncpa);
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(any_ncpa){//apply NCPA DM command.
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					const real hs=parms->wfs[iwfs].hs;
					const real thetax=parms->wfs[iwfs].thetax;
					const real thetay=parms->wfs[iwfs].thetay;
					const real theta=RSS(thetax, thetay);

					for(int idm=0; idm<parms->ndm; idm++){
						if(!P(recon->dm_ncpa, idm)||P(recon->dm_ncpa, idm)->nx==0) continue;
						real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
						real scale=1.-ht/hs;
						real dispx=ht*thetax;
						real dispy=ht*thetay;
						real alpha=-cos(theta*parms->dm[idm].dratio);
						prop_nongrid(P(recon->aloc, idm), P(P(recon->dm_ncpa, idm)),
							powfs[ipowfs].loc, P(P(powfs[ipowfs].opdbias, jwfs)),
							alpha, dispx, dispy, scale, 0, 0);
					}
				}
			}
			if( parms->ncpa.ttr&&powfs[ipowfs].opdbias){
			/*remove average tilt from opdbias and same amount from
			  opdadd. Does not need to be very accurate.*/
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					if(P(powfs[ipowfs].opdbias, jwfs)){
						real ptt[3]={0,0,0};
						loc_calc_ptt(NULL, ptt, powfs[ipowfs].loc, 0, NULL,
							P(P(powfs[ipowfs].amp, jwfs)), P(P(powfs[ipowfs].opdbias, jwfs)));
						loc_sub_ptt(P(powfs[ipowfs].opdbias, jwfs), ptt, powfs[ipowfs].loc);
						loc_sub_ptt(P(powfs[ipowfs].opdadd, jwfs), ptt, powfs[ipowfs].loc);
					}
				}
			}
		}//for ipowfs
	}

	if(parms->save.setup||parms->save.ncpa){
		writebin(aper->opdadd, "evl_opdadd.bin");
		writebin(aper->opdbias, "evl_opdbias.bin");
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			writebin(powfs[ipowfs].opdadd, "powfs%d_opdadd.bin", ipowfs);
			writebin(powfs[ipowfs].opdbias, "powfs%d_opdbias.bin", ipowfs);
		}
	}
	if( parms->ncpa.calib&&parms->ncpa.rmsci&&recon->dm_ncpa&&aper->opdadd){
		//Don't include uncorrectable WFE in science evaluation
		dcellzero(aper->opdadd);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			for(int idm=0; idm<parms->ndm; idm++){
				const real hl=parms->dm[idm].ht;
				const real dispx=P(parms->evl.thetax, ievl)*hl;
				const real dispy=P(parms->evl.thetay, ievl)*hl;
				const real theta=RSS(P(parms->evl.thetax, ievl), P(parms->evl.thetay, ievl));
				const real scale=1-hl/P(parms->evl.hs, ievl);
				const real alpha=cos(theta*parms->dm[idm].dratio);
				prop_nongrid(P(recon->aloc, idm), P(P(recon->dm_ncpa, idm)),
					aper->locs, P(P(aper->opdadd, ievl)),
					alpha, dispx, dispy, scale, 0, 0);
			}
		}
		if(parms->save.setup>1){
			writebin(aper->opdadd, "evl_opdadd_correctable.bin");
		}
	}
	dcellfree(aper->opdbias);
	dcellfree(recon->dm_ncpa);
}
