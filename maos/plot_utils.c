/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "plot_utils.h"
#include "ahst.h"
#include "utils.h"
/*
   A few utility routines
*/
/**
 * Return type of wfs
*/
static const char *
powfs_legend(const parms_t *parms, int ipowfs){
	const char *const legwfs[]={
	"LGS WFS",
	"NGS WFS",
	"PWFS",
	"TTF WFS",
	"TT WFS",
	"Other WFS",
	};
	int ilegwfs=6;
	if(parms->powfs[ipowfs].lo){
		if(parms->powfs[ipowfs].order==1){
			ilegwfs=4;
		} else{
			ilegwfs=3;
		}
	} else{
		if(parms->powfs[ipowfs].trs){
			ilegwfs=0;
		} else{
			if(parms->powfs[ipowfs].type==WFS_PY){
				ilegwfs=2;
			} else{
				ilegwfs=1;
			}
		}
	}
	return legwfs[ilegwfs];
}
/**
   Plot the loc, together with all beams
*/
void plot_loc(const char* fig, const parms_t* parms, int show_recon, int show_evl, 
	loc_t* loc, real ht, const char* format, ...){
	format2fn;
	int ncir=parms->evl.nevl+parms->fit.nfit+parms->nwfs;
	if(parms->ncpa.calib){
		ncir+= parms->ncpa.ndir;
	}
	dmat* cir=dnew(4, ncir);
	const char *legend[ncir+1];
	memset(legend, 0, sizeof(char*)*(ncir+1));
	int count=0;
	for(int ievl=0; show_evl && ievl<parms->evl.nevl; ievl++){
		real hs=P(parms->evl.hs,ievl);
		P(cir, 0, count)=ht*P(parms->evl.thetax,ievl);
		P(cir, 1, count)=ht*P(parms->evl.thetay,ievl);
		P(cir, 2, count)=parms->aper.d*0.5*(1-ht/hs);
		P(cir, 3, count)=0xFF0000;/*rgb color */
		count++;
		if(ievl==0) legend[count]="Evaluation";//after count++ because points are plotted first
	}
	for(int ifit=0; show_recon && ifit<parms->fit.nfit; ifit++){
		real hs=P(parms->fit.hs,ifit);
		P(cir, 0, count)=ht*P(parms->fit.thetax,ifit);
		P(cir, 1, count)=ht*P(parms->fit.thetay,ifit);
		P(cir, 2, count)=parms->aper.d*0.5*(1-ht/hs);
		P(cir, 3, count)=0xFF22DD;/*rgb color */
		count++;
		if(ifit==0) legend[count]="DM Fitting";
	}
	for(int idir=0; show_recon && idir< parms->ncpa.ndir; idir++){
		real hs=P(parms->ncpa.hs,idir);
		P(cir, 0, count)=ht*P(parms->ncpa.thetax,idir);
		P(cir, 1, count)=ht*P(parms->ncpa.thetay,idir);
		P(cir, 2, count)=parms->aper.d*0.5*(1-ht/hs);
		P(cir, 3, count)=0x22FF00;/*rgb color */
		count++;
		if(idir==0) legend[count]="NCPA";
	}

	for(int iwfs=0; show_recon && iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		real hs=parms->wfs[iwfs].hs;
		P(cir, 0, count)=parms->wfs[iwfs].thetax*ht;
		P(cir, 1, count)=parms->wfs[iwfs].thetay*ht;
		P(cir, 2, count)=parms->aper.d*0.5*(1.-ht/hs);
		if(!isinf(hs)){//LGS
			P(cir, 3, count)=0xFF8800;
		} else if(!parms->powfs[ipowfs].lo){//Hi NGS
			P(cir, 3, count)=0xFFFF00;
		} else if(parms->powfs[ipowfs].order>1){//TTF
			P(cir, 3, count)=0x0000FF;//TTF
		} else{
			P(cir, 3, count)=0x0000FF;//TT
		}
		count++;
		if(P(parms->powfs[ipowfs].wfsind,iwfs)==0){
			legend[count]=powfs_legend(parms, ipowfs);
		}
	}
	if(count>ncir){
		error("Overflow\n");
	}else if(count<ncir){
		dresize(cir, NX(cir), count);
	}
	draw(fig, (plot_opts){.ngroup=1, .loc=&loc, .cir=cir, .legend=legend},
		"Coordinate", "x (m)", "y (m)", "%s", fn);
	dfree(cir);
}
/**
   ploted all the different beam directions as points. 
   style is defined as follows:
   bits 1-3:the point style.
   bit  4: whether points are connected.
   bits 5-8: size
   bits 9-32:color
 
   The color follows RGB representation: (Since 2011-02-18)
   bits 32-25: Red. bits 24-17: Green. bits 16-9: Blue.
   */
void plot_dir(const char* fig, const parms_t* parms, real totfov, const char* format, ...){
	format2fn;
	int ncir=1;
	dmat* cir=dnew(4, ncir);
	P(cir, 0, 0)=0;
	P(cir, 1, 0)=0;
	P(cir, 2, 0)=totfov/2;
	P(cir, 3, 0)=0x000000;/*rgb color */
	int ngroup=3+parms->npowfs+parms->nlgswfs;
	const char* legend[ngroup];
	loccell* locs=(loccell*)cellnew(ngroup, 1);
	int32_t* style=mycalloc(ngroup, int32_t);
	int count=0;
	legend[count]="Evaluation";
	style[count]=(0xFF0000<<8)|(4<<4)|3;
	P(locs,count)=locnew(parms->evl.nevl, 0, 0);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		P(locs,count)->locx[ievl]=P(parms->evl.thetax,ievl)*RAD2AS;
		P(locs,count)->locy[ievl]=P(parms->evl.thetay,ievl)*RAD2AS;
	}
	count++;
	legend[count]="DM Fitting";
	style[count]=(0xFF22DD<<8)|(4<<4)|3;
	P(locs,count)=locnew(parms->fit.nfit, 0, 0);
	for(int ifit=0; ifit<parms->fit.nfit; ifit++){
		P(locs,count)->locx[ifit]=P(parms->fit.thetax,ifit)*RAD2AS;
		P(locs,count)->locy[ifit]=P(parms->fit.thetay,ifit)*RAD2AS;
	}
	count++;
	legend[count]="NCPA";
	style[count]=(0x22FF00<<8)|(4<<4)|2;
	P(locs,count)=locnew( parms->ncpa.ndir, 0, 0);
	for(int ifit=0; ifit< parms->ncpa.ndir; ifit++){
		P(locs,count)->locx[ifit]=P(parms->ncpa.thetax,ifit)*RAD2AS;
		P(locs,count)->locy[ifit]=P(parms->ncpa.thetay,ifit)*RAD2AS;
	}
	count++;

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		legend[count]=powfs_legend(parms, ipowfs);
		P(locs,count)=locnew(parms->powfs[ipowfs].nwfs, 0, 0);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
			P(locs,count)->locx[jwfs]=parms->wfs[iwfs].thetax*RAD2AS;
			P(locs,count)->locy[jwfs]=parms->wfs[iwfs].thetay*RAD2AS;
		}
		if(!isinf(parms->powfs[ipowfs].hs)){
			style[count]=(0xFF8800<<8)|(4<<4)|2;/*LGS */
		} else if(!parms->powfs[ipowfs].lo){
			style[count]=(0xFFFF00<<8)|(4<<4)|1;/*Hi NGS*/
		} else if(parms->powfs[ipowfs].order>1){
			style[count]=(0x0000FF<<8)|(4<<4)|4;/*TTF*/
		} else{
			style[count]=(0x0000FF<<8)|(4<<4)|1;/*TT */
		}
		count++;
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
		if(parms->powfs[ipowfs].llt){
			if(jwfs==0) legend[count]="LLT"; else legend[count]=NULL;
			style[count]=(0xFF8800<<8)|(4<<4)|3|8;//connect LLT to LGS
			P(locs,count)=locnew(2, 0, 0);
			const real hs=parms->powfs[ipowfs].hs;
			P(locs,count)->locx[0]=PR(parms->powfs[ipowfs].llt->ox, jwfs)/hs*RAD2AS;
			P(locs,count)->locy[0]=PR(parms->powfs[ipowfs].llt->oy, jwfs)/hs*RAD2AS;
			P(locs,count)->locx[1]=parms->wfs[iwfs].thetax*RAD2AS;
			P(locs,count)->locy[1]=parms->wfs[iwfs].thetay*RAD2AS;
			count++;
		}
	}
	if(count!=ngroup){
		error("count=%d, ngroup=%d. they should equal.\n", count, ngroup);
	}
	real limit[4];
	limit[0]=limit[2]=-totfov/2;
	limit[1]=limit[3]=totfov/2;
	draw(fig, (plot_opts){.ngroup=ngroup, .loc=P(locs), .style=style, .limit=limit, .cir=cir, .legend=legend},
		"Asterism", "x (arcsec)", "y (arcsec)", "%s", fn);
	dfree(cir);
	cellfree(locs);
	free(style);
}
/**
   Plot grid points, amplitude maps and NCPA.
 */
void plot_setup(const parms_t* parms, const powfs_t* powfs,
	const aper_t* aper, const recon_t* recon){
	int draw_single_save=draw_single;
	draw_single=0;
	plot_dir("Aperture", parms, parms->sim.fov*RAD2AS, "fov");/*plot wfs/evaluation direction */
	if(recon){
		if(parms->plot.setup>1){
			plot_loc("Aperture", parms, 1, 0, recon->ploc, 0, "ploc");
			plot_loc("Aperture", parms, 1, 0, recon->floc, 0, "floc");
		}
		for(int idm=0; idm<parms->ndm; idm++){
			real ht=parms->dm[idm].ht;
			plot_loc("Aperture", parms, 1, 1, P(recon->aloc,idm), ht, "aloc %d", idm);
			/*if(recon->actcpl){
				drawopd("Aperture", P(recon->aloc, idm), P(recon->actcpl, idm), 0, 
					"DM Actuator Coupling Factor", "x (m)", "y (m)", "actcpl %d", idm);
			}*/
		}
		if(parms->plot.setup){
			for(int ips=0; ips<recon->npsr; ips++){
				const real ht=P(recon->ht,ips);
				plot_loc("Aperture", parms, 1, 0, P(recon->xloc,ips), ht, "xloc %02d", ips);
			}
		}
	}
	drawopd("Aperture", aper->locs, aper->amp1, 0, "Performance Evaluation Amplitude Map",
		"x (m)", "y (m)", "aper");

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		for(int iamp=0; iamp<PN(powfs[ipowfs].amp); iamp++){
			int iwfs=P(parms->powfs[ipowfs].wfs,iamp);
			if(parms->plot.setup>1){
				drawopd("Aperture", powfs[ipowfs].loc, P(powfs[ipowfs].amp,iamp), 0,
					"WFS Amplitude Map", "x (m)", "y (m)", "amp wfs%d", iwfs);
			}
			if(powfs[ipowfs].saloc->nloc>4){
				drawopd("Aperture", powfs[ipowfs].saloc, P(powfs[ipowfs].saa,iamp), 0,
					"WFS Subaperture Amplitude", "x (m)", "y (m)", "saa wfs%d ", iwfs);
			}
		}
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
			if(powfs[ipowfs].gradoff){
				drawgrad("Goff", powfs[ipowfs].saloc, PR(powfs[ipowfs].saa, jwfs), P(powfs[ipowfs].gradoff, jwfs),
					parms->plot.grad2opd, parms->powfs[ipowfs].trs, 0,
					"WFS Offset", "x (m)", "y (m)", "Gncpa %d", iwfs);
			}
			if(powfs[ipowfs].intstat&&powfs[ipowfs].intstat->cogmask){
				drawints("Gmask", powfs[ipowfs].saloc, powfs[ipowfs].intstat->cogmask, 0,
						"WFS CoG Mask", "x", "y", "WFS %2d", iwfs);
			}
		}
	}
	/*if(recon->floc){
		drawopd("Aperture", recon->floc, recon->W1, 0, "DM Fitting Amplitude Map", "x (m)", "y (m)", "W1");
	}*/
	draw_single=draw_single_save;
}

void plot_dm(const parms_t *parms, const recon_t *recon, const dcell *ac, int modal, const char *title, const char *type){
	if(!ac) return;
	for(int idm=0; idm<NX(ac); idm++){
		if(!draw_current_format("DM", "%s %d", type, idm)) continue;
		dmat *dmc=NULL;
		if(recon->amod && modal){
			dmm(&dmc, 0, P(recon->amod, idm, idm), P(ac,idm), "nn", 1);
		}else{
			dmc=dref(P(ac,idm));
		}
		drawopd("DM", P(recon->aloc, idm), dmc, parms->plot.opdmax, title, "x (m)", "y (m)", "%s %d", type, idm);
		dfree(dmc);
	}
}
void plot_dm_lo(sim_t *simu, dcell *merr, const char *title, const char *type){
	if(!simu||!merr) return;
	int added=0;
	for(int idm=0; idm<NX(simu->dmtmp); idm++){
		if(draw_current_format("DM", "%s %d", type, idm)){
			if(!added){
				added=1;//add only once
				dcellzero(simu->dmtmp);
				addlow2dm(&simu->dmtmp, simu->parms, simu->recon, merr, 1);
				break;
			}
		}
	}
	if(added){
		plot_dm(simu->parms, simu->recon, simu->dmtmp, 1, title, type);
	}
}
void plot_dmreal(sim_t *simu){
    const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
    if(parms->plot.run&&simu->reconisim>=0 && simu->reconisim%parms->plot.run==0){
		if(parms->sim.closeloop){
			plot_dm(parms, recon, P(simu->dmint->mintc, 0), 1, "Deformable Mirror Integrator (Hi)", "Int");
			if(simu->Mint_lo&&!parms->sim.fuseint){
				plot_dm_lo(simu, P(simu->Mint_lo->mintc,0), "Deformable Mirror Integrator (Lo)", "Int Lo");
			}
		}
		if(simu->dmreal){
			//plot_dm(parms, recon, simu->dmcmd, "Deformable Mirror Command", "Cmd");
			plot_dm(parms, recon, simu->dmreal, 0, "Deformable Mirror Command", "Real");

			if(simu->ttmreal&&draw_current("DM", "Real TTM")){
				int idm=0;
				real ptt[3]={0,0,0};
				ptt[1]=P(simu->ttmreal,0);
				ptt[2]=P(simu->ttmreal,1);
				dmat* tmp=dnew(P(recon->aloc,idm)->nloc, 1);
				loc_add_ptt(tmp, ptt, P(recon->aloc,idm));
				drawopd("DM", P(recon->aloc,idm), tmp, parms->plot.opdmax,
					"TTM Command", "x (m)", "y (m)", "Real TTM");
				dfree(tmp);
			}

			if(draw_current("DM", "Real OPD OA")){
				dmat* opd=dnew(simu->aper->locs->nloc, 1);
				for(int idm=0; idm<parms->ndm; idm++){
					int ind=parms->evl.nevl*idm;
					simu->evl_propdata_dm[ind].phiout=P(opd);
					CALL_THREAD(simu->evl_prop_dm[ind], 0);
				}
				dscale(opd, -1);
				if(simu->ttmreal){
					real ptt[]={0,0,0};
					ptt[1]=P(simu->ttmreal,0);
					ptt[2]=P(simu->ttmreal,1);
					loc_add_ptt(opd, ptt, simu->aper->locs);
				}

				drawopd("DM", simu->aper->locs, opd, parms->plot.opdmax,
					"Deformable Mirror OPD On Axis", "x (m)", "y (m)", "Real OPD OA");
				dfree(opd);
			}
		}
	}
}
void plot_recon(sim_t *simu){
    const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(simu->reconisim<0) return;
	if(parms->plot.run&&simu->reconisim%parms->plot.run==0){
		if(simu->dm_wfs){
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				int imoao=parms->powfs[ipowfs].moao;
				if(imoao<0) continue;
				drawopd("DM", P(recon->moao[imoao].aloc, 0), P(simu->dm_wfs, iwfs), parms->plot.opdmax,
					"MOAO DM Command", "x(m)", "y(m)", "WFS %2d", iwfs);
			}
		}
		if(simu->dm_evl){
			int imoao=parms->evl.moao;
			for(int ievl=0; ievl<parms->evl.nevl&&imoao>=0; ievl++){
				drawopd("DM", P(recon->moao[imoao].aloc,0), P(simu->dm_evl,ievl), parms->plot.opdmax,
					"MOAO DM Command", "x(m)", "y(m)", "Evl %d", ievl);
			}
		}
		if(parms->recon.alg==RECON_MVR&&simu->dmrecon){
			plot_dm(parms, recon, simu->dmrecon, 1, "Deformable Mirror Fitting Output", "Fit");
		}
		plot_dm(parms, recon, simu->dmerr, 1, "Deformable Mirror Error Signal (Hi)", "Err Hi");

		if(parms->recon.alg==RECON_MVR&&simu->opdr){
			for(int i=0; i<NX(simu->opdr); i++){
				if(P(simu->opdr,i)){
					drawopd("Opdr", P(recon->xloc,i), P(simu->opdr,i), parms->plot.opdmax,
						"Reconstructed Atmosphere", "x (m)", "y (m)", "Layer %02d", i);
				}
			}
		}
		plot_dm_lo(simu, simu->Merr_lo, "Deformable Mirror Error Signal (Lo)", "Err Lo");
	}
}
void plot_gradoff(sim_t *simu, int iwfs){
	const parms_t *parms=simu->parms;
	if(parms->plot.run&&simu->wfsisim%parms->plot.run==0){
		if(iwfs<0){
			for(iwfs=0; iwfs<parms->nwfs; iwfs++){
				plot_gradoff(simu, iwfs);
			}
		}else{
			int ipowfs=parms->wfs[iwfs].powfs;
			int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
			int draw_single_save=draw_single;
			draw_single=0;
			drawgrad("Goff", simu->powfs[ipowfs].saloc, PR(simu->powfs[ipowfs].saa, jwfs), P(simu->gradoff, iwfs),
				parms->plot.grad2opd, parms->powfs[ipowfs].trs, parms->plot.gmax,
				"WFS Offset", "x (m)", "y (m)", "WFS %2d", iwfs);
			draw_single=draw_single_save;
		}
	}
}
void plot_psf(const_anycell _psf2s, const char* psfname, int type, int ievl, dmat* wvl, int zlog, real psfmin){
	dmat* psftemp=NULL;
	dmat* psfreal=NULL;
	const char* title, * tab;
	for(int iwvl=0; iwvl<NX(_psf2s.c); iwvl++){
		switch(type){
			case 2:
				title="Science Diffraction Limited PSF";
				tab="DL";
				break;
			case 1:
				title="Science Closed Loop PSF";
				tab="CL";
				break;
			case 0:
				title="Science Open Loop PSF";
				tab="OL";
				break;
			default:
				title="PSF";
				tab="PSF";
		}
		char tabname[64];
		snprintf(tabname, sizeof(tabname), "%s%2d %.2f", tab, ievl, P(wvl,iwvl)*1e6);
		if(draw_current(psfname, tabname)){
			cell* c=P(_psf2s.c, iwvl);
			if(dmat_cast(c)){
				psfreal=dmat_cast(c);
			}else if(cmat_cast(c)){
				if(psftemp&&NX(psftemp)!=c->nx){
					dfree(psftemp);
				}
				cabs22d(&psftemp, 0, cmat_cast(c), 1);
				psfreal=psftemp;
			}
			draw(psfname, (plot_opts){.image=psfreal,.zlim={psfmin,1},.zlog=zlog},title, "x", "y", "%s", tabname);
		}
	}
	dfree(psftemp);
}
