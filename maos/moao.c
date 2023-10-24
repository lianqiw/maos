/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "recon_utils.h"
#include "moao.h"
#include "ahst.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file maos/moao.h
   Routings to setup moao and carry out moao DM fitting.
*/

/**
   Free moao_t
*/
void free_recon_moao(recon_t* recon, const parms_t* parms){
	if(!recon||!recon->moao) return;
	for(int imoao=0; imoao<parms->nmoao; imoao++){
		if(!parms->moao[imoao].used) continue;
		cellfree(recon->moao[imoao].aloc);
		cellfree(recon->moao[imoao].amap);
		cellfree(recon->moao[imoao].aimcc);
		dspcellfree(recon->moao[imoao].HA);
		dcellfree(recon->moao[imoao].NW);
		dspcellfree(recon->moao[imoao].actslave);
		dspfree(recon->moao[imoao].W0);
		dfree(recon->moao[imoao].W1);
		lcellfree(recon->moao[imoao].actstuck);
		lcellfree(recon->moao[imoao].actfloat);
	}
	free(recon->moao);
	recon->moao=NULL;
}
/**
   Prepare the propagation H matrix for MOAO and compute the reconstructor. We
   only need a reconstructor for every different MOAO type.  */
void setup_recon_moao(recon_t* recon, const parms_t* parms){
	const int nmoao=parms->nmoao;
	if(nmoao==0) return;
	if(parms->recon.alg!=0){
		error("Moao only works in recon.alg=0 mode MVR\n");
	}
	if(recon->moao){
		free_recon_moao(recon, parms);
	}
	recon->moao=mycalloc(nmoao, moao_t);
	for(int imoao=0; imoao<nmoao; imoao++){
		if(!parms->moao[imoao].used) continue;
		real dxr=parms->moao[imoao].dx;
		real dyr=dxr*parms->moao[imoao].ar;
		map_t* map=0;
		real offset=((int)round(parms->moao[imoao].order)%2)*0.5;
		real guard=parms->moao[imoao].guard*MAX(dxr, dyr);
		create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dxr, dyr, offset, guard, 0, 0, 0, parms->fit.square);
		recon->moao[imoao].aloc=loccellnew(1, 1);
		P(recon->moao[imoao].aloc,0)=map2loc(map, 0);
		P(recon->moao[imoao].aloc,0)->iac=parms->moao[imoao].iac;
		mapfree(map);
		loc_create_map_npad(P(recon->moao[imoao].aloc,0), parms->fit.square?0:1, 0, 0);
		recon->moao[imoao].amap=mapcellnew(1, 1);
		P(recon->moao[imoao].amap,0)=mapref(P(recon->moao[imoao].aloc,0)->map);
		P(recon->moao[imoao].amap,0)->iac=parms->moao[imoao].iac;
		recon->moao[imoao].aimcc=loc_mcc_ptt(P(recon->moao[imoao].aloc,0), NULL);
		dinvspd_inplace(recon->moao[imoao].aimcc);
		recon->moao[imoao].HA=dspcellnew(1, 1);
		P(recon->moao[imoao].HA,0)=mkh(P(recon->moao[imoao].aloc,0), recon->floc, 0, 0, 1);
		if(parms->moao[imoao].actstuck){
			recon->moao[imoao].actstuck=lcellnew(1, 1);
			P(recon->moao[imoao].actstuck,0)=loc_coord2ind(P(recon->moao[imoao].aloc,0), parms->moao[imoao].actstuck);
			if(parms->dbg.recon_stuck){
				act_stuck(recon->moao[imoao].aloc, recon->moao[imoao].HA, recon->moao[imoao].actstuck);
			}
		}
		if(parms->moao[imoao].actfloat){
			recon->moao[imoao].actfloat=lcellnew(1, 1);
			P(recon->moao[imoao].actfloat,0)=loc_coord2ind(P(recon->moao[imoao].aloc,0), parms->moao[imoao].actfloat);
			act_float(recon->moao[imoao].aloc, &recon->moao[imoao].HA, NULL, recon->moao[imoao].actfloat);
		}

		if(parms->moao[imoao].lrt_ptt){
			recon->moao[imoao].NW=dcellnew(1, 1);
			long nloc=P(recon->moao[imoao].aloc,0)->nloc;
			P(recon->moao[imoao].NW,0)=dnew(nloc, 3);
			dmat* pNW=P(recon->moao[imoao].NW,0)/*PDMAT*/;
			real scl=1./nloc;
			real scl2=scl*2./parms->aper.d;
			const real* locx=P(recon->moao[imoao].aloc,0)->locx;
			const real* locy=P(recon->moao[imoao].aloc,0)->locy;
			for(long iloc=0; iloc<nloc; iloc++){
			/*We don't want piston/tip/tilt on the mems. */
				P(pNW, iloc, 0)=scl;/*piston; */
				P(pNW, iloc, 1)=scl2*locx[iloc];/*tip */
				P(pNW, iloc, 2)=scl2*locy[iloc];/*tilt */
			}
		}
		recon->moao[imoao].W0=dspref(recon->W0);
		recon->moao[imoao].W1=dref(recon->W1);
		if(parms->moao[imoao].actslave){
			recon->moao[imoao].actcpl=genactcpl(recon->moao[imoao].HA, recon->moao[imoao].W1);
			recon->moao[imoao].actslave=slaving(recon->moao[imoao].aloc,
				recon->moao[imoao].actcpl,
				parms->dbg.recon_stuck?recon->moao[imoao].actstuck:0,
				recon->moao[imoao].actfloat,
				0.1, 1./recon->floc->nloc, 1);
		}
		if(parms->save.setup){
			writebin(recon->moao[imoao].aloc, "moao%d_aloc", imoao);
			writebin(recon->moao[imoao].HA, "moao%d_HA", imoao);
			writebin(recon->moao[imoao].NW, "moao%d_NW", imoao);
			writebin(recon->moao[imoao].actslave, "moao%d_actslave", imoao);
		}
		if(parms->plot.setup){
			plotloc("FoV", parms, P(recon->moao[imoao].aloc,0), 0, "moao_aloc");
		}
	}/*imoao */
}

/**
   Apply fit right hand side matrix in CG mode without using assembled matrix
   for MOAO. subtract contributions from DMs that are in common path. Be careful
   which time step the dmcommon is. The DM common should use the commands on the
   step that you are going to apply the MOAO command for. That is the integrator
   output after this computation.  */

static void
moao_FitR(dcell** xout, const recon_t* recon, const parms_t* parms, int imoao,
	real thetax, real thetay, real hs,
	const dcell* opdr, const dcell* dmcommon, dcell** rhsout, const real alpha){

	dcell* xp=dcellnew(1, 1);
	P(xp,0)=dnew(recon->floc->nloc, 1);

	for(int ipsr=0; ipsr<recon->npsr; ipsr++){
		const real ht=P(parms->atmr.ht,ipsr);
		real scale=1.-ht/hs;
		if(parms->tomo.square){
			map_t map;
			memcpy(&map, P(recon->xmap,ipsr), sizeof(map_t));
			map.p=P(P(opdr,ipsr));
			prop_grid_stat(&map, recon->floc->stat,
				P(P(xp,0)), 1,
				thetax*ht, thetay*ht, scale, 0, 0, 0);
		} else{
			prop_nongrid(P(recon->xloc,ipsr), P(P(opdr,ipsr)),
				recon->floc, P(P(xp,0)), 1,
				thetax*ht, thetay*ht, scale, 0, 0);
		}
	}
	//static int count=-1; count++;
	//writebin(P(xp,0), "opdfit0_%d", count);
	for(int idm=0; idm<parms->ndm; idm++){
		const real ht=parms->dm[idm].ht;
		real scale=1.-ht/hs;
		prop_nongrid(P(recon->aloc,idm), P(P(dmcommon,idm)),
			recon->floc, P(P(xp,0)), -1,
			thetax*ht, thetay*ht, scale,
			0, 0);
	}
	//writebin(P(xp,0), "opdfit1_%d", count);
	if(rhsout){
		*rhsout=dcelldup(xp);
	}
	real wt=1;
	applyW(xp, recon->moao[imoao].W0, recon->moao[imoao].W1, &wt);
	dcellmm(xout, recon->moao[imoao].HA, xp, "tn", alpha);
	dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode
   without using assembled matrix. Slow. don't
   use. Assembled matridx is faster because of multiple
   directions.
*/

static void
moao_FitL(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const moao_t* moao=(const moao_t*)A;
	dcell* xp=NULL;
	real wt=1;
	dcellmm(&xp, moao->HA, xin, "nn", 1.);
	applyW(xp, moao->W0, moao->W1, &wt);
	dcellmm(xout, moao->HA, xp, "tn", alpha);
	dcellfree(xp);xp=NULL;
	dcellmm(&xp, moao->NW, xin, "tn", 1);
	dcellmm(xout, moao->NW, xp, "nn", alpha);
	dcellfree(xp);
	dcellmm(xout, moao->actslave, xin, "nn", alpha);
}
/**
   moao_recon happens after the common DM fitting and its integrator output
   to take into account the delay in DM commands. there is no close loop
   filtering in MOAO DM commands, but there is still a time delay of 2
   cycles.
*/

void moao_recon(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int nwfs=parms->nwfs;
	const int nevl=parms->evl.nevl;
	dcell* dmcommon=NULL;
	if(1){/*Take High order fitting result */
		dcellcp(&dmcommon, simu->dmrecon);
	} else{/*Take integrator output, remove NGS modes if any. */
		if(parms->sim.closeloop){
			if(parms->sim.fuseint){
				dcellcp(&dmcommon, P(simu->dmint->mintc,0));
				if(parms->recon.split==1){
					ngsmod_remove(simu, dmcommon);
				}
			} else{
				dcellcp(&dmcommon, P(simu->dmint->mintc,0));
			}
		} else{
			dcellcp(&dmcommon, simu->dmerr);
		}
	}
	dcell* rhs=NULL;
	int iy=parms->sim.closeloop?1:0;
	if(simu->dm_wfs){/*There is MOAO DM for WFS */
		dcell* dmmoao=dcellnew(1, 1);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			int imoao=parms->powfs[ipowfs].moao;
			dcell* rhsout=NULL;
			if(imoao<0) continue;
			real hs=parms->wfs[iwfs].hs;
			P(dmmoao,0)=P(simu->dm_wfs,iwfs,iy);
			dcellzero(rhs);
			moao_FitR(&rhs, recon, parms, imoao,
				parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay,
				hs, simu->opdr, dmcommon, parms->plot.run?&rhsout:NULL, 1);
			pcg(&dmmoao, moao_FitL, &recon->moao[imoao], NULL, NULL, rhs,
				parms->fit.cgwarm, parms->fit.maxit);
				/*if(parms->recon.split){//remove the tip/tilt form MEMS DM
				  real ptt[3]={0,0,0};
				  loc_t *aloc=recon->moao[imoao].aloc;
				  dmat *aimcc=recon->moao[imoao].aimcc;
				  loc_calc_ptt(NULL, ptt, aloc, 0, aimcc, NULL, P(P(dmmoao,0)));
				  loc_sub_ptt(P(P(dmmoao,0)), ptt, aloc);
				  }*/
			if(!isinf(parms->moao[imoao].stroke)){
				int nclip=dclip(P(dmmoao,0),
					-parms->moao[imoao].stroke,
					parms->moao[imoao].stroke);
				if(nclip>0){
					dbg("wfs %d: %d actuators clipped\n", iwfs, nclip);
				}
			}
			if(parms->plot.run){
				drawopd("MOAO WFS RHS", recon->floc, P(rhsout,0), P(parms->plot.opdmax),
					"MOAO for WFS", "x (m)", "y(m)", "Wfs rhs %d", iwfs);
				drawopd("MOAO WFS", P(recon->moao[imoao].aloc,0), P(dmmoao,0), P(parms->plot.opdmax),
					"MOAO for WFS", "x (m)", "y(m)", "Wfs %d", iwfs);
			}
			if(parms->save.dm){
				zfarr_push(simu->save->dm_wfs[iwfs], simu->wfsisim, P(dmmoao,0));
			}
			dcellfree(rhsout);
			P(dmmoao,0)=NULL;
		}/*if wfs */
		dcellfree(dmmoao);
	}
	if(simu->dm_evl){/*There is MOAO DM for Science */
		int imoao=parms->evl.moao;
		dcell* dmmoao=dcellnew(1, 1);
		for(int ievl=0; ievl<nevl; ievl++){
			P(dmmoao,0)=P(simu->dm_evl,ievl,iy);
			dcell* rhsout=NULL;
			dcellzero(rhs);
			moao_FitR(&rhs, recon, parms, imoao,
				P(parms->evl.thetax,ievl), P(parms->evl.thetay,ievl),
				INFINITY, simu->opdr, dmcommon, (parms->plot.run||1)?&rhsout:NULL, 1);

			pcg(&dmmoao, moao_FitL, &recon->moao[imoao], NULL, NULL, rhs,
				parms->fit.cgwarm, parms->fit.maxit);
			if(0){
				writebin(P(rhsout,0), "evl_rhs_%d_%d", ievl, simu->perfisim);
				writebin(P(dmmoao,0), "evl_dmrecon_%d_%d", ievl, simu->perfisim);
			}
			/*if(parms->recon.split){//remove the tip/tilt form MEMS DM
			  real ptt[3]={0,0,0};
			  loc_t *aloc=recon->moao[imoao].aloc;
			  dmat *aimcc=recon->moao[imoao].aimcc;
			  loc_calc_ptt(NULL, ptt, aloc, 0, aimcc, NULL, P(P(dmmoao,0)));
			  loc_sub_ptt(P(P(dmmoao,0)), ptt, aloc);
			  }*/
			if(!isinf(parms->moao[imoao].stroke)){
				int nclip=dclip(P(dmmoao,0),
					-parms->moao[imoao].stroke,
					parms->moao[imoao].stroke);
				if(nclip>0){
					dbg("evl %d: %d actuators clipped\n", ievl, nclip);
				}
			}
			if(parms->plot.run){
				drawopd("MOAO EVL RHS", recon->floc, P(rhsout,0), P(parms->plot.opdmax),
					"MOAO for WFS", "x (m)", "y(m)", "Evl %d", ievl);
				drawopd("MOAO EVL", P(recon->moao[imoao].aloc,0), P(dmmoao,0), P(parms->plot.opdmax),
					"MOAO for EVL", "x (m)", "y(m)", "Evl %d", ievl);
			}
			if(parms->save.dm){
				zfarr_push(simu->save->dm_evl[ievl], simu->perfisim, P(dmmoao,0));
			}
			dcellfree(rhsout);
			P(dmmoao,0)=NULL;
		}/*ievl */
		dcellfree(dmmoao);
	}
	dcellfree(dmcommon);
	dcellfree(rhs);
}
