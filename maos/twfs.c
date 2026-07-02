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
#include "twfs.h"
#include "save.h"
#include "dither_utils.h"
#include "plot_utils.h"

/**
   Mapping Truth WFS reconstructed modes (in actuator space) to gradients.
 */
void twfs_setup_GR(recon_t* recon, const parms_t* parms){
	if(parms->itpowfs==-1&&!(parms->ilgspowfs!=-1&&parms->powfs[parms->ilgspowfs].dither==-1&&parms->powfs[parms->ilgspowfs].phytype_sim2==PTYPE_COG)){
		return;
	}
	int nlayer=1;
	//Need two layers when there are multiple TWFS or LGS WFS gradient offset needs projection adjustment
	if((parms->itpowfs!=-1 && parms->powfs[parms->itpowfs].nwfs>1)
		||(parms->ilgspowfs!=-1&&parms->powfs[parms->ilgspowfs].dither==-1&&parms->powfs[parms->ilgspowfs].phytype_sim2==PTYPE_COG&&parms->powfs[parms->ilgspowfs].nwfs>1)){
		nlayer=parms->ndm;
	}
	recon->GRall=dcellnew(parms->nwfs, nlayer);

	const int rmax=parms->recon.twfs_rmax?parms->recon.twfs_rmax:(parms->powfs[parms->itpowfs].order/2);
	const int rmin=parms->recon.twfs_rmin?parms->recon.twfs_rmin:3;
	const int zradonly=parms->recon.twfs_radonly;
	recon->Rmod=dcellnew(parms->ndm, 1);
	info("Truth wfs controls mode %s from order %d to %d on %d layers.\n", zradonly?"radial":"all modes", rmin, rmax, nlayer);
	for(int idm=0; idm<nlayer; idm++){
		const loc_t* loc=P(recon->aloc, idm);
		int rmin2=rmin;
		if(idm>0 && rmin2<3){//don't place those on upper layer.
			rmin2=3;
		}
		dmat* opd=zernike(loc, 0, rmin2, rmax, zradonly);
		if(parms->recon.modal){//conver to modal actuator space
			dmat *opd2=NULL;
			dcellmm(&opd2, P(recon->amodpinv, idm, idm), opd, "nn", 1);
			dfree(opd); opd=opd2;
		}
		OMP_FOR(8)
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].skip==2||(!parms->recon.twfs_offsetdm && !parms->powfs[ipowfs].lo)){
				if(P(recon->GA, iwfs, idm)){
					dcellmm(&P(recon->GRall, iwfs, idm), P(recon->GA, iwfs, idm), opd, "nn", 1);
				}else{
					error("Please implement without GA\n");
				}
			}
		}
		if(parms->recon.twfs_offsetdm){
			P(recon->Rmod, idm)=opd; opd=NULL;
		}else{
			dfree(opd);
		}
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].skip==2){//twfs
			int nlayer2=MIN(parms->powfs[ipowfs].nwfs, nlayer);
			if(parms->powfs[ipowfs].nwfs>1&&nlayer==1 && parms->ndm>1){
				warning("recon.GRwfs should have more than 1 layer when there are multiple twfs and DMs.\n");
			}
			if(rmin<3){
				warning("rmin should be 3 for truth wfs.\n");
			}
			if(!recon->GRtwfs){
				recon->GRtwfs=dcellnew(parms->nwfs, nlayer2);
			}
			for(int ilayer=0; ilayer<nlayer2; ilayer++){
				P(recon->GRtwfs, iwfs, ilayer)=dref(P(recon->GRall, iwfs, ilayer));
			}
		}
		if(!parms->recon.twfs_offsetdm){
			if(parms->powfs[ipowfs].llt&&parms->powfs[ipowfs].dither==-1&&parms->powfs[ipowfs].phytype_sim2==PTYPE_COG){
				int nlayer2=MIN(parms->powfs[ipowfs].nwfs, nlayer);
				if(parms->powfs[ipowfs].nwfs>1&&nlayer==1){
					error("recon.GRwfs should have more than 1 layer for sodium fitting projection.\n");
				}
				if(rmin>2){
					error("rmin should be 1 or 2 for sodium fitting projection\n");
				}
				if(!recon->GRlgs){
					recon->GRlgs=dcellnew(parms->powfs[ipowfs].nwfs, nlayer2);
				}
				int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
				for(int ilayer=0; ilayer<nlayer2; ilayer++){
					P(recon->GRlgs, jwfs, ilayer)=dref(P(recon->GRall, iwfs, ilayer));
				}
			}
		}
	}
	if(recon->GRlgs){
		//2021-10-15: Since we are not selecting modes, there is no need for high threshold
		//to high threshold makes the filtering ill formed
		recon->RRlgs=dcellpinv(recon->GRlgs, NULL);
	}
	if(parms->save.recon){
		writebin(recon->GRall, "twfs_GR");
		if(recon->RRlgs) writebin(recon->RRlgs, "twfs_RRlgs");
		if(parms->recon.twfs_offsetdm) writebin(recon->Rmod, "twfs_Rmod");
	}
}

/**
   Setup reconstructor for TWFS or sodium fit gradient
*/
void twfs_setup_RR(recon_t* recon, const parms_t* parms){
	if(parms->itpowfs==-1){
		return;
	}

	dspcell* neai=NULL;
	//int itwfs=-1;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(parms->powfs[ipowfs].skip==2){//twfs
			if(!neai){
				neai=dspcellnew(parms->nwfs, parms->nwfs);
			}
			P(neai,iwfs,iwfs)=dspref(P(recon->saneai,iwfs,iwfs));
		} 
	}
	//need to set a threshold to avoid other modes reconstruct to spherical modes.
	
	if(recon->GRtwfs){
		cellfree(recon->RRtwfs);
		recon->RRtwfs=dcellpinv2(recon->GRtwfs, neai, 1e-3, 0.1);
	}

	if(parms->save.setup){
		if(recon->RRtwfs) writebin(recon->RRtwfs, "twfs_recon");
	}
	cellfree(neai);
}

/**
   TWFS has output. Accumulate result to simu->gradoff. It is put in wfsgrad.c
   instead of recon.c to avoid race condition because it updates simu->gradoff.
*/
void twfs_recon(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int itpowfs=parms->itpowfs;
	if(simu->wfsflags[itpowfs].gradout){
		gradoff_acc(simu, parms->ilgspowfs);//todo: improve ipowfs index.
		const int nlayer=NY(recon->GRall);
		dcell* Rmod=0;
		//Build radial mode error using closed loop TWFS measurements from this time step.
		dcellmm(&Rmod, recon->RRtwfs, simu->gradcl, "nn", 1);
		if(simu->wfsflags[itpowfs].gradout<5&&parms->itwfssph>-1 && fabs(parms->sim.eptsph/simu->eptwfs-1)>0.01){
			dbg("Step %5d: TWFS output %d spherical mode (%d) gain is boosted from %g to %g\n",
				simu->wfsisim, simu->wfsflags[itpowfs].gradout, parms->itwfssph, simu->eptwfs, parms->sim.eptsph);
			for(int ilayer=0; ilayer<nlayer; ilayer++){
				P(P(Rmod, ilayer), parms->itwfssph)*=(parms->sim.eptsph/simu->eptwfs);
			}
		}
		if(parms->save.extra){
			zfarr_push(simu->save->restwfs, simu->wfsflags[itpowfs].gradout-1, Rmod);
		}
		if(parms->recon.twfs_offsetdm){
			info_once("Step %5d: TWFS output to dmoff with gain %g every %d steps.\n", simu->wfsisim, simu->eptwfs, parms->powfs[itpowfs].dtrat);
			if(!simu->dmoff){
				simu->dmoff=dcellnew(parms->ndm, 1);
			}
			for(int ilayer=0; ilayer<nlayer; ilayer++){
				dmm(&P(simu->dmoff, ilayer), 1, P(recon->Rmod, ilayer), P(Rmod, ilayer), "nn", simu->eptwfs);
			}
			if(parms->plot.run){
				plot_dm(parms, recon, simu->dmoff, 0, "DM Offset", "Offset");
			}
		}else{
			info_once("Step %5d: TWFS output to gradcl with gain %g every %d steps.\n", simu->wfsisim, simu->eptwfs, parms->powfs[itpowfs].dtrat);
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				for(int ilayer=0; ilayer<nlayer; ilayer++){
					if(parms->powfs[ipowfs].skip!=2 && P(recon->GRall, iwfs, ilayer)){
						dmm(&P(simu->gradoff, iwfs), 1, P(recon->GRall, iwfs, ilayer), P(Rmod, ilayer), "nn", -simu->eptwfs);
					}
				}
				if(parms->plot.run){
					plot_gradoff(simu, iwfs);
				}
			}
		}
		if(parms->save.gradoff){
			writebin(simu->gradoff, "extra/gradoff_%d_twfs", simu->wfsisim);
		}
		dcellfree(Rmod);
	}
}