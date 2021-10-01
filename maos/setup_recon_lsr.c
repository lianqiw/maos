/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "setup_recon.h"
#include "ahst.h"

/**
   \file setup_recon_lsr.c

   Contains routines that setup the least squares reconstructor

*/
/**
   Setup the least square reconstruct by directly inverting GA matrix.
   The reconstructor is simply the pseudo inverse of GA matrix:
   \f[\hat{x}=(G_a^TC_g^{-1}G_a)^{-1}G_a^TC_g^{-1}\f]

   This is very close to RR except replacing GX with GA.

   We use the tomograhy parameters for lsr, since lsr is simply "tomography" onto DM directly.
*/
void setup_recon_lsr(recon_t* recon, const parms_t* parms){
	const int ndm=parms->ndm;
	const int nwfs=parms->nwfsr;
	cell* GAlsr;
	cell* GAM=parms->recon.modal?(cell*)recon->GM:(cell*)recon->GA;
	if(parms->recon.split){ //high order wfs only in split mode. 
		GAlsr=parms->recon.modal?(cell*)recon->GMhi:(cell*)recon->GAhi;
	} else{ //all wfs in integrated mode. 
		GAlsr=GAM;
	}
	info("Building recon->LR\n");
	int free_GAlsr=0;
	/*if(P(GAlsr,0)->id!=M_REAL){//Convert low sparsity matrices to full
		dsp* tmp=dsp_cast(P(GAlsr,0));
		if(tmp->nzmax>NX(tmp)*NY(tmp)*0.2){//not very sparse
			dcell* tmp2=0;
			free_GAlsr=1;
			dcelladd(&tmp2, 1, (dspcell*)GAlsr, 1);
			GAlsr=(cell*)tmp2;
		}
	}
	recon->LR.M=dcellmm2(GAlsr, recon->saneai, "tn");*/
	dcellmm_any(&recon->LR.M, GAlsr, CELL(recon->saneai), "tn", 1);
	// Tip/tilt and diff focus removal low rand terms for LGS WFS.
	if(recon->TTF){
		dcellmm_cell(&recon->LR.U, recon->LR.M, recon->TTF, "nn", 1);
		recon->LR.V=dcelltrans(recon->PTTF);
	}

	info("Building recon->LL\n");
	dcellmm_any(&recon->LL.M, recon->LR.M, GAlsr, "nn", 1);
	if(free_GAlsr){
		cellfree(GAlsr);
	}
	if(recon->LR.U){
		recon->LL.U=dcelldup(recon->LR.U);
		dcell* GPTTDF=NULL;
		dcellmm_cell(&GPTTDF, GAM, recon->LR.V, "tn", 1);
		recon->LL.V=dcelldup(GPTTDF);
		dcellfree(GPTTDF);
	}
	real maxeig=pow(recon->neamhi*P(recon->aloc,0)->dx, -2);
	if(parms->recon.modal){
		real strength=1;
		for(int idm=0; idm<ndm; idm++){
			strength*=dnorm(P(recon->amod,idm));
		}
		strength=pow(strength, 2./ndm);
		maxeig*=strength;
	}
	if(fabs(parms->lsr.tikcr)>EPS&&parms->lsr.actslave<=1){
		info2("Adding tikhonov constraint of %g to LLM\n", parms->lsr.tikcr);
		info("The maximum eigen value is estimated to be around %g\n", maxeig);
		dcelladdI_any(recon->LL.M, parms->lsr.tikcr*maxeig);
	}
	dcell* NW=NULL;
	if(!parms->recon.modal){
		if(parms->lsr.alg!=2){
			info("Create piston and check board modes that are in NULL space of GA.\n");
			/* Not SVD, need low rank terms for piston/waffle mode constraint. */
			NW=dcellnew(ndm, 1);
			int nmod=2;/*two modes. */
			for(int idm=0; idm<ndm; idm++){
				const long nloc=P(recon->aloc,idm)->nloc;
				P(NW,idm)=dnew(nloc, ndm*nmod);
				real* p=PCOL(P(NW,idm), idm*nmod);

				//First mode: piston mode.
				const real* cpl=P(P(recon->actcpl,idm));
				for(long iloc=0; iloc<nloc; iloc++){
					if(cpl[iloc]>0.1){
						p[iloc]=1;
					}
				}

				//Second mode: waffle mode
				p=PCOL(P(NW,idm), idm*nmod+1);
				loc_create_map(P(recon->aloc,idm));
				map_t* map=P(recon->aloc,idm)->map;
				for(long iy=0; iy<NY(map); iy++){
					for(long ix=0; ix<NX(map); ix++){
						if(P(map, ix, iy)>0){//Some may be negative due to extend.
							p[(long)P(map, ix, iy)-1]=(real)2*((iy+ix)&1)-1;
						}
					}
				}
			}
			/*scale it to match the magnitude of LL.M */
			dcellscale(NW, sqrt(maxeig));
			if(parms->save.setup){
				writebin(NW, "lsrNW");
			}

			dcellcat2(&recon->LL.U, NW, 2);
			dcellscale(NW, -1);
			dcellcat2(&recon->LL.V, NW, 2);
			dcellfree(NW);
		}

		if(parms->lsr.actslave){
			/*actuator slaving. important. change from 0.5 to 0.1 on 2011-07-14. */
			dspcell* actslave=slaving(recon->aloc, recon->actcpl,
				parms->dbg.recon_stuck?recon->actstuck:0,
				recon->actfloat, parms->lsr.actthres, maxeig, 1);
			if(parms->save.setup){
				writebin(actslave, "lsr_actslave");
			}
			dcelladd_any(&recon->LL.M, 1, CELL(actslave), 1);
			cellfree(actslave);
		}
		if(parms->lsr.actslave>1){
			/*per island regularization to mitigate the island effect*/
			dspcell* actslave=slaving(recon->aloc, recon->actcpl,
				parms->dbg.recon_stuck?recon->actstuck:0,
				recon->actfloat, parms->lsr.actthres2, parms->lsr.tikcr*maxeig, parms->lsr.actslave);
			if(parms->save.setup){
				writebin(actslave, "lsr_actslave2");
			}
			dcelladd_any(&recon->LL.M, 1, CELL(actslave), 1);
			cellfree(actslave);
		}
	}
	if(parms->recon.split==0&&parms->nlowfs){
	/*Low rank terms for low order wfs. Only in Integrated tomography with low order WFS. */
		dcell* ULo=dcellnew(ndm, nwfs);
		dcell* VLo=dcellnew(ndm, nwfs);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip||!parms->powfs[ipowfs].lo){
				continue;
			}
			for(int idm=0; idm<ndm; idm++){
				dspfull(&P(ULo, idm, iwfs), (dsp*)P(recon->LR.M, idm, iwfs), 'n', -1);
				dspfull(&P(VLo, idm, iwfs), (dsp*)P(GAM, iwfs, idm), 't', 1);
			}
		}


		dcellcat2(&recon->LL.U, ULo, 2);
		dcellcat2(&recon->LL.V, VLo, 2);

		dcellfree(ULo);
		dcellfree(VLo);
	}
	if(parms->lsr.fnreg){
		warning("Loading LSR regularization from file %s.\n", parms->lsr.fnreg);
		dspcell* tmp=dspcellread("%s", parms->lsr.fnreg);
		if(tmp) dcelladd_any(&recon->LL.M, 1, CELL(tmp), 1);
		dspcellfree(tmp);
	}
	recon->LL.alg=parms->lsr.alg;
	recon->LL.bgs=parms->lsr.bgs;
	recon->LL.warm=parms->recon.warm_restart;
	recon->LL.maxit=parms->lsr.maxit;
	/*Remove empty cells. */
	dcelldropempty(&recon->LR.U, 2);
	dcelldropempty(&recon->LR.V, 2);
	dcelldropempty(&recon->LL.U, 2);
	dcelldropempty(&recon->LL.V, 2);
	if(parms->save.recon){
		writecell(recon->LR.M, "LRM");
		writebin(recon->LR.U, "LRU");
		writebin(recon->LR.V, "LRV");
		writecell(recon->LL.M, "LLM.bin");/*disable compression */
		writebin(recon->LL.U, "LLU");
		writebin(recon->LL.V, "LLV");
	}
	if(parms->lsr.alg==0||parms->lsr.alg==2){
		if(!parms->lsr.bgs){
			muv_direct_prep(&recon->LL, (parms->lsr.alg==2)*parms->lsr.svdthres);
			if(parms->save.recon){
				if(recon->LL.C)
					chol_save(recon->LL.C, "LLC.bin");
				else
					writebin(recon->LL.MI, "LLMI.bin");
			}
			cellfree(recon->LL.M);
			dcellfree(recon->LL.U);
			dcellfree(recon->LL.V);
		} else{
			muv_direct_diag_prep(&(recon->LL), (parms->lsr.alg==2)*parms->lsr.svdthres);
			if(parms->save.recon){
				for(int ib=0; ib<recon->LL.nb; ib++){
					if(recon->LL.CB)
						chol_save(recon->LL.CB[ib], "LLCB_%d.bin", ib);
					else
						writebin(recon->LL.MI, "LLMIB_%d.bin", ib);
				}
			}
			/*Don't free M, U, V */
		}
	}
}
