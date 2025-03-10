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
#include "recon.h"
#include "recon_utils.h"
/**
   \file recon_fit.c

   Contains routines that setup the DM fitting.

*/
/**
   Setup ray tracing operator HXF from xloc to aperture ploc along DM fiting directions*/
static dspcell*
setup_fit_HXF(const fit_t* fit){
	TIC;tic;
	if(!fit->xloc) return 0;
	const int nfit=NX(fit->thetax);
	const int npsr=NX(fit->xloc);
	dspcell* HXF=dspcellnew(nfit, npsr);
OMP_FOR_COLLAPSE(2, NTHREAD)
	for(int ifit=0; ifit<nfit; ifit++){
		for(int ips=0; ips<npsr; ips++){
			const real hsi=P(fit->hs,ifit);
			const real ht=P(fit->xloc,ips)->ht-fit->floc->ht;
			const real scale=1.-ht/hsi;
			real displace[2];
			displace[0]=P(fit->thetax,ifit)*ht;
			displace[1]=P(fit->thetay,ifit)*ht;
			P(HXF, ifit, ips)=mkh(P(fit->xloc,ips), fit->floc, displace[0], displace[1], scale, 0);
		}
	}
	toc2("HXF");
	return HXF;
}

/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static void
setup_fit_HA(fit_t* fit){
	const int nfit=NX(fit->thetax);
	const int ndm=NX(fit->aloc);
	fit->HA=cellnew(nfit, ndm);
	TIC;tic;
OMP_FOR_COLLAPSE(2, NTHREAD)
	for(int ifit=0; ifit<nfit; ifit++){
		for(int idm=0; idm<ndm; idm++){
			const real hs=P(fit->hs,ifit);
			const real ht=P(fit->aloc,idm)->ht-fit->floc->ht;
			const real scale=1.-ht/hs;
			real displace[2];
			displace[0]=P(fit->thetax,ifit)*ht;
			displace[1]=P(fit->thetay,ifit)*ht;
			loc_t* loc=fit->floc;
			if(fit->misreg&&fit->misreg[ifit+idm*nfit]){
				loc=loctransform(loc, fit->misreg[ifit+idm*nfit]);
			}
			dsp *ha=mkh(P(fit->aloc,idm), loc, displace[0], displace[1], scale, 0);
			if(fit->modal){
				dspmm((dmat**)&P(fit->HA, ifit, idm), ha, P(fit->amod, idm, idm), "nn", 1);
				dspfree(ha);
			}else{
				P(fit->HA, ifit, idm)=(cell*)ha;
			}
			if(loc!=fit->floc){
				locfree(loc);
			}
		}
	}
	toc2("HA");
	if(!fit->modal){
		fit->actcpl=genactcpl(fit->HA, fit->W1);
		if(fit->actfloat){//zero coupling for floating actuators. Do not zero for sutck actuators.
			act_stuck(fit->aloc, fit->actcpl, fit->actfloat);
		}
		if(global->parms->dbg.recon_stuck){
			act_stuck(fit->aloc, fit->HA, fit->actstuck);
		}

		if(fit->flag.actextrap){
			/*
				DM fitting output a is extrapolated to edge actuators by
				actextrap*a. The corresponding ray tracing from DM would be
				HA*actextrap*a. We replace HA by HA*actextrap to take this into
				account during DM fitting. 
				This needs to be used with extrapolation on dmpsol.
			*/
			fit->actextrap=act_extrap(fit->aloc, fit->actcpl, fit->flag.actthres, 0);
			dbg("Replacing HA by HA*fit->actextrap\n");
			cell *HA2=0;
			dcellmm(&HA2, fit->HA, fit->actextrap, "nn", 1);
			cellfree(fit->HA);
			fit->HA=HA2;
		} else if(fit->actfloat){//avoid commanding floating actuators
			act_float(fit->aloc, (dspcell**)&fit->HA, NULL, fit->actfloat);
		}
	}
}
/**
   Setup fitting low rank terms that are in the NULL space of DM fitting
   operator. typically include piston on each DM and tip/tilt on certain
   DMs. Becareful with tip/tilt contraint when using CBS.  */
static void
setup_fit_lrt(fit_t* fit){
	const int ndm=NX(fit->aloc);
	fit->NW=dcellnew(ndm, 1);
	//real fitscl;     /**<strength of fitting FLM low rank terms (vectors)*/
	real fitscl=1./fit->floc->nloc;
	if(fabs(fitscl)<1.e-15){
		error("fit->fitscl is too small\n");
	}
	int nnw=0;
	if(fit->flag.lrt_piston){
		nnw+=ndm;
	}
	if(fit->flag.lrt_tt){
		nnw+=2*(ndm-1);
	}
	if(nnw==0) return;
	dcell* actcpl=dcelldup(fit->actcpl);
	if(global->parms->dbg.recon_stuck){
	//avoid stuck actuators for piston constraint.
		act_stuck(fit->aloc, actcpl, fit->actstuck);
	}
	for(int idm=0; idm<ndm; idm++){
		int nloc=P(fit->aloc,idm)->nloc;
		P(fit->NW,idm)=dnew(nloc, nnw);
	}
	int inw=0;/*current column */
	if(fit->flag.lrt_piston){
		info("Adding piston constraint to fit matrix\n");
		for(int idm=0; idm<ndm; idm++){
			int nloc=P(fit->aloc,idm)->nloc;
			real* p=P(P(fit->NW,idm))+(inw+idm)*nloc;
			const real* cpl=actcpl?P(P(actcpl,idm)):NULL;
			for(int iloc=0; iloc<nloc; iloc++){
				if(!cpl || cpl[iloc]>0.1){ //don't count floating or stuck actuators
					p[iloc]=fitscl;
				}
			}
		}
		inw+=ndm;
	}
	if(fit->flag.lrt_tt){
		real factor=0;
		info("Adding TT constraint on upper DMs to fit matrix.\n");
		factor=fitscl*2./loc_diam(P(fit->aloc,0));
		for(int idm=1; idm<ndm; idm++){
			int nloc=P(fit->aloc,idm)->nloc;
			real* p=P(P(fit->NW,idm))+(inw+(idm-1)*2)*nloc;
			real* p2x=p;
			real* p2y=p+nloc;
			const real* cpl=P(P(actcpl,idm));
			for(int iloc=0; iloc<nloc; iloc++){
				if(cpl[iloc]>0.1){
					p2x[iloc]=P(fit->aloc,idm)->locx[iloc]*factor;/*x tilt */
					p2y[iloc]=P(fit->aloc,idm)->locy[iloc]*factor;/*y tilt */
				}
			}
		}
		inw+=2*(ndm-1);
	}
	if(fit->flag.actslave){
		//TIC;tic;
		fit->actslave=slaving(fit->aloc, fit->actcpl,
			global->parms->dbg.recon_stuck?fit->actstuck:0,
			fit->actfloat, fit->flag.actthres, 1./fit->floc->nloc, 1);

		if(fit->flag.actslave>1){
			dspcell* actslave2=slaving(fit->aloc, fit->actcpl,
				global->parms->dbg.recon_stuck?fit->actstuck:0,
				fit->actfloat, fit->flag.actthres, 1./fit->floc->nloc, 2);
			dcelladd(&fit->actslave, 1, actslave2, 1);
			cellfree(actslave2);
		}
		//toc2("slaving");
	}
	cellfree(actcpl);
}
/**
   Assemble the DM fitting matrix

   The fitting is done by minimizing \f$||H_X x - H_A a||^2_W\f$ where \f$H_X,
   H_A\f$ are ray tracing operator from tomography grid xloc, and deformable
   mirror grid aloc to pupil grid ploc. The norm is weighted using bilinear
   influence functions within the telescope aperture. We have

   \f$a=\left[H_A^T(W_0-W_1 W_1^T)H_A\right]^{-1} H_A^T (W_0-W_1) H_X x\f$

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803
*/
static void
setup_fit_matrix(fit_t* fit){
	const int nfit=NX(fit->thetax);
	const int ndm=NX(fit->aloc);
	if(ndm==0) return;

	cell* HA=fit->HA;
	//print_mem("Before assembling fit matrix");
	/*Assemble Fit matrix. */
	if(!fit->FR.M&&fit->flag.assemble){
		if(fit->HXF){//not idealtomo.
			const int npsr=NX(fit->xloc);
			dbg("Building fit->FR\n");
			fit->FR.M=cellnew(ndm, npsr);
			dspcell* FRM=(dspcell*)fit->FR.M;
			dspcell* HXF=fit->HXF;
			//FRM
			for(int ips=0; ips<npsr; ips++){
				for(int ifit=0; ifit<nfit; ifit++){
					if(fabs(P(fit->wt,ifit))<1.e-12) continue;
					dsp* tmp=dspmulsp(fit->W0, P(HXF, ifit, ips), "nn");
					for(int idm=0; idm<ndm; idm++){
						dcellmm(&P(FRM, idm, ips), P(HA, ifit, idm), tmp, "tn",
							P(fit->wt,ifit));
					}
					dspfree(tmp);
				}
			}
			fit->FR.V=dcellnew(npsr, 1);
			dmat** FRV=P(fit->FR.V);
			dmat frvi={0};
			//FRV
			for(int ips=0; ips<npsr; ips++){
				int nloc=P(fit->xloc,ips)->nloc;
				FRV[ips]=dnew(nloc, nfit);
				for(int ifit=0; ifit<nfit; ifit++){
					/*notice the sqrt. */
					if(fabs(P(fit->wt,ifit))<1.e-12) continue;
					dcols(&frvi, FRV[ips], ifit, 1);
					dspmv(&frvi, P(HXF, ifit, ips), fit->W1, 't', sqrt(P(fit->wt,ifit)));
				}
			}
			cellfree(fit->HXF);
		} else{
			dbg("Avoid building fit->FR.M for idealtomo\n");
			fit->FR.M=NULL;
			fit->FR.V=NULL;
		}
		/*Always need FR.U as it is used to do FL.U, FL.V */
		fit->FR.U=dcellnew(ndm, 1);
		dmat** FRU=P(fit->FR.U);
		dmat frui={0}; dmat *pfrui=&frui;
		for(int idm=0; idm<ndm; idm++){
			int nx=NY(P(HA,0,idm));
			FRU[idm]=dnew(nx, nfit);
			for(int ifit=0; ifit<nfit; ifit++){
			/*notice the sqrt. */
				if(fabs(P(fit->wt,ifit))<1.e-12) continue;
				dcols(&frui, FRU[idm], ifit, 1);
				dcellmm(&pfrui, P(HA, ifit, idm), fit->W1, "tn", sqrt(P(fit->wt,ifit)));
			}
		}
	}

	if(!fit->FL.M){
		dbg("Building fit->FL\n");//TIC;tic;
		fit->FL.M=cellnew(ndm, ndm);
		for(int idm=0; idm<ndm; idm++){
			for(int ifit=0; ifit<nfit; ifit++){
				if(fabs(P(fit->wt,ifit))<1.e-12) continue;
				cell*tmp=NULL;
				dcellmm(&tmp, fit->W0, P(fit->HA, ifit, idm), "nn", 1);
				for(int jdm=0; jdm<ndm; jdm++){
					dcellmm(&P(fit->FL.M, jdm, idm), P(fit->HA, ifit, jdm), tmp, "tn", P(fit->wt,ifit));
				}
				cellfree(tmp);
			}
		}
		//toc("FLM done...");
		if(fabs(fit->flag.tikcr)>1.e-15){
			real tikcr=fit->flag.tikcr;
			/*Estimated from the formula.  1/nloc is due to W0, the other
			  scaling is due to ray tracing between different sampling freq.*/
			int nact=0;
			for(int idm=0; idm<ndm; idm++){
				nact+=P(fit->aloc,idm)->nloc;
			}
			real maxeig=4./nact;
			dbg("Adding tikhonov constraint of %.1e to FLM\n", tikcr);
			dbg("The maximum eigen value is estimated to be around %.1e\n", maxeig);
			dcelladdI(fit->FL.M, tikcr*maxeig);
			//toc("addI done...");
		}

		{/*Low rank terms. */
			fit->FL.U=dcellcat_each(fit->FR.U, fit->NW, 2);
			dcell* tmp=NULL;/*negative NW. */
			dcelladd(&tmp, 1, fit->NW, -1);
			fit->FL.V=dcellcat_each(fit->FR.U, tmp, 2);
			dcellfree(tmp);
			//toc("uv done...");
		}
		if(fit->actslave){
			dcelladd(&fit->FL.M, 1, fit->actslave, 1);
			//toc("slaving done...");
		}
		/*dspcellsym(fit->FL.M); */
		dbg("DM Fit number of Low rank terms: %ld in LHS\n", P(fit->FL.U,0)->ny);
	}
	if(fit->flag.alg==0||fit->flag.alg==2){
		if(fit->flag.alg==0&&fabs(fit->flag.tikcr)<1.e-14){
			warning("tickcr=%g is too small, chol may fail.\n", fit->flag.tikcr);
		}
		if(fit->flag.bgs){
			muv_direct_diag_prep(&(fit->FL), (fit->flag.alg==2)*fit->flag.svdthres);
		} else{
			muv_direct_prep(&(fit->FL), (fit->flag.alg==2)*fit->flag.svdthres);
			if(0){
				writebin(fit->FL.M, "FLM");
				writebin(fit->FL.U, "FLU");
				writebin(fit->FL.V, "FLV");
			}
			cellfree(fit->FL.M);
			dcellfree(fit->FL.U);
			dcellfree(fit->FL.V);
		}
		//print_mem("After cholesky/svd on matrix");
	}

	//print_mem("After assemble fit matrix");
}
/**
   A generic DM fitting routine.
 */
void setup_fit(fit_t* fit, int idealtomo){
	TIC;tic;
	if(!idealtomo&&fit->xloc){
		fit->HXF=setup_fit_HXF(fit);
	}
	setup_fit_HA(fit);
	if(!fit->modal){
		setup_fit_lrt(fit);
	}
	/*always assemble fit matrix, faster if many directions */
	if(fit->flag.assemble||fit->flag.alg!=1){
		setup_fit_matrix(fit);
	}
	/*Fall back function method if FR.M is NULL (!HXF<-idealtomo) */
	fit->FR.Mfun=FitR;
	fit->FR.Mdata=fit;
	/*Fall back function method if FL.M is NULL */
	fit->FL.Mfun=FitL;
	fit->FL.Mdata=fit;
	fit->FL.alg=fit->flag.alg;
	fit->FL.bgs=fit->flag.bgs;
	fit->FL.warm=fit->flag.cgwarm;
	fit->FL.maxit=fit->flag.maxit;
	toc2("Setting up DM Fitting.");
}
void free_fit(fit_t* fit, int nfit){
	if(!fit) return;
	for(int ifit=0; ifit<nfit; ifit++){
		if(!fit[ifit].isref){
			cellfree(fit[ifit].HXF);
			cellfree(fit[ifit].HA);
			cellfree(fit[ifit].actcpl);
			cellfree(fit[ifit].actextrap);
			cellfree(fit[ifit].actslave);
			cellfree(fit[ifit].NW);
			muv_free(&fit[ifit].FR);
			muv_free(&fit[ifit].FL);
		}
	}
	free(fit);
}
/**
   Setup DM fitting parameters
*/
void setup_recon_fit(recon_t* recon, const parms_t* parms){
	fit_t* fit=mycalloc(1, fit_t);
	recon->fit=fit;
	fit->thetax=parms->fit.thetax;
	fit->thetay=parms->fit.thetay;
	fit->wt=parms->fit.wt;
	fit->hs=parms->fit.hs;

	fit->xloc=recon->xloc;
	fit->floc=recon->floc;
	fit->aloc=recon->aloc;
	fit->amod=recon->amod;
	fit->modal=parms->recon.modal;

	fit->W0=recon->W0;
	fit->W1=recon->W1;

	fit->actfloat=recon->actfloat;
	fit->actstuck=recon->actstuck;
	fit->misreg=parms->recon.distortion_dm2sci;
	memcpy(&fit->flag, &parms->fit, sizeof(fit_cfg_t));//use parms->fit.
	if(parms->fit.assemble){
		if(parms->load.fit){
			if(!(zfexist("FRM")&&zfexist("FRU")&&zfexist("FRV"))){
				error("FRM, FRU, FRV (.bin) not all exist\n");
			}
			if(!(zfexist("FLM")&&zfexist("FLU")&&zfexist("FLV"))){
				error("FLM, FLU, FLV (.bin) not all exist\n");
			}
			fit->FR.M=readbin("FRM");
			fit->FR.U=dcellread("FRU");
			fit->FR.V=dcellread("FRV");
			fit->FL.M=readbin("FLM");
			fit->FL.U=dcellread("FLU");
			fit->FL.V=dcellread("FLV");

		}
	}
	setup_fit(fit, parms->sim.idealtomo);
	if(fit->actcpl){
		dcellfree(recon->actcpl);
		recon->actcpl=dcellref(fit->actcpl);
	}
	if(fit->actextrap) {
		dcellfree(recon->actextrap);
		recon->actextrap=dspcellref(fit->actextrap);
	}
	if(parms->save.setup){
		writebin(fit->HA, "HA");
		if(fit->actcpl) writebin(fit->actcpl, "fit_actcpl");
		if(fit->actextrap) writebin(fit->actextrap, "fit_actextrap");
	}
	if(parms->save.recon){
		writebin(fit->FR.M, "FRM");
		writebin(fit->FR.V, "FRV");
		writebin(fit->FR.U, "FRU");

		writebin(fit->FL.M, "FLM");
		writebin(fit->FL.U, "FLU");
		writebin(fit->FL.V, "FLV");
		if(fit->FL.C){
		//chol_convert(fit->FL.C, 1);
			chol_save(fit->FL.C, "FLC.bin");
		}
		if(fit->FL.MI)
			writebin(fit->FL.MI, "FLMI");
		if(fit->FL.Up)
			writebin(fit->FL.Up, "FLUp");
		if(fit->FL.Vp)
			writebin(fit->FL.Vp, "FLVp");
		if(fit->FL.CB){
			for(int ib=0; ib<fit->FL.nb; ib++){
				chol_save(fit->FL.CB[ib], "FLCB_%d.bin", ib);
			}
		}
		if(fit->FL.MIB){
			writebin(fit->FL.MIB, "FLMIB");
		}
	}
	if(fit->flag.alg!=1){
		cellfree(fit->FL.M);
		dcellfree(fit->FL.U);
		dcellfree(fit->FL.V);
	}
}
/**
   Setting fitting parameter for turbulence to WFS lenslet grid.
 */
void setup_powfs_fit(powfs_t* powfs, const recon_t* recon, const parms_t* parms){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].lo) continue;
		int nwfs=parms->powfs[ipowfs].nwfs;
		fit_t* fitall=powfs[ipowfs].fit=mycalloc(nwfs, fit_t);
		loc_t* wfsloc=mkannloc(parms->aper.d+parms->powfs[ipowfs].dsa*2, 0, parms->powfs[ipowfs].dsa, 0);
		wfsloc->ht=0;
		wfsloc->iac=parms->dbg.wfs_iac;//cubic spline better fits the turbulence.
		for(int jwfs=0; jwfs<nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
			fit_t* fit=fitall+jwfs;
			if(jwfs==0){
				memcpy(&fit->flag, &parms->fit, sizeof(fit_cfg_t));//use parms->fit.
				fit->flag.alg=0;
				fit->flag.assemble=0;
				fit->notrecon=1; //not for reconstruction
				fit->wt=dnew(1, 1); P(fit->wt,0)=1;
				fit->hs=dnew(1, 1); P(fit->hs,0)=parms->wfs[iwfs].hs;
				fit->aloc=loccellnew(1, 1); P(fit->aloc,0)=locref(wfsloc);
				fit->floc=locref(recon->floc);
				fit->W0=recon->W0;
				fit->W1=recon->W1;

				fit->thetax=dnew(1, 1);P(fit->thetax,0)=parms->wfs[iwfs].thetax;
				fit->thetay=dnew(1, 1);P(fit->thetay,0)=parms->wfs[iwfs].thetay;
				setup_fit(fit, 1);
				if(fit->flag.alg!=1){
					cellfree(fit->FL.M);
					dcellfree(fit->FL.U);
					dcellfree(fit->FL.V);
				}
			} else{
				memcpy(fitall+jwfs, fitall, sizeof(fit_t));
				fit->isref=1;
				fit->FR.Mdata=fit;
				fit->FL.Mdata=fit;
				fit->thetax=dnew(1, 1);P(fit->thetax,0)=parms->wfs[iwfs].thetax;
				fit->thetay=dnew(1, 1);P(fit->thetay,0)=parms->wfs[iwfs].thetay;
			}
		}
		locfree(wfsloc);
	}
}
/**
   Call this instead of free_fit directly to free powfs fit parameters.
*/
void free_powfs_fit(powfs_t* powfs, const parms_t* parms){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		fit_t* fitall=powfs[ipowfs].fit;
		if(!fitall) continue;
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			fit_t* fit=fitall+jwfs;
			dfree(fit->thetax);
			dfree(fit->thetay);
			if(!fit->isref){
				dfree(fit->wt);
				dfree(fit->hs);
				cellfree(fit->aloc);
				cellfree(fit->floc);
			}
		}
		free_fit(powfs[ipowfs].fit, parms->powfs[ipowfs].nwfs);
	}
}


/**
   Apply fit right hand side matrix in CG mode without using assembled matrix.
   Slow. don't use. Assembled matrix is faster because of multiple directions.
*/
void FitR(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const fit_t* fit=(const fit_t*)A;
	const int nfit=NX(fit->thetax);
	dcell* xp=dcellnew(nfit, 1);

	if(!xin){/*xin is empty. We will trace rays from atmosphere directly */
		const parms_t* parms=global->parms;
		sim_t* simu=global->simu;
		int isim=fit->notrecon?simu->wfsisim:simu->reconisim;
		const real atmscale=simu->atmscale?P(simu->atmscale, isim):1;
		for(int ifit=0; ifit<nfit; ifit++){
			real hs=P(fit->hs, ifit);
			P(xp, ifit)=dnew(fit->floc->nloc, 1);
			for(int ips=0; ips<parms->atm.nps; ips++){
				const real ht=P(parms->atm.ht, ips)-fit->floc->ht;
				real scale=1-ht/hs;
				if(scale<0) continue;
				real dispx=P(fit->thetax, ifit)*ht-P(simu->atm, ips)->vx*isim*parms->sim.dt;
				real dispy=P(fit->thetay, ifit)*ht-P(simu->atm, ips)->vy*isim*parms->sim.dt;
				prop_grid(P(simu->atm, ips), fit->floc, P(P(xp, ifit)),
					atmscale, dispx, dispy, scale, 1, 0, 0);
			}
			/*if(simu->telws){//Wind shake. Enable after cuda code also has it.
				real tmp=P(simu->telws, isim);
				real angle=simu->winddir?P(simu->winddir, 0):0;
				real ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
				loc_add_ptt(P(xp, ifit), ptt, fit->floc);
			}*/
			//No need to utilize ncpa.surf or ncpa.tsurf if ncpa.calib=1
		}
	} else if(fit->HXF){
		dcellmm(&xp, fit->HXF, xin, "nn", 1.);
	} else{/*Do the ray tracing from xloc to ploc */
		const int npsr=NX(fit->xloc);
		for(int ifit=0; ifit<nfit; ifit++){
			real hs=P(fit->hs, ifit);
			P(xp, ifit)=dnew(fit->floc->nloc, 1);
			for(int ips=0; ips<npsr; ips++){
				const real ht=P(fit->xloc, ips)->ht-fit->floc->ht;
				real scale=1-ht/hs;
				if(scale<0) continue;
				real dispx=P(fit->thetax, ifit)*ht;
				real dispy=P(fit->thetay, ifit)*ht;
				prop_nongrid(P(fit->xloc, ips), P(P(xin, ips)), fit->floc,
					P(P(xp, ifit)), 1, dispx, dispy, scale, 0, 0);
			}
		}
	}
	applyW(xp, fit->W0, fit->W1, P(fit->wt));
	dcellmm(xout, fit->HA, xp, "tn", alpha);
	dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode without using assembled
   matrix. Slow. don't use. Assembled matrix is faster because of multiple
   directions.  */
void FitL(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const fit_t* fit=(const fit_t*)A;
	dcell* xp=NULL;
	dcellmm(&xp, fit->HA, xin, "nn", 1.);
	applyW(xp, fit->W0, fit->W1, P(fit->wt));
	dcellmm(xout, fit->HA, xp, "tn", alpha);
	dcellfree(xp);xp=NULL;
	dcellmm(&xp, fit->NW, xin, "tn", 1);
	dcellmm(xout, fit->NW, xp, "nn", alpha);
	dcellfree(xp);
	if(fit->actslave){
		dcellmm(xout, fit->actslave, xin, "nn", alpha);
	}
}
