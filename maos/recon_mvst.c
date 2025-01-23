/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
 * \file recon_mvst.c
 * Setup Minimum variance split tomography. 
 * This function is deprecated.
 */
/**
   compute the MVST split tomography NGS mode reconstructor.

   2010-03-16:
   New Implementation.

   Definition:

   - \f$G_{lgs}\f$ is LGS gradient operator from xloc, contained in recon->GXtomo as dspcell
   - \f$G_{ngs}\f$ is NGS gradient operator from xloc, contained in recon->GXL as dense matrix
   - \f$C_{lgs}, C_{ngs}\f$ is the LGS, and NGS measurement noise covariance matrix.
   - \f$\hat{x}_{lgs}\f$ is LGS tomography output
   - \f$\hat{x}_{ngs}\f$ is NGS minimum variance split tomography output
   - \f$a_{ngs}\f$ is NGS minimum variance DM output.
   - \f$F\f$ is the fitting operator
   - \f$H_A\f$ is the ray tracing operator from aloc to ploc, contained in recon->HA.
   - \f$MVModes\f$ is the MVST NGS modes
   - \f$MVRngs\f$ is the MVST NGS reconstructor

   We have

   \f{eqnarray*}{
   \hat{x}_{lgs}&=&A^{-1}G^{T}_{lgs}C_{lgs}^{-1}s_{lgs}\\
   Uw&=&A^{-1}G_{ngs}^TC_{ngs}^{-1}\\
   \hat{x}_{ngs}&=&Uw(1+G_{ngs}Uw)^{-1}(s_{ngs}-G_{ngs}\hat{x}_{lgs});\\
   a_{NGS}&=&F\hat{x}_{ngs}\\
   MVModes&=&F\cdot Uw\\
   MVRngs&=&(1+G_{ngs}Uw)^{-1}
   \f}

   If we want to orthnormalize the NGS modes. Propagate the NGS modes \f$F\cdot Uw\f$ to
   fitting directions and compute the cross-coupling matrix with weighting on \f$W_0, W_1\f$, we have

   \f{eqnarray*}{
   Q&=&H_A\cdot F \cdot Uw\\
   M_{CC}&=&Q^T(W_0-W_1 W_1^T)Q\\
   \f}

   Do SVD on \f$MCC\f$ we have \f$M_{CC}=u\Sigma v^T \f$ where \f$u{\equiv}v\f$
   because \f$MCC\f$ is symmetric. Redefine the NGS modes and reconstructor as

   \f{eqnarray*}{
   MVModes&=&F\cdot Uw\cdot u \Sigma^{-1/2}\\
   MVRngs&=&\Sigma^{1/2}u^T (1+G_{ngs}A^{-1}G_{ngs}^T C_{ngs}^{-1})^{-1}
   \f}
*/

void
setup_recon_mvst(recon_t* recon, const parms_t* parms){
	TIC;tic;
	/*
	  Notice that: Solve Fitting on Uw and using FUw to form Rngs gives
	  slightly different answer than solve fitting after assemble the
	  reconstructor. 10^-6 relative difference.

	  2010-03-10: Bug found and fixed: The MVST with CBS-CBS method gives worst
	  performance than integrated tomography. The probelm is in PUm
	  computing. I mistakenly called chol_solve, while I should have
	  called muv_direct_solve. The former doesn ot apply the low rank
	  terms.
	*/
	if(parms->recon.split!=2){
		return;
	}
	cellfree(recon->MVRngs);
	cellfree(recon->MVModes);
	cellfree(recon->MVGM);
	cellfree(recon->MVFM);

	dcellfree(recon->GXL);
	dcelladd(&recon->GXL, 1, recon->GXlo, 1);
	//NEA of low order WFS.
	dcell* neailo=dcellnew(parms->nwfsr, parms->nwfsr);
	dcell* nealo=dcellnew(parms->nwfsr, parms->nwfsr);
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(parms->powfs[ipowfs].lo){
			dspfull(&P(neailo,iwfs,iwfs),P(recon->saneai,iwfs,iwfs), 'n', 1);
			dspfull(&P(nealo,iwfs,iwfs), P(recon->sanea,iwfs,iwfs), 'n', 1);
		}
	}
	/* 2012-03-21: Remove focus mode from GL and NEA so no focus is
	   measured/estimated by mvst. It is then estimate separately*/
	if(parms->sim.mffocus){
		dcell* focus=NULL;
		dcellmm(&focus, recon->RFngsg, recon->GXL, "nn", 1);
		dcellmm(&recon->GXL, recon->GFngs, focus, "nn", -1);
		//dcelldropzero(recon->GXL, 1e-12);
		dcellfree(focus);
		dcellmm(&focus, recon->RFngsg, neailo, "nn", 1);
		dcellmm(&neailo, recon->GFngs, focus, "nn", -1);
		//dcelldropzero(neailo, 1e-8);
		dcellfree(focus);
	}

	dcell* U=NULL;
	dcell* FU=NULL;
	if(parms->load.mvst){
		U=dcellread("mvst_U");
		FU=dcellread("mvst_FU");
	} else{
	/*Prepare CBS if not already done */
		if(!recon->RL.C&&!recon->RL.MI){
			muv_direct_prep(&(recon->RL), 0);
		}
		if(!recon->fit->FL.C&&!recon->fit->FL.MI){
			muv_direct_prep(&(recon->fit->FL), 0);
		}
		toc2("MVST: svd prep");
		dcell* GXLT=dcelltrans(recon->GXL);
		muv_direct_solve(&U, &recon->RL, GXLT);
		dcellfree(GXLT);
		dcell* rhs=NULL;
		muv(&rhs, &recon->fit->FR, U, 1);
		muv_direct_solve(&FU, &recon->fit->FL, rhs);
		dcellfree(rhs);
		toc2("MVST: U, FU");

		if(parms->save.mvst||parms->save.setup){
			writebin(U, "mvst_U");
			writebin(FU, "mvst_FU");
		}
	}
	dcell* Uw=NULL;
	dcell* FUw=NULL;

	dcellmm(&Uw, U, neailo, "nn", 1);
	dcellmm(&FUw, FU, neailo, "nn", 1);

	dcell* M=NULL;
	dcellmm(&M, recon->GXL, Uw, "nn", 1);
	dcelladdI(M, 1);
	dcell* Minv=dcellinv(M);
	dcellfree(M);
	if(parms->sim.mffocus){
	//Embed a focus removal. Necessary!
		dcell* focus=NULL;
		dcellmm(&focus, Minv, recon->GFngs, "nn", 1);
		dcellmm(&Minv, focus, recon->RFngsg, "nn", -1);
	}
	if(parms->save.setup){
		writebin(Minv, "mvst_Rngs_0");
		writebin(FUw, "mvst_Modes_0");
	}
	/*Compute the cross coupling matrix of the Modes:
	  FUw'*Ha'*W*Fuw*Ha. Re-verified on 2013-03-24.*/
	dcell* QwQc=NULL;
	{
		dcell* Q=NULL;/*the NGS modes in ploc. */
		dcellmm(&Q, recon->fit->HA, FUw, "nn", 1);
		QwQc=calcWmcc(Q, Q, recon->W0, recon->W1, parms->fit.wt);
		dcellfree(Q);
	}
	/*Compute the wavefront error due to measurement noise. Verified on
	  2013-03-24. The gain optimization is yet a temporary hack because the way
	  PSDs are input.*/
	if(0){
		dcell* RC=NULL;
		dcellmm(&RC, Minv, nealo, "nn", 1);
		dcell* RCRt=NULL;
		dcellmm(&RCRt, RC, Minv, "nt", 1);
		dcell* RCRtQwQ=NULL;
		dcellmm(&RCRtQwQ, RCRt, QwQc, "nn", 1);
		dmat* tmp=dcell2m(RCRtQwQ);
		dmat* ptmp=tmp/*PDMAT*/;
		real rss=0;
		for(int i=0; i<NX(tmp); i++){
			rss+=P(ptmp, i, i);
		}
		dfree(tmp);
		dcellfree(RCRtQwQ);
		dcellfree(RCRt);
		dcellfree(RC);
		recon->sigmanlo=rss;
		dbg("rms=%g nm\n", sqrt(rss)*1e9);

		if(zfexist("../../psd_ngs.bin")){
			warning("Temporary solution for testing\n");
			dmat* psd_ngs=dread("../../psd_ngs.bin");
			if(parms->sim.wspsd){//windshake
				//need to convert from rad to m2.
				dmat* psd_ws_m=ddup(parms->sim.wspsd);
				dmat* psd_ws_y=drefcols(psd_ws_m, 1, 1);
				dscale(psd_ws_y, 4./parms->aper.d); dfree(psd_ws_y);
				add_psd2(&psd_ngs, psd_ws_m, 1); dfree(psd_ws_m);
			}
			writebin(psd_ngs, "psd_ngs_servo");
			dmat* rss2=dnew(1, 1); P(rss2,0)=rss;
			int dtrat=parms->powfs[P(parms->lopowfs,0)].dtrat;
			dcell* res=servo_optim(parms->sim.dt,
				dtrat, parms->sim.allo, M_PI/4, 0, 0, 2, psd_ngs, rss2);
			dfree(rss2);
			dbg("dtrat=%d\n", dtrat);
			dbg("g,a,T was %g,%g,%g\n", P(parms->sim.eplo,0), P(parms->sim.eplo,1), P(parms->sim.eplo,2));
			memcpy(P(parms->sim.eplo), P(P(res,0)), 3*sizeof(real));
			dbg("g,a,T=%g,%g,%g\n", P(parms->sim.eplo,0), P(parms->sim.eplo,1), P(parms->sim.eplo,2));
			dbg("res=%g, resn=%g nm\n", sqrt(P(P(res,0),3))*1e9, sqrt(P(P(res,0),4))*1e9);
			dcellfree(res);
			dfree(psd_ngs);
		}
	}
	if(1){/*Orthnormalize the Modes.*/
	/*
	  Change FUw*Minv -> FUw*(U*sigma^-1/2) * (U*sigma^1/2)'*Minv
	  columes of FUw*(U*sigma^-1/2) are the eigen vectors.

	  U, sigma is the eigen value decomposition of <FUw' HA' W HA FUw>
	*/

		dmat* QSdiag=NULL, * QU=NULL, * QVt=NULL;
		{
			dmat* QwQ=dcell2m(QwQc);
			dsvd(&QU, &QSdiag, &QVt, QwQ);
			dfree(QwQ);
		}
		if(parms->save.setup){
			writebin(QSdiag, "mvst_QSdiag");
		}
		dcwpow_thres(QSdiag, -1./2., 1e-14);
		dmuldiag(QU, QSdiag);/*U*sigma^-1/2 */
		d2cell(&QwQc, QU, NULL);
		dcell* FUw_keep=FUw;FUw=NULL;
		dcellmm(&FUw, FUw_keep, QwQc, "nn", 1);
		dcellfree(FUw_keep);
		dcwpow_thres(QSdiag, -2, 1e-14);
		dmuldiag(QU, QSdiag);/*U*sigma^1/2 (From U*sigma^(-1/2)*sigma) */
		d2cell(&QwQc, QU, NULL);
		dcell* Minv_keep=Minv; Minv=NULL;
		dcellmm(&Minv, QwQc, Minv_keep, "tn", 1);
		dcellfree(Minv_keep);
		dfree(QSdiag);
		dfree(QVt);
		dfree(QU);
	}
	dcellfree(QwQc);

	recon->MVRngs=dcellreduce(Minv, 1);/*1xnwfs cell */
	recon->MVModes=dcellreduce(FUw, 2);/*ndmx1 cell */
	dcellmm(&recon->MVGM, recon->GAlo, recon->MVModes, "nn", 1);
	dcellmm(&recon->MVFM, recon->RFngsg, recon->MVGM, "nn", 1);
	dcellfree(neailo);
	dcellfree(nealo);
	dcellfree(Minv);
	dcellfree(U);
	dcellfree(FU);
	dcellfree(FUw);
	dcellfree(Uw);
	if(parms->save.setup){
		dcell* Qn=NULL;
		dcellmm(&Qn, recon->fit->HA, recon->MVModes, "nn", 1);
		dcell* Qntt=dcellnew(NX(Qn), NY(Qn));
		dmat* TTploc=loc2mat(recon->floc, 1);/*TT mode. need piston mode too! */
		dmat* PTTploc=dpinv(TTploc, recon->W0);/*TT projector. no need w1 since we have piston. */
		dfree(TTploc);
		for(int ix=0; ix<NX(Qn)*NY(Qn); ix++){
			if(!P(Qn,ix)) continue;
			dmm(&P(Qntt,ix), 0, PTTploc, P(Qn,ix), "nn", 1);
		}
		writebin(Qntt, "mvst_modptt");
		dcellfree(Qn);
		dcellfree(Qntt);
		dfree(PTTploc);
		writebin(recon->MVRngs, "mvst_Rngs");
		writebin(recon->MVModes, "mvst_Modes");
	}
	if(parms->dbg.mvstlimit>0){/*limit number of modes used. */
		warning("MVST: Correction is limited to %d modes\n", parms->dbg.mvstlimit);
		dmat* tmp;
		for(int iy=0; iy<NY(recon->MVRngs); iy++){
			tmp=P(recon->MVRngs,iy);
			if(tmp){
				P(recon->MVRngs,iy)=dsub(tmp, 0, parms->dbg.mvstlimit, 0, NY(tmp));
				dfree(tmp);
			}
		}
		for(int ix=0; ix<NX(recon->MVModes); ix++){
			tmp=P(recon->MVModes,ix);
			if(tmp){
				P(recon->MVModes,ix)=dsub(tmp, 0, NX(tmp), 0, parms->dbg.mvstlimit);
				dfree(tmp);
			}
		}
		if(parms->save.setup){
			writebin(recon->MVRngs, "mvst_Rngs_limit");
			writebin(recon->MVModes, "mvst_Modes_limit");
		}
	}
	/*
	if(parms->save.setup){
	dcell *QQ=NULL;
	dcellmm(&QQ, recon->fit.HA, recon->MVModes,"nn", 1);
	dcell *MCC=calcWmcc(QQ,QQ,recon->W0,recon->W1,recon->fitwt);
	writebin(MCC,"mvst_MCC");

	dcellfree(MCC);
	dcellfree(QQ);
	}*/
	toc2("MVST");
}
