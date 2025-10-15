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
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file recon_mvm.c

   Contains routines that setup the matrix for matrix vector multiply

*/
/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
static void
setup_recon_lsr_mvm(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	info("Assembling LSR MVM in CPU\n");
	dcell* MVM=NULL;
	if(recon->LR.Mfun||parms->lsr.alg==1){
	/*
	  First create an identity matrix. then solve each column one by one.
	*/
		const int ndm=parms->ndm;
		const int nwfs=parms->nwfsr;
		int ntotgrad=0;
		long* ngrad=P(recon->ngrad);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			ntotgrad+=powfs[ipowfs].saloc->nloc*2;
		}
		MVM=dcellnew(ndm, nwfs);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(!parms->powfs[ipowfs].skip){
				for(int idm=0; idm<ndm; idm++){
					P(MVM,idm,iwfs)=dnew(P(recon->anloc,idm), powfs[ipowfs].saloc->nloc*2);
				}
			}
		}

		dcell* res=NULL;
		int curg=0, curwfs=0;
		dmat* eye=dnew(ntotgrad, 1);
		dcell* eyec=d2cellref(eye, ngrad, nwfs);
		for(int ig=0; ig<ntotgrad; ig++){
			info_progress(ig, ntotgrad);
			if(ig) P(eye,ig-1)=0;
			P(eye,ig)=1;
			if(!parms->powfs[parms->wfsr[curwfs].powfs].skip){
				dcellzero(res);
				muv_solve(&res, &recon->LL, &recon->LR, eyec);
			}
			for(int idm=0; idm<ndm; idm++){
				dmat* to=P(MVM,idm,curwfs);
				if(to){
					memcpy(PCOL(to, curg), P(P(res,idm)), NX(to)*sizeof(real));
				}
			}
			curg++;
			if(curg>=ngrad[curwfs]){
				curwfs++;
				curg=0;
			}
		}
		dcellfree(res);
		dcellfree(eyec);
		dfree(eye);
	} else{
		dcell* LR=NULL;
		if(P(recon->LR.M,0)->id==M_REAL){
			LR=dcell_cast(recon->LR.M);
		} else{
			dspcellfull(&LR, dspcell_cast(recon->LR.M), 'n', 1);
		}
		if(recon->LR.U&&recon->LR.V){
			dcellmm(&LR, recon->LR.U, recon->LR.V, "nt", -1);
		}
		muv_solve(&MVM, &recon->LL, NULL, LR);
		if(LR!=(dcell*)recon->LR.M){
			dcellfree(LR);
		}
	}
	recon->MVM=dcell2m(MVM);
	dcellfree(MVM);
}


typedef struct{
	const parms_t* parms;
	recon_t* recon;
	dcell* MVMt;
	long(*curp)[2];
	long ntotact;
}MVR_MVM_T;
static void *
setup_recon_mvr_mvm_iact(thread_t* info){
	MVR_MVM_T* data=(MVR_MVM_T*)info->data;
	const parms_t* parms=data->parms;
	recon_t* recon=data->recon;
	const int ndm=parms->ndm;
	const int nwfs=parms->nwfsr;
	const long ntotact=data->ntotact;
	dcell* FLI=NULL;
	dcell* FRT=NULL;
	dcell* RLT=NULL;
	dcell* RRT=NULL;
	dmat* eye=dnew(ntotact, 1);
	dcell* eyec=d2cellref(eye, P(recon->anloc), ndm);
	long(*curp)[2]=data->curp;
	dcell* MVMt=data->MVMt;
	//int nthread=recon->nthread;
	TIC;tic;
	for(long iact=info->start; iact<info->end; iact++){
		int curdm=curp[iact][0];
		int curact=curp[iact][1];
		if(recon->actcpl&&P(P(recon->actcpl,curdm),curact)<EPS){
			continue;
		}
		//info_progress(iact, 64);
		//TIC;tic;
		dcellzero(FRT);
		dcellzero(RRT);
		/*Apply F_L*/
		P(eye,iact)=1;
		muv_solve(&FLI, &recon->fit->FL, NULL, eyec);
		P(eye,iact)=0;
		/*Apply F_R'*/
		muv_trans(&FRT, &recon->fit->FR, FLI, 1);
		//toc2("fit");
		/*Apply R_L*/
		dcellzero(RLT);//warm restart.
		muv_solve(&RLT, &recon->RL, NULL, FRT);
		/*Apply R_R'*/
		muv_trans(&RRT, &recon->RR, RLT, 1);
		//toc2("tomo");
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			dmat* to=P(MVMt,iwfs,curdm);
			if(to){
				memcpy(PCOL(to, curact), P(P(RRT,iwfs)), NX(to)*sizeof(real));
			}
		}
		info("Mean time is %.3f seconds\n", (myclockd()-tk)/(1+iact-info->start));
		fflush(stdout);
		/*{
			writebin(FLI, "cpu_dmrecon_%ld", iact);
			writebin(FRT, "cpu_opdx_%ld", iact);
			writebin(RLT, "cpu_opdr_%ld", iact);
			writebin(RRT, "cpu_grad_%ld", iact);
			}*/
	}
	dcellfree(FLI);
	dcellfree(FRT);
	dcellfree(RLT);
	dcellfree(RRT);
	dfree(eye);
	dcellfree(eyec);
	return NULL;
}


/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
static void
setup_recon_mvr_mvm(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	info("Assembling MVR MVM in CPU\n");
	const int ndm=parms->ndm;
	const int nwfs=parms->nwfsr;
	long ntotact=0;
	for(int idm=0; idm<ndm; idm++){
		ntotact+=P(recon->anloc,idm);
	}
	typedef long long2[2];
	long2* curp=mymalloc(ntotact, long2);
	int nact=0;
	for(int idm=0; idm<ndm; idm++){
		for(int iact=0; iact<P(recon->anloc,idm); iact++){
			curp[nact+iact][0]=idm;
			curp[nact+iact][1]=iact;
		}
		nact+=P(recon->anloc,idm);
	}
	dcell* MVMt=dcellnew(nwfs, ndm);
	for(int idm=0; idm<ndm; idm++){
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(!parms->powfs[ipowfs].skip){
				P(MVMt,iwfs,idm)=dnew(powfs[ipowfs].saloc->nloc*2, P(recon->anloc,idm));
			}
		}
	}
	MVR_MVM_T data={parms, recon, MVMt, curp, ntotact};
	int nthread;
	if(0){//handle one act at a time
		nthread=1;
		recon->nthread=NTHREAD;
	}else{//make each iact a single thread
		nthread=NTHREAD;//recon->nthread;
		recon->nthread=1;//make sure only 1 threaded is used per iact
		extern int fft_disable_threads;
		fft_disable_threads=1;//disable fftw threads in each thread.
	}
	thread_t *tdata=thread_prep(0, ntotact, nthread, setup_recon_mvr_mvm_iact, &data);
	CALL_THREAD(tdata, 1);
	free(tdata);
	dcell* MVM=dcelltrans(MVMt);
	dcellfree(MVMt);
	recon->MVM=dcell2m(MVM);
	dcellfree(MVM);
	free(curp);
}

/**
   Assemble matrix to do matrix vector multiply. Split from setup_recon because GPU may be used.
*/

void setup_recon_mvm(const parms_t* parms, recon_t* recon, const powfs_t* powfs){
	TIC;tic;
	if(!parms->recon.mvm) return;
	if(!recon->MVM){//GPU has not assembled the MVM.
		if(parms->recon.alg==RECON_MVR){
			setup_recon_mvr_mvm(recon, parms, powfs);
		} else{
			setup_recon_lsr_mvm(recon, parms, powfs);
		}
	}
	if(!parms->load.mvm&&(parms->save.setup||parms->save.mvm||parms->save.recon)){
		writebin(recon->MVM, "MVM.bin");
	}
	if(parms->sim.mvmport){
		mvm_client_init(parms->sim.mvmhost, parms->sim.mvmport, recon->MVM, parms->sim.mvmngpu);
	}
	toc2("setup_recon_mvm");
}
