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
#include "fdpcg.h"
#include "recon_utils.h"

/**
 * \file recon_tomo.c
 * @brief Setup tomography parameters
 * 
 */
 

/**
   Regularization terms for tomography. Can be called multiple times for Cn2 update.
*/
void
setup_recon_tomo_reg(recon_t* recon, const parms_t* parms){
	//info("setup_recon_tomo_reg\n");
	/*Free existing struct if already exist.  */
	free_recon_cxx(recon);
	if(parms->tomo.assemble){
	/*We need the old copy of L2 when we update the turbulence profile. */
		dspcellfree(recon->L2save);
		recon->L2save=recon->L2;
	} else{
		dspcellfree(recon->L2);
	}
	recon->L2=NULL;
	/*When layers get a weight less than 1%, we put it at 1% to avoid
	  regularization unstability issues.*/
	dclip(recon->wt, 0.01, 1);
	/*normalize the weights to sum to 1. */
	dnormalize_sumabs(recon->wt, 1);
	const int npsr=recon->npsr;
	recon->cxxalg=parms->tomo.cxxalg;
	/*test_cxx(recon, parms); */
	if(parms->tomo.cxxalg==0){
		if(parms->load.cxx){
			recon->L2=dspcellread("%s", parms->load.cxx);
			if(NX(recon->L2)!=npsr||NY(recon->L2)!=npsr){
				error("Wrong format of loaded L2\n");
			}
		} else{
			recon->L2=dspcellnew(npsr, npsr);
			for(int ips=0; ips<npsr; ips++){
				if(parms->tomo.square){/*periodic bc */
					P(recon->L2,ips,ips)=mklaplacian_map
					(P(recon->xmap,ips)->nx, P(recon->xmap,ips)->nx,
						P(recon->xloc,ips)->dx, recon->r0,
						P(recon->wt,ips));
				} else{/*reflecive bc */
					P(recon->L2,ips,ips)=mklaplacian_loc
					(P(recon->xloc,ips), recon->r0,	P(recon->wt,ips));
				}
			}
		}
		if(parms->save.setup>1){
			writebin(recon->L2, "L2");
		}
		dspcellscale(recon->L2, sqrt(parms->tomo.cxxscale*TOMOSCALE));
	}
	if(parms->tomo.cxxalg==1||(parms->tomo.cxxalg==2&&parms->tomo.precond==1)){
		recon->invpsd=mycalloc(1, invpsd_t);
		if(parms->load.cxx){
			recon->invpsd->invpsd=dcellread("%s", parms->load.cxx);
			if(NX(recon->invpsd->invpsd)!=npsr||NY(recon->invpsd->invpsd)!=1){
				error("Wrong format of loaded invpsd\n");
			}
		} else{
			dcell* invpsd=recon->invpsd->invpsd=dcellnew(npsr, 1);
			for(int ips=0; ips<npsr; ips++){
				long nx=P(recon->xmap,ips)->nx;
				long ny=P(recon->xmap,ips)->ny;
				real r0i=recon->r0*pow(P(recon->wt,ips), -3./5.);
				P(invpsd,ips)=turbpsd(nx, ny, P(recon->xloc,ips)->dx, r0i,
					recon->L0, 0, -1);
				dscale(P(invpsd,ips), pow((real)(nx*ny), -2));
			}
		}
		if(parms->save.setup){
			writebin(recon->invpsd->invpsd, "recon_invpsd");
		}
		dcellscale(recon->invpsd->invpsd, sqrt(parms->tomo.cxxscale*TOMOSCALE));

		ccell* fftxopd=recon->invpsd->fftxopd=ccellnew(recon->npsr, 1);
		for(int ips=0; ips<recon->npsr; ips++){
			P(fftxopd,ips)=cnew(P(recon->xmap,ips)->nx, P(recon->xmap,ips)->ny);
		}
		recon->invpsd->xloc=recon->xloc;
		recon->invpsd->square=parms->tomo.square;
	}
	if(parms->tomo.cxxalg==2){
		recon->fractal=mycalloc(1, fractal_t);
		recon->fractal->xloc=recon->xloc;
		recon->fractal->r0=parms->atmr.r0;
		recon->fractal->L0=parms->atmr.L0;
		recon->fractal->wt=P(parms->atmr.wt);
		recon->fractal->scale=sqrt(parms->tomo.cxxscale*TOMOSCALE);
		recon->fractal->ninit=parms->tomo.ninit;
		dcell* xopd=recon->fractal->xopd=dcellnew(npsr, 1);
		for(int ips=0; ips<npsr; ips++){
			int nn=nextfftsize(MAX(P(recon->xmap,ips)->nx, P(recon->xmap,ips)->ny))+1;
			P(xopd,ips)=dnew(nn, nn);
		}
	}

	if(parms->tomo.piston_cr){
	/*when add constraint, make sure the order of
	  magnitude are at the same range.*/
		dspcellfree(recon->ZZT);
		recon->ZZT=dspcellnew(npsr, npsr);
		for(int ips=0; ips<npsr; ips++){
			real r0=recon->r0;
			real dx=P(recon->xloc,ips)->dx;
			real wt=P(recon->wt,ips);
			real val=pow(laplacian_coef(r0, wt, dx), 2)*1e-6;
			/*dbg("Scaling of ZZT is %g\n",val); */
			/*piston mode eq 47 in Brent 2002 paper */
			int icenter=loccenter(P(recon->xloc,ips));
			int nloc=P(recon->xloc,ips)->nloc;
			dsp* ZZT=P(recon->ZZT, ips, ips)=dspnew(nloc, nloc, 1);
			int icol;
			int count=0;
			for(icol=0; icol<nloc; icol++){
				ZZT->pp[icol]=count;
				if(icol==icenter){
					ZZT->pi[count]=icenter;
					ZZT->px[count]=val;
					count++;
				}
			}
			ZZT->pp[nloc]=count;
		}
		if(parms->save.setup){
			writebin(recon->ZZT, "recon_ZZT");
		}
		dspcellscale(recon->ZZT, parms->tomo.cxxscale*TOMOSCALE);
	}
}

/**
   assemble tomography matrix. In CG mode, this function is not executed if
   tomo.assemble=0, Instead, the algorithm is contained in recon.c. When you
   modify anything, make sure you also do it there.

   For integrated tomograhy:

   \f$\hat{x}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1}+G_{ngs}^{T}C_{ngs}^{-1}
   G_{ngs})^{-1}(G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}+G_{ngs}^{T}C_{ngs}^{-1}s_{ngs}){\equiv}RL^{-1}RR s.\f$

   For split tomography, the terms regarding NGS are dropped:

   \f$\hat{x}_{lgs}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1})^{-1}
   G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}\equiv R_L^{-1} R_R s.\f$

   The left hand side of linear equation (inside the inverse) is stored in RL.
   The right hand side of the linear equation (outside of the inverse) is stored
   in RR. The terms regarding the NGS are handled using low rank terms. The
   gradients from LGS and NGS and concatenated to \f$s\f$.

   In the LGS part, there is a global tip/tilt removal operator because LGS is
   insensitive to global tip/tilts.

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803

*/
void setup_recon_tomo_matrix(recon_t* recon, const parms_t* parms){
	/*if not cg or forced, build explicitly the tomography matrix. */
	int npsr=recon->npsr;
	int nwfs=parms->nwfsr;
	/*Free OLD matrices if any. */
	muv_free(&recon->RR);
	muv_free(&recon->RL);
	print_mem("Before assembling tomo matrix");

	if(parms->load.tomo){
	/*right hand side. */
		warning("Loading saved recon->RR\n");
		recon->RR.M=readbin("RRM");
		if(zfexist("RRU")){
			recon->RR.U=dcellread("RRU");
			recon->RR.V=dcellread("RRV");
		}
		/*Left hand side */
		warning("Loading saved recon->RL\n");
		recon->RL.M=readbin("RLM");
		recon->RL.U=dcellread("RLU");
		recon->RL.V=dcellread("RLV");
		if(parms->tomo.alg==0&&zfexist("RLC")){
			recon->RL.C=chol_read("RLC");
		}
		if(parms->tomo.alg==2&&zfexist("RLMI")){
			recon->RL.MI=dread("RLMI");
		}
	} else{
		info("Building recon->RR\n");
		const dspcell* saneai=recon->saneai;
		/*
		  Reconstruction Right hand side matrix. In split tomography mode, low
		  order NGS are skipped. recon->GXtomo contains GXs that only
		  participate in tomography.
		*/
		dcellmm(&recon->RR.M, recon->GXtomo, saneai, "tn", 1);
		if(recon->TTF){
			//Tip/tilt and diff focus removal low rand terms for LGS WFS.
			dcellmm(&recon->RR.U, recon->RR.M, recon->TTF, "nn", 1);
			recon->RR.V=dcelltrans(recon->PTTF);
		}

		info("Building recon->RL\n"); /*left hand side matrix */
		cell *GtC_hi=cellnew(npsr, parms->nwfsr);
		//RL.M only contains high order WFS even in split=0 case.
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(!parms->powfs[ipowfs].lo && !parms->powfs[ipowfs].skip){
				for(int ipsr=0; ipsr<npsr; ipsr++){
					P(GtC_hi, ipsr, iwfs)=cellref(P(recon->RR.M, ipsr, iwfs));
				}
			}
		}
		dcellmm(&recon->RL.M, GtC_hi, recon->GXhi, "nn", 1);
		cellfree(GtC_hi);
		dspcell* RLM=dspcell_cast(recon->RL.M);
		if(parms->tomo.piston_cr){
			/*single point piston constraint. no need tikholnov.*/
			info("Adding ZZT to RLM\n");
			for(int ips=0; ips<npsr; ips++){
				dspadd(&P(RLM, ips, ips), 1, P(recon->ZZT,ips,ips), 1);
			}
			dspcellfree(recon->ZZT);
		}
		/*Apply tikholnov regularization.*/
		if(fabs(parms->tomo.tikcr)>1.e-15){
			/*Estimated from the Formula */
			real maxeig=pow(recon->neamhi*P(recon->xloc,0)->dx, -2);
			real tikcr=parms->tomo.tikcr;
			info("Adding tikhonov constraint of %.1e to RLM\n", tikcr);
			info("The maximum eigen value is estimated to be around %.1e\n", maxeig);
			dcelladdI(recon->RL.M, tikcr*maxeig);
		}
		/*add L2 and ZZT */
		switch(parms->tomo.cxxalg){
		case 0:/*Add L2'*L2 to RL.M */
			for(int ips=0; ips<npsr; ips++){
				dsp* tmp=dspmulsp(P(recon->L2,ips,ips), P(recon->L2,ips,ips), "tn");
				if(!tmp){
					error("L2 is empty!!\n");
				}
				dspadd(&P(RLM, ips, ips), 1, tmp, 1);
				dspfree(tmp);
			}
			break;
		case 1:/*Need to apply invpsd separately */
			recon->RL.extra=recon->invpsd;
			recon->RL.exfun=apply_invpsd;
			break;
		case 2:/*Need to apply fractal separately */
			recon->RL.extra=recon->fractal;
			recon->RL.exfun=apply_fractal;
		}
		/*Symmetricize, remove values below 1e-15*max and sort RLM (optional). */
		/*dspcellsym(recon->RL.M); */
		if(recon->TTF_LHS || recon->PTTF_LHS){
			dcellmm(&recon->RL.U, recon->RR.M, recon->TTF_LHS, "nn", 1);
			dcellmm(&recon->RL.V, recon->GX, recon->PTTF_LHS, "tt", 1);
		}
		if(!parms->recon.split){
			/*Left hand side low rank terms for low order wfs. Only in left hand side for integrated tomography. */
			dcell* ULo=dcellnew(npsr, nwfs);
			dcell* VLo=dcellnew(npsr, nwfs);
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				int ipowfs=parms->wfsr[iwfs].powfs;
				if(!parms->powfs[ipowfs].skip && parms->powfs[ipowfs].lo){
					for(int ips=0; ips<npsr; ips++){
						dspfull(&P(ULo, ips, iwfs), dsp_cast(P(recon->RR.M, ips, iwfs)), 'n', -1);
						dspfull(&P(VLo, ips, iwfs), P(recon->GX, iwfs, ips), 't', 1);
					}
				}
			}
			dcellcat2(&recon->RL.U, ULo, 2);
			dcellcat2(&recon->RL.V, VLo, 2);
		
			dcellfree(ULo);
			dcellfree(VLo);
		}
		/*Remove empty cells. */
		//dcelldropempty(&recon->RR.U,2);
		//dcelldropempty(&recon->RR.V,2);
		dcelldropempty(&recon->RL.U, 2);
		dcelldropempty(&recon->RL.V, 2);

		if(recon->RL.U){
			/* balance UV. may not be necessary. Just to compare well against laos. */
			real r0=recon->r0;
			real dx=P(recon->xloc,0)->dx;
			real val=laplacian_coef(r0, 1, dx);/*needs to be a constant */
			dcellscale(recon->RL.U, 1./val);
			dcellscale(recon->RL.V, val);
		}
		/*collect statistics.*/
		long nll=0, nlr=0;
		if(recon->RR.U){
			for(int i=0; i<NY(recon->RR.U);i++){
				if(P(recon->RR.U,0,i)){
					nlr+=P(recon->RR.U,0,i)->ny;
				}
			}
		}
		if(recon->RL.U){
			for(int i=0; i<NY(recon->RL.U);i++){
				if(P(recon->RL.U,0,i)){
					nll+=P(recon->RL.U,0,i)->ny;
				}
			}
		}
		info("Tomography number of Low rank terms: %ld in RHS, %ld in LHS\n", nlr, nll);
		if(parms->save.recon){
			writebin(recon->RR.M, "tomo_RRM");
			writebin(recon->RR.U, "tomo_RRU");
			writebin(recon->RR.V, "tomo_RRV");

			writebin(recon->RL.M, "tomo_RLM.bin");/*disable compression */
			writebin(recon->RL.U, "tomo_RLU");
			writebin(recon->RL.V, "tomo_RLV");
		}
			}
	if((parms->tomo.alg==0||parms->tomo.alg==2)&&parms->tomo.bgs){
	/* We need cholesky decomposition in CBS or MVST method. */
		muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	if(((parms->tomo.alg==0||parms->tomo.alg==2)&&!parms->tomo.bgs)||(parms->sim.ecnn&&!parms->recon.mvm)){
		if(parms->load.tomo){
			if(parms->tomo.alg==0&&zfexist("RLC")){
				recon->RL.C=chol_read("RLC");
			}
			if(parms->tomo.alg==2&&zfexist("RLMI")){
				recon->RL.MI=dread("RLMI");
			}
		}
		if(!recon->RL.C&&!recon->RL.MI){
			muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
		}
		print_mem("After cholesky/svd");
	}

	if(parms->save.recon){
		if(recon->RL.C){
		//chol_convert(recon->RL.C, 1);
			chol_save(recon->RL.C, "tomo_RLC.bin");
		}
		if(recon->RL.MI){
			writebin(recon->RL.MI, "tomo_RLMI");
		}
		if(recon->RL.Up){
			writebin(recon->RL.Up, "tomo_RLUp");
			writebin(recon->RL.Vp, "tomo_RLVp");
		}
		if(recon->RL.CB){
			for(int ib=0; ib<recon->RL.nb; ib++){
				chol_save(recon->RL.CB[ib], "recon_RLCB_%d.bin", ib);
			}
		}
		if(recon->RL.MIB){
			writebin(recon->RL.MIB, "tomo_RLMIB");
		}
	}
	/*Don't free PTT. Used in forming LGS uplink err */
	print_mem("After assemble tomo matrix");
}

/**
   Sets up the tomogrpahy turbulence reconstruction structs including wavefront
   reconstructor and DM fitting operator \callgraph

   Calls setup_recon_tomo_matrix()
   and setup_recon_fit_matrix() to setup the tomography and DM
   fitting matrix.

   AHST is handled in setup_ngsmod().

   MVST is handled in setup_recon_mvst().

   MOAO is handled in setup_recon_moao().

*/
void setup_recon_tomo(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	TIC;tic;
	/*setup inverse noise covariance matrix. */
	/*prepare for tomography setup */
	setup_recon_tomo_reg(recon, parms);
	if(!parms->recon.split || parms->tomo.splitlrt){
		//Tomography low rank term
		int full=(!parms->recon.split||parms->tomo.splitlrt==2);
		recon->TTF_LHS=dcellref(full?recon->TTF:recon->FF);
		recon->PTTF_LHS=dcellref(full?recon->PTTF:recon->PFF);
		if(parms->recon.split){
			//To avoid rank deficienciy in AHST, we do not remove the modes in the first wfs.
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				if(parms->powfs[ipowfs].nwfsr>0	&& !parms->powfs[ipowfs].skip){
					int iwfs=P(parms->powfs[ipowfs].wfsr, 0);
					if(recon->TTF_LHS) dfree(P(recon->TTF_LHS, iwfs, iwfs));
					if(recon->PTTF_LHS) dfree(P(recon->PTTF_LHS, iwfs, iwfs));
				}
			}
		}
	}else{
		info("No low rank terms in tomography left hand side.\n");
	}
	if(parms->tomo.assemble||parms->recon.split==2){
	/*assemble the matrix only if not using CG CG apply the
	  individual matrices on fly to speed up and save memory. */
		setup_recon_tomo_matrix(recon, parms);
	}

	if(parms->tomo.precond==1){
		fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
		recon->fdpcg=fdpcg_prepare(parms, recon, powfs, NULL);
	}

	/*Fall back function method if .M is NULL */
	recon->RL.Mfun=TomoL;
	recon->RL.Mdata=recon;
	recon->RR.Mfun=TomoR;
	recon->RR.Mtfun=TomoRt;
	recon->RR.Mdata=recon;
	if(parms->tomo.alg==1){/*CG */
		switch(parms->tomo.precond){
		case 0:/*no preconditioner */
			recon->RL.pfun=NULL;
			recon->RL.pdata=NULL;
			break;
		case 1:
			recon->RL.pfun=fdpcg_precond;
			recon->RL.pdata=recon;
			break;
		default:
			error("Invalid tomo.precond");
		}
	}
	recon->RL.alg=parms->tomo.alg;
	recon->RL.bgs=parms->tomo.bgs;
	recon->RL.warm=parms->tomo.cgwarm;
	recon->RL.maxit=parms->tomo.maxit;

	toc2("setup_recon_tomo");
}

#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_tm
#define tic_tm
#define toc_tm(A)
#else
#define TIC_tm TIC
#define tic_tm tic
#define toc_tm(A) toc2(A);tic
#endif

typedef struct Tomo_T{
	const recon_t* recon;
	const real alpha;
	const dcell* xin;/*input */
	dcell* gg;/*intermediate gradient */
	dcell* xout;/*output */
	int isR;/*Mark that we are right hand side.*/
	unsigned int ips;//counter for CALL
	unsigned int iwfs;//counter for CALL
}Tomo_T;

/**
   Speed up TomoL by gathering the first part of operations (GP*HXW) belonging
   to each WFS to facilitate threading.

   gg  = GP * HXW * xin
*/
static void* Tomo_prop_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t* parms=global->parms;
	sim_t* simu=global->simu;
	const int nps=recon->npsr;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(parms->powfs[ipowfs].skip) continue;
		const int shwfs=parms->powfs[ipowfs].type==WFS_SH;
		const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
		dmat* xx=dnew(recon->ploc->nloc, 1);
		const real hs=parms->wfs[iwfs].hs;
		for(int ips=0; ips<nps; ips++){
			if(parms->tomo.square&&!parms->dbg.tomo_hxw){
			/*Do the ray tracing instead of using HXW. */
				real misregx=(shwfs?parms->wfsr[iwfs].misregx:0);
				real misregy=(shwfs?parms->wfsr[iwfs].misregy:0);
				real shiftx=0;
				real shifty=0;
				if(parms->tomo.predict){
					int ips0=P(parms->atmr.indps, ips);
					shiftx=P(simu->atm, ips0)->vx*delay;
					shifty=P(simu->atm, ips0)->vy*delay;
				}
				prop(&(propdata_t){.mapin=P(recon->xmap, ips), .phiin=P(P(data->xin, ips)),
					.ostat=recon->ploc->stat, .phiout=P(xx), .alpha=1, .hs=hs,
					.thetax=parms->wfsr[iwfs].thetax, .thetay=parms->wfsr[iwfs].thetay,
					.misregx=misregx, .misregy=misregy, .shiftx=shiftx, .shifty=shifty});
			} else{
				dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
				dspmm(&xx, P(HXW, iwfs, ips), P(data->xin, ips), "nn", 1);
			}
		}
		dmat *xxr=NULL;//for rotation misregistration
		if(recon->wfsr[iwfs].misregc){
			if(recon->pmap->nx*recon->pmap->ny==recon->ploc->nloc){//ploc is quare
				xxr=dnew(recon->pmap->nx, recon->pmap->ny);//rotation
				reshape(xx, xxr->nx, xxr->ny);	
				dembed(xxr, xx, -recon->wfsr[iwfs].misregc);
				reshape(xxr, recon->ploc->nloc, 1);
			}else{
				xxr=dnew(recon->ploc->nloc, 1);
				prop(&(propdata_t){.locin=recon->ploc, .phiin=P(xx), .locout=recon->ploc, .phiout=P(xxr), 
					.rot=-recon->wfsr[iwfs].misregc});
			}
		}
		/*Apply the gradient operation */
		dspmm(&P(data->gg, iwfs), P(recon->GP, iwfs), xxr?xxr:xx, "nn", 1);
		dfree(xx);
		dfree(xxr);
		/* For each wfs, Ray tracing takes 1.5 ms.  GP takes 0.7 ms. */
	}
	return NULL;
}
/**
   Wrapper of Tomo_prop_do
*/
static void Tomo_prop(Tomo_T* data, int nthread){
	//thread_t info[nthread];
	//thread_prep(info, 0, NX(data->gg), nthread, Tomo_prop_do, data);
	//CALL_THREAD(info, 1);
	data->iwfs=0;
	nthread=1;
	CALL(Tomo_prop_do, data, nthread, 1);
}
/**
   Speed up TomoL by gathering the second part of operations (GP') belonging to
   each WFS to facilitate threading.

   gg = GP' * NEAI * gg;
*/
static void* Tomo_nea_gpt_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t *parms=global->parms;
	dspcell* NEAI=recon->saneai/*PDSPCELL*/;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		dmat* gg2=NULL;
		/*Apply the gradient operation */
		dspmm(&gg2, P(NEAI, iwfs, iwfs), P(data->gg, iwfs), "nn", 1);
		dfree(P(data->gg, iwfs)); /*We reuse gg. */
		dmat *ggr=NULL;//intermediate for rotation
		dmat **pgg=recon->wfsr[iwfs].misregc?&ggr:&P(data->gg, iwfs);
		dspmm(pgg, P(recon->GP, iwfs), gg2, "tn", data->alpha);
		dfree(gg2);
		if(recon->wfsr[iwfs].misregc){
			if(recon->pmap->nx*recon->pmap->ny==recon->ploc->nloc){//ploc is quare
				reshape(ggr, recon->pmap->nx, recon->pmap->ny);
				P(data->gg, iwfs)=dnew(ggr->nx, ggr->ny);
				dembed(P(data->gg, iwfs), ggr, recon->wfsr[iwfs].misregc);
				reshape(P(data->gg, iwfs), recon->ploc->nloc, 1);
			}else{
				P(data->gg, iwfs)=dnew(recon->ploc->nloc, 1);
				prop(&(propdata_t){.locin=recon->ploc, .phiin=P(ggr), .locout=recon->ploc, .phiout=P(P(data->gg, iwfs)), 
					.rot=recon->wfsr[iwfs].misregc});
			}
			dfree(ggr);
		}
	}
	return NULL;
}
/**
 * gg=NEAI*gg
 */
static void* Tomo_nea_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t *parms=global->parms;
	dspcell* NEAI=recon->saneai/*PDSPCELL*/;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		dmat* gg2=NULL;
		/*Apply the gradient operation */
		dspmm(&gg2, P(NEAI, iwfs, iwfs), P(data->gg, iwfs), "nn", 1);
		dcp(&P(data->gg, iwfs), gg2);
		dfree(gg2);
	}
	return NULL;
}

/**
 * Apply NEA weighting with optional Gp'
 * gg = GP' * NEAI * gg;
 */
void Tomo_nea_gpt(Tomo_T* data, int nthread, int gpt){
	
	data->iwfs=0;
	nthread=1;
	/*Wrapp of Tomo_nea_do for multi-threads*/
	if(gpt){
		CALL(Tomo_nea_gpt_do, data, nthread, 1);
	}else{
		CALL(Tomo_nea_do, data, nthread, 1);
	}
}
/**
   Speed up TomoL by gathering the third part of operations (GP') belonging to
   each WFS to facilitate threading. gg->xout.

   xout = Cxx^-1 * xin + HXW' * gg;
*/
static void* Tomo_iprop_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t* parms=global->parms;
	sim_t* simu=global->simu;
	map_t xmap;
	int ips;
	while((ips=atomic_fetch_add(&data->ips, 1))<NX(data->xout)){
		if(parms->tomo.square&&!parms->dbg.tomo_hxw){
			/*Do the ray tracing instead of using HXW. */
			if(!P(data->xout, ips)){
				P(data->xout, ips)=dnew(P(recon->xloc, ips)->nloc, 1);
			}
			memcpy(&xmap, P(recon->xmap, ips), sizeof(map_t));
			xmap.p=P(P(data->xout, ips));
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				if(!P(data->gg, iwfs)) continue;
				const real hs=parms->wfs[iwfs].hs;
				const int ipowfs=parms->wfs[iwfs].powfs;
				const int shwfs=parms->powfs[ipowfs].type==WFS_SH;
				real misregx=(shwfs?parms->wfsr[iwfs].misregx:0);
				real misregy=(shwfs?parms->wfsr[iwfs].misregy:0);
				real shiftx=0, shifty=0;
				if(parms->tomo.predict){
					const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
					int ips0=P(parms->atmr.indps, ips);
					shiftx=P(simu->atm, ips0)->vx*delay;
					shifty=P(simu->atm, ips0)->vy*delay;
				}
				prop(&(propdata_t){.mapin=&xmap, .ostat=recon->ploc->stat, .phiout=P(P(data->gg, iwfs)), 
					.alpha=1, .hs=hs, .misregx=misregx, .misregy=misregy,
					.shiftx=shiftx, .shifty=shifty});
			}
		} else{
			dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				dspmm(&P(data->xout, ips), P(HXW, iwfs, ips), P(data->gg, iwfs), "tn", 1);
			}
		}
		if(data->xin){/*data->xin is empty when called from TomoR */
			switch(recon->cxxalg){
			case 0:{/*L2 */
				dmat* xx=NULL;
				dspmm(&xx, P(recon->L2, ips, ips), P(data->xin, ips), "nn", 1);
				dspmm(&P(data->xout, ips), P(recon->L2, ips, ips), xx, "tn", data->alpha);
				dfree(xx);
			}
				break;
			case 1:
				apply_invpsd(&data->xout, recon->invpsd, data->xin, data->alpha, ips, ips);
				break;
			case 2:
				apply_fractal(&data->xout, recon->fractal, data->xin, data->alpha, ips, ips);
				break;
			}
			if(recon->ZZT){
				dspmm(&P(data->xout, ips), P(recon->ZZT, ips, ips), P(data->xin, ips), "tn", data->alpha);
			}
		}
	}
	return NULL;
}
/**
   Wrapper of Tomo_iprop_do
 */
static void Tomo_iprop(Tomo_T* data, int nthread){
	data->ips=0;
	nthread=1;
	CALL(Tomo_iprop_do, data, nthread, 1);
}
/**
   Apply tomography right hand operator without using assembled matrix. Fast and
   saves memory.  The operation is the same as the Tomo_nea_gpt and Tomo_iprop in
   TomoL, so merge the implemenations.

   xout=HXW'*GP'*NEAI*(1-TTF*PTTF)*gin.
*/

void TomoR(dcell** xout, const void* A,
	const dcell* gin, const real alpha){
	TIC_tm;tic_tm;
	const recon_t* recon=(const recon_t*)A;
	dcell* gg=NULL;
	dcellcp(&gg, gin);/*copy to gg so we don't touch the input. */
	remove_mode(gg, recon->TTF, recon->PTTF);
	if(!*xout){
		*xout=dcellnew(recon->npsr, 1);
	}
	Tomo_T data={recon, alpha, NULL, gg, *xout, 1};
	Tomo_nea_gpt(&data, recon->nthread, 1);
	Tomo_iprop(&data, recon->nthread);
	dcellfree(gg);
	toc_tm("TomoR");
}

/**
   Transpose operation of TomoRt. From xout -> gin
 */
void TomoRt(dcell** gout, const void* A,
	const dcell* xin, const real alpha){
	const recon_t* recon=(const recon_t*)A;
	if(!*gout){
		*gout=dcellnew(NX(recon->saneai), 1);
	}
	Tomo_T data={recon, alpha, xin, *gout, NULL, 1};
	Tomo_prop(&data, recon->nthread);
	/*Using Tomo_nea_gpt followed by TTFRt is equilvaent as using remove_mode followed by Tomo_nea_gpt.*/
	remove_mode(*gout, recon->TTF, recon->PTTF);
	Tomo_nea_gpt(&data, recon->nthread, 0);
}

/**
   Apply tomography left hand side operator without using assembled matrix. Fast
   and saves memory. Only useful in CG. Accumulates to xout;

   Some timing information:
   Generated in T410s with Intel Core i5 520 M @ 2.4 GHz (1 T is 1 thread, 2 T is 2 thread)
   and Poweredge with dual Xeon X5355 2.66Ghz (1 P is 1 thread, 8 P is 8 thread)
   Time  1 T  2 T  1 P  8 P
   HXW:  7.8  5.4  7.7  5.0
   GP:   4.3  2.6  3.6  1.4
   remove_mode: 0.3  0.3  0.4  0.6
   neai: 0.5  0.3  0.5  0.4
   GP':  3.2  2.1  3.4  1.0
   HXW': 5.5  4.7  7.4  4.9
   cxx:  3.2  2.5  3.5  2.7
*/

void TomoL(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	TIC_tm;tic_tm;
	const recon_t* recon=(const recon_t*)A;
	const parms_t* parms=global->parms;
	assert(NY(xin)==1);/*modify the code for ny>1 case. */
	dcell* gg=dcellnew(parms->nwfsr, 1);
	if(!*xout){
		*xout=dcellnew(recon->npsr, 1);
	}
	Tomo_T data={recon, alpha, xin, gg, *xout, 0};

	Tomo_prop(&data, recon->nthread);

	/*if(parms->recon.split!=1||parms->tomo.splitlrt==2){
		//Remove global Tip/Tilt, focus only in integrated tomography.
		remove_mode(gg, recon->TTF, recon->PTTF);
	} else if(parms->tomo.splitlrt){
		//Remove only focus in split tomography when splitlrt is set
		remove_mode(gg, recon->FF, recon->PFF);
	}*/
	remove_mode(gg, recon->TTF_LHS, recon->PTTF_LHS);
	Tomo_nea_gpt(&data, recon->nthread, 1);
	Tomo_iprop(&data, recon->nthread);
	/*
	   square=1  square=0 (1 thread on T410s)
	   prop:  takes 6 ms 13 ms
	   nea:   takes 4 ms 4 ms
	   iprop: takes 6 ms 9 ms
	*/
	dcellfree(gg);
	/*Tikhonov regularization is not added because it is not necessary in CG mode.*/
	toc_tm("TomoL");
}
