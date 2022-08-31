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
#include "recon.h"
#include "fdpcg.h"
#include "ahst.h"
#include "recon_utils.h"
#include "moao.h"
#include "powfs.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file recon_setup.c

   Contains routines that setup the wavefront reconstructor and DM fitting.

*/
/*
   TOMOSCALE is used for RLM, RRM, in MVST for M, and in CUDA due to
   limited dynamic range of single precision floating point numbers.

   Use parms->wfsr instead of parms->wfs for wfs information,
   which hands GLAO mode correctly.

   2014-04-01: Scale saneai, cxx, zzt by TOMOSCALE.

   All routines in this file depends on saneai and may be called repeatedly during simulation.
*/


/**
   Setup the matrix of the inverse of gradient measurement noise equivalent
   angle covariance matrix. For physical optics wfs, the NEA is computed using
   matched filter output. For geometric optics, the NEA is from input.
*/
static void
setup_recon_saneai(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	const int nwfs=parms->nwfsr;
	dspcellfree(recon->sanea);
	dspcellfree(recon->saneai);
	dspcellfree(recon->saneal);
	dspcell* sanea=recon->sanea=dspcellnew(nwfs, nwfs);//The subaperture NEA
	dspcell* saneal=recon->saneal=dspcellnew(nwfs, nwfs);//The LL' decomposition of sanea
	dspcell* saneai=recon->saneai=dspcellnew(nwfs, nwfs);//The inverse of sanea
	dfree(recon->neam);
	recon->neam=dnew(parms->nwfsr, 1);
	real neam_hi=0;
	int count_hi=0;
	info2("Recon NEA: ");
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].skip==3) continue;
		const int nsa=powfs[ipowfs].saloc->nloc;
		const int ng=parms->powfs[ipowfs].ng;
		int ncol=MAX(ng,3);
		const real pixtheta=parms->powfs[ipowfs].pixtheta;
		dcell* saneac=0;//in unit of rad^2.

		if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy)&&!parms->powfs[ipowfs].phyusenea){
			/*Physical optics use nea from intstat*/
			if(parms->recon.glao && PN(powfs[ipowfs].sanea)>1){//average sanea 
				info("Averaging saneaxy of different WFS for GLAO mode\n");
				saneac=dcellnew(1,1);
				int ni0=PN(powfs[ipowfs].sanea);
				real scale=1./ni0;
				for(int ii0=0; ii0<ni0; ii0++){
					dadd(&P(saneac,0), 1, P(powfs[ipowfs].sanea, ii0), scale);
				}
			}else{
				saneac=dcelldup(powfs[ipowfs].sanea);
			}
		} else{
			if(parms->powfs[ipowfs].neareconfile){
				saneac=dcellread_prefix(parms->powfs[ipowfs].neareconfile, parms, ipowfs);
				for(int i=0; i<NX(saneac)*NY(saneac); i++){
					nea_check(P(saneac,i), nsa, ng);
					nea_mm(&P(saneac,i), P(saneac,i), ng);
				}
			} else{
				saneac=dcellnew(1, 1);
				P(saneac,0)=dnew(nsa, ncol);
				real neamas=parms->powfs[ipowfs].nearecon;
				if(neamas<0.001||neamas > 2000){
					warning("powfs[%d].nearecon=%g mas may have unit incorrect.\n", ipowfs, neamas);
				}
				//convert from mill-arcsec to radian.
				real nearad=pow(neamas/206265000., 2)/(parms->powfs[ipowfs].dtrat);
				for(int isa=0; isa<nsa; isa++){
					P(P(saneac,0), isa, 0)=P(P(saneac,0), isa, 1)=nearad/(P(powfs[ipowfs].saa, isa));
					for(int ig=2; ig<ng; ig++){
						P(P(saneac, 0), isa, ig)=P(P(saneac, 0), isa, 0);
					}
				}
			}
		}
		int do_ref=0;
		if(NX(saneac)==1){
			do_ref=1;
		} else if(parms->recon.glao){
			error("Please average nearecon for GLAO mode.\n");
		}

		const real area_thres=(nsa>4)?0.9*parms->powfs[ipowfs].safill2d:0;
		const real neaextra2=copysign(pow(parms->powfs[ipowfs].neaextra/206265000., 2),
			parms->powfs[ipowfs].neaextra);
		const real neamin2=pow(parms->powfs[ipowfs].neamin/206265000., 2);

		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfsr; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfsr,jwfs);
			int iwfs0=P(parms->powfs[ipowfs].wfsr,0);
			lmat* samask=0;
			if(parms->wfs[iwfs].sabad){
				samask=loc_coord2ind(powfs[ipowfs].saloc, parms->wfs[iwfs].sabad);
			}
			if(!do_ref||iwfs==iwfs0||parms->wfs[iwfs].sabad||parms->wfs[iwfs0].sabad){
				dmat* sanea0=PR(saneac, jwfs, 0);
				dmat* sanea0l=0;
				dmat* sanea0i=0;

				if(!parms->powfs[ipowfs].mtchcpl && NY(sanea0)!=2 && ng==2){
					sanea0->ny=2; //reduce coupling
				}
				real nea2_sum=0;
				long nea2_count=0;
				for(int isa=0; isa<nsa; isa++){
					if(samask && P(samask, isa)){
						warning("wfs %d sa %d is masked\n", iwfs, isa);
						for(int iy=0; iy<ng; iy++){
							P(sanea0, isa, iy)=INFINITY;
						}
					}
					for(int iy=0; iy<ng; iy++){
						if((P(sanea0, isa, iy)+=neaextra2)<neamin2){
							P(sanea0, isa, iy)=neamin2;
						}
					}
					if(P(powfs[ipowfs].saa, isa)>area_thres){
						nea2_sum+=P(sanea0, isa, 0)+P(sanea0, isa, 1);
						nea2_count++;
					}
				}
				
				nea_chol(&sanea0l, sanea0, ng);
				nea_inv(&sanea0i, sanea0, ng, TOMOSCALE);//without TOMOSCALE, it overflows in float mode
				
				real nea_mean=sqrt(nea2_sum/nea2_count*0.5);
				P(recon->neam,iwfs)=nea_mean/(parms->powfs[ipowfs].skip?1:sqrt(TOMOSCALE));
				if(nea_mean>pixtheta*0.33
					&&parms->powfs[ipowfs].usephy
					&&parms->powfs[ipowfs].order<=2
					&&parms->sim.dtrat_lo==parms->sim.dtrat_lo2
					){
					warning("TT WFS %d has too much measurement error: %g mas\". Ignore it\n",
						iwfs, nea_mean*206265000);
					P(sanea,iwfs,iwfs)=dspnewdiag(nsa*ng, NULL, INFINITY);
					P(saneal,iwfs,iwfs)=dspnewdiag(nsa*ng, NULL, 0);
					P(saneai,iwfs,iwfs)=dspnewdiag(nsa*ng, NULL, 0);
				} else{
					P(sanea,iwfs,iwfs)=nea2sp(sanea0, 1, 1, ng);
					P(saneal,iwfs,iwfs)=nea2sp(sanea0l, 1, 0, ng);
					P(saneai,iwfs,iwfs)=nea2sp(sanea0i, 1, 1, ng);
				}
				dfree(sanea0l);
				dfree(sanea0i);
			} else if(do_ref){
				P(sanea, iwfs, iwfs)=dspref(P(sanea, iwfs0, iwfs0));
				P(saneal, iwfs, iwfs)=dspref(P(saneal, iwfs0, iwfs0));
				P(saneai, iwfs, iwfs)=dspref(P(saneai, iwfs0, iwfs0));
				P(recon->neam,iwfs)=P(recon->neam,iwfs0);
			}
			lfree(samask);

			if(!parms->powfs[ipowfs].lo){
				neam_hi+=pow(P(recon->neam,iwfs), 2);
				count_hi++;
			}
		}/*iwfs*/
		dcellfree(saneac);
	}/*ipowfs */

	recon->neamhi=sqrt(neam_hi/count_hi);
	
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].skip==3) continue;
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfsr; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfsr,jwfs);
			const char* neatype;

			if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy)&&
				!parms->powfs[ipowfs].phyusenea){
				switch(parms->powfs[ipowfs].phytype_recon){
				case 1:
					neatype=parms->powfs[ipowfs].mtchcr?"cmf":"mf"; break;
				case 2:
					neatype="cog";break;
				default:
					neatype="unknown";break;
				}
			} else if(parms->powfs[ipowfs].neareconfile){
				neatype="file";
			} else{
				neatype="geom";
			}
			info2("%s(%.2f) ", neatype, P(recon->neam,iwfs)*206265000*(parms->powfs[ipowfs].skip?1:sqrt(TOMOSCALE)));
		}
	}
	info2(" mas\n");
	if(parms->save.setup){
		writebin(recon->sanea, "sanea");
		writebin(recon->saneai, "saneai");
		writebin(recon->saneal, "saneal");
	}
}

/**
   wrapps setup_recon_TTR() and setup_recon_DFR() to removal global tip/tilt and
   differential focus.
*/
static void
setup_recon_TTFR(recon_t* recon, const parms_t* parms){
	cellfree(recon->PTT);
	cellfree(recon->PFF);
	cellfree(recon->PTTF);

	recon->PTT=dcellpinv(recon->TT, recon->saneai);
	recon->PFF=dcellpinv(recon->FF, recon->saneai);
	recon->PTTF=dcellpinv(recon->TTF, recon->saneai);
	if(parms->save.setup){
		writebin(recon->TTF, "TTF");
		writebin(recon->PTT, "PTT");
		writebin(recon->PTTF, "PTTF");
	}
	/*dcellfree(recon->DF);//don't free DF to use in PDF. */
	/*Keep TT, PTT, used in fsm pointing or dithering. */
}
/**
   Frees recon->invpsd or recon->fractal
*/
static void free_cxx(recon_t* recon){
	if(recon->invpsd){
		dcellfree(recon->invpsd->invpsd);
		ccellfree(recon->invpsd->fftxopd);
		free(recon->invpsd);
		recon->invpsd=NULL;
	}
	if(recon->fractal){
		dcellfree(recon->fractal->xopd);
		free(recon->fractal);
		recon->fractal=NULL;
	}
}
/**
   Prepares for tomography. ALlow it to be called multiple times for Cn2 update.
*/
void
setup_recon_tomo_prep(recon_t* recon, const parms_t* parms){
	//info("setup_recon_tomo_prep\n");
	/*Free existing struct if already exist.  */
	free_cxx(recon);
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
	dnormalize_sumabs(P(recon->wt), recon->npsr, 1);
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
		if(parms->save.setup){
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
		dspcell* GX=recon->GX/*PDSPCELL*/;
		const dspcell* saneai=recon->saneai;
		/*
		  Reconstruction Right hand side matrix. In split tomography mode, low
		  order NGS are skipped. recon->GXtomo contains GXs that only
		  participate in tomography.
		*/
		dspcell* GXtomoT=dspcelltrans(recon->GXtomo);
		dcellmm_any(&recon->RR.M,CELL(GXtomoT), CELL(saneai), "nn", 1);
		dspcell* RRM=(dspcell*)recon->RR.M/*PDSPCELL*/;
		/*
		  Tip/tilt and diff focus removal low rand terms for LGS WFS.
		*/
		if(recon->TTF){
			dcellmm_cell(&recon->RR.U, recon->RR.M, recon->TTF, "nn", 1);
			recon->RR.V=dcelltrans(recon->PTTF);
		}

		info("Building recon->RL\n"); /*left hand side matrix */
		dcellmm_any(&recon->RL.M, recon->RR.M, CELL(recon->GXtomo), "nn", 1);
		dspcell* RLM=(dspcell*)recon->RL.M/*PDSPCELL*/;
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
			dcelladdI_any(recon->RL.M, tikcr*maxeig);
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

		/*Low rank terms for low order wfs. Only in Integrated tomography. */
		dcell* ULo=dcellnew(npsr, nwfs);
		dcell* pULo=ULo/*PDELL*/;
		dcell* VLo=dcellnew(npsr, nwfs);
		dcell* pVLo=VLo/*PDELL*/;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip){
				continue;
			}
			if(parms->powfs[ipowfs].lo){
				for(int ips=0; ips<npsr; ips++){
					dspfull(&P(pULo, ips, iwfs), P(RRM, ips, iwfs), 'n', -1);
					dspfull(&P(pVLo, ips, iwfs), P(GX, iwfs, ips), 't', 1);
				}
			}
		}
		if(parms->recon.split!=1||parms->tomo.splitlrt){
			recon->RL.U=dcellcat(recon->RR.U, ULo, 2);
			dcell* GPTTDF=NULL;
			dspcellmm(&GPTTDF, recon->GX, recon->RR.V, "tn", 1);
			recon->RL.V=dcellcat(GPTTDF, VLo, 2);
			dcellfree(GPTTDF);
		} else{
			info("Skipping RL Low rank terms in split tomography\n");
		}
		dcellfree(ULo);
		dcellfree(VLo);
		/*Remove empty cells. */
		//dcelldropempty(&recon->RR.U,2);
		//dcelldropempty(&recon->RR.V,2);
		dcelldropempty(&recon->RL.U, 2);
		dcelldropempty(&recon->RL.V, 2);

		if(recon->RL.U){
			/* balance UV. may not be necessary. Just to compare well against
			   laos. */
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
			writecell(recon->RR.M, "tomo_RRM");
			writebin(recon->RR.U, "tomo_RRU");
			writebin(recon->RR.V, "tomo_RRV");

			writecell(recon->RL.M, "tomo_RLM.bin");/*disable compression */
			writebin(recon->RL.U, "tomo_RLU");
			writebin(recon->RL.V, "tomo_RLV");
		}
		dspcellfree(GXtomoT);
	}
	if((parms->tomo.alg==0||parms->tomo.alg==2)&&parms->tomo.bgs){
	/* We need cholesky decomposition in CBS or MVST method. */
		muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	if(((parms->tomo.alg==0||parms->tomo.alg==2)&&!parms->tomo.bgs)||parms->sim.ecnn){
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

static dcell* setup_recon_ecnn(recon_t* recon, const parms_t* parms, loc_t* locs, lmat* mask){
	/**
	   We compute the wavefront estimation error covariance in science focal
	   plane due to wavefront measurement noise. Basically we compute
	   Hx*E*Cnn*E'*Hx' where E is the tomography operator, and Hx is ray
	   tracing from tomography grid xloc to science focal plane ploc. Since
	   Cnn is symmetrical and sparse, we can decompose it easily into
	   Cnn=Cnl*Cnl'; We first compute L=Hx*E*Cnl, and the result is simply
	   LL'; This is much faster than computing left and right separately,
	   because 1) the number of points in xloc is larger than in Cnn, so
	   after the tomography right hand side vector is applied, the number of
	   rows is larger than number of columns, this causes the right hand
	   side solver to be much slower. 2) Simply real the computation.

	   For HX opeation, build the sparse matrix and do multiply is way
	   slower than doing ray tracing directly.

	   For ad hoc split tomography, we need to remove the five NGS modes
	   from here, as well as in time averaging of estimated turbulence.

	   recon->saneal contains Cnl.
	*/
	TIC;tic;
	read_self_cpu();
	dmat* t1=NULL;
	if(recon->MVM){//MVM
		dspcell* sanealhi=dspcellnew(parms->nwfsr, parms->nwfsr);
		for(int iwfsr=0; iwfsr<parms->nwfsr; iwfsr++){
			int ipowfs=parms->wfsr[iwfsr].powfs;
			if(!parms->powfs[ipowfs].skip){
				P(sanealhi, iwfsr, iwfsr)=dspref(P(recon->saneal, iwfsr, iwfsr));
			}
		}
		//dsp *tmp=dspcell2sp(sanealhi);  dspcellfree(sanealhi);
		dmat* tmp=dspcell2m(sanealhi); dspcellfree(sanealhi);
		dmm(&t1, 0, recon->MVM, tmp, "nn", 1);
		cellfree(tmp);
		toc2("MVM ");tic;
	} else if(parms->recon.alg==0){//MV
		dcell* tmp2=NULL;
		dcell* tmp=NULL;
		dspcellfull(&tmp2, recon->saneal, 'n', 1);
		muv(&tmp, &recon->RR, tmp2, 1);
		toc2("RR ");tic;
		cellfree(tmp2);
		muv_direct_solve(&tmp2, &recon->RL, tmp); dcellfree(tmp);
		toc2("RL ");tic;
		//2015-06-29: Put in ommited DM fitting operation
		muv(&tmp, &recon->fit->FR, tmp2, 1); dcellfree(tmp2);
		toc2("FR ");tic;
		muv_direct_solve(&tmp2, &recon->fit->FL, tmp); dcellfree(tmp);
		toc2("FL ");tic;
		t1=dcell2m(tmp2); dcellfree(tmp2);
	} else{//LSR
		dcell* tmp2=NULL;
		dcell* tmp=NULL;
		dspcellfull(&tmp2, recon->saneal, 'n', 1);
		muv(&tmp, &recon->LR, tmp2, 1); cellfree(tmp2);
		toc2("LR ");tic;
		muv_direct_solve(&tmp2, &recon->LL, tmp); dcellfree(tmp);
		toc2("LL ");tic;
		t1=dcell2m(tmp2); dcellfree(tmp2);
	}
	dcell* ecnn=dcellnew(parms->evl.nevl, 1);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(mask&&!P(mask,ievl)) continue;
		tic;
		/*Build HX for science directions that need ecov.*/
		dmat* x1=dnew(locs->nloc, NY(t1));
		real hs=P(parms->evl.hs,ievl);
		int offset=0;
		for(int idm=0; idm<parms->ndm; idm++){
			const real ht=parms->dm[idm].ht;
			const real scale=1.-ht/hs;
			const real dispx=P(parms->evl.thetax,ievl)*ht;
			const real dispy=P(parms->evl.thetay,ievl)*ht;
			for(int icol=0; icol<NY(t1); icol++){
				prop_nongrid(P(recon->aloc,idm), PCOL(t1, icol)+offset,
					locs, PCOL(x1, icol), 1, dispx, dispy, scale, 0, 0);
			}
			offset+=P(recon->aloc,idm)->nloc;
		}
		toc2("Prop ");tic;
		dmm(&P(ecnn,ievl), 0, x1, x1, "nt", 1);
		dfree(x1);
		toc2("MM ");
	}
	dcellfree(t1);
	return ecnn;
}
/**
   Update assembled tomography matrix with new L2. Called from cn2est when new
   profiles are available.
*/
void setup_recon_update_cn2(recon_t* recon, const parms_t* parms){
	if(parms->sim.evlol) return;
	setup_recon_tomo_prep(recon, parms); /*redo L2, invpsd */
#if USE_CUDA
	if(parms->gpu.tomo){
		gpu_update_recon_cn2(parms, recon);
	}
#endif
	if(parms->tomo.alg==1&&!parms->tomo.assemble){/*no need to do anything */
		return;
	}
	if(parms->tomo.cxxalg==0&&recon->L2save){
	/*Need to adjust RLM with the new L2. */
		dspcell* RLM=(dspcell*)recon->RL.M/*PDSPCELL*/;
		const int npsr=recon->npsr;
		for(int ips=0; ips<npsr; ips++){
			dsp* LL=dspmulsp(P(recon->L2,ips,ips),
				P(recon->L2,ips,ips), "tn");
			dsp* LLold=dspmulsp(P(recon->L2save,ips,ips),
				P(recon->L2save,ips,ips), "tn");
			if(!LL){
				error("L2 is empty!!\n");
			}
			dsp* LLdiff=dspadd2(LL, 1, LLold, -1);/*adjustment to RLM */
			dspadd(&P(RLM, ips, ips), 1, LLdiff, 1);
			dspfree(LLdiff);
			dspfree(LL);
			dspfree(LLold);
		}
	}

	if(parms->tomo.alg==0||parms->tomo.alg==2){
	/*We need cholesky decomposition in CBS or MVST method. */
		if(!parms->tomo.bgs){/*Full Matrix */
			muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
		} else{/*BGS */
			muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
		}
		print_mem("After cholesky/svd");
	}

	if(parms->recon.split==2){
		setup_recon_mvst(recon, parms);
	}
#if USE_CUDA
	if(parms->gpu.tomo){
		gpu_update_recon_control(parms, recon);
	}
#endif
}


/**
   Create the reconstructor to reconstruct the residual focus error due to LGS
   sodium tracking error. Need to average out the focus error caused by
   atmosphere when applying (a low pass filter is applied to the output).  */
static void
setup_recon_focus(recon_t* recon, const parms_t* parms){
	if(parms->nlgspowfs){
		if(parms->recon.split==2&&parms->sim.mffocus){//For MVST.
			dmat* GMGngs=NULL;
			dcell* GMngs=dcellnew(1, parms->nwfs);
			/*Compute focus reconstructor from NGS Grads. fuse grads
			  together to construct a single focus measurement*/
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				if(parms->powfs[ipowfs].trs==0&&parms->powfs[ipowfs].order>1&&parms->powfs[ipowfs].skip!=2){
					info("wfs %d will be used to track focus\n", iwfs);
				} else{
					continue;
				}
				dspmm(&P(GMngs,iwfs), P(recon->saneai,iwfs,iwfs),
					P(recon->GFall,iwfs), "nn", 1);
				dmm(&GMGngs, 1, P(recon->GFall,iwfs), P(GMngs,iwfs), "tn", 1);
			}
			dinvspd_inplace(GMGngs);
			/*A focus reconstructor from all NGS measurements.*/
			dcell* RFngsg=recon->RFngsg=dcellnew(1, parms->nwfs);

			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				if(!P(recon->GFall,iwfs)) continue;
				//NGS gradient to Focus mode reconstructor.
				dmm(&P(RFngsg,iwfs), 0, GMGngs, P(GMngs,iwfs), "nt", 1);
			}
			dfree(GMGngs);
			dcellfree(GMngs);
		}
		/*
		  Compute focus constructor from LGS grads. A constructor for each LGS
		  because each LGS may have different range error. Applyes to parms->wfs, not parms->wfs.
		*/
		cellfree(recon->RFlgsg);
		recon->RFlgsg=dcellnew(parms->nwfs, parms->nwfs);
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(!parms->powfs[ipowfs].llt) continue;
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
				int iwfs0=P(parms->powfs[ipowfs].wfs,0);
				if(iwfs==iwfs0||!parms->recon.glao){
					P(recon->RFlgsg, iwfs, iwfs)=dpinv(P(recon->GFall,iwfs), CELL(P(recon->saneai, iwfs, iwfs)));
				} else{
					P(recon->RFlgsg, iwfs, iwfs)=dref(P(recon->RFlgsg, iwfs0, iwfs0));
				}
			}
		}

		if(parms->save.setup){
			writebin(recon->RFngsg, "focus_RFngsg");
			writebin(recon->RFlgsg, "focus_RFlgsg");
		}
	}
	if(parms->sim.focus2tel){
		dcell* Fdm=dcellnew(parms->ndm, 1);
		recon->RFdm=dcellnew(1, parms->ndm);
		for(int idm=0; idm<parms->ndm; idm++){
			if(idm!=parms->idmground) continue;
			P(Fdm,idm)=dnew(P(recon->anloc,idm), 1);
			loc_add_focus(P(Fdm,idm), P(recon->aloc,idm), 1);
			P(recon->RFdm,idm)=dpinv(P(Fdm,idm), 0);
		}
		dcellfree(Fdm);
		if(parms->save.setup){
			writebin(recon->RFdm, "focus_RFdm");
		}
	}
}

/**
   Setup reconstructor for TWFS or sodium fit gradient
*/
static void
setup_recon_twfs(recon_t* recon, const parms_t* parms){
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
		real thres=1e-10;
		info("RRtwfs svd threshold is %g\n", thres);
		cellfree(recon->RRtwfs);
		recon->RRtwfs=dcellpinv2(recon->GRtwfs, CELL(neai), thres, 0);
	}

	if(parms->save.setup){
		if(recon->RRtwfs) writebin(recon->RRtwfs, "twfs_recon");
	}
	cellfree(neai);
}

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
	dcelladdsp(&recon->GXL, 1, recon->GXlo, 1);
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
		dspcellmm(&Q, recon->fit->HA, FUw, "nn", 1);
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
	dcellmm_cell(&recon->MVGM, recon->GAlo, recon->MVModes, "nn", 1);
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
		dspcellmm(&Qn, recon->fit->HA, recon->MVModes, "nn", 1);
		dcell* Qntt=dcellnew(NX(Qn), NY(Qn));
		dmat* TTploc=loc2mat(recon->floc, 1);/*TT mode. need piston mode too! */
		dmat* PTTploc=dpinv(TTploc, CELL(recon->W0));/*TT projector. no need w1 since we have piston. */
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
void setup_recon_tomo(recon_t* recon, const parms_t* parms, powfs_t* powfs){
	TIC;tic;
	/*setup inverse noise covariance matrix. */
	/*prepare for tomography setup */
	setup_recon_tomo_prep(recon, parms);
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
			recon->RL.pdata=(void*)recon;
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


/**
   Setup either the control matrix using either minimum variance reconstructor by calling setup_recon_mvr()
   or least square reconstructor by calling setup_recon_lsr() 
   It can be called repeatedly to update the reconstructor when NEA updates.

   The results are measurement noise dependent and may be updated during simulation.
   */
void setup_recon_control(recon_t* recon, const parms_t* parms, powfs_t* powfs){
	info("Setup or update control matrix parameters.\n");
	TIC;tic;
	/*assemble noise equiva angle inverse from powfs information */
	setup_recon_saneai(recon, parms, powfs);
	/*setup LGS tip/tilt/diff focus removal */
	setup_recon_TTFR(recon, parms);
	/*mvst uses information here*/
	setup_recon_focus(recon, parms);
	/*setup Truth wfs*/
	setup_recon_twfs(recon, parms);
	
	if(!parms->sim.idealfit&&!parms->sim.idealtomo){
		if(parms->recon.mvm&&parms->load.mvm){
			recon->MVM=dread("%s", parms->load.mvm);
		} else{
			switch(parms->recon.alg){
			case 0:
				setup_recon_tomo(recon, parms, powfs);
				break;
			case 1:
				setup_recon_lsr(recon, parms);
				break;
			default:
				error("recon.alg=%d is not recognized\n", parms->recon.alg);
			}
		}
	}
	toc2("setup_recon_control");
}

/**
   PSD computation for gain update
 */
void setup_recon_psd(recon_t* recon, const parms_t* parms){
	if(!parms->recon.psd) return;
	recon->Herr=dspcellnew(parms->evl.nevl, parms->ndm);
	real d1=parms->aper.d-parms->dm[0].dx;
	real d2=parms->aper.din+parms->dm[0].dx;
	real dx=(d1-d2)*0.1;
	if(dx<parms->dm[0].dx){
		dx=parms->dm[0].dx;
	}
	loc_t* eloc=mkannloc(d1, d2, dx, 0.8);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		for(int idm=0; idm<parms->ndm; idm++){
			real ht=parms->dm[idm].ht;
			real dispx=P(parms->evl.thetax,ievl)*ht;
			real dispy=P(parms->evl.thetay,ievl)*ht;
			real scale=1-ht/P(parms->evl.hs,ievl);
			P(recon->Herr, ievl, idm)=mkh(P(recon->aloc,idm), eloc, dispx, dispy, scale);
		}
	}
	if(parms->recon.psd==2){//don't use signanhi by default
		dcell* ecnn=setup_recon_ecnn(recon, parms, eloc, 0);
		real sigma2e=0;
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			real sigma2i=0;
			for(int iloc=0; iloc<eloc->nloc; iloc++){
				sigma2i+=P(P(ecnn,ievl), iloc, iloc);
			}
			sigma2e+=P(parms->evl.wt,ievl)*(sigma2i/eloc->nloc);
		}
		recon->sigmanhi=sigma2e;
		info("High order WFS mean noise propagation is %g nm\n", sqrt(sigma2e)*1e9);
		if(parms->save.setup){
			writebin(ecnn, "psd_ecnn");
		}
		dcellfree(ecnn);
	}
	if(parms->save.setup){
		writebin(recon->Herr, "psd_Herr");
		writebin(eloc, "psd_eloc.bin");
	}
	locfree(eloc);
}

/**
   A few further operations that needs MVM.
 */
void setup_recon_post(recon_t* recon, const parms_t* parms, const aper_t* aper){
	TIC;tic;
	if(parms->sim.ecnn){
		recon->ecnn=setup_recon_ecnn(recon, parms, aper->locs, parms->evl.psfr);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			char strht[24];
			if(!isinf(P(parms->evl.hs,ievl))){
				snprintf(strht, 24, "_%g", P(parms->evl.hs,ievl));
			} else{
				strht[0]='\0';
			}
			writebin(P(recon->ecnn,ievl), "ecnn_x%g_y%g%s.bin",
				P(parms->evl.thetax,ievl)*206265,
				P(parms->evl.thetay,ievl)*206265, strht);
		}
	}
	if(parms->recon.psd){
		setup_recon_psd(recon, parms);
	}
	toc2("setup_recon_post");
}

/**
   Free unused object in recon struct after preparation is done.
 */
void free_recon_unused(const parms_t* parms, recon_t* recon){
	if(!recon) return;
	/* Free arrays that will no longer be used after reconstruction setup is done. */
	dspcellfree(recon->sanea);
	dspcellfree(recon->saneal);
	if(!(parms->tomo.assemble&&parms->tomo.alg==1)&&!parms->cn2.tomo&&!parms->tomo.bgs){
	/*We no longer need RL.M,U,V */
		cellfree(recon->RL.M);
		dcellfree(recon->RL.U);
		dcellfree(recon->RL.V);
		cellfree(recon->RR.M);
		dcellfree(recon->RR.U);
		dcellfree(recon->RR.V);
	} else{
		dcellfree(recon->TTF);
		dcellfree(recon->PTTF);
	}

	if(parms->tomo.alg==1){
		muv_direct_free(&recon->RL);
	}

	if(recon->RR.M){
		dspcellfree(recon->GP);
	}

	/*
	  The following arrys are not used after preparation is done.
	*/
	dspcellfree(recon->GX);
	dspcellfree(recon->GXtomo);/*we use HXWtomo instead. faster */
	if(!(parms->cn2.tomo&&parms->recon.split==2)){/*mvst needs GXlo when updating. */
		dspcellfree(recon->GXlo);
	}
	if(parms->tomo.square&&!parms->dbg.tomo_hxw&&recon->RR.M){
		dspcellfree(recon->HXWtomo);
	}
	if(parms->recon.mvm){
		muv_free(&recon->RR);
		muv_free(&recon->RL);
		muv_free(&recon->LR);
		muv_free(&recon->LL);
		free_fit(recon->fit, 1); recon->fit=NULL;
		fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
		if(parms->gpu.tomo&&parms->gpu.fit){
			dfree(recon->MVM);//keep GPU copy.
		}
	}
}
/**
   Free the recon struct.
*/
void free_recon(const parms_t* parms, recon_t* recon){
	if(!recon) return;
	ngsmod_free(recon->ngsmod); recon->ngsmod=0;
	free_recon_unused(parms, recon);
	free_recon_moao(recon, parms);
	free_fit(recon->fit, 1); recon->fit=NULL;
	dfree(recon->ht);
	dfree(recon->os);
	dfree(recon->wt);
	dfree(recon->dx);
	dcellfree(recon->MVRngs);
	dcellfree(recon->MVGM);
	dcellfree(recon->MVFM);
	dcellfree(recon->MVModes);
	dspcellfree(recon->GX);
	dspcellfree(recon->GXlo);
	dspcellfree(recon->GXtomo);
	dspcellfree(recon->GP);
	cellfree(recon->GA);
	cellfree(recon->GAlo);
	cellfree(recon->GAhi);
	cellfree(recon->GM);
	cellfree(recon->GMhi);
	dcellfree(recon->GXL);
	dspcellfree(recon->L2);
	dspcellfree(recon->L2save);
	free_cxx(recon);
	dcellfree(recon->TT);
	dcellfree(recon->PTT);
	dcellfree(recon->FF);
	dcellfree(recon->PFF);
	dcellfree(recon->TTF);
	dcellfree(recon->PTTF);
	dcellfree(recon->GFngs);
	dcellfree(recon->GFall);
	dcellfree(recon->RFlgsg);
	dcellfree(recon->RFngsg);
	dcellfree(recon->RFdm);
	dcellfree(recon->GRall);
	dcellfree(recon->RRtwfs);
	dspcellfree(recon->ZZT);
	dspcellfree(recon->HXW);
	dspcellfree(recon->HXWtomo);
	cellfree(recon->Herr);
	dspfree(recon->W0);
	dfree(recon->W1);

	cellfree(recon->xloc);
	cellfree(recon->xmap);
	cellfree(recon->xcmap);
	lfree(recon->xnx);
	lfree(recon->xny);
	lfree(recon->xnloc);
	dcellfree(recon->xmcc);

	lfree(recon->anx);
	lfree(recon->any);
	lfree(recon->anloc);
	lfree(recon->ngrad);
	locfree(recon->floc);
	locfree(recon->ploc);
	cellfree(recon->ploc_tel);
	free(P(recon->amap));free(recon->amap);//data is referenced
	cellfree(recon->amod);
	cellfree(recon->anmod);
	cellfree(recon->acmap);
	cellfree(recon->aloc);
	cellfree(recon->actstuck);
	cellfree(recon->actfloat);
	cellfree(recon->actextrap);
	cellfree(recon->actcpl);
	cellfree(recon->aimcc);/*used in filter.c */
	muv_free(&recon->RR);
	muv_free(&recon->RL);
	muv_free(&recon->LR);
	muv_free(&recon->LL);
	dfree(recon->MVM);
	dspcellfree(recon->sanea);
	dspcellfree(recon->saneal);
	dspcellfree(recon->saneai);
	dfree(recon->neam);
	fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
	cn2est_free(recon->cn2est);
	cellfree(recon->saloc);
	cellfree(recon->DMTT);
	cellfree(recon->DMPTT);
	cellfree(recon->dither_m);
	cellfree(recon->dither_rg);
	cellfree(recon->dither_ra);
	//cellfree(recon->GSF);
	//cellfree(recon->RSF);
	free(recon);
}


/*
  some tips.
  1) UV balance need to be done carefully. The scaling must be a single number
  2) When add regularizations, the number has to be in the same order as the matrix.
  This is supper important!!!
  3) single point piston constraint in Tomography is not implemented correctly.

  2009-11-22
  Found a difference: Ha in LAOS computes the weighting even if not enough points.
  HA in MAOS only computes the weighting only if 4 coupled points all exist.
  AOS used laos geometry that are not enough to cover all points. Will modify mkH/accphi.
  H modified. Performance now agree with laos with AOS W0/W1 and HA

*/
