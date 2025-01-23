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
#include "fdpcg.h"
#include "ahst.h"
#include "recon_utils.h"
#include "moao.h"
#include "powfs.h"
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
void setup_recon_saneai(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	const int nwfs=parms->nwfsr;
	dspcellfree(recon->sanea);
	dspcellfree(recon->saneai);
	dspcellfree(recon->saneal);
	dspcell* sanea=recon->sanea=dspcellnew(nwfs, nwfs);//The subaperture NEA
	dspcell* saneal=recon->saneal=dspcellnew(nwfs, nwfs);//The LL' decomposition of sanea
	dspcell* saneai=recon->saneai=dspcellnew(nwfs, nwfs);//The inverse of sanea
	dfree(recon->neam);
	recon->neam=dnew(parms->nwfsr, 1);
	if(parms->load.saneai){//load from file.
		recon->saneai=dspcellread("%s", parms->load.saneai);
		int mismatch=0;
		if(NX(recon->saneai)!=parms->nwfs){
			mismatch=1;
		}else{
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				const int ipowfs=parms->wfsr[iwfs].powfs;
				const int nsa=powfs[ipowfs].saloc->nloc;
				if(NX(recon->saneai, iwfs, iwfs)!=nsa*2){
					mismatch=1;
				}
				real avg=dsptrace(P(recon->saneai, iwfs, iwfs), -1)/(2*nsa);
				P(recon->neam, iwfs)=sqrt(avg);
			}
		}
		if(mismatch){
			error("load.saneai has wrong dimensions\n");
		}
	}else{//compute
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].skip==3) continue;
		const int nsa=powfs[ipowfs].saloc->nloc;
		const int ng=parms->powfs[ipowfs].ng;
		
		const real neathres=0.33*(parms->powfs[ipowfs].type==WFS_SH?parms->powfs[ipowfs].pixtheta:parms->powfs[ipowfs].fieldstop);
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
			for(int jwfs=0; jwfs<PN(saneac); jwfs++){
				real saat=parms->powfs[ipowfs].saat; if(!saat) saat=0.3;
				dmat *saa=PR(powfs[ipowfs].saa, jwfs);
				for(int isa=0; isa<nsa; isa++){
					real area=P(saa, isa);
					if(area<saat){
						for(int ig=0; ig<ng; ig++){
							P(P(saneac, jwfs), isa, ig)*=1e4;
						}
					}
				}
			}
		} else{
			if(parms->powfs[ipowfs].neareconfile){
				saneac=readwfs(parms->powfs[ipowfs].neareconfile, parms, ipowfs);
				for(int i=0; i<NX(saneac)*NY(saneac); i++){
					nea_check(P(saneac,i), nsa, ng);
					nea_mm(&P(saneac,i), P(saneac,i), ng);
				}
			} else{
				saneac=dcellnew(PN(powfs[ipowfs].saa), 1);
				real neamas=parms->powfs[ipowfs].nearecon;
				if(neamas<0.001||neamas > 2000){
					warning("powfs[%d].nearecon=%g mas may have unit incorrect.\n", ipowfs, neamas);
				}
				//convert from mill-arcsec to radian.
				real nearad=pow(neamas*MAS2RAD, 2)/(parms->powfs[ipowfs].dtrat*parms->sim.dt/parms->sim.dtref);
				for(int jwfs=0; jwfs<PN(powfs[ipowfs].saa); jwfs++){
					P(saneac,jwfs)=dnew(nsa, ng);//in this mode, no x/y cross-coupling.
					real saat=parms->powfs[ipowfs].saat; if(!saat) saat=0.3;
					for(int isa=0; isa<nsa; isa++){
						real area=P(P(powfs[ipowfs].saa, jwfs), isa);
						real nearad2=nearad/(area>saat?area:1e-4);
						for(int ig=0; ig<ng; ig++){
							P(P(saneac, jwfs), isa, ig)=nearad2;
						}
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
		const real neaextra2=copysign(pow(parms->powfs[ipowfs].neaextra*MAS2RAD, 2),
			parms->powfs[ipowfs].neaextra);
		const real neamin2=pow(parms->powfs[ipowfs].neamin*MAS2RAD, 2);

		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfsr; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfsr,jwfs);
			int iwfs0=P(parms->powfs[ipowfs].wfsr,0);
			dmat *saa=PR(powfs[ipowfs].saa, jwfs);
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
					if(P(saa, isa)>area_thres){
						for(int iy=0; iy<ng; iy++){
							nea2_sum+=P(sanea0, isa, iy);
							nea2_count++;
						}
					}
				}
				
				nea_chol(&sanea0l, sanea0, ng);
				nea_inv(&sanea0i, sanea0, ng, TOMOSCALE);//without TOMOSCALE, it overflows in float mode
				
				real nea_mean=sqrt(nea2_sum/nea2_count);
				P(recon->neam,iwfs)=nea_mean/sqrt(TOMOSCALE);//saneai is scaled by TOMOSALE. 
				if(nea_mean>neathres
					&&parms->powfs[ipowfs].usephy
					&&parms->powfs[ipowfs].order<=2
					&&parms->sim.dtrat_lo==parms->sim.dtrat_lo2
					){
					warning("TT WFS %d has too much measurement error: %g mas\". Ignore it\n",
						iwfs, nea_mean*RAD2MAS);
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
		}/*iwfs*/
		dcellfree(saneac);
	}/*ipowfs */
	}
	info2("Recon NEA: ");
	real neam_hi=0;
	int count_hi=0;
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
			info2("%s(%.2f) ", neatype, P(recon->neam,iwfs)*RAD2MAS*sqrt(TOMOSCALE));
			if(!parms->powfs[ipowfs].lo){
				neam_hi+=pow(P(recon->neam,iwfs), 2);
				count_hi++;
			}
		}
	}
	recon->neamhi=sqrt(neam_hi/count_hi);
	info2(" mas\n");
	if(parms->save.setup){
		//writebin(recon->sanea, "sanea");//used for mvst
		writebin(recon->saneai, "saneai");//used for reconstruction
		//writebin(recon->saneal, "saneal");//used for ecnn
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
	ecnn is the wavefront estimation error covariance in science focal
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
static dcell *setup_recon_ecnn(const recon_t *recon, const parms_t *parms, const loc_t *locs, const lmat *mask){
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
		dmat* tmp=dcell2m(sanealhi); dspcellfree(sanealhi);
		dmm(&t1, 0, recon->MVM, tmp, "nn", 1);
		cellfree(tmp);
		toc2("MVM ");tic;
	} else if(parms->recon.alg==RECON_MVR){//MV
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
	setup_recon_tomo_reg(recon, parms); /*redo L2, invpsd */
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
			recon->RFngsg=dcellnew(1, parms->nwfs);

			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				if(!P(recon->GFall,iwfs)) continue;
				//NGS gradient to Focus mode reconstructor.
				dmm(&P(recon->RFngsg,iwfs), 0, GMGngs, P(GMngs,iwfs), "nt", 1);
			}
			dfree(GMGngs);
			dcellfree(GMngs);
			if(parms->save.setup){
				writebin(recon->RFngsg, "focus_RFngsg");
			}
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
					P(recon->RFlgsg, iwfs, iwfs)=dpinv(P(recon->GFall,iwfs), P(recon->saneai, iwfs, iwfs));
				} else{
					P(recon->RFlgsg, iwfs, iwfs)=dref(P(recon->RFlgsg, iwfs0, iwfs0));
				}
			}
		}
		if(parms->save.setup){
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
		cellfree(recon->RRtwfs);
		recon->RRtwfs=dcellpinv(recon->GRtwfs, neai);
	}

	if(parms->save.setup){
		if(recon->RRtwfs) writebin(recon->RRtwfs, "twfs_recon");
	}
	cellfree(neai);
}


/**
   Setup either the control matrix using either minimum variance reconstructor by calling setup_recon_mvr()
   or least square reconstructor by calling setup_recon_lsr() 
   It can be called repeatedly to update the reconstructor when NEA updates.

   The results are measurement noise dependent and may be updated during simulation.
   */
void setup_recon_control(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	info("Compute control matrix parameters.\n");
	TIC;tic;
	/*setup LGS tip/tilt/diff focus removal */
	setup_recon_TTFR(recon, parms);
	/*mvst uses information here*/
	setup_recon_focus(recon, parms);
	/*setup Truth wfs*/
	setup_recon_twfs(recon, parms);
	
	if(!parms->sim.idealtomo){
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

	if(parms->recon.alg==RECON_MVR||parms->sim.dmproj){
		setup_recon_fit(recon, parms);
	}
	if(parms->recon.split){/*split tomography */
		ngsmod_setup(parms, recon);
		if(parms->recon.split==2&&parms->recon.alg==RECON_MVR){/*Need to be after fit */
			setup_recon_mvst(recon, parms);
		}
	}

	toc2("setup_recon_control");
}

/**
   PSD computation for gain update
 */
void setup_recon_psd(recon_t* recon, const parms_t* parms){
	if(!parms->recon.psd) return;
	recon->Herr=cellnew(parms->evl.nevl, parms->ndm);
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
			dsp* Htmp=mkh(P(recon->aloc,idm), eloc, dispx, dispy, scale, 0);
			
			if(parms->recon.modal){
				dmat* Hdtmp=NULL;
				dspmm(&Hdtmp, Htmp, P(recon->amod, idm, idm), "nn", 1);
				dspfree(Htmp);
				P(recon->Herr, ievl, idm)=(cell*)Hdtmp;
			}else{
				P(recon->Herr, ievl, idm)=(cell*)Htmp;
			}
		}
	}
	if(parms->recon.psd==2){//don't use sigmanhi by default
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
 * Petaling mode control
*/
void setup_recon_petal(recon_t *recon, const parms_t *parms, const powfs_t *powfs){
	if(!parms->recon.petal) return;
	recon->petal=mycalloc(2, petal_t*);
	for(int ir=0; ir<2; ir++){
		int ipowfs=0;
		if(ir==0&&parms->ittfpowfs!=-1){
			ipowfs=parms->ittfpowfs;
		} else if(ir==1&&parms->ittpowfs!=-1){
			ipowfs=parms->ittpowfs;
		} else{
			continue;
		}
		real nembed=2;
		real dsa=powfs[ipowfs].pts->dsa;
		real dtheta=parms->powfs[ipowfs].wvlmean/(nembed*dsa);
		real pdtheta=parms->powfs[ipowfs].pixtheta/dtheta;
		dbg("powfs[%d].pdtheta=%g\n", ipowfs, pdtheta);
		if(fabs(pdtheta-1)>0.01){
			warning("TODO: pdtheta!=1 requries resampling PSFs\n");
		}
		//only withtt only for t/t oiwfs unless petaltt>1. 
		//enable it for TTF OIWFS sometimes results in a clocking gradient pattern.
		int withtt=(parms->powfs[ipowfs].order==1||parms->recon.petaltt>1)?parms->recon.petaltt:0;
		recon->petal[ir]=petal_setup(powfs[ipowfs].pts->loc, powfs[ipowfs].loc->dx, P(powfs[ipowfs].amp, 0),
			pdtheta, parms->powfs[ipowfs].pixblur, parms->aper.rot, parms->recon.petalnpsf, withtt);
		if(parms->save.setup){
			petal_save(recon->petal[ir], "petal_%d", ir);
		}
	}
}
/**
   Reconstruction optimization background loops setup.
 */
void setup_recon_misc(recon_t* recon, const parms_t* parms, loc_t* locs, const powfs_t *powfs){
	TIC;tic;
	if(parms->sim.ecnn){
		recon->ecnn=setup_recon_ecnn(recon, parms, locs, parms->evl.psfr);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			char strht[24];
			if(!isinf(P(parms->evl.hs,ievl))){
				snprintf(strht, 24, "_%g", P(parms->evl.hs,ievl));
			} else{
				strht[0]='\0';
			}
			writebin(P(recon->ecnn,ievl), "ecnn_x%g_y%g%s.bin",
				P(parms->evl.thetax,ievl)*RAD2AS,
				P(parms->evl.thetay,ievl)*RAD2AS, strht);
		}
	}
	if(parms->recon.psd){
		setup_recon_psd(recon, parms);
	}
	if(parms->recon.petal){
		int idm=parms->idmground;
		recon->apetal=petal_mkh_loc(P(recon->aloc,idm), 6, parms->aper.rot);
		setup_recon_petal(recon, parms, powfs);
		if(parms->save.setup){
			writebin(recon->apetal, "apetal");
		}
	}
	toc2("setup_recon_misc");
}
/**
   Frees recon->invpsd or recon->fractal
*/
void 
free_recon_cxx(recon_t* recon){
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
	cellfree(recon->amodpinv);
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
	dcellfree(recon->GXL);
	dspcellfree(recon->L2);
	dspcellfree(recon->L2save);
	free_recon_cxx(recon);
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
	dcellfree(recon->GRtwfs);
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
	free(P(recon->amap));free(recon->amap);//data is referenced
	cellfree(recon->amod);
	cellfree(recon->amodpinv);
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
	//cellfree(recon->dither_rm);
	cellfree(recon->dither_ra);
	//cellfree(recon->GSF);
	//cellfree(recon->RSF);
	cellfree(recon->apetal);
	petal_free_arr(recon->petal, 2);
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
