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
#include "dither_utils.h"
#include "powfs_utils.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

/**
 * @brief Calculate expected averaged position of dithering signal during WFS integration. Called when (isim+1)%dtrat=0
 * This function calculates the dither position based on the given parameters.
 * @param cs The cosine of the dither angle.
 * @param ss The sine of the dither angle.
 * @param alfsm The additional latency due to uplink.
 * @param dtrat The dither rate.
 * @param npoint The number of points.
 * @param isim The simulation index.
 * @param deltam The dither delta.
*/
void dither_position(real* cs, real* ss, int alfsm, int dtrat, int npoint, int isim, real deltam){
	//adjust for delay due to propagation, and computation delay if delay is not 2 frame.
	const int adjust=alfsm+1-dtrat;
	//adjust to get delay at beginning of integration
	const int adjust2=dtrat-1;
	const real anglei=(2*M_PI/npoint);
	const real angle0=M_PI*0.5;
	const real angle=angle0+((isim-adjust-adjust2)/dtrat)*anglei+deltam;
	const real angle2=angle0+((isim-adjust)/dtrat)*anglei+deltam;
	const real delay=(real)adjust/dtrat;
	const real beta=1+delay+floor(-delay);
	const real scale=1./(beta*beta+(1-beta)*(1-beta));
	//use average of two places during accumulation and scale
	*cs=(beta*cos(angle)+(1-beta)*cos(angle2))*scale;
	*ss=(beta*sin(angle)+(1-beta)*sin(angle2))*scale;
}
void gradoff_acc(sim_t* simu, int ipowfs){
	(void)ipowfs;
	if(simu->parms->dbg.gradoff_reset==2&&simu->gradoffisim0>0){//accumulate gradoff before updating it.
		int nsim=(simu->wfsisim-simu->gradoffisim);
		if(nsim){
			dcelladd(&simu->gradoffacc, 1, simu->gradoff, nsim);
			info("step %d: gradoff is accumulated with factor %d\n", simu->wfsisim, nsim);
		}
		simu->gradoffisim=simu->wfsisim;
	}
}
/**
   Demodulate the dithering signal to determine the amplitude. Remove trend (detrending) if detrend is set.
 * @param res The result matrix. nmod*1 if not combine.
 * @param signal The signal matrix. nmod*nsim *
 * @param dtrat skip columns due to wfs/sim dt ratio
 * @param npoint The number of points during dithering
 * @param detrend Whether to detrend the signal.
 * @param combine Whether to combine the results. for tip/tilt only
 * @return The calculated dither amplitude.
*/
real dither_amp_calc(dmat **res, const dmat *signal, long dtrat, long npoint, int detrend, int combine){
	long nmod=NY(signal)==1?1:NX(signal);
	long nframe=((NY(signal)==1?signal->nx:signal->ny)-1)/dtrat+1;//number of actual frames
	dmat *slope=dnew(nmod, 1);//for detrending
	long offset=(nframe/npoint-1)*npoint;//number of WFS frame separations between first and last cycle
	if(detrend&&offset){//detrending
		for(long imod=0; imod<nmod; imod++){
			for(long ip=0; ip<npoint; ip++){
				long i0=ip*dtrat*nmod+imod;
				long i1=(ip+offset)*dtrat*nmod+imod;
				P(slope,imod)+=P(signal, i1)-P(signal, i0);
			}
			P(slope,imod)/=(npoint*offset);
		}
		//dbg("slope=%g. npoint=%ld, nmod=%ld, nframe=%ld, offset=%ld\n", slope, npoint, nmod, nframe, offset);
	}
	real anglei=M_PI*2/npoint;
	real angle0=M_PI*0.5;//2023-05-29: bias added to work with npoint==2
	real ipv=0, qdv=0, a2m=0;
	if(combine){//tip and tilt dithering
		if(nmod!=2){
			error("combine only support nmod=2\n");
		}
		for(int iframe=0; iframe<nframe; iframe++){
			real angle=angle0+anglei*iframe;//position of dithering
			real cs=cos(angle);
			real ss=sin(angle);
			real ttx=P(signal, iframe*dtrat*nmod)-P(slope,0)*iframe;
			real tty=P(signal, iframe*dtrat*nmod+1)-P(slope,1)*iframe;
			ipv+=(ttx*cs+tty*ss);
			qdv+=(ttx*ss-tty*cs);
		}
		a2m=sqrt(ipv*ipv+qdv*qdv)/nframe;
		
	}else{//independent mode
		if(nmod>1&&!res){
			error("when there are more than 1 mode. res must be used to return the results\n");
		}
		if(res){
			dinit(res, nmod, 1);
		}
		dmat *ipq=dnew(2,nmod);
		for(int iframe=0; iframe<nframe; iframe++){
			real angle=angle0+anglei*iframe;//position of dithering
			real cs=cos(angle);
			real ss=sin(angle);
			for(int imod=0; imod<nmod; imod++){
				real mod=P(signal, iframe*dtrat*nmod+imod)-P(slope,imod)*iframe;
				P(ipq,0,imod)+=(mod*cs);//ipv
				P(ipq,1,imod)+=(mod*ss);//iqv
			}
		}
		for(int imod=0; imod<nmod; imod++){
			a2m=sqrt(P(ipq,0,imod)*P(ipq,0,imod)+P(ipq,1,imod)*P(ipq,1,imod))/nframe*2.;
			if(res) P(*res, imod)=a2m;
		}
		if(res){
			a2m=P(*res,0);//first mode
		}
		dfree(ipq);
	}
	dfree(slope);
	return a2m;
}
static void wfsgrad_tt_drift(dmat* grad, sim_t* simu, real gain, int iwfs, int remove){
	//gain can be set to 1 if the rate is slower than the main tip/tilt and focus control rate.
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int iwfsr=parms->recon.glao?ipowfs:iwfs;

	if(parms->powfs[ipowfs].trs){
		/* Update T/T drift signal to prevent matched filter from drifting*/
		dmat* tt=dnew(2, 1);
		const dmat* PTT=P(recon->PTT, iwfsr, iwfsr);
		dmm(&tt, 0, PTT, grad, "nn", 1);
		if(remove){
			dmat* TT=P(recon->TT, iwfsr, iwfsr);
			dmm(&grad, 1, TT, tt, "nn", -1);
		}
		//need to use integrator here and then use as offset
		if(gain){
			P(P(simu->fsmerr_drift, iwfs), 0)+=P(tt, 0)*gain;
			P(P(simu->fsmerr_drift, iwfs), 1)+=P(tt, 1)*gain;
			if(iwfs==P(parms->powfs[ipowfs].wfs, 0)){
				dbg("Step %5d: wfs %d uplink drift control error is (%.3f, %.3f) mas output is (%.3f, %.3f) mas.\n",
				simu->wfsisim, iwfs, P(tt, 0)*RAD2MAS, P(tt, 1)*RAD2MAS,
				P(P(simu->fsmerr_drift, iwfs), 0)*RAD2MAS, P(P(simu->fsmerr_drift, iwfs), 1)*RAD2MAS);
			}
		}
		dfree(tt);
	}
}
static void wfsgrad_focus_drift(dmat* grad, sim_t* simu, real gain, int iwfs, int remove){
	//gain can be set to 1 as zoomgain is applied at zoomint.
	const parms_t* parms=simu->parms;
	const int ipowfs=parms->wfs[iwfs].powfs;
	//Output focus error in ib to trombone error signal.
	if(parms->powfs[ipowfs].llt){
		//here we don't use RFlgsg which is noise weighted
		real focus=loc_remove_focus_grad(grad, simu->powfs[ipowfs].saloc, remove?1:0);
		if(gain){
			P(simu->zoomdrift, iwfs)+=gain*focus;//accumulate to the averager
			P(simu->zoomdrift_count, iwfs)++;
		}
	}
}

/**
	Controls gradient drift of individual subapertures due to dithering for matched filter
*/
static void wfsgrad_sa_drift(sim_t* simu, int ipowfs){
	const parms_t* parms=simu->parms;
	const powfs_t* powfs=simu->powfs;
	if(parms->powfs[ipowfs].dither_gdrift==0||parms->powfs[ipowfs].phytype_sim!=PTYPE_MF) return;
	dmat* goff=0;
	intstat_t* intstat=simu->powfs[ipowfs].intstat;
	const int isim=simu->wfsisim;
	gradoff_acc(simu, ipowfs);

	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
		const int nsa=NX(intstat->i0);
		//outer loop to prevent i0 from drifting
		//Compute CoG of i0 + goff and drive it toward gradoff with low gain (0.1)
		if(1){
			dzero(goff);
			shwfs_grad(&goff, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);//force cog
			dadd(&goff, 1, P(simu->gradoff, iwfs), 1);
			if(simu->gradoffdrift){//cog boot strapped
				dadd(&goff, 1, P(simu->gradoffdrift, iwfs), -1);
			} else if(simu->powfs[ipowfs].gradoff){//cmf boot strapped, gradoff is cog of i0
				dadd(&goff, 1, P(simu->powfs[ipowfs].gradoff, jwfs), -1);
			}
			dadd(&P(simu->gradoff, iwfs), 1, goff, -parms->powfs[ipowfs].dither_gdrift);
			if(jwfs==0){
				dbg("Step %5d: powfs %d gradient drift control\n", simu->wfsisim, ipowfs);
			}
		}
		if(1){
			//outer loop to prevent gx/gy direction from drifting.
			//It computes CoG of shifted images (i0+gx/gy) and make sure the angle stays the same.
			//may not be necessary.
			dmat* i0sx=0, * i0sy=0;
			real theta=0;
			const real gyoff=M_PI*0.5;
			const real gshift=parms->powfs[ipowfs].pixtheta*0.1;
			const real cogthres=parms->powfs[ipowfs].cogthres;
			const real cogoff=parms->powfs[ipowfs].cogoff;

			for(int isa=0; isa<nsa; isa++){
				real g0[2], gx[2], gy[2];
				dcp(&i0sx, PR(intstat->i0, isa, jwfs));
				dcog(g0, i0sx, 0., 0., cogthres, cogoff, 0, NULL);
				dcp(&i0sy, PR(intstat->i0, isa, jwfs));
				dadd(&i0sx, 1, PR(intstat->gx, isa, jwfs), gshift);
				dadd(&i0sy, 1, PR(intstat->gy, isa, jwfs), gshift);
				dcog(gx, i0sx, 0., 0., cogthres, cogoff, 0, NULL);
				dcog(gy, i0sy, 0., 0., cogthres, cogoff, 0, NULL);

				//Works in both x/y and r/a coordinate.
				theta+=(atan2(gx[1]-g0[1], gx[0]-g0[0])+atan2(gy[1]-g0[1], gy[0]-g0[0])-gyoff);
			}
			theta*=0.5/nsa;
			simu->dither[iwfs]->deltao=-theta;
			dfree(i0sx);
			dfree(i0sy);
			if(jwfs==0){
				dbg("Step %5d: wfs %d angle drift control deltao is %g.\n", isim, iwfs, simu->dither[iwfs]->deltao);
			}
		}
	}
	dfree(goff);
	if(parms->save.gradoff||parms->save.dither){
		writebin(simu->gradoff, "extra/gradoff_%d_drift", isim);
	}
}
void dither_acc_i0(const parms_t *parms, dither_t* pd, dcell *ints, int iwfs, int isim){
	/*Collect statistics with dithering*/
	const int ipowfs=parms->wfs[iwfs].powfs;
	dcelladd(&pd->imb, 1, ints, 1.);
	if(parms->powfs[ipowfs].dither==1){
		real cs, ss;
		dither_position(&cs, &ss, parms->powfs[ipowfs].alfsm, parms->powfs[ipowfs].dtrat,
			parms->powfs[ipowfs].dither_npoint, isim, pd->deltam);
		//accumulate for matched filter

		dcelladd(&pd->imx, 1, ints, cs);
		dcelladd(&pd->imy, 1, ints, ss);
	}
}
/**
   Accumulate dithering parameters
   - Every step: accumulate signal for phase detection.
   - At PLL output: determine input/output amplitude of dithering signal.
   - At Gain output:determine matched filter i0, gx, gy, or CoG gain.
   - Subtract t/t from gradients for non-comon-path (TT) dithering.

*/
void dither_acc(sim_t* simu, int iwfs){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const powfs_t* powfs=simu->powfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int iwfsr=parms->recon.glao?ipowfs:iwfs;
	const int isim=simu->wfsisim;
	const int pllrat=parms->powfs[ipowfs].dither_pllrat;
	if(!parms->powfs[ipowfs].dither||isim<parms->powfs[ipowfs].dither_pllskip){
		return;
	}
	real cs, ss; //Current phase of tip/tilt dithering signal
	dither_t* pd=simu->dither[iwfs];
	if(parms->powfs[ipowfs].dither==1){ //T/T dithering.
		//Current dithering signal phase
		dither_position(&cs, &ss, parms->powfs[ipowfs].alfsm, parms->powfs[ipowfs].dtrat,
			parms->powfs[ipowfs].dither_npoint, isim, pd->deltam);

		/* Use delay locked loop to determine the phase of actual
		dithering signal (position of LGS spot averaged over a WFS
		integration period) using measured signal (WFS global
		tip/tilt). In actual system, the LGS uplink propagation
		needs to be accounted.
		*/
		real err;
		err=(-ss*(P(P(simu->fsmerr, iwfs), 0))
			+cs*(P(P(simu->fsmerr, iwfs), 1)))/(parms->powfs[ipowfs].dither_amp);
		pd->delta+=parms->powfs[ipowfs].dither_gpll*(err/pllrat);

		//For SHWFS CoG gaim update.
		if(parms->powfs[ipowfs].type==WFS_SH&&parms->powfs[ipowfs].phytype_sim2!=PTYPE_MF&&isim>=parms->powfs[ipowfs].dither_ogskip){
			const int nsa=powfs[ipowfs].saloc->nloc;
			if(!pd->ggm){
				pd->ggm=dnew(nsa*2, 1);
			}
			for(int isa=0; isa<nsa; isa++){
				P(pd->ggm, isa)+=cs*P(P(simu->gradcl, iwfs), isa);
				P(pd->ggm, isa+nsa)+=ss*P(P(simu->gradcl, iwfs), isa+nsa);
			}
		}
	} else if(parms->powfs[ipowfs].dither>1){ //DM dithering.

		const int idm=parms->idmground;
		//Input dither signal
		dmat *mr_in=drefcols(P(pd->mr, 0), isim, 1);
		dmm(&mr_in, 0, P(recon->dither_ra, iwfs, idm), P(simu->dmreal, idm), "nn", 1);
		dfree(mr_in);
		/*P(P(pd->mr, 0), isim)=P(tmp, 0);
		if(PN(tmp)>1){//2 mode dithering
			P(P(pd->mr,2),isim)=P(tmp,1);
		}*/

		dmat *mr_out=drefcols(P(pd->mr, 1), isim, 1);
		//Measured dither signal from gradients
		dmm(&mr_out, 0, P(recon->dither_rg, iwfs, iwfs), P(simu->gradcl, iwfs), "nn", 1);

		/*P(P(pd->mr, 1), isim)=P(mr_out, 0);
		if(PN(mr_out)>1){
			P(P(pd->mr, 3), isim)=P(mr_out, 1);
		}*/
		dfree(mr_out);
	}
	if(simu->wfsflags[ipowfs].pllout&&parms->powfs[ipowfs].dither>0){
		//Synchronous detection of dither signal amplitude in input (DM) and output (gradients).
		//The ratio between the two is used for (optical) gain adjustment.
		const int npoint=parms->powfs[ipowfs].dither_npoint;
		const int ncol=(pllrat-1)*parms->powfs[ipowfs].dtrat+1;
		if(parms->powfs[ipowfs].dither==1){//TT
			//dbg("deltam=%g is updated to %g+%g=%g\n", pd->deltam, pd->delta, pd->deltao, pd->delta+pd->deltao);
			pd->deltam=pd->delta+(pd->deltao*parms->powfs[ipowfs].dither_gdrift);//output PLL
			dmat* tmp=0;
			const int detrend=parms->powfs[ipowfs].llt?0:1;
			tmp=drefcols(P(simu->save->fsmcmds, iwfs), simu->wfsisim-ncol+1, ncol);
			pd->a2m=dither_amp_calc(NULL, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 1);
			dfree(tmp);
			tmp=drefcols(P(simu->save->fsmerrs, iwfs), simu->wfsisim-ncol+1, ncol);
			pd->a2me=dither_amp_calc(NULL, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 1);
			dfree(tmp);
		} else if(parms->powfs[ipowfs].dither>1){//DM
			dmat* tmp=0;
			int detrend=1;//1: default.
			tmp=drefcols(P(pd->mr, 0), simu->wfsisim-ncol+1, ncol);//DM
			pd->a2m=dither_amp_calc(&pd->a2mv, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 0);
			dfree(tmp);
			tmp=drefcols(P(pd->mr, 1), simu->wfsisim-ncol+1, ncol);//Grad
			pd->a2me=dither_amp_calc(&pd->a2mev, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 0);
			dfree(tmp);
			/*if(PN(pd->a2mv)>1){//multi-mode
				warning_once("Use the last dithering mode\n");
				pd->a2m=P(pd->a2mv, PN(pd->a2mv)-1);
				pd->a2me=P(pd->a2mev, PN(pd->a2mev)-1);
			}*/
			/*if(PN(pd->mr)>3){
				tmp=drefcols(P(pd->mr, 2), simu->wfsisim-ncol+1, ncol);//DM
				pd->a2m2=dither_amp_calc(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
				dfree(tmp);
				tmp=drefcols(P(pd->mr, 3), simu->wfsisim-ncol+1, ncol);//Grad
				pd->a2me2=dither_amp_calc(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
				dfree(tmp);
			}*/
		}

		/*
		//Print PLL phase. Moved to OG gain print out
		if(iwfs==P(parms->powfs[ipowfs].wfs, 0)){
			const real anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
			const real scale=parms->powfs[ipowfs].dither==1?(1./parms->powfs[ipowfs].dither_amp):1;
			if(pd->a2m2){
				info2("Step %5d wfs %d PLL: delay=%.2f frame, dither amplitude=(%.2fx %.2fx), estimate=(%.2fx, %.2fx)\n",
					isim, iwfs, pd->deltam/anglei, pd->a2m*scale, pd->a2m2*scale, pd->a2me*scale, pd->a2me2*scale);
			}else{
				info2("Step %5d wfs %d PLL: delay=%.2f frame, dither amplitude=%.2fx, estimate=%.2fx\n",
					isim, iwfs, pd->deltam/anglei, pd->a2m*scale, pd->a2me*scale);
			}
		}*/
		if(simu->resdither){
			int ic=simu->wfsflags[ipowfs].pllout-1;
			P(P(simu->resdither, iwfs), 0, ic)=pd->deltam;
			P(P(simu->resdither, iwfs), 1, ic)=pd->a2m;
			P(P(simu->resdither, iwfs), 2, ic)=pd->a2me;
			if(PN(pd->a2mv)>1){
				int nm=PN(pd->a2mv);
				memcpy(&P(P(simu->resdither, iwfs), 4, ic), P(pd->a2mv), nm*sizeof(real));
				memcpy(&P(P(simu->resdither, iwfs), 4+nm, ic), P(pd->a2mev), nm*sizeof(real));
				/*P(P(simu->resdither, iwfs), 4, ic)=pd->a2m2;
				P(P(simu->resdither, iwfs), 5, ic)=pd->a2me2;*/
			}
		}
	}

	if(simu->wfsflags[ipowfs].ogacc){//Gain update statistics
		/*if(parms->dbg.gradoff_reset==2 && simu->gradoffisim0<=0
			&& parms->powfs[ipowfs].phytype_sim2==PTYPE_MF){//trigger accumulation of gradoff
			simu->gradoffisim0=simu->wfsisim;
			simu->gradoffisim=simu->wfsisim;
		}*/
		if(abs(parms->powfs[ipowfs].dither)==1){//TT Dither or i0 collection
			real scale1=1./pllrat;
			real amp=pd->a2m;
			real scale2=amp?(scale1*2./(amp)):0;
			if(pd->imb){//computer i0, gx, gy for matched filter
				dcellscale(pd->imb, scale1);
				//Accumulate data for matched filter
				dcelladd(&pd->i0, 1, pd->imb, 1);//imb was already scaled
				dcellzero(pd->imb);
				if(parms->powfs[ipowfs].dither==1){
					dcelladd(&pd->gx, 1, pd->imx, scale2);
					dcelladd(&pd->gy, 1, pd->imy, scale2);

					dcellzero(pd->imx);
					dcellzero(pd->imy);
				}
			} else if(pd->ggm&&parms->powfs[ipowfs].dither==1){//cog
				dadd(&pd->gg0, 1, pd->ggm, scale2);
				dzero(pd->ggm);
			}
		}
	}

	if(parms->powfs[ipowfs].dither==1){
		/* subtract estimated tip/tilt dithering signal to avoid perturbing the loop or dithering pattern.*/
		real amp=pd->a2me;
		real tt[2]={-cs*amp, -ss*amp};
		if(parms->powfs[ipowfs].trs){
			//info("fsmerr: %g %g %g %g\n", P(P(simu->fsmerr,iwfs),0), P(P(simu->fsmerr,iwfs),1), -tt[0], -tt[1]);
			if(!amp){//no estimate yet, do not close up FSM loop.
				P(P(simu->fsmerr, iwfs), 0)=0;
				P(P(simu->fsmerr, iwfs), 1)=0;
			} else{
				P(P(simu->fsmerr, iwfs), 0)+=tt[0];
				P(P(simu->fsmerr, iwfs), 1)+=tt[1];
			}
		}
		//also remove from gradient measurements.
		dmulvec(P(P(simu->gradcl, iwfs)), P(recon->TT, iwfsr, iwfsr), tt, 1);
	}
}

/**
   Dither update: zoom corrector, matched filter, gain ajustment, TWFS.
*/
void dither_post(sim_t* simu){
	powfs_t* powfs=simu->powfs;//not const to set instat. //todo: move intstat to simu.
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int isim=simu->wfsisim;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(!parms->powfs[ipowfs].dither || isim<parms->powfs[ipowfs].step || (isim+1)%parms->powfs[ipowfs].dtrat!=0) continue;
		const int nwfs=parms->powfs[ipowfs].nwfs;

		if(simu->wfsflags[ipowfs].ogout){//This is matched filter or cog update
			const int nsa=powfs[ipowfs].saloc->nloc;
			const real scale1=(real)parms->powfs[ipowfs].dither_pllrat/(real)parms->powfs[ipowfs].dither_ograt;
			int ptype2=parms->powfs[ipowfs].phytype_sim2;
			if(parms->powfs[ipowfs].dither==-1||(ptype2==PTYPE_MF)){
				if(parms->powfs[ipowfs].dither==-1){
					if(parms->powfs[ipowfs].llt){
						info2("Step %5d: Update sodium fit for powfs %d\n", isim, ipowfs);
					}else{
						info2("Step %5d: Update i0 processing for powfs %d\n", isim, ipowfs);
					}
				} else{
					info2("Step %5d: Update matched filter for powfs %d\n", isim, ipowfs);
				}
				//For matched filter
				if(!powfs[ipowfs].intstat){
					powfs[ipowfs].intstat=mycalloc(1, intstat_t);
				}
				intstat_t* intstat=powfs[ipowfs].intstat;
				parms->powfs[ipowfs].radgx=0;//ensure derivate is interpreted as along x/y.
				int pixpsax=powfs[ipowfs].pixpsax;
				int pixpsay=powfs[ipowfs].pixpsay;
				if(!intstat->i0||NY(intstat->i0)!=nwfs){
					dcellfree(intstat->i0);
					intstat->i0=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
				}
				if(!intstat->gx||NY(intstat->gx)!=nwfs){
					dcellfree(intstat->gx);
					dcellfree(intstat->gy);
					intstat->gx=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
					intstat->gy=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
				}
				real g2=parms->powfs[ipowfs].dither_glpf;
				if(simu->wfsflags[ipowfs].ogout*g2<1){//not enough accumulations yet.
					g2=1./(simu->wfsflags[ipowfs].ogout);
				}
				gradoff_acc(simu, ipowfs);
				if(g2<1){
					info("Applying LPF with gain %.2f to i0/gx/gy update at update cycle %d\n", g2, simu->wfsflags[ipowfs].ogout);
				}
				for(int jwfs=0; jwfs<nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dither_t* pd=simu->dither[iwfs];
					//Scale the output due to accumulation
					//TODO: remove the LPF which is not useful.
					//TODO: combine pd->i0 with intstat->i0
					//TODO: accumulate directly to intstat->i0 instead of to pd->imx
					for(int isa=0; isa<nsa; isa++){
						dadd(&P(intstat->i0, isa, jwfs), 1-g2, P(pd->i0, isa), scale1*g2);
						if(parms->powfs[ipowfs].dither==1){
							dadd(&P(intstat->gx, isa, jwfs), 1-g2, P(pd->gx, isa), scale1*g2);
							dadd(&P(intstat->gy, isa, jwfs), 1-g2, P(pd->gy, isa), scale1*g2);
						}
					}
					dcellzero(pd->i0);
					dcellzero(pd->gx);
					dcellzero(pd->gy);
					if(parms->powfs[ipowfs].dither==1){
						if(parms->dbg.gradoff_reset==0){
							if(jwfs==0) info("Step %5d: powfs%d reducing gradoff by grad of i0.\n", isim, ipowfs);
							dmat* goff=0;
							/*Compute the gradient of i0 using old gradient
							algorithm and subtract from the gradient offset to
							prevent sudden jump of gradient measurement.*/
							shwfs_grad(&goff, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
							dadd(&P(simu->gradoff, iwfs), 1, goff, -1);
							dfree(goff);
							if(parms->powfs[ipowfs].dither_glpf!=1){
								warning("when dbg.gradoff_reset is enabled, dither_glpf should be 1.\n");
							}
						} else if(parms->dbg.gradoff_reset==1){
							if(jwfs==0) info("Step %5d: powfs%d resetting gradoff to 0.\n", isim, ipowfs);
							dzero(P(simu->gradoff, iwfs));
						} else if(parms->dbg.gradoff_reset==2){
							if(jwfs==0) info("Step %5d: powfs%d reducing gradoff by its average.\n", isim, ipowfs);
							int nacc=simu->gradoffisim-simu->gradoffisim0;
							if(jwfs==0) info("Step %5d: powfs%d gradoffacc is scaled by 1/%d\n", isim, ipowfs, nacc);
							dscale(P(simu->gradoffacc, iwfs), 1./nacc);
							dadd(&P(simu->gradoff, iwfs), 1, P(simu->gradoffacc, iwfs), -1);
						}
					}
				}
				if(ptype2==PTYPE_MF&&parms->powfs[ipowfs].llt){
					//for LGS only. tip/tilt and focus drift control is needed for matched filter with either dithering or sodium fitting
					//2022-07-12: moved before the next block because sodium_fit_wrap() modifies intstat->i0 in place.
					//obsolete: use sodium fit for LGS instead. It works much better and avoids cog bootstrap problem and can work without TWFS.
					dmat *i0grad=0;
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
						shwfs_grad(&i0grad, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);
						if(parms->save.dither){
							writebin(i0grad, "extra/wfs%d_i0grad_%d", iwfs, isim);
						}
						wfsgrad_tt_drift(i0grad, simu, P(parms->sim.eplo, 0), iwfs, 0);
						wfsgrad_focus_drift(i0grad, simu, 1, iwfs, 0);
					}
					dfree(i0grad);
				}
				//i0 is collected
				if(parms->powfs[ipowfs].dither==-1){
					if(parms->powfs[ipowfs].llt){//LGS, require fiting sodium profile
						dmat* sodium=0;
						dcell* grad=0;
						//don't need gradient output for matched filter
						dcell** pgrad=(ptype2==PTYPE_COG)?&grad:NULL;
						dcell** pi0=(ptype2==PTYPE_MF)?&intstat->i0:NULL;
						dcell** pgx=(ptype2==PTYPE_MF)?&intstat->gx:NULL;
						dcell** pgy=(ptype2==PTYPE_MF)?&intstat->gy:NULL;
						//1 iteration of cog/mtch is necessary to get sodium profile, 3 iterations of mtche is necessary for gradient of i0
						int niter=parms->powfs[ipowfs].llt->na_fit_maxit?parms->powfs[ipowfs].llt->na_fit_maxit:(ptype2==PTYPE_MF?1:3);

						if(parms->save.dither){
							writebin(intstat->i0, "extra/powfs%d_i0i_%d", ipowfs, isim);
						}
						sodium_fit_wrap(&sodium, pgrad, pi0, pgx, pgy, intstat->i0, parms, powfs, ipowfs,
							recon->r0, recon->L0, niter, 1);
						if(parms->save.extra){
							writebin(sodium, "extra/powfs%d_fit_sodium_%d", ipowfs, isim);
							if(grad) writebin(grad, "extra/powfs%d_fit_grad_%d", ipowfs, isim);
						}
						if(parms->save.dither){
							if(pi0) writebin(intstat->i0, "extra/powfs%d_i0o_%d", ipowfs, isim);
							if(pgx) writebin(intstat->gx, "extra/powfs%d_gxo_%d", ipowfs, isim);
							if(pgy) writebin(intstat->gy, "extra/powfs%d_gyo_%d", ipowfs, isim);
						}

						if(ptype2==PTYPE_COG){//project pgrad to TWFS corrected modes
							dcelladd(pgrad, 1, powfs[ipowfs].gradoff, -1);
							dbg("project pgrad to TWFS corrected modes\n");
							dcell *zm=0;
							dcellmm(&zm, recon->RRlgs, *pgrad, "nn", 1);
							dcellzero(*pgrad);
							dcellmm(pgrad, recon->GRlgs, zm, "nn", 1);
							if(parms->save.dither){
								writebin(zm, "extra/powfs%d_zm_%d", ipowfs, isim);
								writebin(*pgrad, "extra/powfs%d_fit_grad_proj_%d", ipowfs, isim);
							}
							dcellfree(zm);
						}

						for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
							int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
							if(ptype2==PTYPE_COG){
								if(jwfs==0) dbg("in cog mode, gradoff+=(g_ncpa-grad)\n");
								dadd(&P(simu->gradoff, iwfs), 1, P(grad, jwfs), -1);
								//prevent gradoff from accumulating tip/tilt or focus mode if any. no need to do drift control.
								wfsgrad_tt_drift(P(simu->gradoff, iwfs), simu, 0, iwfs, 1);
								wfsgrad_focus_drift(P(simu->gradoff, iwfs), simu, 0, iwfs, 1);
							} else if(ptype2==PTYPE_MF){
								if(jwfs==0) dbg("in cmf mode, gradoff is reset to 0, and ncpa is used to create i0 with new sodium profile\n");
								//since we are building ideal matched filter with
								//the correct gradoff and sodium profile. no need
								//to use gradient reference vector.
								dzero(P(simu->gradoff, iwfs));
							}
						}
						dcellfree(grad);
						dfree(sodium);
					}else{//NGS, just use i0
						if(ptype2==PTYPE_MF){
							if(simu->gradoff) {
								for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
									int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
									dzero(P(simu->gradoff, iwfs));
								}
							}
						}else{
							warning("Not handled: powfs%d: dither=-1, phytype_sim2=%d\n", ipowfs, ptype2);
						}
					}
				}

				//there is no need to reset trombone error signal
				if((parms->save.gradoff||parms->save.dither)&&parms->dbg.gradoff_reset!=1){
					writebin(simu->gradoff, "extra/gradoff_%d_dither", isim);
				}
				if(parms->dbg.gradoff_reset==2){
					dcellzero(simu->gradoffacc);
					simu->gradoffisim0=isim;
				}
				if(!simu->gradoffdrift&&parms->powfs[ipowfs].dither==1){
					//Use gradoff before adjustment is not good. There are difference between i and i0.
					if(parms->powfs[ipowfs].phytype_sim==PTYPE_MF){
						info2("Step %5d: powfs%d set gradoffdrift to cog of initial i0\n", isim, ipowfs);
					} else{
						info2("Step %5d: powfs%d set gradoffdrift to cog of created i0 + gradoff\n", isim, ipowfs);
					}
					simu->gradoffdrift=dcellnew(nwfs, 1);
					for(int jwfs=0; jwfs<nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
						shwfs_grad(&P(simu->gradoffdrift, iwfs), PCOLR(powfs[ipowfs].intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);
						if(parms->powfs[ipowfs].phytype_sim!=PTYPE_MF){
							dadd(&P(simu->gradoffdrift, iwfs), 1, P(simu->gradoff, iwfs), 1);//added on 2021-08-11
						}
					}
				}

				if(parms->powfs[ipowfs].phytype_sim!=ptype2){
					//the following parms changes need to be moved to simu. It affects the next seed.
					parms->powfs[ipowfs].phytype_sim=ptype2;
					parms->powfs[ipowfs].phytype_recon=parms->powfs[ipowfs].phytype_sim;
					info2("Step %5d: powfs %d changed to %s\n", isim, ipowfs,
						parms->powfs[ipowfs].phytype_sim==PTYPE_MF?"matched filter":"CoG");
				}
				//Generating matched filter
				if(parms->powfs[ipowfs].neareconfile||parms->powfs[ipowfs].phyusenea){
					dbg("Step %5d: powfs %d disable neareconfile and phyusenea\n", isim, ipowfs);
					parms->powfs[ipowfs].neareconfile=NULL;
					parms->powfs[ipowfs].phyusenea=0;
				}
				if(ptype2==PTYPE_MF){
					parms->powfs[ipowfs].phytype_recon=PTYPE_MF;//Make sure nea is used for reconstruction.
					mtch_wrap(parms, powfs, ipowfs);
#if USE_CUDA
					if(parms->gpu.wfs){
						gpu_wfsgrad_update_ref(parms, powfs, ipowfs);
					}
#endif
				}
				if(parms->save.dither>1){
					writebin(intstat->i0, "extra/powfs%d_i0_%d", ipowfs, isim);
					if(ptype2==PTYPE_MF){
						writebin(intstat->gx, "extra/powfs%d_gx_%d", ipowfs, isim);
						writebin(intstat->gy, "extra/powfs%d_gy_%d", ipowfs, isim);
						writebin(intstat->mtche, "extra/powfs%d_mtche_%d", ipowfs, isim);
						writebin(powfs[ipowfs].sanea, "extra/powfs%d_sanea_%d", ipowfs, isim);
					}
				}

				if(!parms->powfs[ipowfs].lo&&parms->recon.alg==RECON_MVR){//no need to update LSR.
					simu->tomo_update=2;
				}
				if(parms->powfs[ipowfs].dither_gdrift>0&&parms->powfs[ipowfs].dither==1){
					wfsgrad_sa_drift(simu, ipowfs);
				}
			} else {
				if(parms->powfs[ipowfs].phytype_sim!=ptype2){
					error("Does not support switching to CoG.\n");
				}
				//For CoG gain
				for(int jwfs=0; jwfs<nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dither_t* pd=simu->dither[iwfs];
					const int ogsingle=!pd->gg0||parms->powfs[ipowfs].dither_ogsingle;
					const int ng=ogsingle?1:(powfs[ipowfs].saloc->nloc*parms->powfs[ipowfs].ng);
					const int nd=PN(pd->a2mv);
					if(!P(simu->gradscale, iwfs)){
						P(simu->gradscale, iwfs)=dnew(ng, 1);
						dset(P(simu->gradscale, iwfs), parms->powfs[ipowfs].gradscale);
					}
					if(nd>1 && !P(simu->gradscale2, iwfs)){
						P(simu->gradscale2, iwfs)=dnew(nd, 1);
						dset(P(simu->gradscale2, iwfs), 1);
					}
					real mgold=dsum(P(simu->gradscale, iwfs))/ng;
					real mgnew=0;//,mgnew2=0;
					const char* ogtype=0;
					//gg0 is output/input of dither dithersig.
					if(ogsingle){//single gain for all subapertures. For Pyramid WFS
						ogtype="globally";
						real gerr=pd->a2m/pd->a2me;
#define HIA_G_UPDATA 0
#if HIA_G_UPDATE //HIA method.
						real adj=parms->powfs[ipowfs].dither_gog*mgold*(gerr-1);
						if(fabs(adj)>0.1) adj*=0.1/fabs(adj);
						while(adj+mgold<0){//prevent negative value
							adj*=0.5;
						}
						dadds(P(simu->gradscale, iwfs), adj);
						mgnew=mgold+adj;
#else
						real adj=pow(gerr, parms->powfs[ipowfs].dither_gog);
						CLIP(adj, 0.7, 1.5);//0.5, 1.5 //clip the adjustment to prevent divergence.
						dscale(P(simu->gradscale, iwfs), adj);
						mgnew=mgold*adj;
#endif
						if(nd>1){//multi-mode dithering
							dmat *gs2=P(simu->gradscale2, iwfs);
							for(int id=0; id<nd; id++){
								const real g2err=P(pd->a2mv, id)/P(pd->a2mev, id);
								/**
								 * A low pass filter is used to update gs2
								 * instead of an integrator like gradscale
								 * because the gradient and a2mev is not
								 * affected by its value.
								*/
								const real g=parms->powfs[ipowfs].dither_gog;
#if HIA_G_UPDATE //HIA method.
									P(gs2, id)+=(g2err-P(gs2, id))*g;
#else
									P(gs2, id)*=pow(g2err/P(gs2, id), g);
#endif
							}
							//mgnew2=P(gs2, nd-1);
							//dshow(gs2, "gradscale2");
							//info2("Scale ratios are High: %g Low: %g\n", P(simu->gradscale2, iwfs), P(P(simu->gradscale, iwfs), 0));
						}

					} else{//separate gain for each gradient. For shwfs.
						ogtype="on average";
						dscale(pd->gg0, scale1); //Scale value at end of accumulation
						for(long ig=0; ig<ng; ig++){
							if(P(pd->gg0, ig)>0.01){//skip weakly determined subapertures.
#if HIA_G_UPDATE //HIA method.
								real adj=parms->powfs[ipowfs].dither_gog*mgold*(1.-P(pd->gg0, ig));
								while(adj+P(P(simu->gradscale, iwfs), ig)<0){
									adj*=0.5;
								}
								P(P(simu->gradscale, iwfs), ig)+=adj;
#else
								real adj=pow(P(pd->gg0, ig), -parms->powfs[ipowfs].dither_gog);
								if(adj>1.5) adj=1.5; else if(adj<0.5) adj=0.5;
								P(P(simu->gradscale, iwfs), ig)*=adj;
#endif
							}
						}
						mgnew=dsum(P(simu->gradscale, iwfs))/ng;
						dzero(pd->gg0);
					}
					/*if(mgnew2){
						info2("Step %5d: wfs %d estimate/dither=%.2f, updated CoG gain=(%5.2f, %5.2f) %s\n",
							isim, iwfs, pd->a2me/pd->a2m, mgnew, mgnew2, ogtype);
					}else*/{
						info2("Step %5d: wfs %d estimate/dither=%.2f, updated CoG gain=%5.2f %s\n",
							isim, iwfs, pd->a2me/pd->a2m, mgnew, ogtype);
					}
					if(simu->resdither){
						int ic=simu->wfsflags[ipowfs].pllout-1;
						P(P(simu->resdither, iwfs), 3, ic)=mgnew;
					}
					//adjust WFS measurement dither dithersig by gain adjustment. used for dither t/t removal from gradients.
					pd->a2me*=(mgnew/mgold);//Adjust for updated gain
					dcellscale(powfs[ipowfs].sanea, pow(mgnew/mgold, 2));
					if(parms->save.dither){
						writebin(P(simu->gradscale, iwfs), "extra/gradscale_wfs%d_%d", iwfs, isim);
						if(P(simu->gradscale2, iwfs)){
							writebin(P(simu->gradscale2, iwfs), "extra/gradscale2_wfs%d_%d", iwfs, isim);
						}
					}
				}
			}
		}
	}
}
struct fit_cache{
	dccell *sepsf;
	dccell *i0m;
	dcell *grad;
	dcell *i0mv;
}fit_cache={0};
void fit_cache_free(){
	cellfree(fit_cache.sepsf);
	cellfree(fit_cache.i0m);
	cellfree(fit_cache.i0mv);
}

/**
 * Fit i0 to sodium profile using iterative algorithm. The steps are as follows
 * 1. Create sub images for each sodium profile bin
 * 2. Fit such sub-images against i0 to determine the profile
 * 3. Determine subaperture tip/tilt comparing i0 and fitted i0 (using matched filter)
 * 4. Repeat 1-4.
 * */
void sodium_fit(
	dmat** sodium, /**<The sodium profile determined by fit*/
	dcell** pgrad, /**<The gradients determined by fit.*/
	dcell** pi0,   /**<The output i0*/
	dcell** pgx,   /**<The output gx*/
	dcell** pgy,   /**<The output gy*/
	const dcell* i0i, /**<The input i0*/
	const dccell* sepsf,   /**<Short exposure PSF*/
	const dtf_t* dtf,     /**<Detector transfer function*/
	const loc_t* saloc,   /**<Saloc*/
	const dcell* saa,      /**<Subaperture area. */
	const dcell* srsa,    /**<Subaperture to LLT distance*/
	const dcell* srot,    /**<Subaperture to LLT clocking*/
	const dmat* siglevs,  /**<Subaperture signal level*/
	const dmat* wvlwts,    /**<Wavelength weights*/
	const dcell* gradoff,/**<NCPA gradient to be used for pi0,pgx,pgy output.*/
	real dh,      /**<The sodium profile sampling in meters*/
	real hs,      /**<LGS focusing height*/
	real htel,    /**<Telescope hegith*/
	real za,      /**<Telescope zenith angle*/
	real svdthres, /**<SVD threshold*/
	int nrep,     /**<Number of iterations*/
	int save,      /**<Save results to file*/
	int use_cache  /**<Use cache*/
){
	static int count=-1; count++;
	//const real ht=25000;//total thickness
	const real hmin=80000;
	const real hmax=105000;//wrapp around happens at 10000 for 10 pixel if alined along x/y.

	long nh=(long)floor((hmax-hmin)/dh)+1;
	const int radgx=0;
	dmat *nai=(sodium&&*sodium)?*sodium:NULL;
	if(!nai){
		nai=dnew(nh, 2);
		if(sodium) *sodium=nai;
	}else if(NX(nai)!=nh || NY(nai)!=2){
		dresize(nai, nh, 2);
	}

	for(long ix=0; ix<nh; ix++){
		P(nai, ix, 0)=hmin+ix*dh;
	}
	const long nsa=NX(i0i);
	const long ni0=NY(i0i);
	real pixthetax=dtf[0].pixthetax;
	real pixthetay=dtf[0].pixthetay;

	dcell* i0tmp=0;
	dcell* gxtmp=0;
	dcell* gytmp=0;
	dcell* gradtmp=0;
	
	//avoid overriding i0 input. Array for gensei full output to build temporal matched filter
	dcell **pi0tmp=(pi0 && (nrep==1 || *pi0 != i0i))?pi0:&i0tmp;
	dcell **pgxtmp=pgx?pgx:&gxtmp;//ok to override input
	dcell **pgytmp=pgy?pgy:&gytmp;
	if(!pgrad) pgrad=&gradtmp;
	if(!*pgrad && (nrep>1 || pgrad!=&gradtmp)){//need output
		*pgrad=dcellnew_same(ni0, 1, nsa*2, 1);
	}
	print_mem("grad");
	TIC;tic;
	dcell* i02=dcellnew(1, 1);
	if(!i0i->m){
		error("i0i->m is not set\n");
	}
	P(i02, 0)=dref(i0i->m);
	dcell* mtche=0;
	dcell* res=0;
	//etf_t** etfs=use_cache?fit_cache.etfs:NULL;
	dccell* i0m=use_cache?fit_cache.i0m:NULL;//sa image for each sodium bin
	dcell* i0mv=use_cache?fit_cache.i0mv:NULL;//vectorized i0m
	dcell *grad=use_cache?fit_cache.grad:*pgrad;
	int skip_first=0;

	if(!i0m){//prepare to create the sub layer image model
		i0mv=dcellnew(1, nh);
		i0m=dccellnew(nh, 1);
		if(gradoff){
			dcellcp(&grad, gradoff);
		}else{
			grad=dcellnew_same(ni0, 1, nsa*2, 1);
		}
		dbg("Initial gradient uses %s.\n", gradoff?"gradoff":"zero");
		
		if(use_cache){
			fit_cache.i0m=i0m;
			fit_cache.i0mv=i0mv;
			fit_cache.grad=grad;
		}
		if(save){
			writebin(grad, "sodium_grad_in");
		}
	}else{
		skip_first=1;//cache available.
		dbg("reuse previous grad, i0m, i0mv\n");
	}

	dbg("svdthres=%g, nrep=%d\n", svdthres, nrep);
	//don't try to cachc fotf. It is per WFS and uses too much storage.
	dcell* ata=0, * atb=0;
	etf_t *etf_full=0;
	dcell *na2s=dcellnew_same(nh, 1, 1, 2);
	for(int irep=0; irep<nrep; irep++){
		dbg("repeat %d of %d\n", irep+1, nrep);
		if(irep>0 || !skip_first){//Compute subaperture sublayer imaging model
			if(irep>0 && &grad!=pgrad){
				dcellcp(&grad, *pgrad);
			}
			OMP_FOR(NTHREAD)
			for(long ix=0; ix<nh; ix++){
				dmat *na2i=P(na2s, ix);
				P(na2i, 0, 0)=P(nai, ix, 0);
				P(na2i, 0, 1)=1;
				//ETF takes a lot of storage but is inexpensive to build. So we choose to build it on the fly
				etf_t *etf_i=mketf(dtf, na2i, 0, srot, srsa, hs, htel, za, 1);
				gensei(&P(i0m, ix), NULL, NULL, NULL, sepsf, dtf, etf_i, saa, radgx?srot:NULL, siglevs, wvlwts, grad, 0, 0);
				etf_free(etf_i);
				if(!P(i0mv, 0, ix)||P(P(i0mv, 0, ix))!=P(P(i0m, ix)->m)){
					dfree(P(i0mv, 0, ix));
					P(i0mv, 0, ix)=dref(P(i0m, ix)->m);
				}

			}
			toc2("gensei each");tic;
			print_mem("gensei");
		}
		dcellzero(ata);
		dcellzero(atb);
		dcellmm(&ata, i0mv, i0mv, "tn", 1);
		dcellmm(&atb, i0mv, i02, "tn", 1);
		if(save){
			writebin(ata, "sodium_ata_%d_%d", count, irep);
			writebin(atb, "sodium_atb_%d_%d", count, irep);
			if(count==0&&irep==0) writebin(i0m, "sodium_i0m_%d_%d", count, irep);
		}
		dcellsvd_pow2(ata, -1, svdthres, 1e-3);
		dcellzero(res);
		dcellmm(&res, ata, atb, "nn", 1);
		real scale=1./dcellsum(res);//make sure nai sum to 1.
		for(long ix=0; ix<nh; ix++){
			P(nai, ix, 1)=P(P(res, ix), 0)*scale;
		}
		if(nrep>1 || pgrad!=&gradtmp){//need to determine error in applying gradient offset
			if(etf_full) etf_free(etf_full);
			//mketf for full profile must use the same no_interp flag 
			etf_full=mketf(dtf, nai, 0, srot, srsa, hs, htel, za, 1);
			toc2("mketf full"); tic;
			gensei(pi0tmp, pgxtmp, pgytmp, NULL, sepsf, dtf, etf_full, saa, radgx?srot:NULL, siglevs, wvlwts, grad, 0, 0);
			toc2("gensei full"); tic;
			mtch_cell(&mtche, NULL, NULL, NULL, *pi0tmp, *pgxtmp, *pgytmp, NULL, NULL, NULL, 0, 0, 3,
				pixthetax, pixthetay, NULL, radgx, 1, 1);
			toc2("mtche create"); tic;

			for(long ii0=0; ii0<ni0; ii0++){
				dmat* grad1=P(grad, ii0);//model
				dmat* grad2=P(*pgrad, ii0);//output
OMP_FOR(NTHREAD)
				for(long isa=0; isa<nsa; isa++){
					real g[2]={0,0};
					dmulvec(g, P(mtche, isa, ii0), P(P(i0i, isa, ii0)), 1);
					P(grad2, isa)    =P(grad1, isa    )+g[0];
					P(grad2, isa+nsa)=P(grad1, isa+nsa)+g[1];
				}
			}
			toc2("mtche apply"); tic;
			
OMP_FOR(NTHREAD)
			for(long ii0=0; ii0<ni0; ii0++){
				//Remove focus mode from the gradients as it degenerates with sodium profile shift.
				loc_remove_focus_grad(P(*pgrad, ii0), saloc, 1);
			}
			if(save){
				writebin(nai, "sodium_prof_%d_%d", count, irep);
				writebin(*pgrad, "sodium_grad_%d_%d", count, irep);
				//writebin(mtche, "sodium_mtche_%d_%d", count, irep);
				//writebin(*pi0tmp, "sodium_i0_%d_%d", count, irep);
			}
		}
	}

	//output is desired. build final i0, gx, gy with the final gradient or ncpa gradient.
	if(pi0 || pgx || pgy){
		dbg("Replacing i0, gx, gy with fitted value\n");
		if(!etf_full){
			etf_full=mketf(dtf, nai, 0, srot, srsa, hs, htel, za, 1);
			toc2("mketf final");tic;
		}
		const dcell *gradf=gradoff?gradoff:(pgrad?(*pgrad):grad);
		gensei(pi0, pgx, pgy, NULL, sepsf, dtf, etf_full, saa, radgx?srot:NULL, siglevs, wvlwts, gradf, 0, 0);
		toc2("gensei final");tic;
	}

	if(etf_full) etf_free(etf_full);
	dcellfree(ata);
	dcellfree(atb);
	if(!use_cache){
		cellfree(i0m);
		cellfree(i0mv);
		cellfree(grad);
	}
	cellfree(mtche);
	cellfree(i02);
	cellfree(res);
	if(!sodium) dfree(nai);
	cellfree(na2s);
	cellfree(i0tmp);
	cellfree(gxtmp);
	cellfree(gytmp);
	cellfree(gradtmp);
}

/**
 * Fit i0 to sodium profile and replace i0, gx, gy with derived parameters
 * */
void sodium_fit_wrap(dmat** psodium, /**<[out] sodium profile*/
	dcell** pgrad, /**<[out] estimated actual gradient*/
	dcell** pi0,   /**<[out] The output i0*/
	dcell** pgx,   /**<[out] The output gx*/
	dcell** pgy,   /**<[out] The output gy*/
	const dcell* i0in, /**<[in]The input sa intensities. may equal to *pi0 */
	const parms_t* parms,/**<[in]parms*/
	powfs_t* powfs, /**<[in]powfs*/
	int ipowfs, /**<[in] ipowfs*/
	real r0,  /**<[in] Fried parameter*/
	real L0,  /**<[in] outer scale*/
	int nrep, /**<[in] Number of iterations. 1 for mtche, 3 for cog*/
	int use_cache /**<[in] cache intermediate results.*/
	){
	dccell* sepsf=use_cache?fit_cache.sepsf:NULL;
	if(!sepsf){
		cccell* otf=NULL, * lotf=NULL;
		otf=genseotf(powfs[ipowfs].pts, powfs[ipowfs].amp,
			NULL, powfs[ipowfs].saa, parms->powfs[ipowfs].wvl, r0, L0,
			parms->powfs[ipowfs].embfac);
		if(parms->powfs[ipowfs].llt){
						//genselotf(parms, powfs, ipowfs);
			lotf=genseotf(powfs[ipowfs].llt->pts, powfs[ipowfs].llt->amp,
				NULL, NULL, parms->powfs[ipowfs].wvl, r0, L0, 
				parms->powfs[ipowfs].embfac);
		}
		gensepsf(&sepsf, otf, lotf, powfs[ipowfs].saa,
			parms->powfs[ipowfs].wvl, powfs[ipowfs].notfx, powfs[ipowfs].notfy);
		cellfree(otf);
		cellfree(lotf);
		if(use_cache){
			fit_cache.sepsf=sepsf;
		}
		if(parms->save.dither){
			static int count=0;
			writebin(sepsf, "sodium_sepsf_%d", count);
			count++;
		}
	}

	sodium_fit(psodium, pgrad, pi0, pgx, pgy, i0in,
		sepsf, powfs[ipowfs].dtf, powfs[ipowfs].saloc, powfs[ipowfs].saa,
		powfs[ipowfs].srsa, powfs[ipowfs].srot,
		parms->powfs[ipowfs].siglevs, parms->powfs[ipowfs].wvlwts, powfs[ipowfs].gradoff,
		parms->powfs[ipowfs].llt->na_fit_dh, parms->powfs[ipowfs].hs, parms->sim.htel, parms->sim.za, 
		parms->powfs[ipowfs].llt->na_fit_svdthres, nrep, parms->save.dither>1, use_cache);//parms->save.setup);
	
	/*if(parms->save.setup){
		if(pi0) writebin(*pi0, "powfs%d_i0_fit", ipowfs);
		if(pgx) writebin(*pgx, "powfs%d_gx_fit", ipowfs);
		if(pgy) writebin(*pgy, "powfs%d_gy_fit", ipowfs);
		if(psodium) writebin(*psodium, "powfs%d_sodium_fit", ipowfs);
	}*/
	if(!use_cache){
		dcellfree(sepsf);
	}else{
		static int registered=0;
		if(!registered){
			registered=1;
			register_deinit(fit_cache_free, NULL);
		}
	}
}
