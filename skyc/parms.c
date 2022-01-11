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
#include <unistd.h>

#include "skyc.h"
#include "parms.h"
#include "nafocus.h"
#undef malloc
#undef calloc
#undef free
#undef realloc

#define READ_INT(A) parms->A = readcfg_int(#A) /*read a key with int value. */
#define READ_DBL(A) parms->A = readcfg_dbl(#A) /*read a key with real value */
#define READ_STR(A) parms->A = readcfg_str(#A) /*read a key with string value. */
#define READ_MAT(A) parms->A = readcfg_dmat(#A) /*read a key with real array. */

static void setup_parms_skyc(PARMS_S* parms){
	READ_INT(skyc.dbg);
	READ_INT(skyc.dbgsky);
	READ_INT(skyc.dbgaster);
	READ_INT(skyc.keeporder);
	READ_INT(skyc.interpg);
	READ_INT(skyc.verbose);
	//READ_INT(skyc.reest);
	READ_INT(skyc.save);
	READ_INT(skyc.start);
	READ_INT(skyc.nsky);
	READ_DBL(skyc.lat);
	READ_DBL(skyc.lon);
	READ_DBL(skyc.catscl);
	READ_DBL(skyc.patfov);
	READ_DBL(skyc.minrad);
	READ_DBL(skyc.keepout);
	READ_INT(skyc.ngsalign);
	READ_INT(skyc.npowfs);
	READ_INT(skyc.psd_scale);
	READ_INT(skyc.noisy);
	READ_INT(skyc.ttfbrightest);
	READ_INT(skyc.ttffastest);
	READ_INT(skyc.bspstrehl);
	READ_INT(skyc.maxaster);
	READ_INT(skyc.maxdtrat);
	READ_INT(skyc.maxstar);
	readcfg_intarr_n(&parms->skyc.nwfsmax, parms->skyc.npowfs, "skyc.nwfsmax");
	readcfg_dblarr_n(&parms->skyc.pixtheta, parms->skyc.npowfs, "skyc.pixtheta");
	readcfg_dblarr_n(&parms->skyc.pixblur, parms->skyc.npowfs, "skyc.pixblur");
	readcfg_intarr_n(&parms->skyc.pixpsa, parms->skyc.npowfs, "skyc.pixpsa");
	readcfg_intarr_n(&parms->skyc.pixguard, parms->skyc.npowfs, "skyc.pixguard");
	readcfg_dblarr_n(&parms->skyc.pixoffx, parms->skyc.npowfs, "skyc.pixoffx");
	readcfg_dblarr_n(&parms->skyc.pixoffy, parms->skyc.npowfs, "skyc.pixoffy");
	readcfg_strarr_nmax(&parms->skyc.fnpsf1, parms->skyc.npowfs, "skyc.fnpsf1");
	READ_INT(skyc.limitnstep);
	READ_DBL(skyc.rne);
	READ_DBL(skyc.excess);
	READ_DBL(skyc.imperrnm);
	READ_DBL(skyc.imperrnmb);
	READ_INT(skyc.mtchcr);
	READ_INT(skyc.mtchfft);
	READ_INT(skyc.phytype);
	READ_MAT(skyc.zb.wvl);
	READ_MAT(skyc.zb.qe);
	READ_MAT(skyc.zb.thruput);
	READ_MAT(skyc.zb.excessbkgrnd);
	READ_MAT(skyc.zb.Z);
	READ_MAT(skyc.zb.B);
	/*readcfg_dblarr_nmax(&parms->skyc.zb.qe, parms->maos.nwvl, "skyc.zb.qe");
	readcfg_dblarr_nmax(&parms->skyc.zb.thruput, parms->maos.nwvl, "skyc.zb.thruput");
	readcfg_dblarr_nmax(&parms->skyc.zb.excessbkgrnd, parms->maos.nwvl, "skyc.zb.excessbkgrnd");*/
	parms->skyc.dtrats=readcfg_dmat("skyc.dtrats");
	parms->skyc.dtrats_mr=readcfg_dmat("skyc.dtrats_mr");
	parms->skyc.snrmin_mr=readcfg_dmat_nmax(PN(parms->skyc.dtrats_mr),"skyc.snrmin_mr");
	READ_INT(skyc.seed);
	READ_INT(skyc.navg);
	READ_INT(skyc.servo);
	READ_INT(skyc.gsplit);
	READ_INT(skyc.evlstart);
	READ_INT(skyc.phystart);
	READ_INT(skyc.neaaniso);
	READ_INT(skyc.neanonlin);
	char* temp;
	temp=readcfg_str("skyc.psd_ws");
	if(temp){
		parms->skyc.psd_ws=dread("%s", temp);
	}
	free(temp);

	READ_STR(skyc.stars);
	READ_INT(skyc.addws);
	READ_DBL(skyc.pmargin);
	READ_INT(skyc.psdcalc);
	READ_DBL(skyc.sdetmax);

	READ_INT(skyc.multirate);
	READ_DBL(skyc.snrmin);
	READ_INT(skyc.usephygrad);

	READ_INT(skyc.estimate);
}
/**
   Setup infromation output from maos run.
*/
static void setup_parms_maos(PARMS_S* parms){
	READ_DBL(maos.r0z);
	READ_DBL(maos.dt);
	READ_DBL(maos.zadeg); parms->maos.za=parms->maos.zadeg*M_PI/180.;
	READ_DBL(maos.hc);
	READ_DBL(maos.hs);
	READ_DBL(maos.D);
	READ_INT(maos.nmod);
	parms->maos.nseed=readcfg_intarr(&parms->maos.seeds, "maos.seeds");

	parms->maos.nwvl=readcfg_dblarr(&parms->maos.wvl, "maos.wvl");
	READ_INT(maos.npowfs);
	readcfg_intarr_n(&parms->maos.msa, parms->maos.npowfs, "maos.msa");
	readcfg_intarr_n(&parms->maos.nsa, parms->maos.npowfs, "maos.nsa");

	readcfg_intarr_n(&parms->maos.ncomp, parms->maos.npowfs, "maos.ncomp");
	readcfg_intarr_n(&parms->maos.embfac, parms->maos.npowfs, "maos.embfac");
	readcfg_dblarr_n(&parms->maos.dxsa, parms->maos.npowfs, "maos.dxsa");

	readcfg_strarr_n(&parms->maos.fnwfsloc, parms->maos.npowfs, "maos.fnwfsloc");
	readcfg_strarr_n(&parms->maos.fnwfsamp, parms->maos.npowfs, "maos.fnwfsamp");
	readcfg_strarr_n(&parms->maos.fnsaloc, parms->maos.npowfs, "maos.fnsaloc");
	char* temp;
	temp=readcfg_str("maos.fnmcc");
	parms->maos.mcc=dread("%s", temp); free(temp);
	parms->maos.mcc_tt=dsub(parms->maos.mcc, 0, 2, 0, 2);

	temp=readcfg_str("maos.fnmcc_oa");
	parms->maos.mcc_oa=dread("%s", temp); free(temp);
	parms->maos.mcc_oa_tt=dsub(parms->maos.mcc_oa, 0, 2, 0, 2);
	parms->maos.mcc_oa_tt2=dsub(parms->maos.mcc_oa, 0, 2, 0, 0);
	if(parms->maos.nmod>5&&parms->maos.mcc->nx<=5){
		warning("Resizing mcc (compatible mode)\n");
		dresize(parms->maos.mcc, parms->maos.nmod, parms->maos.nmod);
		dresize(parms->maos.mcc_oa, parms->maos.nmod, parms->maos.nmod);
	}

	READ_STR(maos.fnmideal);
	READ_STR(maos.fnmidealp);
	READ_INT(maos.evlindoa);
	READ_DBL(maos.ngsgrid);
	READ_INT(maos.nstep);
	READ_INT(maos.ahstfocus);
	READ_INT(maos.mffocus);
	READ_STR(maos.fnrange);
	info("maos.ahstofocus=%d, maos.mffocus=%d\n", parms->maos.ahstfocus, parms->maos.mffocus);

	if(readcfg_peek("maos.wddeg")){
		parms->maos.nwddeg=readcfg_dblarr(&parms->maos.wddeg, "maos.wddeg");
	} else{
		parms->maos.nwddeg=0;
		parms->maos.wddeg=NULL;
	}
	if(readcfg_peek("maos.indps")){
		READ_INT(maos.indps);
		READ_INT(maos.indastig);
		READ_INT(maos.indfocus);
	} else{//old configuration.
		parms->maos.indastig=0;//not implemented
		if(parms->maos.nmod>=5){
			parms->maos.indps=2;
			if(parms->maos.nmod>5){
				parms->maos.indfocus=5;
			}
		} else{
			parms->maos.indps=0;
			if(parms->maos.nmod>2){
				parms->maos.indfocus=2;
			}
		}
	}
}

PARMS_S* setup_parms(const ARG_S* arg){
	/*Setup PATH and result directory */
	char* config_path=find_config("skyc");
	addpath(config_path);
	free(config_path);
	open_config(arg->conf, NULL, 0);
	open_config(arg->confcmd, NULL, 1);
	remove(arg->confcmd);
	PARMS_S* parms=mycalloc(1, PARMS_S);
	parms->skyc.nthread=arg->nthread;
	setup_parms_maos(parms);
	setup_parms_skyc(parms);
	{
		real sum=0;
		real wtsum=0;
		real wt=0;
		parms->skyc.wvlwt=dnew(parms->maos.nwvl, 1);
		for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
			for(int jwvl=0; jwvl<NX(parms->skyc.zb.qe); jwvl++){
				if(fabs(parms->maos.wvl[iwvl]-P(parms->skyc.zb.wvl, jwvl))<0.1e-6){
					wt=P(parms->skyc.zb.qe,iwvl)*P(parms->skyc.zb.thruput,iwvl);
					break;
				}
			}
			P(parms->skyc.wvlwt, iwvl)=wt;
			sum+=parms->maos.wvl[iwvl]*wt;
			wtsum+=wt;
		}
		parms->skyc.wvlmean=sum/wtsum;
		dnormalize_sumabs(P(parms->skyc.wvlwt), parms->maos.nwvl, 1);
	}
	if(parms->skyc.maxdtrat<=0){
		parms->skyc.maxdtrat=parms->skyc.ndtrat;
	}
	if(parms->maos.ahstfocus){
		if(parms->skyc.addws==-1){//auto
			parms->skyc.addws=1;
		}
		if(parms->maos.nmod<=5){
			error("Conflicted parameters: maos.nmod should be >5 when maos.ahstfocus=1\n");
		}
	}
	if(parms->skyc.servo<0){
		parms->skyc.addws=1;
	}
	if(parms->skyc.multirate){
		dfree(parms->skyc.dtrats);
		parms->skyc.dtrats=dref(parms->skyc.dtrats_mr);
		parms->skyc.snrmin=dmin(parms->skyc.snrmin_mr);
	}
	dfree(parms->skyc.dtrats_mr);
	parms->skyc.ndtrat=parms->skyc.dtrats->nx*parms->skyc.dtrats->ny;
	if(parms->skyc.addws==-1){
		parms->skyc.addws=1;
	}
	if(!parms->skyc.psd_ws){
		parms->skyc.addws=0;
	}
	switch(parms->skyc.servo){
	case -1://LQG
		parms->skyc.ngain=0;
		break;
	case 1:
		parms->skyc.ngain=1;
		break;
	case 2:
		parms->skyc.ngain=3;
		break;
	default:
		error("Invalid skyc.servo=%d\n", parms->skyc.servo);
	}
	dbg("maos.mffocus=%d\n", parms->maos.mffocus);
	dbg("skyc.addws=%d\n", parms->skyc.addws);
	if(!parms->skyc.stars){
		if(parms->skyc.keeporder){
			error("When skyc.keeporder is set, skyc.stars need to be set\n");
		}
	} else{
		if(parms->skyc.interpg==1){
			info("disable interpg when skyc.stars is set\n");
			parms->skyc.interpg=0;
		}
	}


	for(int ipowfs=0; ipowfs<parms->skyc.npowfs; ipowfs++){
		parms->skyc.pixtheta[ipowfs]/=206265.;//input is in arcsec
		info("powfs %d, pixtheta=%g mas\n", ipowfs, parms->skyc.pixtheta[ipowfs]*206265000);
	}
	for(int idtrat=1; idtrat<parms->skyc.ndtrat; idtrat++){
		if(P(parms->skyc.dtrats,idtrat)>=P(parms->skyc.dtrats,idtrat-1)){
			error("skyc.dtrats must be specified in descending order\n");
		}
	}
	parms->skyc.fss=dnew(parms->skyc.ndtrat, 1);
	for(long idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
		int dtrat=P(parms->skyc.dtrats,idtrat);
		P(parms->skyc.fss,idtrat)=1./(parms->maos.dt*dtrat);
	}
	parms->skyc.rnefs=dnew(parms->skyc.ndtrat, parms->maos.npowfs);
	if(parms->skyc.rne<0){
	/* Uses the new noise model based on Roger's spread sheet and document
	   TMT.AOS.TEC.13.009.DRF01 H2RG Noise Model */
		dmat* rnefs=parms->skyc.rnefs;
		info("Using frame rate dependent read out noise:\n");
		if(fabs(parms->skyc.rne+1)<EPS){
			const real pixeltime=6.04;//time to read a pixel in us
			const real linetime=2;//time to read a line in us
			const real frametime=3;//time to read a frame in us

			for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
				int N=parms->skyc.pixpsa[ipowfs];
				int Nb=parms->skyc.pixguard[ipowfs];
				int nsa=parms->maos.nsa[ipowfs];
				real t1=nsa*(pixeltime*N*N+linetime*N+frametime);
				real t2=(pixeltime*Nb*Nb+linetime*Nb+frametime);
				for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
					int dtrat=P(parms->skyc.dtrats,idtrat);
					real dt=parms->maos.dt*dtrat*1e6;
					int coadd=floor((dt-t2*nsa)/t1);//number of coadds possible.
					if(coadd<=32&&dtrat<=10){//at high frame rate, read out only 1 subaperture's guard window each time.
						coadd=floor((dt-t2)/t1);
					}
					if(coadd>177){
						coadd=177;//more coadds increases dark noise.
					}
					if(coadd<1) coadd=1;
					real rne_white=10.9*pow(coadd, -0.47);
					real rne_floor=2.4;
					real rne_dark=0.0037*coadd;
					//0.85 is a factor expected for newer detectors
					P(rnefs, idtrat, ipowfs)=0.85*sqrt(rne_white*rne_white+rne_floor*rne_floor+rne_dark*rne_dark);
					info("powfs[%d] %5.1f Hz: %5.1f \n", ipowfs,
						P(parms->skyc.fss,idtrat), P(rnefs, idtrat, ipowfs));
				}
			}
		} else if(fabs(parms->skyc.rne+2)<EPS){//older model.
			for(long idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
				int dtrat=P(parms->skyc.dtrats,idtrat);
				real fs=1./(parms->maos.dt*dtrat);
				info("%5.1f Hz: ", fs);
				for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
					int N=parms->skyc.pixpsa[ipowfs];
					real raw_frame_time=((0.00604*N+0.01033)*N+0.00528)*1e-3;
					real K=1./(raw_frame_time*fs);
					P(rnefs, idtrat, ipowfs)=sqrt(pow(10.6*pow(K, -0.45), 2)+2.7*2.7+0.0034*K);
					info("%5.1f ", P(rnefs, idtrat, ipowfs));
				}
				info("\n");
			}
		} else{
			error("Invalid skyc.rne\n");
		}
	} else{
		info("Using constant read out noise of %g\n", parms->skyc.rne);
		dset(parms->skyc.rnefs, parms->skyc.rne);
	}

	if(parms->skyc.npowfs!=parms->maos.npowfs){
		error("skyc.npowfs should match maos.npowfs\n");
	}
	if(parms->skyc.ngsalign){
		warning("NGS are aligned to grid spaced by %g\"\n", parms->maos.ngsgrid);
	}

	if(parms->skyc.psdcalc){
		dbg("Calculating PSDs from time series\n");
	} else{
		error("Please implement PSD load\n");
	}

	{
		char fn[PATH_MAX];
		snprintf(fn, PATH_MAX, "skyc_%s_%ld.conf", HOST, (long)getpid());
		close_config("%s", fn);
		mysymlink(fn, "skyc_recent.conf");
	}

	if(parms->skyc.neaaniso){
		info("Variance of the gradients in stored PSF is added to NEA\n");
	}
	info("Maximum number of asterisms in each star field is %d\n", parms->skyc.maxaster);

	if(parms->skyc.dbg||parms->skyc.dbgsky>-1){
		dbg("skyc.dbg=%d, skyc.dbgsky=%d, disable multithreading\n", parms->skyc.dbg, parms->skyc.dbgsky);
		if(parms->skyc.verbose<1){
			parms->skyc.verbose=1;
		}
		parms->skyc.dbg=1;
		parms->skyc.nthread=1;
		parms->skyc.interpg=0;
	}
	if(parms->skyc.nsky<20){
		parms->skyc.interpg=0;
	}
	if(arg->detach){
		parms->skyc.verbose=0;
	}
	parms->skyc.nwfstot=0;
	for(int ipowfs=0; ipowfs<parms->skyc.npowfs; ipowfs++){
		parms->skyc.nwfstot+=parms->skyc.nwfsmax[ipowfs];
	}
	info("There are %d MAOS PSF seeds: ", parms->maos.nseed);
	for(int i=0; i<parms->maos.nseed; i++){
		info(" %d", parms->maos.seeds[i]);
	}
	info("\n");
	if(parms->maos.nseed){
		//skip unwanted seeds
		parms->fdlock=mycalloc(parms->maos.nseed, int);
		char fn[81];
		int nseed=0;
		for(int i=0; i<parms->maos.nseed; i++){
			long maos_seed=parms->maos.seeds[i];
			snprintf(fn, 80, "Res%ld_%d.done", maos_seed, parms->skyc.seed);
			if(exist(fn)&&!arg->override){
				parms->fdlock[i]=-1;
				warning("Will skip seed %ld because %s exist.\n", maos_seed, fn);
			} else{
				snprintf(fn, 80, "Res%ld_%d.lock", maos_seed, parms->skyc.seed);
				parms->fdlock[i]=lock_file(fn, 0, 0);
				if(parms->fdlock[i]<0){
					warning("Will skip seed %ld because it is already running.\n", maos_seed);
				} else{
					warning("Lock seed %ld success.\n", maos_seed);
					cloexec(parms->fdlock[i]);
					if(nseed!=i){
						parms->maos.seeds[nseed]=parms->maos.seeds[i];
						parms->fdlock[nseed]=parms->fdlock[i];
					}
					nseed++;
				}
			}
		}
		if(nseed!=parms->maos.nseed){
			info("Skip %d seeds.\n", parms->maos.nseed-nseed);
			parms->maos.nseed=nseed;
		}
		if(!nseed){
			warning("There are no seed to run. Use -O to override.\n");
			//remove log and conf files
			char fnpid[PATH_MAX];
			snprintf(fnpid, PATH_MAX, "run_%d.log", (int)getpid());
			remove(fnpid);
		}
	}
	return parms;
}
/**
   Free the data in parms.
*/
void free_parms(PARMS_S* parms){
	dfree(parms->skyc.psd_ngs);
	dfree(parms->skyc.psd_tt);
	dfree(parms->skyc.psd_ws);
	dfree(parms->skyc.psd_ps);
	dfree(parms->skyc.dtrats);
	free(parms->skyc.stars);parms->skyc.stars=NULL;
	dfree(parms->skyc.fss);
	dfree(parms->skyc.wvlwt);
	dfree(parms->skyc.rnefs);
	dfree(parms->maos.mcc);
	dfree(parms->maos.mcc_oa);
	dfree(parms->maos.mcc_tt);
	dfree(parms->maos.mcc_oa_tt);
	dfree(parms->maos.mcc_oa_tt2);
	free(parms->maos.msa);parms->maos.msa=NULL;
	free(parms->maos.nsa);parms->maos.nsa=NULL;
	free(parms->maos.ncomp);parms->maos.ncomp=NULL;
	free(parms->maos.embfac);parms->maos.embfac=NULL;
	free(parms->maos.dxsa);parms->maos.dxsa=NULL;
	free_strarr(parms->maos.fnwfsloc, parms->maos.npowfs);parms->maos.fnwfsloc=NULL;
	free_strarr(parms->maos.fnwfsamp, parms->maos.npowfs);parms->maos.fnwfsamp=NULL;
	free_strarr(parms->maos.fnsaloc, parms->maos.npowfs);parms->maos.fnsaloc=NULL;
	free(parms);
}
