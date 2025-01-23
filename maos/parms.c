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


#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include "parms.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   Read in configuration parameters and check for consistence.
*/
extern int use_cuda;
/*
  Don't include common.h or types.h, so that we don't have to recompile
  setup_parms.c when these files change.

  This file contains necessary routines to read parametes for
  WFS, DM and wavefront reconstruction.  */

void free_powfs_cfg(powfs_cfg_t *powfscfg){
	dfree(powfscfg->wvl);
	dfree(powfscfg->wvlwts);
	dfree(powfscfg->ncpa);
	if(powfscfg->llt){
		free(powfscfg->llt->fnrange);
		free(powfscfg->llt->fnprof);
		free(powfscfg->llt->fnprep);
		free(powfscfg->llt->fnamp);
		free(powfscfg->llt->fnsurf);
		lfree(powfscfg->llt->i);
		dfree(powfscfg->llt->ox);
		dfree(powfscfg->llt->oy);
		dfree(powfscfg->llt->misreg);
		free(powfscfg->llt);
		powfscfg->llt=NULL;
	}
	if(powfscfg->pywfs){
		pycfg_free(powfscfg->pycfg);
	}
	lfree(powfscfg->wfs);
	lfree(powfscfg->wfsr);
	lfree(powfscfg->wfsind);
	free(powfscfg->pywfs);
	free(powfscfg->fnllt);
	free(powfscfg->piinfile);
	free(powfscfg->sninfile);
	free(powfscfg->neareconfile);
	free(powfscfg->neasimfile);
	free(powfscfg->bkgrndfn);
	free(powfscfg->qe);
	dfree(powfscfg->siglevs);
}
void free_dm_cfg(dm_cfg_t *dmcfg){
	dfree(dmcfg->actstuck);
	dfree(dmcfg->actfloat);
	cellfree(dmcfg->strokescale);
	dfree(dmcfg->stroke);
}
/**
   Free the parms struct.
*/
void free_parms(parms_t *parms){
	dfree(parms->atm.ht);
	dfree(parms->atm.wt);
	dfree(parms->atm.ws);
	dfree(parms->atm.wddeg);
	dfree(parms->atm.L0);
	dfree(parms->atmr.ht);
	dfree(parms->atmr.wt);
	lfree(parms->atmr.os);
	dfree(parms->atm.size);
	lfree(parms->atm.ipsr);
	lfree(parms->atm.overx);
	lfree(parms->atm.overy);
	lfree(parms->atm.nxn);
	lfree(parms->atmr.indps);
	dfree(parms->atm.r0logpsdt);
	dfree(parms->atm.r0logpsds);
	dfree(parms->evl.thetax);
	dfree(parms->evl.thetay);
	dfree(parms->evl.wvl);
	free(parms->evl.wvlname);
	dfree(parms->evl.wt);
	dfree(parms->evl.hs);
	lfree(parms->evl.psf);
	lfree(parms->evl.psfr);
	lfree(parms->evl.psfgridsize);
	lfree(parms->evl.psfsize);
	lfree(parms->evl.pttr);

	dfree(parms->fit.thetax);
	dfree(parms->fit.thetay);
	dfree(parms->fit.wt);
	dfree(parms->fit.hs);

	dfree(parms->sim.aphi);
	dfree(parms->sim.ephi);
	dfree(parms->sim.aplo);
	dfree(parms->sim.eplo);

	lfree(parms->sim.seeds);
	dfree(parms->sim.wspsd);

	dfree(parms->ncpa.thetax);
	dfree(parms->ncpa.thetay);
	dfree(parms->ncpa.wt);
	dfree(parms->ncpa.hs);
	free(parms->sim.mvmhost);
	dfree(parms->cn2.pair);
	lfree(parms->save.gcov);
	for(int isurf=0; isurf<parms->ncpa.nsurf; isurf++){
		free(parms->ncpa.surf[isurf]);
	}
	free(parms->ncpa.surf); parms->ncpa.nsurf=0;
	for(int isurf=0; isurf<parms->ncpa.ntsurf; isurf++){
		free(parms->ncpa.tsurf[isurf]);
	}
	free(parms->ncpa.tsurf); parms->ncpa.ntsurf=0;

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		free_powfs_cfg(&parms->powfs[ipowfs]);
	}
	free(parms->powfs); parms->npowfs=0;
	free(parms->wfsr); parms->nwfsr=0;
	
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		dfree(parms->wfs[iwfs].wvlwts);
		dfree(parms->wfs[iwfs].sabad);
	}
	free(parms->wfs); parms->nwfs=0;

	for(int idm=0; idm<parms->ndm; idm++){
		free_dm_cfg(&parms->dm[idm]);
	}
	free(parms->dm); parms->ndm=0;
	for(int imoao=0; imoao<parms->nmoao; imoao++){
		dfree(parms->moao[imoao].actstuck);
		dfree(parms->moao[imoao].actfloat);
	}
	free(parms->moao); parms->nmoao=0;
	free(parms->aper.fnamp);
	free(parms->aper.pupmask);
	lfree(parms->save.ints);
	lfree(parms->save.wfsopd);
	lfree(parms->save.grad);
	lfree(parms->save.gradnf);
	lfree(parms->save.gradpsol);
	lfree(parms->save.gradgeom);
	free(parms->load.mvm);
	free(parms->load.mvmi);
	free(parms->load.mvmf);
	free(parms->load.ncpa);
	free(parms->load.saneai);
	free(parms->fdlock);
	free(parms->fnlock);
	lfree(parms->hipowfs);
	lfree(parms->lopowfs);

	dfree(parms->aper.misreg);
	free_strarr(parms->distortion.tel2wfs,parms->nwfs);
	free_strarr(parms->distortion.dm2wfs,parms->ndm*parms->nwfs);
	free_strarr(parms->distortion.dm2sci,parms->ndm*parms->evl.nevl);
	free_strarr(parms->recon.distortion_dm2wfs,parms->ndm*parms->nwfsr);
	free_strarr(parms->recon.distortion_dm2sci,parms->ndm*parms->fit.nfit);
	free_strarr(parms->recon.distortion_tel2wfs,parms->nwfsr);
	dfree(parms->dirs);
	lfree(parms->dbg.tomo_maxit);
	dcellfree(parms->dbg.dmoff);
	dcellfree(parms->dbg.gradoff);
	dfree(parms->dbg.atm);
	free(parms);
}
/*static inline int sum_intarr(int n, long *a){
	int sum=0;
	for(int i=0; i<n; i++){
	sum+=(a[i]!=0);
	}
	return sum;
}
static inline int sum_dblarr(int n, real *a){
	real sum=0;
	for(int i=0; i<n; i++){
	sum+=(a[i]!=0);
	}
	return sum;
	}*/

#define MAX_STRLEN 80
#define READ_INT(A) parms->A = readcfg_int(#A) /*read a key with int value. */
#define READ_DBL(A) parms->A = readcfg_dbl(#A) /*read a key with real value */
#define READ_DBL_SCALE(A,B,C) parms->A = readcfg_dbl(#B)*(C) /*read a key with real value */
#define READ_STR(A) parms->A = readcfg_str(#A) /*read a key with string value. */
#define READ_DMAT(A) parms->A= readcfg_dmat(0,0,#A) /*read a key with dmat. */
#define READ_DMAT_N(A,n) parms->A= readcfg_dmat(n,0,#A) /*read a key with dmat. */
#define READ_DMAT_NMAX(A,n) parms->A= readcfg_dmat(n,1,#A) /*read a key with dmat. */
#define READ_DCELL(A) parms->A= readcfg_dcell(#A) /*read a key with dmat. */
#define READ_LMAT(A) parms->A= readcfg_lmat(0,0,#A) /*read a key with lmat. */
#define READ_LMAT_NMAX(A,n) parms->A= readcfg_lmat(n,1,#A) /*read a key with lmat. */

#define READ_POWFS(A,B)						\
    if(readcfg_##A##arr((&A##tmp), npowfs,0, "powfs."#B)==npowfs){\
    for(int i=0; i<npowfs; i++){					\
		parms->powfs[i].B = A##tmp[i];/*doesn't need ## in B*/	\
    }}//else{dbg("Empty array for powfs." #B"\n");}
#define READ_POWFS_RELAX(A,B)					\
    if(readcfg_##A##arr((&A##tmp), npowfs,1, "powfs."#B)==npowfs){\
    for(int i=0; i<npowfs; i++){					\
		parms->powfs[i].B = A##tmp[i];/*doesn't need ## in B*/	\
    }}//else{dbg("Empty array for powfs." #B"\n");}
#define READ_POWFS_MAT(A,B)						\
    if(readcfg_strarr((&strtmp), npowfs, 1,"powfs."#B)==npowfs){\
    for(int i=0; i<npowfs; i++){						\
		parms->powfs[i].B = readstr_##A##mat(0,0,"powfs."#B, strtmp[i]);/*doesn't need ## in B*/ \
		free(strtmp[i]); strtmp[i]=NULL;\
    }}//else{dbg("Empty array for powfs." #B"\n");}

static void convert_theta(real *theta, const char *name, real wvl, real dsa){
	/*Convert pixtheta to radian and do senity check*/
	real val=*theta;
	const char *tmp=NULL;
	if(*theta<=0){//minus means ratio to lambda/dsa
		tmp="-lambda/D";
		*theta=fabs(*theta)*wvl/dsa;
	} else if(*theta<1e-4){
		tmp="radian";
	} else if(*theta>10){
		tmp="mas";
		*theta*=AS2RAD/1000;/*convert form mas to radian. */
	} else{//input is arcsecond.
		tmp="arcsecond";
		*theta*=AS2RAD;/*convert form arcsec to radian. */
	}
	dbg2("Assume %s=%g is in %s. \n", name, val, tmp);
}
/**
   Read wfs geometry. powfs stands for physical optics wfs,
   it is used to represent the types of WFS.
*/
static void readcfg_powfs(parms_t *parms){
	int   npowfs;
	parms->npowfs=npowfs=readcfg_peek_n("powfs.dsa");
	parms->powfs=mycalloc(parms->npowfs,powfs_cfg_t);
	int *inttmp=NULL;
	real *dbltmp=NULL;
	char **strtmp=NULL;
	READ_POWFS(dbl,dsa);
	READ_POWFS(int,nwvl);
	dmat *wvllist=readcfg_dmat(0,0,"powfs.wvl");

	int count=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		int nwvl=parms->powfs[ipowfs].nwvl;
		parms->powfs[ipowfs].wvl=dnew(nwvl,1);
		real wvlmean=0;
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			real wvl=P(wvllist,count+iwvl);
			if(wvl>1e-3){
				wvl=wvl*1e-6;
			}
			wvlmean+=wvl;
			P(parms->powfs[ipowfs].wvl,iwvl)=wvl;
		}
		count+=nwvl;
		parms->powfs[ipowfs].wvlmean=wvlmean/nwvl;
	}
	if(count!=NX(wvllist)){
		error("powfs.wvl has wrong value\n");
	}
	dfree(wvllist);

	READ_POWFS_RELAX(dbl,siglev);
	READ_POWFS_RELAX(dbl,sigrecon);

	READ_POWFS_RELAX(str,saloc);
	READ_POWFS_RELAX(dbl,misregx);
	READ_POWFS_RELAX(dbl,misregy);
	READ_POWFS_RELAX(dbl,misregc);
	READ_POWFS_RELAX(str,amp);
	READ_POWFS_RELAX(str,piinfile);
	READ_POWFS_RELAX(str,sninfile);
	READ_POWFS_RELAX(dbl,saat);
	READ_POWFS_RELAX(dbl,safill2d);
	READ_POWFS_RELAX(dbl,saspherical);
	READ_POWFS_RELAX(dbl,safocuspv);
	READ_POWFS_RELAX(int,neaphy);
	READ_POWFS_RELAX(str,neareconfile);
	READ_POWFS_RELAX(str,neasimfile);
	READ_POWFS_RELAX(dbl,neasim);
	READ_POWFS_RELAX(dbl,neaextra);
	READ_POWFS_RELAX(dbl,neamin);
	READ_POWFS_RELAX(dbl,bkgrnd);
	READ_POWFS_RELAX(dbl,bkgrndc);
	READ_POWFS_RELAX(str,bkgrndfn);
	READ_POWFS_RELAX(str,bkgrndfnc);
	READ_POWFS_RELAX(dbl,pixblur);
	READ_POWFS_RELAX(dbl,radpixtheta);
	READ_POWFS_RELAX(dbl,fieldstop);
	READ_POWFS_RELAX(dbl,astscale);
	READ_POWFS_RELAX(dbl,pixoffx);
	READ_POWFS_RELAX(dbl,pixoffy);
	READ_POWFS_RELAX(int,phyusenea);
	READ_POWFS_RELAX(int,radpix);
	READ_POWFS_RELAX(int,radgx);
	READ_POWFS_RELAX(int,embfac);
	READ_POWFS_RELAX(int,notf);
	READ_POWFS_RELAX(int,psfout);
	READ_POWFS_RELAX(int,pistatout);
	READ_POWFS_RELAX(int,pistatstart);
	READ_POWFS_RELAX(int,pistatstc);
	READ_POWFS_RELAX(int,gtype_sim);
	READ_POWFS_RELAX(int,gtype_recon);
	READ_POWFS_RELAX(int,phytype_recon);
	READ_POWFS_RELAX(int,phytype_sim);
	READ_POWFS_RELAX(int,phytype_sim2);
	READ_POWFS_RELAX(dbl,r0);
	READ_POWFS_RELAX(dbl,L0);
	READ_POWFS_RELAX(int,mtchcpl);
	READ_POWFS_RELAX(int,sigmatch);
	READ_POWFS_RELAX(int,mtchadp);
	READ_POWFS_RELAX(int,mtchfft);
	READ_POWFS_RELAX(dbl,cogthres);
	READ_POWFS_RELAX(dbl,cogoff);
	READ_POWFS_MAT(d,ncpa);
	READ_POWFS_RELAX(int,ncpa_method);
	READ_POWFS_RELAX(int,i0scale);
	READ_POWFS_RELAX(int,i0save);
	READ_POWFS_RELAX(str,i0load);
	READ_POWFS_RELAX(dbl,sigscale);
	READ_POWFS_RELAX(dbl,gradscale);
	READ_POWFS_RELAX(int,moao);
	READ_POWFS_RELAX(int,dither);
	READ_POWFS_RELAX(dbl,dither_amp);
	READ_POWFS_RELAX(int,dither_npoint);
	READ_POWFS_RELAX(int,dither_pllskip);
	READ_POWFS_RELAX(int,dither_pllrat);
	READ_POWFS_RELAX(dbl,dither_gpll);
	READ_POWFS_RELAX(int,dither_ogskip);
	READ_POWFS_RELAX(int,dither_ograt);
	READ_POWFS_RELAX(int,dither_ogsingle);
	READ_POWFS_RELAX(dbl,dither_gog);
	READ_POWFS_RELAX(dbl,dither_gdrift);
	READ_POWFS_RELAX(dbl,dither_glpf);
	//READ_POWFS_RELAX(int, zoomdtrat);
	READ_POWFS_RELAX(int,zoomshare);
	READ_POWFS_RELAX(dbl,zoomgain);
	READ_POWFS_RELAX(dbl,zoomgain_drift);
	READ_POWFS_RELAX(int,zoomset);
	READ_POWFS_RELAX(dbl,apfsm);
	READ_POWFS_RELAX(dbl,epfsm);
	READ_POWFS_RELAX(dbl,alfsm);
	READ_POWFS_RELAX(dbl,zetafsm);
	READ_POWFS_RELAX(dbl,f0fsm);
	READ_POWFS_RELAX(int,idealfsm);
	READ_POWFS_RELAX(int,commonfsm);
	READ_POWFS(dbl,hs);
	READ_POWFS_RELAX(dbl,hc);
	READ_POWFS_RELAX(dbl,nearecon);
	READ_POWFS_RELAX(dbl,rne);
	READ_POWFS_MAT(d,qe);
	READ_POWFS_RELAX(dbl,dx);
	READ_POWFS(dbl,pixtheta);
	READ_POWFS_RELAX(int,trs);
	READ_POWFS_RELAX(int,frs);
	READ_POWFS(int,lo);
	READ_POWFS(int,pixpsa);
	READ_POWFS_RELAX(int,mtchcr);
	READ_POWFS_RELAX(int,mtchstc);
	READ_POWFS_RELAX(int,phystep);
	READ_POWFS_RELAX(int,noisy);
	READ_POWFS_RELAX(int,dtrat);
	READ_POWFS_RELAX(int,skip);
	READ_POWFS(int,type);
	READ_POWFS_RELAX(int,step);
	READ_POWFS(int,nwfs);
	int nllt=0;
	int npywfs=0;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		powfs_cfg_t *powfsi=&parms->powfs[ipowfs];
		if(!isinf(powfsi->hs)){
			nllt++;
		}
		if(powfsi->type==WFS_PY){
			npywfs++;
		}
	}

	if(npywfs){
		READ_POWFS_RELAX(str, pywfs);
	} else{
		readcfg_ignore("powfs.pywfs");
	}
	if(nllt){
		READ_POWFS_RELAX(str, fnllt);
	} else{
		readcfg_ignore("powfs.fnllt");
	}
	dbg("There are %d LGS powfs. \n", nllt);

	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		powfs_cfg_t *powfsi=&parms->powfs[ipowfs];
		real wvlmax=dmax(powfsi->wvl);
		if(powfsi->dsa<=-1){//Order
			powfsi->order=(int)(-powfsi->dsa);
			powfsi->dsa=parms->aper.d/powfsi->order;
		}else if(powfsi->dsa<0){//In unit of aper.d
			//dbg("powfs.dsa=[-1/order] is deprecated. Please use powfs.dsa=[-order] instead.\n");
			powfsi->order=ceil(-1./powfsi->dsa);
			powfsi->dsa*=-parms->aper.d;
		}else{
			if(!powfsi->dsa){
				if(powfsi->lo || !parms->ndm){
					error("powfs[%d].dsa must be set for LO powfs.\n", ipowfs);
				} else{//Follow ground DM.
					if(parms->ndm){
						powfsi->dsa=parms->dm[0].dx;
					} else{
						error("powfs[%d].dsa must be set when there is no DM.\n", ipowfs);
					}
				}
			}
			powfsi->order=ceil(parms->aper.d/powfsi->dsa);
		}
		if(powfsi->fieldstop){
			convert_theta(&powfsi->fieldstop, "fieldstop", powfsi->wvlmean, parms->aper.d);
		}
		if(isfinite(powfsi->hs)){//LGS
			if(parms->sim.htel){
				if(fabs(powfsi->hs-90e3)<2e3){//Only adjust Sodium LGS if it set to 90km+/- 2km
					info("Adjust Sodium LGS range %g by telescope altitude %gm.\n", powfsi->hs, parms->sim.htel);
					powfsi->hs-=parms->sim.htel;
				}else{
					dbg("Telescope altitude is not used.\n");
					parms->sim.htel=0;
				}
			}
			if(!powfsi->fnllt){
				error("powfs%d is at finity range but LLT is not specified.\n", ipowfs);
			}
			char prefix[60]={0};
			snprintf(prefix, 60, "powfs%d_", ipowfs);
#define READ_LLT(T,key) llt->key=readcfg_##T("%sllt."#key, prefix)
#define READ_LLT_ARR(T,key) llt->key=readcfg_##T(0,0,"%sllt."#key, prefix)
			open_config_prefix(powfsi->fnllt, prefix, "powfs.fnllt");
			llt_cfg_t *llt=powfsi->llt=mycalloc(1, llt_cfg_t);
			READ_LLT(dbl, d);
			READ_LLT(dbl, widthp);
			READ_LLT(dbl, focus);
			READ_LLT(dbl, ttrat);
			READ_LLT(dbl, fcfsm);
			READ_LLT(dbl, dhs);
			READ_LLT(str, ttpsd);
			READ_LLT(str, fnrange);
			READ_LLT(str, fnprof);
			READ_LLT(str, fnprep);
			READ_LLT(str, fnsurf);
			READ_LLT(str, fnamp);
			READ_LLT(int, na_smooth);
			READ_LLT(int, na_interp);
			READ_LLT(dbl, na_thres);
			READ_LLT(dbl, na_fit_dh);
			READ_LLT(dbl, na_fit_svdthres);
			READ_LLT(int, na_fit_maxit);
			READ_LLT_ARR(dmat, ox);
			READ_LLT_ARR(dmat, oy);
			READ_LLT_ARR(dmat, misreg);

			READ_LLT(int, ttfr);
			READ_LLT(int, colprep);
			READ_LLT(int, colsim);
			READ_LLT(int, coldtrat);
			READ_LLT(int, nhs);
			llt->nllt=NX(llt->ox);
			if(llt->fcfsm!=0&&llt->nllt>1){
				error("FSM to common LLT FSM offload is only supported for single LLT.\n");
			}
			int mwfs=powfsi->nwfs;
			llt->i=lnew(mwfs, 1);/*default to zero. */
			if(llt->nllt>1){
				/*this is single llt for this powfs. */
				if(llt->nllt!=mwfs)
					error("# of llts should either be 1 or match nwfs for this powfs");
				for(int iwfs=0; iwfs<llt->nllt; iwfs++){
					P(llt->i, iwfs)=iwfs;
				}
			}

		} else{//NGS
			if(powfsi->fnllt){
				warning("powfs%d is NGS but LLT is specified which will be ignored.\n", ipowfs);
				free(powfsi->fnllt);
				powfsi->fnllt=NULL;
			}
			if(powfsi->radpix){
				warning("powfs%d is NGS but radpix is set which will be ignored\n", ipowfs);
				powfsi->radpix=0;
			}
		}
		pywfs_cfg_t *pycfg=NULL;
		if(powfsi->pywfs||powfsi->type==WFS_PY){
			if(powfsi->type==WFS_SH||!powfsi->pywfs){
				error("powfs%d: pywfs must and must only be supplied (%s) for Pyramid WFS.\n", ipowfs, powfsi->pywfs);
			}
			char prefix[60]={0};
			snprintf(prefix, 60, "powfs%d_", ipowfs);
			open_config_prefix(powfsi->pywfs, prefix, "powfs.pywfs");
			pycfg=powfsi->pycfg=mycalloc(1, pywfs_cfg_t);
#define READ_PYWFS(T,key) pycfg->key=readcfg_##T("%spywfs."#key, prefix)
#define READ_PYWFS_MAT(T,key) pycfg->key=readcfg_##T##mat(0,0,"%spywfs."#key, prefix)
			READ_PYWFS(dbl, modulate);pycfg->modulate*=wvlmax/parms->aper.d;
			READ_PYWFS(int, modulpos);
			READ_PYWFS(int, modulring);
			READ_PYWFS(int, nside);
			READ_PYWFS(int, raw);
			READ_PYWFS(dbl, flate);
			READ_PYWFS(dbl, flatv);
			READ_PYWFS(dbl, pupelong);
			READ_PYWFS_MAT(d, psx);
			READ_PYWFS_MAT(d, psy);
			for(int i=0; i<npywfs; i++){
				pycfg->flate*=MAS2RAD;
				pycfg->flatv*=MAS2RAD;
			}
			pycfg->modulpos=pycfg->modulate>0?(pycfg->modulpos/pycfg->nside*pycfg->nside):1;
			pycfg->modulring=pycfg->modulate>0?MAX(1, pycfg->modulring):1;
			pycfg->hs=powfsi->hs;
			pycfg->hc=powfsi->hc;
			pycfg->dx=parms->powfs[ipowfs].dx;
			powfsi->phytype_recon=powfsi->phytype_sim=powfsi->phytype_sim2=2;//like quad cell cog
			powfsi->pixpsa=2;//always 2x2 pixels by definition.
			if(powfsi->phyusenea==-1){
				powfsi->phyusenea=1;
			} else if(powfsi->phyusenea!=1){
				error("PWFS must have phyusenea=1;\n");
			}
			powfsi->pixtheta=0;
			powfsi->ng=pywfs_ng(pycfg); //number of gradients per subaperture
			if(powfsi->sigmatch==-1){
				powfsi->sigmatch=2;//global match
			}
			pycfg->sigmatch=powfsi->sigmatch;
		}else{//SHWFS
			powfsi->ng=2;
			/*Adjust dx if the subaperture does not contain integer, even number of points.*/
			{
				int nx=2*(int)round(0.5*powfsi->dsa/powfsi->dx);
				if(nx<2) nx=2;
				real dx=powfsi->dsa/nx;/*adjust dx. */
				if(fabs(powfsi->dx-dx)>EPS){
					info("powfs%d: Adjusting dx from %g to %g. \n", ipowfs, powfsi->dx, dx);
				}
				powfsi->dx=dx;
			}
			convert_theta(&powfsi->pixtheta, "pixtheta", powfsi->wvlmean, powfsi->dsa);
			if(!powfsi->radpixtheta){
				powfsi->radpixtheta=powfsi->pixtheta;
			} else{
				convert_theta(&powfsi->radpixtheta, "radpixtheta", powfsi->wvlmean, powfsi->dsa);
			}
			if(powfsi->phytype_sim==-1){
				powfsi->phytype_sim=powfsi->phytype_recon;
			}
			if(powfsi->phytype_sim2==-1){
				powfsi->phytype_sim2=powfsi->phytype_sim;
			}
			if(powfsi->phyusenea==-1){
				if(powfsi->phytype_recon==PTYPE_COG){
					powfsi->phyusenea=1;//COG use NEA by default
				} else{
					powfsi->phyusenea=0;
				}
			}
			long pixpsay=powfsi->pixpsa;
			long pixpsax=powfsi->radpix;
			if(!pixpsax) pixpsax=pixpsay;
			if(pixpsax*pixpsay<4){
				powfsi->mtchcr=0;//cannot do constraint.
			}

			if(powfsi->radgx&&!powfsi->radpix){
				powfsi->radgx=0;
			}
			if(powfsi->cogthres<0){
				powfsi->cogthres*=-powfsi->rne;
			}
			if(powfsi->cogoff<0){
				powfsi->cogoff*=-powfsi->rne;
			}
			if(powfsi->sigmatch==-1){
				powfsi->sigmatch=1;
			}
			if(powfsi->phytype_sim==2||powfsi->phytype_sim2==2){//COG
				if((powfsi->cogthres||powfsi->cogoff)&&powfsi->sigmatch!=1){
					error("When cogthres or cogoff is set, only sigmatch==1 is supported but is %d\n", powfsi->sigmatch);
				}
			}
		}


		/*
		if(powfsi->fndither){
			open_config(powfsi->fndither, prefix, 0);
			powfsi->dither=mycalloc(1, dither_cfg_t);
#define READ_DITHER(type,key) powfsi->dither->key=readcfg_##type("%sdither."#key, prefix)
			READ_DITHER(int, mode);
			READ_DITHER(dbl, amp);
			READ_DITHER(int, npoint);
			READ_DITHER(int, pllskip);
			READ_DITHER(int, pllrat);
			READ_DITHER(dbl, gpll);
			READ_DITHER(int, ogskip);
			READ_DITHER(int, ograt);
			READ_DITHER(int, ogsingle);
			READ_DITHER(dbl, gog);
			READ_DITHER(dbl, gdrift);
		}*/
	}

	free(inttmp);
	free(dbltmp);
	free(strtmp);
}
#define READ_WFS(A,B)					\
    if(readcfg_##A##arr((&A##tmp),nwfs,0,"wfs."#B)==nwfs){\
    for(i=0; i<nwfs; i++){				\
		parms->wfs[i].B = A##tmp[i];			\
    }}//else{dbg("Empty array for wfs." #B"\n");}
#define READ_WFS_RELAX(A,B)				\
    if(readcfg_##A##arr((&A##tmp),nwfs,1,"wfs."#B)==nwfs){\
    for(i=0; i<nwfs; i++){				\
		parms->wfs[i].B = A##tmp[i];			\
    }}//else{dbg("Empty array for wfs." #B"\n");}
#define READ_WFS_RELAX_SCALE(A,B,C,D)				\
    if(readcfg_##A##arr((&A##tmp),nwfs,1,"wfs."#B)==nwfs){\
    for(i=0; i<nwfs; i++){				\
		parms->wfs[i].C = D*(A##tmp[i]);			\
    }}//else{dbg("Empty array for wfs." #B"\n");}
#define READ_WFS_DELTA(A,B,BD)				\
    if(readcfg_##A##arr((&A##tmp),nwfs,1,"wfs."#BD)==nwfs){\
    for(i=0; i<nwfs; i++){				\
		parms->wfs[i].B = parms->powfs[parms->wfs[i].powfs].B+A##tmp[i];	\
    }}//else{dbg("Empty array for wfs." #B"\n");}
#define READ_WFS_MAT(A,B)						\
    if(readcfg_strarr((&strtmp), nwfs, 1,"wfs."#B)==nwfs){\
    for(i=0; i<nwfs; i++){						\
		parms->wfs[i].B = readstr_##A##mat(0,0,"wfs."#B,strtmp[i]); \
		free(strtmp[i]); strtmp[i]=NULL;\
    }}//else{dbg("Empty array for wfs." #B"\n");}
/**
   Read in parameters of wfs, excluding signal level.
*/
static void readcfg_wfs(parms_t *parms){
	int i;
	int nwfs=parms->nwfs=readcfg_peek_n("wfs.thetax");
	parms->wfs=mycalloc(parms->nwfs,struct wfs_cfg_t);
	real *dbltmp=NULL;
	int *inttmp=NULL;
	char **strtmp=NULL;
	READ_WFS(dbl,thetax);
	READ_WFS(dbl,thetay);
	for(i=0; i<parms->nwfs; i++){
		parms->wfs[i].thetax*=AS2RAD;
		parms->wfs[i].thetay*=AS2RAD;
	}
	READ_WFS_RELAX(dbl,fitwt);
	READ_WFS_MAT(d, sabad);	
	READ_WFS_RELAX(dbl, misregx);
	READ_WFS_RELAX(dbl, misregy);
	READ_WFS_RELAX(dbl, misregc);//read in as degree.

	/*link wfs with powfs*/
	int wfscount=0;
	int ipowfs=0;
	
	for(int kpowfs=0; kpowfs<parms->npowfs; kpowfs++,ipowfs++){
		if(parms->powfs[kpowfs].nwfs==0){//no stars, remove
			free_powfs_cfg(&parms->powfs[kpowfs]);
			ipowfs--;
			continue;
		} else{
			if(ipowfs<kpowfs){
				memcpy(parms->powfs+ipowfs,parms->powfs+kpowfs,sizeof(powfs_cfg_t));
			}
		}
		int mwfs=parms->powfs[ipowfs].nwfs;
		parms->powfs[ipowfs].wfs=lnew(mwfs,1);
		parms->powfs[ipowfs].wfsind=lnew(parms->nwfs,1);
		int count=0;
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			if(iwfs>=wfscount&&iwfs<wfscount+mwfs){
				parms->wfs[iwfs].powfs=ipowfs;
				P(parms->powfs[ipowfs].wfs,count)=iwfs;
				P(parms->powfs[ipowfs].wfsind,iwfs)=count;
				count++;
			} else{
				P(parms->powfs[ipowfs].wfsind,iwfs)=-1;/*not belong */
			}
		}
		wfscount+=mwfs;
	}
	parms->npowfs=ipowfs;
	if(parms->npowfs==0){
		info("No wfs is specified\n");
		if(!parms->sim.idealtomo&&!parms->sim.evlol){
			error("Cannot proceed\n");
		}
	} else if(parms->nwfs!=wfscount){
		error("parms->nwfs=%d and sum(parms->powfs[*].nwfs)=%d mismatch\n",
			parms->nwfs,wfscount);
	}
	
	READ_WFS_DELTA(dbl,hs,delta_hs);
	READ_WFS_DELTA(dbl,hc,delta_hc);
	free(dbltmp);
	free(inttmp);
	free(strtmp);
	rand_t stat;
	int MISREG_SEQ=1;
	READ_ENV_INT(MISREG_SEQ, 0, INFINITY);
	seed_rand(&stat, MISREG_SEQ);
	
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].astscale==0){
			dbg_once("powfs.astscale=0 is ignored.\n");
		}
		if(parms->powfs[ipowfs].astscale!=0 && parms->powfs[ipowfs].astscale!=1){//scale asterism
			parms->wfs[iwfs].thetax*=parms->powfs[ipowfs].astscale;
			parms->wfs[iwfs].thetay*=parms->powfs[ipowfs].astscale;
		}
		if(parms->powfs[ipowfs].misregx||parms->powfs[ipowfs].misregy||parms->powfs[ipowfs].misregc){
			if(parms->wfs[iwfs].misregx||parms->wfs[iwfs].misregy||parms->wfs[iwfs].misregc){
				warning_once("When both wfs.misreg and powfs.misreg are set, the former is used\n");
			} else{
				if(iwfs==P(parms->powfs[ipowfs].wfs,0)){
					dbg("powfs[%d].misreg=%6.2f %6.2f %7.2f is converted to\n", ipowfs, 
					parms->powfs[ipowfs].misregx, parms->powfs[ipowfs].misregy, parms->powfs[ipowfs].misregc);
				}
				int do_rand=parms->powfs[ipowfs].nwfs>1 && MISREG_SEQ!=0;
				parms->wfs[iwfs].misregx=(do_rand?(2*randu(&stat)-1):1)*parms->powfs[ipowfs].misregx;
				parms->wfs[iwfs].misregy=(do_rand?(2*randu(&stat)-1):1)*parms->powfs[ipowfs].misregy;
				parms->wfs[iwfs].misregc=(do_rand?(2*randu(&stat)-1):1)*parms->powfs[ipowfs].misregc;
			}
		
			dbg("  wfs[%d].misreg=%6.2f %6.2f %7.2f\n",
				iwfs, parms->wfs[iwfs].misregx, parms->wfs[iwfs].misregy, parms->wfs[iwfs].misregc);
		}
		if((parms->wfs[iwfs].misregx||parms->wfs[iwfs].misregy||parms->wfs[iwfs].misregc)&&parms->powfs[ipowfs].type!=WFS_SH){
			warning("wfs.misreg is only used for Shack Hartmann WFS\n");
		}
		//parms->wfs[iwfs].misregx=wrap2range(parms->wfs[iwfs].misregx, -0.5, 0.5);
		//parms->wfs[iwfs].misregy=wrap2range(parms->wfs[iwfs].misregy, -0.5, 0.5);
		parms->wfs[iwfs].misregx*=parms->powfs[ipowfs].dsa;
		parms->wfs[iwfs].misregy*=parms->powfs[ipowfs].dsa;
		parms->wfs[iwfs].misregc*=M_PI/180;//convert degree to radian.
		real rmax=RSS(parms->wfs[iwfs].misregx, parms->wfs[iwfs].misregy);
		if(parms->powfs[ipowfs].misregrmax<rmax){
			parms->powfs[ipowfs].misregrmax=rmax;
		}
		/*
		 * for simulation: raytracing of wfs.misreg and distortion.tel2fs is handled by powfs.loc_tel. amplitude effect is included in powfs.amp. 
		 * for reconstruction: wfsr.misregx,y,c is handled by tomography ray trace and a rotation step. amplitude effect is included in gamp in GP. ploc aligns with valid subapertures. 
		 * ploc (together with saloc, powfs.loc), is defined in the WFS reference frame. 
		 * powfs.loc_tel is defined in the telescope pupil reference frame.
		 * 
		 * notice that misregistration is affected by the rotation.
		 * notice that for ray tracing, misregistration need to be scaled by the cone effect.
		 */
	}
}
/**
 * Compute photon flux for a given magnitude.
*/
static real calc_flux(const dmat *bzero, real wvl, real thruput, real dsa, real dtref, real mag){
	if(NX(bzero)!=2){
		error("bzero dimension is %ldx%ld, nx must be 2.\n", NX(bzero), NY(bzero));
	}
	int iwvl2=0;
	real wvldiff2=INFINITY;
	for(int iwvl=0; iwvl<NY(bzero); iwvl++){
		real wvli=P(bzero, 0, iwvl);
		if(wvli>1e-3){
			wvli*=1e-6;//convert from um to m.
		}
		real wvldiff=fabs(wvli-wvl);
		if(wvldiff<wvldiff2){
			wvldiff2=wvldiff;
			iwvl2=iwvl;
		}
	}
	if(wvldiff2>fabs(wvl)*0.5){
		error("wvl=%g is not found in bzero\n", wvl);
	}

	return pow(10, -0.4*mag)*P(bzero, 1, iwvl2)*thruput*dsa*dsa*dtref;
}
/**
   Read in signal level parameters of wfs.
*/
static void readcfg_siglev(parms_t *parms){
	int wfs_wvl_tot=0;
	int powfs_wvl_tot=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		const int nwvl=parms->powfs[ipowfs].nwvl;
		wfs_wvl_tot+=parms->powfs[ipowfs].nwfs*nwvl;
		powfs_wvl_tot+=nwvl;
	}

	dmat *powfs_wvlwts=readcfg_dmat(powfs_wvl_tot, 1, "powfs.wvlwts");
	dmat *powfs_mag =readcfg_dmat(powfs_wvl_tot, 1, "powfs.mag");
	dmat *powfs_magb=readcfg_dmat(powfs_wvl_tot, 1, "powfs.magbkgrnd");
	dmat *telthruput=readcfg_dmat(powfs_wvl_tot, 1, "powfs.telthruput");
	dmat *atmthruput=readcfg_dmat(powfs_wvl_tot, 1, "powfs.atmthruput");
	//specified for wfs overrides numbers set in powfs.
	dmat *wfs_siglev=readcfg_dmat(0, 0, "wfs.siglev");
	dmat *wfs_wvlwts=readcfg_dmat(0, 0, "wfs.wvlwts");
	dmat *wfs_mag   =readcfg_dmat(0, 0, "wfs.mag");
	dmat *bzero     =readcfg_dmat(0, 0, "sim.bzero");//zero magnitude flux at each wavelength.
	//dmat *wfs_telt=readcfg_dmat(0, 0, "wfs.telthruput");
	//dmat *wfs_atmt=readcfg_dmat(0, 0, "wfs.atmthruput");

	const real cosz=cos(parms->sim.za);
	const real secz=1./cosz;

	int siglev_shared=1;//all wfs in this power share the same siglev and wvlwts.
#define check_dimension(name, arr, size)\
	if(NX(arr)!=0 && NX(arr)!=size){\
		error(name " must be either empty or a vector of %d entries.\n", size);\
	}
	check_dimension("wfs.siglev", wfs_siglev, parms->nwfs);
	check_dimension("wfs.wvlwts", wfs_wvlwts, wfs_wvl_tot);
	check_dimension("wfs.mag",    wfs_mag,    wfs_wvl_tot);
#undef check_dimension
	/*if(dsum(wfs_mag)>0||dsum(powfs_mag)>0){
		info("When both wfs.mag and powfs.mag are specified, wfs.mag takes precedence\n");
	} else{
		info("Computing siglev from %s\n", PN(wfs_mag)?"wfs.mag":"powfs.mag");
	}*/
	int wfs_wvl_count=0;
	int powfs_wvl_count=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		const int nwvl=parms->powfs[ipowfs].nwvl;
		int wvlwts_diff=0;//mark that wvlwts are different
		const int iwfs0=P(parms->powfs[ipowfs].wfs, 0);
		real siglev_sum=0;
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			const int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
			parms->wfs[iwfs].wvlwts=dnew(nwvl,1);
			real siglev_mag=0;
			real bkgrnd_mag=0;
			if(PN(wfs_mag) || PN(powfs_mag) || PN(powfs_magb)){
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					const real mag=(PN(wfs_mag)&&P(wfs_mag, wfs_wvl_count+iwvl))?P(wfs_mag, wfs_wvl_count+iwvl):P(powfs_mag, powfs_wvl_count+iwvl);//do not scale by cosz. It is scaled later in postproc_za.
					const real magb=PN(powfs_magb)?P(powfs_magb, powfs_wvl_count+iwvl):0;
					if(!mag && !magb) continue;//mag==0 is ignored
					const real telt=P(telthruput, powfs_wvl_count+iwvl);
					const real atmt=P(atmthruput, powfs_wvl_count+iwvl);
					const real thruput=telt*pow(atmt, secz); //total throughput
					const real dsa=parms->powfs[ipowfs].dsa;
					const real wvl=P(parms->powfs[ipowfs].wvl, iwvl);
					if(mag){
						P(parms->wfs[iwfs].wvlwts, iwvl)=calc_flux(bzero, wvl, thruput, dsa, parms->sim.dtref, mag);
						siglev_mag+=P(parms->wfs[iwfs].wvlwts, iwvl);
					}
					if(magb){
						real as2=0;//arcsecond squared.
						switch(parms->powfs[ipowfs].type){
						case WFS_PY:
							as2=M_PI*pow(parms->powfs[ipowfs].fieldstop*RAD2AS, 2)/4/parms->powfs[ipowfs].pycfg->nside;
							break;
						case WFS_SH:
							as2=pow(parms->powfs[ipowfs].pixtheta*RAD2AS, 2);
							break;
						default:
							error("Please implement this powfs.type.\n");
						}
						if(as2>1000){
							error("incorrect input: Pixel equivalent area is %g arcsec^2 on sky.\n", as2);
						}
						bkgrnd_mag+=calc_flux(bzero, wvl, thruput, dsa, parms->sim.dtref, magb)*as2;
					}
				}
			}
			if(bkgrnd_mag){
				parms->powfs[ipowfs].bkgrnd=bkgrnd_mag;
			}
			if(siglev_mag){//use siglev and wvlwts from magnitude
				parms->wfs[iwfs].siglev=siglev_mag;
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(parms->wfs[iwfs].wvlwts, iwvl)/=siglev_mag;
				}
				dbg("wfs %d siglev %g is from magnitude\n", iwfs, parms->wfs[iwfs].siglev);
			}else{//use supplied siglev and wvlwts
				if(NX(wfs_siglev)==0||!P(wfs_siglev, iwfs)){
					parms->wfs[iwfs].siglev=parms->powfs[ipowfs].siglev;
					dbg("wfs %d siglev %g is from powfs.siglev\n", iwfs, parms->wfs[iwfs].siglev);
				} else{
					parms->wfs[iwfs].siglev=P(wfs_siglev,iwfs);
					dbg("wfs %d siglev %g is from wfs.siglev\n", iwfs, parms->wfs[iwfs].siglev);
				}
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					real wvlwti=1;
					if(nwvl>1){
						if(NX(wfs_wvlwts)){
							wvlwti=P(wfs_wvlwts, wfs_wvl_count+iwvl);
						}
						if(!wvlwti && NX(powfs_wvlwts)){
							wvlwti=P(powfs_wvlwts, powfs_wvl_count+iwvl);
						}
						if(!wvlwti){
							wvlwti=1./nwvl;//default to 1/nwvl;
						}
					}
					P(parms->wfs[iwfs].wvlwts, iwvl)=wvlwti;
					if(fabs(P(parms->wfs[iwfs0].wvlwts, iwvl)-P(parms->wfs[iwfs].wvlwts, iwvl))>1e-3){
						wvlwts_diff=1;
						siglev_shared=0;
					}
				}
				dnormalize_sum(parms->wfs[iwfs].wvlwts, 1);
			}
			if(fabs(parms->wfs[iwfs].siglev-parms->wfs[iwfs0].siglev)>EPS){
				siglev_shared=0;
			}
			siglev_sum+=parms->wfs[iwfs].siglev;
			wfs_wvl_count+=nwvl;
		}//for jwfs
		if(wvlwts_diff){
			parms->powfs[ipowfs].wvlwts=dnew(nwvl, parms->powfs[ipowfs].nwfs);
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				const int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(parms->powfs[ipowfs].wvlwts, iwvl, jwfs)=P(parms->wfs[iwfs].wvlwts, iwvl);
				}
			}
		}else{
			dcp(&parms->powfs[ipowfs].wvlwts, parms->wfs[iwfs0].wvlwts);
		}
		
		//parms->powfs[ipowfs].wvlmean=ddot(parms->powfs[ipowfs].wvl, parms->wfs[iwfs0].wvlwts)/dsum(parms->wfs[iwfs0].wvlwts);
		if(parms->powfs[ipowfs].lo){
			real nembed=2;
			real dsa=parms->powfs[ipowfs].dsa;
			real wvlpix=parms->powfs[ipowfs].pixtheta*(nembed*dsa);
			if(fabs(wvlpix/parms->powfs[ipowfs].wvlmean-1)<0.5){
				parms->powfs[ipowfs].wvlmean=wvlpix;//change to pixel size equivalent wavelength assuming Nyquist sampling.
				dbg("powfs[%d].wvlmean is set to %g to match detector pixel scale.\n", ipowfs, parms->powfs[ipowfs].wvlmean);
			}else{
				//todo: this weighting does not take into account Strehl ratio and may be too short.
				dbg("powfs[%d].wvlmean is set to %g based on average.\n", ipowfs, parms->powfs[ipowfs].wvlmean);
			}
		}
		
		powfs_wvl_count+=nwvl;
		dbg3("powfs[%d].wvlwts is %ldx%ld\n",ipowfs,NX(parms->powfs[ipowfs].wvlwts),NY(parms->powfs[ipowfs].wvlwts));
		parms->powfs[ipowfs].siglevs=dnew(siglev_shared?1:parms->powfs[ipowfs].nwfs, 1);
		for(int jwfs=0; jwfs<NX(parms->powfs[ipowfs].siglevs); jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
			P(parms->powfs[ipowfs].siglevs,jwfs)=parms->wfs[iwfs].siglev*MAX(1,parms->powfs[ipowfs].dtrat);
		}
		if(!siglev_shared){
			parms->powfs[ipowfs].siglev=siglev_sum/parms->powfs[ipowfs].nwfs;//update to the average of all wfs in this powfs.
		}
		dbg3("powfs[%d].siglevs is %ldx1\n",ipowfs,NX(parms->powfs[ipowfs].siglevs));
	}

	dfree(wfs_siglev);
	dfree(wfs_wvlwts);
	dfree(powfs_wvlwts);

	dfree(wfs_mag);

	dfree(powfs_mag);
	dfree(powfs_magb);

	dfree(telthruput);
	dfree(atmthruput);

	dfree(bzero);
	//dfree(wfs_telt);
	//dfree(wfs_atmt);


}

#define READ_DM(A,B)				\
    readcfg_##A##arr((&A##tmp),ndm,0,"dm."#B);	\
    for(i=0; i<ndm; i++){			\
	parms->dm[i].B = A##tmp[i];		\
    }

#define READ_DM_RELAX(A,B)				\
    readcfg_##A##arr((&A##tmp),ndm,1,"dm."#B);	\
    for(i=0; i<ndm; i++){				\
	parms->dm[i].B = A##tmp[i];			\
    }
#define READ_DM_MAT(A,B)						\
    readcfg_strarr((&strtmp), ndm, 1,"dm."#B);	\
    for(i=0; i<ndm; i++){						\
	parms->dm[i].B = readstr_##A##mat(0,0,"dm."#B,strtmp[i]); \
	free(strtmp[i]); strtmp[i]=NULL;\
    }
/**
   Read in deformable mirror parameters.
*/
static void readcfg_dm(parms_t *parms){
	int ndm,i;
	ndm=parms->ndm=readcfg_peek_n("dm.ht");
	parms->dm=mycalloc(parms->ndm,struct dm_cfg_t);
	int *inttmp=NULL;
	real *dbltmp=NULL;
	char **strtmp=NULL;
	dmat **dmattmp=NULL;
	READ_DM(dbl,ht);
	READ_DM_RELAX(dbl,offset);
	READ_DM_RELAX(dbl,dx);
	READ_DM_RELAX(dbl,ar);
	for(int idm=0; idm<ndm; idm++){
		if(parms->dm[idm].dx<=-1){//this is the order.
			parms->dm[idm].dx=-parms->aper.d/parms->dm[idm].dx;
		} else if(parms->dm[idm].dx<0){
			parms->dm[idm].dx*=-parms->aper.d;
		}
		parms->dm[idm].order=ceil(parms->aper.d/parms->dm[idm].dx);
		parms->dm[idm].dy=parms->dm[idm].dx*parms->dm[idm].ar;
		if(parms->dm[idm].ar<=0){
			error("dm.ar must be positive\n");
		}
	}
	READ_DM_RELAX(dbl,guard);
	{
		char **tmp=0;

		int nstroke=readcfg_strarr(&tmp,0,0,"dm.stroke");
		if(nstroke==1&&!zfexist("%s",tmp[0])){//is number array without separation
			real *ret=0;
			readstr_numarr((void **)&ret,NULL,NULL,ndm,1,M_REAL,"dm.stroke",tmp[0]);
			for(int idm=0; idm<ndm; idm++){
				parms->dm[idm].stroke=dnew(1,1);
				P(parms->dm[idm].stroke,0)=ret[idm];
			}
			free(ret);
			free(tmp[0]); tmp[0]=0;
		} else{
			for(int idm=0; idm<ndm; idm++){
				if(idm==0||nstroke==ndm){
					char *stmp;
					real x=strtod(tmp[idm],&stmp);
					if(stmp==tmp[idm]){//unable to parse number
						parms->dm[idm].stroke=readstr_dmat(0, 0, "dm.stroke", tmp[idm]);
					} else{
						parms->dm[idm].stroke=dnew(1,1);
						P(parms->dm[idm].stroke,0)=x;
					}
					free(tmp[idm]); tmp[idm]=NULL;
				} else if(nstroke==1){
					parms->dm[idm].stroke=dref(parms->dm[0].stroke);
				} else{
					error("dm.stroke is in wrong format\n");
				}
			}
		}
		free(tmp);
	}
	READ_DM_RELAX(dbl,iastroke);
	{
		char **tmp=0;
		int nstroke=readcfg_strarr(&tmp,0,0,"dm.strokescale");
		if(nstroke){
			for(int idm=0; idm<ndm; idm++){
				if(idm==0||nstroke==ndm){
					if(tmp[idm]){
						parms->dm[idm].strokescale=dcellread("%s",tmp[idm]);
						free(tmp[idm]); tmp[idm]=NULL;
					}
				} else if(nstroke==1){
					parms->dm[idm].strokescale=dcellref(parms->dm[0].strokescale);
				} else{
					error("dm.stroke is in wrong format\n");
				}
			}
		}
		free(tmp);
	}
	READ_DM_RELAX(dbl,dratio);
	READ_DM_RELAX(dbl,vmisreg);
	READ_DM_RELAX(dbl,histbin);
	READ_DM_RELAX(int,histn);
	READ_DM_RELAX(int,hist);
	READ_DM_RELAX(dbl,iac);
	READ_DM_RELAX(dbl,hyst);
	READ_DM_RELAX(dbl,hyst_alpha);
	READ_DM_RELAX(dbl,hyst_stroke);
	READ_DM_MAT(d,actfloat);
	READ_DM_MAT(d,actstuck);
	READ_DM_RELAX(dbl, nmod)
	free(strtmp);
	free(inttmp);
	free(dbltmp);
	free(dmattmp);
}
#define READ_MOAO(A,B)					\
    readcfg_##A##arr((&A##tmp),nmoao,0,"moao."#B);	\
    for(i=0; i<nmoao; i++){				\
	parms->moao[i].B = A##tmp[i];			\
    }
#define READ_MOAO_RELAX(A,B)				\
    readcfg_##A##arr((&A##tmp),nmoao,1,"moao."#B);	\
    for(i=0; i<nmoao; i++){				\
	parms->moao[i].B = A##tmp[i];			\
    }
#define READ_MOAO_MAT(A,B)						\
    readcfg_strarr((&strtmp), nmoao, 1,"moao."#B);	\
    for(i=0; i<nmoao; i++){						\
	parms->moao[i].B = readstr_##A##mat(0,0,"moao."#B,strtmp[i]); \
	free(strtmp[i]); strtmp[i]=NULL;\
    }
/**
   Read in MOAO parameters.
*/
static void readcfg_moao(parms_t *parms){
	int nmoao=readcfg_peek_n("moao.dx");
	int i;
	parms->nmoao=nmoao;
	parms->moao=mycalloc(nmoao,moao_cfg_t);
	int *inttmp=NULL;
	real *dbltmp=NULL;
	char **strtmp=NULL;
	READ_MOAO_RELAX(dbl,dx);
	for(int imoao=0; imoao<nmoao; imoao++){
		if(parms->moao[imoao].dx<=-1){
			parms->moao[imoao].dx=parms->aper.d/(-parms->moao[imoao].dx);
		} else if(parms->moao[imoao].dx<0){
			parms->moao[imoao].dx*=-parms->aper.d;
		}
		parms->moao[imoao].order=ceil(parms->aper.d/parms->moao[imoao].dx);
	}
	READ_MOAO_RELAX(dbl,iac);
	READ_MOAO_RELAX(dbl,gdm);
	READ_MOAO_RELAX(dbl,stroke);
	READ_MOAO_RELAX(dbl,ar);
	READ_MOAO_RELAX(int,actslave);
	READ_MOAO_RELAX(int,lrt_ptt);
	READ_MOAO_RELAX(dbl,guard);
	READ_MOAO_MAT(d,actstuck);
	READ_MOAO_MAT(d,actfloat);
	free(inttmp);
	free(dbltmp);
	free(strtmp);
}
/**
   Read in atmosphere parameters.
*/
static void readcfg_atm(parms_t *parms){
	READ_DBL(atm.r0z);
	READ_DBL(atm.dx);
	READ_INT(atm.wdrand);
	READ_INT(atm.method);
	READ_INT(atm.frozenflow);
	READ_INT(atm.ninit);
	READ_INT(atm.share);
	READ_INT(atm.r0evolve);
	READ_INT(atm.dtrat);
	READ_INT(atm.interp);
	READ_DMAT(atm.r0logpsdt);
	READ_DMAT(atm.r0logpsds);
	READ_DMAT_N(atm.size, 2);
	READ_DMAT(atm.ht);
	parms->atm.nps=NX(parms->atm.ht);
	READ_DMAT_N(atm.wt, parms->atm.nps);
	READ_DMAT_N(atm.ws, parms->atm.nps);
	READ_DMAT_NMAX(atm.L0, parms->atm.nps);
	if(!parms->atm.wdrand){
		READ_DMAT_NMAX(atm.wddeg, parms->atm.nps);
	}else{
		readcfg_ignore("atm.wddeg");
	}
	if(parms->atm.dtrat){
		int nps=1;
		if(parms->atm.interp>0 && parms->atm.interp<=2){
			nps=2;//1: do not interpolate, 2: interpolate between two frames.
		}else if(parms->atm.interp==3){
			nps=4;//p-chip interpolation with four frames
		}else{
			error("atm.interp=%d is not supported.\n", parms->atm.interp);
		}
		if(parms->atm.nps!=nps){
			parms->atm.nps=nps;
			dresize(parms->atm.ht, parms->atm.nps, 1);
			dresize(parms->atm.wt, parms->atm.nps, 1);
			dresize(parms->atm.ws, parms->atm.nps, 1);
			dresize(parms->atm.L0, parms->atm.nps, 1);
			if(!parms->atm.wdrand)dresize(parms->atm.wddeg, parms->atm.nps, 1);
			dset(parms->atm.ht, P(parms->atm.ht,0));
		}
	}
}
/**
   Read in atmosphere reconstruction parameters.
*/
static void readcfg_atmr(parms_t *parms){
	READ_DBL(atmr.r0z);
	if(parms->atmr.r0z<=0){
		parms->atmr.r0z=parms->atm.r0z;
	}
	READ_DBL(atmr.L0);
	if(parms->atmr.L0<=0){
		parms->atmr.L0=dsum(parms->atm.L0)/parms->atm.nps;
	}
	READ_DMAT(atmr.ht);
	READ_DMAT(atmr.wt);
	if(NX(parms->atmr.ht)==0){
		dcp(&parms->atmr.ht,parms->atm.ht);
	}
	if(NX(parms->atmr.wt)==0){
		dcp(&parms->atmr.wt,parms->atm.wt);
	} else{
		if(NX(parms->atmr.wt)!=NX(parms->atmr.ht)){
			error("atmr.wt length must match atmr.ht\n");
		}
	}
	parms->atmr.nps=NX(parms->atmr.ht);
	parms->atmr.os=readcfg_lmat(parms->atmr.nps,1,"atmr.os");
	READ_DBL(atmr.dx);
}

/**
   Read in aperture definitions.
*/
static void readcfg_aper(parms_t *parms){
	real *dtmp;
	/*aper.d may contain one for [d] or two numbers for [d din] */
	int nd=readcfg_dblarr(&dtmp,0,0,"aper.d");
	switch(nd){
	case 2:
		parms->aper.din=dtmp[1];
		//fallthrough
	case 1:
		parms->aper.d=dtmp[0];
		break;
	default:
		error("aper.d contains %d elements. But only 1 or 2 elements are supported.\n",nd);
	}
	free(dtmp);

	if(parms->aper.d<=parms->aper.din){
		error("Inner dimeter(%g) should be less than Outer Diameter(%g).\n",parms->aper.din,parms->aper.d);
	} else if(parms->aper.d<0||parms->aper.din<0){
		error("Inner (%g) and outer (%g) diameters should be positive.\n",parms->aper.din,parms->aper.d);
	}
	READ_DBL_SCALE(aper.rot, aper.rotdeg, M_PI/180.);
	parms->aper.misreg=readcfg_dmat(2, 0, "aper.misreg");
	READ_STR(aper.fnamp);
	READ_STR(aper.pupmask);
}
/**<
 * check and scale thetax and thetay to match fov in diameter. If thetax and
 * thetay have only a single on axis point, will make a grid of points
 * sufficient pacing. The unit of thetax, thetay, fov are in arcsec.
 * */
static void scale_fov(dmat *thetax,dmat *thetay,dmat *wt,real fov){
	if(PN(thetax)!=PN(thetay)||PN(thetax)!=PN(wt)){
		error("thetax, thetay, and wt mismatch in dimensions:%ld, %ld, %ld\n", PN(thetax), PN(thetay), PN(wt));
	}
	if(fov==0 || fov==1){//when fov is 0 or 1. will not do anything. use thetax and thetay as it.
		return;
	}
	real maxxy=0;
	for(int i=0; i<PN(thetax); i++){
		if(fabs(P(thetax, i))>maxxy) maxxy=fabs(P(thetax, i));
		if(fabs(P(thetay, i))>maxxy) maxxy=fabs(P(thetay, i));
	}
	if(maxxy==0){//need to make a field
		int nx=round(fov/30.)*2+1; //need sufficient sampling of the focal plane
		if(nx<3) nx=3;
		int np=nx*nx;
		warning("maxxy=0, will make a square field with %dx%d points.\n", nx, nx);

		dresize(thetax, nx, nx);
		dresize(thetay, nx, nx);
		dresize(wt, nx, nx);
		real cx=(nx-1)*0.5;
		real rfov=fov/(2*cx);
		for(int iy=0; iy<nx; iy++){
			//sympson weighting is 1 4 2 4 2 ... 2 4 1 for odd number of points
#define SIMPSON_1D(ix,nx) (((ix)==0||(ix)+1==(nx))?1:(((ix)%2==0)?2:4))
			for(int ix=0; ix<nx; ix++){
				P(thetax, ix, iy)=(ix-cx)*rfov;
				P(thetay, ix, iy)=(iy-cx)*rfov;
				P(wt, ix, iy)=SIMPSON_1D(ix,nx)*SIMPSON_1D(iy,nx);
			}
		}
		//the weighting is very different from simpson weighting
		reshape(thetax, np, 1);
		reshape(thetay, np, 1);
		reshape(wt, np, 1);
	}else{
		if(fabs(maxxy-0.5)>EPS){
			warning("maxxy=%g is not 0.5, adjust the scaling properly.\n", maxxy);
			fov=fov/(2*maxxy);
		}
		dscale(thetax,fov);
		dscale(thetay,fov);
	}
}
const char *wvl2name(real wvl){
	if(wvl<0.01){
		wvl*=1e6; //m -> um
	}else if(wvl>100){
		wvl*=1e-3; //nm->um
	}
	real wc[]={0.365,0.445,0.464,0.551,0.658,0.8,0.9,1.0,1.25,1.65,2.15,3.45};
	const char *n[]={"U","B","G","V","R","I","Z","Y","J","H","K","L"};
	size_t j=0;
	real dist=1;
	for(size_t i=0; i<sizeof(wc)/sizeof(real);i++){
		if(fabs(wc[i]-wvl)<dist){
			dist=fabs(wc[i]-wvl);
			j=i;
		}
	}
	return n[j];
}
/**
   Read in performance evaluation science point parameters.
*/
static void readcfg_evl(parms_t *parms){
	READ_DMAT(evl.thetax);
	READ_DMAT(evl.thetay);
	READ_DMAT_NMAX(evl.wt,NX(parms->evl.thetax));
	real evl_fov=readcfg_dbl("evl.fov");
	scale_fov(parms->evl.thetax,parms->evl.thetay,parms->evl.wt,evl_fov);
	parms->evl.nevl=NX(parms->evl.thetax);//maybe changed by scale_fov
	dnormalize_sumabs(parms->evl.wt, 1);
	READ_DMAT_NMAX(evl.hs,parms->evl.nevl);
	READ_LMAT_NMAX(evl.psf,parms->evl.nevl);
	READ_LMAT_NMAX(evl.psfr,parms->evl.nevl);
	READ_LMAT_NMAX(evl.pttr,parms->evl.nevl);
	READ_DMAT(evl.wvl);
	parms->evl.nwvl=NX(parms->evl.wvl);
	parms->evl.wvlname=mycalloc(parms->evl.nwvl, const char*);
	for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		if(P(parms->evl.wvl,iwvl)>0.1){
			//info("Assume evl.wvl[%d]=%g is supplied in micron.\n", iwvl, P(parms->evl.wvl,iwvl));
			P(parms->evl.wvl,iwvl)*=1e-6;
		}
		parms->evl.wvlname[iwvl]=wvl2name(P(parms->evl.wvl, iwvl));
	}
	parms->evl.psfgridsize=readcfg_lmat(parms->evl.nwvl,1,"evl.psfgridsize");
	parms->evl.psfsize=readcfg_lmat(parms->evl.nwvl,1,"evl.psfsize");
	int ievl;
	real ramin=INFINITY;
	for(ievl=0; ievl<parms->evl.nevl; ievl++){
	/*First Convert theta to radian from arcsec. */
		P(parms->evl.thetax,ievl)*=AS2RAD;
		P(parms->evl.thetay,ievl)*=AS2RAD;
		real ra2=pow(P(parms->evl.thetax,ievl),2)+pow(P(parms->evl.thetay,ievl),2);
		if(ra2<ramin){
			parms->evl.indoa=ievl;
			ramin=ra2;
		}
	}
	READ_DBL(evl.dx); if(parms->evl.dx<=0) parms->evl.dx=parms->atm.dx;
	READ_INT(evl.rmax);
	READ_INT(evl.psfol);
	READ_INT(evl.psfisim);

	READ_INT(evl.psfmean);
	READ_INT(evl.psfhist);
	READ_INT(evl.cov);/*Science OPD covariance. */
	READ_INT(evl.opdmean);/*Science OPD time average.*/
	if(parms->evl.cov){
		parms->evl.opdmean=parms->evl.cov;
	}
	READ_INT(evl.tomo);
	READ_INT(evl.moao);
	READ_INT(evl.split);
	/*it is never good to parallelize the evl ray tracing because it is already so fast */
	parms->evl.nmod=(parms->evl.rmax+1)*(parms->evl.rmax+2)/2;
}
/**
   Read in turbulence tomography parameters.
*/
static void readcfg_tomo(parms_t *parms){
	READ_INT(tomo.pos);
	READ_INT(tomo.cone);
	READ_INT(tomo.square);
	READ_INT(tomo.cxxalg);
	READ_DBL(tomo.cxxscale);
	READ_INT(tomo.guard);
	READ_DBL(tomo.tikcr);
	READ_DBL(tomo.svdthres);
	READ_INT(tomo.piston_cr);
	READ_INT(tomo.ahst_wt);
	READ_INT(tomo.ahst_idealngs);
	READ_INT(tomo.ahst_focus);
	READ_INT(tomo.alg);
	READ_INT(tomo.bgs);
	READ_INT(tomo.precond);
	READ_INT(tomo.maxit);
	READ_INT(tomo.assemble);
	READ_INT(tomo.predict);
	READ_DBL(tomo.minwt);
	READ_DBL(tomo.iac);
	READ_INT(tomo.ninit);
	READ_INT(tomo.nxbase);
	READ_INT(tomo.splitlrt);
}

/**
   Read in DM fit parameters. MOAO is specified elsewhere in readcfg_moao() */
static void readcfg_fit(parms_t *parms){
	READ_DMAT(fit.thetax);
	READ_DMAT(fit.thetay);
	READ_DMAT(fit.wt);
	real fit_fov=readcfg_dbl("fit.fov");
	scale_fov(parms->fit.thetax,parms->fit.thetay,parms->fit.wt,fit_fov);
	parms->fit.nfit=NX(parms->fit.thetax);//maybe changed by scale_fov
	parms->fit.hs=readcfg_dmat(parms->fit.nfit,1,"fit.hs");
	real ramin=INFINITY;
	for(int ifit=0; ifit<parms->fit.nfit; ifit++){
		P(parms->fit.thetax,ifit)*=AS2RAD;
		P(parms->fit.thetay,ifit)*=AS2RAD;
		real ra2=pow(P(parms->fit.thetax,ifit),2)+pow(P(parms->fit.thetay,ifit),2);
		if(ra2<ramin){
			parms->fit.indoa=ifit;
			ramin=ra2;
		}
	}
	//expand fitting directions
	{
		int nwfsfit=0;
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			if(parms->wfs[iwfs].fitwt>EPS){
				nwfsfit++;
			}
		}
		if(nwfsfit){
			int nfit2=parms->fit.nfit+nwfsfit;
			dresize(parms->fit.thetax, nfit2, 1);
			dresize(parms->fit.thetay, nfit2, 1);
			dresize(parms->fit.wt, nfit2, 1);
			dresize(parms->fit.hs, nfit2, 1);
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				if(parms->wfs[iwfs].fitwt<EPS){
					continue;
				}
				info("Including wfs %d with wt %g in fitting\n", iwfs, parms->wfs[iwfs].fitwt);
				P(parms->fit.thetax, parms->fit.nfit)=parms->wfs[iwfs].thetax;
				P(parms->fit.thetay, parms->fit.nfit)=parms->wfs[iwfs].thetay;
				P(parms->fit.hs, parms->fit.nfit)=parms->wfs[iwfs].hs;
				P(parms->fit.wt, parms->fit.nfit)=parms->wfs[iwfs].fitwt;
				parms->fit.nfit++;
			}
			if(nfit2!=parms->fit.nfit){
				error("unexpected: nfit2=%d, parms->fit.nfit=%d\n", nfit2, parms->fit.nfit);
			}
		}
		if(parms->dbg.dmfullfov>1){//Fit over the entire fov. Hurts performance.
			if(!parms->sim.fov){
				error("sim.fov is not specified\n");
			}
			int ndir=8;
			int nfit2=parms->fit.nfit+ndir;
			real rfov=parms->sim.fov*0.5;
			real fitwt=0.001;
			dresize(parms->fit.thetax, nfit2, 1);
			dresize(parms->fit.thetay, nfit2, 1);
			dresize(parms->fit.wt, nfit2, 1);
			dresize(parms->fit.hs, nfit2, 1);
			for(int idir=0; idir<ndir; idir++){
				real theta=(2.*M_PI*idir)/ndir;
				P(parms->fit.thetax, parms->fit.nfit)=rfov*cos(theta);
				P(parms->fit.thetay, parms->fit.nfit)=rfov*sin(theta);
				P(parms->fit.hs, parms->fit.nfit)=P(parms->fit.hs,0);
				P(parms->fit.wt, parms->fit.nfit)=fitwt;
				info("Including virtual direction (%g, %g) with wt %g in fitting\n",
					P(parms->fit.thetax, parms->fit.nfit),P(parms->fit.thetay, parms->fit.nfit),
					P(parms->fit.wt, parms->fit.nfit));
				parms->fit.nfit++;
			}
			if(nfit2!=parms->fit.nfit){
				error("unexpected: nfit2=%d, parms->fit.nfit=%d\n", nfit2, parms->fit.nfit);
			}
		}
	}
	dnormalize_sumabs(parms->fit.wt, 1);
	READ_DBL(fit.tikcr);
	READ_DBL(fit.svdthres);
	READ_INT(fit.actslave);
	READ_DBL(fit.actthres);
	READ_DBL(fit.actthres2);
	READ_INT(fit.actextrap);
	READ_INT(fit.lrt_piston);
	READ_INT(fit.lrt_tt);
	READ_INT(fit.alg);
	READ_INT(fit.bgs);
	READ_INT(fit.precond);
	READ_INT(fit.maxit);
	READ_INT(fit.guard);
	READ_INT(fit.square);
	READ_INT(fit.assemble);
	READ_INT(fit.pos);
	READ_INT(fit.cachedm);
	READ_INT(fit.cachex);
}
/**
   Read LSR parameters.
*/
static void readcfg_lsr(parms_t *parms){
	READ_DBL(lsr.tikcr);
	READ_DBL(lsr.svdthres);
	READ_STR(lsr.fnreg);
	READ_INT(lsr.alg);
	READ_INT(lsr.actslave);
	READ_INT(lsr.bgs);
	READ_INT(lsr.maxit);
	READ_DBL(lsr.actthres);
	READ_DBL(lsr.actthres2);
	READ_INT(lsr.actextrap);
	READ_INT(lsr.splitlrt);
}
/**
   Read general reconstruction parameters
*/
static void readcfg_recon(parms_t *parms){
	READ_INT(recon.alg);
	READ_INT(recon.glao);
	READ_INT(recon.split);
	READ_INT(recon.mvm);
	READ_INT(recon.modal);
	READ_INT(recon.psol);
	READ_DBL(recon.poke);
	parms->nwfsr=parms->recon.glao?parms->npowfs:parms->nwfs;
	readcfg_strarr(&parms->recon.distortion_tel2wfs,parms->nwfsr,1,"recon.distortion_tel2wfs");
	readcfg_strarr(&parms->recon.distortion_dm2wfs,parms->ndm*parms->nwfsr,1,"recon.distortion_dm2wfs");
	readcfg_strarr(&parms->recon.distortion_dm2sci,parms->ndm*parms->fit.nfit,1,"recon.distortion_dm2sci");
	READ_INT(recon.psd);
	READ_INT(recon.psddtrat_hi);
	READ_INT(recon.psddtrat_lo);
	READ_DBL(recon.psdservo_gain);
	READ_INT(recon.psdnseg);
	READ_INT(recon.twfs_rmin);
	READ_INT(recon.twfs_rmax);
	READ_INT(recon.twfs_radonly);
	READ_INT(recon.petal);
	READ_INT(recon.petaldtrat);
	READ_INT(recon.petalstep);
	READ_INT(recon.petalnpsf);
	READ_INT(recon.petaltt);
}
/**
   Read in simulation parameters
*/
static void readcfg_sim(parms_t *parms){
	READ_DBL(sim.fcfocus);
	READ_DBL(sim.fcttm);
	READ_INT(sim.mffocus);
	READ_DBL(sim.epfocus2tel);
	READ_INT(sim.focus2tel);
	READ_DMAT(sim.aphi);
	READ_DMAT(sim.ephi);
	READ_DMAT(sim.aplo);
	READ_DMAT(sim.eplo);
	READ_DBL(sim.f0dm);
	READ_DBL(sim.zetadm);
	//READ_DBL(sim.f0tt);
	//READ_DBL(sim.zetatt);
	READ_DBL(sim.aptwfs);
	READ_DBL(sim.eptwfs);
	READ_DBL(sim.eptsph);
	READ_DBL(sim.alhi);
	READ_DBL(sim.allo);
	/*We append a 0 so that we keep a time history of the integrator. */
	if(NX(parms->sim.aphi)==1){
		dresize(parms->sim.aphi,2,1);
	}
	if(NX(parms->sim.aphi)>2||NY(parms->sim.aphi)!=1){
		error("sim.aphi is invalid. Must have 1 or 2 entries.\n");
	}
	if(NX(parms->sim.aplo)==1){
		dresize(parms->sim.aplo,2,1);
	}
	if(NX(parms->sim.aplo)>2||NY(parms->sim.aplo)!=1){
		error("sim.apli is invalid. Must have 1 or 2 entries.\n");
	}
	parms->sim.seeds=readcfg_lmat(0,0,"sim.seeds");
	parms->sim.nseed=NX(parms->sim.seeds);
	READ_DBL(sim.dt);
	READ_DBL(sim.dtref);
	READ_INT(sim.dtrat_skip);
	READ_INT(sim.start);
	READ_INT(sim.end);
	READ_INT(sim.pause);
	READ_DMAT(sim.wspsd);
	READ_INT(sim.wsseq);
	READ_INT(sim.cachedm);
	READ_INT(sim.fuseint);
	READ_INT(sim.closeloop);
	READ_INT(sim.skysim);
	READ_DBL_SCALE(sim.fov,sim.fov,AS2RAD);
	READ_DBL_SCALE(sim.za,sim.zadeg,M_PI/180);
	READ_INT(sim.htel);
	READ_INT(sim.evlol);
	READ_INT(sim.noatm);
	READ_INT(sim.idealtomo);
	READ_INT(sim.psfr);
	READ_INT(sim.ecnn);
	READ_INT(sim.wfsalias);
	READ_INT(sim.idealwfs);
	READ_INT(sim.idealevl);
	READ_STR(sim.mvmhost);
	READ_INT(sim.mvmport);
	READ_INT(sim.mvmsize);
	READ_INT(sim.mvmngpu);
	READ_STR(sim.dmadd);
}
static void readcfg_ncpa(parms_t *parms){
	READ_INT(ncpa.calib);
	READ_INT(ncpa.ttr);
	READ_DMAT(ncpa.thetax);
	parms->ncpa.ndir=NX(parms->ncpa.thetax);
	READ_DMAT_NMAX(ncpa.thetay, parms->ncpa.ndir);
	READ_DMAT_NMAX(ncpa.wt, parms->ncpa.ndir);
	READ_DMAT_NMAX(ncpa.hs, parms->ncpa.ndir);
	dscale(parms->ncpa.thetax,1.*AS2RAD);
	dscale(parms->ncpa.thetay,1.*AS2RAD);
	if(parms->ncpa.wt){
		dnormalize_sumabs(parms->ncpa.wt, 1);
	}
	READ_INT(ncpa.preload);
	READ_INT(ncpa.rmsci);
	parms->ncpa.nsurf=readcfg_strarr(&parms->ncpa.surf, 0, 0, "ncpa.surf");
	parms->ncpa.ntsurf=readcfg_strarr(&parms->ncpa.tsurf, 0, 0, "ncpa.tsurf");

}
/**
   Read in parameters for Cn2 estimation.
*/
static void readcfg_cn2(parms_t *parms){
/*for Cn2 Estimation. */
	READ_DMAT(cn2.pair);
	READ_INT(cn2.step);
	READ_INT(cn2.reset);
	READ_INT(cn2.tomo);
	READ_INT(cn2.verbose);
	READ_INT(cn2.keepht);
	READ_INT(cn2.nhtomo);
	READ_INT(cn2.moveht);
	READ_INT(cn2.psol);
	READ_DBL(cn2.hmax);
	READ_DBL(cn2.saat);
	if(parms->cn2.pair&&(!NX(parms->cn2.pair)||parms->sim.noatm==1||parms->atm.dtrat)){
		dfree(parms->cn2.pair);
	}
	if(!parms->cn2.pair){/*we are not doing cn2 estimation. */
		parms->cn2.tomo=0;
	}
}

/**
   Specify which variables to plot
*/
static void readcfg_plot(parms_t *parms){
	READ_INT(plot.setup);
	READ_INT(plot.atm);
	READ_INT(plot.run);
	READ_INT(plot.opdx);
	READ_INT(plot.psf);
	READ_INT(plot.all);
	READ_INT(plot.grad2opd);
	READ_DBL(plot.opdmax);
	READ_DBL(plot.gmax);
	READ_DBL(plot.psfmin);
	if(parms->plot.all){
		parms->plot.setup=parms->plot.all;
		if(!parms->plot.run) parms->plot.run=parms->plot.all;
		if(!parms->plot.psf) parms->plot.psf=parms->plot.all;
	}
	/*if(parms->plot.setup||parms->plot.atm||parms->plot.run||parms->plot.opdx||parms->plot.all||parms->plot.psf){
		draw_helper();
	}*/
}

/**
   Read in debugging parameters
*/
static void readcfg_dbg(parms_t *parms){
	READ_DMAT(dbg.atm); if(dsumabs(parms->dbg.atm)==0){ dfree(parms->dbg.atm); parms->dbg.atm=NULL; }
	READ_INT(dbg.mvstlimit);
	READ_INT(dbg.annular_W);
	READ_LMAT(dbg.tomo_maxit);
	READ_INT(dbg.tomo_hxw);
	READ_INT(dbg.ecovxx);
	READ_INT(dbg.useopdr);
	READ_INT(dbg.cmpgpu);
	READ_INT(dbg.wfslinearity);
	READ_INT(dbg.nocgwarm);
	if(readcfg_peek("dbg.test")){
		READ_INT(dbg.test);
	}
	READ_INT(dbg.dmfullfov);
	READ_INT(dbg.tomo);
	READ_INT(dbg.fit);

	READ_INT(dbg.gp_noamp);
	READ_DBL(dbg.gradoff_scale);
	READ_INT(dbg.gradoff_reset);
	READ_DCELL(dbg.dmoff);
	READ_DCELL(dbg.gradoff);
	READ_INT(dbg.wfs_iac);
	READ_INT(dbg.fullatm);
	READ_INT(dbg.lo_blend);
	READ_DBL(dbg.eploscale);
	READ_INT(dbg.ahst_keepfocus);
	READ_INT(dbg.recon_stuck);
}
/**
   Read in GPU options
*/
static void readcfg_gpu(parms_t *parms){
	READ_INT(gpu.wfs);
	READ_INT(gpu.evl);
	READ_INT(gpu.tomo);
	READ_INT(gpu.fit);
	READ_INT(gpu.lsr);
	READ_INT(gpu.psf);
	READ_INT(gpu.moao);
}
/**
   Specify which variables to save
*/
static void readcfg_save(parms_t *parms){
	READ_INT(save.extra);
	READ_INT(save.all);
	READ_INT(save.setup);
	READ_INT(save.recon);
	READ_INT(save.mvst);
	READ_INT(save.ncpa);
	READ_INT(save.fdpcg);
	READ_INT(save.atm);/*Save atmosphere */
	READ_INT(save.run);
	READ_INT(save.opdr);/*reconstructed OPD on XLOC */
	READ_INT(save.opdx);/*ATM propagated to XLOC */
	READ_INT(save.evlopd);/*Science OPD */
	READ_INT(save.dm);/*save DM commands */
	READ_INT(save.dither);
	READ_INT(save.gradoff);//save gradient reference vector
	parms->save.ints=readcfg_lmat(parms->nwfs,1,"save.ints");
	parms->save.wfsopd=readcfg_lmat(parms->nwfs,1,"save.wfsopd");
	parms->save.grad=readcfg_lmat(parms->nwfs,1,"save.grad");
	parms->save.gradnf=readcfg_lmat(parms->nwfs,1,"save.gradnf");
	parms->save.gradpsol=readcfg_lmat(parms->nwfs,1,"save.gradpsol");
	parms->save.gradgeom=readcfg_lmat(parms->nwfs,1,"save.gradgeom");
	if(disable_save){
		parms->save.extra=0;
		if(parms->save.all||parms->save.setup||parms->save.setup){
			error("Save is enabled. Please specify output directory.\n");
		}
	}
	
	if(parms->save.all){
		if(!parms->save.setup) parms->save.setup=parms->save.all;
		if(!parms->save.run) parms->save.run=parms->save.all;
	}
	if(parms->save.setup>1){
		if(!parms->save.recon) parms->save.recon=parms->save.setup;
		if(!parms->save.mvst) parms->save.mvst=parms->save.setup;
		if(!parms->save.ncpa) parms->save.ncpa=parms->save.setup;
		if(!parms->save.fdpcg) parms->save.fdpcg=parms->save.setup;
	}
	if(parms->save.recon||parms->save.mvst||parms->save.ncpa||parms->save.fdpcg){
		if(!parms->save.setup) parms->save.setup=1;
	}
	if(parms->save.run==1){
		info("Saving RTC telemetry.\n");
	}else if(parms->save.run>1){
		info("Saving RTC telemetry and OPDs.\n");
	}
	if(parms->save.run){
		if(!parms->save.gradoff) parms->save.gradoff=parms->save.run;
		if(!parms->save.extra) parms->save.extra=parms->save.run;
		if(!parms->save.dither) parms->save.dither=parms->save.run;
		if(!parms->save.dm) parms->save.dm=parms->save.run;
		lset(parms->save.ints, parms->save.run);
		lset(parms->save.grad, parms->save.run);
		lset(parms->save.gradnf, parms->save.run);
		lset(parms->save.gradpsol, parms->save.run);
		lset(parms->save.gradgeom, parms->save.run);
	}
	if(parms->save.run>1){
		/*The following are run time information that are only enabled with
		save.run>1 or save.all>1 because they take a lot of disk space and slows
		down the simulation dramatically.*/
		if(!parms->save.opdr) parms->save.opdr=parms->save.run;
		if(!parms->save.evlopd) parms->save.evlopd=parms->save.run;
		if(!parms->recon.glao) parms->save.opdx=parms->save.run;
		lset(parms->save.wfsopd, parms->save.run);
	}

	READ_LMAT(save.gcov);
	parms->save.ngcov=NX(parms->save.gcov)/2;
	READ_INT(save.gcovp);
	READ_INT(save.ecov);
	READ_INT(save.mvmi);
	READ_INT(save.mvmf);
	READ_INT(save.mvm);
}
static void readcfg_distortion(parms_t *parms){
	readcfg_strarr(&parms->distortion.tel2wfs,parms->nwfs,1,"distortion.tel2wfs");
	readcfg_strarr(&parms->distortion.dm2wfs,parms->ndm*parms->nwfs,1,"distortion.dm2wfs");
	readcfg_strarr(&parms->distortion.dm2sci,parms->ndm*parms->evl.nevl,1,"distortion.dm2sci");
}
/**
   Specify which variables to load from saved files (Usually from LAOS
   simulations) */
static void readcfg_load(parms_t *parms){

	READ_STR(load.atm);
	READ_STR(load.locs);
	READ_STR(load.aloc);
	READ_STR(load.xloc);
	READ_STR(load.ploc);
	READ_STR(load.floc);
	READ_STR(load.cxx);
	//READ_STR(load.HXF);
	READ_STR(load.HXW);
	//READ_STR(load.HA);
	READ_STR(load.GP);
	READ_STR(load.GA);
	READ_STR(load.mvm);
	READ_STR(load.mvmi);
	READ_STR(load.mvmf);
	READ_INT(load.mvst);
	READ_INT(load.GS0);
	READ_INT(load.tomo);
	READ_INT(load.fit);
	READ_INT(load.W);
	READ_STR(load.ncpa);
	READ_STR(load.saneai);
}
/**
   Scaling necessary values for non-zero zenith angle (za).
   -# r0
   -# turbulence height
   -# WFS height
   -# WFS siglev
   -# reconstruction height.

   DM conjugation range is not changed!
   2012-04-07: Relocated to beginning
*/
static void setup_parms_postproc_za(parms_t *parms){
	real cosz=cos(parms->sim.za);
	real secz=1./cosz;
	
	if(fabs(parms->sim.za)>1.e-14){
		info("Zenith angle is %g degree.\n",parms->sim.za*180./M_PI);
		dbg("    Scaling turbulence height and LGS hs by sec(za).\n"
			"    Scaling r0 by cos(za)^(3/5).\n"
			"    Scaling wind speed and LGS signal level by cos(za).\n");

	}
	/*
	  The input r0z is the r0 at zenith. Scale it if off zenith
	*/
	parms->atm.r0=parms->atm.r0z*pow(cosz,3./5.);
	parms->atmr.r0=parms->atmr.r0z*pow(cosz,3./5.);

	//Adjust turbulence only by zenith angle.
	for(int ips=0; ips<parms->atm.nps; ips++){
		P(parms->atm.ht,ips)*=secz;/*scale atmospheric height */
		P(parms->atm.ws,ips)*=cosz;/*Wind velocity reduced due to line of sight*/
	}

	for(int ips=0; ips<parms->atmr.nps; ips++){
		P(parms->atmr.ht,ips)*=secz;/*scale reconstructed atmospheric height. */
	}
	parms->cn2.hmax*=secz;

	//Adjust LGS height by telescope altitude and sec(za)
	//Adjust LGS signal level by cos(za)
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(!isinf(parms->powfs[ipowfs].hs)){
			
			parms->powfs[ipowfs].hs*=secz;/*scale GS height. */
			parms->powfs[ipowfs].siglev*=cosz;
			dscale(parms->powfs[ipowfs].siglevs,cosz);
			for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs,indwfs);		
				parms->wfs[iwfs].hs*=secz;
				parms->wfs[iwfs].siglev*=cosz;
			}
		}
	}
	dscale(parms->evl.hs,secz);
	dscale(parms->ncpa.hs,secz);
	dscale(parms->fit.hs,secz);
}
/**
   Process simulation parameters to find incompatibility.
*/
static void setup_parms_postproc_sim(parms_t *parms){
	if(parms->sim.skysim){
		if(disable_save){
			error("sim.skysim requires saving. Please specify output folder\n");
		}
		/*if(parms->recon.alg!=RECON_MVR){
		  error("skysim need MVR");
		  }*/
		if(!parms->tomo.ahst_idealngs){
			parms->tomo.ahst_idealngs=1;
		}
		if(parms->ndm>0&&parms->recon.split!=1){
			info("Can only do skysim in split tomography mode 1. Changed\n");
			parms->recon.split=1;
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].lo){
				parms->powfs[ipowfs].psfout=1;
			}
		}
		parms->save.extra=1;
	}
	if(parms->tomo.ahst_idealngs){
		if(parms->sim.fuseint){
			info("Disabling sim.fuseint\n");
			parms->sim.fuseint=0;
		}
	}
	if(NX(parms->dbg.tomo_maxit)){
		warning("dbg.tomo_maxit is set. Will run in open loop mode\n to repeat the simulations"
			" with different values of tomo.maxit.\n");
		parms->sim.closeloop=0;
		parms->atm.frozenflow=1;
		for(int ips=0; ips<parms->atm.nps; ips++){
			P(parms->atm.ws,ips)=0;/*set windspeed to zero. */
		}
		parms->sim.end=NX(parms->dbg.tomo_maxit);
	}
	if(parms->sim.idealtomo){
		if(parms->recon.glao){
			error("idealtomo and recon.glao conflicts\n");
		}
		if(parms->recon.alg!=RECON_MVR){
			error("idealtomo only works in minimum variance reconstruction (recon.alg=0) mode.\n");
			parms->recon.alg=RECON_MVR;
		}
		if(parms->recon.split){
			dbg2("idealtomo only works in integrated tomo mode. changed.\n");
			parms->recon.split=0;
		}
		if(parms->recon.mvm){
			dbg2("idealtomo cannot be used with recon.mvm. changed.\n");
			parms->recon.mvm=0;
		}
		if(parms->sim.closeloop==1){
			dbg2("idealtomo works in open loop only. changed.\n");
			parms->sim.closeloop=0;
		}
		if(parms->recon.modal){
			dbg2("idealtomo only works with zonal reconstruction. changed.\n");
			parms->recon.modal=0;
		}
		if(parms->sim.wfsalias){
			error("wfsalias and idealtomo conflicts\n");
			parms->sim.wfsalias=0;
		}
		if(parms->sim.idealwfs){
			error("idealwfs and idealtomo conflicts\n");
			parms->sim.idealwfs=0;
		}
		if((parms->ncpa.nsurf||parms->ncpa.ntsurf)&&!parms->ncpa.calib){
			warning("idealtomo require ncpa.calib to be enabled when there are surfaces. changed.\n");
			parms->ncpa.calib=1;
		}
		if(parms->evl.tomo){
			error("idealtomo and evl.tomo cannot be used together.\n");
			parms->evl.tomo=0;
		}
		parms->gpu.tomo=0;
	}

	if(parms->sim.wfsalias){
		if(parms->sim.idealwfs){
			error("sim.wfsalias conflicts with sim.idealwfs. Do not enable both.\n");
		}
		if(parms->sim.idealevl){
			error("sim.wfsalias conflicts with sim.idealevl. Do not enable both.\n");
		}
	}

	if(parms->sim.wfsalias||parms->sim.idealwfs||parms->sim.idealevl){
		parms->sim.dmproj=1;/*need dmproj */
	}
	/*if(parms->ncpa.calib && !(parms->nsurf || parms->ntsurf || parms->load.ncpa)){
	info2("No surface found. ncpa.calib is reset to 0.\n");
	parms->ncpa.calib=0;
	}*/
	if(parms->ncpa.calib&&!parms->ncpa.ndir){
		dbg("Using evaluation directions as ncpa calibration directions.\n");
		int ndir=parms->ncpa.ndir=parms->evl.nevl;
		dfree(parms->ncpa.thetax);
		dfree(parms->ncpa.thetay);
		dfree(parms->ncpa.wt);
		dfree(parms->ncpa.hs);
		parms->ncpa.thetax=dnew(ndir,1);
		parms->ncpa.thetay=dnew(ndir,1);
		parms->ncpa.wt=dnew(ndir,1);
		parms->ncpa.hs=dnew(ndir,1);
		dcp(&parms->ncpa.thetax,parms->evl.thetax);
		dcp(&parms->ncpa.thetay,parms->evl.thetay);
		dcp(&parms->ncpa.wt,parms->evl.wt);
		dcp(&parms->ncpa.hs,parms->evl.hs);
	}
}

/**
   postproc various WFS parameters based on other input information

   -# pixtheta if automatic
   -# whether doing GS0.
   -# whether doing phy.
   -# which wfs belongs to which powfs.
   -# necessary adjustments if outputing WFS PSF.
*/
static void setup_parms_postproc_wfs(parms_t *parms){
	if(parms->sim.evlol||parms->sim.idealtomo){
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			free_powfs_cfg(&parms->powfs[ipowfs]);
		}
		free(parms->powfs); parms->powfs=NULL;parms->npowfs=0;
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			dfree(parms->wfs[iwfs].wvlwts);
			free(parms->wfs[iwfs].sabad);
		}
		free(parms->wfs); parms->wfs=NULL;parms->nwfs=0;
		dfree(parms->cn2.pair);
	}
	if(parms->sim.evlol){
		for(int idm=0; idm<parms->ndm; idm++){
			free_dm_cfg(&parms->dm[idm]);
		}
		free(parms->dm); parms->dm=NULL; parms->ndm=0;
	}

	//Note that empty powfs have been removed in readcfg_wfs.
	parms->ittfpowfs=-1;
	parms->ittpowfs=-1;
	parms->itpowfs=-1;
	parms->ilgspowfs=-1;
	parms->nlgspowfs=0;
	parms->hipowfs_hsmin=INFINITY;
	parms->hipowfs_hsmax=INFINITY;
	parms->sim.dtrat_hi=-1;
	parms->sim.dtrat_lo=-1;//maximmum of all lo wfs
	parms->sim.dtrat_lo2=-1;//minimum of all lo wfs
	parms->sim.dtrat_lof=-1;
	parms->step_lo=-1;
	parms->step_hi=-1;

	if(!parms->nwfs||!parms->npowfs){
		return;
	}
	parms->hipowfs_hsmax=0;
	parms->hipowfs=lnew(parms->npowfs, 1);
	parms->lopowfs=lnew(parms->npowfs, 1);
	//first: process types of powfs that are not TWFS
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		powfs_cfg_t *powfsi=&parms->powfs[ipowfs];
		if(powfsi->skip==2) continue;
		if(powfsi->lo){
			P(parms->lopowfs, parms->nlopowfs)=ipowfs;
			parms->nlopowfs++;
			parms->nlowfs+=powfsi->nwfs;
			if(powfsi->trs==1){
				error("Low order wfs should not be tilt removed\n");
			}
			if(powfsi->gtype_sim==GTYPE_G&&powfsi->type==WFS_SH){
				warning("Low order powfs%d is using gtilt instead of ztilt in simulation. "
					"This is not recommended.\n", ipowfs);
			}
			if(parms->powfs[ipowfs].order==2&&parms->ittfpowfs==-1){
				parms->ittfpowfs=ipowfs;
			}
			if(parms->powfs[ipowfs].order==1&&parms->ittpowfs==-1){
				parms->ittpowfs=ipowfs;
			}
		} else{
			P(parms->hipowfs, parms->nhipowfs)=ipowfs;
			parms->nhipowfs++;
			parms->nhiwfs+=powfsi->nwfs;
			if(powfsi->hs<parms->hipowfs_hsmin){
				parms->hipowfs_hsmin=powfsi->hs;
			}
			if(powfsi->hs>parms->hipowfs_hsmax){
				parms->hipowfs_hsmax=powfsi->hs;
			}
		}
		if(powfsi->trs){
			if(!powfsi->llt){
				warning("powfs%d has tip/tilt removed but is not LGS.\n", ipowfs);
			}
			if(powfsi->lo){
				warning("powfs%d has tip/tilt removed but is not high order.\n", ipowfs);
			}
			parms->ntrpowfs++;
		} else{
			if(powfsi->llt){
				warning("powfs%d controls tip/tilt but is an LGS.\n", ipowfs);
			}
			parms->ntipowfs++;
		}
	}
	//2nd: process TWFS and link to LGS WFS.
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		powfs_cfg_t *powfsi=&parms->powfs[ipowfs];
		if(powfsi->skip!=2) continue;
		if(powfsi->dtrat==1){
			error("powfs[%d].dtrat=%d is invalid: truth wfs must have dtrat>1\n", ipowfs, powfsi->dtrat);
		}
		if(parms->itpowfs==-1){
			parms->itpowfs=ipowfs;
			if(parms->ilgspowfs!=-1){
				info("powfs%d is Truth WFS for powfs%d\n", ipowfs, parms->ilgspowfs);
			}
		} else{
			warning("powfs%d: Ignore additional TWFS\n", ipowfs);
		}
	}
	//3rd: determine WFS and control loop parameters
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		powfs_cfg_t *powfsi=&parms->powfs[ipowfs];
		pywfs_cfg_t *pycfg=powfsi->pycfg;
		llt_cfg_t *lltcfg=powfsi->llt;

		if(powfsi->dtrat<=0){
			error("powfs[%d].dtrat=%d is invalid.\n", ipowfs, powfsi->dtrat);
		}else if(powfsi->step<parms->sim.start 
			|| (parms->sim.end>0 && powfsi->step>=parms->sim.end)){//powfs is disabled.
			info("powfs %d is not used.\n", ipowfs);
		}else if(powfsi->step>0){/*round step to be multiple of dtrat. */
			powfsi->step=((powfsi->step+powfsi->dtrat-1)/powfsi->dtrat)*powfsi->dtrat;
		}
		if(parms->sim.wfsalias){
			powfsi->noisy=0;
			powfsi->phystep=-1;
			if(powfsi->type!=WFS_SH){
				error("sim.wfsalias is only supported for SHWFS\n");
			}
		}
		if(powfsi->dither||(parms->recon.petal&&parms->powfs[ipowfs].lo)||powfsi->type==WFS_PY){
			if(powfsi->phystep==-1){
				error("powfs%d: Physical optics mode is required for dithering or petaling control.\n", ipowfs);
			}
		}
		if(powfsi->phystep!=-1){/*round phystep to be multiple of dtrat. */
			if(powfsi->phystep<powfsi->step){
				powfsi->phystep=powfsi->step;
			}else{
				powfsi->phystep=((powfsi->phystep+powfsi->dtrat-1)/powfsi->dtrat)*powfsi->dtrat;
			}
		}
		/*Do we ever do physical optics.*/
		if(powfsi->phystep>=0&&(powfsi->phystep<parms->sim.end||parms->sim.end==0)){
			powfsi->usephy=1;
			parms->nphypowfs++;
			if(!parms->sim.closeloop){
				warning("powfs%d: Physical optics mode does not work well in open loop simulations.\n", ipowfs);
			}
		} else{
			powfsi->usephy=0;
		}
		if(powfsi->phystep>powfsi->step||parms->save.gradgeom){
			powfsi->needGS0=1;
		} else{
			powfsi->needGS0=0;
		}

		if(powfsi->qe&&powfsi->phystep>=0){
			//Check rne input.
			long pixpsay=powfsi->pixpsa;
			long pixpsax=powfsi->radpix;
			if(!pixpsax) pixpsax=pixpsay;

			if(NX(powfsi->qe)*NY(powfsi->qe)!=pixpsax*pixpsay){
				error("Input qe [%ldx%ld] does not match subaperture pixel [%ldx%ld]\n.",
					NX(powfsi->qe),NY(powfsi->qe),pixpsax,pixpsay);
			}
		}
		if(powfsi->dither){
			if(powfsi->dither<-1){
				powfsi->dither*=-1;
				if(powfsi->type==WFS_PY&&parms->recon.modal){
					info("Enable 2nd dithering mode when dither is <-1 for pywfs.\n");
					powfsi->dither_mmd=1;
				}
			}
		
			parms->dither=1;
			if(powfsi->dither==-1){//no dithering, just collect i0
				powfsi->dither_amp=0;
			} else if(powfsi->dither==1){//tip/tilt/arcsec->radian
				powfsi->dither_amp*=AS2RAD;
			} else {//zernike modes. micron-->meter
				powfsi->dither_amp/=1e6;
				if(powfsi->phytype_sim2==PTYPE_MF){
					error("Cannot build matched filter with dither>1.");
				}
			}
			//Convert all rate in unit of WFS frame rate
			//pllrat was already in WFS frame rate.
			powfsi->dither_npoint*=powfsi->dtrat;
			powfsi->dither_pllrat*=powfsi->dtrat;
			powfsi->dither_pllskip*=powfsi->dtrat;
			if(powfsi->dtrat>1){
				info("dither_npoint, pllrat, pllskip scaled by dtrat=%d",powfsi->dtrat);
			}
			powfsi->dither_ograt*=powfsi->dither_pllrat;
			powfsi->dither_ogskip=powfsi->dither_ogskip*powfsi->dither_pllrat+powfsi->dither_pllskip;
			//Convert all in simulation rate (sim.dt).
			if(powfsi->dither_ograt<=0||powfsi->dither_pllrat<=0){
				error("dither_ograt or _pllrat must be positive\n");
			}
			if(powfsi->dither_glpf==0) powfsi->dither_glpf=1;
		}

		/*link LLT with iwfs*/
		if(lltcfg){
			parms->nlgspowfs++;
			if(parms->ilgspowfs==-1){
				parms->ilgspowfs=ipowfs;
			} else{
				dbg_once("There are multiple LGS type. parms->ilgspowfs points to the first one\n");
			}
			if(lltcfg->fcfsm){
				lltcfg->epfsm=fc2lp(lltcfg->fcfsm, parms->sim.dt*powfsi->dtrat);
				if(lltcfg->fcfsm>0) info("powfs%d.llt.fsm: f0=%g Hz, LPF ep=%g\n", ipowfs, lltcfg->fcfsm, lltcfg->epfsm);
			}
		}
		if(lltcfg||powfsi->dither==1){//has FSM
			if(powfsi->epfsm<=0){
				real g=servo_optim_margin(parms->sim.dt, powfsi->dtrat, powfsi->alfsm,
					M_PI/4, powfsi->f0fsm, powfsi->zetafsm);
				powfsi->epfsm=g;
				info("powfs%d.epfsm is set to %g (auto)\n", ipowfs, g);
			}
		}

		if(powfsi->usephy){
			if(powfsi->neaextra){
				info("powfs%d: Adding extra NEA of %.2f mas\n",ipowfs,powfsi->neaextra);
			}
			if(powfsi->neamin){
				info("powfs%d: Limit minimum NEA to %.2f mas\n",ipowfs,powfsi->neamin);
			}
			if(powfsi->siglev<=0 || dmin(powfsi->siglevs)<=0){
				error("powfs%d: siglev must be positive\n", ipowfs);
			}
		}else if(powfsi->bkgrndfn){
			warning("powfs%d: there is sky background, but wfs is in geometric mode. The background is ignored.\n",ipowfs);
		}
		
		int calib_i0=(powfsi->type==WFS_SH&&powfsi->phytype_sim==1&&powfsi->usephy)
			&&!(parms->tomo.ahst_idealngs&&powfsi->skip);
		if(powfsi->ncpa_method==-1){//auto
			if(calib_i0){//mtch
				powfsi->ncpa_method=NCPA_I0;//default to 2
			} else{
				powfsi->ncpa_method=NCPA_G;
			}
		}
		if(powfsi->ncpa_method==NCPA_I0){
			if(!calib_i0){
				dbg("powfs%d: ncpa_method changed from 2 to 1 for non-matched filter mode.\n",ipowfs);
				powfsi->ncpa_method=NCPA_G;
			} else{
				powfsi->mtchstc=0;
			}
		}
	
		powfsi->phytype_sim1=powfsi->phytype_sim;//save value

		if(powfsi->skip!=2){//Only for none Truth WFS
			if(powfsi->type==WFS_PY&&lltcfg){
				error("Pyramid WFS is not available for LGS WFS\n");
			}
			if(powfsi->lo||!powfsi->trs){//has t/t measurement.
				if(parms->sim.dtrat_lo<0){
					parms->sim.dtrat_lo=powfsi->dtrat;
				} else if(parms->sim.dtrat_lo<powfsi->dtrat){
					parms->sim.dtrat_lo=powfsi->dtrat;
				}
				if(parms->sim.dtrat_lo2<0){
					parms->sim.dtrat_lo2=powfsi->dtrat;
				} else if(parms->sim.dtrat_lo2>powfsi->dtrat){
					parms->sim.dtrat_lo2=powfsi->dtrat;
				}
				if(powfsi->order>1){
					if(parms->sim.dtrat_lof<0||parms->sim.dtrat_lof>powfsi->dtrat){
						parms->sim.dtrat_lof=powfsi->dtrat;
					}
				}
				if(parms->step_lo<0||parms->step_lo>powfsi->step){
					parms->step_lo=powfsi->step;
				}
				if(powfsi->noisy){
					parms->sim.noisy_lo=1;
				}
			}
			if(!powfsi->lo && !powfsi->skip){//participate in high order recon
				if(parms->sim.dtrat_hi<0){
					parms->sim.dtrat_hi=powfsi->dtrat;
				} else if(parms->sim.dtrat_hi!=powfsi->dtrat){
					error("powfs.dtrat is invalid: all high order WFS must have the same dtrat.\n");
				}
				if(parms->step_hi<0){
					parms->step_hi=powfsi->step;
				} else if(parms->step_hi!=powfsi->step){
					error("powfs.step is invalid: all high order WFS must have the same enabling step\n");
				}
				if(powfsi->noisy){
					parms->sim.noisy_hi=1;
				}
			}
		}
		if(pycfg){
			pycfg->siglev=powfsi->siglev;
			pycfg->wvl=dref(powfsi->wvl);
			pycfg->wvlwts=dref(powfsi->wvlwts);
			pycfg->poke=parms->recon.poke;//How many meters to poke
			pycfg->D=parms->aper.d;
			pycfg->dsa=powfsi->dsa;
			pycfg->order=powfsi->order;
			pycfg->pixblur=powfsi->pixblur;
			pycfg->fieldstop=powfsi->fieldstop;
			pycfg->saat=powfsi->saat;
			if(pycfg->poke>1e-5||pycfg->poke<1e-10){
				warning("Pyramid WFS poke=%g m is out of the recommended range\n", pycfg->poke);
			}
		}
	}//for ipowfs

	lresize(parms->hipowfs, parms->nhipowfs, 1);
	lresize(parms->lopowfs, parms->nlopowfs, 1);
	if(parms->npowfs&&!parms->nhipowfs){
		warning("There is no high order WFS.\n");
	}
	if(parms->sim.dtrat_lo%parms->sim.dtrat_lo2!=0){
		error("Slower dtrat=%d has to be multiple of %d\n",parms->sim.dtrat_lo,parms->sim.dtrat_lo2);
	}
	dbg("dtrat_lo=%d, dtrat_lo2=%d, dtrat_lof=%d\n",parms->sim.dtrat_lo,parms->sim.dtrat_lo2,parms->sim.dtrat_lof);
	parms->sim.dtlo=parms->sim.dtrat_lo*parms->sim.dt;
	parms->sim.dthi=parms->sim.dtrat_hi*parms->sim.dt;
	if(parms->sim.fcfocus<0){
		parms->sim.fcfocus=0.1/parms->sim.dtlo;
	}

	parms->sim.lpfocushi=fc2lp(parms->sim.fcfocus,parms->sim.dthi);//active only when wfs has output.
	parms->sim.lpfocuslo=fc2lp(parms->sim.fcfocus,parms->sim.dt*parms->sim.dtrat_lof);
}

/**
   The siglev is always specified in sim.dtref. If sim.dt is different, rescale the siglev.
*/
static void setup_parms_postproc_dtref(parms_t *parms){
	real sigscale=parms->sim.dt>0?(parms->sim.dt/parms->sim.dtref):1;
	if(fabs(sigscale-1.)>EPS){
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			parms->wfs[iwfs].siglev*=sigscale;
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			parms->powfs[ipowfs].siglev*=sigscale;
			dscale(parms->powfs[ipowfs].siglevs,sigscale);
			parms->powfs[ipowfs].bkgrnd*=sigscale;
		}
	}

	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		sigscale=parms->powfs[ipowfs].sigscale;
		parms->wfs[iwfs].sigsim=parms->wfs[iwfs].siglev*sigscale;
		if(fabs(sigscale-1)>1.e-12){
			warning("wfs%d: siglev is scaled by %g to %g for simulation (not pixel processing).\n",
				iwfs,sigscale,parms->wfs[iwfs].sigsim);
		}
	}
}
/**
   postproc atmosphere parameters.
   1) drop weak layers.
   2) find ground layer
*/
static void setup_parms_postproc_atm(parms_t *parms){
	/*
	  Drop weak turbulence layers in simulation.

	  int jps=0;
	  for(int ips=0; ips<parms->atm.nps; ips++){
	  if(P(parms->atm.wt,ips)>1.e-4){
	  if(ips!=jps){
	  P(parms->atm.ht,jps)=P(parms->atm.ht,ips);
	  P(parms->atm.wt,jps)=P(parms->atm.wt,ips);
	  P(parms->atm.ws,jps)=P(parms->atm.ws,ips);
	  P(parms->atm.wddeg,jps)=P(parms->atm.wddeg,ips);
	  }
	  jps++;
	  }else{
	  warning("Layer %d has very small weight of %g, will drop it.\n",
	  ips, P(parms->atm.wt,ips));
	  }
	  }
	  if(jps==0){
	  error("There are no valid atmosphere layer\n");
	  }
	  if(parms->atm.nps!=jps){
	  warning("nlayer changed from %d to %d\n", parms->atm.nps, jps);
	  parms->atm.nps=jps;
	  dresize(parms->atm.ht, jps, 1);
	  dresize(parms->atm.wt, jps, 1);
	  dresize(parms->atm.ws, jps, 1);
	  dresize(parms->atm.wddeg, jps, 1);
	  }*/
	if(parms->sim.idealtomo){/*For ideal tomography, we using downsampled atm directly for atmr. */
		dbg("Changing atmr.ht,wt to atm.ht,wt since we are doing fit only\n");
		int nps=parms->atm.nps;
		dresize(parms->atmr.ht,nps,1);
		dresize(parms->atmr.wt,nps,1);
		dcp(&parms->atmr.ht,parms->atm.ht);
		dcp(&parms->atmr.wt,parms->atm.wt);
		lresize(parms->atmr.os,nps,1);
		for(int ips=parms->atmr.nps; ips<nps; ips++){
			P(parms->atmr.os,ips)=P(parms->atmr.os,parms->atmr.nps-1);
		}
		parms->atmr.nps=nps;
	}else if((parms->recon.glao)
		&&parms->recon.alg==RECON_MVR&&NX(parms->atmr.ht)>1&&!parms->sim.idealtomo){
   		//GLAO or single high wfs mode. reconstruct only a single layer near the DM.
		//Comment out because reconstructing into multiple layers can have better performance.
		dbg("In GLAO or single high wfs Mode, use 1 tomography grid near the ground dm.\n");
		dresize(parms->atmr.ht,1,1);
		dresize(parms->atmr.wt,1,1);
		lresize(parms->atmr.os,1,1);
		P(parms->atmr.ht,0)=parms->dm[0].ht;
		P(parms->atmr.wt,0)=1;
		parms->atmr.nps=1;
	}
	else if(parms->npowfs>0){
		int ipsr2=0;
		for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
			if(P(parms->atmr.ht,ipsr)>=parms->hipowfs_hsmax){
				dbg("Tomography Layer %d is above high order guide star and therefore dropped.\n",ipsr);
			} else{
				P(parms->atmr.ht,ipsr2)=P(parms->atmr.ht,ipsr);
				P(parms->atmr.wt,ipsr2)=P(parms->atmr.wt,ipsr);
				ipsr2++;
			}
		}
		if(ipsr2!=parms->atmr.nps){
			parms->atmr.nps=ipsr2;
			dresize(parms->atmr.ht,ipsr2,1);
			dresize(parms->atmr.wt,ipsr2,1);
		}
	}
	dnormalize_sumabs(parms->atm.wt, 1);
	dnormalize_sumabs(parms->atmr.wt, 1);

	/*
	  We don't drop weak turbulence layers in reconstruction. Instead, we make
	  it as least parms->tomo.minwt in setup_recon_tomo_reg
	*/
	if(!parms->recon.glao){
	/*Assign each turbulence layer to a corresponding reconstructon layer. Used
	  to compute opdx in a simple minded way.*/
		parms->atm.ipsr=lnew(parms->atm.nps,1);
		for(int ips=0; ips<parms->atm.nps; ips++){
			real dist=INFINITY;
			int kpsr=-1;
			real ht=P(parms->atm.ht,ips);
			for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
				real htr=P(parms->atmr.ht,ipsr);
				real dist2=fabs(ht-htr);
				if(dist2<dist){
					dist=dist2;
					kpsr=ipsr;
				}
			}
			P(parms->atm.ipsr,ips)=kpsr;
			dbg3("atm layer %d is maped to atmr %d\n",ips,kpsr);
		}

		/* Map reconstructed layers to input layers. for testing tomo.predict*/
		parms->atmr.indps=lnew(parms->atmr.nps,1);
		for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
			P(parms->atmr.indps,ipsr)=-1;
			for(int ips=0; ips<parms->atm.nps; ips++){
				if(fabs(P(parms->atmr.ht,ipsr)-P(parms->atm.ht,ips))<1e-3){
					if(P(parms->atmr.indps, ipsr)==-1){
						P(parms->atmr.indps, ipsr)=ips;
					}else{
						if(!parms->atm.dtrat){
							warning("One ipsr is mapped to multiple ips\n");
						}
					}
				}
			}
		}
	}
	/*
	  Find ground turbulence layer. The ray tracing can be shared between different directions.
	*/
	parms->atm.iground=-1;
	parms->atm.hmax=0;
	for(int ips=0; ips<parms->atm.nps; ips++){
		if(fabs(P(parms->atm.ht,ips))<1.e-10){
			if(parms->atm.iground==-1){
				parms->atm.iground=ips;
			} else if(!parms->atm.dtrat){
				dbg("There are multiple ground atm layers (OK).\n");
			}
		}
		if(P(parms->atm.ht,ips)<0){
			warning("Layer %d height %g is below ground (OK).\n",ips,P(parms->atm.ht,ips));
		}
		if(ips>0 && !parms->atm.dtrat && fabs(P(parms->atm.ht,ips)-P(parms->atm.ht,ips-1))<10){
			warning("Layer %d at %gm is very close to layer %d at %gm (OK).\n",
				ips,P(parms->atm.ht,ips),ips-1,P(parms->atm.ht,ips-1));
		}
		if(parms->atm.hmax<P(parms->atm.ht,ips)){
			parms->atm.hmax=P(parms->atm.ht,ips);
		}
	}
	parms->atmr.hmax=0;
	for(int ips=0; ips<parms->atmr.nps; ips++){
		if(parms->atmr.hmax<P(parms->atmr.ht,ips)){
			parms->atmr.hmax=P(parms->atmr.ht,ips);
		}
	}
	if(parms->sim.closeloop && !parms->atm.frozenflow){
		error("sim.closeloo=1 requires atm.frozenflow=1.\n");
		parms->atm.frozenflow=1;
	}

	if(!parms->atm.frozenflow||parms->dbg.atm){
		parms->atm.r0evolve=0;/*disable r0 evolution*/
	}

	if(!parms->atm.frozenflow){
		if(parms->sim.end>parms->sim.start+10){
			dbg("Disable turbulence file based sharing in open loop nonfrozenflow simulation\n");
			parms->atm.share=0;
		}
		parms->sim.dt=0;
	}

	if(!parms->atmr.hs){
		real hs=NAN;
		/*find out the height to setup cone coordinate. */
		if(parms->tomo.cone){
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			/*skip wfs that does not participate in tomography*/
				if(parms->powfs[ipowfs].lo||parms->powfs[ipowfs].skip){
					continue;
				}
				/*isinf and !isinf both return 0 on inf in FreeBSD 9.0.*/
				if(isnan(hs)){
					hs=parms->powfs[ipowfs].hs;
				} else{
					if(!isinf(hs)||!isinf(parms->powfs[ipowfs].hs)){
						if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
							warning("Two high order POWFS with different hs found: %g and %g. Please double check the results.\n",
								hs,parms->powfs[ipowfs].hs);
							if(parms->powfs[ipowfs].hs>hs){
								hs=parms->powfs[ipowfs].hs;
							}
						}
					}
				}
			}
		}
		if(isnan(hs)) hs=INFINITY;
		parms->atmr.hs=hs;
	}
	if(parms->atmr.dx<EPS){
	/*find out the sampling to setup tomography grid using the maximum order of the wfs and DMs. */
		real mindsa=INFINITY;
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].lo||parms->powfs[ipowfs].skip){
				continue;
			}
			if(parms->powfs[ipowfs].dsa<mindsa){
				mindsa=parms->powfs[ipowfs].dsa;
			}
		}
		if(isinf(mindsa)){
			for(int idm=0; idm<parms->ndm; idm++){
				if(parms->dm[idm].dx<mindsa){
					mindsa=parms->dm[idm].dx;
				}
			}
		}
		if(isinf(mindsa)){
			for(int imoao=0; imoao<parms->nmoao; imoao++){
				if(parms->moao[imoao].dx<mindsa){
					mindsa=parms->moao[imoao].dx;
				}
			}
		}
		if(isinf(mindsa)){
			error("Unable to determine atmr.dx. Please specify manually atmr.dx\n");
		}
		parms->atmr.dx=mindsa;
	}



}
static void setup_parms_postproc_dirs(parms_t *parms){
	//Collect all beam directions
	const int ndir=parms->nwfs+parms->evl.nevl+parms->fit.nfit+(parms->ncpa.calib?parms->ncpa.ndir:0);
	parms->dirs=dnew(4,ndir);
	dmat *pdir=parms->dirs/*PDMAT*/;
	int count=0;

	for(int i=0; i<parms->nwfs; i++){
		P(pdir,0,count)=parms->wfs[i].thetax;
		P(pdir,1,count)=parms->wfs[i].thetay;
		P(pdir,2,count)=parms->wfs[i].hs;
		P(pdir,3,count)=parms->wfs[i].hc;
		count++;
	}

	for(int i=0; i<parms->evl.nevl; i++){
		P(pdir,0,count)=P(parms->evl.thetax,i);
		P(pdir,1,count)=P(parms->evl.thetay,i);
		P(pdir,2,count)=P(parms->evl.hs,i);
		count++;
	}
	for(int i=0; i<parms->fit.nfit; i++){
		P(pdir,0,count)=P(parms->fit.thetax,i);
		P(pdir,1,count)=P(parms->fit.thetay,i);
		P(pdir,2,count)=P(parms->fit.hs,i);
		count++;
	}
	if(parms->ncpa.calib){
		for(int i=0; i<parms->ncpa.ndir; i++){
			P(pdir,0,count)=P(parms->ncpa.thetax,i);
			P(pdir,1,count)=P(parms->ncpa.thetay,i);
			P(pdir,2,count)=P(parms->ncpa.hs,i);
			count++;
		}
	}
	if(count<ndir){
		warning("count=%d, ndir=%d\n",count,ndir);
	} else if(count>ndir){
		error("count=%d, ndir=%d\n",count,ndir);
	}
	dresize(parms->dirs,4,count);
	real rmax=0;
	for(int ic=0; ic<count; ic++){
		real x=P(parms->dirs,0,ic);
		real y=P(parms->dirs,1,ic);
		real r=sqrt(x*x+y*y);
		if(isfinite(P(parms->dirs, 2, ic))){//cone effect
			r-=parms->aper.d/(P(parms->dirs, 2, ic)*2);
	}
		if(r>rmax) rmax=r;
	}
	real fov=2*rmax;
	if(parms->sim.fov<fov){
			dbg("sim.fov=%g is less than actual fov=%g. Changed\n",parms->sim.fov*RAD2AS,fov*RAD2AS);
		parms->sim.fov=fov;
	}
}
/**
   compute minimum size of atm screen to cover all the beam path. same for
   all layers.  todo:may need to consider L0 Must be after
   setup_parms_postproc_za.
*/
static void setup_parms_postproc_atm_size(parms_t *parms){
	const int nps=parms->atm.nps;
	int Nmax=0;
	long nxout[nps],nyout[nps];
	parms->atm.nxn=lnew(nps,1);
	for(int ips=0; ips<nps; ips++){
		double guard=parms->atm.dx*3;
		create_metapupil(0,&nxout[ips],&nyout[ips],parms->dirs,parms->aper.d,P(parms->atm.ht,ips),
			parms->atm.dx,parms->atm.dx,0.5,guard,0,0,0,1);
		P(parms->atm.nxn,ips)=MAX(nxout[ips],nyout[ips]);
		Nmax=MAX(Nmax,P(parms->atm.nxn,ips));
	}
	/*Minimum screen size required. Used to transport atm to GPU. */
	parms->atm.nxnmax=Nmax;
	Nmax=nextpow2(Nmax);
	if(!P(parms->atm.size,0)||!P(parms->atm.size,1)){
		parms->atm.nx=Nmax;
		parms->atm.ny=Nmax;
	} else{/*user specified.*/
		parms->atm.nx=2*(int)round(0.5*P(parms->atm.size,0)/parms->atm.dx);
		parms->atm.ny=2*(int)round(0.5*P(parms->atm.size,1)/parms->atm.dx);
		if(parms->atm.nx<Nmax) parms->atm.nx=Nmax;
		if(parms->atm.ny<Nmax) parms->atm.ny=Nmax;
	}
	if(parms->atm.method==1){/*must be square and 1+power of 2 */
		int nn=parms->atm.nx>parms->atm.ny?parms->atm.nx:parms->atm.ny;
		parms->atm.nx=1+nextpow2(nn);
		parms->atm.ny=parms->atm.nx;
	}
	/*record the size of the atmosphere. */
	real atm_size_x=P(parms->atm.size,0);
	real atm_size_y=P(parms->atm.size,1);
	P(parms->atm.size,0)=parms->atm.nx*parms->atm.dx;
	P(parms->atm.size,1)=parms->atm.ny*parms->atm.dx;
	if(P(parms->atm.size,0)>atm_size_x||P(parms->atm.size,1)>atm_size_y){
		info("Atmospheric size is increased from (%g, %g) to (%g, %g) m.\n",
			atm_size_x,atm_size_y,P(parms->atm.size,0),P(parms->atm.size,1));
	}
	if(P(parms->atm.L0,0)>P(parms->atm.size,0)){
		info("Atmospheric size is smaller than outer scale.\n");
	}
	/*for screen evolving. */
	parms->atm.overx=lnew(parms->atm.nps,1);
	parms->atm.overy=lnew(parms->atm.nps,1);
	for(int ips=0; ips<parms->atm.nps; ips++){
		P(parms->atm.overx,ips)=nxout[ips];
		P(parms->atm.overy,ips)=nyout[ips];
	}
}
/*
  Find entry  that equals to val in array of length n. Append if not exist yet.
*/
/*
static int arrind(real *arr, int *n, real val){
	for(long i=0; i<(*n); i++){
	if(fabs(arr[i]-val)<EPS){
		return i;
	}
	}
	arr[*n]=val;
	(*n)++;
	return (*n)-1;
}
*/
/**
   Setting up DM parameters in order to do DM caching during simulation. High
   resolution metapupils are created for each DM at each magnitude level to
   match the science field or WFS.

   For a MCAO system like NFIRAOS, the WFS and science all have 1/64
   sampling. The ground DM has only 1 metapupil matched to the science field and
   WFS. The upper DM has two metapupils, one for science and NGS and another one
   to match LGS WFS in a reduce sampling depending on the DM altitude and guide
   star range..

*/
static void setup_parms_postproc_dm(parms_t *parms){
	/*disable cache for low order systems. */
	if(parms->sim.cachedm){
		if(parms->evl.nevl<2&&parms->nwfs<2){
			parms->sim.cachedm=0;
			warning("sim.cachedm disabled for SCAO\n");
		}
		if(parms->dbg.cmpgpu){
			parms->sim.cachedm=0;
			warning("sim.cachedm disabled when comparing CPU against GPU\n");
		}
	}
	for(int i=0; i<parms->ndm; i++){
		real ht=parms->dm[i].ht+parms->dm[i].vmisreg;
		if(fabs(ht)<1.e-10){
			parms->dm[i].isground=1;
			parms->idmground=i;
		}
		if(!isinf(P(parms->dm[i].stroke,0))){
			real strokemicron=fabs(P(parms->dm[i].stroke,0))*1e6;
			if(strokemicron>1000){
				dscale(parms->dm[i].stroke,1e-6);
				//info("dm %d: assume stroke=%g is in micron.\n",i,strokemicron*1e-6);
			}
			if(parms->sim.fcttm==0){
				warning("Please set sim.fcttm for offloading to tip/tilt mirror when DM stroke is limited.\n");
			}
		}
		if(!isinf(parms->dm[i].iastroke)&&!parms->dm[i].strokescale){
			real strokemicron=parms->dm[i].iastroke*1e6;
			if(strokemicron>1000){
				parms->dm[i].iastroke*=1e-6;
				//info("dm %d: assume ia stroke=%g is in micron.\n",i,strokemicron*1e-6);
			}
		}
	}
}

/**
   Setting up the cone coordinate for MCAO LGS
   simulation. First find out the guide star conjugate. Only 1
   altitude is allowed.
*/
static void setup_parms_postproc_recon(parms_t *parms){
	{
		//moved from postproc_wfs to ensure initialization
		parms->sim.lpttm=fc2lp(parms->sim.fcttm, parms->sim.dt);//active at every time step. use dt

		if(P(parms->sim.ephi, 0)<=0){
			real g=0.5;
			if(parms->npowfs>0){
				g=servo_optim_margin(parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi,
					M_PI/4, parms->sim.f0dm, parms->sim.zetadm);
			}
			P(parms->sim.ephi, 0)=g;
			info("sim.ephi is set to %g (auto)\n", g);
		}
		if(P(parms->sim.eplo, 0)<=0){
			real g=0.5;
			if(parms->npowfs>0){
				servo_optim_margin(parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo,
				M_PI/4, parms->sim.f0dm, parms->sim.zetadm);//dm is used for tweeter t/t control.
			}
			P(parms->sim.eplo, 0)=g;
			info("sim.eplo is set to %g (auto)\n", g);
		}
	}
	if(parms->nmoao>0){//remove unused moao configurations
		int count=0;
		for(int imoao=0; imoao<parms->nmoao; imoao++){
			int used=0;
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				if(parms->powfs[ipowfs].moao>=parms->nmoao){
					error("invalid powfs[%d].moao=%d",ipowfs,parms->powfs[ipowfs].moao);
				}
				if(parms->powfs[ipowfs].moao==imoao){
					used=1;
					break;
				}
			}
			if(parms->evl.moao>=parms->nmoao){
				error("invalid evl.moao=%d\n",parms->evl.moao);
			}
			if(parms->evl.moao==imoao){
				used=1;
			}
			if(used){
				parms->moao[imoao].used=1;
				count++;
			}
		}
		if(count==0){//no moao
			parms->nmoao=0;
			free(parms->moao);
			parms->moao=NULL;
		}
	}
	if(parms->recon.split){
		if(!parms->nwfs||!parms->sim.closeloop || !parms->ndm || parms->evl.tomo||parms->sim.evlol){
			parms->recon.split=0;
			dbg("Split tomography is not support or needed in current configuration. Changed.\n");
		}
	}
	if(parms->recon.split){
		parms->evl.split=1;//force 1
	}else if(parms->sim.evlol){
		parms->evl.split=0;//force 0
	}
	if(parms->recon.glao&&parms->ndm!=1){
		error("GLAO only works with 1 dm\n");
	}
	/*if(parms->recon.alg==RECON_MVR&&parms->recon.modal){
		error("Modal control is not supported yet with MV reconstructor. Consider change to LSR reconstructor or disable modal control.\n");
		parms->recon.modal=0;
	}*/
	if(parms->recon.alg==RECON_LSR){
		if(parms->recon.split==2){
			error("MVST does not work with least square reconstructor.\n");
		}
		if(parms->lsr.alg==2){
			parms->recon.mvm=1;
		}
		if(parms->lsr.actextrap==-1){
			if(parms->recon.modal){
				//no need in modal lsr control
				parms->lsr.actextrap=0;
			} else{
				parms->lsr.actextrap=1;
			}
		}
	}

	if(parms->fit.square&&parms->load.aloc){
		warning("load.aloc contradicts with fit.square. disable fit.square\n");
		parms->fit.square=0;
	}
	if(!parms->sim.closeloop){
		parms->recon.psol=0;//open loop do not need psol
		parms->cn2.psol=0;//open loop do not need psol
	} else if(parms->recon.psol==-1){
		if(parms->sim.idealtomo){
			parms->recon.psol=1;
		} else if(parms->recon.alg==RECON_MVR){//MV perfers psol
			parms->recon.psol=1;
		} else{//LSR perfers cl
			parms->recon.psol=0;
		}
	}
	
	if(parms->recon.split){
		if(parms->nlopowfs==0){
			if(parms->ntrpowfs>=parms->nhipowfs){
				warning("There is no WFS controlling tip/tilt.\n");
			} else{
				info("Split reconstruction is enabled when there is no low order WFS."
					" Will split the tip/tilt modes from high order wfs\n");
			}
		}
		if(parms->sim.skysim&&(parms->nhipowfs==0||parms->nlopowfs==0)){
			error("There is only high or low order WFS. can not do skycoverage presimulation\n");
		}
		if(parms->tomo.ahst_wt==1){
			if(parms->tomo.ahst_idealngs||!parms->ntipowfs||P(parms->sim.eplo,0)<0.01){
				dbg("Change tomo.ahst_wt from 1 to 3 when there is no NGS WFS control.\n");
				parms->tomo.ahst_wt=3;
			}
		}
		/*if(parms->ndm>2 && parms->tomo.ahst_wt!=4){
			dbg("Change tomo.ahst_wt from %d to 4 when there are more than 2 DMs.\n", parms->tomo.ahst_wt);
			parms->tomo.ahst_wt=4;
		}*/
	}else if(parms->evl.split){
		if(parms->tomo.ahst_wt!=3){
			dbg("Change tomo.ahst_wt from %d to 3 when there is no ahst control.\n", parms->tomo.ahst_wt);
			parms->tomo.ahst_wt=3;
		}
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].r0<=0){
			parms->powfs[ipowfs].r0=parms->atm.r0;
		}
		if(parms->powfs[ipowfs].L0<=0){
			parms->powfs[ipowfs].L0=dsum(parms->atm.L0)/parms->atm.nps;
		}
		if(parms->recon.split&&parms->powfs[ipowfs].lo){
			parms->powfs[ipowfs].skip=1;
		}

		if(parms->save.ngcov>0||(parms->cn2.pair&&!parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].skip)){
			/*focus tracking or cn2 estimation, or save gradient covariance.  */
			parms->powfs[ipowfs].psol=1;
		} else if(parms->recon.psol){//PSOL reconstruction
			/*low order wfs in ahst mode does not need psol. */
			if((parms->recon.split==1&&parms->powfs[ipowfs].skip)){
				parms->powfs[ipowfs].psol=0;
			} else{
				parms->powfs[ipowfs].psol=1;
			}
		}
	}

	if(parms->recon.glao){
		parms->wfsr=mycalloc(parms->npowfs,wfsr_cfg_t);
		parms->nwfsr=parms->npowfs;
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			parms->wfsr[ipowfs].thetax=0;
			parms->wfsr[ipowfs].thetay=0;
			parms->wfsr[ipowfs].hs=parms->powfs[ipowfs].hs;
			//parms->wfsr[ipowfs].hc=parms->powfs[ipowfs].hc;
			parms->wfsr[ipowfs].powfs=ipowfs;
			parms->powfs[ipowfs].nwfsr=1;
			parms->powfs[ipowfs].wfsr=lnew(1,1);
			P(parms->powfs[ipowfs].wfsr,0)=ipowfs;
		}
		/*
		  parms->fit.nfit=1;
		  dresize(parms->fit.thetax, 1, 1);
		  dresize(parms->fit.thetay, 1, 1);
		  dresize(parms->fit.wt, 1, 1);
		  dresize(parms->fit.hs, 1, 1);
		  P(parms->fit.thetax,0)=0;
		  P(parms->fit.thetay,0)=0;
		  P(parms->fit.wt,0)=1;
		*/
	} else{/*Use same information as wfs. */
		parms->wfsr=mycalloc(parms->nwfs, wfsr_cfg_t);
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			parms->wfsr[iwfs].thetax=parms->wfs[iwfs].thetax;
			parms->wfsr[iwfs].thetay=parms->wfs[iwfs].thetay;
			parms->wfsr[iwfs].hs=parms->wfs[iwfs].hs;
			parms->wfsr[iwfs].misregx=parms->wfs[iwfs].misregx;
			parms->wfsr[iwfs].misregy=parms->wfs[iwfs].misregy;
			parms->wfsr[iwfs].misregc=parms->wfs[iwfs].misregc;
			parms->wfsr[iwfs].powfs=parms->wfs[iwfs].powfs;
		}
		parms->nwfsr=parms->nwfs;
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			parms->powfs[ipowfs].nwfsr=parms->powfs[ipowfs].nwfs;
			parms->powfs[ipowfs].wfsr=lref(parms->powfs[ipowfs].wfs);
		}
	}
	if(parms->recon.alg==RECON_MVR){//MVM: tomo+fit
		if(parms->tomo.alg==-1){//default to CG
			parms->tomo.alg=ALG_CG;
		}
		if(parms->tomo.alg==ALG_CG){
			if(parms->nhipowfs>1){
				if(parms->tomo.precond==PCG_FD){
					info("Disable FDPCG when there are multiple high order powfs.\n");
					parms->tomo.precond=PCG_NONE;
				}
			}else if(parms->tomo.precond==-1){
				parms->tomo.precond=PCG_FD;
			}
		}else{
			parms->tomo.precond=PCG_NONE;
		}
		if(parms->fit.alg==-1){
			if(parms->recon.modal){
				parms->fit.alg=ALG_CG;
			}else{
				parms->fit.alg=parms->recon.mvm?ALG_CBS:ALG_CG;//MVM is only good with CBS or SVD.
			}
		}else if(parms->recon.modal && parms->fit.alg==ALG_CBS){
			warning("recon.modal cannot work with CBS. Change to SVD.\n");
			parms->fit.alg=ALG_SVD;
		}
	}
	if(parms->atm.frozenflow&&!parms->dbg.nocgwarm){
		parms->fit.cgwarm=1;
		parms->lsr.cgwarm=1;
		if(!NX(parms->dbg.tomo_maxit)){
			parms->tomo.cgwarm=1;
		}
	}//else: no warm

	if(parms->recon.split==1&&!parms->sim.closeloop&&parms->ndm>1){
		warning("ahst split tomography does not have good NGS correction in open loop.\n");
	}
	if(parms->recon.split==2&&parms->sim.fuseint){
		warning("MVST Mode can only use separate integrator for the moment. Changed.\n");
		parms->sim.fuseint=0;
	}
	if(!parms->recon.split&&!parms->sim.fuseint){
		parms->sim.fuseint=1;/*integrated tomo. only 1 integrator. */
	}

	if(parms->recon.alg==RECON_MVR){
		if(parms->recon.mvm&&parms->fit.alg==ALG_CG){
			warning("CG based fit is not suitable for building MVM\n");
		}
		if(parms->sim.ecnn||parms->load.tomo||parms->tomo.alg!=ALG_CG||parms->tomo.bgs){
			parms->tomo.assemble=1;
		}
		if((parms->tomo.bgs||parms->tomo.alg!=ALG_CG)&&parms->tomo.cxxalg!=0){
			error("Only CG work with non L2 cxx.\n");
			parms->tomo.cxxalg=0;
		}
		if(parms->tomo.predict==1&&parms->tomo.alg!=ALG_CG){
			error("Predictive tomography only works with CG. need to redo CBS/MVM after wind velocity is know.\n");
		}
		if(parms->tomo.alg==ALG_CG){/*MVR with CG*/
			if(parms->tomo.precond>PCG_TOT){
				error("Invalid preconditoner\n");
			}
		}

		/*Fitting tip/tilt constraint is only intended for multi DM*/
		if(parms->ndm<2&&parms->fit.lrt_tt){
			parms->fit.lrt_tt=0;
		}
		if(parms->ndm>2&&parms->fit.lrt_tt==2){
			warning("When there are more than 2 DMs, lrt_tt has to be 1 instead of 2. changed\n");
			parms->fit.lrt_tt=1;
		}
		if(parms->fit.lrt_tt<0||parms->fit.lrt_tt>2){
			error("parms->fit.lrt_tt=%d is invalid\n",parms->fit.lrt_tt);
		}
		/*Assign CG interations*/
		if(parms->tomo.alg==ALG_CG&&parms->tomo.maxit<=0){
			int maxit=4;//minimal 4 iterations is needed
			if(parms->recon.mvm){
				maxit*=parms->load.mvmi?1:25;//assembly mvm needs more steps
			} else{
				maxit*=parms->tomo.cgwarm?1:10;
			}
			if(parms->tomo.precond==PCG_NONE){
				maxit*=10;//non-precond CG needs more steps
			}
			if(!parms->recon.split){
				maxit*=4;//integrated tomo needs more steps
			}
			if(parms->ndm>1){
				//if meta pupil is much larger than the aperture, needs more iterations. 
				//if atmr.dx is smaller than 0.5, also more iterations
				real ratio=pow(1+parms->sim.fov*parms->atmr.hmax/parms->aper.d,2);//meta pupil to pupil area ratio
				if(parms->atmr.dx<0.5){
					ratio*=0.5/parms->atmr.dx;
				}
				if(ratio>10) ratio=10;//limit maximum value.
				maxit=ceil(maxit*ratio);
			}
			parms->tomo.maxit=maxit;
			if(parms->recon.mvm==1&&parms->recon.split&&parms->tomo.splitlrt){
				warning("recon.mvm==1 require tomo.splitlrt=0 due to stability issue. Changed\n");
				parms->tomo.splitlrt=0;
			}
		}
		if(parms->tomo.bgs&&parms->tomo.precond){
			error("Please implement the preconditioner for each block for BGS.\n");
		}
	}

	/*DM Fitting related. fit parameters are also used for dmproj.*/
	if(parms->fit.alg==ALG_CG&&parms->fit.maxit<=0){
		int factor;
		factor=parms->fit.cgwarm?1:10;
		parms->fit.maxit=10*factor;
	}

	if(parms->load.fit||parms->fit.alg!=ALG_CG||parms->fit.bgs){
		parms->fit.assemble=1;
	}
	if(parms->fit.bgs&&parms->fit.precond){
		error("Please implement the preconditioner for each block for BGS.\n");
	}
	if(parms->fit.pos<=0) parms->fit.pos=parms->tomo.pos;
	if(parms->recon.alg==RECON_LSR){
	/*if(parms->lsr.actslave>1 && parms->lsr.tikcr>0){
		info2("lsr.actslave>1 disables lsr.tikcr\n");
		parms->lsr.tikcr=0;
		}*/
		if(parms->lsr.alg==1&&parms->lsr.maxit<=0){
			int factor;
			factor=parms->lsr.cgwarm?1:10;
			parms->lsr.maxit=30*factor;
		}
	}
	if(parms->sim.mffocus==-1){
		parms->sim.mffocus=(parms->nlgspowfs)?1:0;
	}
	if(parms->sim.mffocus>0){
		if(!parms->recon.split||!parms->nlgspowfs){
			info("Focus blending is only implemented for LGS in split tomography. Changed.\n");
			parms->sim.mffocus=0;
		}
	}

	if(parms->sim.mffocus<0||parms->sim.mffocus>2){
		error("parms->sim.mffocus=%d is invalid\n",parms->sim.mffocus);
	}
	if(parms->tomo.ahst_focus){
		if(parms->recon.split!=1||!parms->sim.mffocus){
			parms->tomo.ahst_focus=0;//no need ahst_focus
			dbg("Disable tomo.ahst_focus.\n");
		}
	}

	if(!parms->recon.mvm){
		if(parms->tomo.alg!=ALG_CG&&parms->load.mvmi){
			free(parms->load.mvmi);
			parms->load.mvmi=NULL;
		}
		if(parms->load.mvmf){
			free(parms->load.mvmf);
			parms->load.mvmf=NULL;
		}
	}

	for(int idm=0; idm<parms->ndm; idm++){
		if(!isinf(P(parms->dm[idm].stroke,0))){
			parms->sim.dmclip=1;
		}
		if(!isinf(parms->dm[idm].iastroke)&&parms->dm[idm].iastroke>0){
			parms->sim.dmclipia=1;
		}
	}
	if(parms->sim.psfr){
		int fnd=lsum(parms->evl.psfr);
		if(fnd==0){
			error("sim.psfr is specified, but evl.psfr are all zero\n");
		} else{
			info("Output PSF reconstruction telemetry for %d directions\n",fnd);
		}
		if(!parms->evl.psfmean){
			parms->evl.psfmean=1;/*Saves psfmean for verification. */
		}
		if(!parms->save.ecov){
			parms->save.ecov=1;
		}
		/*required memory to hold memory. */
		long covmem=(long)round(pow(parms->aper.d/parms->evl.dx,4))*8*fnd;
		if(covmem>MAX(NMEM,LONG_MAX/2)&&parms->evl.dx>parms->atmr.dx*0.25+EPS){/*4G or actual */
			warning("parms->evl.dx=%g is probably too large to save ecxx. Recommend parms->evl.dx=%g\n",parms->evl.dx,parms->atmr.dx*0.25);
		}
	}

	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		if(parms->recon.distortion_tel2wfs&&parms->recon.distortion_tel2wfs[iwfs]&&!parms->dbg.tomo_hxw){
			warning_once("Without dbg.tomo_hxw, only pure shift between telescope and LGS WFS is calibrated.\n");
		}
	}
	if(!parms->nwfs||parms->sim.noatm){
		parms->recon.psd=0;
	}

	/*
		It is OK to tune gain without noise. the algorithm doesnot need noise information.
	if(parms->recon.psd){
		if(parms->recon.psddtrat_hi&&!parms->sim.noisy_hi){
			parms->recon.psddtrat_hi=0;
		}
		if(parms->recon.psddtrat_lo&&!parms->sim.noisy_lo){
			parms->recon.psddtrat_lo=0;
		}
		if(!parms->recon.psddtrat_hi&&!parms->recon.psddtrat_lo){
			parms->recon.psd=0;
		}
	}*/
	switch(parms->recon.twfs_rmin){
	case 1:
		parms->itwfssph=parms->recon.twfs_radonly?1:9; break;
	case 2:
		parms->itwfssph=parms->recon.twfs_radonly?1:7; break;
	case 3:
		parms->itwfssph=parms->recon.twfs_radonly?0:4; break;
	default:
		parms->itwfssph=-1;
		error("Please implement\n");
	}
}

static void setup_parms_postproc_seeds(parms_t *parms, int override){
	if(disable_save) return;//do not check if not saved
	if(parms->sim.end>parms->sim.start){
		/*Remove seeds that are already done. */
		int iseed=0;
		int jseed=0;
		parms->fdlock=mycalloc(parms->sim.nseed, int);
		parms->fnlock=mycalloc(parms->sim.nseed, char *);
		char cwd[PATH_MAX];
		if(!getcwd(cwd, PATH_MAX)){
			cwd[0]='.'; cwd[1]='0';
		} else{
			for(char *p=cwd; p[0]; p++){
				if(*p=='/') *p='!';
			}
		}

		for(iseed=0; iseed<parms->sim.nseed; iseed++){
			char fn[PATH_MAX];//may not be NFS shared
			snprintf(fn, sizeof(fn), "Res_%ld.done", P(parms->sim.seeds, iseed));
			if(exist(fn)){
				if(override){
					remove(fn);
				} else{
					parms->fdlock[iseed]=-1;
					warning("Skip seed %ld because %s exists.\n", P(parms->sim.seeds, iseed), fn);
				}
			}
			if(!parms->fdlock[iseed]){
				snprintf(fn, sizeof(fn), "%s/%s_maos_%ld.lock", DIRLOCK, cwd, P(parms->sim.seeds, iseed));
				parms->fdlock[iseed]=lock_file(fn, 0);
				if(parms->fdlock[iseed]<0){
					warning("Skip seed %ld because it is already running.\n",
						P(parms->sim.seeds, iseed));
				} else{
					cloexec(parms->fdlock[iseed]);
					parms->fnlock[iseed]=mystrdup(fn);
					if(jseed!=iseed){//remove gap in array.
						P(parms->sim.seeds, jseed)=P(parms->sim.seeds, iseed);
						parms->fdlock[jseed]=parms->fdlock[iseed];
					}
					jseed++;
				}
			}
		}
		if(jseed!=parms->sim.nseed){
			parms->sim.nseed=jseed;
		}
	}//if sim.end>sim.start
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
	if(parms->sim.nseed<1){
		remove(fn);
		info2("There are no seed to run. Use -O to override. Exit\n");
	} else{
		mysymlink(fn, "run_recent.log");
		info2("There are %d valid simulation seeds: ", parms->sim.nseed);
		for(int i=0; i<parms->sim.nseed; i++){
			info2(" %ld", P(parms->sim.seeds, i));
		}
		info2("\n");
	}
}
/**
   postproc misc parameters.
*/
static void setup_parms_postproc_misc(parms_t *parms){
	if(disable_save){
		if(parms->save.setup||parms->save.all||parms->sim.skysim||parms->evl.psfmean||parms->evl.psfhist){
			error("Please specify -o DIR to enable saving to disk\n");
		}
	}
	if(parms->save.ngcov>0&&parms->save.gcovp<10){
		warning("parms->save.gcovp=%d is too small. It may fill your disk!\n", parms->save.gcovp);
	}
	if(parms->save.gcovp>parms->sim.end){
		parms->save.gcovp=parms->sim.end;
	}
	if(!parms->tomo.square){
		if(parms->dbg.cmpgpu){
			info("Make tomo.square=1 to compare cpu code against gpu implementations.\n");
			parms->tomo.square=1;
		}
	}
	if(!parms->atm.frozenflow){
		info("psfisim is set from %d to %d in openloop mode\n",parms->evl.psfisim,parms->sim.start);
		parms->evl.psfisim=parms->sim.start;
	}
	if(parms->evl.psfisim<parms->sim.start){
		parms->evl.psfisim=parms->sim.start;
	}

	if(parms->evl.psfmean||parms->evl.psfhist){
		int fnd=lsum(parms->evl.psf);
		if(fnd==0){
			warning("Required to output PSF, but evl.psf are all zero\n");
		} else{
			info("Output PSF for %d directions\n",fnd);
		}
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		parms->evl.npsf+=(P(parms->evl.psf,ievl)>0);
		if((P(parms->evl.psf,ievl)&2)){//request ngs mode removed PSF
			if(!parms->recon.split||!isinf(P(parms->evl.hs, ievl))){
				P(parms->evl.psf, ievl)^=2;
			}
		}
		if(parms->tomo.ahst_idealngs==1){
			//Output NGS mode removed PSF as there is no CL control of NGS mode
			if(!(P(parms->evl.psf,ievl)&2)){
				P(parms->evl.psf,ievl)^=2;
			}
		}
	}
	if(parms->dbg.dmoff){
		if(NX(parms->dbg.dmoff)!=parms->ndm){
			dcellresize(parms->dbg.dmoff,parms->ndm,NY(parms->dbg.dmoff));
		}
	}
	if(parms->dbg.gradoff){
		if(NX(parms->dbg.gradoff)!=parms->nwfs){
			dcellresize(parms->dbg.gradoff,parms->nwfs,NY(parms->dbg.gradoff));
		}
	}
}
static void print_alg(int bgs, int alg, int maxit, int precond, real svdthres){
	if(bgs){
		info2("Block Gauss Seidel with ");
	}
	switch(alg){
	case ALG_CBS:
		info2("Cholesky back solve ");
		break;
	case ALG_CG:
		info2("CG%d ", maxit);
		switch(precond){
			case 0:	break;
			case 1:	info2("with Fourier Domain preconditioner "); break;
			default: info2("Unknown preconditioner "); break;
		}
		break;
	case ALG_SVD:
		info2("SVD with threshold %g ", svdthres);
		break;
	default:
		info2("Invalid algorithm ");
	}
}

/**
   Selectively print out parameters for easy diagnose of possible mistakes.
*/
static void print_parms(const parms_t *parms){

	int i;
	const char *const phytype[]={
	"Skip",
	"matched filter",
	"CoG",
	"Maximum a posteriori tracing (MAP)",
	"correlation (peak first)",
	"correlation (sum first)",
	"Invalid"
	};
	const char *const closeloop[]={
	"open",
	"close"
	};

	info2("Aperture is %g m with sampling 1/%g m\n",
		parms->aper.d,1/parms->evl.dx);
	if(!parms->sim.noatm){
		real fgreen=calc_greenwood(parms->atm.r0z, parms->atm.nps, P(parms->atm.ws), P(parms->atm.wt));
		real theta0z=calc_aniso(parms->atm.r0z,parms->atm.nps,P(parms->atm.ht),P(parms->atm.wt));
	
		info2("Turbulence at %g degree zenith angle: r0=%gm, L0=%gm, %d layers.\n",
			parms->sim.za*180./M_PI,parms->atm.r0,P(parms->atm.L0,0),parms->atm.nps);
		info("    Greenwood freq is %.1fHz, anisoplanatic angle is %.2f as",
			fgreen,theta0z*RAD2AS);
		if(parms->ndm>1 && isfinite(theta0z)){
			real H1=parms->dm[0].ht;
			real H2=parms->dm[1].ht;
			real theta2z=calc_aniso2(parms->atm.r0z,parms->atm.nps,P(parms->atm.ht),P(parms->atm.wt),H1,H2);
			info(", generalized aa is %.2f as\n",theta2z*RAD2AS);
		} else{
			info("\n");
		}
		info("    Sampled %dx%d at 1/%gm. wind dir is%s randomized.\n",
			parms->atm.nx,parms->atm.ny,1./parms->atm.dx,
			(parms->atm.wdrand?"":" not"));
		if(parms->atm.nps>1&&theta0z*RAD2AS>4){
			warning("Atmosphere theta0 maybe wrong\n");
		}
		for(int ips=0; ips<parms->atm.nps; ips++){
			info("    layer %2d: height is %6.0f m, weight is %5.3f, wind speed is %4.1f m/s\n",
				ips,P(parms->atm.ht,ips),P(parms->atm.wt,ips),P(parms->atm.ws,ips));
		}
	}
	if(parms->recon.alg==RECON_MVR){
		info2("Reconstruction: r0=%gm L0=%gm. %d layers.%s\n",
			parms->atmr.r0,parms->atmr.L0,
			parms->atmr.nps,(parms->tomo.cone?" use cone coordinate.":""));

		for(int ips=0; ips<parms->atmr.nps; ips++){
			info("    layer %2d: height is %6.0f m, weight is %5.3f\n",
				ips,P(parms->atmr.ht,ips),P(parms->atmr.wt,ips));
		}
	}
	info2("There are %d powfs\n",parms->npowfs);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		info("  powfs%d: Order %2d, %sGS at %3.3g km. Sampling 1/%g m. Thres %g%%. ",
			ipowfs,parms->powfs[ipowfs].order,(parms->powfs[ipowfs].llt?"L":"N"),
			parms->powfs[ipowfs].hs/1000,1./parms->powfs[ipowfs].dx,parms->powfs[ipowfs].saat*100);
		int lrt=(parms->recon.split&&parms->tomo.splitlrt);
		if(parms->powfs[ipowfs].trs){
			info("%s%sis removed in %s side in tomography.",
				parms->powfs[ipowfs].trs?"T/T ":"",parms->powfs[ipowfs].frs?"Focus ":"",
				lrt?"both":"right hand");
		}
		info("\n");
		if(parms->powfs[ipowfs].usephy){
			if(parms->powfs[ipowfs].type==WFS_SH){
				info("    CCD image is %dx%d @ %gx%g mas, blur %g%% (sigma), %gHz, ",
					(parms->powfs[ipowfs].radpix?parms->powfs[ipowfs].radpix:parms->powfs[ipowfs].pixpsa),
					parms->powfs[ipowfs].pixpsa,
					parms->powfs[ipowfs].radpixtheta*RAD2MAS,parms->powfs[ipowfs].pixtheta*RAD2MAS,
					parms->powfs[ipowfs].pixblur*100,
					1./parms->sim.dt/parms->powfs[ipowfs].dtrat);
			} else{
				info("    PWFS, %gHz, ",1./parms->sim.dt/parms->powfs[ipowfs].dtrat);
			}
			info("wvl: [");
			for(int iwvl=0; iwvl<parms->powfs[ipowfs].nwvl; iwvl++){
				info(" %g",P(parms->powfs[ipowfs].wvl,iwvl));
			}
			info("]\n");
		}
		info("    %s in reconstruction. ",
			parms->powfs[ipowfs].gtype_recon==GTYPE_G?"Gtilt":"Ztilt");

		if(parms->powfs[ipowfs].step!=parms->powfs[ipowfs].phystep){
			info("Geomtric optics start at %d with %s ",
				parms->powfs[ipowfs].step,
				parms->powfs[ipowfs].gtype_sim==GTYPE_G?"gtilt":"ztilt");
		}
		if(parms->powfs[ipowfs].phystep>-1&&parms->powfs[ipowfs].phystep<parms->sim.end){
			info("Physical optics start at %d with '%s' ",
				parms->powfs[ipowfs].phystep, phytype[parms->powfs[ipowfs].phytype_sim]);
		}

		if(parms->powfs[ipowfs].noisy){
			info("(noisy)\n");
		} else{
			info("(noise free)\n");
		}
		if(parms->powfs[ipowfs].dither){
			info("    Delay locked loop starts at step %d and outputs every %d WFS frames.\n",
				parms->powfs[ipowfs].dither_pllskip,parms->powfs[ipowfs].dither_pllrat);
			info("    Pixel processing update starts at step %d and outputs every %d WFS frames.\n",
				parms->powfs[ipowfs].dither_ogskip,parms->powfs[ipowfs].dither_ograt);
		}
	}
	info2("There are %d wfs\n",parms->nwfs);
	for(i=0; i<parms->nwfs; i++){
		const int ipowfs=parms->wfs[i].powfs;
		const real rho=RSS(parms->wfs[i].thetax, parms->wfs[i].thetay)*RAD2AS;
		real th=rho==0?0:atan2(parms->wfs[i].thetay, parms->wfs[i].thetax)*180/M_PI;
		//if(th<0) th+=360;
		info("    wfs %d: powfs %d, at (%7.2f, %7.2f) (%5.1f, %4.0f°) arcsec, %3.0f km, siglev is %7.1f", i,
			parms->wfs[i].powfs,parms->wfs[i].thetax*RAD2AS,
			parms->wfs[i].thetay*RAD2AS, rho, th, 
			parms->wfs[i].hs*1e-3,parms->wfs[i].siglev*parms->powfs[ipowfs].dtrat);
		if((parms->wfs[i].siglev-parms->wfs[i].sigsim)>EPS){
			info(" (%g in simulation)",parms->wfs[i].sigsim);
		}
		info(" bkgrnd is %g",parms->powfs[ipowfs].bkgrnd);
		info("\n");
		if(fabs(parms->wfs[i].thetax)>1||fabs(parms->wfs[i].thetay)>1){
			warning("wfs %d thetax or thetay appears too large\n", i);
		}
	}
	info2("There are %d DMs\n",parms->ndm);
	for(i=0; i<parms->ndm; i++){
		info("  DM %d: at %g km, pitch %g m, offset %g, %g micron stroke, %g micron inter-actuator stroke.\n",
			i, parms->dm[i].ht/1000, parms->dm[i].dx/parms->dm[i].dratio,
			parms->dm[i].offset,
			fabs(P(parms->dm[i].stroke, 0))*1e6, fabs(parms->dm[i].iastroke)*1e6);
		if(parms->dm[i].iac){
			info("    Normalized cubic influence function with inter-actuator coupling of %g\n",
				parms->dm[i].iac);
		} else{
			info("    Bilinear influence function.\n");
		}
	}
	if(parms->recon.alg==RECON_MVR){
		if(!parms->sim.idealtomo){
			info2("Tomography is using ");
			print_alg(parms->tomo.bgs, parms->tomo.alg, parms->tomo.maxit, parms->tomo.precond, parms->tomo.svdthres);
			switch(parms->recon.split){
			case 0:
				info2(", integrated tomo.\n");break;
			case 1:
				info2(", ad hoc split tomo.\n"); break;
			case 2:
				info2(", minimum variance split tomo\n"); break;
			default:
				error(", Invalid\n");
			}
		}
		info2("DM Fitting is using ");
		
		print_alg(parms->fit.bgs,parms->fit.alg, parms->fit.maxit, parms->fit.precond, parms->fit.svdthres);
		info2("\nThere are %d DM fitting directions\n",parms->fit.nfit);
		for(i=0; i<parms->fit.nfit; i++){
			info("    Fit %d: weight is %5.3f, at (%7.2f, %7.2f) arcsec\n",
				i,P(parms->fit.wt,i),P(parms->fit.thetax,i)*RAD2AS,
				P(parms->fit.thetay,i)*RAD2AS);
			if(fabs(P(parms->fit.thetax,i))>1||fabs(P(parms->fit.thetay,i))>1){
				warning("fit %d thetax or thetay appears too large\n", i);
			}
		}
		if(parms->fit.nfit==1 && parms->ndm>1 && parms->evl.nevl>1){
			warning("There are multiple DMs and science evaluation directions but only one DM fitting direction.\n");
		}
	} else if(parms->recon.alg==RECON_LSR){
		info2("Least square reconstructor is using ");
		print_alg(parms->lsr.bgs, parms->lsr.alg, parms->lsr.maxit, 0, parms->lsr.svdthres);
		info2("\n");
	} else{
		error("parms->recon.alg=%d is not supported.\n",parms->recon.alg);
	}

	info2("There are %d evaluation directions at sampling 1/%g m.\n",
		parms->evl.nevl,1./parms->evl.dx);
	for(i=0; i<parms->evl.nevl; i++){
		info("    Evl %d: weight is %5.3f, at (%7.2f, %7.2f) arcsec\n",
			i,P(parms->evl.wt,i),P(parms->evl.thetax,i)*RAD2AS,
			P(parms->evl.thetay,i)*RAD2AS);
		if(fabs(P(parms->evl.thetax,i))>1||fabs(P(parms->evl.thetay,i))>1){
			warning("evl %d thetax or thetay appears too large\n", i);
		}
	}

	info2("Simulation start at step %d, end at step %d, "
		"with time step 1/%gs, %s loop.\n",
		parms->sim.start,parms->sim.end,1./parms->sim.dt, closeloop[parms->sim.closeloop]);
}

/**
   This routine calles other routines in this file to setup the parms parameter
   struct parms and check for possible errors. parms is kept constant after
   returned from setup_parms. */
parms_t *setup_parms(const char *mainconf,const char *extraconf,int override){
	char *config_path=find_config("maos");
	/*Setup PATH and result directory so that the config_path is in the back of path */
	if(!config_path||!exist(config_path)){
		error("Unable to find usable configuration file\n");
	}else{
		/*info2("Using config files found in %s\n", config_path); */
		addpath(config_path);
		addpath2(0, "%s/%s", config_path, "bin");
		addpath2(0, "%s/%s", config_path, "atm");
		addpath2(0, "%s/%s", config_path, "examples");
		free(config_path); config_path=NULL;
	}
	if(mainconf){
		open_config(mainconf, 0);/*user supplied main .conf file. */
	}
	if(extraconf && strlen(extraconf)){
		open_config(extraconf, 1);/*overriding .conf or lines*/
	}
	if(!mainconf && //sanity check for completeness.
		(!readcfg_peek("sim.skysim")||!readcfg_peek("save.evlopd")||
			!readcfg_peek("atm.r0z")||!readcfg_peek("recon.psdnseg"))){
		open_config("default.conf", 0);/*updates to this file is not tracked by git. */
	}
	parms_t *parms=mycalloc(1,parms_t);
	/*
		Conversion of input (e.g., in units) should be done in readcfg_* not in postproc_* as the values might be used early on.
	*/
	readcfg_dbg(parms);//2022-08-26: moved to front to use flags here.
	readcfg_sim(parms);
	readcfg_ncpa(parms);
	readcfg_aper(parms);
	readcfg_atm(parms);
	readcfg_dm(parms);
	readcfg_powfs(parms);//depends on readcfg_dm results
	readcfg_wfs(parms);
	readcfg_siglev(parms);
	readcfg_moao(parms);
	readcfg_atmr(parms);
	readcfg_tomo(parms);
	readcfg_fit(parms);
	readcfg_lsr(parms);
	readcfg_recon(parms);
	readcfg_evl(parms);
	readcfg_cn2(parms);
	readcfg_plot(parms);
	readcfg_gpu(parms);
	readcfg_save(parms);
	readcfg_distortion(parms);
	readcfg_load(parms);
	
	setup_parms_postproc_seeds(parms, override);
	/*
	  Output all the readed parms to a single file that can be used to reproduce
	  the same simulation.
	*/
	if(disable_save||parms->sim.nseed==0){
		close_config(NULL);
	} else{
		char fn[PATH_MAX];
		snprintf(fn, PATH_MAX, "maos_%s_%ld.conf", HOST, (long)getpid());
		close_config("%s", fn);
		mysymlink(fn, "maos_recent.conf");
	}
	/*
	  Postprocess the parameters for integrity. The ordering of the following
	  routines are critical.
	*/
	setup_parms_postproc_za(parms);
	setup_parms_postproc_sim(parms);
	setup_parms_postproc_wfs(parms);
	setup_parms_postproc_dtref(parms);
	setup_parms_postproc_dirs(parms);
	setup_parms_postproc_atm(parms);
	setup_parms_postproc_atm_size(parms);
	setup_parms_postproc_dm(parms);
	setup_parms_postproc_recon(parms);
	setup_parms_postproc_misc(parms);

	if(parms->sim.nseed>0){
		print_parms(parms);
		print_mem("After setup_parms");
	}
	return parms;
}
/**
   Additional setup_parms code to run when maos is running. It only contains GPU
   initialization code for the moment.
*/
void setup_parms_gpu(parms_t *parms,int *gpus,int ngpu){
#if USE_CUDA
	if(parms->sim.end==0){
		use_cuda=0;
	} else{
		use_cuda=1;
	}
	if(use_cuda){
		if(parms->sim.evlol){
			parms->gpu.tomo=0;
			parms->gpu.fit=0;
			parms->gpu.wfs=0;
			parms->gpu.lsr=0;
		}
		if(parms->sim.idealtomo){
			parms->gpu.wfs=0;
			parms->gpu.tomo=0;/*no need tomo.*/
			parms->fit.cachex=0;
		}
		if(parms->evl.rmax>1&&parms->gpu.evl){
			warning("evl.rmax>1 is not implemented in gpu. disable gpu.evl\n");
			parms->gpu.evl=0;
		}
		if(parms->recon.alg==RECON_MVR){/*MV*/
			parms->gpu.lsr=0;
			if(parms->gpu.tomo&&parms->dbg.tomo_hxw){
				warning("Disable gpu.tomo when dbg.tomo_hxw=1\n");
				parms->gpu.tomo=0;
			}
			if(parms->gpu.tomo&&parms->tomo.cxxalg!=0){
				parms->gpu.tomo=0;
				warning("\n\nGPU reconstruction is only available for tomo.cxxalg==0. Disable GPU Tomography.\n");
			}
			if(parms->gpu.tomo&&parms->tomo.alg>2){
				parms->gpu.tomo=0;
				warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Tomography.\n");
			}
			if(parms->gpu.fit&&parms->fit.alg>2){
				warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Fitting.\n");
				parms->gpu.fit=0;
			}
			if(parms->sim.idealtomo&&parms->gpu.fit){
				parms->gpu.fit=2;//in idealtomo FR is not assembled.
			}

		} else if(parms->recon.alg==RECON_LSR){
			parms->gpu.tomo=0;
			parms->gpu.fit=0;
		}
		parms->gpu.recon=(parms->gpu.tomo||parms->gpu.fit||parms->gpu.lsr);
		if(!parms->atm.frozenflow){
			warning("Atm is not frozen flow. Disabled gpu.evl and gpu.wfs.\n");
			parms->gpu.evl=0;
			parms->gpu.wfs=0;
		}
		if(parms->plot.run){
			//Do not accumulate PSF in GPU in order to plot individual PSF frame.
			parms->gpu.psf=0;
		}
	}
	/*use a max of one gpu if there is only 1 wfs.*/
	if(parms->nwfs==1&&ngpu==0) ngpu=1;
	if(use_cuda) use_cuda=gpu_init(parms,gpus,ngpu);
#else
	use_cuda=0; (void)gpus; (void)ngpu;
#endif
	//Other flags that depends on GPU enabling flags
	if(use_cuda){
		if(parms->recon.alg==RECON_MVR){/*MV*/
			if(parms->gpu.tomo||parms->gpu.fit==2){
				/*Tomography RHS in cuda always requrie full grid.*/
				parms->tomo.square=1;
			} 
			if(parms->evl.tomo && !parms->tomo.square){
				warning("evl.tomo without tomo.square is not implemented in gpu. disable gpu.evl\n");
				parms->gpu.evl=0;
			}
			if(parms->gpu.fit==1&&!parms->fit.assemble){
				info("\n\nGPU fitting=1 requries fit.assemble. Changed\n");
				parms->fit.assemble=1;
			}
			if(parms->gpu.fit==2&&!parms->fit.square){
				info("GPU fitting=2 requires fit.square=1. Changed\n");
				parms->fit.square=1;
			}
			if(parms->gpu.tomo&&parms->tomo.bgs){
				error("BGS in GPU is not implemented yet\n");
			}
			if(parms->gpu.fit!=2){//cache in gpu. only for grid based fit.
				parms->fit.cachedm=0;
				parms->fit.cachex=0;
			}
			if(parms->nmoao>0){
				if(parms->gpu.moao||parms->gpu.fit){
					if(!parms->fit.square){
						info("GPU moao=1 requires fit.square=1. Changed\n");
						parms->fit.square=1;
					}
					if(!parms->tomo.square){
						info("GPU moao=1 requires tomo.square=1. Changed\n");
						parms->tomo.square=1;
					}
				}
			}
		} else{
			parms->fit.square=0;
		}
		if((parms->gpu.evl||!parms->evl.nevl)&&(parms->gpu.wfs||!parms->nwfs)){
			parms->sim.cachedm=0; /*No need in CUDA. */
		}
	} else{
		memset(&(parms->gpu),0,sizeof(gpu_cfg_t));
	}
}
