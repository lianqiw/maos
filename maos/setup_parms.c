/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include "../lib/aos.h"
#include "parms.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
extern int use_cuda;
/*
  Don't include common.h or types.h, so that we don't have to recompile
  setup_parms.c when these files change.
 */
/**
   \file maos/setup_parms.c
   This file contains necessary routines to read parametes for
   WFS, DM and wavefront reconstruction.  */

void free_powfs_cfg(POWFS_CFG_T *powfscfg){
    dfree(powfscfg->wvl);
    if(powfscfg->wvlwts){
	dfree(powfscfg->wvlwts);
    }
    if(powfscfg->llt){
	free(powfscfg->llt->fnrange);
	free(powfscfg->llt->fnprof);
	free(powfscfg->llt->fnamp);
	free(powfscfg->llt->fnsurf);
	lfree(powfscfg->llt->i);
	dfree(powfscfg->llt->ox);
	dfree(powfscfg->llt->oy);
	dfree(powfscfg->llt->misreg);
	free(powfscfg->llt);
    }
    lfree(powfscfg->wfs);
    lfree(powfscfg->wfsind);
    free(powfscfg->fnllt);
    free(powfscfg->piinfile);
    free(powfscfg->sninfile);
    free(powfscfg->neareconfile);
    free(powfscfg->neasimfile);
    free(powfscfg->bkgrndfn);
}
void free_strarr(char **str, int n){
    if(str){
	for(int i=0; i<n; i++){
	    free(str[i]);
	}
	free(str);
    }
}
/**
   Create first order low pass filter coeffcient from cross over frequency and sampling rate.
*/
static double fc2lp(double fc, double dt){
    double lp=2*M_PI*fc*dt;
    if(lp>1){
	lp=1;
    }
    return lp;
}
/**
   Free the parms struct.
 */
void free_parms(PARMS_T *parms){
    dfree(parms->atm.ht);
    dfree(parms->atm.wt);
    dfree(parms->atm.ws);
    dfree(parms->atm.wddeg);
    dfree(parms->atmr.ht);
    dfree(parms->atmr.wt);
    lfree(parms->atmr.os);
    dfree(parms->atm.size);
    lfree(parms->atm.ipsr);
    lfree(parms->atm.overx);
    lfree(parms->atm.overy);
    lfree(parms->atmr.indps);
    dfree(parms->atm.r0logpsdt);
    dfree(parms->atm.r0logpsds);
    dfree(parms->evl.thetax);
    dfree(parms->evl.thetay);
    dfree(parms->evl.wvl);
    dfree(parms->evl.wt);
    dfree(parms->evl.hs);
    lfree(parms->evl.psf);
    lfree(parms->evl.psfr);
    lfree(parms->evl.psfgridsize);
    lfree(parms->evl.psfsize);
    lfree(parms->evl.pttr);
    lfree(parms->evl.psfngsr);

    dfree(parms->fit.thetax);
    dfree(parms->fit.thetay);
    dfree(parms->fit.wt);
    dfree(parms->fit.hs);

    dfree(parms->sim.apdm);
    dfree(parms->sim.epdm);
    dfree(parms->sim.aplo);
    dfree(parms->sim.eplo);
    dfree(parms->sim.apfsm);
    dfree(parms->sim.epfsm);
    lfree(parms->sim.seeds);
    free(parms->sim.wspsd);

    dfree(parms->sim.ncpa_thetax);
    dfree(parms->sim.ncpa_thetay);
    dfree(parms->sim.ncpa_wt);
    dfree(parms->sim.ncpa_hs);
    free(parms->sim.mvmhost);
    dfree(parms->cn2.pair);
    lfree(parms->save.gcov);
    for(int isurf=0; isurf<parms->nsurf; isurf++){
	free(parms->surf[isurf]);
    }
    free(parms->surf);
    for(int isurf=0; isurf<parms->ntsurf; isurf++){
	free(parms->tsurf[isurf]);
    }
    free(parms->tsurf);

    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	free_powfs_cfg(&parms->powfs[ipowfs]);
    }
    free(parms->powfs);
    if(parms->wfs!=parms->wfsr){
	free(parms->wfsr);
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	dfree(parms->wfs[iwfs].wvlwts);
	free(parms->wfs[iwfs].sabad);
    }
    free(parms->wfs);
    for(int idm=0; idm<parms->ndm; idm++){
	free(parms->dm[idm].hyst);
	free(parms->dm[idm].actstuck);
	free(parms->dm[idm].actfloat);
	free(parms->dm[idm].iastrokefn);
	cellfree(parms->dm[idm].iastrokescale);
	dfree(parms->dm[idm].stroke);
    }
    free(parms->dm);
    for(int imoao=0; imoao<parms->nmoao; imoao++){
	free(parms->moao[imoao].actstuck);
	free(parms->moao[imoao].actfloat);
    }
    free(parms->moao);
    free(parms->aper.fnamp);
    free(parms->aper.pupmask);
    lfree(parms->save.ints);
    lfree(parms->save.wfsopd);
    lfree(parms->save.grad);
    lfree(parms->save.gradgeom);
    free(parms->load.mvm);
    free(parms->load.mvmi);
    free(parms->load.mvmf);
    free(parms->load.ncpa);
    lfree(parms->fdlock);
    lfree(parms->hipowfs);
    lfree(parms->lopowfs);

    dfree(parms->misreg.pupil);
    free_strarr(parms->misreg.tel2wfs, parms->nwfs);
    free_strarr(parms->misreg.dm2wfs, parms->ndm*parms->nwfs);
    free_strarr(parms->misreg.dm2sci, parms->ndm*parms->evl.nevl);
    free_strarr(parms->recon.misreg_dm2wfs, parms->ndm*parms->nwfsr);
    free_strarr(parms->recon.misreg_dm2sci, parms->ndm*parms->fit.nfit);
    free_strarr(parms->recon.misreg_tel2wfs,parms->nwfsr);
    dfree(parms->dirs);
    lfree(parms->dbg.tomo_maxit);
    free(parms);
}
static inline int sum_intarr(int n, long *a){
    int sum=0;
    for(int i=0; i<n; i++){
	sum+=(a[i]!=0);
    }
    return sum;
}
static inline int sum_dblarr(int n, double *a){
    double sum=0;
    for(int i=0; i<n; i++){
	sum+=(a[i]!=0);
    }
    return sum;
}

#define MAX_STRLEN 80
#define READ_INT(A) parms->A = readcfg_int(#A) /*read a key with int value. */
#define READ_DBL(A) parms->A = readcfg_dbl(#A) /*read a key with double value */
#define READ_STR(A) parms->A = readcfg_str(#A) /*read a key with string value. */
#define READ_DMAT(A) parms->A= readcfg_dmat(#A) /*read a key with dmat. */

#define READ_POWFS(A,B)						\
    readcfg_##A##arr_n((void*)(&A##tmp), npowfs, "powfs."#B);	\
    for(i=0; i<npowfs; i++){					\
	parms->powfs[i].B = A##tmp[i];/*doesn't need ## in B*/	\
    }								
#define READ_POWFS_RELAX(A,B)						\
    readcfg_##A##arr_nmax((void*)(&A##tmp), npowfs, "powfs."#B);	\
    for(i=0; i<npowfs; i++){					\
	parms->powfs[i].B = A##tmp[i];/*doesn't need ## in B*/	\
    }								

/**
   Read wfs geometry. powfs stands for physical optics wfs,
   it is used to represent the types of WFS.
*/
static void readcfg_powfs(PARMS_T *parms){
    int     npowfs,i;
    parms->npowfs=npowfs=readcfg_peek_n("powfs.dsa");
    parms->powfs=calloc(parms->npowfs,sizeof(POWFS_CFG_T));
    int    *inttmp=NULL;
    double *dbltmp=NULL;
    char  **strtmp=NULL;
    READ_POWFS(dbl,dsa);
    READ_POWFS(int,nwvl);
    dmat* wvllist=readcfg_dmat("powfs.wvl");
    dmat* wvlwts=readcfg_dmat("powfs.wvlwts");
    dmat* siglev=readcfg_dmat("powfs.siglev");
    if(wvllist->nx != wvlwts->nx && wvlwts->nx != 0){
	error("powfs.wvl is not empty and does not match powfs.wvlwts\n");
    }
    if(siglev->nx!=0 && siglev->nx!=parms->npowfs){
	error("powfs.siglev is not empty and does not match npowfs");
    }
    int count=0;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int nwvl=parms->powfs[ipowfs].nwvl;
	parms->powfs[ipowfs].wvl=dnew(nwvl, 1);
	double wvlm=0;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    double wvl=wvllist->p[count+iwvl];
	    if(wvl>1e-3){
		wvl=wvl*1e-6;
	    }
	    if(wvl<=wvlm){
		error("Wavelength must in ascend order\n");
	    }
	    wvlm=wvl;
	    parms->powfs[ipowfs].wvl->p[iwvl]=wvl;
	}
	if(wvlwts->nx){
	    parms->powfs[ipowfs].wvlwts=dnew(nwvl, 1);
	    memcpy(parms->powfs[ipowfs].wvlwts->p, wvlwts->p+count,sizeof(double)*nwvl);
	    normalize_sum(parms->powfs[ipowfs].wvlwts->p, nwvl, 1);
	}
	if(siglev->nx){
	    parms->powfs[ipowfs].siglev=siglev->p[ipowfs];
	}else{
	    parms->powfs[ipowfs].siglev=-1;
	}
	count+=nwvl;
    }
    if(count!=wvllist->nx){
	error("powfs.wvl has wrong value\n");
    }
    dfree(wvllist);
    dfree(wvlwts);
    dfree(siglev);
    READ_POWFS_RELAX(str,saloc);
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
    READ_POWFS_RELAX(dbl,neaspeckle);
    READ_POWFS_RELAX(dbl,bkgrnd);
    READ_POWFS_RELAX(dbl,bkgrndc);
    READ_POWFS_RELAX(str,bkgrndfn);
    READ_POWFS_RELAX(str,bkgrndfnc);
    READ_POWFS_RELAX(dbl,pixblur);
    READ_POWFS_RELAX(dbl,radpixtheta);
    READ_POWFS_RELAX(int,radgx);
    READ_POWFS_RELAX(dbl,fieldstop);
    READ_POWFS_RELAX(dbl,pixoffx);
    READ_POWFS_RELAX(dbl,pixoffy);
    READ_POWFS_RELAX(int,phyusenea);
    READ_POWFS_RELAX(int,radpix);
    READ_POWFS_RELAX(int,radrot);
    READ_POWFS_RELAX(int,embfac);
    READ_POWFS_RELAX(int,ncomp);
    READ_POWFS_RELAX(int,psfout);
    READ_POWFS_RELAX(int,pistatout);
    READ_POWFS_RELAX(int,pistatstart);
    READ_POWFS_RELAX(int,pistatstc);
    READ_POWFS_RELAX(int,gtype_sim);
    READ_POWFS_RELAX(int,gtype_recon);
    READ_POWFS_RELAX(int,phytype);
    READ_POWFS_RELAX(int,phytypesim);
    READ_POWFS_RELAX(int,phytypesim2);
    READ_POWFS_RELAX(dbl,r0);
    READ_POWFS_RELAX(dbl,L0);
    READ_POWFS_RELAX(dbl,mtchcra);
    READ_POWFS_RELAX(int,mtchcpl);
    READ_POWFS_RELAX(int,mtchscl);
    READ_POWFS_RELAX(int,mtchadp);
    READ_POWFS_RELAX(dbl,cogthres);
    READ_POWFS_RELAX(dbl,cogoff);
    READ_POWFS_RELAX(int,ncpa_method);
    READ_POWFS_RELAX(int,i0scale);
    READ_POWFS_RELAX(dbl,sigscale);
    READ_POWFS_RELAX(int,moao);
    READ_POWFS_RELAX(int,dither);
    READ_POWFS_RELAX(dbl,gradscale);
    READ_POWFS_RELAX(dbl,dither_amp);
    READ_POWFS_RELAX(int,dither_npoint);
    READ_POWFS_RELAX(int,dither_pllskip);
    READ_POWFS_RELAX(int,dither_pllrat);
    READ_POWFS_RELAX(dbl,dither_gpll);
    READ_POWFS_RELAX(int,dither_ogskip);
    READ_POWFS_RELAX(int,dither_ograt);
    READ_POWFS_RELAX(dbl,dither_gog);

    READ_POWFS_RELAX(int, zoomdtrat);
    READ_POWFS_RELAX(int, zoomshare);
    READ_POWFS_RELAX(dbl, zoomgain);
    READ_POWFS_RELAX(int, zoomset);
    READ_POWFS(dbl,hs);
    READ_POWFS(dbl,nearecon);
    READ_POWFS(dbl,rne);
    READ_POWFS_RELAX(dbl,dx);
    READ_POWFS(dbl,pixtheta);
    READ_POWFS(str,fnllt);
    READ_POWFS(int,trs);
    READ_POWFS(int,dfrs);
    READ_POWFS(int,lo);
    READ_POWFS(int,pixpsa);
    READ_POWFS_RELAX(dbl,mtchcr);
    READ_POWFS_RELAX(int,mtchstc);
    READ_POWFS_RELAX(int,phystep);
    READ_POWFS_RELAX(int,noisy);
    READ_POWFS_RELAX(int,dtrat);
    READ_POWFS_RELAX(int,skip); 
    READ_POWFS_RELAX(int,type);
    READ_POWFS_RELAX(int,step);
    READ_POWFS_RELAX(dbl,modulate);
    READ_POWFS_RELAX(int,modulpos);
    READ_POWFS(int,nwfs);
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	POWFS_CFG_T *powfsi=&parms->powfs[ipowfs];
	if(!isfinite(powfsi->hs) && powfsi->fnllt){
	    warning2("powfs %d is at infinity, disable LLT\n", ipowfs);
	    free(powfsi->fnllt);
	    powfsi->fnllt=NULL;
	}
	if(powfsi->fnllt){
	    char prefix[60];
	    snprintf(prefix,60,"powfs%d_",ipowfs);
	    open_config(powfsi->fnllt,prefix,0);
	    powfsi->llt=calloc(1, sizeof(LLT_CFG_T));
	    powfsi->llt->d=readcfg_dbl("%sllt.d",prefix);
	    powfsi->llt->widthp=readcfg_dbl("%sllt.widthp",prefix);
	    powfsi->llt->ttrat=readcfg_dbl("%sllt.ttrat",prefix);
	    powfsi->llt->ttpsd=readcfg_str("%sllt.ttpsd",prefix);
	    powfsi->llt->fnrange=readcfg_str("%sllt.fnrange",prefix);
	    powfsi->llt->fnprof=readcfg_str("%sllt.fnprof",prefix);
	    powfsi->llt->fnamp=readcfg_str("%sllt.fnamp",prefix);
	    powfsi->llt->fnsurf=readcfg_str("%sllt.fnsurf",prefix);
	    powfsi->llt->ttfr=readcfg_int("%sllt.ttfr",prefix);
	    powfsi->llt->colprep=readcfg_int("%sllt.colprep",prefix); 
	    powfsi->llt->colsim=readcfg_int("%sllt.colsim",prefix);
	    powfsi->llt->colsimdtrat=readcfg_int("%sllt.colsimdtrat",prefix);
	    powfsi->llt->misreg=readcfg_dmat_n(2, "%sllt.misreg",prefix);
	    powfsi->llt->ox=readcfg_dmat("%sllt.ox",prefix);
	    powfsi->llt->oy=readcfg_dmat("%sllt.oy",prefix);
	    powfsi->llt->n=powfsi->llt->ox->nx;
	}else{/*there is no LLT. */
	    powfsi->llt=NULL;
	    if(isfinite(powfsi->hs)){
		warning2("powfs%d has finite hs at %g but no llt specified\n",
			ipowfs, powfsi->hs);
	    }
	    if(powfsi->radpix){
		warning2("powfs%d has no LLT, disable radial coordinate.\n", ipowfs);
		powfsi->radpix=0;
	    }
	}
	if(powfsi->radrot && !powfsi->radpix){
	    powfsi->radrot=0;
	    warning2("powfs%d does not have polar ccd. radrot should be zero. changed\n",ipowfs);
	}
	if(powfsi->llt && !powfsi->radpix && !powfsi->mtchcpl){
	    powfsi->mtchcpl=1;
	    warning2("powfs%d has llt, but no polar ccd or mtchrot=1, we need mtchcpl to be 1. changed\n",ipowfs);
	}
	double wvlmax=dmax(powfsi->wvl);
	if(powfsi->type==0){//shwfs
	    if(powfsi->phytypesim==-1){
		powfsi->phytypesim=powfsi->phytype;
	    }
	    if(powfsi->phytypesim2==-1){
		powfsi->phytypesim2=powfsi->phytypesim;
	    }
	    if(powfsi->mtchcra==-1){
		powfsi->mtchcra=powfsi->mtchcr;
	    }
	    int pixpsay=powfsi->pixpsa;
	    int pixpsax=powfsi->radpix;
	    if(!pixpsax) pixpsax=pixpsay;
	    if(pixpsax*pixpsay<4 && (powfsi->mtchcr>0 || powfsi->mtchcra>0)){
		powfsi->mtchcr=0;
		powfsi->mtchcra=0;
	    }
	    /*Senity check pixtheta*/
	    if(powfsi->pixtheta<0){
		powfsi->pixtheta=fabs(powfsi->pixtheta)*wvlmax/powfsi->dsa;
	    }else if(powfsi->pixtheta<1e-4){
		warning("powfs%d: pixtheta should be supplied in arcsec\n", ipowfs);
	    }else{
		powfsi->pixtheta/=206265.;/*convert form arcsec to radian. */
	    }
	}else if(powfsi->type==1){//pywfs only uses cog for the moment
	    powfsi->phytype=powfsi->phytypesim=powfsi->phytypesim2=2;//like quad cell cog
	    if(powfsi->phystep!=0){
		warning("PWFS must run in physical optics mode, changed.\n");
		powfsi->phystep=0;
	    }
	    powfsi->pixpsa=2;//always 2x2 pixels by definition.
	    //Input of modulate is in unit of wvl/D. Convert to radian
	    powfsi->modulate*=wvlmax/parms->aper.d;
	}
	if(powfsi->dither && powfsi->phystep!=0){
	    warning("Dither requrie physical optics mode from the beginning, changed.\n");
	    powfsi->phystep=0;
	}else if(powfsi->phystep>0){
	    /*round phystep to be multiple of dtrat. */
	    powfsi->phystep=((powfsi->phystep+powfsi->dtrat-1)/powfsi->dtrat)*powfsi->dtrat;
	}

	if(powfsi->fieldstop>0 && (powfsi->fieldstop>10 || powfsi->fieldstop<1e-4)){
	    error("powfs%d: fieldstop=%g. probably wrong unit. (arcsec)\n", ipowfs, powfsi->fieldstop);
	}
	powfsi->fieldstop/=206265.;

	if(powfsi->dither){
	    parms->dither=1;
	    if(powfsi->dither==1){//tip/tilt/arcsec->radian
		powfsi->dither_amp/=206265.;
	    }else{//zernike modes. micron-->meter
		powfsi->dither_amp/=1e6;
	    }
	    //Convert all rate in unit of WFS frame rate
	    //pllrat was already in WFS frame rate.
	    powfsi->dither_ograt*=powfsi->dither_pllrat;
	    
	    powfsi->dither_ogskip=powfsi->dither_ogskip*powfsi->dither_pllrat+powfsi->dither_pllskip;
	    //Convert all in simulation rate (sim.dt).
	    powfsi->dither_pllskip*=powfsi->dtrat;
	    powfsi->dither_ogskip*=powfsi->dtrat;
	    if(powfsi->dither_ograt<=0 || powfsi->dither_pllrat<=0){
		error("dither_ograt or _pllrat must be positive\n");
	    }
	}
    }/*ipowfs */
    free(inttmp);
    free(dbltmp);
    free(strtmp);
}
#define READ_WFS(A,B)							\
    readcfg_##A##arr_n((void*)(&A##tmp),nwfs,"wfs."#B);			\
    for(i=0; i<nwfs; i++){						\
	parms->wfs[i].B = A##tmp[i];					\
    }									
#define READ_WFS_RELAX(A,B)						\
    readcfg_##A##arr_nmax((void*)(&A##tmp),nwfs,"wfs."#B);		\
    for(i=0; i<nwfs; i++){						\
	parms->wfs[i].B = A##tmp[i];					\
    }									

/**
   Read in parameters of wfs, including GS direction, signal level, wvlwts, etc.
*/
static void readcfg_wfs(PARMS_T *parms){
    int i;
    int nwfs=parms->nwfs=readcfg_peek_n("wfs.thetax");
    parms->wfs=calloc(parms->nwfs,sizeof(struct WFS_CFG_T));
    double *dbltmp=NULL;
    int    *inttmp=NULL;
    char  **strtmp=NULL;
    READ_WFS(dbl,thetax);
    READ_WFS(dbl,thetay);
    for(i=0; i<parms->nwfs; i++){
	parms->wfs[i].thetax/=206265.;
	parms->wfs[i].thetay/=206265.;
    }
    READ_WFS_RELAX(dbl,hs);
    READ_WFS_RELAX(dbl,fitwt);
    READ_WFS_RELAX(str,sabad);

    /*link wfs with powfs*/
    int wfscount=0;
    int ipowfs=0;
    for(int kpowfs=0; kpowfs<parms->npowfs; kpowfs++, ipowfs++){
	if(parms->powfs[kpowfs].nwfs==0){//no stars.
	    free_powfs_cfg(&parms->powfs[kpowfs]);
	    parms->nwfs-=parms->powfs[kpowfs].nwfs;
	    ipowfs--;
	    continue;
	}else{
	    if(ipowfs<kpowfs){
		memcpy(parms->powfs+ipowfs, parms->powfs+kpowfs, sizeof(POWFS_CFG_T));
	    }
	}
	int mwfs=parms->powfs[ipowfs].nwfs;
	parms->powfs[ipowfs].wfs=lnew(mwfs, 1);
	parms->powfs[ipowfs].wfsind=lnew(parms->nwfs, 1);
	int count=0;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(iwfs>=wfscount && iwfs<wfscount+mwfs){
		parms->wfs[iwfs].powfs=ipowfs;
		parms->powfs[ipowfs].wfs->p[count]=iwfs;
		parms->powfs[ipowfs].wfsind->p[iwfs]=count;
		count++;
	    }else{
		parms->powfs[ipowfs].wfsind->p[iwfs]=-1;/*not belong */
	    }
	}
	wfscount+=mwfs;
    }
    parms->npowfs=ipowfs;
    if(parms->npowfs==0){
	warning("No wfs is found\n");
	if(!parms->sim.idealfit && !parms->sim.evlol){
	    error("Cannot proceed\n");
	}
    }else if(parms->nwfs!=wfscount){
	error("parms->nwfs=%d and sum(parms->powfs[*].nwfs)=%d mismatch\n", 
	      parms->nwfs, wfscount);
    }

    dmat *wvlwts=readcfg_dmat("wfs.wvlwts");
    dmat *siglev=readcfg_dmat("wfs.siglev");
    int powfs_siglev_override=readcfg_peek_override("powfs.siglev");
    int powfs_wvlwts_override=readcfg_peek_override("powfs.wvlwts");
    int count=0;
    if(siglev->nx!=0 && siglev->nx!=parms->nwfs){
	error("wfs.siglev can be either empty or %d\n",parms->nwfs);
    }
    for(i=0; i<parms->nwfs; i++){
	ipowfs=parms->wfs[i].powfs;
	int nwvl=parms->powfs[ipowfs].nwvl;
	parms->wfs[i].wvlwts=dnew(nwvl, 1);
	if(wvlwts->nx==0){
	    memcpy(parms->wfs[i].wvlwts->p,parms->powfs[ipowfs].wvlwts->p,sizeof(double)*nwvl);
	}else{
	    memcpy(parms->wfs[i].wvlwts->p,wvlwts->p+count,sizeof(double)*nwvl);
	    count+=nwvl;
	    if(parms->powfs[ipowfs].wvlwts && powfs_wvlwts_override){
		error("when both powfs.wvlwts and wfs.wvlwts are overriden "
		      "must set powfs.wvlwts=[]\n");
	    }
	}
	if(siglev->nx==0){
	    parms->wfs[i].siglev=parms->powfs[ipowfs].siglev;
	}else{
	    parms->wfs[i].siglev=siglev->p[i];
	    if(parms->powfs[ipowfs].siglev>0 && powfs_siglev_override){
		error("when both powfs.siglev and wfs.siglev are overriden "
		      "must set powfs.siglev=[]\n");
	    }
	}
    }
    if(count!=wvlwts->nx){
	error("Supplied %ld wvlwts but need %d for all wfs.\n",wvlwts->nx,count);
    }
    dfree(siglev);
    dfree(wvlwts);
    free(dbltmp);
    free(inttmp);
    free(strtmp);
}
#define READ_DM(A,B)					     \
    readcfg_##A##arr_n((void*)(&A##tmp),ndm,"dm."#B);	     \
    for(i=0; i<ndm; i++){				     \
	parms->dm[i].B = A##tmp[i];			     \
    }							     

#define READ_DM_RELAX(A,B) \
    readcfg_##A##arr_nmax((void*)(&A##tmp),ndm,"dm."#B);     \
    for(i=0; i<ndm; i++){				     \
	parms->dm[i].B = A##tmp[i];			     \
    }							     

/**
   Read in deformable mirror parameters.
*/
static void readcfg_dm(PARMS_T *parms){
    int ndm,i;
    ndm=parms->ndm=readcfg_peek_n("dm.ht");
    parms->dm=calloc(parms->ndm,sizeof(struct DM_CFG_T));
    int* inttmp=NULL;
    double *dbltmp=NULL;
    char **strtmp=NULL;
    dmat **dmattmp=NULL;
    READ_DM(dbl,ht);
    READ_DM(dbl,offset);
    READ_DM_RELAX(dbl,dx);
    READ_DM_RELAX(dbl,ar);
    for(int idm=0; idm<ndm; idm++){
	parms->dm[idm].order=parms->aper.d/parms->dm[idm].dx;
	parms->dm[idm].dy=parms->dm[idm].dx*parms->dm[idm].ar;
	if(parms->dm[idm].ar<=0){
	    error("ar must be positive\n");
	}
    }
    READ_DM_RELAX(dbl,guard);
    {
	char **tmp=0;
	int nstroke=readcfg_strarr(&tmp, "dm.stroke");
	for(int idm=0; idm<ndm; idm++){
	    if(nstroke==ndm){
		parms->dm[idm].stroke=readstr_dmat(tmp[idm]);
		free(tmp[idm]); tmp[idm]=0;
	    }else if(nstroke==1){
		parms->dm[idm].stroke=readstr_dmat(tmp[0]);
	    }else{
		error("dm.stroke is in wrong format\n");
	    }
	}
	free(tmp[0]);
	free(tmp);
    }
    READ_DM_RELAX(dbl,iastroke);
    READ_DM_RELAX(str,iastrokefn);
    READ_DM_RELAX(dbl,vmisreg);
    READ_DM_RELAX(dbl,histbin);
    READ_DM_RELAX(int,histn);
    READ_DM_RELAX(int,hist); 
    READ_DM_RELAX(dbl,iac);
    READ_DM_RELAX(str,hyst);
    READ_DM_RELAX(str,actfloat);
    READ_DM_RELAX(str,actstuck);
    free(strtmp);
    free(inttmp);
    free(dbltmp);
    free(dmattmp);
}
#define READ_MOAO(A,B)					      \
    readcfg_##A##arr_n((void*)(&A##tmp),nmoao,"moao."#B);     \
    for(i=0; i<nmoao; i++){				      \
	parms->moao[i].B = A##tmp[i];			      \
    }							      
#define READ_MOAO_RELAX(A,B)				      \
    readcfg_##A##arr_nmax((void*)(&A##tmp),nmoao,"moao."#B);  \
    for(i=0; i<nmoao; i++){				      \
	parms->moao[i].B = A##tmp[i];			      \
    }							      

/**
   Read in MOAO parameters.
*/
static void readcfg_moao(PARMS_T *parms){
    int nmoao=readcfg_peek_n("moao.dx");
    int i;
    parms->nmoao=nmoao;
    parms->moao=calloc(nmoao, sizeof(MOAO_CFG_T));
    int *inttmp=NULL;
    double *dbltmp=NULL;
    char **strtmp=NULL;
    READ_MOAO_RELAX(dbl,dx);
    for(int imoao=0; imoao<nmoao; imoao++){
	parms->moao[imoao].order=parms->aper.d/parms->moao[imoao].dx;
    }
    READ_MOAO_RELAX(dbl,iac);
    READ_MOAO_RELAX(dbl,gdm);
    READ_MOAO_RELAX(dbl,stroke);
    READ_MOAO_RELAX(dbl,ar);
    READ_MOAO_RELAX(int,actslave);
    READ_MOAO_RELAX(int,lrt_ptt);
    READ_MOAO_RELAX(dbl,guard);
    READ_MOAO_RELAX(str,actstuck);
    READ_MOAO_RELAX(str,actfloat);
    free(inttmp);
    free(dbltmp);
    free(strtmp);
}
/**
   Read in atmosphere parameters.
*/
static void readcfg_atm(PARMS_T *parms){
    READ_DBL(atm.r0z);
    READ_DBL(atm.L0);
    READ_DBL(atm.dx);
    READ_INT(atm.wdrand);
    READ_INT(atm.fractal);
    READ_INT(atm.frozenflow);
    READ_INT(atm.ninit);
    READ_INT(atm.share);
    READ_INT(atm.r0evolve);
    READ_DMAT(atm.r0logpsdt);
    READ_DMAT(atm.r0logpsds);
    READ_DMAT(atm.ht);
    //parms->atm.r0logpsdt=readcfg_dmat("atm.r0logpsdt");
    //parms->atm.r0logpsds=readcfg_dmat("atm.r0logpsds");
    parms->atm.size=readcfg_dmat_n(2, "atm.size");
    //parms->atm.ht=readcfg_dmat("atm.ht");
    parms->atm.nps=parms->atm.ht->nx;
    parms->atm.wt=readcfg_dmat_n(parms->atm.nps,"atm.wt");
    parms->atm.ws=readcfg_dmat_n(parms->atm.nps,"atm.ws");
    parms->atm.wddeg=readcfg_dmat_nmax(parms->atm.nps,"atm.wddeg");
    for(int ih=0; ih<parms->atm.nps; ih++){
	if(fabs(parms->atm.wddeg->p[ih])>1){
	    warning("wddeg is not zero. Disable wdrand\n");
	    parms->atm.wdrand=0;
	    break;
	}
    }
}
/**
   Read in atmosphere reconstruction parameters.
*/
static void readcfg_atmr(PARMS_T *parms){
    READ_DBL(atmr.r0z);
    if(parms->atmr.r0z<=0){
	parms->atmr.r0z=parms->atm.r0z;
    }
    READ_DBL(atmr.L0);
    if(parms->atmr.L0<=0){
	parms->atmr.L0=parms->atm.L0;
    }
    READ_DMAT(atmr.ht);
    //parms->atmr.ht=readcfg_dmat("atmr.ht");
    parms->atmr.nps=parms->atmr.ht->nx;
    parms->atmr.wt=readcfg_dmat_n(parms->atmr.nps, "atmr.wt");
    parms->atmr.os=readcfg_lmat_nmax(parms->atmr.nps, "atmr.os");
    READ_DBL(atmr.dx);
}

/**
   Read in aperture definitions.
*/
static void readcfg_aper(PARMS_T *parms){
    double *dtmp;
    /*aper.d may contain one for [d] or two numbers for [d din] */
    int nd=readcfg_dblarr(&dtmp, "aper.d");
    switch(nd){
    case 2:
	parms->aper.din=dtmp[1];/*don't break here. */
    case 1:
	parms->aper.d=dtmp[0];
	break;
    default:
	error("aper.d contains %d elements. But only 1 or 2 elements expected.\n", nd);
    }
    free(dtmp);

    if(parms->aper.d <= parms->aper.din){
	error("Inner dimeter: %g, Outer Diameter: %g. Illegal\n", parms->aper.din, parms->aper.d);
    }
    READ_DBL(aper.rotdeg);
    parms->aper.fnampuser=readcfg_peek_override("aper.fnamp");
    READ_STR(aper.fnamp);
    READ_STR(aper.pupmask);
}

/**
   Read in performance evaluation science point parameters.
*/
static void readcfg_evl(PARMS_T *parms){
    READ_DMAT(evl.thetax);
    //parms->evl.thetax=readcfg_dmat("evl.thetax");
    parms->evl.nevl=parms->evl.thetax->nx;
    parms->evl.thetay=readcfg_dmat_n(parms->evl.nevl, "evl.thetay");
    parms->evl.wt=readcfg_dmat_nmax(parms->evl.nevl, "evl.wt");
    parms->evl.hs=readcfg_dmat_nmax(parms->evl.nevl, "evl.hs");
    normalize_sum(parms->evl.wt->p, parms->evl.nevl, 1);
    parms->evl.psf=readcfg_lmat_nmax(parms->evl.nevl, "evl.psf");
    parms->evl.psfr=readcfg_lmat_nmax(parms->evl.nevl, "evl.psfr");
    READ_DMAT(evl.wvl);
    //parms->evl.wvl=readcfg_dmat("evl.wvl");
    parms->evl.nwvl=parms->evl.wvl->nx;
    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
	if(parms->evl.wvl->p[iwvl]>0.1){
	    warning("wvl should be supplied in unit of meter. scale %g by 1e-6\n",
		    parms->evl.wvl->p[iwvl]);
	    parms->evl.wvl->p[iwvl]*=1e-6;
	}
    }
    parms->evl.psfgridsize=readcfg_lmat_nmax(parms->evl.nwvl, "evl.psfgridsize");
    parms->evl.psfsize=readcfg_lmat_nmax(parms->evl.nwvl, "evl.psfsize");
    int ievl;
    double ramin=INFINITY;
    for(ievl=0; ievl<parms->evl.nevl; ievl++){
	/*First Convert theta to radian from arcsec. */
	parms->evl.thetax->p[ievl]/=206265.;
	parms->evl.thetay->p[ievl]/=206265.;
	double ra2=pow(parms->evl.thetax->p[ievl], 2)+pow(parms->evl.thetay->p[ievl], 2);
	if(ra2<ramin){
	    parms->evl.indoa=ievl;
	    ramin=ra2;
	}
    }
    READ_DBL(evl.dx);
    READ_INT(evl.rmax);
    READ_INT(evl.psfol);
    READ_INT(evl.psfisim);
    parms->evl.pttr=readcfg_lmat_nmax(parms->evl.nevl, "evl.pttr");
    parms->evl.psfngsr=readcfg_lmat_nmax(parms->evl.nevl, "evl.psfngsr");
    READ_INT(evl.psfmean); 
    READ_INT(evl.psfhist); 
    READ_INT(evl.cov);/*Science OPD covariance. */
    READ_INT(evl.tomo);
    READ_INT(evl.moao);
    /*it is never good to parallelize the evl ray tracing because it is already so fast */
    parms->evl.nmod=(parms->evl.rmax+1)*(parms->evl.rmax+2)/2;
}
/**
   Read in turbulence tomography parameters.
*/
static void readcfg_tomo(PARMS_T *parms){
    READ_INT(tomo.pos);
    READ_INT(tomo.cone);
    READ_INT(tomo.square);
    READ_INT(tomo.cxx);
    READ_DBL(tomo.cxxscale);
    READ_INT(tomo.guard);
    READ_DBL(tomo.tikcr);
    READ_DBL(tomo.svdthres);
    READ_INT(tomo.piston_cr);
    READ_INT(tomo.ahst_wt);
    READ_INT(tomo.ahst_idealngs);
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
    READ_DBL(tomo.cgthres);
    READ_INT(tomo.splitlrt);
}

/**
   Read in DM fit parameters. MOAO is specified elsewhere in readcfg_moao() */
static void readcfg_fit(PARMS_T *parms){
    READ_DMAT(fit.thetax);
    //parms->fit.thetax=readcfg_dmat("fit.thetax");
    parms->fit.nfit=parms->fit.thetax->nx;
    parms->fit.thetay=readcfg_dmat_n(parms->fit.nfit, "fit.thetay");
    parms->fit.wt=readcfg_dmat_nmax(parms->fit.nfit, "fit.wt");
    parms->fit.hs=readcfg_dmat_nmax(parms->fit.nfit, "fit.hs");
    double ramin=INFINITY;
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	parms->fit.thetax->p[ifit]/=206265.;
	parms->fit.thetay->p[ifit]/=206265.;
	double ra2=pow(parms->fit.thetax->p[ifit], 2)+pow(parms->fit.thetay->p[ifit], 2);
	if(ra2<ramin){
	    parms->fit.indoa=ifit;
	    ramin=ra2;
	}
    }

    READ_DBL(fit.tikcr);
    READ_DBL(fit.svdthres);
    READ_INT(fit.actslave);
    READ_DBL(fit.actthres);
    READ_INT(fit.actinterp);
    READ_INT(fit.lrt_piston);
    READ_INT(fit.lrt_tt);
    READ_INT(fit.alg);
    READ_INT(fit.bgs);
    READ_INT(fit.precond);
    READ_INT(fit.maxit);
    READ_INT(fit.square);
    READ_INT(fit.assemble);
    READ_INT(fit.pos);
    READ_INT(fit.cachedm);
    READ_INT(fit.cachex);
}
/**
   Read LSR parameters.
*/
static void readcfg_lsr(PARMS_T *parms){
    READ_DBL(lsr.tikcr);
    READ_DBL(lsr.svdthres);
    READ_STR(lsr.fnreg);
    READ_INT(lsr.alg);
    READ_INT(lsr.actslave);
    READ_INT(lsr.bgs);
    READ_INT(lsr.maxit);
    READ_DBL(lsr.actthres);
    READ_INT(lsr.actinterp);
}
/**
   Read general reconstruction parameters
 */
static void readcfg_recon(PARMS_T *parms){
    READ_INT(recon.alg);
    READ_INT(recon.glao);
    READ_INT(recon.split);
    READ_INT(recon.mvm);
    READ_INT(recon.modal);
    READ_INT(recon.nmod);
    READ_INT(recon.psol);
    READ_DBL(recon.poke);
    parms->nwfsr=parms->recon.glao?parms->npowfs:parms->nwfs;
    readcfg_strarr_nmax(&parms->recon.misreg_tel2wfs,parms->nwfsr, "recon.misreg_tel2wfs");  
    readcfg_strarr_nmax(&parms->recon.misreg_dm2wfs,parms->ndm*parms->nwfsr, "recon.misreg_dm2wfs");  
    readcfg_strarr_nmax(&parms->recon.misreg_dm2sci,parms->ndm*parms->fit.nfit, "recon.misreg_dm2sci");
    READ_INT(recon.psd);
    READ_INT(recon.psddtrat);
    READ_INT(recon.psddtrat_lo); 
    READ_INT(recon.psddtrat_twfs);
    READ_INT(recon.psdnseg);
    READ_STR(recon.fnsphpsd);
}
/**
   Read in simulation parameters
*/
static void readcfg_sim(PARMS_T *parms){
    READ_DBL(sim.fcfocus);
    READ_DBL(sim.fcttm);
    READ_INT(sim.mffocus);
    READ_DBL(sim.epfocus2tel);
    READ_INT(sim.focus2tel);
    READ_INT(sim.fsmideal);
    READ_DMAT(sim.apdm);
    READ_DMAT(sim.epdm);
    READ_DMAT(sim.aplo);
    READ_DMAT(sim.eplo);
    READ_DMAT(sim.apfsm);
    READ_DMAT(sim.epfsm);
    READ_DBL(sim.zetafsm);
    READ_DBL(sim.f0fsm);
    READ_DBL(sim.aptwfs);
    READ_DBL(sim.eptwfs);
    READ_INT(sim.aldm);
    READ_INT(sim.allo);
    READ_INT(sim.alfsm);
    /*We append a 0 so that we keep a time history of the integrator. */
    if(parms->sim.apdm->nx==1){
	dresize(parms->sim.apdm, 2, 1);
    }
    if(parms->sim.apdm->nx>2 || parms->sim.apdm->ny!=1){
	error("Invalid use of apdm\n");
    }
    if(parms->sim.aplo->nx==1){
	dresize(parms->sim.aplo, 2, 1);
    }
    if(parms->sim.aplo->nx>2 || parms->sim.aplo->ny!=1){
	error("Invalid use of aplo\n");
    }
    parms->sim.seeds=readcfg_lmat("sim.seeds");
    parms->sim.nseed=parms->sim.seeds->nx;
    READ_DBL(sim.dt);
    READ_INT(sim.dtrat_skip);
    READ_INT(sim.start);
    READ_INT(sim.end);
    READ_INT(sim.pause);
    READ_STR(sim.wspsd);
    READ_INT(sim.wsseq);
    READ_INT(sim.cachedm);
    READ_INT(sim.fuseint);
    READ_INT(sim.closeloop);
    READ_INT(sim.skysim);
    READ_DBL(sim.fov);parms->sim.fov/=206265;
    parms->sim.za = readcfg_dbl("sim.zadeg")*M_PI/180.;
    READ_INT(sim.evlol);
    READ_INT(sim.noatm);
    READ_INT(sim.idealfit);
    READ_INT(sim.idealtomo);
    READ_INT(sim.psfr);
    READ_INT(sim.ecnn);
    READ_INT(sim.wfsalias);
    READ_INT(sim.idealwfs);
    READ_INT(sim.idealevl);
    READ_INT(sim.ahstfocus);

    READ_STR(sim.mvmhost);
    READ_INT(sim.mvmport);
    READ_INT(sim.mvmsize);
    READ_INT(sim.mvmngpu);

    READ_INT(sim.ncpa_calib);
    READ_INT(sim.ncpa_ttr);
    parms->sim.ncpa_thetax=readcfg_dmat("sim.ncpa_thetax");
    parms->sim.ncpa_ndir=parms->sim.ncpa_thetax->nx;
    parms->sim.ncpa_thetay=readcfg_dmat_n(parms->sim.ncpa_ndir, "sim.ncpa_thetay");
    parms->sim.ncpa_wt=readcfg_dmat_n(parms->sim.ncpa_ndir, "sim.ncpa_wt");
    parms->sim.ncpa_hs=readcfg_dmat_nmax(parms->sim.ncpa_ndir, "sim.ncpa_hs");
    dscale(parms->sim.ncpa_thetax, 1./206265);
    dscale(parms->sim.ncpa_thetay, 1./206265);
    normalize_sum(parms->sim.ncpa_wt->p, parms->sim.ncpa_ndir, 1);
    READ_STR(sim.dmadd);
}
/**
   Read in parameters for Cn2 estimation.
 */
static void readcfg_cn2(PARMS_T *parms){
/*for Cn2 Estimation. */
    READ_DMAT(cn2.pair);
    //parms->cn2.pair = readcfg_dmat("cn2.pair");
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
    if(parms->cn2.pair && !parms->cn2.pair->nx){
	dfree(parms->cn2.pair);
    }
    if(!parms->cn2.pair){/*we are not doing cn2 estimation. */
	parms->cn2.tomo=0;
    }
}
/**
   Specify which variables to plot
*/
static void readcfg_plot(PARMS_T *parms){
  
    READ_INT(plot.setup);
    READ_INT(plot.atm);
    READ_INT(plot.run);
    READ_INT(plot.opdx);
    READ_INT(plot.all);
    if(parms->plot.all){
	parms->plot.setup=parms->plot.all;
	parms->plot.run=parms->plot.all;
    }
    if(parms->plot.setup || parms->plot.atm || parms->plot.run || parms->plot.opdx || parms->plot.all){
	draw_helper();
    }
}
/**
   Read in debugging parameters
*/
static void readcfg_dbg(PARMS_T *parms){
    READ_INT(dbg.wamethod);
    READ_INT(dbg.atm);
    READ_INT(dbg.mvstlimit);
    READ_INT(dbg.annular_W);
    parms->dbg.tomo_maxit=readcfg_lmat("dbg.tomo_maxit");
    READ_INT(dbg.tomo_hxw);
    READ_INT(dbg.ecovxx);
    READ_INT(dbg.useopdr);
    READ_INT(dbg.cmpgpu);
    READ_INT(dbg.pupmask);
    READ_INT(dbg.wfslinearity);
    READ_INT(dbg.nocgwarm);
    if(readcfg_peek("dbg.test")){
	READ_INT(dbg.test);
    }
    READ_INT(dbg.dmfullfov);
    READ_INT(dbg.tomo);
    READ_INT(dbg.fit);
    READ_INT(dbg.na_smooth);
    READ_INT(dbg.na_interp);
    READ_INT(dbg.ncpa_preload);
    READ_INT(dbg.ncpa_nouncorr);
    READ_INT(dbg.i0drift);
    READ_DBL(dbg.gradoff_scale);
}
/**
   Read in GPU options
*/
static void readcfg_gpu(PARMS_T *parms){
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
static void readcfg_save(PARMS_T *parms){
    READ_INT(save.extra);
    READ_INT(save.all);
    READ_INT(save.setup);
    READ_INT(save.recon);
    READ_INT(save.mvst);
    READ_INT(save.ncpa);
    READ_INT(save.atm);/*Save atmosphere */
    READ_INT(save.run);
    READ_INT(save.opdr);/*reconstructed OPD on XLOC */
    READ_INT(save.opdx);/*ATM propagated to XLOC */
    READ_INT(save.evlopd);/*Science OPD */
    READ_INT(save.dm);/*save DM commands */
    READ_INT(save.dither);
    parms->save.ints=readcfg_lmat_nmax(parms->nwfs, "save.ints");
    parms->save.wfsopd=readcfg_lmat_nmax(parms->nwfs, "save.wfsopd");
    parms->save.grad=readcfg_lmat_nmax(parms->nwfs, "save.grad");
    parms->save.gradgeom=readcfg_lmat_nmax(parms->nwfs, "save.gradgeom");
    if(disable_save){
	parms->save.extra=0;
    }
    if(parms->save.all){/*enables everything */
	if(disable_save){
	    error("please specify output directory\n");
	}
	warning("Enabling saving everything.\n");
	/*The following 3 are for setup. */
	if(!parms->save.setup) parms->save.setup=parms->save.all;
	if(!parms->save.recon) parms->save.recon=parms->save.all;
	parms->save.mvst=parms->save.all;
	/*The following are run time information that are not enabled by
	  save.run because they take a lot of space*/
	parms->save.opdr=parms->save.all;
	parms->save.opdx=parms->save.all;
	parms->save.evlopd=parms->save.all;
	parms->save.run=parms->save.all;/*see following */
	parms->save.ncpa=parms->save.all;
	parms->save.dither=parms->save.all;
    }

    if(parms->save.run){
	parms->save.dm=1;
	if(!parms->sim.idealfit){
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		parms->save.ints->p[iwfs]=1;
		parms->save.wfsopd->p[iwfs]=1;
		parms->save.grad->p[iwfs]=1;
		parms->save.gradgeom->p[iwfs]=1;
	    }
	}
    }
    parms->save.gcov=readcfg_lmat("save.gcov");
    parms->save.ngcov=parms->save.gcov->nx/2;
    READ_INT(save.gcovp);
    READ_INT(save.mvmi);
    READ_INT(save.mvmf);
    READ_INT(save.mvm);
}
static void readcfg_misreg(PARMS_T *parms){
    parms->misreg.pupil=readcfg_dmat_nmax(2, "misreg.pupil");
    readcfg_strarr_nmax(&parms->misreg.tel2wfs, parms->nwfs, "misreg.tel2wfs");
    readcfg_strarr_nmax(&parms->misreg.dm2wfs, parms->ndm*parms->nwfs, "misreg.dm2wfs");
    readcfg_strarr_nmax(&parms->misreg.dm2sci, parms->ndm*parms->evl.nevl, "misreg.dm2sci");
}
/**
   Specify which variables to load from saved files (Usually from LAOS
   simulations) */
static void readcfg_load(PARMS_T *parms){
  
    READ_STR(load.atm);
    READ_STR(load.locs);
    READ_STR(load.aloc);
    READ_STR(load.xloc);
    READ_STR(load.ploc);
    READ_STR(load.floc);
    READ_STR(load.cxx);
    READ_STR(load.HXF);
    READ_STR(load.HXW);
    READ_STR(load.HA);
    READ_STR(load.GP);
    READ_STR(load.GA);
    READ_STR(load.mvm);
    READ_STR(load.mvmi);
    READ_STR(load.mvmf);
    READ_STR(load.i0);
    READ_INT(load.mvst);
    READ_INT(load.GS0);
    READ_INT(load.tomo);
    READ_INT(load.fit);
    READ_INT(load.W);
    READ_STR(load.ncpa);
}
/**
   Process simulation parameters to find incompatibility.
*/
static void setup_parms_postproc_sim(PARMS_T *parms){
    if(parms->sim.skysim){
	if(disable_save){
	    error("sim.skysim requires saving. Please specify output folder\n");
	}
	/*if(parms->recon.alg!=0){
	    error("skysim need MVR");
	}*/
	parms->tomo.ahst_idealngs=1;
	if(parms->tomo.ahst_wt==1){//gradient weighting not available.
	    /*2013-1-30: ahst_wt=2 is not good. It resulted in higher NGS mode than ahst_wt=3*/
	    warning("in skycoverage presimulation, ahst_wt need to be 3. Changed\n");
	    parms->tomo.ahst_wt=3;
	}
	if(parms->ndm>0 && parms->recon.split!=1){
	    warning("Can only do skysim in split tomography mode 1. Changed\n");
	    parms->recon.split=1;
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].lo){
		parms->powfs[ipowfs].psfout=1;
	    }
	}
	parms->save.extra=1;
    }
    if(parms->dbg.tomo_maxit->nx){
	warning("dbg.tomo_maxit is set. Will run in open loop mode\n to repeat the simulations"
		" with different values of tomo.maxit.\n");
	parms->sim.closeloop=0;
	parms->atm.frozenflow=1;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    parms->atm.ws->p[ips]=0;/*set windspeed to zero. */
	}
	parms->sim.end=parms->dbg.tomo_maxit->nx;
    }
    if(parms->sim.idealfit || parms->sim.idealtomo){
	if(parms->sim.idealfit && parms->sim.idealtomo){
	    warning("idealfit takes precedence over idealtomo\n");
	}
	if(parms->recon.alg!=0){
	    warning("idealfit only works in recon.alg=0 mode. changed\n");
	    parms->recon.alg=0;
	}
	if(parms->recon.split){
	    warning("idealfit only works in integrated tomo mode. changed\n");
	    parms->recon.split=0;
	}
	if(parms->recon.mvm){
	    warning("idealfit cannot be used with recon.mvm. changed\n");
	    parms->recon.mvm=0;
	}
	if(parms->sim.wfsalias){
	    error("wfsalias and idealtomo/idealfit conflicts\n");
	}
	if(parms->sim.idealwfs){
	    error("idealwfs and idealtomo/idealfit conflicts\n");
	}
    }
    if(parms->sim.evlol){
	parms->recon.split=0;
    }
    if(parms->recon.glao && parms->ndm!=1){
	error("GLAO only works with 1 dm\n");
    }
    if(parms->recon.alg==1 && parms->recon.split==2){
	error("MVST does not work with least square reconstructor.\n");
    }
    if(parms->ndm>1 && parms->recon.modal){
	warning("Modal control is not supported for multiple DMs\n");
	parms->recon.modal=0;
    }
    if(parms->recon.alg==0 && parms->recon.modal){
	warning("Modal control is not supported yet with MV reconstructor. Disabled.\n");
	parms->recon.modal=0;
    }
    if(parms->lsr.actinterp==-1){
	if(parms->recon.alg==1 && parms->recon.modal){
	    //no need in modal lsr control
	    parms->lsr.actinterp=0;
	}else{
	    parms->lsr.actinterp=1;
	}
    }
    if(parms->recon.alg==1 && parms->lsr.alg==2){
	parms->recon.mvm=1;
    }
    if(parms->sim.wfsalias){
	if(parms->sim.idealwfs){
	    error("sim.wfsalias conflicts with sim.idealwfs. Do not enable both.\n");
	}
	if(parms->sim.idealevl){
	    error("sim.wfsalias conflicts with sim.idealevl. Do not enable both.\n");
	}
    }
    
    if(parms->sim.wfsalias || parms->sim.idealwfs || parms->sim.idealevl){
	parms->sim.dmproj=1;/*need dmproj */
    }
    if(parms->sim.ahstfocus){
	if(parms->recon.split!=1 || !parms->sim.mffocus){
	    parms->sim.ahstfocus=0;
	    warning("Disable ahstfocus\n");
	}
    }
    if(parms->sim.ncpa_calib && !parms->sim.ncpa_ndir){
	info2("Using evaluation directions as ncpa calibration directions if needed.\n");
	int ndir=parms->sim.ncpa_ndir=parms->evl.nevl;
	dfree(parms->sim.ncpa_thetax);
	dfree(parms->sim.ncpa_thetay);
	dfree(parms->sim.ncpa_wt);
	dfree(parms->sim.ncpa_hs);
	parms->sim.ncpa_thetax=dnew(ndir, 1);
 	parms->sim.ncpa_thetay=dnew(ndir, 1);
 	parms->sim.ncpa_wt=dnew(ndir, 1);
 	parms->sim.ncpa_hs=dnew(ndir, 1);
	dcp(&parms->sim.ncpa_thetax, parms->evl.thetax);
	dcp(&parms->sim.ncpa_thetay, parms->evl.thetay);
	dcp(&parms->sim.ncpa_wt, parms->evl.wt);
	dcp(&parms->sim.ncpa_hs, parms->evl.hs);
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
static void setup_parms_postproc_wfs(PARMS_T *parms){
    if(parms->sim.evlol){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    free_powfs_cfg(&parms->powfs[ipowfs]);
	}
	free(parms->powfs); parms->powfs=NULL;
	parms->npowfs=0;
	free(parms->wfs); parms->wfs=NULL;
	parms->nwfs=0;
    }
    //Check powfs.dsa
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].dsa<0){
	    parms->powfs[ipowfs].dsa*=-parms->aper.d;
	}else if(parms->powfs[ipowfs].dsa==0){
	    if(parms->ndm){
		parms->powfs[ipowfs].dsa=parms->dm[0].dx;
	    }else{
		parms->powfs[ipowfs].dsa=0.5;
	    }
	}
	parms->powfs[ipowfs].order=round(parms->aper.d/parms->powfs[ipowfs].dsa);
    }
    /*link wfs with powfs*/
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int mwfs=parms->powfs[ipowfs].nwfs;
	if(parms->powfs[ipowfs].llt){
	    parms->powfs[ipowfs].llt->i=lnew(mwfs, 1);/*default to zero. */
	    if(parms->powfs[ipowfs].llt->n>1){
		/*this is single llt for this powfs. */
		if(parms->powfs[ipowfs].llt->n!=mwfs)
		    error("# of llts should either be 1 or match nwfs for this powfs");
		for(int iwfs=0; iwfs<parms->powfs[ipowfs].llt->n; iwfs++){
		    parms->powfs[ipowfs].llt->i->p[iwfs]=iwfs;
		}
	    }
	}

	double wfs_hs=0;
	for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[indwfs];
	    if(parms->wfs[iwfs].hs<=EPS){
		parms->wfs[iwfs].hs=parms->powfs[ipowfs].hs;
	    }else{
		parms->sim.cachedm=0;
		wfs_hs+=parms->wfs[iwfs].hs;
	    }
	}
	wfs_hs/=parms->powfs[ipowfs].nwfs;
	if(parms->powfs[ipowfs].hs<=EPS){
	    if(wfs_hs>EPS){
		parms->powfs[ipowfs].hs=wfs_hs;
		warning2("powfs[%d].hs is set to %g\n", ipowfs, parms->powfs[ipowfs].hs);
	    }else{
		error("either wfs.wfs or powfs.hs has to be specified\n");
	    }
	}else if(wfs_hs>0 && fabs(wfs_hs-parms->powfs[ipowfs].hs)>10000){
	    warning2("powfs[%d].hs is %g, but wfs average hs is %g\n", ipowfs, parms->powfs[ipowfs].hs, wfs_hs);
	}
    }//for ipowfs
    
    //Match TWFS to LGS POWFS
    parms->itpowfs=-1;
    parms->ilgspowfs=-1;
    parms->nlgspowfs=0;
    for(int lgspowfs=0; lgspowfs<parms->npowfs; lgspowfs++){
	if(parms->powfs[lgspowfs].llt){
	    parms->nlgspowfs++;
	    if(parms->ilgspowfs==-1){
		parms->ilgspowfs=lgspowfs;
	    }else{
		warning("There are multiple LGS type. parms->ilgspowfs points to the first one\n");
	    }
	}
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].skip==2){//TWFS
	    parms->itpowfs=ipowfs;
	    int lgspowfs=parms->ilgspowfs;
	    if(lgspowfs!=-1){
		warning("powfs %d is TWFS for powfs %d\n", ipowfs, lgspowfs);
		if(parms->powfs[ipowfs].dtrat<1){
		    int mtchdtrat=parms->powfs[lgspowfs].dtrat*parms->powfs[lgspowfs].dither_ograt;
		    parms->powfs[ipowfs].dtrat=mtchdtrat;
		    warning("powfs %d dtrat is set to %d\n", ipowfs, mtchdtrat);
		    //Set TWFS integration start time to LGS matched filter acc step
		    parms->powfs[ipowfs].step=parms->powfs[lgspowfs].dither_ogskip;
		}else if(parms->powfs[lgspowfs].dither && parms->powfs[ipowfs].step<parms->powfs[lgspowfs].dither_pllskip){
		    //Set TWFS integration start time to pll start time to synchronize with matched filter.
		    parms->powfs[ipowfs].step=parms->powfs[lgspowfs].dither_pllskip;
		}
		//floor to multiple of dtrat.
		const int dtrat=parms->powfs[ipowfs].dtrat;
		parms->powfs[ipowfs].step=((parms->powfs[ipowfs].step+dtrat-1)/dtrat)*dtrat;
		warning("powfs %d step is set to %d, dtrat=%d\n", ipowfs, parms->powfs[ipowfs].step, dtrat);
	    }
	}
    }

    parms->hipowfs=lnew(parms->npowfs, 1);
    parms->lopowfs=lnew(parms->npowfs, 1);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs>0){
	    if(parms->powfs[ipowfs].lo){
		parms->lopowfs->p[parms->nlopowfs]=ipowfs;
		parms->nlopowfs++;
		parms->nlowfs+=parms->powfs[ipowfs].nwfs;
		if(parms->powfs[ipowfs].trs==1){
		    error("Low order wfs should not be tilt removed\n");
		}
		if(parms->powfs[ipowfs].gtype_sim==0 && parms->powfs[ipowfs].type==0){
		    warning("Low order POWFS %d is using gtilt in simulation. "
			    "This is not recommended\n", ipowfs);
		}
	    }else{
		parms->hipowfs->p[parms->nhipowfs]=ipowfs;
		parms->nhipowfs++;
		parms->nhiwfs+=parms->powfs[ipowfs].nwfs;
	    }
	    if(parms->powfs[ipowfs].trs){
		if(!parms->powfs[ipowfs].llt){
		    warning("WFS with tip/tilt removed should be LGS\n");
		}
		if(parms->powfs[ipowfs].lo){
		    warning("WFS with tip/tilt removed should be high order\n");
		}
		parms->ntrpowfs++;
	    }else{
		if(parms->powfs[ipowfs].llt){
		    warning("WFS with tip/tilt include should not be LGS\n");
		}
		parms->ntipowfs++;
	    }
	    if(parms->powfs[ipowfs].llt){
		if(!isfinite(parms->powfs[ipowfs].hs)){
		    warning("powfs with llt should have finite hs\n");
		}
	    }else{
		if(isfinite(parms->powfs[ipowfs].hs)){
		    warning("powfs without llt should infinite hs\n");
		}
	    }
	}
	if(fabs(parms->powfs[ipowfs].radpixtheta)<EPS){
	    parms->powfs[ipowfs].radpixtheta=parms->powfs[ipowfs].pixtheta;
	}else{
	    if(parms->powfs[ipowfs].radpixtheta>1e-4){
		parms->powfs[ipowfs].radpixtheta/=206265.;
	    }else if(parms->powfs[ipowfs].radpixtheta<0){
		error("powfs %d radpixtheta<0\n", ipowfs);
	    }
	}
	if (parms->powfs[ipowfs].phystep!=0 || parms->save.gradgeom){
	    parms->powfs[ipowfs].needGS0=1;
	}else{
	    parms->powfs[ipowfs].needGS0=0;
	}
	/*Do we ever do physical optics.*/
	if(parms->powfs[ipowfs].phystep>=0&&(parms->powfs[ipowfs].phystep<parms->sim.end||parms->sim.end==0)){
	    parms->powfs[ipowfs].usephy=1;
	    parms->nphypowfs++;
	}else{
	    parms->powfs[ipowfs].usephy=0;
	}
	if(!parms->powfs[ipowfs].usephy && parms->powfs[ipowfs].bkgrndfn){
	    warning("powfs %d: there is sky background, but is using geometric wfs. "
		    "background won't be effective.\n", ipowfs);
	} 
	if(parms->sim.ncpa_calib && (parms->nsurf || parms->ntsurf)){
	    if(parms->powfs[ipowfs].ncpa_method==2 && parms->powfs[ipowfs].mtchstc){
		warning("powfs %d: Disabling shifting i0 to center in the presence of NCPA.\n", ipowfs);
		parms->powfs[ipowfs].mtchstc=0;
	    }
	    if((!parms->powfs[ipowfs].usephy || parms->powfs[ipowfs].phytypesim2!=1)
	       && parms->powfs[ipowfs].ncpa_method==2){
		warning("powfs %d: ncpa_method changed from 2 to 1 in geometric wfs or CoG mode\n", ipowfs);
		parms->powfs[ipowfs].ncpa_method=1;
	    }
	    if(parms->tomo.ahst_idealngs && parms->powfs[ipowfs].ncpa_method==2 && parms->powfs[ipowfs].skip){
		warning("powfs %d: ncpa_method changed from 2 to 1 in idealngs mode\n", ipowfs);
		parms->powfs[ipowfs].ncpa_method=1;
	    }
	}
	{
	    /*Adjust dx if the subaperture does not contain integer, even number of points.*/
	    const double dsa=parms->powfs[ipowfs].dsa;
	    int nx = 2*(int)round(0.5*dsa/parms->powfs[ipowfs].dx);
	    if(nx<2) nx=2;
	    double dx=dsa/nx;/*adjust dx. */
	    if(fabs(parms->powfs[ipowfs].dx-dx)>EPS){
		warning("powfs %d: Adjusting dx from %g to %g. \n",
			ipowfs,parms->powfs[ipowfs].dx, dx);
	    }
	    parms->powfs[ipowfs].dx=dx;
	}
	if(!parms->sim.closeloop && parms->powfs[ipowfs].dtrat!=1){
	    warning("powfs %d: in open loop mode, only dtrat=1 is supported. Changed\n", ipowfs);
	    parms->powfs[ipowfs].dtrat=1;
	}
	if(parms->sim.wfsalias){
	    parms->powfs[ipowfs].noisy=0;
	    parms->powfs[ipowfs].phystep=-1;
	}
	if(parms->powfs[ipowfs].mtchscl==-1){
	    if(fabs(parms->powfs[ipowfs].sigscale-1)>EPS){
		parms->powfs[ipowfs].mtchscl=1;
	    }else{
		parms->powfs[ipowfs].mtchscl=0;
	    }
	}
    }
    parms->hipowfs->nx=parms->nhipowfs;
    parms->lopowfs->nx=parms->nlopowfs;
    if(!parms->nhipowfs){
	warning("There is no high order WFS!!!\n");
    }
    parms->sim.dtrat_hi=-1;
    parms->sim.dtrat_lo=-1;
    parms->step_lo=-1;
    parms->step_hi=-1;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].skip==2) continue;
	if(parms->powfs[ipowfs].type==1 && parms->powfs[ipowfs].llt){
	    error("Pyramid WFS is not available for LGS WFS\n");
	}
	if(parms->powfs[ipowfs].lo || (parms->nlopowfs==0 && !parms->powfs[ipowfs].trs)){//has t/t measurement.
	    if(parms->sim.dtrat_lo<0){
		parms->sim.dtrat_lo=parms->powfs[ipowfs].dtrat;
	    }else if(parms->sim.dtrat_lo!=parms->powfs[ipowfs].dtrat){
		error("We don't handle multiple framerate of the Tilt included WFS yet\n");
	    }
	    if(parms->step_lo<0){
		parms->step_lo=parms->powfs[ipowfs].step;
	    }else if(parms->step_lo!=parms->powfs[ipowfs].step){
		error("Different low order WFS has different enabling step\n");
	    }
	}
	if(!parms->powfs[ipowfs].lo){
	    if(!parms->powfs[ipowfs].skip){//participate in high order recon
		if(parms->sim.dtrat_hi<0){
		    parms->sim.dtrat_hi=parms->powfs[ipowfs].dtrat;
		}else if(parms->sim.dtrat_hi!=parms->powfs[ipowfs].dtrat){
		    error("We don't handle multiple framerate of the LO WFS yet\n");
		}
		if(parms->step_hi<0){
		    parms->step_hi=parms->powfs[ipowfs].step;
		}else if(parms->step_hi!=parms->powfs[ipowfs].step){
		    error("Different high order WFS has different enabling step\n");
		}
	    }
	}
    }
    parms->sim.dtlo=parms->sim.dtrat_lo*parms->sim.dt;
    parms->sim.dthi=parms->sim.dtrat_hi*parms->sim.dt;
    if(parms->sim.fcfocus<=0){
	parms->sim.fcfocus=1./parms->sim.dtlo/10;
	if(parms->sim.fcfocus<10){
	    parms->sim.fcfocus=10;
	}
    }
    parms->sim.lpfocushi=fc2lp(parms->sim.fcfocus, parms->sim.dthi);
    parms->sim.lpfocuslo=fc2lp(parms->sim.fcfocus, parms->sim.dtlo);
    parms->sim.lpttm=fc2lp(parms->sim.fcttm, parms->sim.dthi);
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
static void setup_parms_postproc_za(PARMS_T *parms){
    /*
      The input r0z is the r0 at zenith. Scale it if off zenith
    */
    parms->atm.r0=parms->atm.r0z*pow(cos(parms->sim.za),3./5.);
    parms->atmr.r0=parms->atmr.r0z*pow(cos(parms->sim.za),3./5.);

    if(fabs(parms->sim.za)>1.e-14){
	warning("Scaling turbulence height and LGS hs to zenith angle %gdeg\n",
		parms->sim.za*180./M_PI);
	double cosz=cos(parms->sim.za);
	double secz=1./cosz;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    parms->atm.ht->p[ips] *= secz;/*scale atmospheric height */
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    parms->powfs[ipowfs].hs *= secz;/*scale GS height. */
	    for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[indwfs];
		parms->wfs[iwfs].hs *= secz;
		
		if(isfinite(parms->wfs[iwfs].hs)){
		    double siglev=parms->wfs[iwfs].siglev;
		    parms->wfs[iwfs].siglev=siglev*cosz;/*scale signal level. */
		    info2("iwfs%d: siglev scaled from %g to %g\n", iwfs,siglev,parms->wfs[iwfs].siglev);
		    warning2("Need to update to account for the transmittance\n");
		}
	    }
	}
	warning("Scaling reconstruction height to zenith angle %gdeg\n",parms->sim.za*180./M_PI);
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    parms->atmr.ht->p[ips] *= secz;/*scale reconstructed atmospheric height. */
	}
	parms->cn2.hmax*=secz;
    }
}
/**
   The siglev is always specified in 800 Hz. If sim.dt is not 1/800, rescale the siglev.
*/
static void setup_parms_postproc_siglev(PARMS_T *parms){
    double sigscale=parms->sim.dt*800;
    if(fabs(sigscale-1.)>EPS){
	info2("sim.dt is 1/%g, need to scale siglev.\n",1/parms->sim.dt);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    double siglev=parms->wfs[iwfs].siglev;
	    parms->wfs[iwfs].siglev=siglev*sigscale;
	    info2("wfs%d: siglev scaled from %g to %g\n", iwfs,siglev,parms->wfs[iwfs].siglev);
	} 
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    double bkgrnd=parms->powfs[ipowfs].bkgrnd;
	    if(bkgrnd>0){
		parms->powfs[ipowfs].bkgrnd=bkgrnd*sigscale;
		info2("powfs%d: bkgrnd scaled from %g to %g\n", 
		      ipowfs,bkgrnd,parms->powfs[ipowfs].bkgrnd);
	    }
	}
    }
    
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	sigscale=parms->powfs[ipowfs].sigscale;
	if(fabs(sigscale-1)>0.5){
	    warning("powfs[%d].sigscale=%g\n",ipowfs,sigscale);
	}
	double siglev=parms->wfs[iwfs].siglev;
	parms->wfs[iwfs].siglevsim=siglev*sigscale;
	if(fabs(sigscale-1)>1.e-12){
	    warning("wfs%d: siglev in simulation scaled from %g to %g\n", 
		    iwfs,siglev,parms->wfs[iwfs].siglevsim);
	}
    }
}
/**
   postproc atmosphere parameters.
   1) drop weak layers.
   2) find ground layer
*/
static void setup_parms_postproc_atm(PARMS_T *parms){
    /*
      Drop weak turbulence layers in simulation.
    */
    int jps=0;
    for(int ips=0; ips<parms->atm.nps; ips++){
	if(parms->atm.wt->p[ips]>1.e-3){
	    if(ips!=jps){
		parms->atm.ht->p[jps]=parms->atm.ht->p[ips];
		parms->atm.wt->p[jps]=parms->atm.wt->p[ips];
		parms->atm.ws->p[jps]=parms->atm.ws->p[ips];
		parms->atm.wddeg->p[jps]=parms->atm.wddeg->p[ips];
	    }
	    jps++;
	}else{
	    warning("Layer %d has very small weight of %g, will drop it.\n", 
		    ips, parms->atm.wt->p[ips]);
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
    }
    if(parms->sim.idealfit){/*If fit only, we using atm for atmr. */
	warning("Changing atmr.ht,wt to atm.ht,wt since we are doing fit only\n");
	int nps=parms->atm.nps;
	dresize(parms->atmr.ht, nps, 1);
	dresize(parms->atmr.wt, nps, 1);
	dcp(&parms->atmr.ht, parms->atm.ht);
	dcp(&parms->atmr.wt, parms->atm.wt);
	lresize(parms->atmr.os, nps, 1);
	for(int ips=parms->atmr.nps; ips<nps; ips++){
	    parms->atmr.os->p[ips]=parms->atmr.os->p[parms->atmr.nps-1];
	}
	parms->atmr.nps=nps;
    }
    if((parms->recon.glao || parms->nhiwfs<=1) && parms->recon.alg==0){
	/*GLAO or single high wfs mode. reconstruct only a single layer near the DM.*/
	warning2("In GLAO or single high wfs Mode, use 1 tomography grid near the ground dm.\n");
	dresize(parms->atmr.ht, 1, 1);
	dresize(parms->atmr.wt, 1, 1);
	parms->atmr.os->nx=1;
	parms->atmr.ht->p[0]=parms->dm[0].ht;
	parms->atmr.wt->p[0]=1;
	parms->atmr.nps=1;
    }
    normalize_sum(parms->atm.wt->p, parms->atm.nps, 1);
    normalize_sum(parms->atmr.wt->p, parms->atmr.nps, 1);
    normalize_sum(parms->fit.wt->p, parms->fit.nfit, 1);
    /*
      We don't drop weak turbulence layers in reconstruction. Instead, we make
      it as least parms->tomo.minwt in setup_recon_tomo_prep
    */
    if(!parms->recon.glao){
	/*Assign each turbulence layer to a corresponding reconstructon layer. Used
	  to compute opdx in a simple minded way.*/
	parms->atm.ipsr=lnew(parms->atm.nps, 1);
	for(int ips=0; ips<parms->atm.nps; ips++){
	    double dist=INFINITY;
	    int kpsr=-1;
	    double ht=parms->atm.ht->p[ips];
	    for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
		double htr=parms->atmr.ht->p[ipsr];
		double dist2=fabs(ht-htr);
		if(dist2<dist){
		    dist=dist2;
		    kpsr=ipsr;
		}
	    }
	    parms->atm.ipsr->p[ips]=kpsr;
	    /*info("atm layer %d is maped to atmr %d\n", ips,kpsr); */
	}
    
	/* Map reconstructed layers to input layers. for testing tomo.predict*/
	parms->atmr.indps=lnew(parms->atmr.nps, 1);
	for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
	    parms->atmr.indps->p[ipsr]=-1;
	    for(int ips=0; ips<parms->atm.nps; ips++){
		if(fabs(parms->atmr.ht->p[ipsr]-parms->atm.ht->p[ips])<1e-3){
		    if(parms->atmr.indps->p[ipsr]>-1){
			warning("One ipsr is mapped to multiple ips\n");
		    }
		    parms->atmr.indps->p[ipsr]=ips;
		}
	    }
	}
    }
    /*
      Find ground turbulence layer. The ray tracing can be shared between different directions.
    */
    parms->atm.iground=-1;
    parms->atm.hmax=-INFINITY;
    for(int ips=0; ips<parms->atm.nps; ips++){
	if(fabs(parms->atm.ht->p[ips])<1.e-10){
	    if(parms->atm.iground==-1){
		parms->atm.iground=ips;
	    }else{
		error("Multiple grounds atm. Please combine them together.\n");
	    }
	}
	if(parms->atm.ht->p[ips]<0){
	    warning("Layer %d height %g is below ground\n",ips,parms->atm.ht->p[ips]);
	}
	if(ips > 0 && fabs(parms->atm.ht->p[ips]-parms->atm.ht->p[ips-1])<20){
	    warning("Layer %d at %gm is too close to layer %d at %gm\n",
		    ips,parms->atm.ht->p[ips], ips-1, parms->atm.ht->p[ips-1]);
	}
	if(parms->atm.hmax<parms->atm.ht->p[ips]){
	    parms->atm.hmax=parms->atm.ht->p[ips];
	}
    }
    parms->atmr.hmax=-INFINITY;
    for(int ips=0; ips<parms->atmr.nps; ips++){
	if(parms->atmr.hmax<parms->atmr.ht->p[ips]){
	    parms->atmr.hmax=parms->atmr.ht->p[ips];
	}
    }
    if(parms->atm.iground==-1){
	warning("There is no ground layer\n");
    }
    parms->atm.frozenflow = (parms->atm.frozenflow || parms->sim.closeloop);
    if(!parms->atm.frozenflow || parms->dbg.atm>0){
	parms->atm.r0evolve=0;/*disable r0 evolution*/
    }
    if(parms->dbg.atm==0){
	if(parms->atm.fractal){
	    parms->atm.fun=fractal_screen;
	}else{
	    parms->atm.fun=vonkarman_screen;
	}
    }else if(parms->dbg.atm==-1){
	info2("Generating Biharmonic Atmospheric Screen...");
	parms->atm.fun=biharmonic_screen;
    }else{
	parms->atm.fun=NULL;
    }
    if(!parms->atm.frozenflow && parms->sim.end>parms->sim.start+10){
	warning("Disable turbulence file based sharing in open loop nonfrozenflow simulation\n");
	parms->atm.share=0;
    }
}
static void setup_parms_postproc_dirs(PARMS_T *parms){
    //Collect all beam directions 
    const int ndir=parms->nwfs+parms->evl.nevl+parms->fit.nfit+(parms->sim.ncpa_calib?parms->sim.ncpa_ndir:0);
    parms->dirs=dnew(3, ndir);
    PDMAT(parms->dirs, pdir);
    int count=0;

    for(int i=0; i<parms->nwfs; i++){
	pdir[count][0]=parms->wfs[i].thetax;
	pdir[count][1]=parms->wfs[i].thetay;
	pdir[count][2]=parms->wfs[i].hs;
	count++;
    }
    
    for(int i=0; i<parms->evl.nevl; i++){
	pdir[count][0]=parms->evl.thetax->p[i];
	pdir[count][1]=parms->evl.thetay->p[i];
	pdir[count][2]=parms->evl.hs->p[i];
	count++;
    }
    for(int i=0; i<parms->fit.nfit; i++){
	pdir[count][0]=parms->fit.thetax->p[i];
	pdir[count][1]=parms->fit.thetay->p[i];
	pdir[count][2]=parms->fit.hs->p[i];
	count++;
    }
    if(parms->sim.ncpa_calib){
	for(int i=0; i<parms->sim.ncpa_ndir; i++){
	    pdir[count][0]=parms->sim.ncpa_thetax->p[i];
	    pdir[count][1]=parms->sim.ncpa_thetay->p[i];
	    pdir[count][2]=parms->sim.ncpa_hs->p[i];
	    count++;
	}
    }
    if(count<ndir){
	warning("count=%d, ndir=%d\n", count, ndir);
    }else if(count>ndir){
	error("count=%d, ndir=%d\n", count, ndir);
    }
    dresize(parms->dirs, 3, count);
}
/**
   compute minimum size of atm screen to cover all the beam path. same for
   all layers.  todo:may need to consider L0 Must be after
   setup_parms_postproc_za.
*/
static void setup_parms_postproc_atm_size(PARMS_T *parms){
    const int nps=parms->atm.nps;
    int Nmax=0;
    long nxout[nps],nyout[nps];
    for(int ips=0; ips<nps; ips++){
	create_metapupil(0,&nxout[ips],&nyout[ips],parms->dirs, parms->aper.d,parms->atm.ht->p[ips],
			 parms->atm.dx,parms->atm.dx,0.5,
			 parms->atm.dx*3,0,0,0,1);
	if(nxout[ips]>Nmax) Nmax=nxout[ips];
	if(nyout[ips]>Nmax) Nmax=nyout[ips];
    }
    /*Minimum screen size required. Used to transport atm to GPU. */
    parms->atm.nxn=Nmax;
    parms->atm.nyn=Nmax;
    parms->atm.nxm=nextpow2(Nmax);
    parms->atm.nym=parms->atm.nxm;
    if(fabs(parms->atm.size->p[0])<EPS ||fabs(parms->atm.size->p[1])<EPS){
	parms->atm.nx=parms->atm.nxm;
	parms->atm.ny=parms->atm.nym;
    }else{/*user specified.*/
	parms->atm.nx=2*(int)round(0.5*parms->atm.size->p[0]/parms->atm.dx);
	parms->atm.ny=2*(int)round(0.5*parms->atm.size->p[1]/parms->atm.dx);
	if(parms->atm.nx<parms->atm.nxm) parms->atm.nx=parms->atm.nxm;
	if(parms->atm.ny<parms->atm.nym) parms->atm.ny=parms->atm.nym;
    }
    if(parms->atm.fractal){/*must be square and 1+power of 2 */
	int nn=parms->atm.nx>parms->atm.ny?parms->atm.nx:parms->atm.ny;
	parms->atm.nx=1+nextpow2(nn);
	parms->atm.ny=parms->atm.nx;
    }
    /*record the size of the atmosphere. */
    parms->atm.size->p[0]=parms->atm.nx*parms->atm.dx;
    parms->atm.size->p[1]=parms->atm.ny*parms->atm.dx;
    if(parms->atm.L0 > parms->atm.size->p[0]){
	warning("Atmospheric size is smaller than outer scale!\n");
    }
    /*for screen evolving. */
    parms->atm.overx = lnew(parms->atm.nps, 1);
    parms->atm.overy = lnew(parms->atm.nps, 1);
    for(int ips=0; ips<parms->atm.nps; ips++){
	parms->atm.overx->p[ips] = nxout[ips];
	parms->atm.overy->p[ips] = nyout[ips];
    }
}
/*
  Find entry  that equals to val in array of length n. Append if not exist yet.
*/
static int arrind(double *arr, int *n, double val){
    for(long i=0; i<(*n); i++){
	if(fabs(arr[i]-val)<EPS){
	    return i;
	}
    }
    arr[*n]=val;
    (*n)++;
    return (*n)-1;
}

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
static void setup_parms_postproc_dm(PARMS_T *parms){
    /*disable cache for low order systems. */
    if(parms->sim.cachedm){
	if(parms->evl.nevl<2 && parms->nwfs<2){
	    parms->sim.cachedm=0;
	    warning("cachedm disabled for SCAO\n");
	}
	if(parms->dbg.cmpgpu){
	    parms->sim.cachedm=0;
	    warning("cachedm disabled when comparing CPU against GPU\n");
	}
    }
    for(int i=0; i<parms->ndm; i++){
	double ht=parms->dm[i].ht+parms->dm[i].vmisreg;
	if(fabs(ht)<1.e-10){
	    parms->dm[i].isground=1;
	    parms->idmground=i;
	}
	if(isfinite(parms->dm[i].stroke->p[0])){
	    double strokemicron=fabs(parms->dm[i].stroke->p[0])*1e6;
	    if(strokemicron<1 || strokemicron>50){
		warning("dm %d: stroke %g um is probably wrong\n",
			i, strokemicron);
	    }
	}
	if(isfinite(parms->dm[i].iastroke) && !parms->dm[i].iastrokefn){
	    double strokemicron=parms->dm[i].iastroke*1e6;
	    if(strokemicron<.1 || strokemicron>50){
		warning("dm %d: iastroke %g um is probably wrong\n",
			i, strokemicron);
	    }
	}
    }
}

/**
   Setting up the cone coordinate for MCAO LGS
   simulation. First find out the guide star conjugate. Only 1
   altitude is allowed.
*/
static void setup_parms_postproc_recon(PARMS_T *parms){
   if(parms->nmoao>0){//remove unused moao configurations
	int count=0;
	for(int imoao=0; imoao<parms->nmoao; imoao++){
	    int used=0;
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].moao>=parms->nmoao){
		    error("invalid powfs[%d].moao=%d", ipowfs, parms->powfs[ipowfs].moao);
		}
		if(parms->powfs[ipowfs].moao==imoao){
		    used=1;
		    break;
		}
	    }
	    if(parms->evl.moao>=parms->nmoao){
		error("invalid evl.moao=%d\n", parms->evl.moao);
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
    if((parms->recon.split) && parms->ndm==0){
	warning("Disable split tomography since there is no common DM\n");
	parms->recon.split=0;
    }
    if(parms->fit.square && parms->load.aloc){
	warning("load.aloc contradicts with fit.square. disable fit.square\n");
	parms->fit.square=0;
    }
    if(!parms->sim.closeloop){
	parms->recon.psol=0;//open loop does not need psol
    }else if(parms->recon.psol==-1){//automatic.
	if(parms->recon.alg==0){//MV perfers psol
	    parms->recon.psol=1;
	}else{//LSR perfers cl
	    parms->recon.psol=0;
	}
    }
    if(parms->recon.split){
	if(parms->nlopowfs==0){
	    if(parms->ntrpowfs>=parms->nhipowfs){
		warning("There is no WFS controlling tip/tilt.\n");
	    }else{
		warning("Split reconstruction is enabled when there is no low order WFS."
			" Will split the tip/tilt modes from high order wfs\n");
	    }
	}
	if(parms->sim.skysim && (parms->nhipowfs==0 || parms->nlopowfs==0)){
	    error("There is only high or low order WFS. can not do skycoverage presimulation\n");
	}
	if(!parms->nlopowfs && parms->tomo.ahst_wt==1){
	    warning("When there is no lowfs. Change ahst_wt from 1 to 3\n");
	    parms->tomo.ahst_wt=3;
	}
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].r0<=0){
	    parms->powfs[ipowfs].r0=parms->atm.r0;
	}
	if(parms->powfs[ipowfs].L0<=0){
	    parms->powfs[ipowfs].L0=parms->atm.L0;
	}
	if(parms->recon.split && parms->powfs[ipowfs].lo){
	    parms->powfs[ipowfs].skip=1;
	}
	if(parms->save.ngcov>0 || (parms->cn2.pair && !parms->powfs[ipowfs].lo && !parms->powfs[ipowfs].skip)){
	    /*focus tracking or cn2 estimation, or save gradient covariance.  */
	    parms->powfs[ipowfs].psol=1;
	}else if(parms->recon.psol){//PSOL reconstruction
	    /*low order wfs in ahst mode does not need psol. */
	    if((parms->recon.split==1 && parms->powfs[ipowfs].skip)){
		parms->powfs[ipowfs].psol=0;
	    }else{
		parms->powfs[ipowfs].psol=1;
	    }
	}
    }
  
    if(parms->recon.glao){
	parms->wfsr=calloc(parms->npowfs, sizeof(WFS_CFG_T));
	parms->nwfsr=parms->npowfs;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    parms->wfsr[ipowfs].thetax=0;
	    parms->wfsr[ipowfs].thetay=0;
	    parms->wfsr[ipowfs].hs=parms->powfs[ipowfs].hs;
	    parms->wfsr[ipowfs].powfs=ipowfs;
	    parms->powfs[ipowfs].nwfsr=1;
	}
	parms->fit.nfit=1;
	dresize(parms->fit.thetax, 1, 1);
	dresize(parms->fit.thetay, 1, 1);
	dresize(parms->fit.wt, 1, 1);
	dresize(parms->fit.hs, 1, 1);
	parms->fit.thetax->p[0]=0;
	parms->fit.thetay->p[0]=0;
	parms->fit.wt->p[0]=1;
    }else{/*Use same information as wfs. */
	parms->wfsr = parms->wfs;
	parms->nwfsr= parms->nwfs;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    parms->powfs[ipowfs].nwfsr=parms->powfs[ipowfs].nwfs;
	}
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
		warning2("Including wfs %d with wt %g in fitting\n", iwfs, parms->wfs[iwfs].fitwt);
		parms->fit.thetax->p[parms->fit.nfit]=parms->wfs[iwfs].thetax;
		parms->fit.thetay->p[parms->fit.nfit]=parms->wfs[iwfs].thetay;
		parms->fit.hs->p[parms->fit.nfit]=parms->wfs[iwfs].hs;
		parms->fit.wt->p[parms->fit.nfit]=parms->wfs[iwfs].fitwt;
		parms->fit.nfit++;
	    }
	    if(nfit2!=parms->fit.nfit){
		error("nfit2=%d, parms->fit.nfit=%d\n", nfit2, parms->fit.nfit);
	    }
	}
    }

    parms->recon.warm_restart = !parms->dbg.nocgwarm && parms->atm.frozenflow && !(parms->dbg.tomo_maxit && parms->dbg.tomo_maxit->nx>0);
    {
	double hs=NAN;
	/*find out the height to setup cone coordinate. */
	if(parms->tomo.cone){
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		/*skip wfs that does not participate in tomography*/
		if (parms->powfs[ipowfs].lo || parms->powfs[ipowfs].skip){
		    continue;
		}
		/*isinf and isfinite both return 0 on inf in FreeBSD 9.0.*/
		if(isnan(hs)){
		    hs=parms->powfs[ipowfs].hs;
		}else{
		    if(isfinite(hs) || isfinite(parms->powfs[ipowfs].hs)){
			if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
			    warning("Two high order POWFS with different hs found: %g and %g\n", 
				    hs, parms->powfs[ipowfs].hs);
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
	double mindsa=INFINITY;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if (parms->powfs[ipowfs].lo || parms->powfs[ipowfs].skip){
		continue;
	    }
	    if(parms->powfs[ipowfs].dsa<mindsa){
		mindsa=parms->powfs[ipowfs].dsa;
	    }
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    if(parms->dm[idm].dx<mindsa){
		mindsa=parms->dm[idm].dx;
	    }
	}
	for(int imoao=0; imoao<parms->nmoao; imoao++){
	    if(parms->moao[imoao].dx<mindsa){
		mindsa=parms->moao[imoao].dx;
	    }
	}
	parms->atmr.dx=mindsa;
    }


    if(parms->recon.split==1 && !parms->sim.closeloop){
	warning("ahst split tomography does not have good NGS correction in open loop\n");
    }
    if(parms->recon.split==2 && parms->sim.fuseint==1){
	warning("MVST Mode can only use separated integrator for the moment. Changed\n");
	parms->sim.fuseint=0;
    }
    if(!parms->recon.split && !parms->sim.fuseint){
	parms->sim.fuseint=1;/*integrated tomo. only 1 integrator. */
    }

    /*Tomography related*/
    if(parms->sim.closeloop && parms->evl.tomo){
	warning("Evaluating tomography performance is best done in open loop\n");
    }
    if(parms->recon.split && parms->evl.tomo){
	warning("Evaluating tomography performance is best done with integrated tomography.\n");
    }
    if(parms->sim.ecnn || parms->load.tomo || parms->tomo.alg!=1 || parms->tomo.bgs){
	parms->tomo.assemble=1;
    }
    if((parms->tomo.bgs || parms->tomo.alg != 1) && parms->tomo.cxx !=0){
	error("Only CG work with non L2 cxx.\n");
	parms->tomo.cxx=0;
    }
    if(parms->recon.alg==0 && parms->tomo.predict==1 && parms->tomo.alg!=1){
	error("Predictive tomography only works with CG. need to redo CBS/MVM after wind velocity is know.\n");
    }
    if(parms->tomo.alg==1 && parms->recon.alg==0){/*MVR with CG*/
	if(parms->tomo.precond>1){
	    error("Invalid preconditoner\n");
	}
    }else{
	parms->tomo.precond=0;/*No Preconditioner is available*/
    }
    /*Assign CG interations*/
    if(parms->tomo.alg==1 && parms->tomo.maxit==0){
	int factor;
	if(parms->recon.mvm){
	    factor=parms->load.mvmi?1:10;
	}else{
	    factor=parms->recon.warm_restart?1:10;
	}
	if(!parms->tomo.precond){
	    factor*=10;
	}
	parms->tomo.maxit=3*factor;
	if(parms->recon.mvm==1 && parms->tomo.splitlrt){
	    warning("recon.mvm==1 require tomo.splitlrt=0 due to stability issue. Changed\n");
	    parms->tomo.splitlrt=0;
	}
    }
    if(parms->tomo.bgs && parms->tomo.precond){
	error("Please implement the preconditioner for each block for BGS.\n");
    }

    /*DM Fitting related*/
    if(parms->fit.alg==1 && parms->fit.maxit==0){
	int factor;
	factor=parms->recon.warm_restart?1:10;
	parms->fit.maxit=4*factor;
    }
    /*Fitting tip/tilt constraint is only intended for multi DM*/
    if(parms->ndm<2 && parms->fit.lrt_tt){
	parms->fit.lrt_tt=0;
    }
    if(parms->ndm>2 && parms->fit.lrt_tt==2){
	warning("When there are more than 2 DMs, lrt_tt has to be 1 instead of 2. changed\n");
	parms->fit.lrt_tt=1;
    }
    if(parms->fit.lrt_tt<0 || parms->fit.lrt_tt>2){
	error("parms->fit.lrt_tt=%d is invalid\n", parms->fit.lrt_tt);
    }
    if(parms->load.fit || parms->fit.alg!=1 || parms->fit.bgs){
	parms->fit.assemble=1;
    }
    if(parms->fit.bgs && parms->fit.precond){
	error("Please implement the preconditioner for each block for BGS.\n");
    }
    if(parms->lsr.alg==1 && parms->lsr.maxit==0){
	int factor;
	factor=parms->recon.warm_restart?1:10;
	parms->lsr.maxit=30*factor;
    }
 
    if(parms->sim.mffocus){
	if(!parms->sim.closeloop || parms->sim.idealfit){
	    warning("mffocus is set, but we are in open loop mode or doing fitting only. disable\n");
	    parms->sim.mffocus=0;
	}
	if(parms->sim.mffocus<0 || parms->sim.mffocus>2){
	    error("parms->sim.mffocus=%d is invalid\n", parms->sim.mffocus);
	}
    }
    if(!parms->recon.mvm){
	if(parms->tomo.alg!=1 && parms->load.mvmi){
	    free(parms->load.mvmi);
	    parms->load.mvmi=NULL;
	}
	if(parms->load.mvmf){
	    free(parms->load.mvmf);
	    parms->load.mvmf=NULL;
	}
    }
 
    for(int idm=0; idm<parms->ndm; idm++){
	if(isfinite(parms->dm[idm].stroke->p[0])){
	    parms->sim.dmclip=1;
	}
	if(isfinite(parms->dm[idm].iastroke) && parms->dm[idm].iastroke>0){
	    parms->sim.dmclipia=1;
	    if(parms->dm[idm].iastrokefn){
		parms->dm[idm].iastrokescale=dcellread("%s", parms->dm[idm].iastrokefn);
	    }
	}
    }
    if(parms->sim.psfr){
	int fnd=sum_intarr(parms->evl.nevl, parms->evl.psfr->p);
	if(fnd==0){
	    error("sim.psfr is specified, but evl.psfr are all zero\n");
	}else{
	    info2("Output PSF reconstruction telemetry for %d directions\n", fnd);
	}
	parms->evl.psfmean=1;/*Saves psfmean for verification. */
	/*required memory to hold memory. */
	long covmem=(long)round(pow(parms->aper.d/parms->evl.dx,4))*8*fnd;
	if(covmem>MAX(NMEM, LONG_MAX/2) && parms->evl.dx > parms->atmr.dx*0.25+EPS){/*4G or actual */
	    error("parms->evl.dx=%g is probably too large to save ecxx. Recommend parms->evl.dx=%g\n", parms->evl.dx, parms->atmr.dx*0.25);
	}
    }
    if(parms->fit.pos<=0) parms->fit.pos=parms->tomo.pos;
 
    if(parms->recon.misreg_tel2wfs){
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    if(parms->recon.misreg_tel2wfs[iwfs]){
		warning("Set dbg.tomo_hxw=1\n");
		parms->dbg.tomo_hxw=1;
		parms->gpu.tomo=0;
	    }
	}
    }
}


/**
   postproc misc parameters.
*/
static void setup_parms_postproc_misc(PARMS_T *parms, int override){
    if(!disable_save){
	/*Remove seeds that are already done. */
	char fn[80];
	int iseed=0; 
	int jseed=0;
	parms->fdlock=lnew(parms->sim.nseed, 1);
	for(iseed=0; iseed<parms->sim.nseed; iseed++){
	    snprintf(fn, 80, "Res_%ld.done",parms->sim.seeds->p[iseed]);
	    if(exist(fn) && !override){
		parms->fdlock->p[iseed]=-1;
		warning2("Skip seed %ld because %s exist.\n", parms->sim.seeds->p[iseed], fn);
	    }else{
		remove(fn);
	    	snprintf(fn, 80, "Res_%ld.lock",parms->sim.seeds->p[iseed]);
		parms->fdlock->p[iseed]=lock_file(fn, 0, 0);
		if(parms->fdlock->p[iseed]<0){
		    warning2("Skip seed %ld because it is already running.\n",
			     parms->sim.seeds->p[iseed]);
		}else{
		    cloexec(parms->fdlock->p[iseed]);
		    if(jseed!=iseed){
			parms->sim.seeds->p[jseed]=parms->sim.seeds->p[iseed];
			parms->fdlock->p[jseed]=parms->fdlock->p[iseed];
		    }
		    jseed++;
		}
	    }
	}
	if(jseed!=parms->sim.nseed){
	    parms->sim.nseed=jseed;
	}
	if(parms->sim.nseed<1){
	    scheduler_finish(0);
	    info2("There are no seed to run. Use -O to override. Exit\n");
	    quit();
	}
    }
    info2("There are %d valid simulation seeds: ",parms->sim.nseed);
    for(int i=0; i<parms->sim.nseed; i++){
	info2(" %ld", parms->sim.seeds->p[i]);
    }
    info2("\n");
    if(parms->sim.nseed>1 && parms->dither){
	warning("Some of the dither mode updates parameters still persist for different seeds.\n");
    }
    if(parms->save.ngcov>0 && parms->save.gcovp<10){
	warning("parms->save.gcovp=%d is too small. It may fill your disk!\n",
		parms->save.gcovp);
    }
    if(parms->save.gcovp>parms->sim.end){
	parms->save.gcovp=parms->sim.end;
    }

    if(parms->dbg.cmpgpu){
	warning("Make cpu code follows gpu implementations.\n");
	parms->tomo.square=1;
    }
  
    if(!parms->atm.frozenflow){
	warning2("psfisim is set from %d to %d in openloop mode\n", parms->evl.psfisim, parms->sim.start);
	parms->evl.psfisim=parms->sim.start;
    }
    if(parms->evl.psfisim<parms->sim.start){
	parms->evl.psfisim=parms->sim.start;
    }
    if(parms->sim.closeloop==0 || parms->evl.tomo){
	/*disable parallelizing the big loop. */
	extern int PARALLEL;
	PARALLEL=0;
    }
    if(parms->evl.psfmean || parms->evl.psfhist){
	int fnd=sum_intarr(parms->evl.nevl, parms->evl.psf->p);
	if(fnd==0){
	    error("Required to output PSF, but evl.psf are all zero\n");
	}else{
	    info2("Output PSF for %d directions\n", fnd);
	}
    }
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	parms->evl.npsf+=(parms->evl.psf->p[ievl]>0);
	if(!parms->recon.split){
	    parms->evl.psfngsr->p[ievl]=0;
	}
	if(isfinite(parms->evl.hs->p[ievl]) && parms->evl.psfngsr->p[ievl]){
	    parms->evl.psfngsr->p[ievl]=0;
	    if(parms->evl.psfmean || parms->evl.psfhist || parms->evl.cov){
		warning("evl %d: star is not at infinity. disable NGS mode removal for it\n", ievl);
	    }
	}
    }
}

    /**
   Selectively print out parameters for easy diagnose of possible mistakes.
*/
static void print_parms(const PARMS_T *parms){
    
    int i;
    const char *phytype[]={
	"Skip",
	"\033[0;32mmatched filter\033[0;0m",
	"\033[0;31mthresholded center of gravity\033[0;0m",
	"\033[0;31mMaximum A Priori Tracing (MAP)\033[0;0m"
    };
    const char* tomo_precond[]={
	"\033[0;32mNo\033[0;0m",
	"\033[0;32mFourer Domain\033[0;0m",
	"\033[0;32mInvalid\033[0;0m"
    };
    const char *closeloop[]={
	"\033[0;31mopen\033[0;0m",
	"\033[0;32mclose\033[0;0m"
    };

    info2("\033[0;32mAperture\033[0;0m is %g m with sampling 1/%g m\n",
	  parms->aper.d, 1/parms->evl.dx);
    double fgreen=calc_greenwood(parms->atm.r0z, parms->atm.nps, parms->atm.ws->p, parms->atm.wt->p);
    double theta0z=calc_aniso(parms->atm.r0z, parms->atm.nps, parms->atm.ht->p, parms->atm.wt->p);
    
    info2("\033[0;32mTurbulence at zenith:\033[0;0m\n"
	  "Fried parameter r0 is %gm, Outer scale is %gm Greenwood freq is %.1fHz\n"
	  "Anisoplanatic angle is %.2f\"",
	  parms->atm.r0, parms->atm.L0, fgreen, theta0z*206265);
    if(parms->ndm==2){
	double H1=parms->dm[0].ht;
	double H2=parms->dm[1].ht;
	double theta2z=calc_aniso2(parms->atm.r0z, parms->atm.nps, parms->atm.ht->p, parms->atm.wt->p, H1, H2);
	info2(", generalized is %.2f\"", theta2z*206265);
    }
    info2("\n");
    info2("There are %d layers, sampled %dx%d at 1/%gm. ZA is %g deg. wind dir is%s randomized.\n",
	  parms->atm.nps, parms->atm.nx, parms->atm.ny,  1./parms->atm.dx,  
	  parms->sim.za*180/M_PI, (parms->atm.wdrand?"":" not"));
    if(parms->atm.nps>1 && theta0z*206265>4){
	warning("Atmosphere theta0 maybe wrong\n");
    }
    for(int ips=0; ips<parms->atm.nps; ips++){
	info2("layer %d: ht= %6.0f m, wt= %5.3f, ws= %4.1f m/s\n",
	      ips,parms->atm.ht->p[ips],parms->atm.wt->p[ips],parms->atm.ws->p[ips]);
    }
    if(parms->recon.alg==0){
	info2("\033[0;32mReconstruction\033[0;0m: r0=%gm l0=%gm "
	      "ZA is %g deg. %d layers.%s\n", 
	      parms->atmr.r0, parms->atmr.L0,  
	      parms->sim.za*180/M_PI, 
	      parms->atmr.nps,(parms->tomo.cone?" use cone coordinate.":""));
  
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    info2("layer %d: ht= %6.0f m, wt= %5.3f\n",
		  ips,parms->atmr.ht->p[ips],parms->atmr.wt->p[ips]);
	}
    }
    info2("\033[0;32mThere are %d powfs\033[0;0m\n", parms->npowfs);
    for(i=0; i<parms->npowfs; i++){
	info2("powfs %d: Order %2d, %sGS at %3.3g km. Thres %g%%",
	      i,parms->powfs[i].order, (parms->powfs[i].llt?"L":"N"),
	      parms->powfs[i].hs/1000,parms->powfs[i].saat*100);
	int lrt=(parms->recon.split && parms->tomo.splitlrt);
	if(parms->powfs[i].trs){
	    info2("\033[0;32m Tip/tilt is removed in %s side in tomography.\033[0;0m", lrt?"both":"right hand");
	    if(!parms->powfs[i].llt){
		warning("\n\ntrs=1, but this powfs doesn't have LLT!\n\n");
	    }
	}
	if(parms->powfs[i].dfrs){
	    if(parms->powfs[i].nwfs<2){
		parms->powfs[i].dfrs=0;
	    }else{
		info2("\033[0;32m Diff focus is removed in %s side in tomography.\033[0;0m",
		      lrt?"both":"right hand");
	    }
	    if(!parms->powfs[i].llt){
		warning("\n\ndfrs=1, but this powfs doesn't have LLT!\n\n");
	    }
	}
	if(parms->powfs[i].pixblur>1.e-12){
	    info2("\033[0;32m Pixel is blurred by %g.\033[0;0m", parms->powfs[i].pixblur);
	}
	info2("\n");
	if(parms->powfs[i].type==0){
	    info2("    CCD image is %dx%d @ %gx%gmas, %gHz, ", 
		  (parms->powfs[i].radpix?parms->powfs[i].radpix:parms->powfs[i].pixpsa), 
		  parms->powfs[i].pixpsa, 
		  parms->powfs[i].radpixtheta*206265000,parms->powfs[i].pixtheta*206265000,
		  1./parms->sim.dt/parms->powfs[i].dtrat);
	}else{
	    info2("    PWFS, %gHz, ", 1./parms->sim.dt/parms->powfs[i].dtrat);
	}
	info2("wvl: [");
	for(int iwvl=0; iwvl<parms->powfs[i].nwvl; iwvl++){
	    info2(" %g",parms->powfs[i].wvl->p[iwvl]);
	}
	info2("]\n");
	info2("    %s in reconstruction. ", 
	      parms->powfs[i].gtype_recon==0?"Gtilt":"Ztilt");
	if(parms->powfs[i].phystep>-1){
	    info2("Physical optics start at %d with '%s' %s",
		  parms->powfs[i].phystep, 
		  phytype[parms->powfs[i].phytypesim],
		  parms->powfs[i].mtchscl?"scaled":"");
	}else{
	    info2("Geomtric optics uses %s ",
		  parms->powfs[i].gtype_sim==0?"gtilt":"ztilt");
	}
	
	if(parms->powfs[i].noisy){
	    info2("\033[0;32m(noisy)\033[0;0m\n");
	}else{
	    info2("\033[0;31m(noise free)\033[0;0m\n");
	}
	if(parms->powfs[i].dither){
	    info2("    Delay locked loop starts at step %d and outputs every %d WFS frames.\n",
		  parms->powfs[i].dither_pllskip, parms->powfs[i].dither_pllrat);
	    info2("    Pixel processing update starts at step %d and outputs every %d WFS frames.\n",
		  parms->powfs[i].dither_ogskip, parms->powfs[i].dither_ograt);
	}
    }
    info2("\033[0;32mThere are %d wfs\033[0;0m\n", parms->nwfs);
    for(i=0; i<parms->nwfs; i++){
	info2("wfs %d: type is %d, at (%7.2f, %7.2f) arcsec, %g km, siglev is %g",
	      i,parms->wfs[i].powfs,parms->wfs[i].thetax*206265,
	      parms->wfs[i].thetay*206265, parms->wfs[i].hs*1e-3, parms->wfs[i].siglev);
	if((parms->wfs[i].siglev-parms->wfs[i].siglevsim)>EPS){
	    info2(" (%g in simulation)", parms->wfs[i].siglevsim);
	}
	const int ipowfs=parms->wfs[i].powfs;
	info2(" bkgrnd is %g", parms->powfs[ipowfs].bkgrnd);
	info2("\n");
	if(fabs(parms->wfs[i].thetax)>1 || fabs(parms->wfs[i].thetay)>1){
	    error("wfs thetax or thetay is too large\n");
	}
    }
    info2("\033[0;32mThere are %d DMs\033[0;0m\n",parms->ndm);
    for(i=0; i<parms->ndm; i++){
	info2("DM %d: Order %g, at %4gkm, actuator pitch %gm, offset %3g, with %f micron stroke.\n",
	      i, parms->dm[i].order,
	      parms->dm[i].ht/1000, parms->dm[i].dx,
	      parms->dm[i].offset, 
	      fabs(parms->dm[i].stroke->p[0])*1e6);
	if(parms->dm[i].iac){
	    info2("     Normalized cubic influence function with inter-actuator coupling of %g\n",
		  parms->dm[i].iac);
	}else{
	    info2("     Bilinear influence function.\n");
	}
    }
    if(parms->recon.alg==0){
	info2("\033[0;32mTomography\033[0;0m is using ");
	if(parms->tomo.bgs){
	    info2("Block Gauss Seidel with ");
	}
	switch(parms->tomo.alg){
	case 0:
	    info2("Cholesky back solve ");
	    break;
	case 1:
	    info2("CG, with %s preconditioner, \033[0;32m%d\033[0;0m iterations, ",
		  tomo_precond[parms->tomo.precond], parms->tomo.maxit);
	    break;
	case 2:
	    info2("SVD direct solve ");
	    break;
	case 3:
	    info2("Block Gauss Seidel ");
	    break;
	default:
	    error("Invalid\n");
	}
	switch(parms->recon.split){
	case 0:
	    info2(" integrated tomo.\n");break;
	case 1:
	    info2(" ad hoc split tomo.\n"); break;
	case 2:
	    info2(" minimum variance split tomo\n"); break;
	default:
	    error(" Invalid\n");
	}
	info2("\033[0;32mDM Fitting\033[0;0m is using ");
	if(parms->fit.bgs){
	    info2("Block Gauss Seidel with ");
	}
	switch(parms->fit.alg){
	case 0:
	    info2("Cholesky back solve");
	    break;
	case 1:
	    info2("CG, with %s preconditioner, \033[0;32m%d\033[0;0m iterations, ",
		  tomo_precond[parms->fit.precond], parms->fit.maxit);
	    break;
	case 2:
	    info2("SVD direct solve ");
	    break;
	case 3:
	    info2("Block Gauss Seidel ");
	    break;
	default:
	    error("Invalid");
	}
    }else if(parms->recon.alg==1){
	info2("\033[0;32mLeast square reconstructor\033[0;0m is using ");
	if(parms->tomo.bgs){
	    info2("Block Gauss Seidel with ");
	}
	switch(parms->lsr.alg){
	case 0:
	    info2("Cholesky back solve ");
	    break;
	case 1:
	    info2("CG%d", parms->tomo.maxit);
	    break;
	case 2:
	    info2("SVD direct solve ");
	    break;
	default:
	    error("Invalid\n");
	}
    }else{
	error("parms->recon.alg=%d is illegal\n", parms->recon.alg);
    }
    info2("\n");
    info2("\033[0;32mSimulation\033[0;0m start at step %d, end at step %d, "
	  "with time step 1/%gs, %s loop \n", 
	  parms->sim.start, parms->sim.end, 1./parms->sim.dt, 
	  closeloop[parms->sim.closeloop]);
    info2("\033[0;32mThere are %d fit directions\033[0;0m\n", parms->fit.nfit);
    for(i=0; i<parms->fit.nfit; i++){
	info2("Fit %d: wt is %5.3f, at (%7.2f, %7.2f) arcsec\n",
	      i,parms->fit.wt->p[i],parms->fit.thetax->p[i]*206265, 
	      parms->fit.thetay->p[i]*206265);
	if(fabs(parms->fit.thetax->p[i])>1 || fabs(parms->fit.thetay->p[i])>1){
	    error("fit thetax or thetay is too large\n");
	}
    }
    info2("\033[0;32mThere are %d evaluation directions\033[0;0m\n", parms->evl.nevl);
    for(i=0; i<parms->evl.nevl; i++){
	info2("Eval %d: wt is %5.3f, at (%7.2f, %7.2f) arcsec\n",
	      i,parms->evl.wt->p[i],parms->evl.thetax->p[i]*206265, 
	      parms->evl.thetay->p[i]*206265);
	if(fabs(parms->evl.thetax->p[i])>1 || fabs(parms->evl.thetay->p[i])>1){
	    error("evl thetax or thetay is too large\n");
	}
    }
}

/**
   This routine calles other routines in this file to setup the parms parameter
   struct parms and check for possible errors. parms is kept constant after
   returned from setup_parms. */
PARMS_T * setup_parms(char *mainconf, char *extraconf, int override){
    if(!mainconf){
	mainconf="default.conf";
    }
    info2("Main config file is %s\n", mainconf);

    /*Setup PATH and result directory so that the config_path is in the back of path */
    char *config_path=find_config("maos");
    if(!config_path || !exist(config_path)){
	error("Unable to find usable configuration file\n");
    }
    /*info2("Using config files found in %s\n", config_path); */
    char *bin_path=stradd(config_path, "/bin", NULL);
    addpath(config_path);
    addpath(bin_path);
    free(bin_path);
    free(config_path);
    open_config(mainconf,NULL,0);/*main .conf file. */
    open_config(extraconf, NULL, 1);
    PARMS_T* parms=calloc(1, sizeof(PARMS_T));
    readcfg_sim(parms);
    readcfg_aper(parms);
    readcfg_atm(parms);
    readcfg_powfs(parms);
    readcfg_wfs(parms);
    readcfg_dm(parms);
    readcfg_moao(parms);
    readcfg_atmr(parms);
    readcfg_tomo(parms);
    readcfg_fit(parms);
    readcfg_lsr(parms);
    readcfg_recon(parms);
    readcfg_evl(parms);
    readcfg_cn2(parms);
    readcfg_plot(parms);
    readcfg_dbg(parms);
    readcfg_gpu(parms);
    readcfg_save(parms);
    readcfg_misreg(parms);
    readcfg_load(parms);
    parms->nsurf=readcfg_strarr(&parms->surf, "surf");
    parms->ntsurf=readcfg_strarr(&parms->tsurf,"tsurf");
    /*
      Output all the readed parms to a single file that can be used to reproduce
      the same simulation.
    */

    if(disable_save){
	close_config(NULL);
    }else{
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "maos_%s_%ld.conf", HOST, (long)getpid());
	close_config("%s", fn);
	remove("maos_recent.conf");
	mysymlink(fn, "maos_recent.conf");
    }
    /*
      Postprocess the parameters for integrity. The ordering of the following
      routines are critical.
    */
    if(disable_save){
	if(parms->save.setup || parms->save.all || parms->sim.skysim || parms->evl.psfmean || parms->evl.psfhist){
	    error("Please specify -o to enable saving to disk\n");
	}
    }
    setup_parms_postproc_sim(parms);
    setup_parms_postproc_wfs(parms);
    setup_parms_postproc_za(parms);
    setup_parms_postproc_siglev(parms);
    setup_parms_postproc_dirs(parms);
    setup_parms_postproc_atm(parms);
    setup_parms_postproc_atm_size(parms);
    setup_parms_postproc_dm(parms);
    setup_parms_postproc_recon(parms);
    setup_parms_postproc_misc(parms, override);
    print_parms(parms);
    return parms;
}
/**
   Additional setup_parms code to run when maos is running. It only contains GPU
   initialization code for the moment.
 */
void setup_parms_gpu(PARMS_T *parms, int *gpus, int ngpu){
#if USE_CUDA 
    if(parms->sim.end==0){
	use_cuda=0;
    }else{
	use_cuda=1;
    }
    if(use_cuda){
	if(parms->evl.tomo){
	    parms->gpu.evl=0;
	}
	if(parms->recon.alg==0){/*MV*/
	    parms->gpu.lsr=0;
	    if(parms->gpu.tomo && parms->tomo.cxx!=0){
		parms->gpu.tomo=0;
		warning("\n\nGPU reconstruction is only available for tomo.cxx==0. Disable GPU Tomography.\n");
	    }
	    if(parms->gpu.tomo && parms->tomo.alg >2){
		parms->gpu.tomo=0;
		warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Tomography.\n");
	    }
	    if(parms->gpu.fit && parms->fit.alg > 2){
		warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Fitting.\n");
		parms->gpu.fit=0;
	    }
	}else if(parms->recon.alg==1){
	    parms->gpu.tomo=0;
	    parms->gpu.fit=0;
	}
    }
    /*use a max of one gpu if there is only 1 wfs.*/
    if(parms->nwfs==1 && ngpu==0) ngpu=1;
    if(use_cuda) use_cuda=gpu_init(parms, gpus, ngpu);
#else
    use_cuda=0;
#endif
    if(use_cuda){
	if(parms->sim.evlol){
	    memset(&parms->gpu, 0, sizeof(GPU_CFG_T));
	}
	if(parms->sim.idealfit){
	    parms->gpu.tomo=0;/*no need tomo.*/
	    parms->fit.cachex=0;
	}
	if(parms->sim.idealtomo){
	    parms->gpu.tomo=0;
	}
	if(parms->recon.alg==0){/*MV*/
	    if(parms->sim.idealfit && parms->gpu.fit){
		parms->gpu.fit=2;//In idealfit, FR is not assembled.
	    }
	    if(parms->gpu.fit==1 && !parms->fit.assemble){
		warning("\n\nGPU fitting=1 requries fit.assemble. Changed\n");
		parms->fit.assemble=1;
	    }
	    if(parms->gpu.fit==2 && !parms->fit.square){
		warning("GPU fitting=2 requires fit.square=1. Changed\n");
		parms->fit.square=1;
	    }
	    if(parms->nmoao>0){
		if(parms->gpu.moao || parms->gpu.fit){
		    if(!parms->fit.square){
			info2("GPU moao=1 requires fit.square=1. Changed\n");
			parms->fit.square=1;
		    }
		    if(!parms->tomo.square){
			info2("GPU moao=1 requires tomo.square=1. Changed\n");
			parms->tomo.square=1;
		    }
		}
	    }
	}else{
	    parms->fit.square=0;
	}
	if(!parms->atm.frozenflow){
	    warning("Atm is not frozen flow. Disable gpu.evl and gpu.wfs\n");
	    parms->gpu.evl=0;
	    parms->gpu.wfs=0;
	}
	if(parms->gpu.evl && parms->gpu.wfs){
	    parms->sim.cachedm=0; /*No need in CUDA. */
	}
	if(parms->gpu.tomo || parms->gpu.fit==2){
	    /*Tomography RHS in cuda requries full grid.*/
	    parms->tomo.square=1;
	}
	if(parms->gpu.tomo && parms->tomo.bgs){
	    error("BGS in GPU is not implemented yet\n");
	}
	if(parms->gpu.fit!=2){
	    parms->fit.cachedm=0;
	    parms->fit.cachex=0;
	}
    }else{
	memset(&(parms->gpu), 0, sizeof(GPU_CFG_T));
    }
}
