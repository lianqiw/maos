/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include "../lib/aos.h"
#include "parms.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
extern int use_cuda;
/*
  Don't include maos.h or types.h, so that we don't have to recompile
  setup_parms.c when these files change.
 */
/**
   \file maos/setup_parms.c
   This file contains necessary routines to read parametes for
   WFS, DM and wavefront reconstruction.  */

void free_powfs_cfg(POWFS_CFG_T *powfscfg){
    free(powfscfg->wvl);
    if(powfscfg->wvlwts){
	free(powfscfg->wvlwts);
    }
    if(powfscfg->llt){
	free(powfscfg->llt->fnrange);
	free(powfscfg->llt->fnprof);
	free(powfscfg->llt->fnamp);
	free(powfscfg->llt->fnsurf);
	free(powfscfg->llt->i);
	free(powfscfg->llt->ox);
	free(powfscfg->llt->oy);
	free(powfscfg->llt->misreg);
	free(powfscfg->llt);
    }
    free(powfscfg->wfs);
    free(powfscfg->wfsind);
    free(powfscfg->scalegroup);
    free(powfscfg->fnllt);
    free(powfscfg->piinfile);
    free(powfscfg->sninfile);
    free(powfscfg->neareconfile);
    free(powfscfg->neasimfile);
    free(powfscfg->bkgrndfn);
    free(powfscfg->misreg);
    free(powfscfg->ncpa);
}
/**
   Free the parms struct.
 */
void free_parms(PARMS_T *parms){
    free(parms->atm.ht);
    free(parms->atm.wt);
    free(parms->atm.ws);
    free(parms->atm.wddeg);
    free(parms->atm.ipsr);
    free(parms->atm.overx);
    free(parms->atm.overy);
    free(parms->atmr.ht);
    free(parms->atmr.os);
    free(parms->atmr.wt);
    free(parms->atmr.indps);
    free(parms->atm.size);

    free(parms->evl.thetax);
    free(parms->evl.thetay);
    free(parms->evl.wvl);
    free(parms->evl.wt);
    free(parms->evl.hs);
    free(parms->evl.psf);
    free(parms->evl.psfr);
    free(parms->evl.psfgridsize);
    free(parms->evl.psfsize);
    free(parms->evl.misreg);
    free(parms->evl.psfpttr);
    free(parms->evl.psfngsr);

    free(parms->fit.thetax);
    free(parms->fit.thetay);
    free(parms->fit.wt);
    free(parms->fit.hs);

    dfree(parms->sim.apdm);
    dfree(parms->sim.epdm);
    dfree(parms->sim.aplo);
    dfree(parms->sim.eplo);
    dfree(parms->sim.apupt);
    dfree(parms->sim.epupt);
    free(parms->sim.seeds);
    free(parms->sim.gtypeII_lo);
    free(parms->sim.wspsd);

    free(parms->cn2.pair);
    free(parms->save.gcov);
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
	free(parms->wfs[iwfs].wvlwts);
    }
    free(parms->wfs);
    for(int idm=0; idm<parms->ndm; idm++){
	free(parms->dm[idm].dxcache);
	free(parms->dm[idm].misreg);
	free(parms->dm[idm].hyst);
	free(parms->dm[idm].actstuck);
	free(parms->dm[idm].actfloat);
    }
    free(parms->dm);
    for(int imoao=0; imoao<parms->nmoao; imoao++){
	free(parms->moao[imoao].actstuck);
	free(parms->moao[imoao].actfloat);
    }
    free(parms->moao);
    free(parms->evl.scalegroup);
    free(parms->aper.fnamp);
    free(parms->aper.pupmask);
    free(parms->aper.misreg);
    free(parms->save.ints);
    free(parms->save.wfsopd);
    free(parms->save.grad);
    free(parms->save.gradgeom);
    free(parms->fdlock);
    free(parms);
}
static inline int sum_intarr(int n, int *a){
    int sum=0;
    for(int i=0; i<n; i++){
	sum+=(a[i]!=0);
    }
    return sum;
}

#define MAX_STRLEN 80
#define READ_INT(A) parms->A = readcfg_int(#A) /*read a key with int value. */
#define READ_DBL(A) parms->A = readcfg_dbl(#A) /*read a key with double value */
#define READ_STR(A) parms->A = readcfg_str(#A) /*read a key with string value. */


#define READ_POWFS(A,B)						\
    readcfg_##A##arr_n((void*)(&A##tmp), npowfs, "powfs."#B);	\
    for(i=0; i<npowfs; i++){					\
	parms->powfs[i].B = A##tmp[i];/*doesn't need ## in B*/	\
    }								

/**
   Read wfs geometry. powfs stands for physical optics wfs,
   it is used to represent the types of WFS.
*/
static void readcfg_powfs(PARMS_T *parms){
    int     npowfs,i;
    parms->npowfs=npowfs=readcfg_peek_n("powfs.order");
    parms->powfs=calloc(parms->npowfs,sizeof(POWFS_CFG_T));
    int    *inttmp=NULL;
    double *dbltmp=NULL;
    char  **strtmp=NULL;
    READ_POWFS(int,order);
    READ_POWFS(int,nwvl);
    double *wvllist=NULL;
    int nwvllist=readcfg_dblarr(&wvllist, "powfs.wvl");
    double *wvlwts=NULL;
    int nwvlwts=readcfg_dblarr(&wvlwts, "powfs.wvlwts");
    double *siglev=NULL;
    int nsiglev=readcfg_dblarr(&siglev, "powfs.siglev");
    if(nwvllist != nwvlwts && nwvlwts != 0){
	error("powfs.wvlwts is not empty and does not match powfs.wvlwts\n");
    }
    if(nsiglev!=0 && nsiglev!=parms->npowfs){
	error("powfs.siglev is not empty and does not match npowfs");
    }
    int count=0;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int nwvl=parms->powfs[ipowfs].nwvl;
	parms->powfs[ipowfs].wvl=calloc(nwvl, sizeof(double));
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    double wvl=wvllist[count+iwvl];
	    if(wvl>1e-3){
		wvl=wvl*1e-6;
	    }
	    parms->powfs[ipowfs].wvl[iwvl]=wvl;
	}
	if(nwvlwts){
	    parms->powfs[ipowfs].wvlwts=calloc(nwvl, sizeof(double));
	    memcpy(parms->powfs[ipowfs].wvlwts, wvlwts+count,sizeof(double)*nwvl);
	    normalize_sum(parms->powfs[ipowfs].wvlwts, nwvl, 1);
	}
	if(nsiglev){
	    parms->powfs[ipowfs].siglev=siglev[ipowfs];
	}else{
	    parms->powfs[ipowfs].siglev=-1;
	}
	count+=nwvl;
    }
    if(count!=nwvllist){
	error("powfs.wvl has wrong value\n");
    }
    free(wvllist);
    if(nwvlwts){
	free(wvlwts);
    }
    if(nsiglev){
	free(siglev);
    }
    READ_POWFS(str,piinfile);
    READ_POWFS(str,sninfile);
    READ_POWFS(dbl,hs);
    READ_POWFS(dbl,saat);
    READ_POWFS(int,neaphy);
    READ_POWFS(str,neareconfile);
    READ_POWFS(str,neasimfile);
    READ_POWFS(dbl,nearecon);
    READ_POWFS(dbl,neasim);
    READ_POWFS(dbl,neaspeckle);
    READ_POWFS(dbl,bkgrnd);
    READ_POWFS(dbl,bkgrndc);
    READ_POWFS(str,bkgrndfn);
    READ_POWFS(str,bkgrndfnc);
    READ_POWFS(dbl,pixblur);
    READ_POWFS(dbl,rne);
    READ_POWFS(dbl,dx);
    READ_POWFS(dbl,pixtheta);
    READ_POWFS(dbl,radpixtheta);
    READ_POWFS(dbl,fieldstop);
    READ_POWFS(dbl,pixoffx);
    READ_POWFS(dbl,pixoffy);
    READ_POWFS(int,phyusenea);
    READ_POWFS(str,fnllt);
    READ_POWFS(int,trs);
    READ_POWFS(int,dfrs);
    READ_POWFS(int,lo);
    READ_POWFS(int,pixpsa);
    READ_POWFS(int,radpix);
    READ_POWFS(int,radrot);
    READ_POWFS(int,ncomp);
    READ_POWFS(int,embfac);
    READ_POWFS(int,psfout);
    READ_POWFS(int,pistatout);
    READ_POWFS(int,pistatstart);
    READ_POWFS(int,pistatstc);
    READ_POWFS(int,gtype_sim);
    READ_POWFS(int,gtype_recon);
    READ_POWFS(int,phytype);
    READ_POWFS(int,phytypesim);
    READ_POWFS(dbl,mtchcr);
    READ_POWFS(dbl,mtchcra);
    READ_POWFS(int,mtchcpl);
    READ_POWFS(int,mtchstc);
    READ_POWFS(int,mtchscl);
    READ_POWFS(int,mtchadp);
    READ_POWFS(dbl,cogthres);
    READ_POWFS(dbl,cogoff);
    READ_POWFS(int,phystep);
    READ_POWFS(int,noisy);
    READ_POWFS(str,misreg);
    READ_POWFS(str,ncpa);
    READ_POWFS(int,ncpa_method);
    READ_POWFS(int,dtrat);
    READ_POWFS(int,i0scale);
    READ_POWFS(dbl,sigscale);
    READ_POWFS(int,moao);

    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	if(!isfinite(parms->powfs[ipowfs].hs) && parms->powfs[ipowfs].fnllt){
	    warning2("powfs %d is at infinity, disable LLT\n", ipowfs);
	    free(parms->powfs[ipowfs].fnllt);
	    parms->powfs[ipowfs].fnllt=NULL;
	}
	if(parms->powfs[ipowfs].fnllt){
	    char prefix[60];
	    snprintf(prefix,60,"powfs%d_",ipowfs);
	    open_config(parms->powfs[ipowfs].fnllt,prefix,0);
	    parms->powfs[ipowfs].llt=calloc(1, sizeof(LLT_CFG_T));
	    parms->powfs[ipowfs].llt->d=readcfg_dbl("%sllt.d",prefix);
	    parms->powfs[ipowfs].llt->widthp=readcfg_dbl("%sllt.widthp",prefix);
	    parms->powfs[ipowfs].llt->fnrange=readcfg_str("%sllt.fnrange",prefix);
	    parms->powfs[ipowfs].llt->fnprof=readcfg_str("%sllt.fnprof",prefix);
	    parms->powfs[ipowfs].llt->fnamp=readcfg_str("%sllt.fnamp",prefix);
	    parms->powfs[ipowfs].llt->fnsurf=readcfg_str("%sllt.fnsurf",prefix);
	    parms->powfs[ipowfs].llt->smooth=readcfg_int("%sllt.smooth",prefix);
	    parms->powfs[ipowfs].llt->colprep=readcfg_int("%sllt.colprep",prefix);
	    parms->powfs[ipowfs].llt->colsim=readcfg_int("%sllt.colsim",prefix);
	    parms->powfs[ipowfs].llt->colsimdtrat=readcfg_int("%sllt.colsimdtrat",prefix);
	    readcfg_dblarr_n(&parms->powfs[ipowfs].llt->misreg,2,"%sllt.misreg",prefix);
	    parms->powfs[ipowfs].llt->n=readcfg_dblarr(&(parms->powfs[ipowfs].llt->ox),"%sllt.ox",prefix);
	    readcfg_dblarr_n(&(parms->powfs[ipowfs].llt->oy),parms->powfs[ipowfs].llt->n,"%sllt.oy",prefix);

	}else{/*there is no LLT. */
	    parms->powfs[ipowfs].llt=NULL;
	    if(isfinite(parms->powfs[ipowfs].hs)){
		warning2("powfs%d has finite hs at %g but no llt specified\n",
			ipowfs, parms->powfs[ipowfs].hs);
	    }
	    if(parms->powfs[ipowfs].radpix){
		warning2("powfs%d has no LLT, disable radial coordinate.\n", ipowfs);
		parms->powfs[ipowfs].radpix=0;
	    }
	}
	if(parms->powfs[ipowfs].radrot && !parms->powfs[ipowfs].radpix){
	    parms->powfs[ipowfs].radrot=0;
	    warning2("powfs%d does not have polar ccd. radrot should be zero. changed\n",ipowfs);
	}
	if(parms->powfs[ipowfs].llt && !parms->powfs[ipowfs].radpix && !parms->powfs[ipowfs].mtchcpl){
	    parms->powfs[ipowfs].mtchcpl=1;
	    warning2("powfs%d has llt, but no polar ccd or mtchrot=1, we need mtchcpl to be 1. changed\n",ipowfs);
	}
	if(parms->powfs[ipowfs].phytypesim==-1){
	    parms->powfs[ipowfs].phytypesim=parms->powfs[ipowfs].phytype;
	}
	/*round phystep to be multiple of dtrat. */
	if(parms->powfs[ipowfs].phystep>0){
	    parms->powfs[ipowfs].phystep=(parms->powfs[ipowfs].phystep/parms->powfs[ipowfs].dtrat)
		*parms->powfs[ipowfs].dtrat;
	}
	if(parms->powfs[ipowfs].mtchcra==-1){
	    parms->powfs[ipowfs].mtchcra=parms->powfs[ipowfs].mtchcr;
	}
	if(parms->powfs[ipowfs].phytypesim==1 || parms->powfs[ipowfs].phytype==1){
	    int pixpsay=parms->powfs[ipowfs].pixpsa;
	    int pixpsax=parms->powfs[ipowfs].radpix;
	    if(!pixpsax) pixpsax=pixpsay;
	    if(pixpsax*pixpsay<4 && (parms->powfs[ipowfs].mtchcr>0
				     || parms->powfs[ipowfs].mtchcra>0)){
		warning("Disable constraint matched filte for quadcell\n");
		parms->powfs[ipowfs].mtchcr=0;
		parms->powfs[ipowfs].mtchcra=0;
	    }
	}
	if(parms->powfs[ipowfs].fieldstop>0 && (parms->powfs[ipowfs].fieldstop>10 || parms->powfs[ipowfs].fieldstop<1e-4)){
	    error("powfs%d: fieldstop=%g. probably wrong unit. (arcsec)\n", ipowfs, parms->powfs[ipowfs].fieldstop);
	}
	parms->powfs[ipowfs].fieldstop=parms->powfs[ipowfs].fieldstop/206265.;
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
    READ_WFS(dbl,thetax);
    READ_WFS(dbl,thetay);
    for(i=0; i<parms->nwfs; i++){
	parms->wfs[i].thetax/=206265.;
	parms->wfs[i].thetay/=206265.;
    }
    READ_WFS(int,powfs);
    double *wvlwts=0;
    int nwvlwts=readcfg_dblarr(&wvlwts,"wfs.wvlwts");
    double *siglev=0;
    int nsiglev=readcfg_dblarr(&siglev,"wfs.siglev");
    int powfs_siglev_override=readcfg_peek_override("powfs.siglev");
    int powfs_wvlwts_override=readcfg_peek_override("powfs.wvlwts");
    int count=0;
    if(nsiglev!=0 && nsiglev!=parms->nwfs){
	error("wfs.siglev can be either empty or %d\n",parms->nwfs);
    }
    for(i=0; i<parms->nwfs; i++){
	int ipowfs=parms->wfs[i].powfs;
	int nwvl=parms->powfs[ipowfs].nwvl;
	parms->wfs[i].wvlwts=malloc(nwvl*sizeof(double));
	if(nwvlwts==0){
	    memcpy(parms->wfs[i].wvlwts,parms->powfs[ipowfs].wvlwts,sizeof(double)*nwvl);
	}else{
	    memcpy(parms->wfs[i].wvlwts,wvlwts+count,sizeof(double)*nwvl);
	    count+=nwvl;
	    if(parms->powfs[ipowfs].wvlwts && powfs_wvlwts_override){
		error("when both powfs.wvlwts and wfs.wvlwts are overriden "
		      "must set powfs.wvlwts=[]\n");
	    }
	}
	if(nsiglev==0){
	    parms->wfs[i].siglev=parms->powfs[ipowfs].siglev;
	}else{
	    parms->wfs[i].siglev=siglev[i];
	    if(parms->powfs[ipowfs].siglev>0 && powfs_siglev_override){
		error("when both powfs.siglev and wfs.siglev are overriden "
		      "must set powfs.siglev=[]\n");
	    }
	}
    }
    if(nsiglev>0){
	free(siglev);
    }
    free(wvlwts);
    if(count!=nwvlwts){
	error("Supplied %d wvlwts but need %d for all wfs.\n",nwvlwts,count);
    }
    free(dbltmp);
    free(inttmp);
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
    ndm=parms->ndm=readcfg_peek_n("dm.order");
    parms->dm=calloc(parms->ndm,sizeof(struct DM_CFG_T));
    int* inttmp=NULL;
    double *dbltmp=NULL;
    char **strtmp=NULL;
    READ_DM(int,order);
    READ_DM(dbl,ht);
    READ_DM(dbl,offset);
    READ_DM_RELAX(dbl,guard);
    READ_DM_RELAX(dbl,stroke);
    READ_DM_RELAX(dbl,vmisreg);
    READ_DM_RELAX(dbl,histbin);
    READ_DM_RELAX(int,histn);
    READ_DM_RELAX(int,hist); 
    READ_DM_RELAX(int,cubic);
    READ_DM_RELAX(dbl,iac);
    READ_DM_RELAX(str,misreg);
    READ_DM_RELAX(str,hyst);
    READ_DM_RELAX(str,actfloat);
    READ_DM_RELAX(str,actstuck);
    free(strtmp);
    free(inttmp);
    free(dbltmp);
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
    int nmoao,i;
    nmoao=parms->nmoao=readcfg_peek_n("moao.order");
    parms->moao=calloc(nmoao, sizeof(MOAO_CFG_T));
    int *inttmp=NULL;
    double *dbltmp=NULL;
    char **strtmp=NULL;
    READ_MOAO(int,order);
    READ_MOAO_RELAX(int,cubic);
    READ_MOAO_RELAX(dbl,iac);
    READ_MOAO_RELAX(dbl,gdm);
    READ_MOAO_RELAX(dbl,stroke);
    READ_MOAO_RELAX(int,actslave);
    READ_MOAO_RELAX(int,lrt_ptt);
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
    READ_DBL(atm.l0);
    READ_DBL(atm.dx);
    READ_INT(atm.wdrand);
    READ_INT(atm.fractal);
    READ_INT(atm.evolve);
    READ_INT(atm.frozenflow);
    READ_INT(atm.ninit);
    READ_INT(atm.share);
    readcfg_dblarr_n(&(parms->atm.size),2,"atm.size");
    parms->atm.nps=readcfg_dblarr(&(parms->atm.ht),"atm.ht");
    readcfg_dblarr_n(&(parms->atm.wt),parms->atm.nps,"atm.wt");
    readcfg_dblarr_n(&(parms->atm.ws),parms->atm.nps,"atm.ws");
    readcfg_dblarr_n(&(parms->atm.wddeg),parms->atm.nps,"atm.wddeg");
    for(int ih=0; ih<parms->atm.nps; ih++){
	if(fabs(parms->atm.wddeg[ih])>1){
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
    READ_DBL(atmr.l0);
    if(parms->atmr.l0<=0){
	parms->atmr.l0=parms->atm.l0;
    }
    parms->atmr.nps=readcfg_dblarr(&(parms->atmr.ht),"atmr.ht");
    readcfg_dblarr_n(&(parms->atmr.wt), parms->atmr.nps, "atmr.wt");
    readcfg_intarr_n(&(parms->atmr.os), parms->atmr.nps, "atmr.os");
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
    READ_STR(aper.fnamp);
    READ_STR(aper.pupmask);
    readcfg_dblarr_n(&parms->aper.misreg, 2, "aper.misreg");
    if(fabs(parms->aper.misreg[0])>EPS || fabs(parms->aper.misreg[1])>EPS){
	parms->aper.ismisreg=1;
    }
}

/**
   Read in performance evaluation science point parameters.
*/
static void readcfg_evl(PARMS_T *parms){
    parms->evl.nevl=readcfg_peek_n("evl.thetax");
    readcfg_dblarr_n(&(parms->evl.thetax),parms->evl.nevl, "evl.thetax");
    readcfg_dblarr_n(&(parms->evl.thetay),parms->evl.nevl, "evl.thetay");
    readcfg_dblarr_n(&(parms->evl.wt),parms->evl.nevl, "evl.wt");
    readcfg_dblarr_nmax(&(parms->evl.hs), parms->evl.nevl, "evl.hs");
    readcfg_dblarr_n(&(parms->evl.misreg),2, "evl.misreg");
    if(fabs(parms->evl.misreg[0])>EPS || fabs(parms->evl.misreg[1])>EPS){
	parms->evl.ismisreg=1;
    }
    normalize_sum(parms->evl.wt, parms->evl.nevl, 1);
    readcfg_intarr_nmax(&(parms->evl.psf), parms->evl.nevl, "evl.psf");
    readcfg_intarr_nmax(&(parms->evl.psfr), parms->evl.nevl, "evl.psfr");
    parms->evl.nwvl = readcfg_dblarr(&(parms->evl.wvl), "evl.wvl");
    readcfg_intarr_nmax(&(parms->evl.psfgridsize), parms->evl.nwvl, "evl.psfgridsize");
    readcfg_intarr_nmax(&(parms->evl.psfsize), parms->evl.nwvl, "evl.psfsize");
    int ievl;
    double ramin=INFINITY;
    for(ievl=0; ievl<parms->evl.nevl; ievl++){
	/*First Convert theta to radian from arcsec. */
	parms->evl.thetax[ievl]/=206265.;
	parms->evl.thetay[ievl]/=206265.;
	double ra2=pow(parms->evl.thetax[ievl], 2)+pow(parms->evl.thetay[ievl], 2);
	if(ra2<ramin){
	    parms->evl.indoa=ievl;
	    ramin=ra2;
	}
    }
    READ_DBL(evl.dx);
    READ_INT(evl.rmax);
    READ_INT(evl.psfol);
    READ_INT(evl.psfisim);
    readcfg_intarr_nmax(&parms->evl.psfpttr, parms->evl.nevl, "evl.psfpttr");
    readcfg_intarr_nmax(&parms->evl.psfngsr, parms->evl.nevl, "evl.psfngsr");
    READ_INT(evl.psfmean); 
    READ_INT(evl.psfhist); 
    READ_INT(evl.opdcov);/*Science OPD covariance. */
    READ_INT(evl.tomo);
    READ_INT(evl.moao);
    /*it is never good to parallelize the evl ray tracing because it is already so fast */
    parms->evl.nthread=1;/*parms->sim.nthread; */
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
    READ_INT(tomo.ahst_rtt);
    READ_INT(tomo.alg);
    READ_INT(tomo.bgs);
    READ_INT(tomo.precond);
    READ_INT(tomo.maxit);
    READ_INT(tomo.assemble);
    READ_INT(tomo.predict);
    READ_DBL(tomo.minwt);
    READ_INT(tomo.cubic);
    READ_DBL(tomo.iac);
    READ_INT(tomo.ninit);
    READ_INT(tomo.psol);
    READ_INT(tomo.nxbase);
}

/**
   Read in DM fit parameters. MOAO is specified elsewhere in readcfg_moao() */
static void readcfg_fit(PARMS_T *parms){
    parms->fit.nfit=readcfg_peek_n("fit.thetax");
    readcfg_dblarr_n(&(parms->fit.thetax), parms->fit.nfit, "fit.thetax");
    readcfg_dblarr_n(&(parms->fit.thetay), parms->fit.nfit, "fit.thetay");
    readcfg_dblarr_n(&(parms->fit.wt), parms->fit.nfit, "fit.wt");
    readcfg_dblarr_nmax(&(parms->fit.hs), parms->fit.nfit, "fit.hs");
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	parms->fit.thetax[ifit]/=206265.;
	parms->fit.thetay[ifit]/=206265.;
    }

    READ_DBL(fit.tikcr);
    READ_DBL(fit.svdthres);
    READ_INT(fit.actslave);
    READ_INT(fit.lrt_piston);
    READ_INT(fit.lrt_tt);
    READ_INT(fit.alg);
    READ_INT(fit.bgs);
    READ_INT(fit.precond);
    READ_INT(fit.maxit);
    READ_INT(fit.square);
    READ_INT(fit.assemble);
    READ_INT(fit.pos);
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
}
/**
   Read general reconstruction parameters
 */
static void readcfg_recon(PARMS_T *parms){
    READ_INT(recon.alg);
    READ_INT(recon.glao);
    READ_INT(recon.split);
    READ_INT(recon.mvm);
}
/**
   Read in simulation parameters
*/
static void readcfg_sim(PARMS_T *parms){
    parms->sim.apupt=readcfg_dmat("sim.apupt");
    parms->sim.epupt=readcfg_dmat("sim.epupt");
    READ_DBL(sim.epfocus);
    READ_DBL(sim.lpfocus);
    READ_INT(sim.mffocus);
    READ_INT(sim.uptideal);

    parms->sim.apdm=readcfg_dmat("sim.apdm");
    parms->sim.epdm=readcfg_dmat("sim.epdm");
    parms->sim.aplo=readcfg_dmat("sim.aplo");
    parms->sim.eplo=readcfg_dmat("sim.eplo");
    /*We append a 0 so that we keep a time history of the integrator. */
    if(parms->sim.apdm->nx==1){
	dresize(parms->sim.apdm, 2, 1);
    }
    if(parms->sim.aplo->nx==1){
	dresize(parms->sim.aplo, 2, 1);
    }
    if(fabs(dsum(parms->sim.apdm)-1)>1.e-10){
	warning("sum(sim.apdm)=%g. Should be 1.\n", dsum(parms->sim.apdm));
    }
    if(fabs(dsum(parms->sim.aplo)-1)>1.e-10){
	warning("sum(sim.aplo)=%g. Should be 1.\n", dsum(parms->sim.aplo));
    }
    if(fabs(dsum(parms->sim.apupt)-1)>1.e-10){
	warning("sum(sim.apupt)=%g. Should be 1.\n", dsum(parms->sim.apupt));
    }
    parms->sim.nseed=readcfg_intarr(&parms->sim.seeds,"sim.seeds");
    READ_DBL(sim.dt);
    READ_INT(sim.start);
    READ_INT(sim.end);
    READ_STR(sim.wspsd);
    READ_INT(sim.wsseq);
    READ_INT(sim.cachedm);
    READ_INT(sim.fuseint);
    READ_INT(sim.closeloop);
    READ_INT(sim.skysim);
    READ_DBL(sim.fov);
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
    READ_INT(sim.parallel);
}
/**
   Read in parameters for Cn2 estimation.
 */
static void readcfg_cn2(PARMS_T *parms){
/*for Cn2 Estimation. */
    parms->cn2.npair = readcfg_intarr(&parms->cn2.pair,"cn2.pair");
    READ_INT(cn2.step);
    READ_INT(cn2.reset);
    READ_INT(cn2.tomo);
    READ_INT(cn2.verbose);
    READ_INT(cn2.keepht);
    READ_INT(cn2.nhtomo);
    READ_INT(cn2.moveht);
    READ_DBL(cn2.hmax);
    READ_DBL(cn2.saat);
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
	parms->plot.setup=1;
	parms->plot.atm=1;
	parms->plot.run=1;
    }
    if(parms->plot.setup || parms->plot.atm || parms->plot.run || parms->plot.opdx || parms->plot.all){
	draw_helper();
    }
}
/**
   Read in debugging parameters
*/
static void readcfg_dbg(PARMS_T *parms){
  
    READ_INT(dbg.psol);
    READ_INT(dbg.wamethod);
    READ_INT(dbg.atm);
    READ_INT(dbg.mvstlimit);
    READ_INT(dbg.annular_W);
    parms->dbg.ntomo_maxit=readcfg_intarr(&parms->dbg.tomo_maxit, "dbg.tomo_maxit");
    READ_INT(dbg.tomo_hxw);
    READ_INT(dbg.splitlrt);
    READ_INT(dbg.ecovxx);
    READ_INT(dbg.useopdr);
    READ_INT(dbg.usegwr);
    READ_INT(dbg.dxonedge);
    READ_INT(dbg.cmpgpu);
    READ_INT(dbg.pupmask);
    READ_INT(dbg.wfslinearity);
    READ_INT(dbg.nocgwarm);
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
    READ_INT(save.all);
    READ_INT(save.setup);
    READ_INT(save.recon);
    READ_INT(save.mvst);

    READ_INT(save.atm);/*Save atmosphere */
    READ_INT(save.run);
    READ_INT(save.opdr);/*reconstructed OPD on XLOC */
    READ_INT(save.opdx);/*ATM propagated to XLOC */
    READ_INT(save.evlopd);/*Science OPD */
    READ_INT(save.dm);/*save DM commands */
    READ_INT(save.dmpttr);
    readcfg_intarr_nmax(&parms->save.ints, parms->nwfs, "save.ints");
    readcfg_intarr_nmax(&parms->save.wfsopd, parms->nwfs, "save.wfsopd");
    readcfg_intarr_nmax(&parms->save.grad, parms->nwfs, "save.grad");
    readcfg_intarr_nmax(&parms->save.gradgeom, parms->nwfs, "save.gradgeom");
  
    if(parms->save.all){/*enables everything */
	warning("Enabling saving everything.\n");
	/*The following 3 are for setup. */
	parms->save.setup=1;
	parms->save.recon=1;
	parms->save.mvst=1;
	/*The following are run time information that are not enabled by
	  save.run because they take a lot of space*/
	parms->save.opdr=1;
	parms->save.opdx=1;
	parms->save.evlopd=1;
	parms->save.dmpttr=1;
	
	parms->save.run=1;/*see following */
    }

    if(parms->save.run){
	parms->save.dm=1;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    parms->save.ints[iwfs]=1;
	    parms->save.wfsopd[iwfs]=1;
	    parms->save.grad[iwfs]=1;
	    parms->save.gradgeom[iwfs]=1;
	}
    }
    parms->save.ngcov=readcfg_intarr(&parms->save.gcov,"save.gcov")/2;
    READ_INT(save.gcovp);
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
    READ_STR(load.MVM);
    READ_INT(load.mvst);
    READ_INT(load.GS0);
    READ_INT(load.tomo);
    READ_INT(load.fit);
    READ_INT(load.W);
    READ_INT(load.i0);
}
/**
   Process simulation parameters to find incompatibility.
*/
static void setup_parms_postproc_sim(PARMS_T *parms){
    if(parms->sim.skysim){
	if(parms->recon.alg!=0){
	    error("skysim need MVR");
	}
	parms->tomo.ahst_idealngs=1;
	if(parms->tomo.ahst_wt!=3){
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
		parms->powfs[ipowfs].pistatout=1;
	    }
	}
    }
    if(parms->dbg.ntomo_maxit){
	warning("dbg.tomo_maxit is set. Will run in open loop mode\n to repeat the simulations"
		" with different values of tomo.maxit.\n");
	parms->sim.closeloop=0;
	parms->atm.frozenflow=1;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    parms->atm.ws[ips]=0;/*set windspeed to zero. */
	}
	parms->sim.end=parms->dbg.ntomo_maxit;
    }
    if(parms->sim.idealfit){
	if(parms->recon.alg!=0){
	    warning("idealfit only works in recon.alg=0 mode. changed\n");
	    parms->recon.alg=0;
	}
	if(parms->recon.split){
	    warning("idealfit only works in integrated tomo mode. changed\n");
	    parms->recon.split=0;
	}
	if(parms->recon.mvm){
	    parms->recon.mvm=0;
	}
    }
    if(parms->recon.glao && parms->ndm!=1){
	error("GLAO only works with 1 dm\n");
    }
    if(parms->recon.alg==1 && parms->recon.split==2){
	info2("MVST does not work with least square reconstructor. Changed to AHST");
	parms->recon.split=1;
    }
    if(parms->sim.wfsalias){
	if(parms->sim.idealwfs){
	    error("sim.wfsalias conflicts with sim.idealwfs. Do not enable both.\n");
	}
	if(parms->sim.idealevl){
	    error("sim.wfsalias conflicts with sim.idealevl. Do not enable both.\n");
	}
    }
    if(parms->sim.idealtomo){
	if(parms->sim.wfsalias){
	    error("wfsalias and idealtomo conflicts\n");
	}
	if(parms->sim.idealfit){
	    error("idealfit and idealtomo conflicts\n");
	}
	if(parms->sim.idealwfs){
	    error("idealwfs and idealtomo conflicts\n");
	}
    }
    if(parms->sim.wfsalias || parms->sim.idealwfs || parms->sim.idealevl){
	parms->sim.dmproj=1;/*need dmproj */
    }
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
	    parms->atm.ht[ips] *= secz;/*scale atmospheric height */
	}
	//parms->atm.hmax*=secz;
	//parms->atmr.hmax*=secz;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(isfinite(parms->powfs[ipowfs].hs)){
		parms->powfs[ipowfs].hs *= secz;/*scale GS height. */
		for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs[indwfs];
		    double siglev=parms->wfs[iwfs].siglev;
		    parms->wfs[iwfs].siglev=siglev*cosz;/*scale signal level. */
		    info("iwfs%d: siglev scaled from %g to %g\n", 
			 iwfs,siglev,parms->wfs[iwfs].siglev);
		    warning("Need to update to account for the transmittance\n");
		}
	    }
	}
	warning("Scaling reconstruction height to zenith angle %gdeg\n",parms->sim.za*180./M_PI);
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    parms->atmr.ht[ips] *= secz;/*scale reconstructed atmospheric height. */
	}
	parms->cn2.hmax*=secz;
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
    /*link wfs with powfs*/
    int jpowfs=0;/*records last powfs.*/
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	parms->powfs[ipowfs].wfs=calloc(parms->nwfs, sizeof(int));
	parms->powfs[ipowfs].wfsind=calloc(parms->nwfs, sizeof(int));
	int count=0;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int kpowfs=parms->wfs[iwfs].powfs;
	    if(kpowfs==ipowfs){
		parms->powfs[ipowfs].wfs[count]=iwfs;
		parms->powfs[ipowfs].wfsind[iwfs]=count;
		count++;
		parms->wfs[iwfs].powfs=jpowfs; /*update the pointer*/
	    }else{
		parms->powfs[ipowfs].wfsind[iwfs]=-1;/*not belong */
	    }
	}
	parms->powfs[ipowfs].nwfs=count;
	if(count>0){/*non empty powfs*/
	    if(jpowfs<ipowfs){
		memcpy(parms->powfs+jpowfs, parms->powfs+ipowfs, sizeof(POWFS_CFG_T));
	    }
	    jpowfs++;
	}else{
	    info("Removing powfs %d\n", ipowfs);
	    free_powfs_cfg(parms->powfs+ipowfs);
	    continue;
	}
	if(parms->powfs[ipowfs].llt){
	    parms->powfs[ipowfs].llt->i=calloc(count, sizeof(int));/*default to zero. */
	    if(parms->powfs[ipowfs].llt->n>1){
		/*this is single llt for this powfs. */
		if(parms->powfs[ipowfs].llt->n!=count)
		    error("# of llts should either be 1 or match nwfs for this powfs");
		for(int iwfs=0; iwfs<parms->powfs[ipowfs].llt->n; iwfs++){
		    parms->powfs[ipowfs].llt->i[iwfs]=iwfs;
		}
	    }
	}
    }
    parms->npowfs=jpowfs;
    if(jpowfs==0){
	warning("No wfs is found\n");
	if(!parms->sim.idealfit && !parms->sim.evlol){
	    error("Cannot proceed\n");
	}
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	/*Figure out order of High order WFS if not specified.*/
	if(parms->powfs[ipowfs].order==0){
	    if(parms->ndm>0){
		parms->powfs[ipowfs].order=parms->dm[0].order;
	    }else{
		error("Please specify powfs[%d].order in MOAO mode\n", ipowfs);
	    }
	}
	if(parms->powfs[ipowfs].nwfs>0){
	    if(parms->powfs[ipowfs].lo){
		parms->nlopowfs++;
		if(parms->powfs[ipowfs].trs==1){
		    error("Low order wfs should not be tilt removed\n");
		}
	    }else{
		parms->nhipowfs++;
	    }
	    if(parms->powfs[ipowfs].trs){
		if(parms->powfs[ipowfs].lo){
		    error("WFS with tip/tilt removed should be high order\n");
		}
		parms->ntrspowfs++;
	    }
	}
	/* 
	   Figure out pixtheta if specified to be auto (<0).
	   -pixtheta is the ratio to nominal value.
	*/
	const double dxsa 
	    = parms->aper.d/(double)parms->powfs[ipowfs].order;
	double wvl=0;
	for(int iwvl=0; iwvl<parms->powfs[ipowfs].nwvl; iwvl++){
	    if(parms->powfs[ipowfs].wvl[iwvl]>wvl)
		wvl=parms->powfs[ipowfs].wvl[iwvl];
	}
	/*should supply in arcsec. Was supplied in radian pre 2011-02-17. So we
	  test magnitude and then apply conversion*/
	if(parms->powfs[ipowfs].pixtheta>1e-4){
	    parms->powfs[ipowfs].pixtheta/=206265.;/*convert form arcsec to radian. */
	}else if(parms->powfs[ipowfs].pixtheta>0){
	    warning2("ipowfs %d:pixtheta should be supplied in unit of arcsec.\n", ipowfs);
	}
	if(parms->powfs[ipowfs].pixtheta<0){
	    parms->powfs[ipowfs].dl=1;/*mark as diffraction limited. */
	    double ratio=(-parms->powfs[ipowfs].pixtheta);
	    parms->powfs[ipowfs].pixtheta=ratio*wvl/dxsa;
	    info2("powfs %d pixtheta set to %.1fx %g/%g: %g mas\n",
		  ipowfs, ratio, wvl,dxsa,parms->powfs[ipowfs].pixtheta*206265000);
	}else if(parms->powfs[ipowfs].pixtheta<wvl/dxsa*1.22){
	    parms->powfs[ipowfs].dl=1;/*mark as diffraction limited. */
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
	    parms->powfs[ipowfs].hasGS0=1;
	}else{
	    parms->powfs[ipowfs].hasGS0=0;
	}
	/*Do we ever do physical optics.*/
	if(parms->powfs[ipowfs].phystep>=0
	   &&(parms->powfs[ipowfs].phystep<parms->sim.end||parms->sim.end==0)){
	    parms->powfs[ipowfs].usephy=1;
	}else{
	    parms->powfs[ipowfs].usephy=0;
	}
	if(!parms->powfs[ipowfs].usephy && parms->powfs[ipowfs].bkgrndfn){
	    warning("powfs %d: there is sky background, but is using geometric wfs. background won't be effective.\n", ipowfs);
	} 
	if(parms->powfs[ipowfs].ncpa && parms->powfs[ipowfs].mtchstc){
	    warning("powfs %d: Disabling shifting i0 to center in the presence of NCPA.\n", ipowfs);
	    parms->powfs[ipowfs].mtchstc=0;
	}
	if(parms->powfs[ipowfs].ncpa && !parms->powfs[ipowfs].usephy 
	   && parms->powfs[ipowfs].ncpa_method==2){
	    warning("powfs %d: ncpa_method changed from 2 to 1 in geometric wfs mdoe\n", ipowfs);
	    parms->powfs[ipowfs].ncpa_method=1;
	}
	{
	    /*Adjust dx if the subaperture does not contain integer, even number of points.*/
	    int nx = 2*(int)round(0.5*dxsa/parms->powfs[ipowfs].dx);
	    double dx=dxsa/nx;/*adjust dx. */
	    if(fabs(parms->powfs[ipowfs].dx-dx)>EPS){
		warning("powfs %d: Adjusting dx from %g to %g. \n"
			"Please adjust evl.dx, powfs.dx, atm.dx to match the new value for best efficiency.\n",
			ipowfs,parms->powfs[ipowfs].dx, dx);
	    }
	    parms->powfs[ipowfs].dx=dx;
	}
	if(!parms->sim.closeloop && parms->powfs[ipowfs].dtrat!=1){
	    warning("powfs %d: in open loop mode, only dtrat=1 is supported. Changed\n", ipowfs);
	    parms->powfs[ipowfs].dtrat=1;
	}
	if(parms->sim.wfsalias){
	    if(parms->powfs[ipowfs].noisy){
		warning("powfs%d: wfsalias is set. make it noise free\n", ipowfs);
		parms->powfs[ipowfs].noisy=0;
	    }
	    if(parms->powfs[ipowfs].phystep>=parms->sim.start || parms->powfs[ipowfs].phystep<parms->sim.end){
		warning("powfs%d: wfsalias is set. make it geometric wfs since the mtched fitler won't work well\n", ipowfs);
		parms->powfs[ipowfs].phystep=-1;
	    }
	}
	if(parms->powfs[ipowfs].mtchscl==-1){
	    if(fabs(parms->powfs[ipowfs].sigscale-1)>EPS){
		parms->powfs[ipowfs].mtchscl=1;
	    }else{
		parms->powfs[ipowfs].mtchscl=0;
	    }
	}
    }

    if(parms->recon.split){
	if(parms->nlopowfs==0){
	    if(parms->ntrspowfs>=parms->nhipowfs){
		warning("There is no WFS controlling tip/tilt.\n");
	    }else{
		warning("Split reconstruction is enabled when there is no low order WFS. Will split the tip/tilt modes from high order wfs\n");
	    }
	}
	int hi_found=0;
	int hi_hastt=0;
	int lo_found=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs>0){
		if(parms->powfs[ipowfs].lo){
		    lo_found=1;
		    if(parms->powfs[ipowfs].trs==1){
			error("Low order wfs should not be tilt removed\n");
		    }
		}else{
		    hi_found=1;
		    if(parms->powfs[ipowfs].trs!=1){
			hi_hastt=1;
		    }
		}
	    }
	}
	if(parms->sim.skysim && (parms->nhipowfs==0 || parms->nlopowfs==0)){
	    error("There is only high or low order WFS. can not do skycoverage presimulation\n");
	}
	if(!parms->nhipowfs){
	    warning("There is no high order WFS!!!\n");
	}
    }
   
    if((parms->recon.split) && parms->ndm==0){
	warning("Disable split tomography since there is no common DM\n");
	parms->recon.split=0;
    }

    if(parms->recon.glao){
	parms->wfsr=calloc(parms->npowfs, sizeof(WFS_CFG_T));
	parms->nwfsr=parms->npowfs;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    parms->wfsr[ipowfs].thetax=0;
	    parms->wfsr[ipowfs].thetay=0;
	    parms->wfsr[ipowfs].powfs=ipowfs;
	}
    }else{/*Use same information as wfs. */
	parms->wfsr = parms->wfs;
	parms->nwfsr= parms->nwfs;
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].lo && parms->powfs[ipowfs].gtype_sim==0)
	    warning("Low order POWFS %d is using gtilt in simulation. "
		    "This is not recommended\n",ipowfs);
    }
    /*
    parms->sim.ndtrat=0;
    parms->sim.dtrats=calloc(parms->npowfs, sizeof(int));
    parms->sim.idtrat_hi=-1;
    parms->sim.idtrat_lo=-1;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int idtrat;
	for(idtrat=0; idtrat<parms->sim.ndtrat; idtrat++){
	    if(parms->powfs[ipowfs].dtrat==parms->sim.dtrats[idtrat]){
		parms->powfs[ipowfs].idtrat=idtrat;
		break;
	    }
	}
	if(idtrat==parms->sim.ndtrat){
	    parms->powfs[ipowfs].idtrat=idtrat;
	    parms->sim.dtrats[idtrat]=parms->powfs[ipowfs].dtrat;
	    parms->sim.ndtrat++;
	}
	if(!parms->powfs[ipowfs].lo){
	    if(parms->sim.idtrat_hi==-1){
		parms->sim.idtrat_hi=parms->powfs[ipowfs].idtrat;
	    }else if(parms->sim.idtrat_hi>-1){
		if(parms->sim.idtrat_hi!=parms->powfs[ipowfs].idtrat){
		    warning("High order WFS has multiple dtrats\n");
		    parms->sim.idtrat_hi=-2;
		}
	    }
	}else{
	    if(parms->sim.idtrat_lo==-1){
		parms->sim.idtrat_lo=parms->powfs[ipowfs].idtrat;
	    }else if(parms->sim.idtrat_lo>-1){
		if(parms->sim.idtrat_lo!=parms->powfs[ipowfs].idtrat){
		    warning("High order WFS has multiple dtrats\n");
		    parms->sim.idtrat_lo=-2;
		}
	    }
	}
	}*/
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
	if(parms->atm.wt[ips]>1.e-3){
	    if(ips!=jps){
		parms->atm.ht[jps]=parms->atm.ht[ips];
		parms->atm.wt[jps]=parms->atm.wt[ips];
		parms->atm.ws[jps]=parms->atm.ws[ips];
		parms->atm.wddeg[jps]=parms->atm.wddeg[ips];
	    }
	    jps++;
	}else{
	    warning("Layer %d has very small weight of %g, will drop it.\n", 
		    ips, parms->atm.wt[ips]);
	}
    }
    if(jps==0){
	error("There are no valid atmosphere layer\n");
    }
    if(parms->atm.nps!=jps){
	warning("nlayer changed from %d to %d\n", parms->atm.nps, jps);
	parms->atm.nps=jps;
	parms->atm.ht=realloc(parms->atm.ht, sizeof(double)*jps);
	parms->atm.wt=realloc(parms->atm.wt, sizeof(double)*jps);
	parms->atm.ws=realloc(parms->atm.ws, sizeof(double)*jps);
	parms->atm.wddeg=realloc(parms->atm.wddeg, sizeof(double)*jps);
    }
    if(parms->sim.idealfit){/*If fit only, we using atm for atmr. */
	warning("Changing atmr.ht,wt to atm.ht,wt since we are doing fit only\n");
	int nps=parms->atm.nps;
	parms->atmr.ht=realloc(parms->atmr.ht, sizeof(double)*nps);
	parms->atmr.wt=realloc(parms->atmr.wt, sizeof(double)*nps);
	memcpy(parms->atmr.ht, parms->atm.ht, sizeof(double)*nps);
	memcpy(parms->atmr.wt, parms->atm.wt, sizeof(double)*nps);
	parms->atmr.os=realloc(parms->atmr.os, sizeof(int)*nps);
	for(int ips=parms->atmr.nps; ips<nps; ips++){
	    parms->atmr.os[ips]=parms->atmr.os[parms->atmr.nps-1];
	}
	parms->atmr.nps=nps;
    }
    if(parms->recon.glao && parms->recon.alg==0){
	/*GLAO mode. reconstruct only a single layer near the DM. Using only 1 fitting direction on axis.*/
	warning2("In GLAO Mode, use 1 tomography grid near the ground dm and 1 on axis fitting direction\n");
	parms->atmr.ht=realloc(parms->atmr.ht, sizeof(double));
	parms->atmr.wt=realloc(parms->atmr.wt, sizeof(double));
	parms->atmr.os=realloc(parms->atmr.os, sizeof(double));
	parms->atmr.ht[0]=parms->dm[0].ht;
	parms->atmr.wt[0]=1;
	parms->atmr.nps=1;
	
	parms->fit.nfit=1;
	parms->fit.thetax=realloc(parms->fit.thetax, sizeof(double));
	parms->fit.thetay=realloc(parms->fit.thetay, sizeof(double));
	parms->fit.wt=realloc(parms->fit.wt, sizeof(double));
	parms->fit.hs=realloc(parms->fit.hs, sizeof(double));
	parms->fit.thetax[0]=0;
	parms->fit.thetay[0]=0;
	parms->fit.wt[0]=1;
    }
    normalize_sum(parms->atm.wt, parms->atm.nps, 1);
    normalize_sum(parms->atmr.wt, parms->atmr.nps, 1);
    normalize_sum(parms->fit.wt, parms->fit.nfit, 1);
    /*
      We don't drop weak turbulence layers in reconstruction. Instead, we make
      it as least parms->tomo.minwt in setup_recon_tomo_prep
    */
 
    /*
      Find ground turbulence layer. The ray tracing can be shared between different directions.
    */
    parms->atm.iground=-1;
    parms->atm.hmax=-INFINITY;
    for(int ips=0; ips<parms->atm.nps; ips++){
	if(fabs(parms->atm.ht[ips])<1.e-10){
	    if(parms->atm.iground==-1){
		parms->atm.iground=ips;
	    }else{
		error("Multiple grounds atm. Please combine them together.\n");
	    }
	}
	if(parms->atm.ht[ips]<0){
	    warning("Layer %d height %g is below ground\n",ips,parms->atm.ht[ips]);
	}
	if(ips > 0 && fabs(parms->atm.ht[ips]-parms->atm.ht[ips-1])<20){
	    warning("Layer %d at %gm is too close to layer %d at %gm\n",
		    ips,parms->atm.ht[ips], ips-1, parms->atm.ht[ips-1]);
	}
	if(parms->atm.hmax<parms->atm.ht[ips]){
	    parms->atm.hmax=parms->atm.ht[ips];
	}
    }
    parms->atmr.hmax=-INFINITY;
    for(int ips=0; ips<parms->atmr.nps; ips++){
	if(parms->atmr.hmax<parms->atmr.ht[ips]){
	    parms->atmr.hmax=parms->atmr.ht[ips];
	}
    }
    if(parms->atm.iground==-1){
	warning("There is no ground layer\n");
    }
    parms->atm.frozenflow = (parms->atm.frozenflow || parms->sim.closeloop);
    if(parms->atm.fractal && !parms->atm.evolve && parms->atm.frozenflow){
	warning("atm.fractal requires atm.evolve=1, changed\n");
	parms->atm.evolve=1;
    }
    if((!parms->atm.frozenflow || parms->dbg.atm>0) && parms->atm.evolve){
	warning("Disable atm.evolve since there is no frozenflow or in debug atm mode.\n");
	parms->atm.evolve=0;/*disable evolve in open loop mode. */
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

/**
   compute minimum size of atm screen to cover all the beam path. same for
   all layers.  todo:may need to consider l0 Must be after
   setup_parms_postproc_za.
*/
static void setup_parms_postproc_atm_size(PARMS_T *parms){
    const int nps=parms->atm.nps;
    int Nmax=0;
    long nxout[nps],nyout[nps];
    for(int ips=0; ips<nps; ips++){
	create_metapupil(parms,parms->atm.ht[ips],parms->atm.dx,0.5,
			 &nxout[ips],&nyout[ips],NULL,NULL,NULL,parms->atm.dx*3,0,0,0,1);
	if(nxout[ips]>Nmax) Nmax=nxout[ips];
	if(nyout[ips]>Nmax) Nmax=nyout[ips];
    }
    /*Minimum screen size required. Used to transport atm to GPU. */
    parms->atm.nxn=Nmax;
    parms->atm.nyn=Nmax;
    parms->atm.nxm=nextpow2(Nmax);
    parms->atm.nym=parms->atm.nxm;
    if(parms->atm.l0 > parms->atm.nxm * parms->atm.dx){
	warning("Atmospheric size is smaller than outer scale!\n");
    }
    if(fabs(parms->atm.size[0])<EPS ||fabs(parms->atm.size[1])<EPS){
	parms->atm.nx=parms->atm.nxm;
	parms->atm.ny=parms->atm.nym;
    }else{/*user specified.*/
	parms->atm.nx=2*(int)round(0.5*parms->atm.size[0]/parms->atm.dx);
	parms->atm.ny=2*(int)round(0.5*parms->atm.size[1]/parms->atm.dx);
	if(parms->atm.nx<parms->atm.nxm) parms->atm.nx=parms->atm.nxm;
	if(parms->atm.ny<parms->atm.nym) parms->atm.ny=parms->atm.nym;
    }
    if(parms->atm.fractal){/*must be square and 1+power of 2 */
	int nn=parms->atm.nx>parms->atm.ny?parms->atm.nx:parms->atm.ny;
	parms->atm.nx=1+nextpow2(nn);
	parms->atm.ny=parms->atm.nx;
    }
    /*record the size of the atmosphere. */
    parms->atm.size[0]=parms->atm.nx*parms->atm.dx;
    parms->atm.size[1]=parms->atm.ny*parms->atm.dx;
    /*for screen evolving. */
    parms->atm.overx = calloc(parms->atm.nps, sizeof(long));
    parms->atm.overy = calloc(parms->atm.nps, sizeof(long));
    for(int ips=0; ips<parms->atm.nps; ips++){
	parms->atm.overx[ips] = nxout[ips];
	parms->atm.overy[ips] = nyout[ips];
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
    int ndm=parms->ndm;
    for(int idm=0; idm<ndm; idm++){
	parms->dm[idm].dx=parms->aper.d/parms->dm[idm].order;
    }
    /*
      Setup the parameters used to do DM caching on a finer grid.
    */
    for(int idm=0; idm<ndm && parms->sim.cachedm; idm++){
	double ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
	if(fabs(ht)<1.e-10){
	    parms->dm[idm].isground=1;
	}
	int nscale=0;
	int nscalemax=parms->npowfs+parms->evl.nevl;/*maximum number of possible scalings */
	double scale[nscalemax];
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    double dxscl=(1.-ht/parms->powfs[ipowfs].hs)*parms->powfs[ipowfs].dx;
	    if(idm==0){
		parms->powfs[ipowfs].scalegroup=calloc(ndm,sizeof(int));
	    }
	    if(parms->dm[idm].dx>parms->powfs[ipowfs].dx*10){
		parms->powfs[ipowfs].scalegroup[idm]=arrind(scale, &nscale, dxscl);
	    }else{
		parms->powfs[ipowfs].scalegroup[idm]=-1;
	    }
	}
	/*evl; */
	if(parms->dm[idm].dx>parms->evl.dx*10){
	    if(idm==0){
		parms->evl.scalegroup=calloc(ndm*parms->evl.nevl, sizeof(int)); 
	    }
	    for(int ievl=0; ievl<parms->evl.nevl ;ievl++){
		double dxscl=(1. - ht/parms->evl.hs[ievl])*parms->evl.dx;
		parms->evl.scalegroup[idm+ievl*ndm]=arrind(scale, &nscale, dxscl);
	    }
	}
	parms->dm[idm].ncache=nscale;
	parms->dm[idm].dxcache=calloc(1, nscale*sizeof(double));
	for(int iscale=0; iscale<nscale; iscale++){
	    parms->dm[idm].dxcache[iscale]=scale[iscale];
	}
    }
}

/**
   Setting up the cone coordinate for MCAO LGS
   simulation. First find out the guide star conjugate. Only 1
   altitude is allowed.
*/
static void setup_parms_postproc_recon(PARMS_T *parms){    
    parms->recon.warm_restart = !parms->dbg.nocgwarm && parms->atm.frozenflow && !parms->dbg.ntomo_maxit;
    {
	double hs=INFINITY;
	/*find out the height to setup cone coordinate. */
	if(parms->tomo.cone){
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if (parms->powfs[ipowfs].lo){/*skip low order wfs */
		    continue;
		}
		/*isinf and isfinite both return 0 on inf in FreeBSD 9.0.*/
		if(isfinite(parms->powfs[ipowfs].hs)){/*at finite. */
		    if(isfinite(hs) && fabs(hs-parms->powfs[ipowfs].hs)>1.e-6){
			error("Two high order POWFS with different hs found");
		    }else{
			hs = parms->powfs[ipowfs].hs;
		    }
		}
	    }
	}
	parms->atmr.hs=hs;
    }
    {
	/*find out the sampling to setup tomography grid using the maximum order of the wfs and DMs. */
	double maxorder=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].order>maxorder){
		maxorder=parms->powfs[ipowfs].order;
	    }
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    if(parms->dm[idm].order > maxorder){
		maxorder=parms->dm[idm].order;
	    }
	}
	for(int imoao=0; imoao<parms->nmoao; imoao++){
	    if(parms->moao[imoao].order>maxorder){
		maxorder=parms->moao[imoao].order;
	    }
	}
	parms->atmr.dx=parms->aper.d/maxorder;
    }
    if(parms->sim.ecnn){
	parms->tomo.assemble=1;
    }
    if((parms->tomo.bgs || parms->tomo.alg != 1) && parms->tomo.cxx !=0){
	error("Only CG work with non L2 cxx.\n");
	parms->tomo.cxx=0;
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
    if(parms->sim.closeloop && parms->evl.tomo){
	warning("Evaluating tomography performance is best done in open loop\n");
    }
    if(parms->recon.split && parms->evl.tomo){
	warning("Evaluating tomography performance is best done with integrated tomography.\n");
    }
    if(parms->tomo.alg==1 && parms->recon.alg==0){
	if(parms->tomo.precond>1){
	    error("Invalid preconditoner\n");
	}
	/*if(parms->tomo.precond==1 && parms->tomo.piston_cr){
	    warning("FDPCG does not perform well with piston_cr=1. Disabled piston_cr\n");
	    parms->tomo.piston_cr=0;
	    }*/
	if(parms->tomo.precond==1 && parms->tomo.square!=1){
	    warning("FDPCG prefers square XLOC.\n");
	}
    }else{
	parms->tomo.precond=0;
    }
    /*check cg iterations*/
    if(parms->recon.mvm || !parms->recon.warm_restart){
	if(parms->recon.alg==0){
	    if(parms->tomo.alg==1 && parms->tomo.maxit<10){
		warning("In MVM or non warm restart mode, CG in TOMO need more iterations than %d. Make 10x bigger\n",
			parms->tomo.maxit);
		//parms->tomo.maxit*=10;
	    }
	    if(parms->fit.alg==1 && parms->fit.maxit<10){
		warning("In MVM or non warm restart mode, CG in FIT need more iterations than %d. Make 10x bigger\n",
			parms->fit.maxit);
		//parms->fit.maxit*=10;
	    }
	}else{
	    if(parms->lsr.alg==1 && parms->lsr.maxit<10){
		warning("In MVM or non warm restart mode, CG in LSR need more iterations than %d. Make 10x bigger\n",
			parms->lsr.maxit);
		//parms->lsr.maxit*=10;
	    }
	}	
    }
    if(parms->sim.mffocus && (!parms->sim.closeloop || parms->sim.idealfit)){
	warning("mffocus is set, but we are in open loop mode or doing fitting only. disable\n");
	parms->sim.mffocus=0;
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->recon.split && parms->powfs[ipowfs].lo){
	    parms->powfs[ipowfs].skip=1;
	}else{
	    parms->powfs[ipowfs].skip=0;
	}
	if(parms->sim.mffocus || parms->save.ngcov>0 || (parms->cn2.pair && !parms->powfs[ipowfs].lo)){
	    /*focus tracking or cn2 estimation, or save gradient covariance.  */
	    parms->powfs[ipowfs].psol=1;
	}else{/*no focus tracking */
	    if(parms->recon.alg==0){/*MV */
		/*low order wfs in ahst mode does not need psol. */
		if((parms->recon.split==1 && parms->powfs[ipowfs].skip) || !parms->tomo.psol){
		    parms->powfs[ipowfs].psol=0;
		}else{
		    parms->powfs[ipowfs].psol=1;
		}
	    }else{
		parms->powfs[ipowfs].psol=0;
	    }
	}
    }
    for(int idm=0; idm<parms->ndm; idm++){
	if(parms->dm[idm].hist){
	    parms->sim.dmttcast=1;
	}
	if(isfinite(parms->dm[idm].stroke)){
	    parms->sim.dmclip=1;
	}
    }
    if(parms->save.dmpttr || parms->sim.dmclip || parms->sim.dmttcast){
	parms->sim.dmttcast=1;
	if(!parms->sim.fuseint){
	    error("Sorry, clipping only works in fuseint=1 mode\n");
	}
    }
    if(parms->sim.psfr){
	int fnd=sum_intarr(parms->evl.nevl, parms->evl.psfr);
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
	    if(fabs(bkgrnd)>1.e-50){
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
   postproc misc parameters.
*/
static void setup_parms_postproc_misc(PARMS_T *parms, ARG_T *arg){
    parms->pause=arg->pause;
    {
	/*Remove seeds that are already done. */
	char fn[80];
	int iseed=0; 
	int jseed=0;
	parms->fdlock=calloc(parms->sim.nseed, sizeof(int));
	for(iseed=0; iseed<parms->sim.nseed; iseed++){
	    snprintf(fn, 80, "Res_%d.done",parms->sim.seeds[iseed]);
	    if(exist(fn) && !arg->override){
		parms->fdlock[iseed]=-1;
		warning2("Skip seed %d because %s exist.\n", parms->sim.seeds[iseed], fn);
	    }else{
	    	snprintf(fn, 80, "Res_%d.lock",parms->sim.seeds[iseed]);
		parms->fdlock[iseed]=lock_file(fn, 0, 0);
		if(parms->fdlock[iseed]<0){
		    warning2("Skip seed %d because it is already running.\n",
			     parms->sim.seeds[iseed]);
		}else{
		    cloexec(parms->fdlock[iseed]);
		    if(jseed!=iseed){
			parms->sim.seeds[jseed]=parms->sim.seeds[iseed];
			parms->fdlock[jseed]=parms->fdlock[iseed];
		    }
		    jseed++;
		}
	    }
	}
	if(jseed!=parms->sim.nseed){
	    info2("Skip %d seeds.\n", parms->sim.nseed - jseed);
	}
	parms->sim.nseed=jseed;
	if(parms->sim.nseed<1){
	    warning("There are no seed to run. Use -O to override. Exit\n");
	    scheduler_finish(0);
	    raise(SIGUSR1);
	    exit(1);
	}
    }
    info2("There are %d valid simulation seeds: ",parms->sim.nseed);
    for(int i=0; i<parms->sim.nseed; i++){
	info2(" %d", parms->sim.seeds[i]);
    }
    info2("\n");
  
    if(parms->save.ngcov>0 && parms->save.gcovp<10){
	warning("parms->save.gcovp=%d is too small. It may fill your disk!\n",
		parms->save.gcovp);
    }
    if(parms->save.gcovp>parms->sim.end){
	parms->save.gcovp=parms->sim.end;
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
    /*load.tomo==1 loads tomography matrix saved by LAOS for debugging purpose.*/
    if(parms->load.tomo || parms->tomo.alg!=1 || parms->tomo.bgs){
	parms->tomo.assemble=1;
    }
    if(parms->load.fit || parms->fit.alg!=1 || parms->fit.bgs){
	parms->fit.assemble=1;
    }
    if(parms->tomo.bgs && parms->tomo.precond){
	error("Please implement the preconditioner for each block for BGS.\n");
    }
    if(parms->fit.bgs && parms->fit.precond){
	error("Please implement the preconditioner for each block for BGS.\n");
    }
    if(parms->recon.alg==0 && parms->tomo.predict==1 && parms->tomo.alg!=1){
	error("Predictive tomography only works with CG. need to redo CBS/MVM after wind velocity is know.\n");
    }
    /*disable cache for low order systems. */
    if(parms->evl.nevl<2){
	if(parms->sim.cachedm==1){
	    parms->sim.cachedm=0;
	    warning("cachedm disabled for SCAO\n");
	}
    }
    if(parms->evl.tomo){
	warning("evl.tomo in cuda not implemented yet.\n");
	use_cuda=0;
    }
    if(use_cuda){
	if(parms->sim.evlol){
	    parms->gpu.lsr=parms->gpu.tomo=parms->gpu.fit=parms->gpu.moao=parms->gpu.wfs=0;
	}
	if(parms->sim.idealfit){
	    parms->gpu.tomo=0;/*no need tomo.*/
	}
	if(parms->recon.alg==0){/*MV*/
	    if(parms->gpu.tomo && parms->tomo.pos !=2){
		parms->gpu.tomo=0;
		warning("\n\nGPU reconstruction is only available for CG with tomo.pos=2 for the moment.\n");
	    }
	    if(parms->gpu.tomo && parms->tomo.alg >2){
		parms->gpu.tomo=0;
		warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Tomography.\n");
	    }
	    if(parms->gpu.fit && parms->fit.alg > 2){
		warning("\n\nGPU reconstruction is only available for CBS/CG. Disable GPU Fitting.\n");
		parms->gpu.fit=0;
	    }
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
	    int moao_used=0;
	    if(parms->evl.moao>=0) moao_used++;
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].moao>=0) moao_used++;
	    }
	    if(moao_used>0){
		if(parms->gpu.moao){
		    if(!parms->fit.square){
			warning("GPU moao=1 requires fit.square=1. Changed\n");
			parms->fit.square=1;
		    }
		    if(!parms->tomo.square){
			warning("GPU moao=1 requires tomo.square=1. Changed\n");
			parms->tomo.square=1;
		    }
		}
	    }
	}else if(parms->recon.alg==1){
	    parms->gpu.tomo=0;
	    parms->gpu.fit=0;
	    if(parms->gpu.lsr){
		warning("\n\nGPU reconstruction for LSR is not available yet\n");
		parms->gpu.lsr=0;
	    }
	}
	if(!parms->recon.glao){
	    int has_dfr=0;
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs<=1) continue;
		if(parms->powfs[ipowfs].dfrs){
		    has_dfr=1;
		}
	    }
	    if(has_dfr){
		warning("\n\nGPU reocnstruction is not yet available for differential focus removal\n\n");
		parms->gpu.tomo=0;
	    }
	}

	if(!parms->atm.frozenflow){
	    warning("Atm is not frozen flow. Disable gpu.evl and gpu.wfs\n");
	    parms->gpu.evl=0;
	    parms->gpu.wfs=0;
	}
	if(parms->gpu.evl && parms->gpu.wfs){
	    parms->sim.cachedm=0; /*Done in CUDA. */
	}
	if(parms->gpu.tomo){
	    parms->tomo.square=1;
	    parms->dbg.dxonedge=1;
	    /*parms->dbg.splitlrt=0;*//*need extensity comparison. */
	}
    }else{
	memset(&(parms->gpu), 0, sizeof(GPU_CFG_T));
    }
    if(parms->dbg.cmpgpu){
	warning("Make cpu code follows gpu implementations.\n");
	parms->sim.cachedm=0;
	parms->tomo.square=1;
	parms->dbg.dxonedge=1;
	/*parms->dbg.splitlrt=0;*//*need extensity comparison. */
    }
    /*Assign each turbulence layer to a corresponding reconstructon layer. Used
      to compute opdx in a simple minded way.*/
    parms->atm.ipsr=calloc(parms->atm.nps, sizeof(int));
    for(int ips=0; ips<parms->atm.nps; ips++){
	double dist=INFINITY;
	int kpsr=-1;
	double ht=parms->atm.ht[ips];
	for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
	    double htr=parms->atmr.ht[ipsr];
	    double dist2=fabs(ht-htr);
	    if(dist2<dist){
		dist=dist2;
		kpsr=ipsr;
	    }
	}
	parms->atm.ipsr[ips]=kpsr;
	/*info("atm layer %d is maped to atmr %d\n", ips,kpsr); */
    }
    /* Map reconstructed layers to input layers. for testing tomo.predict*/
    parms->atmr.indps=calloc(parms->atmr.nps, sizeof(int));
    for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
	parms->atmr.indps[ipsr]=-1;
	for(int ips=0; ips<parms->atm.nps; ips++){
	    if(fabs(parms->atmr.ht[ipsr]-parms->atm.ht[ips])<1e-3){
		if(parms->atmr.indps[ipsr]>-1){
		    warning("One ipsr is mapped to multiple ips\n");
		}
		parms->atmr.indps[ipsr]=ips;
	    }
	}
    }
    if(!parms->sim.closeloop){
	warning2("psfisim is set from %d to %d in openloop mode\n", parms->evl.psfisim, parms->sim.start);
	parms->evl.psfisim=parms->sim.start;
    }
    if(parms->evl.psfisim<parms->sim.start){
	parms->evl.psfisim=parms->sim.start;
    }
    if(NCPU==1 || parms->sim.closeloop==0 || parms->evl.tomo || parms->sim.nthread==1){
	/*disable parallelizing the big loop. */
	parms->sim.parallel=0;
    }
    if(parms->evl.psfmean || parms->evl.psfhist){
	int fnd=sum_intarr(parms->evl.nevl, parms->evl.psf);
	if(fnd==0){
	    error("Required to output PSF, but evl.psf are all zero\n");
	}else{
	    info2("Output PSF for %d directions\n", fnd);
	}
    }
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	parms->evl.npsf+=(parms->evl.psf[ievl]>0);
	if(!parms->recon.split){
	    parms->evl.psfngsr[ievl]=0;
	}
	if(isfinite(parms->evl.hs[ievl]) && parms->evl.psfngsr[ievl]){
	    parms->evl.psfngsr[ievl]=0;
	    if(parms->evl.psfmean || parms->evl.psfhist || parms->evl.opdcov){
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
    double wh=0;
    double wv=0;
    for(int ips=0; ips<parms->atm.nps; ips++){
	wh+=pow(fabs(parms->atm.ht[ips]),5./3.)*parms->atm.wt[ips];
	wv+=pow(fabs(parms->atm.ws[ips]),5./3.)*parms->atm.wt[ips];
    }
    double theta0z=0.3144*parms->atm.r0z*pow(wh,-3./5.);
    double fgreen=0.426/parms->atm.r0z*pow(wv,3./5.);
    double theta2z=0;
    
    info2("\033[0;32mTurbulence at zenith:\033[0;0m\n"
	  "Fried parameter r0 is %gm, Outer scale is %gm Greenwood freq is %.1fHz\n"
	  "Anisoplanatic angle is %.2f\"",
	  parms->atm.r0, parms->atm.l0, fgreen, theta0z*206265);
    if(parms->ndm==2){
	double wf=0;
	double H1=parms->dm[0].ht;
	double H2=parms->dm[1].ht;
	double HH=pow(H2-H1,5./3.);
	for(int ips=0; ips<parms->atm.nps; ips++){
	    double ht=parms->atm.ht[ips];
	    double t1=0.5*pow(fabs(ht-H1),5./3.)+0.5*pow(fabs(ht-H2),5./3.);
	    double t2=-0.25*HH-0.25/HH*pow(pow(fabs(ht-H1),5./3.)-pow(fabs(ht-H2),5./3.),2);
	    wf+=parms->atm.wt[ips]*(t1+t2);
	}
	theta2z=0.3144*parms->atm.r0z*pow(wf,-3./5.);
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
	      ips,parms->atm.ht[ips],parms->atm.wt[ips],parms->atm.ws[ips]);
    }
    if(parms->recon.alg==0){
	info2("\033[0;32mTomography\033[0;0m: r0=%gm l0=%gm "
	      "ZA is %g deg. %d layers.%s\n", 
	      parms->atmr.r0, parms->atmr.l0,  
	      parms->sim.za*180/M_PI, 
	      parms->atmr.nps,(parms->tomo.cone?" use cone coordinate.":""));
  
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    info2("layer %d: ht= %6.0f m, wt= %5.3f\n",
		  ips,parms->atmr.ht[ips],parms->atmr.wt[ips]);
	}
    }
    info2("\033[0;32mThere are %d powfs\033[0;0m\n", parms->npowfs);
    for(i=0; i<parms->npowfs; i++){
	info2("powfs %d: Order %2d, %sGS at %3.3g km. Thres %g%%",
	      i,parms->powfs[i].order, (parms->powfs[i].llt?"L":"N"),
	      parms->powfs[i].hs/1000,parms->powfs[i].saat*100);
	if(parms->powfs[i].trs){
	    info2("\033[0;32m Tip/tilt is removed.\033[0;0m");
	}
	if(parms->powfs[i].dfrs){
	    info2("\033[0;32m Diff focus is removed.\033[0;0m");
	    if(!parms->powfs[i].llt){
		warning("\n\ndfrs=1, but this powfs doesn't have LLT!\n\n");
	    }
	}
	if(parms->powfs[i].pixblur>1.e-12){
	    info2("\033[0;32m Pixel is blurred by %g.\033[0;0m", parms->powfs[i].pixblur);
	}
	info2("\n");
	info2("    CCD image is %dx%d @ %gx%gmas, %gHz, ", 
	      (parms->powfs[i].radpix?parms->powfs[i].radpix:parms->powfs[i].pixpsa), 
	      parms->powfs[i].pixpsa, 
	      parms->powfs[i].radpixtheta*206265000,parms->powfs[i].pixtheta*206265000,
	      1./parms->sim.dt/parms->powfs[i].dtrat);
	info2("wvl: [");
	for(int iwvl=0; iwvl<parms->powfs[i].nwvl; iwvl++){
	    info2(" %g",parms->powfs[i].wvl[iwvl]);
	}
	info2("]\n");
	info2("    %s in reconstruction. ", 
	      parms->powfs[i].gtype_recon==0?"Gtilt":"Ztilt");
	if(parms->powfs[i].phystep>-1){
	    info2("Physical optics start at %d with '%s' %s",
		  parms->powfs[i].phystep, 
		  phytype[parms->powfs[i].phytypesim-1],
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
	info2("wvl: [");
	for(int iwvl=0; iwvl<parms->powfs[i].nwvl; iwvl++){
	    info2(" %g",parms->powfs[i].wvl[iwvl]);
	}
	info2("]\n");
	
    }
    info2("\033[0;32mThere are %d wfs\033[0;0m\n", parms->nwfs);
    for(i=0; i<parms->nwfs; i++){
	info2("wfs %d: type is %d, at (%7.2f, %7.2f) arcsec\n",
	      i,parms->wfs[i].powfs,parms->wfs[i].thetax*206265,
	      parms->wfs[i].thetay*206265);
	if(fabs(parms->wfs[i].thetax)>1 || fabs(parms->wfs[i].thetay)>1){
	    error("wfs thetax or thetay is too large\n");
	}
    }
    info2("\033[0;32mThere are %d DMs\033[0;0m\n",parms->ndm);
    for(i=0; i<parms->ndm; i++){
	info2("DM %d: Order %d, at %4gkm, actuator pitch %gm, offset %3g, with %f micron stroke.\n",
	      i, parms->dm[i].order,
	      parms->dm[i].ht/1000, parms->dm[i].dx,
	      parms->dm[i].offset, 
	      parms->dm[i].stroke*1e6);
	if(parms->dm[i].cubic){
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
	switch(parms->tomo.alg){
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
	info2("Fit %d: Fit  wt is %5.3f, at (%7.2f, %7.2f) arcsec\n",
	      i,parms->fit.wt[i],parms->fit.thetax[i]*206265, 
	      parms->fit.thetay[i]*206265);
	if(fabs(parms->fit.thetax[i])>1 || fabs(parms->fit.thetay[i])>1){
	    error("fit thetax or thetay is too large\n");
	}
    }
    info2("\033[0;32mThere are %d evaluation directions\033[0;0m\n", parms->evl.nevl);
    for(i=0; i<parms->evl.nevl; i++){
	info2("Evl %d: Eval wt is %5.3f, at (%7.2f, %7.2f) arcsec\n",
	      i,parms->evl.wt[i],parms->evl.thetax[i]*206265, 
	      parms->evl.thetay[i]*206265);
	if(fabs(parms->evl.thetax[i])>1 || fabs(parms->evl.thetay[i])>1){
	    error("evl thetax or thetay is too large\n");
	}
    }
    if(parms->plot.setup){
	plotdir("FoV",parms,parms->sim.fov,"fov");/*plot wfs/evaluation direction */
    }
}
/**
   Limited sanity check of the parameters to prevent obvious mistakes.
*/
static void check_parms(const PARMS_T *parms){
    int i;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].llt){
	    if(!isfinite(parms->powfs[ipowfs].hs)){
		warning("powfs with llt should have finite hs\n");
	    }
	    if(parms->powfs[ipowfs].trs==0){
		warning("powfs with llt should have trs=1. Will disable uplink ray tracing, but keep ETF.\n");
	    }
	}else{
	    if(isfinite(parms->powfs[ipowfs].hs)){
		warning("powfs without llt should infinite hs\n");
	    }
	    if(parms->powfs[ipowfs].trs==1){
		warning("powfs without llt should have trs=0\n");
	    }
	}
    }

    for(i=0; i<parms->ndm; i++){
	if(isfinite(parms->dm[i].stroke)){
	    double strokemicron=parms->dm[i].stroke*1e6;
	    if(strokemicron<1 || strokemicron>50){
		warning("dm %d: stroke %g m is probably wrong\n",
			i,parms->dm[i].stroke);
	    }
	}
    }
 
    if(parms->tomo.alg<0 || parms->tomo.alg>2){
	error("parms->tomo.alg=%d is invalid\n", parms->tomo.alg);
    }
    if(parms->fit.alg<0 || parms->fit.alg>2){
	error("parms->fit.alg=%d is invalid\n", parms->tomo.alg);
    }
    if(parms->lsr.alg<0 || parms->lsr.alg>2){
	error("parms->fit.alg=%d is invalid\n", parms->tomo.alg);
    }
}

/**
   This routine calles other routines in this file to setup the parms parameter
   struct parms and check for possible errors. parms is kept constant after
   returned from setup_parms. */
PARMS_T * setup_parms(ARG_T *arg){
    info2("Main config file is %s\n",arg->conf);
    open_config(arg->conf,NULL,0);/*main .conf file. */
    open_config(arg->confcmd, NULL, 1);
    remove(arg->confcmd);
    free(arg->conf);
    free(arg->confcmd);
    PARMS_T* parms=calloc(1, sizeof(PARMS_T));
    parms->sim.nthread=arg->nthread;
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
    readcfg_sim(parms);
    readcfg_cn2(parms);
    readcfg_plot(parms);
    readcfg_dbg(parms);
    readcfg_gpu(parms);
    readcfg_save(parms);
    readcfg_load(parms);
    parms->nsurf=readcfg_strarr(&parms->surf, "surf");
    parms->ntsurf=readcfg_strarr(&parms->tsurf,"tsurf");
    /*
      Output all the readed parms to a single file that can be used to reproduce
      the same simulation.
    */
    char fnconf[PATH_MAX];
    snprintf(fnconf, PATH_MAX, "maos_%ld.conf", (long)getpid());
    close_config("%s",fnconf);
    mysymlink(fnconf, "maos_recent.conf");
#if USE_CUDA 
    if(parms->nwfs==1 && arg->ngpu==0) arg->ngpu=1;/*use a single gpu is there is only 1 wfs.*/
    use_cuda=gpu_init(arg->gpus, arg->ngpu);
#else
    use_cuda=0;
#endif
    /*
      Postprocess the parameters for integrity. The ordering of the following
      routines are critical.
    */
    setup_parms_postproc_sim(parms);
    setup_parms_postproc_za(parms);
    setup_parms_postproc_wfs(parms);
    setup_parms_postproc_atm(parms);
    setup_parms_postproc_atm_size(parms);
    setup_parms_postproc_dm(parms);
    setup_parms_postproc_recon(parms);
    setup_parms_postproc_siglev(parms);
    setup_parms_postproc_misc(parms, arg);
    check_parms(parms);
    print_parms(parms);
    return parms;
}
