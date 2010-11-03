/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "maos.h"
#include "mtch.h"
/**
   \file maos/mtch.c
   Setting up matched filter
*/
/**
   shift i0 without wraping into i0x1 (+1) and i0x2 (-1)
*/
static void mki0shx(double *i0x1, double *i0x2, dmat *i0, double scale){
    int nx=i0->nx;
    double (*i0p)[nx]=(void*)i0->p;
    double (*i0x1p)[nx]=(void*)i0x1;
    double (*i0x2p)[nx]=(void*)i0x2;
    for(int iy=0; iy<i0->ny; iy++){
	for(int ix=0; ix<i0->nx-1; ix++){
	    i0x1p[iy][ix+1]=i0p[iy][ix]*scale;
	    i0x2p[iy][ix]=i0p[iy][ix+1]*scale;
	}
    }
}

/**
  shift i0 without wraping into i0y1 (+1) and i0y2 (-1)
*/
static void mki0shy(double *i0y1, double *i0y2, dmat *i0, double scale){
    int nx=i0->nx;
    double (*i0p)[nx]=(void*)i0->p;
    double (*i0y1p)[nx]=(void*)i0y1;
    double (*i0y2p)[nx]=(void*)i0y2;
    for(int iy=0; iy<i0->ny-1; iy++){
	for(int ix=0; ix<i0->nx; ix++){
	    i0y1p[iy+1][ix]=i0p[iy][ix]*scale;
	    i0y2p[iy][ix]=i0p[iy+1][ix]*scale;
	}
    }
}
/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const PARMS_T *parms, POWFS_T *powfs, const int ipowfs){
    if(!powfs[ipowfs].intstat){
	error("Please create intstat before calling mtch");
    }
    const double pixtheta=parms->powfs[ipowfs].pixtheta;
    const double kp=1./pixtheta;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrndfn_res=(1.-parms->powfs[ipowfs].bkgrndrm);
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const int sub_i0=1;//doesn't make any difference.
    int ni0=powfs[ipowfs].intstat->i0->ny;
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
    }
    const int nsa=powfs[ipowfs].pts->nsa;
    if(powfs[ipowfs].intstat->mtche){
	warning("powfs %d: matched filter already exists. free them\n",ipowfs);
	dcellfree(powfs[ipowfs].intstat->mtche);
	dcellfree(powfs[ipowfs].intstat->sanea);
	dcellfree(powfs[ipowfs].intstat->saneara);
	dcellfree(powfs[ipowfs].intstat->saneaxy);
	dcellfree(powfs[ipowfs].intstat->saneaixy);
	dfree(powfs[ipowfs].intstat->i0sum);
    }
    powfs[ipowfs].intstat->mtche=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->sanea=dcellnew(ni0,1);
    powfs[ipowfs].intstat->saneara=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->saneaxy=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->saneaixy=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->i0sum=dnew(nsa,ni0);
    PDCELL(powfs[ipowfs].intstat->i0,i0s);
    PDCELL(powfs[ipowfs].intstat->gx,gxs);
    PDCELL(powfs[ipowfs].intstat->gy,gys);
    PDMAT(powfs[ipowfs].intstat->i0sum,i0sum);
    dmat *(*mtche)[nsa]=
	(dmat*(*)[nsa])powfs[ipowfs].intstat->mtche->p;
    dmat *(*saneaxy)[nsa]
	=(void*)powfs[ipowfs].intstat->saneaxy->p;
    dmat *(*saneaixy)[nsa]
	=(void*)powfs[ipowfs].intstat->saneaixy->p;
    dmat *(*saneara)[nsa]=NULL;
    if(parms->powfs[ipowfs].radpix){
	saneara=(void*)powfs[ipowfs].intstat->saneara->p;
    }
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    int irot_multiplier=nllt>1?1:0;

    int nmod=3;
    int mtchcrx=0;
    int mtchcry=0;
    double shiftx=0, shifty=0;
    //always use saved i0. cyclic shift is not good
    //because of the wrapped ring.
    if(fabs(parms->powfs[ipowfs].mtchcrx)>1.e-10){
	shiftx=parms->powfs[ipowfs].mtchcrx;
	mtchcrx=nmod;
	nmod+=2;
	if(fabs(shiftx-1)>1.e-10){
	    error("Only constraint of 1 pixel is implemented\n");
	}
    }
    if(fabs(parms->powfs[ipowfs].mtchcry)>1.e-10){
	shifty=parms->powfs[ipowfs].mtchcry;
	mtchcry=nmod;
	nmod+=2;
	if(fabs(shifty-1)>1.e-10){
	    error("Only constraint of 1 pixel is implemented\n");
	}
    }
    const int i0n=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;

    dmat *i0m=dnew(2,nmod);
    dmat *i0g=dnew(i0n,nmod);
    PDMAT(i0g, pi0g);
    PDMAT(i0m, pi0m);
    dmat *i0x1=NULL, *i0x2=NULL, *i0y1=NULL, *i0y2=NULL;
    dmat *wt=dnew(i0n,1);
    for(int ii0=0; ii0<ni0; ii0++){
	int iwfs=parms->powfs[ipowfs].wfs[ii0];
	double *srot=NULL;
	if(powfs[ipowfs].srot){
	    int irot=ii0*irot_multiplier;
	    srot=powfs[ipowfs].srot->p[irot]->p;
	}
	powfs[ipowfs].intstat->sanea->p[ii0]=dnew(nsa,2);
	PDMAT(powfs[ipowfs].intstat->sanea->p[ii0], psanea);
	pi0m[0][0]=1;
	pi0m[1][1]=1;

	if(mtchcrx){//constrained x(radial)
	    double shift=parms->powfs[ipowfs].pixtheta*shiftx;
	    pi0m[mtchcrx][0]=shift*kp;//k is here to ensure good conditioning
	    pi0m[mtchcrx+1][0]=-shift*kp;
	}
	if(mtchcry){//constrained y(azimuthal).
	    double shift=parms->powfs[ipowfs].pixtheta*shifty;
	    pi0m[mtchcry][1]=shift*kp;
	    pi0m[mtchcry+1][1]=-shift*kp;
	}

	for(int isa=0; isa<nsa; isa++){
	    i0sum[ii0][isa]=dsum(i0s[ii0][isa]);
	    double siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
	    if(i0sum[ii0][isa]<siglev*0.05 || i0sum[ii0][isa]>siglev){
		warning("powfs %d: i0 sum to %g, but wfs %d has siglev of %g\n",
			ipowfs, i0sum[ii0][isa], iwfs, siglev);
	    }
	    double* bkgrnd2=NULL;
	    if(powfs[ipowfs].bkgrnd && powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]){
		bkgrnd2= powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]->p; 
	    }
	    dzero(i0g);//don't forget to zero out
	    adddbl(pi0g[0], 1,gxs[ii0][isa]->p, i0n, 1);
	    adddbl(pi0g[1], 1,gys[ii0][isa]->p, i0n, 1);
	    adddbl(pi0g[2], 1,i0s[ii0][isa]->p, i0n, kp);
	    if(powfs[ipowfs].bkgrnd && powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]){
		//notice that the bkgrnd is already scaled by sim.dt and dtrat properly.
		adddbl(pi0g[2], 1, bkgrnd2, i0n, bkgrndfn_res);
	    }
	    if(mtchcrx){
		/*
		  constrained matched filter. compute separately for each wfs.
		*/
		mki0shx(pi0g[mtchcrx],pi0g[mtchcrx+1],i0s[ii0][isa],kp);
		if(sub_i0){
		    adddbl(pi0g[mtchcrx],1,i0s[ii0][isa]->p, i0n, -kp);
		    adddbl(pi0g[mtchcrx+1],1,i0s[ii0][isa]->p, i0n, -kp);
		}
		if(bkgrnd2){
		    adddbl(pi0g[mtchcrx], 1, bkgrnd2, i0n, bkgrndfn_res);
		    adddbl(pi0g[mtchcrx+1],1,bkgrnd2,i0n,  bkgrndfn_res);
		}
	    }
	    if(mtchcry){
		mki0shy(pi0g[mtchcry],pi0g[mtchcry+1],i0s[ii0][isa],kp);
		if(sub_i0){
		    adddbl(pi0g[mtchcry],1,i0s[ii0][isa]->p, i0n, -kp);
		    adddbl(pi0g[mtchcry+1],1,i0s[ii0][isa]->p,i0n, -kp);
		}
		 if(bkgrnd2){ 
		     adddbl(pi0g[mtchcry], 1, bkgrnd2, i0n, bkgrndfn_res);
		     adddbl(pi0g[mtchcry+1],1, bkgrnd2,i0n,  bkgrndfn_res);
		}
	    }
	  
	    if(bkgrnd2){
		//adding rayleigh backscatter poisson noise.
		for(int i=0; i<i0n; i++){//noise weighting.
		    wt->p[i]=1./(rne*rne+bkgrnd+i0s[ii0][isa]->p[i]+bkgrnd2[i]);
		}	
	    }else{
		for(int i=0; i<i0n; i++){//noise weighting.
		    wt->p[i]=1./(rne*rne+bkgrnd+i0s[ii0][isa]->p[i]);
		}
	    }
	    dmat *tmp=dpinv(i0g, wt, NULL);
	    

	    dmm(&mtche[ii0][isa],i0m, tmp, "nn", 1);
	    dfree(tmp);
	    for(int i=0; i<i0n; i++){//noise weighting.
		wt->p[i]=1./wt->p[i];
	    }
	    dmat *nea2=dtmcc(mtche[ii0][isa], wt);

	    if(parms->powfs[ipowfs].mtchcpl==0){
		//remove coupling between r/a measurements.
		nea2->p[1]=nea2->p[2]=0;
	    }
	    psanea[0][isa]=nea2->p[0];
	    psanea[1][isa]=nea2->p[3];
		
	    if(parms->powfs[ipowfs].radpix){
		saneara[ii0][isa]=nea2;
		double theta=srot[isa]; 
		drotvecnn(&saneaxy[ii0][isa], saneara[ii0][isa], theta);
	    }else{
		saneaxy[ii0][isa]=nea2;
	    }
	    saneaixy[ii0][isa]=dinvspd(saneaxy[ii0][isa]);
	}//isa 
    }//ii0
    info2("powfs %d: matched filter sanea:\n",ipowfs);
    if(powfs[ipowfs].sprint){
	if(parms->powfs[ipowfs].llt->n!=ni0){
	    warning("nllt!=ni0\n");
	}
	for(int ii0=0; ii0<ni0; ii0++){
	    int illt=0;
	    if(ni0==parms->powfs[ipowfs].llt->n){
		illt=ii0;
	    }else if(ni0==parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].llt->n==1){
		illt=0;
	    }else{
		error("Invalid combination\n");
	    }
	    info2("llt %d:\n",illt);
	    PDMAT(powfs[ipowfs].intstat->sanea->p[ii0], psanea);
	    for(int ksa=0; ksa<powfs[ipowfs].sprint->p[illt]->nx; ksa++){
		int isa=(int)powfs[ipowfs].sprint->p[illt]->p[ksa];
		if(isa>0){
		    info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			  isa, powfs[ipowfs].srsa->p[illt]->p[isa], 
			  sqrt(psanea[0][isa])*206265000,
			  sqrt(psanea[1][isa])*206265000);
		}
	    }
	}
    }else{
	for(int ii0=0; ii0<ni0; ii0++){
	    info("ii0=%d:\n",ii0);
	    PDMAT(powfs[ipowfs].intstat->sanea->p[ii0], psanea);
	    double dsa=powfs[ipowfs].saloc->dx;
	    double llimit=-dsa/2;
	    double ulimit=dsa/2;
	    for(int isa=0; isa<nsa; isa++){
		double locx=powfs[ipowfs].saloc->locx[isa];
		double locy=powfs[ipowfs].saloc->locy[isa];
		if(nsa<10 || (locx>0&&locy>llimit&&locy<ulimit)){
		    info2("sa %5d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			  isa, locx, sqrt(psanea[0][isa])*206265000,
			  sqrt(psanea[1][isa])*206265000);
		}
	    }//isa 
	}//ii0
    }
    if(parms->save.setup){
	dcellwrite(powfs[ipowfs].intstat->mtche,
		   "%s/powfs%d_mtche",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].intstat->saneara,
		   "%s/powfs%d_saneara",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].intstat->saneaxy,
		   "%s/powfs%d_saneaxy",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].intstat->saneaixy,
		   "%s/powfs%d_saneaixy",dirsetup,ipowfs);
	dcellwrite(powfs[ipowfs].intstat->sanea,
		   "%s/powfs%d_sanea",dirsetup,ipowfs);
    }
    dfree(i0m);
    dfree(i0g);
    dfree(i0x1); dfree(i0x2); dfree(i0y1); dfree(i0y2);
    dfree(wt);
}
