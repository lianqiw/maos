/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "mtch.h"
/**
   \file maos/mtch.c
   Setting up matched filter
*/

/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const PARMS_T *parms, POWFS_T *powfs, const int ipowfs){
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    const real pixthetax=parms->powfs[ipowfs].radpixtheta;
    const real pixthetay=parms->powfs[ipowfs].pixtheta;
    const real rne=parms->powfs[ipowfs].rne;
    const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
    int ni0=intstat->i0->ny;
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
    }
    const int nsa=powfs[ipowfs].saloc->nloc;
    //Prevent printing of NEA during recomputing of matched filter
    const int print_nea=intstat->mtche?0:1;
    dcellfree(intstat->mtche);
    dfree(intstat->i0sum);

    intstat->mtche=dcellnew(nsa,ni0);
    dcellfree(powfs[ipowfs].sanea);
    dcell *sanea=powfs[ipowfs].sanea=dcellnew(ni0,1);
    intstat->i0sum=dnew(nsa,ni0);
    intstat->i0sumsum=dnew(ni0, 1);

    dcell *i0s=intstat->i0;
    dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx/*PDELL*/;
    dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy/*PDELL*/;
    dmat *i0sum=intstat->i0sum;
    dcell *mtche=intstat->mtche;
  
    //dcell *saneaxy=powfs[ipowfs].saneaxy;
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    int irot_multiplier=nllt>1?1:0;
    const int mtchadp=parms->powfs[ipowfs].mtchadp;
  
    for(int ii0=0; ii0<ni0; ii0++){
	int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	const real siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
	real i0thres=MAX(0.1*siglev, rne*10);
	real nea2thres=pixthetax*pixthetay*100;
	real *srot=NULL;
	if(powfs[ipowfs].srot){
	    int irot=ii0*irot_multiplier;
	    srot=powfs[ipowfs].srot->p[irot]->p;
	}
	sanea->p[ii0]=dnew(nsa,3);
	dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
	real i0sumsum=0;
	int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
	int ncrdisable=0;
	const int radgx=parms->powfs[ipowfs].radgx;
	dmat *nea2=0;
	for(int isa=0; isa<nsa; isa++){
	    real pixrot=0;//pixel rotation
	    if(srot && parms->powfs[ipowfs].radpix){
		pixrot=srot[isa]; 
	    }
	    if(mtchadp){
		long fwhm=dfwhm(P(i0s,isa,ii0));
		if(fwhm>4){
		    crdisable=0;
		}else{
		    crdisable=1;
		    ncrdisable++;
		}
	    }
	    dmat* bkgrnd2=NULL;
	    dmat* bkgrnd2c=NULL;
	    if(powfs[ipowfs].bkgrnd){
		bkgrnd2= powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]; 
	    }
	    if(powfs[ipowfs].bkgrndc){
		bkgrnd2c= powfs[ipowfs].bkgrndc->p[ii0*nsa+isa]; 
	    }
	    P(mtche,isa,ii0)=mtch(&nea2, P(i0s,isa,ii0),
				    gxs?P(gxs,isa,ii0):0, gys?P(gys,isa,ii0):0, 
				    parms->powfs[ipowfs].qe,
				    bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				    pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
	    
	    P(i0sum,isa,ii0)=dsum(P(i0s,isa,ii0));
	    i0sumsum+=P(i0sum,isa,ii0);

	    if(P(i0sum,isa,ii0)<i0thres || nea2->p[0]>nea2thres || nea2->p[3]>nea2thres){
		//Signal level too low or error to high.
		nea2->p[0]=nea2->p[3]=nea2thres;
		nea2->p[1]=nea2->p[2]=0;
		dset(P(mtche,isa,ii0), 0);
	    }
	    if(parms->powfs[ipowfs].mtchcpl==0 
	       && (!parms->powfs[ipowfs].radpix || parms->powfs[ipowfs].radgx)){
		/*remove coupling between r/a (x/y) measurements. */
		nea2->p[1]=nea2->p[2]=0;
	    }
	    P(psanea,isa,0)=nea2->p[0];
	    P(psanea,isa,1)=nea2->p[3];
	    P(psanea,isa,2)=nea2->p[1];
	}/*isa  */
	dfree(nea2);

	if(mtchadp){
	    info("Mtched filter contraint are disabled for %d subaps out of %d.\n",
		  ncrdisable, nsa);
	}
	intstat->i0sumsum->p[ii0]=i0sumsum;
    }/*ii0 */
    if(print_nea){
	info("Matched filter sanea:\n");
	if(powfs[ipowfs].sprint){/*print nea for select subapertures.*/
	    for(int ii0=0; ii0<ni0; ii0++){
		int illt=0;
		if(ni0==parms->powfs[ipowfs].llt->n){
		    illt=ii0;
		}else if(ni0==parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].llt->n==1){
		    illt=0;
		}else{
		    error("Invalid combination\n");
		}
		info("ii0 %d, llt %d.\n", ii0, illt);
		info("sa index   dist   noise equivalent angle\n");
		dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
		for(int ksa=0; ksa<powfs[ipowfs].sprint->p[illt]->nx; ksa++){
		    int isa=(int)powfs[ipowfs].sprint->p[illt]->p[ksa];
		    if(isa>0){
			info("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			      isa, powfs[ipowfs].srsa->p[illt]->p[isa], 
			      sqrt(P(psanea,isa,0))*206265000,
			      sqrt(P(psanea,isa,1))*206265000);
		    }
		}
	    }
	}else{
	    real dsa=powfs[ipowfs].saloc->dx;
	    real llimit=-dsa/2;
	    real ulimit=dsa/2;
	    info("index: position noise equivalent angle\n");
	    for(int isa=0; isa<nsa; isa++){
		real locx=powfs[ipowfs].saloc->locx[isa];
		real locy=powfs[ipowfs].saloc->locy[isa];
		if((parms->powfs[ipowfs].llt && (nsa<10 || (locx>0&&locy>llimit&&locy<ulimit)))
		   ||(!parms->powfs[ipowfs].llt && locx>=0 && locx<dsa*0.6 && locy>=0 && locy<dsa*0.6)
		    || nsa<=4){
		    info("sa%4d:%6.1fm",isa, locx);
		    for(int ii0=0; ii0<ni0; ii0++){
			info(" (%4.1f,%4.1f)", 
			      sqrt(P(sanea->p[ii0],isa,0))*206265000,
			      sqrt(P(sanea->p[ii0],isa,1))*206265000);
		    }//for ii0
		    info(" mas\n");
		}
	    }/*isa  */
	}
    }
    if(parms->save.setup){
	writebin(sanea, "powfs%d_sanea", ipowfs);
    }
    if(parms->powfs[ipowfs].phytype_recon==1 && parms->recon.glao && ni0>0){
	info("Averaging saneaxy of different WFS for GLAO mode\n");
	dmat *sanea2=0;
	real scale=1./ni0;
	for(int ii0=0; ii0<ni0; ii0++){
	    dadd(&sanea2, 1, sanea->p[ii0], scale);
	}
	dcellfree(powfs[ipowfs].sanea);
	powfs[ipowfs].sanea=dcellnew(1,1);
	powfs[ipowfs].sanea->p[0]=sanea2;
    }
}
