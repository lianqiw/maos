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
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const double bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;

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
    dcell *sanea=dcellnew(ni0,1);
    
    intstat->i0sum=dnew(nsa,ni0);
    intstat->i0sumsum=dnew(ni0, 1);

    dcell *i0s=intstat->i0;
    dcell* gxs=intstat->gx/*PDELL*/;
    dcell* gys=intstat->gy/*PDELL*/;
    dmat *i0sum=intstat->i0sum;
    dcell *mtche=intstat->mtche;
    if(parms->powfs[ipowfs].phytype==1){//use MF nea for recon
	dcellfree(powfs[ipowfs].saneaxy);
	powfs[ipowfs].saneaxy=dcellnew(nsa,ni0);
    }
    dcell *saneaxy=powfs[ipowfs].saneaxy;
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    int irot_multiplier=nllt>1?1:0;
    const int mtchadp=parms->powfs[ipowfs].mtchadp;
    double neaspeckle=parms->powfs[ipowfs].neaspeckle/206265000.;
    if(neaspeckle>pixthetax){
	error("parms->powfs[%d].neaspeckle=%g is bigger than pixel size\n",
	      ipowfs, neaspeckle);
    }
    if(neaspeckle>0){
	warning2("powfs%d: Adding speckle noise of %.2f mas\n", ipowfs, neaspeckle*206265000);
    }
    double neaspeckle2=pow(neaspeckle,2);
    for(int ii0=0; ii0<ni0; ii0++){
	int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	double *srot=NULL;
	if(powfs[ipowfs].srot){
	    int irot=ii0*irot_multiplier;
	    srot=powfs[ipowfs].srot->p[irot]->p;
	}
	sanea->p[ii0]=dnew(nsa,2);
	dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
	double i0summax=0;
	double i0sumsum=0;
	int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
	int ncrdisable=0;
	const int radgx=parms->powfs[ipowfs].radgx;
	for(int isa=0; isa<nsa; isa++){
	    double pixrot=0;//pixel rotation
	    if(srot && parms->powfs[ipowfs].radpix){
		pixrot=srot[isa]; 
	    }
	    if(mtchadp){
		long fwhm=dfwhm(IND(i0s,isa,ii0));
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
	    dmat *nea2=0;
	    IND(mtche,isa,ii0)=mtch(&nea2, IND(i0s,isa,ii0), IND(gxs,isa,ii0), IND(gys,isa,ii0), 
				    parms->powfs[ipowfs].qe,
				    bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				    pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
	    
	    IND(i0sum,isa,ii0)=dsum(IND(i0s,isa,ii0));
	    i0sumsum+=IND(i0sum,isa,ii0);
	    if(IND(i0sum,isa,ii0)>i0summax){
		i0summax=IND(i0sum,isa,ii0);
	    }

	    nea2->p[0]+=neaspeckle2;
	    nea2->p[3]+=neaspeckle2;
	    if(IND(i0sum,isa,ii0)<EPS){//zero flux
		nea2->p[0]=nea2->p[3]=pixthetax*10;
	    }
	    if(parms->powfs[ipowfs].mtchcpl==0 
	       && (!parms->powfs[ipowfs].radpix || parms->powfs[ipowfs].radgx)){
		/*remove coupling between r/a (x/y) measurements. */
		nea2->p[1]=nea2->p[2]=0;
	    }
	    IND(psanea,isa,0)=nea2->p[0];
	    IND(psanea,isa,1)=nea2->p[3];
	    IND(saneaxy, isa, ii0)=nea2;
	}/*isa  */
	double siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
	if(i0summax<siglev*0.1 || i0summax>siglev*1.1){
	    warning("i0 sum to maximum of %g, wfs %d has siglev of %g\n",
		    i0summax, iwfs, siglev);
	}
	if(mtchadp){
	    info2("Mtched filter contraint are disabled for %d subaps out of %d.\n",
		  ncrdisable, nsa);
	}
	intstat->i0sumsum->p[ii0]=i0sumsum;
    }/*ii0 */
    if(print_nea){
	info2("Matched filter sanea:\n");
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
		info2("ii0 %d, llt %d.\n", ii0, illt);
		info2("sa index   dist   noise equivalent angle\n");
		dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
		for(int ksa=0; ksa<powfs[ipowfs].sprint->p[illt]->nx; ksa++){
		    int isa=(int)powfs[ipowfs].sprint->p[illt]->p[ksa];
		    if(isa>0){
			info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			      isa, powfs[ipowfs].srsa->p[illt]->p[isa], 
			      sqrt(IND(psanea,isa,0))*206265000,
			      sqrt(IND(psanea,isa,1))*206265000);
		    }
		}
	    }
	}else{
	    double dsa=powfs[ipowfs].saloc->dx;
	    double llimit=-dsa/2;
	    double ulimit=dsa/2;
	    info2("sa index: radius   noise equivalent angle\n");
	    for(int isa=0; isa<nsa; isa++){
		double locx=powfs[ipowfs].saloc->locx[isa];
		double locy=powfs[ipowfs].saloc->locy[isa];
		if((parms->powfs[ipowfs].llt && (nsa<10 || (locx>0&&locy>llimit&&locy<ulimit)))
		   ||(!parms->powfs[ipowfs].llt && locx>=0 && locx<dsa*0.6 && locy>=0 && locy<dsa*0.6)
		    ){
		    info2("sa%4d:%4.1fm",isa, locx);
		    for(int ii0=0; ii0<ni0; ii0++){
			info2(" (%4.1f,%4.1f)", 
			      sqrt(IND(sanea->p[ii0],isa,0))*206265000,
			      sqrt(IND(sanea->p[ii0],isa,1))*206265000);
		    }//for ii0
		    info2("mas\n");
		}
	    }/*isa  */
	}
    }
    if(parms->powfs[ipowfs].phytype==1 && parms->save.setup){
	writebin(sanea, "powfs%d_sanea", ipowfs);
    }
    if(parms->powfs[ipowfs].phytype==1 && parms->recon.glao && ni0>0){
	info2("Averaging saneaxy of different WFS for GLAO mode\n");
	dcell *saneaxy2=dcellnew(nsa, 1);
	double scale=1./ni0;
	for(int isa=0; isa<nsa; isa++){
	    for(int ii0=0; ii0<ni0; ii0++){
		dadd(&saneaxy2->p[isa], 1, IND(saneaxy, isa, ii0), scale);
	    }
	}
	dcellfree(powfs[ipowfs].saneaxy);
	powfs[ipowfs].saneaxy=saneaxy2;
    }
    dcellfree(sanea);
}
