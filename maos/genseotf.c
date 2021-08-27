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

/*
  2009-11-26: changed to rotate OTF instead of psf to comply with the approach
  in wfsints.  this gives slightly larger intensity because max value of OTF is
  preserved which corresponds to the sum of psf.


*/

#include "common.h"
#include "genseotf.h"

/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs){
	intstat_t* intstat=powfs[ipowfs].intstat;
	const real pixthetax=parms->powfs[ipowfs].radpixtheta;
	const real pixthetay=parms->powfs[ipowfs].pixtheta;
	const real rne=parms->powfs[ipowfs].rne;
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
	const int radgx=parms->powfs[ipowfs].radgx;
	int ni0=intstat->i0->ny;

	if(powfs[ipowfs].bkgrnd&&powfs[ipowfs].bkgrnd->ny>ni0){
		error("Please implement the case when bkgrnd has more columns\n");
	}
	info("Generating matched filter for %d\n", ipowfs);
	if(ni0!=1&&ni0!=parms->powfs[ipowfs].nwfs){
		error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
	}
	const int nsa=powfs[ipowfs].saloc->nloc;
	//Prevent printing of NEA during recomputing of matched filter
	//const int print_nea=intstat->mtche?0:1;
		
	dcellfree(powfs[ipowfs].sanea);
	dcell* sanea=powfs[ipowfs].sanea=dcellnew_same(ni0, 1, nsa, 3);
	dfree(intstat->i0sum);
	dmat* i0sum=intstat->i0sum=dnew(nsa, ni0);
	dfree(intstat->i0sumsum);
	intstat->i0sumsum=dnew(ni0, 1);

	const dcell* i0s=intstat->i0;
	const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx;
	const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy;
	
	long npix=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
	dcellfree(intstat->mtche);
	dcell* mtche=intstat->mtche=dcellnew_same(nsa, ni0, 2, npix);

	//dcell *saneaxy=powfs[ipowfs].saneaxy;
	int nllt;
	if(parms->powfs[ipowfs].llt){
		nllt=parms->powfs[ipowfs].llt->n;
	} else{
		nllt=0;
	}
	int irot_multiplier=nllt>1?1:0;
	const int mtchadp=parms->powfs[ipowfs].mtchadp;
	real sigratio=parms->powfs[ipowfs].sigrecon>0?(parms->powfs[ipowfs].sigrecon/parms->powfs[ipowfs].siglev):1;
	real sigratior=1./sigratio;
	if(sigratio<0){
		error("sigratio cannot be negative\n");
	}

	for(int ii0=0; ii0<ni0; ii0++){
		int iwfs=P(parms->powfs[ipowfs].wfs,ii0);
		const real siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
		real i0thres=MAX(0.1*siglev, rne*10);
		real nea2thres=pixthetax*pixthetay*100;
		real* srot=NULL;
		if(powfs[ipowfs].srot){
			int irot=ii0*irot_multiplier;
			srot=P(powfs[ipowfs].srot,irot)->p;
		}
		//P(sanea,ii0)=dnew(nsa, 3);
		dmat* psanea=P(sanea,ii0)/*PDMAT*/;
		real i0sumsum=0;
		int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
		int ncrdisable=0;
		dmat* nea2=0;
		for(int isa=0; isa<nsa; isa++){
			real pixrot=0;//pixel rotation
			if(srot&&parms->powfs[ipowfs].radpix){
				pixrot=srot[isa];
			}
			if(mtchadp){
				long fwhm=dfwhm_gauss(P(i0s, isa, ii0));
				if(fwhm>4){
					crdisable=0;
				} else{
					crdisable=1;
					ncrdisable++;
				}
			}
			dmat* bkgrnd2=NULL;
			dmat* bkgrnd2c=NULL;
			if(powfs[ipowfs].bkgrnd){
				bkgrnd2=P(powfs[ipowfs].bkgrnd, isa, ii0);
			}
			if(powfs[ipowfs].bkgrndc){
				bkgrnd2c=P(powfs[ipowfs].bkgrndc, isa, ii0);
			}
			mtch(&P(mtche, isa, ii0), &nea2, P(i0s, isa, ii0),
				gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
				parms->powfs[ipowfs].qe,
				bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
			if(fabs(sigratio-1)>1e-5){
				dscale(P(i0s, isa, ii0), sigratio);
				if(gxs) dscale(P(gxs, isa, ii0), sigratio);
				if(gys) dscale(P(gys, isa, ii0), sigratio);
				mtch(NULL, &nea2, P(i0s, isa, ii0),
					gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
					parms->powfs[ipowfs].qe,
					bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
					pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
				dscale(P(i0s, isa, ii0), sigratior);
				if(gxs) dscale(P(gxs, isa, ii0), sigratior);
				if(gys) dscale(P(gys, isa, ii0), sigratior);
			}
			P(i0sum, isa, ii0)=dsum(P(i0s, isa, ii0));
			i0sumsum+=P(i0sum, isa, ii0);

			if(P(i0sum, isa, ii0)<i0thres||P(nea2,0)>nea2thres||P(nea2,3)>nea2thres){
			//Signal level too low or error to high.
				P(nea2,0)=P(nea2,3)=nea2thres;
				P(nea2,1)=P(nea2,2)=0;
				dset(P(mtche, isa, ii0), 0);
			}
			if(parms->powfs[ipowfs].mtchcpl==0
				&&(!parms->powfs[ipowfs].radpix||radgx)){
			 /*remove coupling between r/a (x/y) measurements. */
				P(nea2,1)=P(nea2,2)=0;
			}
			P(psanea, isa, 0)=P(nea2,0);
			P(psanea, isa, 1)=P(nea2,3);
			P(psanea, isa, 2)=P(nea2,1);
		}/*isa  */
		dfree(nea2);

		if(mtchadp){
			info("Mtched filter contraint are disabled for %d subaps out of %d.\n",
				ncrdisable, nsa);
		}
		P(intstat->i0sumsum,ii0)=i0sumsum;
	}/*ii0 */
	if(1 /*print_nea*/){
		info2("Matched filter sanea:\n");
		if(powfs[ipowfs].sprint){/*print nea for select subapertures.*/
			for(int ii0=0; ii0<ni0; ii0++){
				int illt=0;
				if(ni0==parms->powfs[ipowfs].llt->n){
					illt=ii0;
				} else if(ni0==parms->powfs[ipowfs].nwfs&&parms->powfs[ipowfs].llt->n==1){
					illt=0;
				} else{
					error("Invalid combination\n");
				}
				info2("ii0 %d, llt %d.\n", ii0, illt);
				info2("sa index   dist   noise equivalent angle\n");
				dmat* psanea=P(sanea,ii0)/*PDMAT*/;
				for(int ksa=0; ksa<P(powfs[ipowfs].sprint,illt)->nx; ksa++){
					int isa=(int)P(P(powfs[ipowfs].sprint,illt),ksa);
					if(isa>0){
						info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n",
							isa, P(P(powfs[ipowfs].srsa,illt),isa),
							sqrt(P(psanea, isa, 0))*206265000,
							sqrt(P(psanea, isa, 1))*206265000);
					}
				}
			}
		} else{
			real dsa=powfs[ipowfs].saloc->dx;
			real llimit=-dsa/2;
			real ulimit=dsa/2;
			info2("index: position noise equivalent angle\n");
			for(int isa=0; isa<nsa; isa++){
				real locx=powfs[ipowfs].saloc->locx[isa];
				real locy=powfs[ipowfs].saloc->locy[isa];
				if((parms->powfs[ipowfs].llt&&(nsa<10||(locx>0&&locy>llimit&&locy<ulimit)))
					||(!parms->powfs[ipowfs].llt&&locx>=0&&locx<dsa*0.6&&locy>=0&&locy<dsa*0.6)
					||nsa<=4){
					info2("sa%4d:%6.1fm", isa, locx);
					for(int ii0=0; ii0<ni0; ii0++){
						info2(" (%4.1f,%4.1f)",
							sqrt(P(P(sanea,ii0), isa, 0))*206265000,
							sqrt(P(P(sanea,ii0), isa, 1))*206265000);
					}//for ii0
					info2(" mas\n");
				}
			}/*isa  */
		}
	}
}
