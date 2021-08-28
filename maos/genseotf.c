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
void mtch_wrap(
	dcell** pmtche,  /**<[out] matched filter. */
	dcell** psanea,  /**<[out] subaperture noise equivalent angle*/
	dmat** pi0sum,   /**<[out] sum of subaperture intensity*/
	dmat** pi0sumsum,/**<[out] sum of all subaperture intensity*/
	const dcell* i0s, /**<Subaperture intensities*/
	const dcell* gxs, /**<Subaperture intensity gradients along x (radial)*/
	const dcell* gys, /**<Subaperture intensity gradients along y (azimuthal)*/
	real pixthetax, /**<Pixel size along x (radial)*/
	real pixthetay, /**<Pixel size along y (azimuthal)*/
	real rne,     /**<Read out noise per pixel*/
	const dmat* qe, /**<Quantum efficiency of each pixel*/
	const dmat* siglevs, /**<Signal levels at dtrat*/
	real bkgrnd,  /**<bkgrnd per pixel*/
	real bkgrndc, /**<bkgrnd correction per pixel*/
	const dcell* bkgrnd2, /**<bkgrnd image*/
	const dcell* bkgrnd2c,/**<bkgrnd correction image*/
	real sigratio, /**<scale signal level to increase NEA*/
	dcell* srot,   /**<subaperture rotation*/
	int radpix,    /**<Radial coordinate pixel*/
	int radgx,     /**<Leave gradients at radial direction */
	int mtchcr     /**<constraint. -1: auto*/
	){
	int ni0=NY(i0s);
	const int nsa=NX(i0s);
	if(psanea){
		dcellfree(*psanea);
		*psanea=dcellnew_same(ni0, 1, nsa, 3);
	}
	if(pi0sum){
		dfree(*pi0sum);
		*pi0sum=dnew(nsa, ni0);
	}
	if(pi0sumsum){
		dfree(*pi0sumsum);
		*pi0sumsum=dnew(ni0, 1);
	}

	const int npix=PN(P(i0s, 0, 0));
	if(pmtche){
		dcellfree(*pmtche);
		*pmtche=dcellnew_same(nsa, ni0, 2, npix);
	}
	real sigratior=1./sigratio;
	for(int ii0=0; ii0<ni0; ii0++){
		const real siglev=PR(siglevs, ii0, 0);
		real i0thres=MAX(0.1*siglev, rne*10);
		real nea2thres=pixthetax*pixthetay*100;
		//P(sanea,ii0)=dnew(nsa, 3);
		real i0sumsum=0;
		int ncrdisable=0;
		dmat* nea2=0;
		for(int isa=0; isa<nsa; isa++){
			real pixrot=0;//pixel rotation
			if(srot&&radpix){
				pixrot=P(PR(srot, ii0, 0), isa);
			}
			int cr=mtchcr;
			if(cr==-1){
				long fwhm=dfwhm_gauss(P(i0s, isa, ii0));
				if(fwhm>4){
					cr=1;
				} else{
					cr=0;
					ncrdisable++;
				}
			}
			dmat* bkgrnd2i=bkgrnd2?PR(bkgrnd2, isa, ii0):NULL;
			dmat* bkgrnd2ci=bkgrnd2c?PR(bkgrnd2c, isa, ii0):NULL;

			mtch(&P(*pmtche, isa, ii0), &nea2, P(i0s, isa, ii0),
				gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0, qe,
				bkgrnd2i, bkgrnd2ci, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				pixrot, radgx, cr);
			if(fabs(sigratio-1)>1e-5){
				dscale(P(i0s, isa, ii0), sigratio);
				if(gxs) dscale(P(gxs, isa, ii0), sigratio);
				if(gys) dscale(P(gys, isa, ii0), sigratio);
				mtch(NULL, &nea2, P(i0s, isa, ii0),
					gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
					qe, bkgrnd2i, bkgrnd2ci, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
					pixrot, radgx, cr);
				dscale(P(i0s, isa, ii0), sigratior);
				if(gxs) dscale(P(gxs, isa, ii0), sigratior);
				if(gys) dscale(P(gys, isa, ii0), sigratior);
			}
			
			real i0sum=dsum(P(i0s, isa, ii0));
			i0sumsum+=i0sum;
			if(pi0sum){
				P(*pi0sum, isa, ii0)=i0sum;
			}
			if(i0sum<i0thres||P(nea2, 0)>nea2thres||P(nea2, 3)>nea2thres){
			//Signal level too low or error to high.
				P(nea2, 0)=P(nea2, 3)=nea2thres;
				P(nea2, 1)=P(nea2, 2)=0;
				dset(P(*pmtche, isa, ii0), 0);
			}
			if(psanea){
				P(P(*psanea, ii0), isa, 0)=P(nea2, 0);
				P(P(*psanea, ii0), isa, 1)=P(nea2, 3);
				P(P(*psanea, ii0), isa, 2)=P(nea2, 1);
			}
		}/*isa  */
		dfree(nea2);

		if(mtchcr==-1){
			info("Mtched filter contraint are disabled for %d subaps out of %d.\n",	ncrdisable, nsa);
		}
		if(pi0sumsum){
			P(*pi0sumsum, ii0)=i0sumsum;
		}
	}/*ii0 */
}
void print_nea(const dcell* sanea, const dcell* sprint, const loc_t* saloc, const dcell *srsa){
	info2("Matched filter sanea:\n");
	if(sprint){/*print nea for select subapertures.*/
		for(int ii0=0; ii0<NX(sanea); ii0++){
			dmat *saindex=PR(sprint,ii0,0);
			info2("ii0 %d\n", ii0);
			info2("sa index   dist   noise equivalent angle\n");
			dmat* psanea=P(sanea, ii0);
			for(int ksa=0; ksa<NX(saindex); ksa++){
				int isa=(int)P(saindex, ksa);
				if(isa>0){
					info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n",
						isa, P(PR(srsa, ii0, 0), isa),
						sqrt(P(psanea, isa, 0))*206265000,
						sqrt(P(psanea, isa, 1))*206265000);
				}
			}
		}
	} else{
		real dsa=saloc->dx;
		real llimit=-dsa*0.6;
		real ulimit=dsa*0.4;
		info2("index: position noise equivalent angle\n");
		for(int isa=0; isa<saloc->nloc; isa++){
			real locx=saloc->locx[isa];
			real locy=saloc->locy[isa];
			if((locx>=0 && locy>llimit && locy<ulimit) || saloc->nloc<=4){
				info2("sa%4d:%6.1fm", isa, locx);
				for(int ii0=0; ii0<NX(sanea); ii0++){
					info2(" (%4.1f,%4.1f)",
						sqrt(P(P(sanea, ii0), isa, 0))*206265000,
						sqrt(P(P(sanea, ii0), isa, 1))*206265000);
				}//for ii0
				info2(" mas\n");
			}
		}/*isa  */
	}
}

void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs){
//	powfs[ipowfs].intstat->i0sum=dcellsum_each(powfs[ipowfs].intstat->i0);
//	powfs[ipowfs].intstat->i0sumsum=dsum_col(powfs[ipowfs].intstat->i0sum);
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
	info("Generating matched filter for %d\n", ipowfs);
	intstat_t* intstat=powfs[ipowfs].intstat;
	const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx;
	const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy;
	dmat* siglevs=parms->powfs[ipowfs].siglevs;
	real sigratio=parms->powfs[ipowfs].sigrecon>0?(parms->powfs[ipowfs].sigrecon/parms->powfs[ipowfs].siglev):1;
	const int mtchadp=parms->powfs[ipowfs].mtchadp;
	int mtchcr=mtchadp?-1:parms->powfs[ipowfs].mtchcr;
	mtch_wrap(&intstat->mtche, &powfs[ipowfs].sanea, &intstat->i0sum, &intstat->i0sumsum,
		intstat->i0, gxs, gys, parms->powfs[ipowfs].radpixtheta, parms->powfs[ipowfs].pixtheta,
		parms->powfs[ipowfs].rne, parms->powfs[ipowfs].qe, siglevs,
		bkgrnd, bkgrndc, powfs[ipowfs].bkgrnd, powfs[ipowfs].bkgrndc,
		sigratio, powfs[ipowfs].srot, parms->powfs[ipowfs].radpix, parms->powfs[ipowfs].radgx, mtchcr
	);
	print_nea(powfs[ipowfs].sanea, powfs[ipowfs].sprint, powfs[ipowfs].saloc, powfs[ipowfs].srsa);
}