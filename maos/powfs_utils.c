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
/*
  2009-11-26: changed to rotate OTF instead of psf to comply with the approach
  in wfsints.  this gives slightly larger intensity because max value of OTF is
  preserved which corresponds to the sum of psf.


*/

#include "common.h"
#include "powfs_utils.h"


void print_nea(const dcell* sanea, const lcell* sprint, const loc_t* saloc){
	//info2("Matched filter sanea:\n");
	if(sprint){/*print nea for select subapertures.*/
		for(int ii0=0; ii0<NX(sanea); ii0++){
			lmat* saindex=PR(sprint, ii0, 0);
			info2("sa index   location       noise equivalent angle\n");
			dmat* psanea=P(sanea, ii0);
			for(int ksa=0; ksa<NX(saindex); ksa++){
				long isa=P(saindex, ksa);
				if(isa>0){
					info2("%8ld: (%5.1f, %5.1f) m, (%6.2f, %6.2f) mas\n",
						isa, saloc->locx[isa], saloc->locy[isa],
						sqrt(P(psanea, isa, 0))*RAD2MAS,
						sqrt(P(psanea, isa, 1))*RAD2MAS);
				}
			}
		}
	} else{
		real dsa=saloc->dx;
		real llimit=-dsa*0.6;
		real ulimit=dsa*0.4;
		info2("sa index   location       noise equivalent angle\n");
		for(int isa=0; isa<saloc->nloc; isa++){
			real locx=saloc->locx[isa];
			real locy=saloc->locy[isa];
			if((locx>=0&&locy>llimit&&locy<ulimit)||saloc->nloc<=4){
				info2("%8d: (%5.1f, %5.1f) m", isa, locx, locy);
				for(int ii0=0; ii0<NX(sanea); ii0++){
					info2(" (%6.2f, %6.2f)",
						sqrt(P(P(sanea, ii0), isa, 0))*RAD2MAS,
						sqrt(P(P(sanea, ii0), isa, 1))*RAD2MAS);
				}//for ii0
				info2(" mas\n");
			}
		}/*isa  */
	}
}

void mtch_wrap(const parms_t* parms, powfs_t* powfs, const int ipowfs){
	info("Generating matched filter for powfs %d\n", ipowfs);
	const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
	const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;

	intstat_t* intstat=powfs[ipowfs].intstat;
	const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx;
	const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy;
	real sigrecon=parms->powfs[ipowfs].sigrecon;
	if(sigrecon<=0){
		if(parms->powfs[ipowfs].siglev>1200){
			info("powfs%d: sigrecon is limited to %g for too high siglev %g\n", ipowfs, sigrecon, parms->powfs[ipowfs].siglev);
			sigrecon=1200;
		}else{
			sigrecon=parms->powfs[ipowfs].siglev;
		}
	}
	real sigratio=sigrecon/parms->powfs[ipowfs].siglev;
	const int mtchadp=parms->powfs[ipowfs].mtchadp;
	int mtchcr=mtchadp?-1:parms->powfs[ipowfs].mtchcr;
	mtch_cell(&intstat->mtche, &powfs[ipowfs].sanea, &intstat->i0sum, &intstat->i0sumsum,
		intstat->i0, gxs, gys, parms->powfs[ipowfs].qe,
		powfs[ipowfs].bkgrnd, powfs[ipowfs].bkgrndc, bkgrnd, bkgrndc,
		parms->powfs[ipowfs].rne, parms->powfs[ipowfs].radpixtheta, parms->powfs[ipowfs].pixtheta,
		parms->powfs[ipowfs].radpix?powfs[ipowfs].srot:NULL, parms->powfs[ipowfs].radgx, mtchcr, sigratio
	);
	if(powfs[ipowfs].srsa){
		print_nea(powfs[ipowfs].sanea, powfs[ipowfs].sprint, powfs[ipowfs].saloc);
	}
}

