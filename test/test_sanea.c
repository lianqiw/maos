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


#include "../aos.h"
char *dirout;
static void test_ints(){
    rand_t init;
    seed_rand(&init,1);
    const int nwfs=6;
    lrand(&init);/*atm */
    rand_t wfs_rand[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	seed_rand(wfs_rand+iwfs,lrand(&init));
    }
    dcell *mtche=dcellread("powfs0_mtche.bin");
    int nsim=500;
    int nsa=2582;
    dmat *nea=dnew(nsim,nsa*2);
    real(*pnea)[nsim]=(void*)P(nea);
    real rne=3;
    real bkgrnd=0;
    real siglev=1000;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	for(int isim=0; isim<nsim; isim++){
	    dbg("iwfs %d isim=%d\n",iwfs,isim);
	    /*dcell *ints=dcellread("ints_%d_wfs%d.bin",isim,iwfs); */
	    dcell *ints=dcellread("ints_%d_wfs%d.bin",isim,iwfs);
	    /*dcell *i0=dcellread("powfs0_i0.bin"); */
	    dmat *im=NULL, *imy=NULL;
	    real gnf[2], gny[2];

	    for(int isa=0; isa<nsa; isa++){
		dcp(&im,P(ints,isa));
		dscale(im,siglev);
		gnf[0]=0; gnf[1]=0;
		dmulvec(gnf, P(mtche,isa), P(im),1.);
		gny[0]=0; gny[1]=0;
		dcp(&imy, im);
		addnoise(imy, &wfs_rand[iwfs], bkgrnd, bkgrnd, 0,0,rne);
		dmulvec(gny, P(mtche,isa), P(imy),1.);
		P(pnea,isim,isa)=gny[0]-gnf[0];
		P(pnea,isim,isa+nsa)=gny[1]-gnf[1];
	    }
	}
	writebin(nea,"test_sanea_wfs%d.bin",iwfs);
    }
}
/*
static void test_i0(){
    rand_t i0rand;
    seed_rand(&i0rand,100);
    dcell *mtche=dcellread("powfs0_mtche.bin");
    int nsim=500;
    int nsa=2582;
    dmat *nea=dnew(nsim,nsa*2);
    real(*pnea)[nsim]=(void*)P(nea);
    real rne=3;
    real bkgrnd=0;
    real siglev=1000;
    dcell *i0=dcellread("powfs0_i0.bin");
    for(int isim=0; isim<nsim; isim++){
	dbg("isim=%d\n",isim);
	dmat *im=NULL, *imy=NULL;
	real gnf[2], gny[2];
	for(int isa=0; isa<nsa; isa++){
	    dcp(&im,P(i0,isa));
	    dscale(im,siglev);
	    gnf[0]=0; gnf[1]=0;
	    dmulvec(gnf, P(mtche,isa), P(im),1.);
	    gny[0]=0; gny[1]=0;
	    dcp(&imy, im);
	    addnoise(imy, &i0rand, bkgrnd, bkgrnd,0,0,rne);
	    dmulvec(gny, P(mtche,isa), P(imy),1.);
	    P(pnea,isim,isa)=gny[0]-gnf[0];
	    P(pnea,isim,isa+nsa)=gny[1]-gnf[1];
	}
    }
    writebin(nea,"test_sanea_i0.bin");
    }*/
int main(){
    test_ints();
}
