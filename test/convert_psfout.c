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


#include "../lib/aos.h"
int main(){
    ccell *psf;
    int nsim=5000;
    for(int iwfs=6; iwfs<64; iwfs++){
	zfarr *ca=zfarr_init(nsim,1,"psfout_wfs%d.bin",iwfs);
	for(int isim=0; isim<nsim; isim++){
	    if(isim%100==0)
		dbg("iwfs=%d, isim=%d\n",iwfs,isim);
	    psf=ccellread("psfout/psfout_wfs%d_isim%d.bin",iwfs,isim);
	    zfarr_push(ca, isim, psf);
	    ccellfree(psf);
	}
	zfarr_close(ca);
    }
}
