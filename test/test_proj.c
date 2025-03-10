/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

/*
  Rectangular grid. (different X/Y Spacing)
 */

static void test_grid_proj(){
    if(!zfexist("M3_p.bin")){
	info("M3_p.bin does not exist\n");
	return;
    }
    dmat *junk=dread("M3_p.bin");
    dmat *X=dread("M3_x.bin");
    dmat *Y=dread("M3_y.bin");
    dmat *tmp=dread("M3_theta.bin");
    real bx=P(tmp,0);
    real by=P(tmp,1);
    dfree(tmp);
    rmap_t *mapin=mycalloc(1,rmap_t);
    P(mapin)=P(junk);
    mapin->ox=P(X,0);
    mapin->oy=P(Y,0);
    mapin->dx=P(X,1)-P(X,0);
    mapin->dy=P(Y,1)-P(Y,0);
    mapin->nx=X->nx*X->ny;
    mapin->ny=Y->nx*Y->ny;
    real d_m3_f=20.;/*from m3 to focus */
    real d_exitpupil_f=46.38661051;
    /*real d_m2_m3=23.59375000; */
    /*real d_m2_f=43.59375000; */
    real d_exitpupil_m3=d_exitpupil_f-d_m3_f;
    real r_exitpupil=1.546220350;
    real r_pupil=15;

    if(P(X,X->nx*X->ny-1)-P(X,0)-mapin->dx*X->nx*X->ny>1.e-10){
	    error("X has non even spacing\n");
    }
    if(P(Y,Y->nx*Y->ny-1)-P(Y,0)-mapin->dy*Y->nx*Y->ny>1.e-10){
	    error("Y has non even spacing\n");
    }
    free(junk);/*don't dfree */
    dfree(X); dfree(Y);
    /*direction of guide star */
    /*loc_t *loc2=mksqloc2(2000,2000,1./64.); */
    loc_t* loc2=locread("aper_locs.bin");
    
    dmat *amp=dread("aper_amp.bin");
    dmat *phi2=dnew(loc2->nloc,1);
    proj_rect_grid(P(phi2),mapin,M_PI*0.75,M_PI*0.5,
		   loc2,-r_exitpupil/r_pupil,r_exitpupil/r_pupil,
		   P(amp),-2,d_exitpupil_f,d_exitpupil_m3,-bx,-by);
    /*drawopdamp("test_proj",loc2,phi2,amp,"phi"); */
    
    writebin(phi2,"phi");
    /*locwrite(loc2,"loc"); */
}
int main(){
    test_grid_proj();
}
