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
/*
  Test the difference in generated atmosphere in taurus and opteron.
*/
#include "../lib/aos.h"
int main(){
    rand_t rstat;
    rand_t bstat;
    seed_rand(&bstat,1);
    seed_rand(&rstat,lrand(&bstat));
    int nlayer=2;
    real wt[7]={0.2887, 0.17795, 0.06602, 0.07833, 0.1405, 0.1216, 0.1269};
    real r0=0.1987; 
    real L0[7]={30,30,30,30,30,30,30};
    real dx=1./64.;
    int m=4096*2;
    NTHREAD=6;
    THREAD_POOL_INIT(NTHREAD);
    GENATM_T data;
    memset(&data, 0, sizeof(GENATM_T));
    data.rstat=&rstat;
    data.nx=m;
    data.ny=m;
    data.dx=dx;
    data.r0=r0;
    data.L0=L0;
    data.wt=wt;
    data.nlayer=nlayer;
    mapcell *map=vonkarman_screen(&data);
    writebin(map, "atm.bin");
}
