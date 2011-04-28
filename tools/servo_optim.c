/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
/**
   standalone routine that does servo filtering.
*/
#include "aos.h"
int main(int argc, char **argv){
    if(argc<6){
	info2("Usage: %s psd.bin dtrat dt sigma gainout.bin\n", argv[0]);
	exit(1);
    }
    dmat *psd=dread("%s",argv[1]);
    long dtrat=strtol(argv[2],NULL,10);
    double dt=strtod(argv[3],NULL);
    double sigma=strtod(argv[4],NULL);//m^2
    dmat *sigma2=dnew(1,1); sigma2->p[0]=sigma;
    dcell *gain=servo_typeII_optim(psd,dtrat,dt,sigma2);
    dwrite(gain->p[0],"%s",argv[5]);
    dfree(sigma2);
    dcellfree(gain);
    dfree(psd);
}

