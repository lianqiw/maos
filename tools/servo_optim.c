/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   \file servo_optim.c

   Standalone routine that does servo filtering.
*/
#include "../lib/aos.h"
int main(int argc, char **argv){
    enum{
	P_EXE,
	P_OUT,
	P_PSD,
	P_DT,
	P_DTRAT,
	P_SIGMAN,
	P_STYPE,
	P_TOT,
    };
    if(argc!=P_TOT){
	info2("Usage: %s gainout.bin psd.bin dt dtrat sigman servotype\n", argv[0]);
	info2(
	      "PSD psd.bin should be in m^2/hz\n"
	      "dt is the AO fundemental sampling period.\n"
	      "dtrat is ratio of sampling period of the WFS over dt.\n"
	      "sigman is wavefront error due to noise in m^2.\n"
	      "servotype is 1 or 2 for typeI, typeII controller.\n"
	      );
	exit(1);
    }
    dmat *psd=dread("%s",argv[P_PSD]);
    double dt=strtod(argv[P_DT],NULL);
    long dtrat=strtol(argv[P_DTRAT],NULL,10);
    double sigma=strtod(argv[P_SIGMAN],NULL);/*m^2 */
    dmat *sigma2=dnew(1,1); sigma2->p[0]=sigma;
    int servotype=(int)strtol(argv[P_STYPE],NULL,10);
    dcell *gain=servo_optim(psd,dt,dtrat,M_PI/4,sigma2,servotype);
    writebin(gain->p[0],"%s",argv[P_OUT]);
    double gain_n;
    double res=servo_residual(&gain_n, psd, dt, dtrat,  gain->p[0], servotype);
    info2("sigma res=%g, noise prop=%g. signal res=%g noise prop=%g\n", 
	  gain->p[0]->p[3], gain->p[0]->p[4], res, sigma*gain_n);
    dfree(sigma2);
    dcellfree(gain);
    dfree(psd);
}

