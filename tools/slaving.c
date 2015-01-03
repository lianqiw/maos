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
#include "../lib/aos.h"
/**
   \file slaving.c
   
   Warps of slaving() for use standalone.

   Usage: slaving loc.bin HA.bin amp.bin
*/
int main(int argc, char **argv){
    if(argc<3){
	info2("Usage:\n"
	      "slaving loc.bin HA.bin [amp.bin] [thres]\n"
	      "loc.bin contains the grid. \n"
	      "HA.bin  contains the influence function from loc to destination\n"
	      "amp.bin contains the amplitude weighting in the destination (optional)\n"
	      "thres   contains the threshold to treat actuator as not coupled [0 1] (optional)\n"
	      "The result will be saved in slaving.bin\n"
	      );				
	exit(0);
    }
    loccell *aloc=cellnew(1,1);
    aloc->p[0]=locread("%s",argv[1]);
    dsp *ha=dspread("%s",argv[2]);
    double thres=argc>3?strtod(argv[3],NULL):0.5;
    dspcell*has=cellnew(1,1); has->p[0]=ha;
    dcell *actcpl=genactcpl(has, NULL);
    dspcell *slave=slaving(aloc, actcpl, NULL, NULL, NULL, thres, 1);
    dcellfree(actcpl);
    dspwrite(slave->p[0],"slaving");
    dspcellfree(slave);
}
