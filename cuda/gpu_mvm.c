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

/**
   \file gpu_mvm.c
  */
#include <stdlib.h>
#include <stdio.h>
void gpu_mvm_daemon(int);
int main(int argc, char *argv[]){
    enum{
	P_EXE,
	P_PORT,
	P_TOT,
    };
    if(argc<P_TOT){
	fprintf(stderr,"Usage: %s port\n", argv[0]);
	exit(0);
    }
    int port=strtol(argv[P_PORT], NULL, 10);
    gpu_mvm_daemon(port);
}
