/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <sched.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include "../lib/aos.h"
#include "../cuda/gpu.h"

int main(int argc, char *argv[]){
    info("HOST is %s\n", HOST);
    int ngpu;
    int gpus[8];
    if(!strcmp(HOST, "cassiopeia") || !strcmp(HOST, "kepler")){
	ngpu=2;
	gpus[0]=0;
	gpus[1]=1;
    }else if(!strcmp(HOST, "orion")){
	ngpu=2;
	gpus[0]=0;
	gpus[1]=1;
    }else if(!strcmp(HOST, "geforce")){
	ngpu=8;
	/*gpus[0]=4;
	  gpus[1]=6;*/
	for(int i=0; i<ngpu; i++){
	    gpus[i]=i;
	}
    }else{
	ngpu=1;
	gpus[0]=0;
    }
    set_realtime(1,-20);

    int testcase=0;
    if(argc>1){
	testcase=strtol(argv[1], NULL, 10);
    }
    int nstep=8000;
    if(argc>2){
	nstep=strtol(argv[2], NULL, 10);
    }
    switch(testcase){
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
	for(int jgpu=2; jgpu<=2; jgpu++){
	    switch(testcase){
	    case 0:
		mvmfull_iwfs(gpus, jgpu, nstep);
		break;
	    case 1:
		mvmfull_real(gpus, jgpu, nstep);
		break;
	    case 2:
		mvm_iwfs(gpus, jgpu, nstep);
		break;
	    case 3:
		mvm_only(gpus, jgpu, nstep);
		break;
	    case 4:
		mvmfull_pipe("mvm1.bin", "mvm2.bin", "pix1.bin", "pix2.bin", "mtch.bin", gpus, jgpu, nstep);
		break;
	    }
	}
	break;
    case 10:
	//mvm_test(gpus[0]);
	break;
    }
}
