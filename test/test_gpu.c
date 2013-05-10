#define _GNU_SOURCE 
#include <sched.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include "../lib/aos.h"
#include "../cuda/gpu.h"

int main(int argc, char *argv[]){
    const char *host=myhostname();
    info2("host is %s\n", host);
    int ngpu;
    int gpus[8];
    if(!strcmp(host, "cassiopeia") || !strcmp(host, "kepler")){
	ngpu=2;
	gpus[0]=0;
	gpus[1]=1;
    }else if(!strcmp(host, "orion")){
	ngpu=2;
	gpus[0]=0;
	gpus[1]=1;
    }else if(!strcmp(host, "geforce")){
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
    set_realtime(0,-20);

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
	for(int jgpu=1; jgpu<=1; jgpu++){
	    switch(testcase){
	    case 0:
		mvmfull_iwfs(gpus, jgpu, nstep);
		break;
	    case 1:
		mvm_iwfs(gpus, jgpu, nstep);
		break;
	    case 2:
		mvm_only(gpus, jgpu, nstep);
		break;
	    case 3:
		mvmfull_pipe("mvm1.bin", "mvm2.bin", "pix1.bin", "pix2.bin", "mtch.bin", gpus, jgpu, nstep);
		break;
	    }
	}
	break;
    case 10:
	mvm_test(gpus[0]);
	break;
    }
}
