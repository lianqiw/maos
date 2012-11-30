#define _GNU_SOURCE 
#include <sched.h>
#include <sys/types.h>

#include "../lib/aos.h"
#include "../cuda/gpu.h"

int main(int argc, char *argv[]){
    const char *host=myhostname();
    info2("host is %s\n", host);
    int ngpu;
    int gpus[8];
    if(!strcmp(host, "cassiopeia")){
	ngpu=2;
	gpus[0]=2;
	gpus[1]=3;
    }else if(!strcmp(host, "orion")){
	ngpu=1;
	gpus[0]=0;
    }else if(!strcmp(host, "geforce")){
	ngpu=2;
	gpus[0]=4;
	gpus[1]=6;
    }else{
	ngpu=1;
	gpus[0]=0;
    }
    
    cpu_set_t cpuset={0};
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
    
    int testcase=0;
    if(argc>1){
	testcase=strtol(argv[1], NULL, 10);
    }
    int nstep=8000;
    if(argc>2){
	nstep=strtol(argv[2], NULL, 10);
    }
    for(int jgpu=2; jgpu<=ngpu; jgpu++){
	if(testcase==0){
	    mvmfull_iwfs("mvm1.bin", "mvm2.bin", "pix1.bin", "pix2.bin", "mtch.bin", gpus, jgpu, nstep);
	}else{
	    mvm_iwfs("mvm1.bin", "mvm2.bin", "grad1.bin", "grad2.bin", gpus, jgpu, nstep);
	}
    }
}
