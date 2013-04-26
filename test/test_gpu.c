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
#ifdef __linux__    
    cpu_set_t cpuset={{0}};
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
#endif
    mlockall(MCL_FUTURE | MCL_CURRENT);
    //fault stack
    struct rlimit rl;
    if(!getrlimit(RLIMIT_STACK, &rl)){
	const int NSTACK=rl.rlim_cur/2;
	char tmp[NSTACK];
	memset(tmp, 0, NSTACK);
    }
    if(getuid()==0){
	info2("Set priority to -18\n");
	setpriority(PRIO_PROCESS, getpid(), -18);
	struct sched_param param;
	sched_getparam(getpid(), &param);
	param.sched_priority=sched_get_priority_max(SCHED_FIFO)-1;
	sched_setscheduler(getpid(), SCHED_FIFO, &param);
    }
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
    case 3:
	for(int jgpu=1; jgpu<=1; jgpu++){
	    switch(testcase){
	    case 0:
		mvmfull_iwfs(gpus, jgpu, nstep);
		break;
	    case 1:
		mvm_iwfs("mvm1.bin", "mvm2.bin", "grad1.bin", "grad2.bin", gpus, jgpu, nstep);
		break;
	    case 3:
		mvmfull_pipe("mvm1.bin", "mvm2.bin", "pix1.bin", "pix2.bin", "mtch.bin", gpus, jgpu, nstep);
		break;
	    }
	}
	break;
    case 2:
	mvm_test(gpus[0]);
	break;
    }
}
