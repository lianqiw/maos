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
#include <sys/file.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <errno.h>
#include "utils.h"
#include "curmat.h"
#include "cucmat.h"
#if !USE_CPP
extern "C"{
#endif
#include "../maos/mvm_client.h"
#if !USE_CPP
}
#endif
#include "cudata.h"

/**
   \file mvm_daemon.cu

   A standalone mvm daemon by gpu_mvm_daemon(). Waiting correction from MAOS,
   read pixels, and send DM commands back.

   This is not part of maos executable.
   
   For RTC benchmarking, use mvmfull_iwfs.cu which does both gradient processing
   and MVM, in a piplined way.
*/

real tic_save;
pthread_t *threads=NULL;
int *cmds=NULL;
int NGPU_MAX=0;

#define READ_CMD(p) stread(sock_mvm, p,N_CMD*sizeof(int))
#define WRITE_CMD(p) stwrite(sock_mvm, p,N_CMD*sizeof(int))

enum{
    CMD_MCOPY=1,
    CMD_GMUL,
    CMD_GMULF,
    CMD_DMCP,
    CMD_DMADD,
};

pthread_t thread_init;
int thread_init_joined=0;
typedef struct mvm_g_mul_t{
    int nact;
    int ngtot;
    int icol;
    int ngpu;
    int k;
    Array<GTYPE,Pinned>g;/*full g*/
    Array<Array<ATYPE,Pinned> > ac;/*a for each gpu.*/
    Array<ATYPE,Pinned> a;/*final a*/
    X(mat) *mvm;
    int *icols;
    int *kcols;
}mvm_t;
mvm_t *mvm_data=NULL;

static int mp_count;/*number multi-processors on each GPU*/

int ndone;
int nleft;
char *start;

int sock_mvm;
__global__ static void mvm_g_mul_do(const Real *restrict mvm, ATYPE *restrict a, const GTYPE *restrict g, int nact, int ng){
    extern __shared__ Real acc[];
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	acc[threadIdx.x]=0;
	for(int ig=0; ig<ng; ig++){
	    register Real mvmi=mvm[nact*ig+iact];
	    acc[threadIdx.x]+=mvmi*(Real)(g[ig]);
	}
	a[iact]+=(ATYPE)acc[threadIdx.x];
    }
}
/*A couple of threads that does jobs upon available.*/
static void mvm_thread(void* ithread0){
    long ithread=(long) ithread0;
    int nact=mvm_data->nact;
    int nact1=(nact+NCPU-1)/NCPU;
    int iact1=nact1*ithread;
    int iact2=MIN(iact1+nact1, nact);
    int naeach=(nact+mp_count-1)/mp_count;
    if(ithread<NGPU){
	gpu_set(ithread);
    }
    while(1){
	switch(cmds[ithread]){
	case 0:
	    usleep(10);
	    break;
	case CMD_MCOPY:{
	    if(ithread<NGPU){
		gpu_set(ithread);
	    }
	    X(mat) *mvm=mvm_data->mvm;
	    cp2gpu(cudata->mvm_m, mvm, cudata->mvm_stream);
	    cudata->mvm_a.init(mvm->nx,1);
	    cudata->mvm_g.init(mvm->ny,1);
	  
	    cudaMemsetAsync(cudata->mvm_a(), 0, mvm->nx*sizeof(ATYPE), cudata->mvm_stream);
	    cudaStreamSynchronize(cudata->mvm_stream);
	    cmds[ithread]=0;
	}
	    break;
	case CMD_GMUL:{
	    int icol=mvm_data->icols[ithread];
	    int k=mvm_data->kcols[ithread];
	    cudaMemcpyAsync(cudata->mvm_g()+icol, mvm_data->g()+icol, k*sizeof(GTYPE), 
			    cudaMemcpyHostToDevice, cudata->mvm_stream);
	    
	    mvm_g_mul_do<<<mp_count, naeach, sizeof(Real)*naeach, cudata->mvm_stream>>>
		(cudata->mvm_m.Col(icol), cudata->mvm_a(), cudata->mvm_g()+icol, nact, k);
	 
	    cmds[ithread]=0;
	}
	    break;
	case CMD_GMULF:{
	    int ki=(mvm_data->k+NGPU-1)/NGPU;
	    int icol=mvm_data->icol+ki*ithread;
	    int ki2=mvm_data->k+mvm_data->icol-icol;
	    int k=MIN(ki, ki2);
	    cudaMemcpyAsync(cudata->mvm_g()+icol, mvm_data->g()+icol, k*sizeof(GTYPE), 
			    cudaMemcpyHostToDevice, cudata->mvm_stream);
	    mvm_g_mul_do<<<mp_count, naeach, sizeof(Real)*naeach, cudata->mvm_stream>>>
		(cudata->mvm_m.Col(icol), cudata->mvm_a(), cudata->mvm_g()+icol, nact, k);

	    cmds[ithread]=0;
	}
	    break;
	case CMD_DMCP:{
	    cudaMemcpyAsync(mvm_data->ac[ithread], cudata->mvm_a, nact*sizeof(ATYPE),
			    cudaMemcpyDeviceToHost, cudata->mvm_stream);
	    cudaStreamSynchronize(cudata->mvm_stream);
	    cudata->mvm_a.zero(cudata->mvm_stream);
	    cmds[ithread]=0;
	}
	    break;
	case CMD_DMADD:{
	    for(int iact=iact1; iact<iact2; iact++){
		register ATYPE temp=0;
		for(int igpu=0; igpu<NGPU; igpu++){
		    temp+=mvm_data->ac[igpu][iact];
		}
		mvm_data->a[iact]=temp;
	    }
	    cmds[ithread]=0;//mark we are done.
	}
	    break;
	    
	}
    }
}
static void mvm_data_free(void){
    free(mvm_data->icols);
    free(mvm_data->kcols);
    free(mvm_data);
}
real tim_cmd=0, tim_gsend=0, tim_gcp=0, tim_dmcp=0, tim_queue=0, tim_dmsum=0, tim_dmsend=0;
static int respond(int sock){
    TIC;tic;
    sock_mvm=sock;
    int cmd[N_CMD]={0};
    READ_CMD(cmd);
    static real tim_gfirst=0;
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	int ngpu=cmd[1];
	int nact=cmd[2];
	int ngtot=cmd[3];
	dbg("Receiving mvm %dx%d\n", nact, ngtot);
	if(mvm_data){
	    mvm_data_free();
	}
	mvm_data=(mvm_t*)calloc(1, sizeof(mvm_t));
	mvm_data->nact=nact;
	mvm_data->ngtot=ngtot;
	mvm_data->mvm=X(new)(nact, ngtot);
	stread(sock_mvm, mvm_data->mvm->p, (nact*ngtot)*sizeof(Real));
	if(!thread_init_joined){
	    pthread_join(thread_init, NULL);
	    thread_init_joined=1;
	}
	if(ngpu>0 && ngpu<=NGPU_MAX){
	    NGPU=ngpu;
	    warning("Using %d GPUs\n", ngpu);
	}
	mvm_data->ngpu=NGPU;
	mvm_data->g.init(ngtot,1);
	mvm_data->a.init(nact,1);
	mvm_data->ac.init(NGPU,1);
	for(int ig=0; ig<NGPU; ig++){
	    mvm_data->ac[ig].init(nact,1);
	}
		
	mvm_data->icols=(int*)calloc(NGPU, sizeof(int));
	mvm_data->kcols=(int*)calloc(NGPU, sizeof(int));

	toc("Read mvm");tic;
	if(!threads){
	    threads=(pthread_t*)calloc(NCPU, sizeof(pthread_t));
	    cmds=(int*)calloc(NCPU, sizeof(int));
	    for(long i=0; i<NCPU; i++){
		cmds[i]=CMD_MCOPY;
		pthread_create(threads+i, NULL, (void*(*)(void*))mvm_thread, (void*)i);
	    }
	}
	//Wait for all copying to be finish.
	for(int i=0; i<NGPU; i++){
	    while(cmds[i]==CMD_MCOPY){
		usleep(1);
	    }
	}
	toc("copy mvm to gpu");
	info("done");
	X(free)(mvm_data->mvm);
    }
	break;
    case GPU_MVM_G:{/*maos sends gradients*/
	assert(cmd[1]==0);
	int ngeach=cmd[2];
	tim_cmd+=toc3; tic;
	int ngtot=mvm_data->ngtot;
	tim_gfirst=myclockd();
	for(int icol=cmd[1]; icol<ngtot; icol+=ngeach){
	    int k=MIN(ngeach, ngtot-icol);
	    stread(sock_mvm, mvm_data->g()+icol, k*sizeof(GTYPE));
	    tim_gsend+=toc3;tic;
	    if(cmd[2]<1800){//Use next available GPU to handle this task.
		int igpu=gpu_next();
		while(cmds[igpu]){//wait until last job is finished.
		    usleep(10);
		}
		mvm_data->icols[igpu]=icol;
		mvm_data->kcols[igpu]=k;
		cmds[igpu]=CMD_GMUL;//Notify new job.
	    }else{ //Let each GPU figure out which part of gradients to copy and compute.
		mvm_data->icol=icol;
		mvm_data->k=k;
		//CALL_THREAD(mvm_data->mvm_g_mul, NGPU, 0);
		for(int i=0; i<NGPU; i++){
		    cmds[i]=CMD_GMULF;
		}
		for(int i=0; i<NGPU; i++){
		    while(cmds[i]==CMD_GMULF){
			usleep(1);
		    }
		}
		tim_queue+=toc3;tic;
	    }
	}
	for(int i=0; i<NGPU; i++){
	    while(cmds[i]){
		usleep(10);
	    }
	    //Copy partial DM commands back.
	    cmds[i]=CMD_DMCP;
	}
	for(int i=0; i<NGPU; i++){
	    while(cmds[i]==CMD_DMCP){
		usleep(10);
	    }
	}
	tim_dmcp+=toc3;tic;
	tic_save=tic;
	//Ask each thread to sum the DM vector.
	for(int i=0; i<NCPU; i++){
	    cmds[i]=CMD_DMADD;
	}
	for(int i=0; i<NCPU; i++){
	    while(cmds[i]==CMD_DMADD){
		usleep(10);
	    }
	}
	tim_dmsum+=toc3;tic;
	stwrite(sock_mvm, mvm_data->a, mvm_data->nact*sizeof(ATYPE));
	tim_dmsend+=toc3;tic;
	info("k=%4d CMD %1.0f, gsend %2.0f, gcp %3.0f, queue %3.0f, sync %3.0f sum %3.0f, send %2.0f, total %4.0f\n", ngeach,
	      tim_cmd*1e6, tim_gsend*1e6, tim_gcp*1e6, tim_queue*1e6, tim_dmcp*1e6, 
	      tim_dmsum*1e6, tim_dmsend*1e6, (myclockd()-tim_gfirst)*1e6);
	tim_cmd=tim_gsend=tim_gcp=tim_dmcp=tim_queue=tim_dmsum=tim_dmsend=0;
	
    }
	break;
    }
    return 0;
}
static void* gpu_mvm_gpu_init(void *A){
    (void)A;
    gpu_init(NULL, 0, 0);
    DO(cudaFuncSetCacheConfig(mvm_g_mul_do, cudaFuncCachePreferShared));
    struct cudaDeviceProp prop;
    DO(cudaGetDeviceProperties(&prop, 0));
    mp_count=prop.multiProcessorCount;
    /*
    switch(prop.major){
    case 1:
	mp_core=8;
	break;
    case 2:{
	switch(prop.minor){
	case 0:
	    mp_core=32;
	    break;
	case 1:
	    mp_core=48;
	    break;
	default:
	    error("Please fill in this");
	}
	break;
    }
    default:
	error("Please fill on this");
    }*/
    /*Creating stream for the first time is slow. So do it here to avoid latency
      later. Do not init stream_t in multithread. It is slower.*/
    TIC;tic;
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	stream_t temp;
	toc("init gpu");
    }
    NCPU=NGPU;
    NGPU_MAX=NGPU;
    return NULL;
}
void gpu_mvm_daemon(int port){
    /* GPU initialization may take a few seconds to finish. Launch in a separate
     * thread and join before using GPU.*/
    pthread_create(&thread_init, NULL, gpu_mvm_gpu_init, NULL);
    info("Starting MVM daemon at port %d\n", port);
    redirect();
    listen_port(port, NULL, respond, -1, NULL, 1);
}
