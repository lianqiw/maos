/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "utils.h"
#include "curmat.h"
#include "cucmat.h"
extern "C"{
#include "../maos/mvm_client.h"
#include <sys/file.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <errno.h>
}
double tic_save;
pthread_t *threads=NULL;
int *cmds=NULL;
int NGPU_MAX=0;
/*Read into double array until we get all*/
#define READ_ARR(p,n,type)					\
    nleft=(n)*sizeof(type);					\
    start=(char*)(p);						\
    do{								\
	int nread=read(sock_mvm, start, nleft);			\
	nleft-=nread;						\
	start+=nread;						\
	if(nread<=0){						\
	    warning("nread=%d, nleft=%d\n", nread, nleft);	\
	    return -1;						\
	}							\
    }while(nleft>0);

#define WRITE_ARR(p,n,type)						\
    if((ndone=write(sock_mvm, p, sizeof(type)*(n)))!=sizeof(type)*(n)){	\
	perror("write");						\
	warning("error writing. want %ld wrote %d\n",(n)*sizeof(type), ndone); \
	return -1;							\
    }

#define READ_CMD(p) READ_ARR(p,N_CMD,int)
#define WRITE_CMD(p) WRITE_ARR(p,N_CMD,int)

enum{
    CMD_MCOPY=1,
    CMD_GMUL,
    CMD_GMULF,
    CMD_DMCP,
    CMD_DMADD,
};

/**
 * 2012-06-22: Serious bug found:
 * Initially was using curnew. But later replaced by new while should be using cudaMalloc
 * This caused misterious bugs, such as no answer or even kernel panic.
 */
pthread_t thread_init;
int thread_init_joined=0;
typedef struct mvm_g_mul_t{
    int nact;
    int ngtot;
    int icol;
    int ngpu;
    int k;
    GTYPE *g;/*full g*/
    ATYPE **ac;/*a for each gpu.*/
    ATYPE *a;/*final a*/
    smat *mvm;
    int *icols;
    int *kcols;
}mvm_t;
mvm_t *mvm_data=NULL;

int mp_count;/*number multi-processors on each GPU*/
int mp_core;/*number of cuda cores per multi-processor*/

int ndone;
int nleft;
char *start;

int sock_mvm;
__global__ static void mvm_g_mul_do(float *restrict mvm, ATYPE *restrict a, const GTYPE *restrict g, int nact, int ng){
    extern __shared__ float acc[];
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	acc[threadIdx.x]=0;
	for(int ig=0; ig<ng; ig++){
	    register float mvmi=mvm[nact*ig+iact];
	    acc[threadIdx.x]+=mvmi*(float)(g[ig]);
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
	    /*info2("thread %d sleeped at %.1f\n", ithread, (myclockd()-tic_start)*1e6);
	    usleep(1);
	    info2("thread %d wake up at %.1f\n", ithread, (myclockd()-tic_start)*1e6);*/
	    break;
	case CMD_MCOPY:{
	    if(ithread<NGPU){
		gpu_set(ithread);
	    }
	    smat *mvm=mvm_data->mvm;
	    cudata->mvm_stream=new stream_t;
	    cp2gpu(&cudata->mvm_m, mvm, *cudata->mvm_stream);
	    if(cudata->mvm_a){
		cudaFree(cudata->mvm_a); cudata->mvm_a=NULL;
		cudaFree(cudata->mvm_g); cudata->mvm_g=NULL;
	    }
	    cudaMalloc(&cudata->mvm_a, mvm->nx*sizeof(ATYPE));
	    cudaMalloc(&cudata->mvm_g, mvm->ny*sizeof(GTYPE));
	    cudaMemsetAsync(cudata->mvm_a, 0, mvm->nx*sizeof(ATYPE), *cudata->mvm_stream);
	    cudaStreamSynchronize(*cudata->mvm_stream);
	    if(cudata->mvm_a2){
		free(cudata->mvm_a2);
	    }
	    cudata->mvm_a2=(ATYPE**)calloc(NGPU, sizeof(ATYPE*));
	    cmds[ithread]=0;
	}
	    break;
	case CMD_GMUL:{
	    int icol=mvm_data->icols[ithread];
	    int k=mvm_data->kcols[ithread];
	    cudaMemcpyAsync(cudata->mvm_g+icol, mvm_data->g+icol, k*sizeof(GTYPE), 
			    cudaMemcpyHostToDevice, *cudata->mvm_stream);
	    
	    mvm_g_mul_do<<<mp_count, naeach, sizeof(float)*naeach, *cudata->mvm_stream>>>
	      (cudata->mvm_m->p+nact*icol, cudata->mvm_a, cudata->mvm_g+icol, nact, k);
	    /*float one=1;
	    DO(cublasSgemv(cudata->mvm_stream[0], CUBLAS_OP_N, nact, k, &one, cudata->mvm_m->p+nact*icol,
	    nact, cudata->mvm_g+icol, 1, &one, cudata->mvm_a, 1));*/
	    cmds[ithread]=0;
	}
	    break;
	case CMD_GMULF:{
	    int ki=(mvm_data->k+NGPU-1)/NGPU;
	    int icol=mvm_data->icol+ki*ithread;
	    int ki2=mvm_data->k+mvm_data->icol-icol;
	    int k=MIN(ki, ki2);
	    cudaMemcpyAsync(cudata->mvm_g+icol, mvm_data->g+icol, k*sizeof(GTYPE), 
			    cudaMemcpyHostToDevice, *cudata->mvm_stream);
	    mvm_g_mul_do<<<mp_count, naeach, sizeof(float)*naeach, *cudata->mvm_stream>>>
	    (cudata->mvm_m->p+nact*icol, cudata->mvm_a, cudata->mvm_g+icol, nact, k);
	    /*float one=1;
	    DO(cublasSgemv(cudata->mvm_stream[0], CUBLAS_OP_N, nact, k, &one, cudata->mvm_m->p+nact*icol,
	    nact, cudata->mvm_g+icol, 1, &one, cudata->mvm_a, 1));*/
	    cmds[ithread]=0;
	}
	    break;
	case CMD_DMCP:{
	    cudaMemcpyAsync(mvm_data->ac[ithread], cudata->mvm_a, nact*sizeof(ATYPE),
			    cudaMemcpyDeviceToHost, *cudata->mvm_stream);
	    cudaStreamSynchronize(*cudata->mvm_stream);
	    cudaMemsetAsync(cudata->mvm_a, 0, nact*sizeof(ATYPE), *cudata->mvm_stream);
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
    cudaFreeHost(mvm_data->g);
    cudaFreeHost(mvm_data->a);
    for(int i=0; i<mvm_data->ngpu; i++){
	cudaFreeHost(mvm_data->ac[i]);
    }
    free(mvm_data->ac);
    free(mvm_data->icols);
    free(mvm_data->kcols);
    free(mvm_data);
}
double tim_cmd=0, tim_gsend=0, tim_gcp=0, tim_dmcp=0, tim_queue=0, tim_dmsum=0, tim_dmsend=0;
static int respond(int sock){
    TIC;tic;
    sock_mvm=sock;
    int cmd[N_CMD];
    READ_CMD(cmd);
    static double tim_gfirst=0;
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	int ngpu=cmd[1];

	int nact=cmd[2];
	int ngtot=cmd[3];
	info("Receiving mvm %dx%d\n", nact, ngtot);
	if(mvm_data){
	    mvm_data_free();
	}
	mvm_data=(mvm_t*)calloc(1, sizeof(mvm_t));
	mvm_data->nact=nact;
	mvm_data->ngtot=ngtot;
	mvm_data->mvm=snew(nact, ngtot);
	READ_ARR(mvm_data->mvm->p, (nact*ngtot), float);
	if(!thread_init_joined){
	    pthread_join(thread_init, NULL);
	    thread_init_joined=1;
	}
	if(ngpu>0 && ngpu<=NGPU_MAX){
	    NGPU=ngpu;
	    warning("Using %d GPUs\n", ngpu);
	}
	mvm_data->ngpu=NGPU;
	mvm_data->ac=new ATYPE*[NGPU];
	cudaMallocHost(&mvm_data->g, sizeof(GTYPE)*ngtot);
	cudaMallocHost(&mvm_data->a, sizeof(ATYPE)*nact);
	memset(mvm_data->a, 0., nact*sizeof(ATYPE));
	for(int ig=0; ig<NGPU; ig++){
	    cudaMallocHost(&mvm_data->ac[ig], sizeof(ATYPE)*nact);
	}
	mvm_data->icols=(int*)calloc(NGPU, sizeof(int));
	mvm_data->kcols=(int*)calloc(NGPU, sizeof(int));

	toc22("Read mvm");tic;
	if(!threads){
	    threads=(pthread_t*)calloc(NCPU, sizeof(pthread_t));
	    cmds=(int*)calloc(NCPU, sizeof(int));
	    for(int i=0; i<NCPU; i++){
		pthread_create(threads+i, NULL, (void*(*)(void*))mvm_thread, (void*)i);
		cmds[i]=0;
	    }
	}
	for(int i=0; i<NGPU; i++){
	    cmds[i]=CMD_MCOPY;
	}
	for(int i=0; i<NGPU; i++){
	    while(cmds[i]==CMD_MCOPY){
		usleep(1);
	    }
	}
	toc22("copy mvm to gpu");
	info2("done");
	sfree(mvm_data->mvm);
    }
	break;
    case GPU_MVM_G:{/*maos sends gradients*/
	int icol=cmd[1];/*starting column*/
	assert(icol==0);
	int ngeach=cmd[2];
	tim_cmd+=toc3; tic;
	int ngtot=mvm_data->ngtot;
	tim_gfirst=myclockd();
	for(int icol=cmd[1]; icol<ngtot; icol+=ngeach){
	    int k=MIN(ngeach, ngtot-icol);
	    READ_ARR(mvm_data->g+icol, k, GTYPE);
	    tim_gsend+=toc3;tic;
	    if(cmd[2]<1800){//part of grads
		int igpu=gpu_next();
		while(cmds[igpu]){
		}
		mvm_data->icols[igpu]=icol;
		mvm_data->kcols[igpu]=k;
		cmds[igpu]=CMD_GMUL;
	    }else{ //Send to different gpus
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
	    }
	    cmds[i]=CMD_DMCP;
	}
	for(int i=0; i<NGPU; i++){
	    while(cmds[i]==CMD_DMCP){
		usleep(1);
	    }
	}
	//CALL_THREAD(mvm_data->mvm_a_cp, NGPU, 0);
	tim_dmcp+=toc3;tic;
	tic_save=tic;
	for(int i=0; i<NCPU; i++){
	    cmds[i]=CMD_DMADD;
	}
	for(int i=0; i<NCPU; i++){
	    while(cmds[i]==CMD_DMADD){
		usleep(1);
	    }
	}
	tim_dmsum+=toc3;tic;
	WRITE_ARR(mvm_data->a, mvm_data->nact, ATYPE);
	tim_dmsend+=toc3;tic;
	info2("k=%4d CMD %1.0f, gsend %2.0f, gcp %3.0f, queue %3.0f, sync %3.0f sum %3.0f, send %2.0f, total %4.0f\n", ngeach,
	      tim_cmd*1e6, tim_gsend*1e6, tim_gcp*1e6, tim_queue*1e6, tim_dmcp*1e6, 
	      tim_dmsum*1e6, tim_dmsend*1e6, (myclockd()-tim_gfirst)*1e6);
	tim_cmd=tim_gsend=tim_gcp=tim_dmcp=tim_queue=tim_dmsum=tim_dmsend=0;
	
    }
	break;
    }
    return 0;
}
void* gpu_mvm_gpu_init(void* A){
    (void)A;
    gpu_init(NULL, 0);
    DO(cudaFuncSetCacheConfig(mvm_g_mul_do, cudaFuncCachePreferShared));
    struct cudaDeviceProp prop;
    DO(cudaGetDeviceProperties(&prop, 0));
    mp_count=prop.multiProcessorCount;
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
    }
    /*Creating stream for the first time is slow. So do it here to avoid latency
      later. Do not init stream_t in multithread. It is slower.*/
    TIC;tic;
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	stream_t temp;
	toc22("init gpu");
    }
    NCPU=NGPU;
    NGPU_MAX=NGPU;
    return NULL;
}
void gpu_mvm_daemon(int port){
    info2("Starting MVM daemon at port %d\n", port);
    pthread_create(&thread_init, NULL, gpu_mvm_gpu_init, NULL);
    redirect();
    listen_port(port, respond, 0, NULL, 1);
}
