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
#include "../lib/sys/scheduler_server.h"
#include <sys/file.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
}
pthread_t thread_init;
typedef struct mvm_g_mul_t{
    smat *g;/*full g*/
    scell *ac;/*a for each gpu.*/
    smat *a;/*final a*/
}mvm_g_mul_t;

typedef struct mvm_t{
    mvm_g_mul_t data;
    thread_t *mvm_g_mul;
    thread_t *mvm_a_cp;
    thread_t *mvm_a_sum;
}mvm_t;
mvm_t *mvm_data=NULL;

int mp_count;/*number multi-processors on each GPU*/
int mp_core;/*number of cuda cores per multi-processor*/
enum{
    GPU_MVM_ZERO,
    GPU_MVM_M,
    GPU_MVM_G,
    GPU_MVM_A,
};
int ndone;
#define WRITE_INTARR(p,n)						\
    if((ndone=write(sock_mvm, p, sizeof(int)*n))!=sizeof(int)*n){	\
	perror("write");						\
	error("error writing. want %d write %ld\n", n, ndone/sizeof(int)); \
    }
#define READ_INTARR(p,n)						\
    if((ndone=read(sock_mvm, p, sizeof(int)*n))!=sizeof(int)*n){	\
	perror("read");							\
	error("error reading. want %d got %ld\n", n, ndone/sizeof(int)); \
    }
int nleft;
char *start;

/*Read into double array until we get all*/
#define READ_ARR(p,n,type)			\
    nleft=(n)*sizeof(type);			\
    start=(char*)(p);				\
    do{						\
	int nread=read(sock_mvm, start, nleft);	\
	if(nread<0) {				\
	    error("Read failed\n");		\
	}else{						\
	    nleft-=nread;				\
	    start+=nread;				\
	}					\
    }while(nleft>0);

#define WRITE_ARR(p,n,type)						\
    if((ndone=write(sock_mvm, p, sizeof(type)*(n)))!=sizeof(type)*(n)){	\
	perror("write");						\
	error("error writing. want %ld wrote %d\n",n*sizeof(type), ndone); \
    }

#define N_CMD 4
int sock_mvm;

__global__ void mvm_g_mul_do(float *mvm, float *a, const float *g, int nact, int ng){
    extern __shared__ float acc[];
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	acc[threadIdx.x]=0;
	for(int ig=0; ig<ng; ig++){
	    register float *mvmi=mvm+nact*ig+iact;
	    acc[threadIdx.x]+=(*mvmi)*g[ig];
	}
	a[iact]=acc[threadIdx.x];
    }
}
__global__ void mvm_g_mulacc_do(float *mvm, float *a, const float *g, int nact, int ng){
    extern __shared__ float acc2[];
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	acc2[threadIdx.x]=0;
	for(int ig=0; ig<ng; ig++){
	    register float *mvmi=mvm+nact*ig+iact;
	    acc2[threadIdx.x]+=*mvmi*g[ig];
	}
	a[iact]+=acc2[threadIdx.x];
    }
}
/* multiply mvm against g for each GPU. each gpu get an equal slice of the whole g.*/
static void mvm_g_mul(thread_t *info){
    //TIC;tic;
    int igpu=info->ithread;
    gpu_set(igpu);
    int icol=info->start;
    int k=info->end-info->start;
    int m=cudata->mvm_m->nx;

    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    smat *g=data->g;
    cp2gpu(&cudata->mvm_g, g->p+icol, k, cudata->mvm_stream[0]);

    int neach=(m+mp_count-1)/mp_count;
    mvm_g_mul_do<<<mp_count, neach, sizeof(float)*neach, *cudata->mvm_stream>>>
	(cudata->mvm_m->p+m*icol, cudata->mvm_a->p, cudata->mvm_g->p, m, k);

    //toc22("gpu queue");
    cp2cpu(&data->ac->p[igpu], cudata->mvm_a, *cudata->mvm_stream);
    cudaStreamSynchronize(*cudata->mvm_stream);
    //toc22("gpu sync");
}
/*copy data from each gpu to cpu*/
static void mvm_a_cp(thread_t *info){
    int igpu=info->ithread;
    gpu_set(igpu);
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    cp2cpu(&data->ac->p[igpu], cudata->mvm_a, *cudata->mvm_stream);
    curzero(cudata->mvm_a, *cudata->mvm_stream);
    cudaStreamSynchronize(*cudata->mvm_stream);
}
/*sum the DM commands from different GPUs together.*/
static void mvm_a_sum(thread_t *info){
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    float *restrict pout=data->a->p;
    for(int igpu=0; igpu<NGPU; igpu++){
	const float *restrict pin=data->ac->p[igpu]->p;
	for(int i=info->start; i<info->end; i++){
	    pout[i]+=pin[i];
	}
    }
}
static void mvm_copy_m(thread_t *info){
    smat *mvm=(smat *)info->data;
    int igpu=info->ithread;
    gpu_set(igpu);
    pthread_mutex_init(&cudata->mvm_mutex, NULL);
    cudata->mvm_stream=new stream_t;
    cp2gpu(&cudata->mvm_m, mvm, *cudata->mvm_stream);
    int nact=mvm->nx;
    cudata->mvm_a=curnew(nact, 1);
    curzero(cudata->mvm_a, *cudata->mvm_stream);
    cudaStreamSynchronize(*cudata->mvm_stream);
}
static int respond(int sock){
    TIC;
    double tk0=tic;
    sock_mvm=sock;
    int cmd[N_CMD];
    READ_INTARR(cmd, N_CMD);
    double tim_cmd=toc3; tic;
    static double tim_gfirst=0;
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	info("Receiving mvm\n");
	smat *mvm=snew(cmd[2], cmd[3]);
	READ_ARR(mvm->p, cmd[2]*cmd[3],float);
	toc22("Read mvm");tic;
	int nact=cmd[2];
	int ngtot=cmd[3];
	thread_t info[NGPU];
	thread_prep(info, 0, NGPU, NGPU, mvm_copy_m, mvm);
	pthread_join(thread_init, NULL);
	CALL_THREAD(info, NGPU, 1);
	toc22("copy mvm to gpu");
	sfree(mvm);
	mvm_data=new mvm_t;
	mvm_data->data.g=snew(ngtot, 1);
	mvm_data->data.ac=scellnew(NGPU, 1);
	mvm_data->data.a=snew(nact, 1);
	for(int ig=0; ig<NGPU; ig++){
	    mvm_data->data.ac->p[ig]=snew(nact, 1);
	}
	mvm_data->mvm_g_mul=new thread_t[NGPU];
	thread_prep(mvm_data->mvm_g_mul, 0, ngtot, NGPU, mvm_g_mul, &mvm_data->data);
	mvm_data->mvm_a_cp=new thread_t[NGPU];
	thread_prep(mvm_data->mvm_a_cp, 0, NGPU, NGPU, mvm_a_cp, &mvm_data->data);
	mvm_data->mvm_a_sum=new thread_t[NCPU];
	thread_prep(mvm_data->mvm_a_sum, 0, nact, NCPU, mvm_a_sum, &mvm_data->data);
	info2("done");
    }
	break;
    case GPU_MVM_G:{/*maos sends gradients*/
	int icol=cmd[1];/*starting column*/
	const int m=cudata->mvm_m->nx;
	const int k=cmd[2];
	smat *g=mvm_data->data.g;
	READ_ARR(g->p+icol, k, float);
	if(cmd[2]<cudata->mvm_m->ny){//part of grads
	    if(icol==0){
		tim_gfirst=myclockd();
	    }
	    gpu_set(gpu_next());//use next GPUs
	    cp2gpu(&cudata->mvm_g, g->p+icol, k, cudata->mvm_stream[0]);
	    int neach=(m+mp_count-1)/mp_count;
	    mvm_g_mulacc_do<<<mp_count, neach, sizeof(float)*neach, *cudata->mvm_stream>>>
		(cudata->mvm_m->p+m*icol, cudata->mvm_a->p, cudata->mvm_g->p, m, k);
	}else{ //all grads. break to different gpus
	    double tim_g=toc3;tic;
	    double tim_gmul=0, tim_dmcp=0, tim_dmsend=0;
	    CALL_THREAD(mvm_data->mvm_g_mul, NGPU, 1);
	    tim_gmul=toc3; tic;
	    szero(mvm_data->data.a);tic;
	    CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
	    tim_dmcp=toc3; tic;
	    WRITE_ARR(mvm_data->data.a->p, m, float);
	    tim_dmsend=toc3; tic;
	    info2("CMD %2.0f, receive g %4.0f, gmul %4.0f, dmcp %4.0f, dmsend %4.0f, total %5.0f\n", 
		  tim_cmd*1e6, tim_g*1e6, tim_gmul*1e6, tim_dmcp*1e6, tim_dmsend*1e6, (myclockd()-tk0)*1e6);
	}
    }
	break;
    case GPU_MVM_A:{/*maos finish sending gradients. Push dm commands back. only
		      when partial gradient sending each time.*/
	double tim_gsend=myclockd()-tim_gfirst; tic;
	CALL_THREAD(mvm_data->mvm_a_cp, NGPU, 1);
	double tim_dmcp=toc3;tic;
	szero(mvm_data->data.a);
	CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
	double tim_dm=toc3;tic;
	WRITE_ARR(mvm_data->data.a->p, cudata->mvm_m->nx, float);
	info2("CMD %2.0f, receive g %4.0f, dmcp %4.0f, dmsum %4.0f, dmsend %4.0f, total %5.0f\n",
	      tim_cmd*1e6, tim_gsend*1e6, tim_dmcp*1e6, tim_dm*1e6, toc3*1e6, (myclockd()-tim_gfirst)*1e6);
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
    THREAD_POOL_INIT(NCPU);
    /*Creating stream for the first time is slow. So do it here to avoid latency
      later. Do not init stream_t in multithread. It is slower.*/
    TIC;tic;
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	stream_t temp;
	toc22("init gpu");
    }
}
void gpu_mvm_daemon(int port){
    info2("Starting MVM daemon at port %d\n", port);
    pthread_create(&thread_init, NULL, gpu_mvm_gpu_init, NULL);
    listen_port(port, respond, 0, NULL);
}