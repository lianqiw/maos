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
    start=(char*)p;				\
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
int port_mvm;

__global__ void mvm_g_mul_do(float *mvm, float *a, const float *g, int nact, int ng){
    extern __shared__ float acc[];
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	acc[threadIdx.x]=0;
	for(int ig=0; ig<ng; ig++){
	    register float *mvmi=mvm+nact*ig+iact;
	    acc[threadIdx.x]+=*mvmi*g[ig];
	}
	a[iact]=acc[threadIdx.x];
    }
}
/* multiply mvm against g for each GPU. each gpu get an equal slice of the whole g.*/
static void mvm_g_mul(thread_t *info){
    TIC;tic;
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

    toc2("gpu %d queue", igpu);
    cp2cpu(&data->ac->p[igpu], cudata->mvm_a, *cudata->mvm_stream);
    cudaStreamSynchronize(*cudata->mvm_stream);
    toc2("gpu %d sync", igpu);
}
/*copy data from each gpu to cpu*/
static void mvm_a_cp(thread_t *info){
    int igpu=info->ithread;
    gpu_set(igpu);
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    cp2cpu(&data->ac->p[igpu], cudata->mvm_a, *cudata->mvm_stream);
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
static int respond(int sock){
    TIC;tic;
    sock_mvm=sock;
    int cmd[N_CMD];
    READ_INTARR(cmd, N_CMD);
    toc2("read cmd");tic;
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	info("Receiving mvm\n");
	smat *mvm=snew(cmd[2], cmd[3]);
	READ_ARR(mvm->p, cmd[2]*cmd[3],float);
	int nact=cmd[2];
	int ngtot=cmd[3];
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    assert(!cudata->mvm_m);
	    cp2gpu(&cudata->mvm_m, mvm);
	    cudata->mvm_a=curnew(nact, 1);
	    cudata->mvm_stream=new stream_t;
	    pthread_mutex_init(&cudata->mvm_mutex, NULL);
	    info2("Initializing gpu %d\n", igpu);
	}
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
	toc2("read grad"); tic;
	const float alpha=1.;
	const float beta=1.;
	const int m=cudata->mvm_m->nx;
	const int n=1;
	const int k=cmd[2];
	smat *g=mvm_data->data.g;
	READ_ARR(g->p+icol, k, float);
	if(cmd[2]<cudata->mvm_m->ny){//part of grads
	    gpu_set(gpu_next());//use next GPUs
	    cp2gpu(&cudata->mvm_g, g->p+icol, k, cudata->mvm_stream[0]);
	    DO(cublasSgemm(*cudata->mvm_stream, CUBLAS_OP_N, CUBLAS_OP_N,
			   m, n, k, &beta, 
			   cudata->mvm_m->p+m*icol, m,
			   cudata->mvm_g->p, k,
			   &alpha, cudata->mvm_a->p, m));
	}else{ //all grads. break to different gpus
	    assert(icol==0);
	    CALL_THREAD(mvm_data->mvm_g_mul, NGPU, 1);
	    szero(mvm_data->data.a);
	    CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
	    WRITE_ARR(mvm_data->data.a->p, m, float);
	}
	toc("MVM tot");
    }
	break;
    case GPU_MVM_A:{/*maos finish sending gradients. Push dm commands back. only
		      when partial gradient sending each time.*/
	CALL_THREAD(mvm_data->mvm_a_cp, NGPU, 1);
	CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
	toc2("to cpu");tic;
	WRITE_ARR(mvm_data->data.a->p, cudata->mvm_m->nx, float);
	toc2("write back");
    }
	break;
    }
    return 0;
}

/** It is called by maos to launch mvm server and forks to background afterwards.*/
void gpu_mvm_init(const char *host, int port){
    port=12347;
    port_mvm=port;
    if((sock_mvm=connect_port(host, port, 1, 0))<0){
	error("Unable to connect to mvm server\n");
    }
    {
	/*Applications that require lower latency on every packet sent should be
	  run on sockets with TCP_NODELAY enabled. It can be enabled through the
	  setsockopt command with the sockets API.  

	  For this to be used effectively, applications must avoid doing small,
	  logically related buffer writes. Because TCP_NODELAY is enabled, these
	  small writes will make TCP send these multiple buffers as individual
	  packets, which can result in poor overall performance.  */
	int one=1;
	//setsockopt(sock_mvm, SOL_TCP, TCP_NODELAY|TCP_QUICKACK|TCP_CORK, &one, sizeof(one));
	setsockopt(sock_mvm, SOL_TCP, TCP_NODELAY, &one, sizeof(one));
    }
    //fcntl(sock_mvm, F_SETFD, O_NONBLOCK);
    /*First we fork*/
    /*
    gpu_cleanup();
    pid_t pid=fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if(pid>0){
	warning("Launching mvm and connect to it\n");
	sleep(1);
	port_mvm=port;
	if((sock_mvm=connect_port("localhost", port, 1, 0))<0){
	    error("Unable to connect to mvm server\n");
	}
	return;//parent continue to run maos
    }else{
	// child of maos
	redirect();
	gpu_cleanup();
	gpu_init(NULL, 0);
	listen_port(port, respond, 0, NULL);
    }*/
}
void gpu_mvm_daemon(int port){
    info2("Starting MVM daemon at port %d\n", port);
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
    THREAD_POOL_INIT(NGPU);
    listen_port(port, respond, 0, NULL);
}
/**
   It is called by maos to send mvm matrix to mvm server
*/
void gpu_mvm_send_m(dcell *mvm){
    info("sending mvm ...");
    dmat *mvmd=dcell2m(mvm);
    int cmd[N_CMD]={GPU_MVM_M, 0, 0, 0};
    cmd[2]=mvmd->nx;
    cmd[3]=mvmd->ny;
    WRITE_INTARR(cmd,N_CMD);
    info2("cmd written");
    float *fmvm=new float[mvmd->nx*mvmd->ny];
    for(int i=0; i<mvmd->nx*mvmd->ny; i++){
	fmvm[i]=(float)mvmd->p[i];
    }
    WRITE_ARR(fmvm, mvmd->nx*mvmd->ny, float);
    dfree(mvmd);
    delete [] fmvm;
    info2("done");
}

void gpu_mvm_recon(dcell *dm, dcell *grad){
    assert(dm);
    /* first move gradients to continuous buffer*/
    int nwfs=grad->nx;
    int ngtot=0;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(grad->p[iwfs]){
	    int n=grad->p[iwfs]->nx;
	    ngtot+=n;
	}
    }
    float *gall=new float[ngtot];
    float *pgall=gall;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(!grad->p[iwfs]) continue;
	int ng=grad->p[iwfs]->nx;
	for(int ig=0; ig<ng; ig++){
	    pgall[ig]=(float)grad->p[iwfs]->p[ig];
	}
	pgall+=grad->p[iwfs]->nx;
    }
    int natot=0;
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    natot+=dm->p[idm]->nx;
	}
    }
    float *dmall=new float[natot];
    float *pdmall=dmall;
    static int neach=50;
    neach=ngtot;//temporary
    while(ngtot%neach){
	neach++;
    }
    int cmd[N_CMD]={GPU_MVM_G, 0, 0, 1};
    cmd[2]=neach;
    TIC;tic;
    for(int i=0; i<ngtot; i+=neach){
	cmd[1]=i;
	WRITE_INTARR(cmd,N_CMD);
	WRITE_ARR(gall+i, neach, float);
    }
    toc2("send grad");tic;
    cmd[0]=GPU_MVM_A;
    if(neach!=ngtot){
	WRITE_INTARR(cmd,N_CMD);
    }
    READ_ARR(dmall, natot, float);
    toc2("read dm");tic;
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    int nact=dm->p[idm]->nx;
	    double *restrict pdm=dm->p[idm]->p;
	    for(int i=0; i<nact; i++){
		pdm[i]=(double)pdmall[i];
	    }
	    //memcpy(dm->p[idm]->p, pdmall, sizeof(double)*nact);
	    pdmall+=nact;
	}
    }
    toc2("copy dm");
    delete[] dmall;
    delete[] gall;
}
