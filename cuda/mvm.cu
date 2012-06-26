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
#include "../maos/mvm_client.h"
#include <sys/file.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <errno.h>
}

/*Read into double array until we get all*/
#define READ_ARR(p,n,type)					\
    nleft=(n)*sizeof(type);					\
    start=(char*)(p);						\
    do{								\
	int nread=read(sock_mvm, start, nleft);			\
	nleft-=nread;						\
	start+=nread;						\
	if(nread<=0){						\
	    warning2("nread=%d, nleft=%d\n", nread, nleft);	\
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



/**
 * 2012-06-22: Serious bug found:
 * Initially was using curnew. But later replaced by new while should be using cudaMalloc
 * This caused misterious bugs, such as no answer or even kernel panic.
 */
pthread_t thread_init;
typedef struct mvm_g_mul_t{
    int nact;
    int ngtot;
    int icol;
    int ngpu;
    int k;
    GTYPE *g;/*full g*/
    ATYPE **ac;/*a for each gpu.*/
    ATYPE *a;/*final a*/
}mvm_g_mul_t;

typedef struct mvm_t{
    mvm_g_mul_t data;
    thread_t *mvm_g_mul;
    thread_t *mvm_a_cp;
    thread_t *mvm_a_sum;
    thread_t *mvm_a_cp_sum;
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
/* multiply mvm against g for each GPU. each gpu get an equal slice of the g.*/
static void mvm_g_mul(thread_t *info){
    int igpu=info->ithread;
    gpu_set(igpu);
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    int ki=(data->k+NGPU-1)/NGPU;
    int icol=data->icol+ki*igpu;
    int ki2=data->k+data->icol-icol;
    int k=MIN(ki, ki2);
    int m=cudata->mvm_m->nx;

    cudaMemcpyAsync(cudata->mvm_g+icol, data->g+icol, k*sizeof(GTYPE), 
		    cudaMemcpyHostToDevice, *cudata->mvm_stream);
    int neach=(m+mp_count-1)/mp_count;
    mvm_g_mul_do<<<mp_count, neach, sizeof(float)*neach, *cudata->mvm_stream>>>
	(cudata->mvm_m->p+m*icol, cudata->mvm_a, cudata->mvm_g+icol, m, k);
}
__global__ static void mvm_a_sum_do(ATYPE *restrict a, ATYPE *restrict a2, int nact){
    int iact=threadIdx.x+blockIdx.x*blockDim.x;
    if(iact<nact){
	a[iact]+=a2[iact];
    }
}
/*
  To reduce the data (a) across all GPUs, use GPU reduction method we used in kernels:
  Do then in GPUs.
*/
static void mvm_a_cp_sum(thread_t *info){
    int igpu=info->ithread;
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    int nact=data->nact;
    int neach=(nact+mp_count-1)/mp_count;

    for(int step=(NGPU+1)>>1; step>0; step>>=1){
	if(igpu<step && igpu+step<NGPU){
	    gpu_set(igpu);
	    DO(cudaStreamSynchronize(cudata_all[igpu+step].mvm_stream[0]));
#define COPY_A2 1
#if COPY_A2
	    if(!cudata->mvm_a2[igpu+step]){
		DO(cudaMalloc(&(cudata->mvm_a2[igpu+step]), sizeof(ATYPE)*nact));
	    }
	    DO(cudaMemcpyAsync(cudata_all[igpu].mvm_a2[igpu+step], 
			       cudata_all[igpu+step].mvm_a, 
			       sizeof(ATYPE)*nact, cudaMemcpyDefault, *cudata->mvm_stream));
	    mvm_a_sum_do<<<mp_count, neach, 0, *cudata->mvm_stream>>>
		(cudata->mvm_a, cudata->mvm_a2[igpu+step], data->nact);
#else
	    mvm_a_sum_do<<<mp_count, neach, 0, *cudata->mvm_stream>>>
		(cudata->mvm_a, cudata_all[igpu+step].mvm_a, data->nact);  
#endif
	}
    }
    if(igpu==0){
	gpu_set(igpu);/*gpu to add the data*/
	cudaMemcpyAsync(data->a, cudata->mvm_a, data->nact*sizeof(ATYPE),
			cudaMemcpyDeviceToHost, *cudata->mvm_stream);
	cudaStreamSynchronize(*cudata->mvm_stream);
    }
}
static void mvm_a_zero(thread_t *info){
    int igpu=info->ithread;
    gpu_set(igpu);
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    int nact=data->nact;
    cudaMemsetAsync(cudata->mvm_a, 0, nact*sizeof(ATYPE), *cudata->mvm_stream);
}
/*copy data from each gpu to cpu and clear accumulation*/
static void mvm_a_cp(thread_t *info){
    int igpu=info->ithread;
    gpu_set(igpu);
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    int nact=data->nact;
    cudaMemcpyAsync(data->ac[igpu], cudata->mvm_a, nact*sizeof(ATYPE),
		    cudaMemcpyDeviceToHost, *cudata->mvm_stream);
    cudaStreamSynchronize(*cudata->mvm_stream);
    cudaMemsetAsync(cudata->mvm_a, 0, nact*sizeof(ATYPE), *cudata->mvm_stream);
}
/*sum the DM commands from different GPUs together.*/
static void mvm_a_sum(thread_t *info){
    mvm_g_mul_t *data=(mvm_g_mul_t *)info->data;
    ATYPE *restrict pout=data->a;
    for(int igpu=0; igpu<NGPU; igpu++){
	const ATYPE *restrict pin=data->ac[igpu];
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
}
static void mvm_data_free(void){
    cudaFreeHost(mvm_data->data.g);
    cudaFreeHost(mvm_data->data.a);
    for(int i=0; i<mvm_data->data.ngpu; i++){
	cudaFreeHost(mvm_data->data.ac[i]);
    }
    free(mvm_data->data.ac);
    free(mvm_data->mvm_g_mul);
    free(mvm_data->mvm_a_cp);
    free(mvm_data->mvm_a_sum);
    free(mvm_data->mvm_a_cp_sum);
    free(mvm_data);
}
static int respond(int sock){
    TIC;
    double tk0=tic;
    sock_mvm=sock;
    int cmd[N_CMD];
    READ_CMD(cmd);
    double tim_cmd=toc3; tic;
    static int ksave=0;
    static double tim_gfirst=0;
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	int nact=cmd[2];
	int ngtot=cmd[3];
	info("Receiving mvm %dx%d\n", nact, ngtot);
	smat *mvm=snew(nact, ngtot);
	READ_ARR(mvm->p, (nact*ngtot), float);
	pthread_join(thread_init, NULL);
	toc22("Read mvm");tic;
	thread_t info[NGPU];
	thread_prep(info, 0, NGPU, NGPU, mvm_copy_m, mvm);
	CALL_THREAD(info, NGPU, 0);
	toc22("copy mvm to gpu");
	if(mvm_data){
	    mvm_data_free();
	}
	mvm_data=(mvm_t*)calloc(1, sizeof(mvm_t));
	mvm_data->data.nact=mvm->nx;
	mvm_data->data.ngtot=mvm->ny;
	mvm_data->data.ngpu=NGPU;
	mvm_data->data.ac=new ATYPE*[NGPU];
	cudaMallocHost(&mvm_data->data.g, sizeof(GTYPE)*ngtot);
	cudaMallocHost(&mvm_data->data.a, sizeof(ATYPE)*nact);
	memset(mvm_data->data.a, 0., nact*sizeof(ATYPE));
	for(int ig=0; ig<NGPU; ig++){
	    cudaMallocHost(&mvm_data->data.ac[ig], sizeof(ATYPE)*nact);
	}
	mvm_data->mvm_g_mul=new thread_t[NGPU];
	thread_prep(mvm_data->mvm_g_mul, 0, ngtot, NGPU, mvm_g_mul, &mvm_data->data);
	mvm_data->mvm_a_cp=new thread_t[NGPU];
	mvm_data->mvm_a_sum=new thread_t[NCPU];
#define USE_CP_SUM 0
#if USE_CP_SUM
	thread_prep(mvm_data->mvm_a_cp, 0, NGPU, NGPU, mvm_a_cp_sum, &mvm_data->data);
	thread_prep(mvm_data->mvm_a_sum, 0, nact, NCPU, mvm_a_zero, &mvm_data->data);
#else
	thread_prep(mvm_data->mvm_a_cp, 0, NGPU, NGPU, mvm_a_cp, &mvm_data->data);
	thread_prep(mvm_data->mvm_a_sum, 0, nact, NCPU, mvm_a_sum, &mvm_data->data);
#endif
	info2("done");
	sfree(mvm);
    }
	break;
    case GPU_MVM_G:{/*maos sends gradients*/
	int icol=cmd[1];/*starting column*/
	static int count;
	const int nact=mvm_data->data.nact;
	const int k=cmd[2];
	if(icol==0){
	    tim_gfirst=myclockd();
	    ksave=k;
	    count=k;
	}else{
	    count+=k;
	}

	READ_ARR(mvm_data->data.g+icol, k, GTYPE);
	if(cmd[2]<1800){//part of grads
	    gpu_set(gpu_next());//use next GPUs
	    cudaMemcpyAsync(cudata->mvm_g+icol, mvm_data->data.g+icol, k*sizeof(GTYPE), 
			    cudaMemcpyHostToDevice, cudata->mvm_stream[0]);
	    int naeach=(nact+mp_count-1)/mp_count;
	    mvm_g_mul_do<<<mp_count, naeach, sizeof(float)*naeach, *cudata->mvm_stream>>>
		(cudata->mvm_m->p+nact*icol, cudata->mvm_a, cudata->mvm_g+icol, nact, k);
	}else{ //Send to different gpus
	    mvm_data->data.icol=icol;
	    mvm_data->data.k=k;
	    CALL_THREAD(mvm_data->mvm_g_mul, NGPU, 0);
	}
	if(count==mvm_data->data.ngtot){/*gradient send is over*/
#if USE_CP_SUM
	    double tim_gsend=myclockd()-tim_gfirst; tic;
	    CALL_THREAD(mvm_data->mvm_a_cp, NGPU, 0);
	    double tim_dmcp=toc3;tic;
	    double tim_dmsum=toc3;tic;
	    WRITE_ARR(mvm_data->data.a, mvm_data->data.nact, ATYPE);
	    CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
#else

	    double tim_gsend=myclockd()-tim_gfirst; tic;
	    CALL_THREAD(mvm_data->mvm_a_cp, NGPU, 0);
	    double tim_dmcp=toc3;tic;
	    CALL_THREAD(mvm_data->mvm_a_sum, NCPU, 1);
	    double tim_dmsum=toc3;tic;
	    WRITE_ARR(mvm_data->data.a, mvm_data->data.nact, ATYPE);
	    memset(mvm_data->data.a, 0., mvm_data->data.nact*sizeof(ATYPE));
#endif
	    info2("k=%4d CMD %1.0f, receive g %4.0f, dmcp %4.0f, dmsum %4.0f, dmsend %4.0f, total %5.0f\n", ksave,
		  tim_cmd*1e6, tim_gsend*1e6, tim_dmcp*1e6, tim_dmsum*1e6, toc3*1e6, (myclockd()-tim_gfirst)*1e6);
	}
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
    return NULL;
}
void gpu_mvm_daemon(int port){
    info2("Starting MVM daemon at port %d\n", port);
    pthread_create(&thread_init, NULL, gpu_mvm_gpu_init, NULL);
    listen_port(port, respond, 0, NULL);
}
