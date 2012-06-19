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
}

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

#define N_CMD 2
int sock_mvm;
int port_mvm;
static int respond(int sock){
    sock_mvm=sock;
    int cmd[N_CMD];
    int dim[2];
    READ_INTARR(cmd, N_CMD);
    switch(cmd[0]){
    case GPU_MVM_M:{/*maos sends M matrix*/
	info("Receiving mvm\n");
	READ_INTARR(dim, 2);
	smat *mvm=snew(dim[0], dim[1]);
	READ_ARR(mvm->p, dim[0]*dim[1],float);
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    assert(!cudata->mvm_m);
	    cp2gpu(&cudata->mvm_m, mvm);
	    cudata->mvm_a=curnew(mvm->nx, 1);
	    cudata->mvm_stream=new stream_t;
	    pthread_mutex_init(&cudata->mvm_mutex, NULL);
	    info2("Initializing gpu %d\n", igpu);
	}
	info2("done");
	sfree(mvm);
    }
	break;
    case GPU_MVM_G:{/*maos sends gradients*/
	int icol=cmd[1];/*starting column*/
	READ_INTARR(dim, 2);
	smat *g=snew(dim[0], dim[1]);
	READ_ARR(g->p, dim[0]*dim[1], float);
	gpu_set(gpu_next());/*use next GPUs*/
	cp2gpu(&cudata->mvm_g, g, cudata->mvm_stream[0]);
	float alpha=1.;
	float beta=1.;
	int m=cudata->mvm_m->nx;
	int n=1;
	int k=dim[0];
	DO(cublasSgemm(*cudata->mvm_stream, CUBLAS_OP_N, CUBLAS_OP_N,
		       m, n, k, &beta, 
		       cudata->mvm_m->p+m*icol, m,
		       cudata->mvm_g->p, k,
		       &alpha, cudata->mvm_a->p, m));
	sfree(g);
    }
	break;
    case GPU_MVM_A:{/*maos finish sending gradients. Push commands back*/
	dmat *a=NULL;
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(gpu_next());
	    add2cpu(&a, cudata->mvm_a, *cudata->mvm_stream, &cudata->mvm_mutex);
	    curzero(cudata->mvm_a, *cudata->mvm_stream);
	}
	dim[0]=a->nx;
	dim[1]=a->ny;
	WRITE_INTARR(dim, 2);
	WRITE_ARR(a->p, dim[0]*dim[1], double);
	dfree(a);
    }
	break;
    }
    return 0;
}

/** It is called by maos to launch mvm server and forks to background afterwards.*/
void gpu_mvm_init(int port){
    port=12347;
    port_mvm=port;
    if((sock_mvm=connect_port("localhost", port, 1, 0))<0){
	error("Unable to connect to mvm server\n");
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
    int gpu[1]={0};
    gpu_init(NULL, 0);
    listen_port(port, respond, 0, NULL);
}
/**
   It is called by maos to send mvm matrix to mvm server
*/
void gpu_mvm_send_m(dcell *mvm){
    info("sending mvm ...");
    dmat *mvmd=dcell2m(mvm);
    int cmd[N_CMD]={GPU_MVM_M, 0};
    int dim[2];
    dim[0]=mvmd->nx;
    dim[1]=mvmd->ny;
    WRITE_INTARR(cmd,N_CMD);
    info2("cmd written");
    WRITE_INTARR(dim,2);
    info2("dim written");
    float *fmvm=new float[mvmd->nx*mvmd->ny];
    for(int i=0; i<mvmd->nx*mvmd->ny; i++){
	fmvm[i]=(float)mvmd->p[i];
    }
    WRITE_ARR(fmvm, dim[0]*dim[1], float);
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
    double *dmall=new double[natot];
    double *pdmall=dmall;
    static int neach=50;
    neach=ngtot/NGPU;
    while(ngtot%neach){
	neach++;
    }
    int cmd[N_CMD]={GPU_MVM_G, 0};
    int dim[2];
    dim[0]=neach;
    dim[1]=1;
    TIC;tic;
    for(int i=0; i<ngtot; i+=neach){
	cmd[1]=i;
	WRITE_INTARR(cmd,N_CMD);
	WRITE_INTARR(dim, 2);
	WRITE_ARR(gall+i, neach, float);
    }
    toc2("send grad");
    cmd[0]=GPU_MVM_A;
    WRITE_INTARR(cmd,N_CMD);
    READ_INTARR(dim,2);
    READ_ARR(dmall, dim[0]*dim[1], double);
    toc2("read dm");
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    int nact=dm->p[idm]->nx;
	    memcpy(dm->p[idm]->p, pdmall, sizeof(double)*nact);
	    pdmall+=nact;
	}
    }
    toc2("copy dm");
    delete[] dmall;
    delete[] gall;
}
