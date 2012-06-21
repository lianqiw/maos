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

#include "mvm_client.h"
#include <sys/file.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
int sock_mvm;
int ndone;
#define N_CMD 4
enum{
    GPU_MVM_ZERO,
    GPU_MVM_M,
    GPU_MVM_G,
    GPU_MVM_A,
};
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

void mvm_client_init(const char *host, int port){
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
}

/**
   It is called by maos to send mvm matrix to mvm server
*/
void mvm_client_send_m(dcell *mvm){
    info2("sending mvm ...");
    dmat *mvmd=dcell2m(mvm);
    int cmd[N_CMD]={GPU_MVM_M, 0, 0, 0};
    cmd[2]=mvmd->nx;
    cmd[3]=mvmd->ny;
    WRITE_INTARR(cmd,N_CMD);
    float *fmvm=malloc(sizeof(float)*(mvmd->nx*mvmd->ny));
    for(int i=0; i<mvmd->nx*mvmd->ny; i++){
	fmvm[i]=(float)mvmd->p[i];
    }
    info2("prep");
    WRITE_ARR(fmvm, mvmd->nx*mvmd->ny, float);
    writeflt(fmvm, mvmd->nx, mvmd->ny, "fmvm");
    writedbl(mvmd->p,  mvmd->nx, mvmd->ny, "dmvm");
    dcellwrite(mvm, "cmvm");
    dwrite(mvmd, "dmvm2");
    dfree(mvmd);
    free(fmvm);
    info2("done\n");
}

void mvm_client_recon(const PARMS_T *parms, dcell *dm, dcell *grad){
    /* first move gradients to continuous buffer*/
    int nwfs=grad->nx;
    int ngtot=0;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(grad->p[iwfs]){
	    int n=grad->p[iwfs]->nx;
	    ngtot+=n;
	}
    }
    float *gall=malloc(ngtot*sizeof(float));
    float *pgall=gall;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(!grad->p[iwfs]) continue;
	int ng=grad->p[iwfs]->nx;
	for(int ig=0; ig<ng; ig++){
	    *(pgall++)=(float)grad->p[iwfs]->p[ig];
	}
    }
    int natot=0;
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    natot+=dm->p[idm]->nx;
	}
    }
    float *dmall=malloc(sizeof(float) *natot);
    float *pdmall=dmall;
    static int neach=0;
    if(neach==0){
	neach=parms->sim.mvmsize;
	if(neach==0){
	    neach=ngtot;
	}else{
	    while(ngtot % neach){
		neach++;
	    }
	}
	info("neach = %d\n", neach);
    }
    //neach=ngtot;//temporary
    //neach=ngtot/7;//temporary.

    TIC;double tk0=tic;
    int cmd[N_CMD]={GPU_MVM_G, 0, 0, 1};
    cmd[2]=neach;
    for(int i=0; i<ngtot; i+=neach){
	cmd[1]=i;
	WRITE_INTARR(cmd,N_CMD);
	WRITE_ARR(gall+i, neach, float);
    }
    double tim_gsend=toc3; tic;
    cmd[0]=GPU_MVM_A;
    if(neach!=ngtot){
	WRITE_INTARR(cmd,N_CMD);
    }
    READ_ARR(dmall, natot, float);
    double tim_aread=toc3; tic;
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    int nact=dm->p[idm]->nx;
	    double *restrict pdm=dm->p[idm]->p;
	    for(int i=0; i<nact; i++){
		pdm[i]=(double) *(pdmall++);
	    }
	    //memcpy(dm->p[idm]->p, pdmall, sizeof(double)*nact);
	}
    }
    double tim_acp=toc3; tic;
    info2("send gradiens %3.0f, get dm %4.0f, copy dm %2.0f, total %4.0f\n",
	  tim_gsend*1e6, tim_aread*1e6, tim_acp*1e6, (myclockd()-tk0)*1e6);
    free(dmall);
    free(gall);
}
