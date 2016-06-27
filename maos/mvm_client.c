/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <sys/socket.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <errno.h>
static int sock_mvm=-1;
/*Read into double array until we get all*/
#define READ_ARR(p,n,type)				\
    {							\
	int nleft;					\
	char *start;					\
	nleft=(n)*sizeof(type);				\
	start=(char*)(p);				\
	do{						\
	    int nread=read(sock_mvm, start, nleft);	\
	    nleft-=nread;				\
	    start+=nread;				\
	    if(nread<=0){				\
		error("Read failed\n");			\
	    }						\
	}while(nleft>0);				\
    }

#define WRITE_ARR(p,n,type)						\
    {									\
	size_t ndone;							\
	if((ndone=write(sock_mvm, p, sizeof(type)*(n)))!=sizeof(type)*(n)){ \
	    perror("write");						\
	    error("error writing. want %ld wrote %zu\n",(long)(n)*sizeof(type), ndone); \
	}								\
    }
#define READ_CMD(p) READ_ARR(p,N_CMD,int)
#define WRITE_CMD(p) WRITE_ARR(p,N_CMD,int)

/**
   It is called by maos to send mvm matrix to mvm server
*/
void mvm_client_send_m(dmat *mvmd, int ngpu){
    int cmd[N_CMD]={GPU_MVM_M, 0, 0, 0};
    cmd[1]=ngpu;
    cmd[2]=mvmd->nx;
    cmd[3]=mvmd->ny;
    info2("sending mvm %dx%d ...", cmd[2], cmd[3]);
    WRITE_CMD(cmd);
    float *fmvm=malloc(sizeof(float)*(mvmd->nx*mvmd->ny));
    for(int i=0; i<mvmd->nx*mvmd->ny; i++){
	fmvm[i]=(float)(mvmd->p[i]*GSCALE1*ASCALE);
    }
    WRITE_ARR(fmvm, mvmd->nx*mvmd->ny, float);
    free(fmvm);
    info2(" done\n");
}
void mvm_client_init(const char *host, int port, dmat *mvm, int ngpu){
    if((sock_mvm=connect_port(host, port, 1, 1))<0){
	error("Unable to connect to mvm server\n");
    }
    mvm_client_send_m(mvm, ngpu);
}
void mvm_client_recon(int mvmsize, dcell *dm, dcell *grad){
    if(sock_mvm==-1){
	error("please call mvm_client_init first\n");
    }
    /* first move gradients to continuous buffer with correct numerical type.*/
    int nwfs=grad->nx;
    int ngtot=0;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(grad->p[iwfs]){
	    int n=grad->p[iwfs]->nx;
	    ngtot+=n;
	}
    }
    GTYPE *gall=malloc(ngtot*sizeof(GTYPE));
    GTYPE *pgall=gall;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	if(!grad->p[iwfs]) continue;
	int ng=grad->p[iwfs]->nx;
	for(int ig=0; ig<ng; ig++){
	    *(pgall++)=(GTYPE)(grad->p[iwfs]->p[ig]*GSCALE);
	}
    }
    int natot=0;
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    natot+=dm->p[idm]->nx;
	}
    }
    ATYPE *dmall=malloc(sizeof(ATYPE) *natot);
    ATYPE *pdmall=dmall;
    int neach=mvmsize;//testing parms->sim.mvmsize;
    if(neach<=0){//scan different sizes.
	static int neach0=10;
	neach=(neach0+=10);
    }
    if(neach>ngtot) neach=ngtot;
    TIC;double tk0=tic;
    int cmd[N_CMD]={GPU_MVM_G, 0, 0, 1};
    cmd[1]=0;//NGPU
    cmd[2]=neach;//Number of grads each time.
    WRITE_CMD(cmd);
#if 0 //write set of data each time.
    for(int i=0; i<ngtot; i+=neach){
	WRITE_ARR(gall+i, MIN(neach, ngtot-i), GTYPE);
    }
#else //write all the data once
    WRITE_ARR(gall, ngtot, GTYPE);
#endif
    double tim_gsend=toc3;
    //Read computed DM command
    READ_ARR(dmall, natot, ATYPE);
    double tim_aread=toc3;
    //Copy DM command to the right place.
    for(int idm=0; idm<dm->nx; idm++){
	if(dm->p[idm]){
	    int nact=dm->p[idm]->nx;
	    double *restrict pdm=dm->p[idm]->p;
	    for(int i=0; i<nact; i++){
		pdm[i]=(double)(*(pdmall++)*ASCALE1);
	    }
	}
    }
    double tim_acp=toc3; tic;
    info2("k=%d, gsend %3.0f, dmread %4.0f, dmcopy %2.0f, total %4.0f\n", neach,
	  tim_gsend*1e6, tim_aread*1e6, tim_acp*1e6, (myclockd()-tk0)*1e6);
    free(dmall);
    free(gall);
}

void mvm_client_close(void){
    close(sock_mvm);
}
static __attribute__((destructor))void deinit(){
    if(sock_mvm!=-1){
	close(sock_mvm);
    }
}
