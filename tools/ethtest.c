/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

/**
   \file ethtest.c
   Testing network throughput.
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 
#endif
#include <time.h>
#include <errno.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include "../lib/aos.h"

int nrep;
int nbuf=65536;

char *buf1, *buf2;
int nstep=1;
int nmax=1000;
/**
   The following two routines overwrite the data when new data comes in. Just for testing jitter.
 */
/*Write long messages with smaller buffer*/
int stwrite3(int sfd, const void *p, size_t len, size_t buflen){
    if(buflen>len) buflen=len;
    size_t nwrite;
    long left=len;//do not use size_t which is unsigned
#ifdef __linux__
    int is_sock=issock(sfd);
#endif
    do{
#ifdef __linux__
	if(is_sock){
	    nwrite=send(sfd, p, len, MSG_NOSIGNAL);
	}else
#endif
	    nwrite=write(sfd, p, buflen);
	
	left-=nwrite;
    }while(nwrite>0 && left>0);
    return left?-1:0;
}
/*Read long messages with smaller buffer*/
int stread3(int sfd, void *p, size_t len, size_t buflen){
    if(buflen>len) buflen=len;
    size_t nread;
    long left=len;
    do{
	nread=read(sfd, p, buflen);
	left-=nread;
    }while(nread>0 && left>0);
    return left?-1:0;
}
int server(int sock){
    if(stread(sock, &nstep, sizeof(int))
       || stread(sock, &nmax, sizeof(int))
       || stread(sock, &nrep, sizeof(int))){
	return -1;
    }
    info2("nstep=%d, nmax=%d, nrep=%d\n", nstep, nmax, nrep);
    buf2=malloc(nmax*nstep);
    //warm up
    for(int i=0; i<10; i++){
	if(stread(sock, buf2, nmax*nstep) || stwrite(sock, buf2, 64)){
	    close(sock);
	    return -1;
	}
    }
    for(int len=1; len<=nmax; len++){
	for(int irep=0; irep<nrep; irep++){
	    if(stread(sock, buf2, len*nstep) || stwrite(sock, buf2, 64)){
		close(sock);
		return -1;
	    }
	}
    }
    return -1;
}

int client(const char *hostname, int port){
    int sock=connect_port(hostname, port, 0, 1);
    if(sock<0 || stwriteint(sock, nstep) 
       || stwriteint(sock, nmax)
       || stwriteint(sock, nrep)) {
	warning("Unable to connecto to %s\n", hostname);
	close(sock); return 1;
    }
    buf1=malloc(nmax*nstep);
    for(int i=0;i<10;i++){//warm up
	stwrite(sock, buf1, nmax*nstep);
	stread(sock, buf1, 64);
	usleep(500);
    }
    double tim1, tim2, tim3;
    dmat *timing=dnew(nrep, nmax);
    dmat *timing2=dnew(nrep, nmax);
    for(int len=1; len<=nmax; len++){
	info2("len=%d\n", len);
	for(int irep=0; irep<nrep; irep++){
	    usleep(500);
	    tim1=myclockd();
	    stwrite(sock, buf1, len*nstep);
	    tim2=myclockd();
	    stread(sock, buf1, 64);
	    tim3=myclockd();
	    timing->p[irep+(len-1)*nrep]=tim3-tim1;
	    timing2->p[irep+(len-1)*nrep]=tim2-tim1;
	}
    }
    close(sock);
    dwrite(timing, "timing_%d_%d", nstep, nmax);
    dwrite(timing2, "timing2_%d_%d", nstep, nmax);
    info("done\n");
    return 0;
}

//server for mvmfull_iwfs
int mvm_server(int sock){
    int nact, nsa, sastep, pixpsa, nstep0;
    int cmd[6];
    if(streadintarr(sock, cmd, 6)){
	return -1;
    }
    nact=cmd[0];
    nsa=cmd[1];
    sastep=cmd[2];
    pixpsa=cmd[3];
    nstep=cmd[4];
    nstep0=cmd[5];
    info("nact=%d, nsa=%d, sastep=%d, pixpsa=%d, nstep=%d\n",
	 nact, nsa, sastep, pixpsa, nstep);
    smat *pix=snew(pixpsa, nsa);
    smat *dmres=snew(nact, 1);
    rand_t rseed;
    seed_rand(&rseed, 1);
    srandn(pix, 20, &rseed);
    int ready;
    struct timespec ct;
    int readtime_ns=500000;//500 micro-second read out time. 
    int frametime_ns=1250000;//frame time.
    int nsend=((nsa+sastep-1)/sastep);//segments sent along read out
    int int1_ns=readtime_ns/nsend;//interval between segment sending
    int int2_ns=frametime_ns-readtime_ns;//interval after last segment.
    streadint(sock, &ready); //wait for client to be ready.
    clock_gettime(CLOCK_MONOTONIC, &ct);
    //double t0=ct.tv_sec+ct.tv_nsec*1e-9;
    TIC;
    for(int istep=-nstep0; istep<nstep; istep++){
	//info2("\rSend trigger ");
	tk=(double)ct.tv_sec+(double)ct.tv_nsec*1.e-9;//start of the frame.
	for(int isa=0; isa<nsa; isa+=sastep){
	    clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ct, NULL);
	    int nleft=(nsa-isa)<sastep?(nsa-isa):sastep;
	    if(stwrite(sock, pix->p+isa*pixpsa, 2*nleft*pixpsa)){//2 byte data.
		warning("failed: %s\n", strerror(errno));
		return -1;
	    }
	    ct.tv_nsec+=int1_ns;
	    while(ct.tv_nsec>=1000000000){
		ct.tv_nsec-=1000000000;
		ct.tv_sec++;
	    }
	    
	}
	if(stread(sock, dmres->p, sizeof(float)*nact)){
	    warning("read dmres failed: %s\n", strerror(errno));
	    return -1;
	}
	ready=(int)(toc3*1e6);//mvm is finished.
	if(stwriteint(sock, ready)){
	    warning("write ready failed: %s\n", strerror(errno));
	    return -1;
	}
	//set next frame start time.
	ct.tv_nsec+=int2_ns;
	while(ct.tv_nsec>=1000000000){
	    ct.tv_nsec-=1000000000;
	    ct.tv_sec++;
	}
	if((istep & 0b11111111) == 0b11111111){
	    info2("%d %d us\n", istep, ready);
	}

    }
    info2("\n");
    return -1;
}

int main(int argc, char *argv[]){
    if(argc<3){
	error2("Usage: \n\t./ethtest client hostname port \n\t./ethtest server port\n\t./ethtest mvm_server port\n");
	_Exit(1);
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
	info2("Set priority to -20\n");
	setpriority(PRIO_PROCESS, getpid(), -20);
	struct sched_param param;
	sched_getparam(getpid(), &param);
	param.sched_priority=sched_get_priority_max(SCHED_FIFO);

	sched_setscheduler(getpid(), SCHED_FIFO, &param);
    }
    char *NBUF=getenv("NBUF");
    if(NBUF){
	nbuf=strtol(NBUF, NULL, 10);
	info("nbuf=%d\n", nbuf);
    }

    if(!strcmp(argv[1], "server")){
	int port=strtol(argv[2], NULL, 10);
	listen_port(port, NULL, server, 0, NULL, 1);
    }else if(!strcmp(argv[1],"mvm_server")){
	int port=strtol(argv[2], NULL, 10);
	listen_port(port, NULL, mvm_server, 0, NULL, 1);	
    }else{
	if(argc<4){
	    error("Usage: ./ethtest client hostname port [nstep] [nmax] [nrep]");
	}
	const char *host=argv[2];
	int port=strtol(argv[3], NULL, 10);
	if(argc>4){
	    nstep=strtol(argv[4], NULL, 10);
	}
	if(argc>5){
	    nmax=strtol(argv[5], NULL, 10);
	}
	if(argc>6){
	    nrep=strtol(argv[6], NULL, 10);
	}
	client(host, port);
    }
}
