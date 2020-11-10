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

/**
   \file ethtest.c
   Testing network throughput.
*/
#include <errno.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include "../lib/aos.h"

int nbuf=65536;

char* buf1, * buf2;
/**
   The following two routines overwrite the data when new data comes in. Just for testing jitter.
 */
/*Write long messages with smaller buffer*/
int stwrite3(int sfd, const void* p, size_t len, size_t buflen){
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
		} else
#endif
			nwrite=write(sfd, p, buflen);

		left-=nwrite;
	} while(nwrite>0&&left>0);
	return left?-1:0;
}
/*Read long messages with smaller buffer*/
int stread3(int sfd, void* p, size_t len, size_t buflen){
	if(buflen>len) buflen=len;
	size_t nread;
	long left=len;
	do{
		nread=read(sfd, p, buflen);
		left-=nread;
	} while(nread>0&&left>0);
	return left?-1:0;
}
int server(int sock){
	int nstep;
	int nmin;
	int nmax;
	int nrep;
	if(stread(sock, &nstep, sizeof(int))
		||stread(sock, &nmin, sizeof(int))
		||stread(sock, &nmax, sizeof(int))
		||stread(sock, &nrep, sizeof(int))){
		return -1;
	}
	info("nstep=%d, nmin=%d, nmax=%d, nrep=%d\n", nstep, nmin, nmax, nrep);
	buf2=(char*)malloc(nmax*nstep);
	//warm up
	for(int i=0; i<10; i++){
		if(stread(sock, buf2, nmax)||stwrite(sock, buf2, 64)){
			close(sock);
			return -1;
		}
	}
	for(int len=nmin; len<=nmax; len+=nstep){
		for(int irep=0; irep<nrep; irep++){
			if(stread(sock, buf2, len)||stwrite(sock, buf2, 64)){
				close(sock);
				return -1;
			}
		}
	}
	return -1;
}

int client(const char* hostname, int port, int nmin, int nmax, int nstep, int nrep){
	int sock=connect_port(hostname, port, 0, 1);
	if(sock<0||stwriteint(sock, nstep)
		||stwriteint(sock, nmin)
		||stwriteint(sock, nmax)
		||stwriteint(sock, nrep)){
		warning("Unable to connecto to %s\n", hostname);
		close(sock); return 1;
	}
	buf1=(char*)malloc(nmax*nstep);
	for(int i=0;i<10;i++){//warm up
		stwrite(sock, buf1, nmax);
		stread(sock, buf1, 64);
		usleep(500);
	}
	double tim1, tim2, tim3;
	int nlen=(nmax-nmin+nstep)/nstep;
	dmat* timing=dnew(nrep, nlen);
	dmat* timing2=dnew(nrep, nlen);
	int ilen=-1;
	for(int len=nmin; len<=nmax; len+=nstep){
		ilen++;
		info("len=%d\n", len);
		for(int irep=0; irep<nrep; irep++){
			if(irep%800==0){
				info("irep=%d of %d\n", irep, nrep);
			}
			usleep(500);
			tim1=myclockd();
			stwrite(sock, buf1, len);
			tim2=myclockd();
			stread(sock, buf1, 64);
			tim3=myclockd();
			timing->p[irep+ilen*nrep]=tim3-tim1;
			timing2->p[irep+ilen*nrep]=tim2-tim1;
		}
	}
	close(sock);
	writebin(timing, "pix_timing_%s_%d_%d_%d", HOST, nmin, nmax, nstep);
	writebin(timing2, "pix_timing2_%s_%d_%d_%d", HOST, nmin, nmax, nstep);
	dbg("done\n");
	return 0;
}
//server for mvmfull_real
int mvm_server(int sock){
	int cmd[7];
	if(streadintarr(sock, cmd, 7)){
		return -1;
	}
	int nact=cmd[0];
	int nsa=cmd[1];
	int sastep=cmd[2];
	int totpix=cmd[3];
	int pixpsa=totpix;
	int nstep=cmd[4];
	int nstep0=cmd[5];
	int type=cmd[6];
	dbg("type=%d, nact=%d, nsa=%d, sastep=%d, %s=%d, nstep=%d\n", type,
		nact, nsa, sastep, type==1?"pixpsa":"totpix", totpix, nstep);
	int* saind=NULL;
	if(type==1){//mvmfull_iwfs
		totpix=pixpsa*nsa;
	} else{//mvmfull_real
		saind=mymalloc((nsa+1), int);
		if(streadintarr(sock, saind, nsa+1)){
			return -1;
		}
	}
	short* pix=mymalloc(totpix, short);
	if(type==1){
		rand_t rseed;
		seed_rand(&rseed, 1);
		for(int i=0; i<totpix; i++){
			pix[i]=(short)randu(&rseed);
		}
	} else{
		if(stread(sock, pix, totpix*sizeof(short))){
			return -1;
		}
	}
	smat* dmres=snew(nact, 1);
	int ready;
	streadint(sock, &ready); //wait for client to be ready.
#if __linux__
	struct timespec ct;
	clock_gettime(CLOCK_MONOTONIC, &ct);
	int readtime_ns=500000;//500 micro-second read out time. 
	int frametime_ns=1250000;//frame time.
	int nsend=((nsa+sastep-1)/sastep);//number of segments sent along read out
	int int1_ns=readtime_ns/nsend;//interval between segment sending
	int int2_ns=frametime_ns-readtime_ns;//interval after last segment.
#endif
	TIC;tic;
	for(int istep=-nstep0; istep<nstep; istep++){
	//info("\rSend trigger ");
#if __linux__
	//scheduled start time of the frame.
		double tk0=(double)ct.tv_sec+(double)ct.tv_nsec*1.e-9;
#else
		tic;
#endif
		for(int isa=0; isa<nsa; isa+=sastep){
#if __linux__
			if(clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ct, NULL)){
				warning("clock_nanosleep is interrupted\n");
			}
			if(isa==0){
				tic;
			}
#endif
			int nleft;
			if(type==1){
				nleft=((nsa-isa)<sastep?(nsa-isa):sastep)*pixpsa;
			} else{
				if(nsa<isa+sastep){//terminate
					nleft=totpix-saind[isa];
				} else{
					nleft=saind[isa+sastep]-saind[isa];
				}
			}
			if(stwrite(sock, pix+(type==1?pixpsa*isa:saind[isa]), 2*nleft)){//2 byte data.
				warning("failed: %s\n", strerror(errno));
				return -1;
			}
#if __linux__
			ct.tv_nsec+=int1_ns;
			while(ct.tv_nsec>=1000000000){
				ct.tv_nsec-=1000000000;
				ct.tv_sec++;
			}
#endif
		}
		if(stread(sock, dmres->p, sizeof(float)*nact)){
			warning("read dmres failed: %s\n", strerror(errno));
			return -1;
		}
		ready=(int)(toc3*1e6);//mvm is finished.
#if __linux__
		if(nstep<100){
			dbg("tk=%.6f tic=%.6f, toc=%.6f, ready=%.6f\n", tk0, tk, myclockd(), ready*1e-6);
		}
#endif
		if(stwriteint(sock, ready)){
			warning("write ready failed: %s\n", strerror(errno));
			return -1;
		}
		//set next frame start time.
#if __linux__
		ct.tv_nsec+=int2_ns;
		while(ct.tv_nsec>=1000000000){
			ct.tv_nsec-=1000000000;
			ct.tv_sec++;
		}
#endif
		if((istep&0xFF)==0xFF){
			info("%d %d us.\n", istep, ready);
		}

	}
	info("\n");
	return -1;
}
int main(int argc, char* argv[]){
	if(argc<3){
		error("Usage: \n\t./ethtest client hostname port \n\t./ethtest server port\n\t./ethtest mvm_server port\n");
		_Exit(1);
	}
	set_realtime(0, -20);


	char* NBUF=getenv("NBUF");
	if(NBUF){
		nbuf=strtol(NBUF, NULL, 10);
		dbg("nbuf=%d\n", nbuf);
	}
	char* ip=0;
	if(argc>3){
		ip=argv[3];
	}

	if(!strcmp(argv[1], "server")){
		int port=strtol(argv[2], NULL, 10);
		listen_port(port, ip, server, 0, NULL, 1);
	} else if(!strcmp(argv[1], "mvm_server")){
		int port=strtol(argv[2], NULL, 10);
		listen_port(port, ip, mvm_server, 0, NULL, 1);
	} else{
		if(argc<4){
			error("Usage: ./ethtest client hostname port [nrep] [nmin] [nmax] [nstep]");
		}
		const char* host=argv[2];
		int port=strtol(argv[3], NULL, 10);
		int nstep=100;
		int nmin=409584;
		int nmax=409584;
		int nrep=1000;
		if(argc>4){
			nrep=strtol(argv[4], NULL, 10);
		}
		if(argc>5){
			nmin=strtol(argv[5], NULL, 10);
		}
		if(argc>6){
			nmax=strtol(argv[6], NULL, 10);
		}
		if(argc>7){
			nstep=strtol(argv[7], NULL, 10);
		}

		client(host, port, nmin, nmax, nstep, nrep);
	}
}
