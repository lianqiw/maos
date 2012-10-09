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

/**
   Testing network throughput.
*/
#include "../lib/aos.h"

int N=10000;
double *temp=NULL; 
int nrw;
int nleft;
char *start;
#define READ(sock,p,n)						\
    nleft=(n);							\
    start=(char*)(p);						\
    do{								\
	int nread=read(sock, start, nleft);			\
	nleft-=nread;						\
	start+=nread;						\
	if(nread<=0){						\
	    warning("nread=%d, nleft=%d\n", nread, nleft);	\
	    return -1;						\
	}							\
    }while(nleft>0);						

#define WRITE(sock, p, n)					\
    if((nrw=write(sock,p,n)!=n)) error2("read failed\n");

int server(int sock){
    double tim1=myclockd();
    long nr;
    READ(sock, &nr, sizeof(long));
    double tim2=myclockd();
    WRITE(sock, temp, sizeof(double));
    double tim3=myclockd();
    READ(sock, temp, sizeof(double)*nr);
    double tim4=myclockd();
    WRITE(sock, temp, sizeof(double));
    double tim5=myclockd();
    info2("read %5.1f, send %5.1f, read2 %5.1f send2 %5.1f\n", (tim2-tim1)*1e6, (tim3-tim2)*1e6, (tim4-tim3)*1e6, (tim5-tim4)*1e6);
    return 0;
}
int client(const char *hostname, int port, int type){
    int sock;
    int nretry=100;
    if(type==1){
	sock=connect_port(hostname, port, 0, 0);
	if(sock<0) exit(1);
	dmat *timing=dnew(nretry, N);
	for(long iN=1; iN<=N; iN++){
	    for(int i=0; i<nretry; i++){
		double tim1=myclockd();
		WRITE(sock, &iN, sizeof(long));
		double tim2=myclockd();
		READ(sock, temp, sizeof(double));
		double tim3=myclockd();
		WRITE(sock, temp, sizeof(double)*iN);
		double tim4=myclockd();
		READ(sock, temp, sizeof(double));
		double tim5=myclockd();
		if(i==0) info2("N=%ld, send %5.1f, read %5.1f, send2 %5.1f read2 %5.1f\n", iN, (tim2-tim1)*1e6, (tim3-tim2)*1e6, (tim4-tim3)*1e6, (tim5-tim4)*1e6);
		timing->p[i+(iN-1)*nretry]=(tim5-tim3)*1e6;
	    }
	}
	dwrite(timing, "timing_%s", hostname);
    }else{
	sock=socket(AF_INET, SOCK_DGRAM, 0);
    }
    return 0;
}
int main(int argc, char *argv[]){
    if(argc<3){
	error("Usage: ./ethtest client hostname port or ./ethtest server port");
    }
    int type=1;
    if(!temp) temp=malloc(N*sizeof(double));
    if(!strcmp(argv[1], "server")){
	int port=strtol(argv[2], NULL, 10);
	if(argc>3){
	    type=strtol(argv[3], NULL, 10);
	}
	listen_port(port, server, 0, NULL);
    }else{
	if(argc<4){
	    error("Usage: ./ethtest client hostname port");
	}
	const char *host=argv[2];
	int port=strtol(argv[3], NULL, 10);
	if(argc>4){
	    type=strtol(argv[4], NULL, 10);
	}
	client(host, port, type);
    }
}
