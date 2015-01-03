/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <stdio.h>
#include <stdlib.h>
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <limits.h>
#include <string.h>
#define error(A...) ({							\
	    fprintf(stderr, "\033[01;31mFatal error\033[00;00m\t");	\
	    fprintf(stderr, A);						\
	    raise(SIGTERM);})
#define info(A...) fprintf(stderr, A)
#define warning(A...) ({					\
	    fprintf(stderr,"\033[00;31m");			\
	    fprintf(stderr,A);fprintf(stderr,"\033[00;00m"); }) 

/**
   make a server port and bind to localhost on all addresses
*/
int bind_socket (uint16_t port, int type){
    int sock;
    struct sockaddr_in name;
    
    /* Create the socket. */
    switch(type){
    case 1:{
	sock = socket(PF_INET, SOCK_STREAM, 0);//tcp
	/*Applications that require lower latency on every packet sent should be
	  run on sockets with TCP_NODELAY enabled. It can be enabled through the
	  setsockopt command with the sockets API.  

	  For this to be used effectively, applications must avoid doing small,
	  logically related buffer writes. Because TCP_NODELAY is enabled, these
	  small writes will make TCP send these multiple buffers as individual
	  packets, which can result in poor overall performance.  */
	int one=1;
	setsockopt(sock, SOL_TCP, TCP_NODELAY, &one, sizeof(one));
    }
	break;
    case 2:
	sock = socket(PF_INET, SOCK_DGRAM, 0); //udp
	break;
    default:
	error("Invalid type");
    }
    if (sock < 0){
	perror ("socket");
	exit (EXIT_FAILURE);
    }
    
    setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,NULL,sizeof(int));
    /* Give the socket a name. */
    name.sin_family = AF_INET;
    name.sin_port = htons(port);
    name.sin_addr.s_addr = htonl(INADDR_ANY);
    int count=0;
    while(bind(sock,(struct sockaddr *)&name, sizeof (name))<0){
	info("errno=%d. port=%d,sock=%d: ",errno,port,sock);
	perror ("bind");
	sleep(10);
	count++;
	if(count>100){
	    error("Failed to bind to port %d\n",port);
	}
    }
    info("binded to port %hd at sock %d\n",port,sock);
    return sock;
}

void listen_port(uint16_t port, int (*respond)(int), double timeout_sec, void (*timeout_fun)()){
    fd_set read_fd_set;
    fd_set active_fd_set;
    struct sockaddr_in clientname;
    int sock = bind_socket (port, 1);
    if (listen (sock, 1) < 0){
	perror("listen");
	exit(EXIT_FAILURE);
    }
    FD_ZERO (&active_fd_set);
    FD_SET (sock, &active_fd_set);

    while(1){
	struct timeval timeout;
	timeout.tv_sec=timeout_sec;
	timeout.tv_usec=0;
	struct timeval *timeout2;
	timeout2=timeout_sec>0?&timeout:NULL;

	read_fd_set = active_fd_set;
	if(select(FD_SETSIZE, &read_fd_set, NULL, NULL, timeout2)<0){
	    perror("select");
	    if(errno!=EINTR){
		break;
	    }else{
		continue;
	    }
	}
	for(int i=0; i<FD_SETSIZE; i++){
	    if(FD_ISSET(i, &read_fd_set)){
		if(i==sock){
		    /* Connection request on original socket. */
		    socklen_t size=sizeof(struct sockaddr_in);;
    		    int port2=accept(sock, (struct sockaddr*)&clientname, &size);
		    if(port2<0){
			perror("accept");
			break;
		    }
		    FD_SET(port2, &active_fd_set);
		}else{
		    /* Data arriving on an already-connected socket. */
		    if(respond(i)<0){
			warning("Close port %d\n", i);
			shutdown(i, SHUT_RD);
			FD_CLR(i, &active_fd_set);
			close(i);
		    }
		}
	    }
	}
	if(timeout_fun){
	    timeout_fun();
	}
    }
    /* Error happened. We close all connections and this server socket.*/
    shutdown(sock, SHUT_RDWR);
    for(int i=0; i<FD_SETSIZE; i++){
	if(FD_ISSET(i, &active_fd_set)){
	    shutdown(i, SHUT_RDWR);
	    usleep(100);
	    close(i);
	    FD_CLR(i, &active_fd_set);
	}
    }
    close(sock);
    usleep(100);
    _Exit(1);
}
int init_sockaddr (struct sockaddr_in *name,
		   const char *hostname, uint16_t port){
    struct hostent *hostinfo;
    
    name->sin_family = AF_INET;
    name->sin_port = htons(port);
    hostinfo = gethostbyname (hostname);
    if (hostinfo == NULL){
	return -1;
    }else{
	struct in_addr *addr = (struct in_addr *) hostinfo->h_addr_list[0];
	if(addr){
	    name->sin_addr = *addr;
	    return 0;
	}else{
	    /*warning("h_addr_list is NULL for host %s\n", hostname); */
	    return -1;
	}
    }
}

/*
  Connect to a host at port.
  mode=0: read/write
  mode=1: read only by the client. the server won't read
  mode=2: write only by the client. the server won't write
*/
int connect_port(const char *hostname, int port, int block, int mode){
    int sock;
    int count=0;
    struct sockaddr_in servername;
    do{
	/* Create the socket. */
	sock = socket (PF_INET, SOCK_STREAM, 0);
	if (sock < 0) {
	    perror ("socket (scheduler)");
	    return sock;
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
	    setsockopt(sock, SOL_TCP, TCP_NODELAY, &one, sizeof(one));
	}
	/* Give the socket a name. */
	init_sockaddr(&servername, hostname, port);
	if(!block){
	    fcntl(sock, F_SETFD, O_NONBLOCK);
	}
	if(connect(sock, (struct sockaddr *)&servername, sizeof (servername))<0){
	    perror("connect");
	    close(sock);
	    if(!block){
		return -1;
	    }
	    sleep(4);
	    count++;
	}else{
	    return sock;
	}
    }while(count<10);
    return -1; 

}

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
    if((nrw=write(sock,p,n)!=n)) error("read failed\n");
static double myclockd(void){
    struct timeval tk;
    gettimeofday(&tk,NULL);
    return (double)tk.tv_sec+(double)tk.tv_usec*1e-6;
}
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
    //info("read %5.1f, send %5.1f, read2 %5.1f send2 %5.1f\n", (tim2-tim1)*1e6, (tim3-tim2)*1e6, (tim4-tim3)*1e6, (tim5-tim4)*1e6);
    return 0;
}
int client(const char *hostname, int port, int type){
    int sock;
    int nretry=10;
    if(type==1){
	sock=connect_port(hostname, port, 0, 0);
	if(sock<0) exit(1);
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
		info("N=%ld, send %5.1f, read %5.1f, send2 %5.1f read2 %5.1f\n", iN, (tim2-tim1)*1e6, (tim3-tim2)*1e6, (tim4-tim3)*1e6, (tim5-tim4)*1e6);
	    }
	}
    }else{
	sock=socket(AF_INET, SOCK_DGRAM, 0);
    }
    return 0;
}
int main(int argc, char *argv[]){
    if(argc<3){
	error("Usage: ./ethtest client servername port or ./ethtest server port");
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
