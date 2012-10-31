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
   \file sock.c
   
   Routines handle socket connection

   use ioctl with SIOCETHTOOL detect whether line is connected.

   Note for TCP socket:

   Applications that require lower latency on every packet sent should be run on
   sockets with TCP_NODELAY enabled. It can be enabled through the
   setsockopt command with the sockets API.
   
   For this to be used effectively, applications must avoid doing small,
   logically related buffer writes. Because TCP_NODELAY is enabled, these
   small writes will make TCP send these multiple buffers as individual
   packets, which can result in poor overall performance. 


   Changelog:
   2012-10-25:
   New creation. Migrated from scheduler_client.c and scheduler_server.c
   Make scheduler multi-threaded so that select() is not blocked due to pending handling requests. 

*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h> 

#include <string.h>
#include <errno.h>
#include <netdb.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#ifndef SOL_TCP
#define SOL_TCP IPPROTO_TCP
#endif
#include "misc.h"
#include "common.h"
#include "sock.h"
#include "thread.h"
/**
   When the keepalive flag is on, the socket will receive notice when the
   connection to remote socket is disrupted.

  \todo Find keepalive options in mac.  */
static void socket_tcp_keepalive(int sock){
    int keeplive=1;
#ifdef __linux__
    int keepidle =1;/*second before try to probe */
    int keepintvl=1;/*wait this seconds before repeat */
    int keepcnt  =2;/*repeat before declare dead */
#endif
    if(!setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, &keeplive, sizeof(int))
#if defined(__linux__) 
       && !setsockopt(sock, SOL_TCP, TCP_KEEPCNT, &keepcnt, sizeof(int))
       && !setsockopt(sock, SOL_TCP, TCP_KEEPIDLE, &keepidle, sizeof(int))
       && !setsockopt(sock, SOL_TCP, TCP_KEEPINTVL, &keepintvl, sizeof(int))
#endif
       ){
    }else{
	perror("setsockopt");
	warning("Keepalive failed. sock=%d\n", sock);
    }
}

#if defined(__INTEL_COMPILER)
/*with htons defined in glibc 2.4, intel compiler complains
  about conversion from in to uint16. THis is an ugly workaround*/
#undef htons
#define htons myhtons
static inline uint16_t myhtons(uint16_t port){
    uint16_t ans;
#if __BYTE_ORDER == __BIG_ENDIAN
    ans=(port);
#else
    ans=(unsigned short int)
	((((port) >> 8) & 0xff) | (unsigned short int)(((port) & 0xff) << 8));
#endif
    return ans;
}
#endif

/**
   Tune socket for more real time.
*/
static void socket_tcp_nodelay(int sock){
    int one=1;
    setsockopt(sock, SOL_TCP, TCP_NODELAY, &one, sizeof(one));
}
/**
   make a server port and bind to localhost on all addresses
*/
int bind_socket (uint16_t port, int type){
    int sock=-1;
    struct sockaddr_in name;
    
    /* Create the socket. */
    switch(type){
    case 1:{//TCP
	sock = socket(PF_INET, SOCK_STREAM, 0);//tcp
	socket_tcp_keepalive(sock);
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
    
    cloexec(sock);
    setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,NULL,sizeof(int));
    /* Give the socket a name. */
    name.sin_family = AF_INET;
    name.sin_port = htons(port);
    name.sin_addr.s_addr = htonl(INADDR_ANY);
    int count=0;
    while(bind(sock,(struct sockaddr *)&name, sizeof (name))<0){
	info3("errno=%d. port=%d,sock=%d: ",errno,port,sock);
	perror ("bind");
	sleep(10);
	count++;
	if(count>100){
	    error("Failed to bind to port %d\n",port);
	}
    }
    info2("binded to port %hd at sock %d\n",port,sock);
    return sock;
}

static int quit_listen=0;
static void signal_handler(int sig){
    quit_listen=1;
}
/**
   Open a port and listen to it. Calls respond(sock) to handle data. If
   timeout_fun is not NULL, it will be called when 1) connection is
   establed/handled, 2) every timeout_sec if timeout_sec>0. This function only
   return at error.
 */
void listen_port(uint16_t port, int (*responder)(int), double timeout_sec, void (*timeout_fun)(), int nodelay){
    register_signal_handler(signal_handler);

    fd_set read_fd_set;
    fd_set active_fd_set;
    struct sockaddr_in clientname;
    int sock = bind_socket (port, 1);//has to be tcp.
    if(nodelay){
	socket_tcp_nodelay(sock);
    }
    if (listen (sock, 1) < 0){
	perror("listen");
	exit(EXIT_FAILURE);
    }
    FD_ZERO (&active_fd_set);
    FD_SET (sock, &active_fd_set);

    while(!quit_listen){
	struct timeval timeout;
	timeout.tv_sec=timeout_sec;
	timeout.tv_usec=0;
	struct timeval *timeout2;
	timeout2=timeout_sec>0?&timeout:NULL;

	read_fd_set = active_fd_set;
	if(select(FD_SETSIZE, &read_fd_set, NULL, NULL, timeout2)<0){
	    perror("select");
	    if(errno!=EINTR){
		warning("select failed\n");
		continue;
	    }else if(errno==EBADF){
		break;//bad file descriptor
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
			warning("accept failed\n");
			break;
		    }
		    info("port %d is connected\n", port2);
		    cloexec(port2);
		    FD_SET(port2, &active_fd_set);
		}else{
		    /* Data arriving on an already-connected socket. Call responder to handle.
		       On return:
		       negative value: Close read of socket. 
		       -1: also close the socket.
		     */
		    int ans=responder(i);
		    if(ans<0){
			shutdown(i, SHUT_RD);
			FD_CLR(i, &active_fd_set);
			if(ans==-1){
			    warning2("close port %d\n", i);
			    close(i);
			}else{
			    warning2("shutdown port %d for reading\n", i);
			}
		    }
		}
	    }
	}
	if(timeout_fun){
	    timeout_fun();
	}
    }
    /* Error happened. We close all connections and this server socket.*/
    warning("listen_port exited\n");
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

int init_sockaddr (struct sockaddr_in *name, const char *hostname, uint16_t port){
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

/**
  Connect to a host at port.
  mode=0: read/write
  mode=1: read only by the client. the server won't read
  mode=2: write only by the client. the server won't write
*/
int connect_port(const char *hostname, int port, int block, int nodelay){
    int sock;
    int count=0;
    struct sockaddr_in servername;
    do{
	count++;
	/* Create the socket. */
	sock = socket (PF_INET, SOCK_STREAM, 0);
	cloexec(sock);
	socket_tcp_keepalive(sock);
	if(nodelay){
	    socket_tcp_nodelay(sock);
	}
	if(!block){
	    fcntl(sock, F_SETFD, O_NONBLOCK);
	}
	/* Give the socket a name. */
	if(init_sockaddr(&servername, hostname, port) || 
	   connect(sock, (struct sockaddr *)&servername, sizeof (servername))<0){
	    close(sock);
	    if(!block){
		return -1;
	    }
	    sleep(4);
	}else{
	    return sock;
	}
    }while(count<10);
    return -1; 
}
