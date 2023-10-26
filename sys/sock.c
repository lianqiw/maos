/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   Routines to establish socket connection

   use ioctl with SIOCETHTOOL detect whether line is connected.

   Note for TCP socket:

   Applications that require lower latency on every packet sent should be run on
   sockets with TCP_NODELAY enabled. It can be enabled through the
   setsockopt command with the sockets API.

   For this to be used effectively, applications must avoid doing small,
   logically related buffer writes. Because TCP_NODELAY is enabled, these
   small writes will make TCP send these multiple buffers as individual
   packets, which can result in poor overall performance.

   Note on TCP server/client connection 3-step establishment process:
   0) Server creates a socket and bind() on a port, and then starts listen() on the port
   1) Client connects to the server by calling connect(), which send a SYN to the server.
   2) The server receives the SYN. It sends SYN+ACK to the client.
   3) The client receives SYN+ACK. It sends back ACK, and then moves to established state
   *) Server receives ACK and moves to estibalished state.

   Note on TCP peer to peer simultaneous open process:
   1) peer 1 and peer 2 send SYN to each other
   2) Each received SYN and send back ACK
   3) Each received ACK and move to estibalished stat.

   Note on TCP 4-step shutdown process:
   1) peer 1 send a FIN using shutdown(socket, SHUT_WR) or close() to active close
   2) peer 2 recieve the FIN and send back ACK. It can still write the data for peer 1 to use.
   3) peer 2 send a FIN using shutdown(socket, SHUT_WR) or close() to passive close
   4) peer 1 receives the FIN and (ACK before or after the FIN) and send back ACK.
   *) peer 2 receive ACK and now connection is closed.
   peer 1 now is in TIME_WAIT state (always) to avoid delayed packets confusing the system.

   In order for the server to close the port and reopen it immediately, it has
   to be peer 2 instead of peer 1. This requires it to first notify the client
   to shutdown the connection.

   shutdown(socket, SHUR_RD) doesn't send anything to its peer. It just blocking reading the socket.
   close() destroys the socket in addition to call shutdown on it.

   Changelog:
   2012-10-25:
   New creation. Migrated from scheduler_client.c and scheduler_server.c
   Make scheduler multi-threaded so that select() is not blocked due to pending handling requests.

*/
#include <unistd.h>
#include <fcntl.h> 
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/un.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>

#ifndef SOL_TCP
#define SOL_TCP IPPROTO_TCP
#endif
#include "misc.h"
#include "common.h"
#include "sock.h"
#include "thread.h"
#include "sockio.h"
#include "scheduler_client.h"
/**
   When the keepalive flag is on, the socket will receive notice when the
   connection to remote socket is disrupted.

  \todo Find keepalive options in mac.  */
static void socket_tcp_keepalive(int sock){
	if(sock==-1) return;
	int keeplive=1;
#ifdef __linux__
	int keepidle=1;/*second before try to probe */
	int keepintvl=1;/*wait this seconds before repeat */
	int keepcnt=2;/*repeat before declare dead */
#endif
	if(!setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, &keeplive, sizeof(int))
#if defined(__linux__) 
		&&!setsockopt(sock, SOL_TCP, TCP_KEEPCNT, &keepcnt, sizeof(int))
		&&!setsockopt(sock, SOL_TCP, TCP_KEEPIDLE, &keepidle, sizeof(int))
		&&!setsockopt(sock, SOL_TCP, TCP_KEEPINTVL, &keepintvl, sizeof(int))
#endif
		){
		//success
	} else{
		warning("Keepalive failed. sock=%d: %s\n", sock, strerror(errno));
	}
}
static void socket_reuse_addr(int sock){
	if(sock==-1) return;
	const int one=1;
	if(setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(int))){
		perror("socket_reuse_addr");
	}
}
int socket_nopipe(int sock){
	if(sock==-1) return -1;
#ifdef SO_NOSIGPIPE
	const int one=1;
	if(setsockopt(sock, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(int))==-1) return -1;
#else
#if ! defined(__linux__)
	if(signal(SIGPIPE, SIG_IGN)==SIG_ERR) return -1;
#endif
	(void)sock;
#endif
	return 0;
}
/**
   Set socket to unblocking mode
*/
int socket_block(int sock, int block){
	if(sock==-1) return -1;
	int arg=fcntl(sock, F_GETFD);
	if(arg==-1){
		perror("fcntl");
		return -1;
	}
	if(block){
		arg&=~O_NONBLOCK;
	} else{
		arg|=O_NONBLOCK;
	}
	if(fcntl(sock, F_SETFD, arg)==-1){
		perror("fcntl");
		return -1;
	}
	return 0;
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
		((((port)>>8)&0xff)|(unsigned short int)(((port)&0xff)<<8));
#endif
	return ans;
}
#endif

/**
   Tune socket for more real time.
*/
static void socket_tcp_nodelay(int sock){
	if(sock==-1) return;
	int one=1;
	setsockopt(sock, SOL_TCP, TCP_NODELAY, &one, sizeof(one));
}
/**
 * Set recv timeout in seconds
 * */
int socket_recv_timeout(int sock, double sec){
	if(sock==-1) return -1;
	struct timeval val={sec,0};
	int ans;
	if((ans=setsockopt(sock, SOL_SOCKET, SO_RCVTIMEO, &val, sizeof(val)))<0){
		warning("socket_recv_timeout: setsockopt failed\n");
	}
	return ans;
}
/**
 * Set recv timeout in seconds
 * */
int socket_send_timeout(int sock, double sec){
	if(sock==-1) return -1;
	struct timeval val={sec,0};
	int ans;
	if((ans=setsockopt(sock, SOL_SOCKET, SO_SNDTIMEO, &val, sizeof(val)))<0){
		warning("socket_send_timeout: setsockopt failed\n");
	}
	return ans;
}
/**
   make a server port and bind to sockpath. AF_UNIX.
 */
static int bind_socket_unix(char* sockpath){
	int sock=socket(AF_UNIX, SOCK_STREAM, 0);
	socket_nopipe(sock);
	struct sockaddr_un addr={0};
	addr.sun_family=AF_UNIX;
	strncpy(addr.sun_path, sockpath, sizeof(addr.sun_path)-1);
	if(bind(sock, (struct sockaddr*)&addr, sizeof(struct sockaddr_un))==-1){
		perror("bind_socket_unix");
		warning_time("bind to %s failed\n", sockpath);
		close(sock);
		sock=-1;
	} else{
		dbg_time("binded to %s at sock %d\n", sockpath, sock);
	}
	return sock;
}
/**
   make a server port and bind to localhost on all addresses.
   protocl is usually SOCK_DGRAM or SOCK_STREAM.
   ip is usually zero.
   if port is 0, will bind an assigned port
*/
int bind_socket(int protocol, char* ip, uint16_t port){
	struct sockaddr_in name;
	/* Create the socket. */
	int sock=socket(PF_INET, protocol, 0);//tcp
	if(sock<0){
		perror("bind_socket");
		exit(EXIT_FAILURE);
	}
	if(protocol==SOCK_STREAM&&port){
		socket_tcp_keepalive(sock);
		socket_reuse_addr(sock);
	}
	socket_nopipe(sock);
	cloexec(sock);

	/* Give the socket a name. */
	name.sin_family=AF_INET;
	name.sin_port=htons(port);
	if(ip){
		name.sin_addr.s_addr=inet_addr(ip);
	} else{
		name.sin_addr.s_addr=htonl(INADDR_ANY);
	}
	int count=0;
	while(bind(sock, (struct sockaddr*)&name, sizeof(name))<0){
		dbg_time("errno=%d. port=%d,sock=%d: %s. Sleep 10 seconds and then retry.", errno, port, sock, strerror(errno));
		sleep(10);
		count++;
		if(count>100){
			warning("Failed to bind to port %d\n", port);
			close(sock);
			sock=-1;
		}
	}
	if(sock!=-1){
		dbg_time("binded to port %hd at socket %d\n", port, sock);
	}
	return sock;
}
/**
 * determine the socket port
 * */
int socket_port(int sock){
	struct sockaddr_in name;
	socklen_t addrlen=sizeof(name);
	if(!getsockname(sock, (struct sockaddr*)&name, &addrlen)){
		return ntohs(name.sin_port);
	} else{
		return -1;
	}
}

/**
 * determine the socket peer IP address
 * */
in_addr_t socket_peer(int sock){
	struct sockaddr_in name;
	socklen_t addrlen=sizeof(name);
	if(!getpeername(sock, (struct sockaddr*)&name, &addrlen)){
		return name.sin_addr.s_addr;
	} else{
		return 0;
	}
}
/**
 * convert in_addr_t to host name for print out
 * */
const char* addr2name(in_addr_t s_addr){
	struct in_addr addr;
	addr.s_addr=s_addr;
	return (const char*)inet_ntoa(addr);
}

static volatile int quit_listen=0;
static int listen_signal_handler(int sig){
	/*quit listening upon signal and do clean up.*/
	dbg_time("scheduler caught: %s(%d)\n", strsignal(sig), sig);
	quit_listen=1;
	return 1;
}
/**
   Open a port and listen to it. Calls respond(sock) to handle data. If
   timeout_fun is not NULL, it will be called when 1) connection is
   establed/handled, 2) every timeout_sec if timeout_sec>0. This function only
   return at error.
   if localpath is not null and start with /, open the local unix port also
   if localpath is not null and contains ip address, only bind to that ip address
 */
void listen_port(uint16_t port, char* localpath, int (*responder)(int),
	double timeout_sec, void (*timeout_fun)(), int nodelay){

	fd_set read_fd_set;
	fd_set active_fd_set;
	FD_ZERO(&active_fd_set);
	char* ip=0;
	int sock_local=-1;
	int sock=-1;
	if(port){
		if((sock=bind_socket(SOCK_STREAM, ip, port))==-1){
			warning("bind_socket failed\n");
			return;
		}
		if(nodelay){//turn off tcp caching.
			dbg("Turn on tcp nodelay\n");
			socket_tcp_nodelay(sock);
		}
		if(listen(sock, 10)<0){
			perror("listen (sock)");
			exit(EXIT_FAILURE);
		}
		FD_SET(sock, &active_fd_set);
		socket_recv_timeout(sock, 60);
		socket_send_timeout(sock, 60);
	}
	if(localpath){
		if(localpath[0]=='/'){
			//also bind to AF_UNIX
			remove(localpath);//if tcp port bind is success, no other process is running
			sock_local=bind_socket_unix(localpath);
			if(sock_local==-1){
				dbg_time("bind to %s failed\n", localpath);
			} else{
				if(!listen(sock_local, 10)){
					FD_SET(sock_local, &active_fd_set);
					socket_send_timeout(sock_local, 60);
					socket_recv_timeout(sock_local, 60);
				} else{
					perror("listen (sock_local)");
					close(sock_local);
					sock_local=-1;
				}
			}
		} else{
			ip=localpath;
		}
	}

	register_signal_handler(listen_signal_handler);
	int nlisten=2;
	while(quit_listen!=2&&nlisten){
		//int new_connection=0;
		if(quit_listen==1){
			/*
			  shutdown(sock, SHUT_WR) sends a FIN to the peer, therefore
			  initialize a active close and the port will remain in TIME_WAIT state.

			  shutdown(sock, SHUT_RD) sends nother, just mark the scoket as not readable.

			  accept() create a new socket on the existing port. A socket is identified by client ip/port and server ip/port.

			  It is the port that remain in TIME_WAIT state.
			 */
			//shutdown(sock, SHUT_RDWR);
			//shutdown(sock_local, SHUT_RDWR);
			/*Notice existing client to shutdown*/
			dbg_time("terminate pending: ask all clients (count=%d) to close.\n", nlisten);
			static double quit_time=0;
			if(!quit_time) quit_time=myclockd();
			//stop listening.
			FD_CLR(sock, &active_fd_set);
			close(sock);
			if(sock_local!=-1){
				FD_CLR(sock_local, &active_fd_set);
				close(sock_local);
			}

			//We ask the client to actively close the connection
			for(int i=0; i<FD_SETSIZE; i++){
				if(FD_ISSET(i, &active_fd_set)){
					int cmd[3]={-1, 0, 0};
					stwrite(i, cmd, 3*sizeof(int));
				}
			}

			if(myclockd()>quit_time+5){//wait for a few seconds before force terminate.
				dbg_time("force terminate while %d clients remaining.\n", nlisten);
				quit_listen=2;//force quit.
			}
			//don't break. Listen for connection close events.
		}

		struct timeval timeout;
		timeout.tv_sec=(time_t)timeout_sec;
		timeout.tv_usec=(timeout_sec-timeout.tv_sec)*1e6;
		read_fd_set=active_fd_set;
		int navail=select(FD_SETSIZE, &read_fd_set, NULL, NULL, timeout_sec>0?&timeout:0);
		if(navail<0){//select failed
			if(errno==EINTR){
				dbg_time("select failed: %s\n", strerror(errno));
				continue;
			} else if(errno==EBADF){
				warning_time("bad file descriptor: %s\n", strerror(errno));
				break;//bad file descriptor
			} else{
				warning_time("unknown error: %s\n", strerror(errno));
				break;
			}
		} else if(navail>0){
			for(int i=0; i<FD_SETSIZE; i++){
				if(FD_ISSET(i, &read_fd_set)){
					if(i==sock||i==sock_local){// Connection request on original socket. 
						union{
							struct sockaddr_un un;
							struct sockaddr_in in;
						} client;
						socklen_t size=sizeof(client);
						int sock2=accept(i, (struct sockaddr*)&client, &size);
						if(sock2<0){
							warning_time("accept failed: %s. \n", strerror(errno));
						} else{
							dbg_time("socket %d is connected\n", sock2);
							FD_SET(sock2, &active_fd_set);
							//socket_recv_timeout(sock2, 5);//do not set recv timeout. The socket may be passed to draw() that does not use select.
							socket_send_timeout(sock2, 60);
							//new_connection=1;
						}
					} else{
						/* Data arriving on an already-connected socket. Call responder to handle.
						   On return:
						   negative value: Close read of socket.
						   -1: also close the socket.
						 */
						int ans=responder(i);
						if(ans<0){
							FD_CLR(i, &active_fd_set);
							close(i);
							dbg_time("close socket %d, ans=%d\n", i, ans);
							//new_connection=-1;
						}
					}
				}
			}
		}else if(timeout_fun){
			//only run timeout_fun when timeout actually happens.
			//if((timeout_sec==0||navail==0)&&!new_connection&&timeout_fun){
			timeout_fun();
		}
		nlisten=0;
		for(int i=0; i<FD_SETSIZE; i++){
			if(FD_ISSET(i, &active_fd_set)){
				nlisten++;
			}
		}
		if(!nlisten){
			dbg_time("all clients have closed.\n");
		}
	}
	/* Error happened. We close all connections and this server socket.*/
	if(sock_local!=-1){
		close(sock_local);
		FD_CLR(sock_local, &active_fd_set);
	}
	if(sock!=-1){
		close(sock);
		FD_CLR(sock, &active_fd_set);
	}
	dbg_time("listen_port exited\n");
	for(int i=0; i<FD_SETSIZE; i++){
		if(FD_ISSET(i, &active_fd_set)){
			dbg_time("sock %d is still connected\n", i);
			close(i);
			FD_CLR(i, &active_fd_set);
		}
	}
	//sync();
}

/**
  Connect to a host at port. If hostname already includes port after :, port is ignored.
  hostname may start with / to indicate a local UNIX port
*/
int connect_port(const char* hostname,/**<The hostname can be just name or name:port.*/
	int port,            /**<The port if hostname does not include it*/
	int block,           /**<Do we block until connection is established*/
	int nodelay			 /**<set TCP nodelay flag*/
	){
	int sock=-1;
	if(hostname[0]=='/'){//connect locally so we can pass fd.
		sock=socket(PF_UNIX, SOCK_STREAM, 0);
		socket_nopipe(sock);
		socket_send_timeout(sock, 60);
		struct sockaddr_un addr={0};
		addr.sun_family=AF_UNIX;
		strncpy(addr.sun_path, hostname, sizeof(addr.sun_path)-1);
		if(connect(sock, (struct sockaddr*)&addr, sizeof(struct sockaddr_un))<0){
			warning("connect locally (%s) failed: %s. \n", hostname, strerror(errno));
			close(sock);
			sock=-1;
		}
		return sock;
		//hostname="localhost";
	} else{
		//if(sock!=-1) return sock;
		struct sockaddr_in servername;
		for(int count=0; count<25; count++){
			sock=socket(PF_INET, SOCK_STREAM, 0);
			socket_nopipe(sock);
			socket_send_timeout(sock, 60);//do not set recv timeout as client may rely on blocking read without using select()
			socket_tcp_keepalive(sock);
			if(nodelay){
				socket_tcp_nodelay(sock);
			}
			if(!block){
				socket_block(sock, 0);
			}
			int res;
			struct addrinfo* result;
			struct addrinfo hints={0};
			hints.ai_family=AF_INET;
			hints.ai_socktype=SOCK_STREAM;
			char hoststr[512];
			char portstr[32];
			hostname=lookup_host(hostname);
			const char* col=strchr(hostname, ':');
			if(col&&strlen(col+1)>0){//port is part of hostname
				size_t nn=col-hostname;
				if(nn+1>sizeof(hoststr)) nn=sizeof(hoststr)-1;
				strncpy(hoststr, hostname, nn); hoststr[nn]='\0';
				nn=strlen(col+1);
				if(nn+1>sizeof(portstr)) nn=sizeof(portstr)-1;
				strncpy(portstr, col+1, nn); portstr[nn]='\0';
			} else{
				snprintf(hoststr, sizeof(hoststr), "%s", hostname);
				snprintf(portstr, sizeof(portstr), "%d", port);
			}
			if((res=getaddrinfo(hoststr, portstr, &hints, &result))){
				dbg("getaddrinfo for %s failed with %d: %s\n", hostname, res, gai_strerror(res));
				return -1;
			}
			res=connect(sock, result->ai_addr, sizeof(servername));
			freeaddrinfo(result);
			/* Give the socket the target hostname. */
			if(res<0){
				//warning_once("connect to %s at %d failed: %s\n", hostname, port, strerror(errno));
				close(sock);
				if(!block){
					return -1;
				}
				sleep(10);
			} else{
				return sock;
			}
		}
	}
	return -1;
}
