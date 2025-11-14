/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
	\file scheduler_proxy.c This file contains proxy code that linkes TCP server
	with websocket client.

	ws_proxy_command: handles command comming from websocket clients. If the
	command is for receiving job list or drawing data, a proxy is setup by
	ws_proxy_add with a helper routine. The data is received by ws_proxy_read
	which calls the helper routine to send the data.

	The built-in websocket server can avoid a proxy and let scheduler directly
	write to the websocket. This is also mediated in this file.
*/
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <sys/socket.h>
#include <poll.h>
#include "../sys/sys.h"
#include "scheduler.h"

static ws_proxy_t *ws_proxy=NULL;
static int nws_proxy=0;
/**
 * @brief Add a new proxy link. fd is websocket fd
 *
 * @param ws information for the helper routine
 */
static void ws_proxy_add(ws_proxy_t ws){
	int jws=-1;
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].fd==ws.fd){
			jws=iws;//replace existing entry
			break;
		}
		if(jws==-1 && ws_proxy[iws].fd==-1){
			jws=iws;//save insert location
		}
	}
	if(jws==-1){
		jws=nws_proxy;
		nws_proxy++;
		ws_proxy_t *tmp=myrealloc(ws_proxy, nws_proxy, ws_proxy_t);
		if(tmp){
			ws_proxy=tmp;
		} else{
			nws_proxy--;
			warning("Unable to allocate memory for ws_proxy\n");
			return;
		}
	}
	ws_proxy[jws]=ws;
	//info_time("ws_proxy_add: fd=%d for %p\n", sock_tcp, userdata);
}
/**
 * @brief Remove an existing proxy link
 *
 * @param userdata 	parameter to identify the proxy
 * @param toclose close the socket if set
 */
void ws_proxy_remove(void *userdata, int toclose){
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata){
			if(ws_proxy[iws].ismonitor){
				monitor_remove((int)(long)(ws_proxy[iws].userdata));
			}
			listen_port_del(ws_proxy[iws].fd, toclose, "ws_proxy_remove");
			//info_time("ws_proxy_remove: fd=%d for %p\n", ws_proxy[iws].fd, userdata);
			if(ws_proxy[iws].fd_remote>0)close(ws_proxy[iws].fd_remote);//notify maos that the proxy is closed
			
			ws_proxy[iws].fd_remote=-1;
			ws_proxy[iws].fd=-1;
			ws_proxy[iws].userdata=NULL;
			return;
		}
	}
	dbg_time("%p not found.\n", userdata);//ok, may not be registered yet.
}
/**
 * @brief Obtain the server socket number for the proxy link
 *
 * @param userdata 	parameter to identify the proxy
 * @return int  	fd for the socket
 */
static int ws_proxy_get_fd(void *userdata){
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata){
			return ws_proxy[iws].fd;
		}
	}
	return -1;
}
/**
 * @brief Read from TCP server and forward to the websocket client
 *
 * @param pfd	File descriptor
 * @param flag	Maybe set to -1 to ask the client to close the connection
 * @return int	Return -1 to indicate error, 0 Otherwise.
 */
static int ws_proxy_read(struct pollfd *pfd, int flag){
	if(!pfd||pfd->fd==-1||(pfd->revents&POLLIN)==0) return 0;
	if(pfd->revents&POLLHUP){
		warning_time("POLLHUP is received. fd=%d\n", pfd->fd);
		return -1;
	}
	if(flag==-1){//ask client to shutdown
		stwriteint(pfd->fd, -1);
	}
	int iws;
	for(iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].fd==pfd->fd){
			break;
		}
	}
	if(iws==nws_proxy){
		warning("unable to find ws_proxy entry.\n");
		return -1;
	}

	int sock=ws_proxy[iws].fd;

	int cmd[3]={0};
	int nlen2=recv(sock, &cmd, sizeof(cmd), MSG_PEEK);//peek at the command
	if(nlen2<(int)sizeof(int)){
		warning_time("failed to peek; nlen2=%d, revents=%d\n", nlen2, pfd->revents);
		return 0;
	}
	if(cmd[0]==DRAW_HEARTBEAT||cmd[2]==DRAW_HEARTBEAT){
		int nlen3;
		if((nlen3=recv(sock, &cmd, nlen2, MSG_DONTWAIT))<nlen2){
			warning_time("Unable to read %d ints, got %d\n", nlen2, nlen3);
			if(nlen3<0 && (errno==EWOULDBLOCK || errno==EAGAIN)){
				return 0;
			}else{
				return -1;
			}
		}
		return 0;//ignored
	}
	if(cmd[0]==DRAW_ENTRY){
		int nlen=cmd[1];
		int mode=cmd[2]==-1?0:1;
		if(cmd[2]==0){//text data, no need for the header. read and drop the header
			if(stread(sock, &cmd, sizeof(cmd))){
				warning_time("Unable to read 3 ints\n");
				return -1;
			}
		} else{
			nlen+=3*sizeof(int);//additional data
		}
		if(ws_proxy[iws].forward(sock, nlen, mode, ws_proxy[iws].userdata)){
			ws_proxy[iws].fd=-1;
			return -1;
		}
	} else{
		warning("cmd=%d is not prefixed by DRAW_ENTRY. close connection.\n", cmd[0]);
		return -1;
	}
	return 0;
}

/**
 * @brief handle requests (in text format) from websockets client. The format is
  `PID&COMMAND_NAME;` with field separated by & and terminated by ; or \n
 *
 * @param in    The command buffer
 * @param len 	The length of the command buffer
 * @param ws 	Information for the helper function to write to websocket client
 */
void ws_proxy_command(char *in, size_t len, ws_proxy_t ws){
	if(in[len]){
		warning_time("buf string is not terminated\n");
		in[len]=0;//make sure string is null terminated
	}
	char *end=NULL;//find and mark end of current command
	while(len>1){//len includes final 0
		if(end){//handle command separation
			if(end<in+len){//locate next command
				len-=(end-in);
				in=end;
			} else{
				len=0;
				break;
			}
		}
		char *tmp1=strchr(in, '\n');
		char *tmp2=strchr(in, ';');
		if(tmp1&&(!tmp2||tmp1<=tmp2)){
			*tmp1=0;
			end=tmp1+1;
		} else if(tmp2&&(!tmp1||tmp2<=tmp1)){
			*tmp2=0;
			end=tmp2+1;
		} else{
			end=in+len+1;
		}
		char *sep=strchr(in, '&');
		if(!sep){
			warning_time("Unable to handle cmd: {%.20s}, len={%zu}\n", in, len);
			continue;
		}
		*sep=0;
		int pid=strtol(in, 0, 10);
		*sep='&';
		sep++;

		char *sep2=strchr(sep, '&');
		if(sep2){
			*sep2=0; sep2++;//next entry
		}

		//LOCK(mutex_sch);
		if(pid>0 && !strcmp(sep, "REMOVE")){
			runned_remove(pid);
		} else if(pid>0 && !strcmp(sep, "KILL")){
			info_time("HTTP client send term signal to %5d term signal.\n", pid);
			running_kill(pid);
		} else if(pid>0 && !strcmp(sep, "MONITOR")){
			if(ws.forward){//forward mode
				if(ws_proxy_get_fd(ws.userdata)>-1){
					dbg_time("Reopen monitor for %p.\n", ws.userdata);
					ws_proxy_remove(ws.userdata, 1);
				}
				int sock=scheduler_connect_self(1);
				int cmd[2]={CMD_MONITOR, 1<<31|1<<2|1<<1};
				if(sock==-1 || stwriteintarr(sock, cmd, 2)){
					warning("Unable to talk to scheduler");
					if(sock!=-1) close(sock);
				}else{
					ws.fd=sock;
					ws.fd_remote=-1;
					ws_proxy_add(ws);//add proxy first to have higher priority in poll()
					listen_port_add(sock, POLLIN, ws_proxy_read, "ws_proxy_read");
				}
			} else if(ws.send){
				monitor_add((int)(long)(ws.userdata), 1<<31|1<<2|1<<1, ws.send, ws.userdata);
				ws.fd=(int)(long)(ws.userdata);
				ws.ismonitor=1;
				ws_proxy_add(ws);//keep reference so that we can disconnect it when browser disconnects
			}
		} else if(!strcmp(sep, "DRAW")){
			int sock=ws_proxy_get_fd(ws.userdata);
			if(sock>-1){
				dbg_time("Reopen drawing for %p.\n", ws.userdata);
				ws_proxy_remove(ws.userdata, 1);
			}
			
			int sv2[2];//socket pair for communication.
			/*one end of sv2 will be passed back to call, the other end of sv2 will be passed to drawdaemon.*/
			if(!socketpair(AF_UNIX, SOCK_STREAM, 0, sv2)){
				int status=send_draw_sock(sv2[1], pid);
				if(!status){
					streadint(sv2[0], &status);
				}
				if(!status){
					ws.fd=sv2[0];
					ws.fd_remote=sv2[1];
					ws_proxy_add(ws);//add proxy first to have higher priority in poll()
					listen_port_add(sv2[0], POLLIN, ws_proxy_read, "ws_proxy_read");
				} else{
					close(sv2[0]);
					warning_time("Failed to request DRAW to maos pid=%d\n", pid);
				}
				close(sv2[1]);
			} else{
				perror("socketpair");
			}
		
		} else if(!strcmp(sep, "DRAW_FIGFN")){
			char *fig=0, *fn=0;
			if(sep2){
				char *sep3=strchr(sep2, '&');
				fig=sep2;
				if(sep3){
					*sep3=0; sep3++;
					fn=sep3;
				}
				sep2=strchr(sep3, '&');
				if(sep2) *sep2=0;//terminate string
			}
			if(fig&&fn){
				int sock=ws_proxy_get_fd(ws.userdata);
				if(sock>-1){
					if(stwriteint(sock, DRAW_FIGFN)||(fig&&stwritestr(sock, fig))||(fn&&stwritestr(sock, fn))){
						info_time("Send %s {%s} {%s} to draw failed.\n", sep, fig?fig:"", fn?fn:"");
					}
				}
			}
		} else if(!strcmp(sep, "DRAW_PAUSE")){
			int sock=ws_proxy_get_fd(ws.userdata);
			if(sock>-1){
				if(stwriteint(sock, DRAW_PAUSE)){
					info_time("Send %s to draw failed.\n", sep);
				}
			}
		} else if(!strcmp(sep, "DRAW_RESUME")){
			int sock=ws_proxy_get_fd(ws.userdata);
			if(sock>-1){
				if(stwriteint(sock, DRAW_RESUME)){
					info_time("Send %s to draw failed.\n", sep);
				}
			}
		} else if(!strcmp(sep, "DRAW_SINGLE")){
			int sock=ws_proxy_get_fd(ws.userdata);
			if(sock>-1){
				if(stwriteint(sock, DRAW_SINGLE)){
					info_time("Send %s to draw failed.\n", sep);
				}
			}
		} else{
			warning_time("Unable to handle cmd: {%s}, {%.20s}, len={%zu}\n", sep, in, len);
		}
	}
	//UNLOCK(mutex_sch);
}
