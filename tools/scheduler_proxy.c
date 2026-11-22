/*
  Copyright 2009-2026 Lianqi Wang
  
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
   Add a new proxy link, or replace the existing one with the same fd.

   @param[in] ws Information for the helper routine. ws.fd is the socket to read from.
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
   Remove the proxy links of a client, removing their sockets from listen_port
   and closing the sockets passed to maos.

   @param[in] userdata Parameter to identify the client.
   @param[in] toclose Close the socket if set.
   @param[in] keep_monitor Keep the link of a direct monitor (whose socket is the client itself).
*/
static void ws_proxy_remove_do(void *userdata, int toclose, int keep_monitor){
	int found=0;
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata && !(keep_monitor && ws_proxy[iws].ismonitor)){
			found=1;
			if(ws_proxy[iws].ismonitor){
				monitor_remove((int)(long)(ws_proxy[iws].userdata));
			}
			if(ws_proxy[iws].fd>=0){//-1 if the link is already closed by ws_proxy_read
				listen_port_del(ws_proxy[iws].fd, toclose, "ws_proxy_remove");
			}
			//info_time("ws_proxy_remove: fd=%d for %p\n", ws_proxy[iws].fd, userdata);
			if(ws_proxy[iws].fd_remote>0)close(ws_proxy[iws].fd_remote);//notify maos that the proxy is closed
			
			ws_proxy[iws].fd_remote=-1;
			ws_proxy[iws].fd=-1;
			ws_proxy[iws].ismonitor=0;
			ws_proxy[iws].userdata=NULL;
		}
	}
	if(!found){
		dbg_time("%p not found.\n", userdata);//ok, may not be registered yet.
	}
}
/**
   Remove all proxy links of a client, removing their sockets from listen_port
   and closing the sockets passed to maos.

   @param[in] userdata Parameter to identify the client.
   @param[in] toclose Close the socket if set.
*/
void ws_proxy_remove(void *userdata, int toclose){
	ws_proxy_remove_do(userdata, toclose, 0);
}
/**
   Obtain the server socket number for the proxy link of a client. The link of
   a direct monitor is skipped as its socket is the client itself.

   @param[in] userdata Parameter to identify the client.
   @return fd for the socket, or -1 if not found.
*/
static int ws_proxy_get_fd(void *userdata){
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata && !ws_proxy[iws].ismonitor && ws_proxy[iws].fd>=0){
			return ws_proxy[iws].fd;
		}
	}
	return -1;
}
/**
   Read from TCP server (scheduler or draw) and forward to the websocket client.

   @param[in] pfd Poll entry of the TCP socket.
   @param[in] flag May be set to -1 to ask the client to close the connection.
   @return -1 to indicate error (the link is closed), 0 otherwise.
*/
static int ws_proxy_read(struct pollfd *pfd, int flag){
	if(!pfd||pfd->fd==-1) return 0;
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
	int ans=0;
	if(flag==-1){//ask client to shutdown
		stwriteint(sock, -1);
	}else if(!(pfd->revents&POLLIN)){//POLLERR, POLLNVAL, or POLLHUP without data. poll() would return immediately again.
		warning_time("revents=%d is received. fd=%d\n", pfd->revents, sock);
		ans=-1;
	}else{
		int cmd[3]={0};
		int nlen2=(int)recv(sock, cmd, sizeof(cmd), MSG_PEEK);//peek at the command. do not consume it as the proxy needs it
		if(nlen2>=(int)sizeof(int) && cmd[0]==DRAW_ENTRY && nlen2<(int)sizeof(cmd)){//wait for the full header
			nlen2=(int)recv(sock, cmd, sizeof(cmd), MSG_PEEK|MSG_WAITALL);
		}
		if(nlen2<(int)sizeof(int)){
			warning_time("failed to peek; nlen2=%d, revents=%d\n", nlen2, pfd->revents);
			ans=-1;//break the link
		}else if(cmd[0]==DRAW_HEARTBEAT){//ignored
			if(stread(sock, cmd, sizeof(int))){
				warning_time("Unable to read 1 int, close connection.\n");
				ans=-1;
			}
		}else if(cmd[0]==DRAW_ENTRY && nlen2==(int)sizeof(cmd)){
			//cmd[1] is the payload length. cmd[2] is -1 for text from the scheduler, or the draw command for binary data.
			int nlen=cmd[1];
			int mode=cmd[2]==-1?0:1;
			if(nlen<0 || nlen>INT_MAX-(int)sizeof(cmd)){
				warning_time("invalid length %d, close connection.\n", nlen);
				ans=-1;
			}else if(mode==0){//text data, the client does not need the header. read and drop the header
				if(stread(sock, cmd, sizeof(cmd))){
					warning_time("Unable to read 3 int, close connection.\n");
					ans=-1;
				}
			} else{//binary data. the client wants the header as well.
				nlen+=(int)sizeof(cmd);
			}
			if(!ans && ws_proxy[iws].forward(sock, nlen, mode, ws_proxy[iws].userdata)){
				ans=-1;
			}
		} else{
			warning("cmd=%d is not prefixed by DRAW_ENTRY. close connection.\n", cmd[0]);
			ans=-1;
		}
	}
	if(ans){//listen_port closes the socket.
		if(ws_proxy[iws].fd_remote>0) close(ws_proxy[iws].fd_remote);
		ws_proxy[iws].fd_remote=-1;
		ws_proxy[iws].fd=-1;
	}
	return ans;
}

/**
   Handle requests (in text format) from websockets client. The format is
   `PID&COMMAND_NAME;` with field separated by & and terminated by ; or \n.
   Supported commands are REMOVE, KILL, RESTART, MONITOR, DRAW, DRAW_FIGFN,
   DRAW_PAUSE, DRAW_RESUME, and DRAW_SINGLE.

   @param[in,out] in The command buffer. Modified in place to separate fields.
   @param[in] len The length of the command buffer.
   @param[in] ws Information for the helper function to write to websocket client.
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
		int pid=(int)strtol(in, 0, 10);
		*sep='&';
		sep++;

		char *sep2=strchr(sep, '&');
		if(sep2){
			*sep2=0; sep2++;//next entry
		}

		//LOCK(mutex_sch);
		if(!strcmp(sep, "REMOVE")){
			runned_remove(pid);
		} else if(!strcmp(sep, "KILL")){
			info_time("HTTP client send term signal to %5d term signal.\n", pid);
			running_kill(pid);
		} else if(!strcmp(sep, "RESTART")){
			dbg_time("HTTP client request restart job %d\n", pid);
			runned_restart(pid);
		} else if(pid>0 && !strcmp(sep, "MONITOR")){
			if(ws.forward){//forward mode
				if(ws_proxy_get_fd(ws.userdata)>-1){
					dbg_time("Reopen monitor for %p.\n", ws.userdata);
					ws_proxy_remove_do(ws.userdata, 1, 1);
				}
				int sock=scheduler_connect_self(1);
				int cmd[2]={CMD_MONITOR, INT_MIN|1<<2|1<<1};
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
				monitor_add((int)(long)(ws.userdata), INT_MIN|1<<2|1<<1, ws.send, ws.userdata);
				ws.fd=(int)(long)(ws.userdata);
				ws.ismonitor=1;
				ws_proxy_add(ws);//keep reference so that we can disconnect it when browser disconnects
			}
		} else if(!strcmp(sep, "DRAW")){
			int sock=ws_proxy_get_fd(ws.userdata);
			if(sock>-1){
				dbg_time("Reopen drawing for %p.\n", ws.userdata);
				ws_proxy_remove_do(ws.userdata, 1, 1);
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
					ws.fd_remote=-1;//sv2[1] is closed below once passed. Do not close it again in ws_proxy_remove.
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
					sep2=strchr(sep3, '&');
					if(sep2) *sep2=0;//terminate string
				}
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
