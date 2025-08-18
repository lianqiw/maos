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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <syslog.h>
#include <sys/time.h>
#include <libwebsockets.h>
#include <poll.h>
//#include "../sys/sys.h"
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include "scheduler_ws.h"
#include "../sys/scheduler.h"
#include "../sys/mem.h" //make sure memory debugging is used 
#include "../sys/misc.h"
#include "../sys/sockio.h"
#include "../sys/sock.h"
static struct lws_context* context=0;
//Calls by scheduler timeout periodically or upon job status change to service lws requests.
int ws_service_fd(struct pollfd *pollfd, int flag){
	if(!context) return 0;
	//info("ws_service_fd with %d (revents=%d, read=%d, write=%d)\n", pollfd->fd, pollfd->revents, pollfd->revents&POLLIN, pollfd->revents&POLLOUT);
	if(flag==-1){
		lws_context_destroy(context);//close context.
		context=0;
		return 0;
	}else{
		return lws_service_fd(context, pollfd);
	}
}
/**
   Contains code for handling websocket connection using libwebsockets. Part of
   the code is copied from test-server.c in libwebsockets

   Upon conection, existing job list is converted to payload list and send to the new client.
   Upon job status change, job info is converted to payload ringbuffer and send to all the clients.
*/
static int
callback_http(struct lws* wsi, enum lws_callback_reasons reason, void* user, void* in, size_t nin){
	switch(reason){
		case LWS_CALLBACK_ADD_POLL_FD:
		{
			struct lws_pollargs *p = (struct lws_pollargs *)in;
			listen_port_add(p->fd, p->events, ws_service_fd);
		}break;
		case LWS_CALLBACK_CHANGE_MODE_POLL_FD:
		{
			struct lws_pollargs *p = (struct lws_pollargs *)in;
			listen_port_add(p->fd, p->events, ws_service_fd);
		}break;
		case LWS_CALLBACK_DEL_POLL_FD:{
			struct lws_pollargs *p = (struct lws_pollargs *)in;
			listen_port_remove(p->fd, 0);
		}break;
		default://let the default handler handle the rest
		return lws_callback_http_dummy(wsi, reason, user, in, nin);
	}
	return 0;
}
struct a_message{
	char* payload;
	size_t len;//valid payload (exclude LWS_PRE)
	size_t size;//total memory allocated
};
typedef struct p_message{
	char *p;
	int len;
	int size;
	int pending;
}p_message;
void ws_proxy_remove(void *userdata, int closed);
#define MAX_MESSAGE_QUEUE 1024
/*The receiver modifies the head of the ring buffer. Each client maintains its
  own tail of the ring buffer so they can proceed that their own speed. */
static struct a_message ringbuffer[MAX_MESSAGE_QUEUE]={0};
static int ringbuffer_head=0;//ring buffer head for all clients
struct per_session_data__maos_monitor{//per client
	struct lws* wsi;
	l_message* lhead;//head of job list for initial connection
	p_message* pdata;
	int ringbuffer_tail;//ring buffer tail per client for job update
};
int nclient=0;//count number of clients
static int
callback_maos_monitor(struct lws* wsi, enum lws_callback_reasons reason, void* user, void* in, size_t nin){
	struct per_session_data__maos_monitor* pss=(struct per_session_data__maos_monitor*)user;
	int pending=0;
	switch(reason){
	case LWS_CALLBACK_ESTABLISHED:
		nclient++;
		lwsl_notice("LWS_CALLBACK_ESTABLISHED, nclient=%d\n", nclient);
		pss->ringbuffer_tail=ringbuffer_head;//initialize to zero length
		pss->wsi=wsi;
		pss->lhead=0;
		scheduler_push_ws(&pss->lhead, LWS_PRE);//send existing job information upon connection
		lws_callback_on_writable(wsi);
		break;
	case LWS_CALLBACK_CLOSED://client closed
		nclient--;
		lwsl_notice("LWS_CALLBACK_CLOSED, nclient=%d\n", nclient);
		ws_proxy_remove(pss, 0);
		break;
		
	case LWS_CALLBACK_PROTOCOL_DESTROY://upon server exits
		lwsl_notice("LWS_CALLBACK_PROTOCOL_DESTROY, nclient=%d\n", nclient);
		for(int n=0; n<(int)(sizeof ringbuffer/sizeof ringbuffer[0]); n++)
			if(ringbuffer[n].payload)
				free(ringbuffer[n].payload);
		break;

	case LWS_CALLBACK_SERVER_WRITEABLE:	
		//Notice that lws_write automatically retries when partial data is sent. No need to handle in user code.
	#define CHECK_WRITE_SUCCESS(payload, len, FLAG) \
				if(lws_write(wsi, (unsigned char*)payload+LWS_PRE, len, FLAG)<0){\
					lwsl_err("ERROR: Failed to writing to client %p\n", wsi);\
					return -1;\
				}else

		do{
			pending=0;//continue writing if set
			if(pss->pdata && pss->pdata->pending){//proxy data is pending
				CHECK_WRITE_SUCCESS(pss->pdata->p, pss->pdata->len, LWS_WRITE_BINARY)
				{
					pss->pdata->pending=0;//mark sent
				}
			}else if(pss->lhead){/*initialization with list*/
				CHECK_WRITE_SUCCESS(pss->lhead->payload, pss->lhead->len, LWS_WRITE_TEXT)
				{
					l_message* tmp=pss->lhead;
					pss->lhead=pss->lhead->next;
					free(tmp->payload);
					free(tmp);
				}
				if(pss->lhead) pending=1;
			}else if(pss->ringbuffer_tail!=ringbuffer_head){//send using ring buffer tail
				CHECK_WRITE_SUCCESS(ringbuffer[pss->ringbuffer_tail].payload, ringbuffer[pss->ringbuffer_tail].len, LWS_WRITE_TEXT)
				{
					if(pss->ringbuffer_tail==(MAX_MESSAGE_QUEUE-1)){
						pss->ringbuffer_tail=0;
					}else{
						pss->ringbuffer_tail++;
					}
				}
				if(pss->ringbuffer_tail!=ringbuffer_head){
					pending=2;//more data to write
				}
			}
		}while(!lws_send_pipe_choked(wsi) && pending);//lws_send_pipe_choked() check whether you can continue write data
		if(pending){//request further writing
			lws_callback_on_writable(wsi);
		}
	
		break;

	case LWS_CALLBACK_RECEIVE://receive client message
		scheduler_receive_ws((char*)in, nin, pss);
		break;

		/*
		 * this just demonstrates how to use the protocol filter. If you won't
		 * study and reject connections based on header content, you don't need
		 * to handle this callback
		 */

	case LWS_CALLBACK_FILTER_PROTOCOL_CONNECTION:
		/* you could return non-zero here and kill the connection */
		break;

	default:
		break;
	}

	return 0;
}

static struct lws_protocols protocols[]={
	/* first protocol must always be HTTP handler */
	{ "http", callback_http, 0, 0, 0, NULL},
	
	{
	"maos-monitor-protocol",
	callback_maos_monitor,
	sizeof(struct per_session_data__maos_monitor),
	1024,                   /* max frame size / rx buffer */
	0,                      //id unused
	0,                      //user unused
	},
	{ NULL, NULL, 0, 0, 0, 0 } /* terminator */
};
static const struct lws_http_mount mount = {
	.mountpoint		= "/",				/* mountpoint URL */
	.origin			= SRCDIR "/tools",	/* serve from dir */
	.def			= "monitor.html",	/* default filename */
	.origin_protocol= LWSMPRO_FILE,		/* serve with files in a dir */
	.mountpoint_len	= 1,				/* char count */
};


int ws_start(short port){
	if(context) return 0;
	struct lws_context_creation_info info;
	memset(&info, 0, sizeof info);
	info.port=port;
	info.protocols=protocols;
	/*this is critical. otherwise UTF-8 error in client.*/
	info.gid=-1;
	info.uid=-1;
	info.mounts = &mount;
	context=lws_create_context(&info);
	lws_set_log_level(LLL_ERR|LLL_WARN, NULL);//default is err, warn, notice.
	if(!context){
		lwsl_err("libwebsocket init failed\n");
		return -1;
	}
	return 0;
}
void ws_end(){
	if(context){
		lws_context_destroy(context);
		context=0;
		lwsl_notice("libwebsockets-test-server exited cleanly\n");
	}
}

//Initiate server to client message with ring buffer.
void ws_append(const char* in, size_t len){
	if(!context || !nclient) return;
	size_t new_size=LWS_PRE+len;
	//make sure buffer is big enough.
	if(ringbuffer[ringbuffer_head].size<new_size){
		ringbuffer[ringbuffer_head].payload=realloc(
			ringbuffer[ringbuffer_head].payload, new_size);
		ringbuffer[ringbuffer_head].size=new_size;
	}
	ringbuffer[ringbuffer_head].len=len;
	memcpy(ringbuffer[ringbuffer_head].payload+LWS_PRE, in, len);
	if(ringbuffer_head+1==MAX_MESSAGE_QUEUE){
		ringbuffer_head=0;/*wrap over*/
	} else{
		ringbuffer_head++;
	}
	lws_callback_on_writable_all_protocol(context, protocols+1);
}

#define CATCH_TO(A,p) if(A) {if(errno==EAGAIN ||errno==EWOULDBLOCK){return 0;}\
else{ ws_proxy[iws].fd=-1; dbg_time("read " #p " failed %s.", strerror(errno)); return -1;}}
typedef struct ws_proxy_t{
	int fd;
	void *userdata;
}ws_proxy_t;
ws_proxy_t* ws_proxy=NULL;
int nws_proxy=0;
void ws_proxy_add(int sock_tcp, void *userdata){
	int iws;
	for(iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].fd==sock_tcp || ws_proxy[iws].fd==-1){
			break;
		}
	}
	if(iws==nws_proxy){
		nws_proxy++;
		ws_proxy=realloc(ws_proxy, nws_proxy*sizeof(ws_proxy_t));
	}
	ws_proxy[iws].fd=sock_tcp;
	ws_proxy[iws].userdata=userdata;
}
void ws_proxy_remove(void *userdata, int closed){
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata){
			if(closed){
				listen_port_remove(ws_proxy[iws].fd, closed);
			}
			ws_proxy[iws].fd=-1;
			ws_proxy[iws].userdata=NULL;
			break;
		}
	}
}
void ws_proxy_write(int cmd, char *fig, char *fn, void *userdata){
	for(int iws=0; iws<nws_proxy; iws++){
		if(ws_proxy[iws].userdata==userdata){
			int sock=ws_proxy[iws].fd;
			if(stwriteint(sock, cmd)||
				(fig && stwritestr(sock, fig))||
				(fn && stwritestr(sock, fn))){
				info_time("Send %d {%s} {%s} to draw failed.\n", cmd, fig?fig:"", fn?fn:"");
			}else{
				//info_time("Send %d {%s} {%s} to draw ok.\n", cmd, fig?fig:"", fn?fn:"");
			}
		}
	}
}
int ws_proxy_read(struct pollfd *pfd, int flag){
	if(!pfd || pfd->fd==-1||(pfd->revents&POLLIN)==0) return 0;
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
	struct per_session_data__maos_monitor* pss=(struct per_session_data__maos_monitor*)ws_proxy[iws].userdata;

	if(!(pss->pdata)){
		pss->pdata=mycalloc(1, p_message);
	}
	p_message *pdata=pss->pdata;
	struct lws*wsi=pss->wsi;
	if(pdata->pending){
		//info_time("buffer is pending\n");
		return 0;
	}
	int sock=ws_proxy[iws].fd;
	int cmd=-1;
	int nlen=0;
	CATCH_TO(streadint(sock, &cmd), CMD);
	if(cmd==DRAW_ENTRY){
		CATCH_TO(streadint(sock, &nlen), nlen);
		CATCH_TO(streadint(sock, &cmd), CMD2);
		if(cmd==DRAW_HEARTBEAT) return 0;
		int nlen2=nlen+3*sizeof(int)+LWS_PRE;
		if(pdata->size<nlen2){
			pdata->size=nlen2;
			pdata->p=realloc(pdata->p, nlen2);
		}
		int *p2=(int*)(pdata->p+LWS_PRE);
		p2[0]=DRAW_ENTRY;
		p2[1]=nlen;
		p2[2]=cmd;
		if(nlen>0){
			CATCH_TO(stread(sock, &p2[3], nlen), p);
		}
		pdata->len=nlen+3*sizeof(int);
		//info("read %d with payload %d\n", cmd, pdata->len);
		pdata->pending=1;
		lws_callback_on_writable(wsi);
		//TODO: append to bhead
	}else if(cmd!=DRAW_HEARTBEAT){
		warning("cmd=%d is not prefixed by DRAW_ENTRY\n", cmd);
		//close(sock);
		//ws_fds[iws].fd=-1;
	}
	return 0;
}
