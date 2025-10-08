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
	\file scheduler_ws.c This file contains routines to use libwebsockets
	handle http connection, serving files, upgrade the connection to websocket,
	and read and write websocket messages.
*/
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
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

#include "../sys/mem.h" //make sure memory debugging is used 
#include "../sys/misc.h"
#include "../sys/sockio.h"
#include "../sys/sock.h"
#include "../sys/dlist.h"
#include "scheduler.h"
static struct lws_context *context=0;
const int tx_buffer_size=65536;
static int lws_forward_tcp(int fd, int nlen, int mode, void *userdata);
/**
 * @brief Register with listen_socket to handle lws clients
 *
 * @param pollfd
 * @param flag
 * @return int
 */
static int
lws_respond(struct pollfd *pollfd, int flag){
	if(!context) return 0;
	//info("lws_respond with %d (revents=%d, read=%d, write=%d)\n", pollfd->fd, pollfd->revents, pollfd->revents&POLLIN, pollfd->revents&POLLOUT);
	if(flag==-1){
		dbg_time("flag==-1, call lws_context_destroy on %d\n", pollfd->fd);
		lws_context_destroy(context);//close context.
		context=0;
		return 0;
	} else{
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
callback_http(struct lws *wsi, enum lws_callback_reasons reason, void *user, void *in, size_t nin){
	switch(reason){
	case LWS_CALLBACK_ADD_POLL_FD:
	{
		struct lws_pollargs *p=(struct lws_pollargs *)in;
		listen_port_add(p->fd, p->events, lws_respond, "ws");
	}break;
	case LWS_CALLBACK_CHANGE_MODE_POLL_FD:
	{//change_mode maybe folloed by del_poll. 
		struct lws_pollargs *p=(struct lws_pollargs *)in;
		listen_port_add(p->fd, p->events, lws_respond, "ws");
	}break;
	case LWS_CALLBACK_DEL_POLL_FD:{
		struct lws_pollargs *p=(struct lws_pollargs *)in;
		listen_port_del(p->fd, 0, "callback_http");//do not close fd. lws owns it.
	}break;
	default://let the default handler handle the rest
		return lws_callback_http_dummy(wsi, reason, user, in, nin);
	}
	return 0;
}

typedef struct msg_t{
	dlist list[1];
	char *p;
	double ck;//timing
	int len;//valid payload (exclude LWS_PRE)
	int size;//total memory allocated
	int offset;//offset from the beginning for continuation
	int mode;//0: text, 1: binary
}msg_t;

struct protocol_t{//per client
	struct lws *wsi;
	dlist *phead;//active proxy data list head
	dlist *ptail;//active proxy data list tail
	int p_count;//messages pending.
	int p_drop;
};

static int
callback_maos(struct lws *wsi, enum lws_callback_reasons reason, void *user, void *in, size_t nin){
	struct protocol_t *pss=(struct protocol_t *)user;
	const char *pname=lws_get_protocol(wsi)->name;
	int pending=0;
	switch(reason){
	case LWS_CALLBACK_ESTABLISHED:
		lwsl_notice("LWS_CALLBACK_ESTABLISHED for %s\n", pname);
		pss->wsi=wsi;
		if(!strcmp(pname, "maos-monitor-protocol")){
			char cmd[11]="1&MONITOR;";
			ws_proxy_command(cmd, strlen(cmd), (ws_proxy_t){ .forward=lws_forward_tcp, .userdata=pss });
		}
		break;
	case LWS_CALLBACK_CLOSED://client closed
		lwsl_notice("LWS_CALLBACK_CLOSED for %s\n", pname);
		while(pss->phead){//remove buffer list
			msg_t *pdata=(msg_t*)remove_head(&pss->phead, &pss->ptail);
			free(pdata->p);
			free(pdata);
		}
		ws_proxy_remove(pss, 1);
		break;

	case LWS_CALLBACK_PROTOCOL_DESTROY://upon server exits
		lwsl_notice("LWS_CALLBACK_PROTOCOL_DESTROY for %s\n", pname);
		break;

	case LWS_CALLBACK_SERVER_WRITEABLE:
		//Notice that lws_write automatically retries when partial data is sent. No need to handle in user code.
#define CHECK_WRITE_SUCCESS(payload, len, FLAG) \
				if(lws_write(wsi, (unsigned char*)payload+LWS_PRE, len, (enum lws_write_protocol)FLAG)<0){\
					lwsl_err("ERROR: Failed to writing to client %p\n", wsi);\
					return -1;\
				}else

		do{
			pending=0;//continue writing if set
			msg_t* msg=(msg_t*)pss->phead;
			if(msg&&msg->len){//proxy data is pending
				int len=MIN(msg->len-msg->offset, tx_buffer_size);//how much to send. 64kB
				int protocol;
				if(msg->offset+len<msg->len){//partial send
					if(msg->offset){//continuation
						protocol=LWS_WRITE_CONTINUATION|LWS_WRITE_NO_FIN;
					} else{//first frame with continuation
						protocol=msg->mode|LWS_WRITE_NO_FIN;
						msg->ck=myclockd();
					}
				} else{//final segment or no fragmentation.
					if(msg->offset){//continuation
						protocol=LWS_WRITE_CONTINUATION; //final fragment
					} else{
						protocol=msg->mode; //no fragmentation
					}
				}
				CHECK_WRITE_SUCCESS(msg->p+msg->offset, len, protocol){
					if(msg->offset+len==msg->len){//finished
						if(msg->offset){
							//info_time("count=%d payload %d %.3f MB/s. (%d dropped)\n", pss->p_count, msg->len, msg->len/(msg->ck*1048576), pss->p_drop);
							pss->p_drop=0;
						}
						msg_t *pdata=(msg_t*)remove_head(&pss->phead, &pss->ptail);
						pdata->len=0;//set data as empty
						pdata->offset=0;//reset offset
						add_tail(&pss->phead, &pss->ptail, pdata->list);
						pss->p_count--;
					} else{
						msg->offset+=len;//wait for next run
					}
				}
				msg=(msg_t*)pss->phead;
				if(msg&&msg->len) pending=1;
			}
		} while(!lws_send_pipe_choked(wsi)&&pending);//lws_send_pipe_choked() check whether you can continue write data
		if(pending){//request further writing
			lws_callback_on_writable(wsi);
		}

		break;

	case LWS_CALLBACK_RECEIVE://receive client message
		ws_proxy_command((char *)in, nin, (ws_proxy_t){ .forward=lws_forward_tcp, .userdata=pss });
		break;

	default:
		break;
	}

	return 0;

}
static struct lws_protocols protocols[]={ /* first protocol must always be HTTP handler */
	{ "http", callback_http, 0, 0, 0, NULL},
	{"maos-monitor-protocol",callback_maos,sizeof(struct protocol_t),1024,0, 0},
	{"maos-drawdaemon-protocol",callback_maos,sizeof(struct protocol_t), 1024, 0, 0},
	{ NULL, NULL, 0, 0, 0, 0 } /* terminator */
};
static const struct lws_http_mount mount={
	.mountpoint="/",				/* mountpoint URL */
	.origin=SRCDIR "/tools",	/* serve from dir */
	.def="monitor.html",	/* default filename */
	.origin_protocol=LWSMPRO_FILE,		/* serve with files in a dir */
	.mountpoint_len=1,				/* char count */
};

/**
 * @brief Called by scheduler to start lws web server
 *
 * @param port
 * @return int
 */
int start_lws(short port){
	if(context) return 0;
	struct lws_context_creation_info info;
	memset(&info, 0, sizeof info);
	info.port=port;
	info.protocols=protocols;
	/*this is critical. otherwise UTF-8 error in client.*/
	info.gid=-1;
	info.uid=-1;
	info.mounts=&mount;
	context=lws_create_context(&info);
	lws_set_log_level(LLL_ERR|LLL_WARN, NULL);//default is err, warn, notice.
	if(!context){
		lwsl_err("libwebsocket init failed\n");
		return -1;
	}
	return 0;
}
/**
 * @brief Read TCP data and forward to websocket clients
 *
 * @param fd 	File descriptor to read from
 * @param nlen 	Length of data to read
 * @param mode 	text (0) or binary (1) mode
 * @param userdata 	lws context
 * @return int 	return -1 to indicate error
 */
static int
lws_forward_tcp(int fd, int nlen, int mode, void *userdata){
	struct protocol_t *pss=(struct protocol_t *)userdata;
	struct lws *wsi=pss->wsi;

	msg_t *pdata=(msg_t*)pss->ptail;
	if(pdata&&(!pdata->len||pss->p_count>100)){//empty slot or too much being buffered.
		pdata=(msg_t*)remove_tail(&pss->phead, &pss->ptail);
	} else{//allocate new data
		pdata=mycalloc(1, msg_t);
	}
	if(!pdata->len){//not replacing old data
		pss->p_count++;
	} else{
		//dbg_time("sock=%d count=%d payload %d is dropped\n", sock, pss->p_count, pdata->len);
		pss->p_drop++;
	}


	int nlen2=nlen+LWS_PRE;
	if(pdata->size<nlen2){
		pdata->size=nlen2;
		pdata->p=(char*)realloc(pdata->p, nlen2);
	}
	pdata->len=nlen;
	pdata->mode=mode;
	int *p2=(int *)(pdata->p+LWS_PRE);
	if(stread(fd, p2, nlen)){
		warning_time("read tcp failed\n");
		return -1;
	}

	msg_t *msg=(msg_t*)pss->phead;
	if(msg&&msg->offset!=0){//do not disturb a partially send frame
		msg=(msg_t*)remove_head(&pss->phead, &pss->ptail);
	}
	add_head(&pss->phead, &pss->ptail, pdata->list);
	if(msg){
		add_head(&pss->phead, &pss->ptail, msg->list); msg=NULL;
	}
	lws_callback_on_writable(wsi);
	return 0;
}
