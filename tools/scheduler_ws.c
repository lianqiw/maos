/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
//#include "../sys/sys.h"
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include "scheduler_ws.h"
#include "../sys/mem.h" //make sure memory debugging is used 
#include "../sys/misc.h"
/**
   Contains code for handling websocket connection using libwebsockets. Part of
   the code is copied from test-server.c in libwebsockets
*/

#define LOCAL_RESOURCE_PATH SRCDIR"/tools"
const char* resource_path=LOCAL_RESOURCE_PATH;
/*
 * We take a strict whitelist approach to stop ../ attacks
 */
struct serveable{
	const char* urlpath;
	const char* mimetype;
};

struct per_session_data__http{
	int fd;
};

/*
 * this is just an example of parsing handshake headers, you don't need this
 * in your code unless you will filter allowing connections by the header
 * content
 */
static void
dump_handshake_dbg(struct lws* wsi){
	static const char* token_names[]={
		/*[WSI_TOKEN_GET_URI]		=*/ "GET URI",
		/*[WSI_TOKEN_POST_URI]		=*/ "POST URI",
		/*[WSI_TOKEN_HOST]		=*/ "Host",
		/*[WSI_TOKEN_CONNECTION]	=*/ "Connection",
		/*[WSI_TOKEN_KEY1]		=*/ "key 1",
		/*[WSI_TOKEN_KEY2]		=*/ "key 2",
		/*[WSI_TOKEN_PROTOCOL]		=*/ "Protocol",
		/*[WSI_TOKEN_UPGRADE]		=*/ "Upgrade",
		/*[WSI_TOKEN_ORIGIN]		=*/ "Origin",
		/*[WSI_TOKEN_DRAFT]		=*/ "Draft",
		/*[WSI_TOKEN_CHALLENGE]		=*/ "Challenge",

		/* new for 04 */
		/*[WSI_TOKEN_KEY]		=*/ "Key",
		/*[WSI_TOKEN_VERSION]		=*/ "Version",
		/*[WSI_TOKEN_SWORIGIN]		=*/ "Sworigin",

		/* new for 05 */
		/*[WSI_TOKEN_EXTENSIONS]	=*/ "Extensions",

		/* client receives these */
		/*[WSI_TOKEN_ACCEPT]		=*/ "Accept",
		/*[WSI_TOKEN_NONCE]		=*/ "Nonce",
		/*[WSI_TOKEN_HTTP]		=*/ "Http",

		"Accept:",
		"If-Modified-Since:",
		"Accept-Encoding:",
		"Accept-Language:",
		"Pragma:",
		"Cache-Control:",
		"Authorization:",
		"Cookie:",
		"Content-Length:",
		"Content-Type:",
		"Date:",
		"Range:",
		"Referer:",
		"Uri-Args:",

		/*[WSI_TOKEN_MUXURL]	=*/ "MuxURL",
	};
	char buf[256];

	for(size_t n=0; n<sizeof(token_names)/sizeof(token_names[0]); n++){
		if(!lws_hdr_total_length(wsi, (enum lws_token_indexes)n))
			continue;

		lws_hdr_copy(wsi, buf, sizeof buf, (enum lws_token_indexes)n);

		fprintf(stderr, "    %s = %s\n", token_names[n], buf);
	}
}

const char* get_mimetype(const char* file){
	int n=strlen(file);

	if(n<5)
		return NULL;

	if(!strcmp(&file[n-4], ".ico"))
		return "image/x-icon";

	if(!strcmp(&file[n-4], ".png"))
		return "image/png";

	if(!strcmp(&file[n-5], ".html"))
		return "text/html";

	return NULL;
}

/* this protocol server (always the first one) just knows how to do HTTP */

static int callback_http(struct lws* wsi,
	enum lws_callback_reasons reason,
	void* user,
	void* in,
	size_t len){
	(void)user;
	char buf[256];
	char leaf_path[1024];
	char b64[64];
	struct timeval tv;
	char* other_headers;
	const char* mimetype;
	unsigned char* p;
	int n;
	switch(reason){
	case LWS_CALLBACK_HTTP:
		dump_handshake_dbg(wsi);
		if(len<1){
			lws_return_http_status(wsi, HTTP_STATUS_BAD_REQUEST, NULL);
			return -1;
		}

		/* this server has no concept of directories */
		if(strchr((const char*)in+1, '/')){
			lws_return_http_status(wsi, HTTP_STATUS_FORBIDDEN, NULL);
			return -1;
		}

		/* if a legal POST URL, let it continue and accept data */
		if(lws_hdr_total_length(wsi, WSI_TOKEN_POST_URI))
			return 0;

		/* if not, send a file the easy way */
		strcpy(buf, resource_path);
		if(strcmp((const char*)in, "/")){
			strcat(buf, "/");
			strncat(buf, (char*)in, sizeof(buf)-strlen(resource_path));
		} else{ /* default file to serve */
			strcat(buf, "/monitor.html");
		}
		buf[sizeof(buf)-1]='\0';

		/* refuse to serve files we don't understand */
		mimetype=get_mimetype(buf);
		if(!mimetype){
			lwsl_err("Unknown mimetype for %s\n", buf);
			lws_return_http_status(wsi, HTTP_STATUS_UNSUPPORTED_MEDIA_TYPE, NULL);
			return -1;
		}

		/* demostrates how to set a cookie on / */

		other_headers=leaf_path;
		p=(unsigned char*)leaf_path;
		if(!strcmp((const char*)in, "/")&&
			!lws_hdr_total_length(wsi, WSI_TOKEN_HTTP_COOKIE)){
			/* this isn't very unguessable but it'll do for us */
			gettimeofday(&tv, NULL);
			n=sprintf(b64, "LWS_%u_%u_COOKIE; Max-Age=360000",
				(unsigned int)tv.tv_sec,
				(unsigned int)tv.tv_usec);
			if(lws_add_http_header_by_name(wsi,
				(unsigned char*)"set-cookie:",
				(unsigned char*)b64, n, &p,
				(unsigned char*)leaf_path+sizeof(leaf_path)))
				return 1;
		}
		n=(char*)p-leaf_path;
		n=lws_serve_http_file(wsi, buf, mimetype, other_headers, n);
		if(n<0||((n>0)&&lws_http_transaction_completed(wsi)))
			return -1; /* error or can't reuse connection: close the socket */

		break;

	case LWS_CALLBACK_HTTP_BODY:
		strncpy(buf, (char*)in, 20);
		buf[20]='\0';
		if(len<20)
			buf[len]='\0';

		lwsl_notice("LWS_CALLBACK_HTTP_BODY: %s... len %d\n",
			(const char*)buf, (int)len);

		break;

	case LWS_CALLBACK_HTTP_BODY_COMPLETION:
		lwsl_notice("LWS_CALLBACK_HTTP_BODY_COMPLETION\n");
		/* the whole of the sent body arried, close the connection */
		lws_return_http_status(wsi, HTTP_STATUS_OK, NULL);

		return -1;

	case LWS_CALLBACK_HTTP_FILE_COMPLETION:
	/* kill the connection after we sent one file */
		return -1;

	case LWS_CALLBACK_HTTP_WRITEABLE:

		return -1;

		/*
		 * callback for confirming to continue with client IP appear in
		 * protocol 0 callback since no websocket protocol has been agreed
		 * yet.  You can just ignore this if you won't filter on client IP
		 * since the default uhandled callback return is 0 meaning let the
		 * connection continue.
		 */

	case LWS_CALLBACK_FILTER_NETWORK_CONNECTION:
#if 0
		lws_get_peer_addresses(context, wsi, (int)(long)in, client_name,
			sizeof(client_name), client_ip, sizeof(client_ip));

		fprintf(stderr, "Received network connect from %s (%s)\n",
			client_name, client_ip);
#endif
	/* if we returned non-zero from here, we kill the connection */
		break;

	case LWS_CALLBACK_GET_THREAD_ID:
	/*
	 * if you will call "lws_callback_on_writable"
	 * from a different thread, return the caller thread ID
	 * here so lws can use this information to work out if it
	 * should signal the poll() loop to exit and restart early
	 */

	/* return pthread_getthreadid_np(); */

		break;

	default:
		break;
	}

	return 0;
}

struct a_message{
	char* payload;
	size_t len;//valid payload
	size_t size;//total memory allocated
};

#define MAX_MESSAGE_QUEUE 1024
/*The receiver modifies the head of the ring buffer. Each client maintains its
  own tail of the ring buffer so they can proceed that their own speed. */
static struct a_message ringbuffer[MAX_MESSAGE_QUEUE]={0};
static int ringbuffer_head=0;//ring buffer head for all clients
struct per_session_data__maos_monitor{//per client
	struct lws* wsi;
	l_message* head;//head of list just for initialization
	l_message* tail;//tail of list just for initialization
	int ringbuffer_tail;//ring buffer tail per client
};
int nclient=0;//count number of clients
int pending=0;//pending message for sending
static int
callback_maos_monitor(struct lws* wsi,
	enum lws_callback_reasons reason,
	void* user, void* in, size_t len){
	int n;
	struct per_session_data__maos_monitor* pss=(struct per_session_data__maos_monitor*)user;

	switch(reason){

	case LWS_CALLBACK_ESTABLISHED:
		nclient++;
		lwsl_info("callback_maos_monitor: LWS_CALLBACK_ESTABLISHED, nclient=%d\n", nclient);
		pss->ringbuffer_tail=ringbuffer_head;//initialize to zero length
		pss->wsi=wsi;
		pss->head=pss->tail=0;
		html_push_all(&pss->head, &pss->tail,
			LWS_SEND_BUFFER_PRE_PADDING, LWS_SEND_BUFFER_POST_PADDING);
		lws_callback_on_writable(wsi);
		break;

	case LWS_CALLBACK_PROTOCOL_DESTROY:
		nclient--;
		lwsl_notice("mirror protocol cleaning up, nclient=%d\n", nclient);
		for(n=0; n<(int)(sizeof ringbuffer/sizeof ringbuffer[0]); n++)
			if(ringbuffer[n].payload)
				free(ringbuffer[n].payload);
		break;

	case LWS_CALLBACK_SERVER_WRITEABLE:
		do{
			pending=0;
			if(pss->head){/*initialization with list*/
				n=lws_write(wsi, (unsigned char*)
					pss->head->payload+
					LWS_SEND_BUFFER_PRE_PADDING,
					pss->head->len,
					LWS_WRITE_TEXT);
				if(n<0){
					lwsl_err("ERROR %d writing to mirror socket\n", n);
					return -1;
				} else if(n<(int)pss->head->len){
					lwsl_err("mirror partial write %d out of %ld, retry.\n", n, (long)pss->head->len);
				}else{
					l_message* tmp=pss->head;
					pss->head=pss->head->next;
					free(tmp->payload);
					free(tmp);
				}
				if(pss->head) pending=1;
			}else if(pss->ringbuffer_tail!=ringbuffer_head){//send using ring buffer tail
				n=lws_write(wsi, (unsigned char*)
					ringbuffer[pss->ringbuffer_tail].payload+
					LWS_SEND_BUFFER_PRE_PADDING,
					ringbuffer[pss->ringbuffer_tail].len,
					LWS_WRITE_TEXT);
				if(n<0){
					lwsl_err("ERROR %d writing to mirror socket\n", n);
					return -1;
				}
				if(n<(int)ringbuffer[pss->ringbuffer_tail].len){
					lwsl_err("mirror partial write %d of %ld, retry.\n",
						n, (long)ringbuffer[pss->ringbuffer_tail].len);
				}else{
					if(pss->ringbuffer_tail==(MAX_MESSAGE_QUEUE-1))
						pss->ringbuffer_tail=0;
					else
						pss->ringbuffer_tail++;
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

	case LWS_CALLBACK_RECEIVE:
		scheduler_handle_ws((char*)in, len);
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
	{
	"http-only",		/* name */
	callback_http,		/* callback */
	sizeof(struct per_session_data__http),	/* per_session_data_size */
	0,			/* max frame size / rx buffer */
	0,                      //id unused
	0,                      //user unused
	},
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

static struct lws_context* context=0;
int ws_start(short port){
	if(context) return 0;
	struct lws_context_creation_info info;
	memset(&info, 0, sizeof info);
	info.port=port;
	lwsl_notice("libwebsockets based server "
		"(C) Copyright 2010-2013 Andy Green <andy@warmcat.com> - "
		"licensed under LGPL2.1\n");
	info.protocols=protocols;
	/*this is critical. otherwise UTF-8 error in client.*/
	info.gid=-1;
	info.uid=-1;
	context=lws_create_context(&info);
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
//Calls by scheduler timeout periodically.
int ws_service(int waiting){
	static int ans=0;
	if(!context){
		return -1;
	}
	if(waiting){
		lws_callback_on_writable_all_protocol(context, protocols+1);
	}
	//int repeat=0;
	pending=0;
	/*returns immediately if no task is pending when timeout is 0 ms.*/
	if(!ans){
		do{
#if LWS_LIBRARY_VERSION_NUMBER > 3002000
			ans=lws_service(context, -1);//>=3.2 stable, need to use -1 for polling
#else		
			ans=lws_service(context, 0);//<3.2 stable, use 0 for polling, but may not work for version close to 3.2
#endif		
			if(pending && waiting){
				mysleep(0.001);//this sleep is necessary 
				//repeat++;
			}
		}while(pending && (waiting--)>0);//service all available requests
		/*if(repeat>1){
			lwsl_notice("repeat %d times.\n", repeat);
		}*/
	}
	if(ans<0){
		ws_end();
	}
	return ans?-1:pending;
}
//Initiate server to client message with ring buffer.
void ws_push(const char* in, size_t len){
	if(!context || !nclient) return;
	size_t new_size=LWS_SEND_BUFFER_PRE_PADDING+len+
		LWS_SEND_BUFFER_POST_PADDING;
	//make sure buffer is big enough.
	if(ringbuffer[ringbuffer_head].size<new_size){
		ringbuffer[ringbuffer_head].payload=realloc(
			ringbuffer[ringbuffer_head].payload, new_size);
		ringbuffer[ringbuffer_head].size=new_size;
	}
	ringbuffer[ringbuffer_head].len=len;
	memcpy((char*)ringbuffer[ringbuffer_head].payload+
		LWS_SEND_BUFFER_PRE_PADDING, in, len);
	if(ringbuffer_head+1==MAX_MESSAGE_QUEUE){
		ringbuffer_head=0;/*wrap over*/
	} else{
		ringbuffer_head++;
	}
	lws_callback_on_writable_all_protocol(context, protocols+1);
}
