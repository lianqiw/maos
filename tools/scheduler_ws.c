


#include <getopt.h>


#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <syslog.h>
#include <sys/time.h>

#include <libwebsockets.h>
#include "../sys/sys.h"
#include "scheduler_ws.h"
/**
   Contains code for handling websocket connection using libwebsockets. Part of
   the code is copied from test-server.c in libwebsockets
*/

#define LOCAL_RESOURCE_PATH SRCDIR"/tools"
const char *resource_path = LOCAL_RESOURCE_PATH;
/*
 * We take a strict whitelist approach to stop ../ attacks
 */
struct serveable {
    const char *urlpath;
    const char *mimetype;
}; 

struct per_session_data__http {
    int fd;
};

/*
 * this is just an example of parsing handshake headers, you don't need this
 * in your code unless you will filter allowing connections by the header
 * content
 */
static void
dump_handshake_info(struct libwebsocket *wsi)
{
    static const char *token_names[] = {
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

    for (size_t n = 0; n < sizeof(token_names) / sizeof(token_names[0]); n++) {
	if (!lws_hdr_total_length(wsi, (lws_token_indexes)n))
	    continue;

	lws_hdr_copy(wsi, buf, sizeof buf, (lws_token_indexes)n);

	fprintf(stderr, "    %s = %s\n", token_names[n], buf);
    }
}

const char * get_mimetype(const char *file)
{
    int n = strlen(file);

    if (n < 5)
	return NULL;

    if (!strcmp(&file[n - 4], ".ico"))
	return "image/x-icon";

    if (!strcmp(&file[n - 4], ".png"))
	return "image/png";

    if (!strcmp(&file[n - 5], ".html"))
	return "text/html";

    return NULL;
}

/* this protocol server (always the first one) just knows how to do HTTP */

static int callback_http(struct libwebsocket_context *context,
			 struct libwebsocket *wsi,
			 enum libwebsocket_callback_reasons reason, 
			 void *user,
			 void *in, 
			 size_t len)
{
    (void)user;
    char buf[256];
    char leaf_path[1024];
    char b64[64];
    struct timeval tv;
    char *other_headers;
    const char *mimetype;

    switch (reason) {
    case LWS_CALLBACK_HTTP:
	dump_handshake_info(wsi);
	if (len < 1) {
	    libwebsockets_return_http_status(context, wsi,
					     HTTP_STATUS_BAD_REQUEST, NULL);
	    return -1;
	}

	/* this server has no concept of directories */
	if (strchr((const char *)in + 1, '/')) {
	    libwebsockets_return_http_status(context, wsi,
					     HTTP_STATUS_FORBIDDEN, NULL);
	    return -1;
	}

	/* if a legal POST URL, let it continue and accept data */
	if (lws_hdr_total_length(wsi, WSI_TOKEN_POST_URI))
	    return 0;

	/* if not, send a file the easy way */
	strcpy(buf, resource_path);
	if (strcmp((const char*)in, "/")) {
	    strcat(buf, "/");
	    strncat(buf, (char*)in, sizeof(buf) - strlen(resource_path));
	}else{ /* default file to serve */
	    strcat(buf, "/monitor.html");
	}
	buf[sizeof(buf) - 1] = '\0';

	/* refuse to serve files we don't understand */
	mimetype = get_mimetype(buf);
	if (!mimetype) {
	    lwsl_err("Unknown mimetype for %s\n", buf);
	    libwebsockets_return_http_status(context, wsi, HTTP_STATUS_UNSUPPORTED_MEDIA_TYPE, NULL);
	    return -1;
	}

	/* demostrates how to set a cookie on / */

	other_headers = NULL;
	if (!strcmp((const char *)in, "/") &&
	    !lws_hdr_total_length(wsi, WSI_TOKEN_HTTP_COOKIE)) {
	    /* this isn't very unguessable but it'll do for us */
	    gettimeofday(&tv, NULL);
	    sprintf(b64, "LWS_%u_%u_COOKIE",
		    (unsigned int)tv.tv_sec,
		    (unsigned int)tv.tv_usec);

	    sprintf(leaf_path,
		    "Set-Cookie: test=LWS_%u_%u_COOKIE;Max-Age=360000\x0d\x0a",
		    (unsigned int)tv.tv_sec, (unsigned int)tv.tv_usec);
	    other_headers = leaf_path;
	    lwsl_err(other_headers);
	}

	if (libwebsockets_serve_http_file(context, wsi, buf, mimetype, other_headers)){
	    return -1; /* through completion or error, close the socket */
	}
	break;

    case LWS_CALLBACK_HTTP_BODY:
	strncpy(buf, (char*)in, 20);
	buf[20] = '\0';
	if (len < 20)
	    buf[len] = '\0';

	lwsl_notice("LWS_CALLBACK_HTTP_BODY: %s... len %d\n",
		    (const char *)buf, (int)len);

	break;

    case LWS_CALLBACK_HTTP_BODY_COMPLETION:
	lwsl_notice("LWS_CALLBACK_HTTP_BODY_COMPLETION\n");
	/* the whole of the sent body arried, close the connection */
	libwebsockets_return_http_status(context, wsi,
					 HTTP_STATUS_OK, NULL);

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
	libwebsockets_get_peer_addresses(context, wsi, (int)(long)in, client_name,
					 sizeof(client_name), client_ip, sizeof(client_ip));

	fprintf(stderr, "Received network connect from %s (%s)\n",
		client_name, client_ip);
#endif
	/* if we returned non-zero from here, we kill the connection */
	break;

    case LWS_CALLBACK_GET_THREAD_ID:
	/*
	 * if you will call "libwebsocket_callback_on_writable"
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

struct a_message {
    char *payload;
    size_t len;
};
#define MAX_MESSAGE_QUEUE 64
/*The receiver modifies the head of the ring buffer. Each client maintains its
  own tail of the ring buffer so they can proceed that their own speed. */
static struct a_message ringbuffer[MAX_MESSAGE_QUEUE];
static int ringbuffer_head;
struct per_session_data__maos_monitor {
    struct libwebsocket *wsi;
    l_message *head;/*for initialization*/
    l_message *tail;
    int ringbuffer_tail;
};

static int
callback_maos_monitor(struct libwebsocket_context *context,
		      struct libwebsocket *wsi,
		      enum libwebsocket_callback_reasons reason,
		      void *user, void *in, size_t len)
{
    int n;
    struct per_session_data__maos_monitor *pss = (struct per_session_data__maos_monitor *)user;

    switch (reason) {

    case LWS_CALLBACK_ESTABLISHED:
	lwsl_info("callback_maos_monitor: LWS_CALLBACK_ESTABLISHED\n");
	pss->ringbuffer_tail = ringbuffer_head;
	pss->wsi = wsi;
	pss->head=pss->tail=0;
	html_convert_all(&pss->head, &pss->tail, 
			 LWS_SEND_BUFFER_PRE_PADDING, LWS_SEND_BUFFER_POST_PADDING);
	libwebsocket_callback_on_writable(context, wsi);
	lwsl_notice("head=%p, tail=%p\n", pss->head, pss->tail);
	break;

    case LWS_CALLBACK_PROTOCOL_DESTROY:
	lwsl_notice("mirror protocol cleaning up\n");
	for (n = 0; n < (int)(sizeof ringbuffer / sizeof ringbuffer[0]); n++)
	    if (ringbuffer[n].payload)
		free(ringbuffer[n].payload);
	break;

    case LWS_CALLBACK_SERVER_WRITEABLE:
	while(pss->head){/*initialization*/
	    n = libwebsocket_write(wsi, (unsigned char *)
				   pss->head->payload +
				   LWS_SEND_BUFFER_PRE_PADDING,
				   pss->head->len,
				   LWS_WRITE_TEXT);
	    if(n<0){
		lwsl_err("ERROR %d writing to mirror socket\n", n);
		return -1;
	    }else if(n<(int)pss->head->len){
		lwsl_err("mirror partial write %d vs %d\n", n, pss->head->len);
	    }
	    l_message *tmp=pss->head;
	    pss->head=pss->head->next;
	    free(tmp->payload);
	    free(tmp);
	    if(lws_send_pipe_choked(wsi)){
		libwebsocket_callback_on_writable(context, wsi);
		break;
	    }
#ifdef _WIN32
	    Sleep(1);
#else
	    usleep(1);
#endif
	}
	while (pss->ringbuffer_tail != ringbuffer_head) {
	    n = libwebsocket_write(wsi, (unsigned char *)
				   ringbuffer[pss->ringbuffer_tail].payload +
				   LWS_SEND_BUFFER_PRE_PADDING,
				   ringbuffer[pss->ringbuffer_tail].len,
				   LWS_WRITE_TEXT);
	    if (n < 0) {
		lwsl_err("ERROR %d writing to mirror socket\n", n);
		return -1;
	    }
	    if (n < (int)ringbuffer[pss->ringbuffer_tail].len)
		lwsl_err("mirror partial write %d vs %d\n",
			 n, ringbuffer[pss->ringbuffer_tail].len);

	    if (pss->ringbuffer_tail == (MAX_MESSAGE_QUEUE - 1))
		pss->ringbuffer_tail = 0;
	    else
		pss->ringbuffer_tail++;


	    if (lws_send_pipe_choked(wsi)) {
		libwebsocket_callback_on_writable(context, wsi);
		break;
	    }
	    /*
	     * for tests with chrome on same machine as client and
	     * server, this is needed to stop chrome choking
	     */
#ifdef _WIN32
	    Sleep(1);
#else
	    usleep(1);
#endif
	}
	break;

    case LWS_CALLBACK_RECEIVE:
	scheduler_handle_ws((char*)in, len);
	//ws_push(in, len);
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

static struct libwebsocket_protocols protocols[] = {
    /* first protocol must always be HTTP handler */
    {
	"http-only",		/* name */
	callback_http,		/* callback */
	sizeof (struct per_session_data__http),	/* per_session_data_size */
	0,			/* max frame size / rx buffer */
    },
    {
	"maos-monitor-protocol",
	callback_maos_monitor,
	sizeof(struct per_session_data__maos_monitor),
	1024,                   /* max frame size / rx buffer */
    },
    { NULL, NULL, 0, 0 } /* terminator */
};
static struct libwebsocket_context *context=0;
int ws_start(short port){
    struct lws_context_creation_info info;
    memset(&info, 0, sizeof info);
    info.port = port;
    lwsl_notice("libwebsockets test server - "
		"(C) Copyright 2010-2013 Andy Green <andy@warmcat.com> - "
		"licensed under LGPL2.1\n");
    info.protocols = protocols;
    info.extensions = libwebsocket_get_internal_extensions();
    /*this is critical. otherwise UTF-8 error in client.*/
    info.gid = -1;
    info.uid = -1;
    info.ka_time=100;
    context = libwebsocket_create_context(&info);
    if (!context){
	lwsl_err("libwebsocket init failed\n");
	return -1;
    }
    return 0;
}
void ws_end(){
    if(context){
	libwebsocket_context_destroy(context);
	context=0;
	lwsl_notice("libwebsockets-test-server exited cleanly\n");
    }
}
int ws_service(){
    /*returns immediately if no task is pending.*/
    if(context){
	int ans=libwebsocket_service(context, 0);
	if(ans<0){
	    ws_end();
	    context=0;
	}
	return ans;
    }else{
	return -1;
    }
}
void ws_push(const char *in, int len){
    if(!context) return;
    free(ringbuffer[ringbuffer_head].payload);
    ringbuffer[ringbuffer_head].payload =(char*)
	malloc(LWS_SEND_BUFFER_PRE_PADDING + len +
	       LWS_SEND_BUFFER_POST_PADDING);
    ringbuffer[ringbuffer_head].len = len;
    memcpy((char *)ringbuffer[ringbuffer_head].payload +
	   LWS_SEND_BUFFER_PRE_PADDING, in, len);
    if (ringbuffer_head == (MAX_MESSAGE_QUEUE - 1)){
	ringbuffer_head = 0;/*wrap over*/
    }else{
	ringbuffer_head++;
    }
    libwebsocket_callback_on_writable_all_protocol(protocols+1);
}
