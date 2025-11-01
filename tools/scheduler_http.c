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
	\file scheduler_http.c This file contains custom web server routines to
	handle http connection, serving files, upgrade the connection to websocket,
	and read and write websocket messages.
*/
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#ifdef HAVE_OPENSSL 
#include <openssl/ssl.h>
#include <openssl/err.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <poll.h>
#include <arpa/inet.h>
#include <stdint.h>
#include <sys/stat.h>
#include "../sys/sys.h"
#include "scheduler.h"
#define BUF_SIZE 4096
/**
 * @brief Context for an HTTP connection.
 *
 * Encapsulates transport-specific state and I/O operations for a single
 * HTTP connection so callers can treat plain and TLS-backed sockets
 * uniformly.
 *
 * Responsibilities and semantics:
 *  - Holds the underlying OS file descriptor for the connection.
 *  - Optionally holds a TLS/SSL handle when built with OpenSSL.
 *  - Provides function pointers for receiving and sending data. These
 *    functions should follow conventional POSIX-style semantics:
 *      - Return the number of bytes read/written on success.
 *      - Return 0 to indicate EOF (for recv) where applicable.
 *      - Return -1 on error and set errno (or use a transport-specific
 *        error reporting mechanism).
 *  - The caller is responsible for managing the lifecycle of the fd and
 *    the ssl pointer (closing/freeing them as appropriate).
 *
 * Thread-safety:
 *  - Access to a single http_context_t instance is not guaranteed to be
 *    thread-safe. Synchronize externally if the same context will be used
 *    concurrently from multiple threads.
 *
 * Members:
 *  - fd: File descriptor for the underlying socket/connection.
 *  - ssl: (present only when HAVE_OPENSSL is enabled) Pointer to an
 *         OpenSSL SSL object or equivalent TLS state.
 *  - recv: Pointer to a receive function with signature:
 *          int (*recv)(struct http_context_t *ctx, char *buf, size_t nbuf);
 *          Must read up to nbuf bytes into buf.
 *  - send: Pointer to a send function with signature:
 *          int (*send)(struct http_context_t *ctx, const char *buf, size_t nbuf);
 *          Must write up to nbuf bytes from buf.
 */
typedef struct http_context_t{
	int fd;
#if HAVE_OPENSSL
	SSL *ssl;
#endif	
	int(*recv)(struct http_context_t*ctx, char *buf, size_t nbuf);
	int(*send)(struct http_context_t*ctx, const char *buf, size_t nbuf);
}http_context_t;
http_context_t *ssl_context=NULL;
int nssl_context=0;
#ifdef HAVE_OPENSSL
static int has_ssl=0;//guard initialization and tell the status. 1: yes. 0 or -1: no. 0: not initialized
SSL_CTX *ssl_ctx=NULL;
void init_ssl(){
	if(has_ssl) return;
	SSL_load_error_strings();
	OpenSSL_add_ssl_algorithms();

	const SSL_METHOD *method = TLS_server_method();  // TLS 1.2/1.3 compatible
	ssl_ctx = SSL_CTX_new(method);
	if (!ssl_ctx) {
		perror("Unable to create SSL context");
		ERR_print_errors_fp(stderr);
		has_ssl=-1;
	}else{
		has_ssl=1;
		do{
			char fn_crt[PATH_MAX];
			char fn_key[PATH_MAX];
			mysnprintf(fn_crt, sizeof(fn_crt), "%s/.aos/%s.crt", HOME, HOST);
			if(!exist(fn_crt)){
				mysnprintf(fn_crt, sizeof(fn_crt), "%s/.aos/Server.crt", HOME);
			}
			if(!exist(fn_crt)){
				warning_time("TLS certificate not found at %s\n", fn_crt);
				has_ssl=-1; break;
			}
			mysnprintf(fn_key, sizeof(fn_key), "%s/.aos/%s.key", HOME, HOST);
			if(!exist(fn_key)){
				mysnprintf(fn_key, sizeof(fn_key), "%s/.aos/Server.key", HOME);
			}
			if(!exist(fn_key)){
				warning_time("TLS key not found at %s\n", fn_key);
				has_ssl=-1; break;
			}

			if (SSL_CTX_use_certificate_file(ssl_ctx, fn_crt, SSL_FILETYPE_PEM) <= 0 ||
				SSL_CTX_use_PrivateKey_file(ssl_ctx, fn_key, SSL_FILETYPE_PEM) <= 0) {
				ERR_print_errors_fp(stderr);
				has_ssl=-1; break;
			}
		}while(0);
		if(has_ssl==-1){
			SSL_CTX_free(ssl_ctx);
			ssl_ctx=NULL;
		}
	}
}
#endif
static void http_context_remove(http_context_t *ctx){
	ctx->fd=-1;
#if HAVE_OPENSSL
	if(ssl_ctx && ctx->ssl){
		SSL_free(ctx->ssl);
		ctx->ssl=NULL;
	}
#endif
}
#if HAVE_OPENSSL
static int https_reader(http_context_t*ctx, char *buf, size_t nbuf){
	return ctx?SSL_read(ctx->ssl, (void*)buf, nbuf):-1;
}
static int https_writer(http_context_t*ctx, const char *buf, size_t nbuf){
	if(!ctx) return -1;
	if(!nbuf) return 0;
	ssize_t nsend;
	do{
		nsend=SSL_write(ctx->ssl, buf, nbuf);
		//dbg_time("fd=%d nsend=%ld, nbuf=%lu\n", ctx->fd, nsend, nbuf);
		if(nsend>0){
			buf=buf+nsend;
			nbuf-=nsend;
		}
	}while(nsend>0 && nbuf);
	if(nbuf){
		http_context_remove(ctx);
		return -1;
	}else{
		return 0;
	}
}
#endif
static int http_reader(http_context_t*ctx, char *buf, size_t nbuf){
	return ctx?recv(ctx->fd, (void*)buf, nbuf, 0):-1;
}
static int http_writer(http_context_t*ctx, const char *buf, size_t nbuf){
	if(!ctx) return -1;
	if(!nbuf) return 0;
	ssize_t nsend;
	do{
		nsend=send(ctx->fd, buf, nbuf, 0);
		//dbg_time("fd=%d nsend=%ld, nbuf=%lu\n", ctx->fd, nsend, nbuf);
		if(nsend>0){
			buf=buf+nsend;
			nbuf-=nsend;
		}
	}while(nsend>0 && nbuf);
	if(nbuf){
		http_context_remove(ctx);
		return -1;
	}else{
		return 0;
	}
}

static http_context_t *http_context_get(int fd){
	for (int ic=0; ic<nssl_context; ic++){
		if(ssl_context[ic].fd==fd){
			return &ssl_context[ic];
		}
	}
	warning_time("http_context not found for fd=%d\n", fd);
	return NULL;
}
/**Create a context for fd. overriding existing one if exists */
static http_context_t *http_context_create(int fd){
	char buf[5]={0};
	int len=recv(fd, buf, 4, MSG_PEEK);
	int ishttp=(len==4 && !strcmp(buf, "GET "));
	int ishttps=(len>=3 && buf[0] == 0x16 && buf[1] == 0x03 &&
    		(buf[2] == 0x00 || buf[2] == 0x01 || buf[2] == 0x02 || buf[2] == 0x03));
	if(!ishttp && !ishttps){
		warning("Connection is not http or https\n");
		return NULL;
	}
	int jc=-1;
	for (int ic=0; ic<nssl_context; ic++){
		if(ssl_context[ic].fd==fd){
			jc=ic;
			break;//found matching
		}else if(ssl_context[ic].fd==-1 && jc==-1){
			jc=ic;//mark for insertion
		}
	}
	if(jc==-1){
		jc=nssl_context;
		nssl_context++;
		ssl_context=myrealloc(ssl_context, nssl_context, http_context_t);
	}
	memset(&ssl_context[jc], 0, sizeof(http_context_t));
	ssl_context[jc].fd=fd;
	if(ishttps){
#if HAVE_OPENSSL
		if(has_ssl==0) init_ssl();
		if(ssl_ctx){
			SSL* ssl=SSL_new(ssl_ctx);
			SSL_set_fd(ssl, fd);
			if (SSL_accept(ssl) > 0) {//success
				ssl_context[jc].ssl=ssl;
				ssl_context[jc].recv=https_reader;
				ssl_context[jc].send=https_writer;
			}else{
				unsigned long err = ERR_peek_error();
				warning_time("SSL_accept failed: %s\n", ERR_reason_error_string(err));
				SSL_free(ssl);
				ssl=NULL;
			}
		}
#endif
		if(!ssl_context[jc].recv){
			warning("HTTPS is not supported\n");
			return NULL;
		}
	}else if(ishttp){
		ssl_context[jc].recv=http_reader;
		ssl_context[jc].send=http_writer;
	}
	return &ssl_context[jc];
}


const int http_padding=16;
/** 
 * Write to websockets. allocated memory must have at least http_padding padding before the pointer*/
static int ws_send(char *buf, int nlen, int mode, void *userdata){
	int offset=2;//no mask field can be present
	if(nlen>UINT16_MAX){
		offset+=8;
	}else if(nlen>125){
		offset+=2;
	}
	char *header=buf-offset;
	header[0] = 0x81+mode;//82 is binary. 81 is text. 
	if(nlen<126){
		header[1]=nlen;
	}else if(nlen>UINT16_MAX){//8 byte length
		header[1]=127;
		*((uint64_t*)&header[2])=swap8bytes((long)nlen);
	}else{//2 byte length
		header[1]=126;
		*((uint16_t*)&header[2])=swap2bytes(nlen);
	}
	int fd_dest=(int)(long)userdata;
	http_context_t* ctx=http_context_get(fd_dest);
	if(!ctx || ctx->send(ctx, header, nlen+offset)){
		return -1;
	}
	return 0;
}
/** 
 * Read tcp data and send to websockets */
static int ws_forward(int fd, int nlen, int mode, void *userdata){
	static char *buf=0;//persistent buffer
	static int size=0;
	int offset=2;//no mask field can be present
	if(nlen>UINT16_MAX){
		offset+=8;
	}else if(nlen>125){
		offset+=2;
	}
	if(size<offset+nlen){
		size=offset+nlen;
		buf=myrealloc(buf, size, char);
	}
	if(stread(fd, buf+offset, nlen)) return -1;
	return ws_send(buf+offset, nlen, mode, userdata);
}
static void ws_write_close(http_context_t* ctx){//send a ws close frame.
	if(!ctx) return;
	char buf[4];
	buf[0]=0x88;//FIN=1; opcode=8
	buf[1]=0x2;//payload=2
	buf[2]=0x3;	buf[3]=0xE8;//1000
	if(ctx->send(ctx, buf, sizeof(buf))){
		dbg_time("Send close failed.\n");
	}
}
static int ws_fail(int fd, const char *format, ...) CHECK_ARG(2);
static int ws_fail(int fd, const char *format, ...) {
	format2fn;
	dbg_time("sock %d: %s\n", fd, fn);
	ws_proxy_remove((void*)(long)fd, 1);
	return -1;
}
/*
	Receive websocket messages and call ws_proxy_command to handle it.
*/
static int ws_receive(struct pollfd *pfd, int flag){
	const int fd=pfd->fd;
	http_context_t* ctx=http_context_get(pfd->fd);
	if(!ctx){
		return -1;
	}
	if(flag==-1){//server ask browser client to close
		ws_write_close(ctx);//send browser client close frame
		shutdown(fd, SHUT_WR);//make sure we don't reply to the close message.
		return ws_fail(fd, "listen_sock requests shutdown");
	}	
	static char *prev=NULL;//save previous message to handle continuation
	static size_t nprev=0;//length of previous message
	static int oprev=0;//opcode of previous message
    char buf0[BUF_SIZE];
	char *buf=buf0;
	const int len = ctx->recv(ctx, buf, BUF_SIZE-1);
	if (len<6){//server received msg must be at least 6 bytes
		return ws_fail(fd, "len=%d", len);
	}
	const int FIN=((unsigned char)(buf[0])>>7);//must be unsigned char for this to be 1.
	int opcode=(unsigned char)(buf[0]) & 0x0F;
	int offset=6;//offset of payload from buf for small message (2 byte opcode+msg size, 4 byte mask)
	size_t msg_len = buf[1] & 0x7F;
	if(msg_len==127){
		msg_len=swap8bytes(*((uint64_t*)&buf[2]));
		offset+=8;//additional offset
	}else if(msg_len==126){
		msg_len=swap2bytes(*((uint16_t*)&buf[2]));
		offset+=2;
	}
	if(msg_len+offset>BUF_SIZE-1){//buffer is too small. We read the remaining data in block mode.
		//dbg_time("Buffer is too small\n");
		buf=(char*)malloc(msg_len+offset+1);
		memcpy(buf, buf0, BUF_SIZE-1);
		if(ctx->recv(ctx, buf+BUF_SIZE-1, msg_len+offset-(BUF_SIZE-1))){
			return ws_fail(fd, "recv");
		}
	}
	buf[msg_len+offset]=0;//always terminate string
	/*Unmask message from client. A client must mask all frames sent to the
	server. A server must not mask any frames sent to the client*/
	char mask[4];
	memcpy(mask, &buf[offset-4], 4);
	for (size_t i = 0; i < msg_len; i++){
		buf[offset + i] ^= mask[i % 4];
	}
	//dbg_time("FIN=%d, opcode=%d, prev=%p, nprev=%ld, msg_len=%ld\n", FIN, opcode, prev, nprev, msg_len);
	if(FIN==0 || (FIN==1 && opcode==0)){//there is segmentation
		//first segment: FIN=0, opcode=1 or 2
		//middle segment: FIN=0, opcode=0
		//final segment: FIN=1, opcode=0
		
		if(opcode!=0){//first segment
			oprev=opcode;//save opcode
			if(prev || nprev){
				warning("unexpected. free prev\n");
				free(prev); prev=NULL; nprev=0;
			}
		}
		prev=myrealloc(prev, (nprev+msg_len+1), char);//save data
		memcpy(prev+nprev, buf+offset, msg_len+1); nprev+=msg_len;//save length
		offset=0;//in prev, no offset
		if(FIN==1 && opcode==0){//assemble the final message
			if(buf!=buf0){
				free(buf);
			}
			opcode=oprev;
			buf=prev;
			msg_len=nprev;
			prev=NULL;
			nprev=0;
		}
	}else if(opcode==0x8){//client replied//request close
		ws_write_close(ctx);
		return ws_fail(fd, "client requests close");
	}else if(opcode==0x9){//client sent a ping. we reply with pong
		buf[offset-2]=0x8A;//pong
		buf[offset-1]=msg_len;//cannot be more than 126
		if(ctx->send(ctx, buf+offset-2, msg_len+2)){
			return ws_fail(fd, "pong");
		}else{
			return 0;
		}
	}
	if(FIN==1){
		if(opcode==0x1){
			ws_proxy_command(buf+offset, msg_len, (ws_proxy_t){.forward=ws_forward, .userdata=(void*)(long)fd});
		}else if(opcode==0x2){
			warning_time("opcode=%d is not handled\n", opcode);
		}
	}
	if(buf!=buf0){//heap memory is allocated.
		free(buf);
	}
    return 0;
}
/*
	Parse keys from http payload. The client needs to free the returned memory.
*/
char* http_parse_key(const char *buf, const char *name){
	const char *key_start = strstr(buf, name);
	if (!key_start) return 0;
	key_start += strlen(name);
	while (*key_start == ' ') key_start++;
	const char *key_end = strchr(key_start, '\r');
	if (!key_end) return 0;
	return strndup(key_start, key_end-key_start);
}
/*
	Upgrade http connection to websocket. it replaces the handler in listen_port.
*/
void http_upgrade_websocket(http_context_t *ctx, char *buf) {
	//Sec-WebSocket-Protocol: maos-monitor-protocol
	//dbg2_time("Upgrading %d to websockets\n", fd);
	char *key=http_parse_key(buf, "Sec-WebSocket-Key:");
	char *protocol=http_parse_key(buf, "Sec-WebSocket-Protocol:");
    char accept_src[128];
	char accept_key[64];
    snprintf(accept_src, sizeof(accept_src), "%s258EAFA5-E914-47DA-95CA-C5AB0DC85B11", key);
	base64_sha1(accept_src, accept_key);
    char resp[256];
    snprintf(resp, sizeof(resp),
             "HTTP/1.1 101 Switching Protocols\r\n"
             "Upgrade: websocket\r\n"
             "Connection: Upgrade\r\n"
             "Sec-WebSocket-Accept: %s\r\n"
			 "Sec-WebSocket-Protocol: %s\r\n\r\n",
             accept_key, protocol?protocol:"chat");
    ctx->send(ctx, resp, strlen(resp));
	listen_port_add(ctx->fd, POLLIN, ws_receive, "ws_receive");//switch to websocket handler
	if(!strcmp(protocol, "maos-monitor-protocol")){
		char cmd[11]="1&MONITOR;";
		ws_proxy_command(cmd, strlen(cmd), (ws_proxy_t){.send=ws_send, .userdata=(void*)(long)(ctx->fd)});
	}
	
	free(protocol);
	free(key);
}
/*
	Send a close request to client which will indicate closure. 
*/
void http_close(http_context_t *ctx){
	if(!ctx) return;
	const char *resp =
		"HTTP/1.1 200 OK\r\n"
		"Content-Type: text/plain\r\n"
		"Content-Length: 12\r\n"
		"Connection: close\r\n"
		"\r\n"
		"Hello world\n";

	ctx->send(ctx, resp, strlen(resp));
	shutdown(ctx->fd, SHUT_WR);  // optional graceful signal
}

/*
	HTTP handler callable by listen_port. 
	It supports sending plain files (folder is disabled for security reasons) and websocket upgrade
*/
int http_handler(struct pollfd *pfd, int flag){
	int fd=pfd->fd;
	http_context_t* ctx=http_context_get(pfd->fd);
	if (!ctx) {
		warning_time("http_context_get returned NULL for fd=%d\n", fd);
		return -1;
	}
	if(flag==-1){
		http_close(ctx);
		return 0;//wait for client to initiate close
	}
	char buf[BUF_SIZE];
	
	int n= ctx->recv(ctx, buf, sizeof(buf) - 1);
	if (n <= 0) return -1;
	buf[n] = 0;
	//info_time("received '%s'\n", buf);
    if (strncmp(buf, "GET ", 4) != 0){
		warning_time("buf shall start with 'GET': '%s'. close connection. fd=%d\n", buf, fd);
		return -1;//error
	}
    if (strstr(buf, "Upgrade: websocket")) {//handle upgrade
        http_upgrade_websocket(ctx, buf);
		return 0;
    } else {//handle http send file
		int ans=0;
		const char *path = buf + 4;
		const char* close=strstr(buf, "Connection: close")?"Connection: close\r\n":"";
		char *sp = strchr(buf+4, ' ');
		if (!sp) return -1;
		*sp = 0;
		if (strcmp(path, "/") == 0){
        	path = "/monitor.html";
		}
		if(path[0]=='/'){
			path++;//ignore leading /
		}
		const char *texttype=NULL;
		char fullpath[256];
		FILE *fp=NULL;
		if(strchr(path, '/')){//do not allow directory. return 404.
			//warning("directory is not supported in path: '%s'", path);
			path=NULL;
		}else{
			if(check_suffix(path, ".html") || check_suffix(path, ".htm")){
				texttype="text/html";
			}else if(check_suffix(path, ".css")){
				texttype="text/css";
			}else if(check_suffix(path, ".js")||check_suffix(path, ".jsx")){
				texttype="text/javascript";
			}else if(check_suffix(path, ".jpg")||check_suffix(path, ".jpeg")){
				texttype="image/jpeg";
			}else if(check_suffix(path, ".png")){
				texttype="image/png";
			}
			if(texttype){
        		snprintf(fullpath, sizeof(fullpath), "%s/%s", SRCDIR "/tools", path);	
				fp = fopen(fullpath, "rb");
			}
		}
		char header[256];

		if (!fp) {//send 404
			snprintf(header, sizeof(header), "HTTP/1.1 404 Not Found\r\n%sContent-Length: 9\r\n\r\nNot Found", close);
			ans=ctx->send(ctx, header, strlen(header));
		}else{//send file
			fseek(fp, 0, SEEK_END);
			long size = ftell(fp);
			rewind(fp);
			
			snprintf(header, sizeof(header), "HTTP/1.1 200 OK\r\n%sContent-Type: %s\r\nContent-Length: %ld\r\n\r\n", close, texttype, size);
			if(ctx->send(ctx, header, strlen(header))){
				ans=-1;
			}else{
				char buf2[BUF_SIZE];
				size_t n2;
				
				while ((n2 = fread(buf2, 1, sizeof(buf2), fp)) > 0 && !ans){
					ans=ctx->send(ctx, buf2, n2);
				}
			}
			fclose(fp);
		}
		if(close[0]!=0){//close connection
			ans=-1;
		}
		return ans;
		
	}
}
int http_handshake(struct pollfd *pfd, int flag){
	http_context_t *ctx=http_context_create(pfd->fd);
	if(!ctx){
		return -1;
	}
	if(flag==-1){
		http_close(ctx);
		return 0;//wait for client to initiate close
	}
	listen_port_add(pfd->fd, POLLIN, http_handler, "http_handler");
	return 0;
}
