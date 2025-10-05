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

/* Minimal SHA-1 Implementation */

typedef struct {
    uint32_t h[5];
    unsigned char block[64];
    uint64_t bitlen;
    size_t curlen;
} SHA1_CTX;

static void sha1_init(SHA1_CTX *ctx) {
    ctx->h[0] = 0x67452301;
    ctx->h[1] = 0xEFCDAB89;
    ctx->h[2] = 0x98BADCFE;
    ctx->h[3] = 0x10325476;
    ctx->h[4] = 0xC3D2E1F0;
    ctx->curlen = 0;
    ctx->bitlen = 0;
}

static void sha1_transform(SHA1_CTX *ctx, const unsigned char *data) {
    uint32_t w[80];
    uint32_t a, b, c, d, e, t;

    for (int i = 0; i < 16; i++)
        w[i] = (data[4*i]<<24)|(data[4*i+1]<<16)|(data[4*i+2]<<8)|data[4*i+3];
    for (int i = 16; i < 80; i++) {
        uint32_t v = w[i-3]^w[i-8]^w[i-14]^w[i-16];
        w[i] = (v << 1) | (v >> 31);
    }

    a = ctx->h[0]; b = ctx->h[1]; c = ctx->h[2]; d = ctx->h[3]; e = ctx->h[4];

    for (int i = 0; i < 80; i++) {
        uint32_t f, k;
        if (i < 20) { f = (b & c) | ((~b) & d); k = 0x5A827999; }
        else if (i < 40) { f = b ^ c ^ d; k = 0x6ED9EBA1; }
        else if (i < 60) { f = (b & c) | (b & d) | (c & d); k = 0x8F1BBCDC; }
        else { f = b ^ c ^ d; k = 0xCA62C1D6; }

        t = ((a << 5) | (a >> 27)) + f + e + k + w[i];
        e = d; d = c; c = (b << 30) | (b >> 2); b = a; a = t;
    }

    ctx->h[0] += a; ctx->h[1] += b; ctx->h[2] += c;
    ctx->h[3] += d; ctx->h[4] += e;
}

static void sha1_update(SHA1_CTX *ctx, const unsigned char *data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        ctx->block[ctx->curlen++] = data[i];
        if (ctx->curlen == 64) {
            sha1_transform(ctx, ctx->block);
            ctx->bitlen += 512;
            ctx->curlen = 0;
        }
    }
}

static void sha1_final(SHA1_CTX *ctx, unsigned char *out) {
    ctx->bitlen += ctx->curlen * 8;
    ctx->block[ctx->curlen++] = 0x80;
    if (ctx->curlen > 56) {
        while (ctx->curlen < 64) ctx->block[ctx->curlen++] = 0;
        sha1_transform(ctx, ctx->block);
        ctx->curlen = 0;
    }
    while (ctx->curlen < 56) ctx->block[ctx->curlen++] = 0;

    for (int i = 7; i >= 0; i--)
        ctx->block[ctx->curlen++] = (ctx->bitlen >> (i * 8)) & 0xFF;

    sha1_transform(ctx, ctx->block);
    for (int i = 0; i < 5; i++) {
        out[i*4] = (ctx->h[i] >> 24) & 0xFF;
        out[i*4+1] = (ctx->h[i] >> 16) & 0xFF;
        out[i*4+2] = (ctx->h[i] >> 8) & 0xFF;
        out[i*4+3] = ctx->h[i] & 0xFF;
    }
}

/**  Base64 Encode  */
static void base64_encode(const unsigned char *in, int in_len, char *out) {
	static const char b64_table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int i, j;
    for (i = 0, j = 0; i < in_len;) {
        uint32_t a = i < in_len ? in[i] : 0; i++;
        uint32_t b = i < in_len ? in[i] : 0; i++;
        uint32_t c = i < in_len ? in[i] : 0; i++;
        uint32_t triple = (a << 16) | (b << 8) | c;

        out[j++] = b64_table[(triple >> 18) & 0x3F];
        out[j++] = b64_table[(triple >> 12) & 0x3F];
        out[j++] = (i > in_len + 1) ? '=' : b64_table[(triple >> 6) & 0x3F];
        out[j++] = (i > in_len) ? '=' : b64_table[triple & 0x3F];
    }
    out[j] = 0;
}

const int http_padding=16;
/** 
 * Read tcp data and send to websockets. allocated memory must have at least http_padding padding before the pointer*/
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
	//read tcp and send to websockets
	if(stwrite(fd_dest, header, nlen+offset)){
		return -1;
	}
	return 0;
}
/** 
 * Read tcp data and send to websockets */
static int ws_forward(int fd, int nlen, int mode, void *userdata){
	static char *buf=0;
	static int size=0;
	if(size<http_padding+nlen){
		size=http_padding+nlen;
		buf=realloc(buf, size);
	}
	int offset=2;//no mask field can be present
	buf[0] = 0x81+mode;//82 is binary. 81 is text. 
	if(nlen<126){
		buf[1]=nlen;
	}else if(nlen>UINT16_MAX){
		buf[1]=127;
		*((uint64_t*)&buf[2])=swap8bytes((long)nlen);
		offset+=8;//8 byte length
	}else{
		buf[1]=126;
		*((uint16_t*)&buf[2])=swap2bytes(nlen);
		offset+=2;//2 byte length
	}
	int fd_dest=(int)(long)userdata;
	//read tcp and send to websockets
	if(stread(fd, buf+offset, nlen) || stwrite(fd_dest, buf, nlen+offset)){
		return -1;
	}
	return 0;
}
static void ws_write_close(int fd){//send a ws close frame.
	char buf[4];
	buf[0]=0x88;//FIN=1; opcode=8
	buf[1]=0x2;//payload=2
	buf[2]=0x3;	buf[3]=0xE8;//1000
	if(stwrite(fd, buf, sizeof(buf))){
		dbg_time("Send close failed.\n");
	}
}
static int ws_fail(int fd, const char *format, ...) CHECK_ARG(2);
static int ws_fail(int fd, const char *format, ...) {
	format2fn;
	info("ws_received failed: %s", fn);
	ws_proxy_remove((void*)(long)fd, 1);
	return -1;
}
/*
	Receive websocket messages and call ws_proxy_command to handle it.
*/
static int ws_receive(struct pollfd *pfd, int flag){
	int fd=pfd->fd;
	if(flag==-1){//server ask client to close
		dbg_time("listen_sock requests shutdown\n");
		ws_write_close(fd);
		shutdown(fd, SHUT_WR);//make sure we don't reply to the close message.
	}
	static char *prev=NULL;//save previous message to handle continuation
	static size_t nprev=0;//length of previous message
	static int oprev=0;//opcode of previous message
    char buf0[BUF_SIZE];
	char *buf=buf0;
	const int len = recv(fd, buf, BUF_SIZE-1, 0);
	if (len<6){//server received msg must be at least 6 bytes
		return ws_fail(fd, "len=%d\n", len);
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
		buf=malloc(msg_len+offset+1);
		memcpy(buf, buf0, BUF_SIZE-1);
		if(stread(fd, buf+BUF_SIZE-1, msg_len+offset-(BUF_SIZE-1))){
			return ws_fail(fd, "stread");
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
		prev=realloc(prev, (nprev+msg_len+1));//save data
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
		//ssize_t len2 = recv(fd, buf+len, sizeof(buf)-len, 0);
	}else if(opcode==0x8){//client replied//request close
		dbg_time("client requests close. fd=%d\n", fd);
		ws_write_close(fd);
		return -1;//we need to close it first.
	}else if(opcode==0x9){//client sent a ping. we reply with pong
		buf[offset-2]=0x8A;//pong
		buf[offset-1]=msg_len;//cannot be more than 126
		if(stwrite(fd, buf+offset-2, msg_len+2)){
			return ws_fail(fd, "pong\n");
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
	char *key_start = strstr(buf, name);
	if (!key_start) return 0;
	key_start += strlen(name);
	while (*key_start == ' ') key_start++;
	char *key_end = strchr(key_start, '\r');
	if (!key_end) return 0;
	return strndup(key_start, key_end-key_start);
}
/*
	Upgrade http connection to websocket. it replaces the handler in listen_port.
*/
void http_upgrade_websocket(int fd, char *buf) {
	//Sec-WebSocket-Protocol: maos-monitor-protocol
	//dbg2_time("Upgrading %d to websockets\n", fd);
	char *key=http_parse_key(buf, "Sec-WebSocket-Key:");
	char *protocol=http_parse_key(buf, "Sec-WebSocket-Protocol:");
    unsigned char sha[20];
    char accept_src[128];
    snprintf(accept_src, sizeof(accept_src), "%s258EAFA5-E914-47DA-95CA-C5AB0DC85B11", key);
    SHA1_CTX ctx;
    sha1_init(&ctx);
    sha1_update(&ctx, (unsigned char*)accept_src, strlen(accept_src));
    sha1_final(&ctx, sha);
    char accept_key[64];
    base64_encode(sha, 20, accept_key);
	
    char resp[256];
    snprintf(resp, sizeof(resp),
             "HTTP/1.1 101 Switching Protocols\r\n"
             "Upgrade: websocket\r\n"
             "Connection: Upgrade\r\n"
             "Sec-WebSocket-Accept: %s\r\n"
			 "Sec-WebSocket-Protocol: %s\r\n\r\n",
             accept_key, protocol?protocol:"chat");
    stwrite(fd, resp, strlen(resp));
	listen_port_add(fd, POLLIN, ws_receive, "ws_receive");//switch to websocket handler
	if(!strcmp(protocol, "maos-monitor-protocol")){
		char cmd[11]="1&MONITOR;";
		ws_proxy_command(cmd, strlen(cmd), (ws_proxy_t){.send=ws_send, .userdata=(void*)(long)(fd)});
	}
	
	free(protocol);
	free(key);
}
/*
	Send a close request to client which will indicate closure. 
*/
void http_close(int client_fd){
	const char *resp =
		"HTTP/1.1 200 OK\r\n"
		"Content-Type: text/plain\r\n"
		"Content-Length: 12\r\n"
		"Connection: close\r\n"
		"\r\n"
		"Hello world\n";

	stwrite(client_fd, resp, strlen(resp));
	shutdown(client_fd, SHUT_WR);  // optional graceful signal
}
/*
	HTTP handler callable by listen_port. 
	It supports sending plain files (folder is disabled for security reasons) and websocket upgrade
*/
int http_handler(struct pollfd *pfd, int flag) {
	int fd=pfd->fd;
	if(flag==-1){
		http_close(fd);
		return 0;//wait for client to initiate close
	}
    char buf[BUF_SIZE];
    int n = recv(fd, buf, sizeof(buf) - 1, 0);
    if (n <= 0) return -1;
    buf[n] = 0;
	//info_time("received '%s'\n", buf);
    if (strncmp(buf, "GET ", 4) != 0){
		warning_time("buf shall start with 'GET': '%s'. close connection. fd=%d\n", buf, fd);
		return -1;//error
	}
    if (strstr(buf, "Upgrade: websocket")) {//handle upgrade
        http_upgrade_websocket(fd, buf);
		return 0;
    } else {//handle http send file
		int ans=0;
		char *path = buf + 4;
		char *sp = strchr(path, ' ');
		char* close=strstr(buf, "Connection: close")?"Connection: close\r\n":"";
		if (!sp) return -1;
		*sp = 0;
		if (strcmp(path, "/") == 0){
        	path = "/monitor.html";
		}
		if(path[0]=='/'){
			path++;//ignore leading /
		}
		char fullpath[256];
		if(strchr(path, '/')){//do not allow directory. return 404.
			//warning("directory is not supported in path: '%s'", path);
			path=NULL;
		}else{
        	snprintf(fullpath, sizeof(fullpath), "%s/%s", SRCDIR "/tools", path);	
		}
		//open File
		FILE *f = path?fopen(fullpath, "rb"):NULL;
		char header[256];

		if (!f) {//send 404
			snprintf(header, sizeof(header), "HTTP/1.1 404 Not Found\r\n%sContent-Length: 9\r\n\r\nNot Found", close);
			stwrite(fd, header, strlen(header));
		}else{//send file
			fseek(f, 0, SEEK_END);
			long size = ftell(f);
			rewind(f);

			snprintf(header, sizeof(header), "HTTP/1.1 200 OK\r\n%sContent-Length: %ld\r\n\r\n", close, size);
			if(stwrite(fd, header, strlen(header))){
				ans=-1;
			}else{
				char buf2[BUF_SIZE];
				size_t n2;
				
				while ((n2 = fread(buf2, 1, sizeof(buf2), f)) > 0 && !ans){
					ans=stwrite(fd, buf2, n2);
				}
			}
			fclose(f);
		}
		if(close[0]!=0){//close connection
			ans=-1;
		}
		return ans;
		
    }
}
