/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include <errno.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include <poll.h>
#include "misc.h"
#include "common.h"
#include "sockio.h"
#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif
int stwrite(int sfd, const void* p, size_t len){
	int is_sock=issock(sfd);
	ssize_t nwrite;
	long left=(long)len;
	do{
		if(is_sock){
			nwrite=send(sfd, p, len, MSG_NOSIGNAL);
		} else{
			nwrite=write(sfd, p, left);
		}
		p=(char*)p+nwrite; left-=nwrite;
	} while(nwrite>0&&left>0);
	return (nwrite<0||left)?-1:0;//-1 indicate error/closed
}
int stread(int sfd, void* p, size_t len){
	int is_sock=issock(sfd);
	ssize_t nread;
	long left=len;
	do{
		if(is_sock){
			nread=recv(sfd, p, left, 0);//MSG_WAITALL);
		} else{
			nread=read(sfd, p, left);
		}
		p=(char*)p+nread; left-=nread;
	} while(nread>0&&left>0);
	//if(left) info("stread failed: nread=%ld, left=%ld\n", nread, left);
	return (nread<0||(nread==0&&len!=0)||left)>0?-1:0; //-1 indicated error/closed
}
/*Write long messages with smaller buffer*/
int stwrite2(int sfd, const void* p, size_t len, size_t nbuf){
	if(nbuf>len) nbuf=len;
	ssize_t nwrite;
	long left=len;//do not use size_t which is unsigned
#ifdef __linux__
	int is_sock=issock(sfd);
#endif
	do{
#ifdef __linux__
		if(is_sock){
			nwrite=send(sfd, p, len, MSG_NOSIGNAL);
		} else
#endif
			nwrite=write(sfd, p, nbuf);

		p=(char*)p+nwrite; left-=nwrite;
	} while(nwrite>0&&left>0);
	return (nwrite<0||left)?-1:0;
}
/*Read long messages with smaller buffer*/
int stread2(int sfd, void* p, size_t len, size_t nbuf){
	if(nbuf>len) nbuf=len;
	ssize_t nread;
	long left=len;
	do{
		nread=read(sfd, p, nbuf);
		p=(char*)p+nread; left-=nread;
	} while(nread>0&&left>0);
	return (nread<0||(nread==0&&len!=0)||left)?-1:0;
}
/**
   Write a string to socket
*/
int stwritestr(int sfd, const char* str){
	if(str){
		int len=strlen(str)+1;
		return stwriteint(sfd, len)||stwrite(sfd, str, len);
	} else{
		return stwriteint(sfd, 0);
	}
}
/**
   Read a string from socket
*/
int streadstr(int sfd, char** str){
	int len;
	int ans=streadint(sfd, &len);
	if(!ans&&len>0){
		*str=(char*)calloc(1, sizeof(char)*len);
		ans=stread(sfd, *str, len);
		if(ans){
			free(*str);
			*str=NULL;
		}
	} else{
		*str=NULL;
	}
	return ans;
}
/**
   Write a string array to socket
 */
int stwritestrarr(int sfd, const char* const* str, int nstr){
	int ans=stwriteint(sfd, nstr);
	for(int i=0; i<nstr&&!ans; i++){
		ans=stwritestr(sfd, str[i]);
	}
	return ans;
}
/**
   Read a string array from socket
*/
int streadstrarr(int sfd, char*** str, int* nstr){
	int ans=streadint(sfd, nstr);
	if(ans) return ans;
	*str=mycalloc(*nstr, char*);
	for(int istr=0; istr<*nstr&&!ans; istr++){
		ans=streadstr(sfd, &(*str)[istr]);
	}
	return ans;
}
/**
   Wrie a file descriptor to socket. Transfer a open file (socket) descriptor to
   another process in the same hosts. sfd has to be a AF_UNIX socket.
*/
int stwritefd(int sfd, int fd){
	int ans=0;
	char buf[1]={0};
	struct iovec iov;
	iov.iov_base=buf;
	iov.iov_len=1;
	struct msghdr msg={0};
	msg.msg_iov=&iov;
	msg.msg_iovlen=1;
	char cms[CMSG_SPACE(sizeof(int))];
	msg.msg_control=(caddr_t)cms;
	msg.msg_controllen=CMSG_LEN(sizeof(int));
	struct cmsghdr* cmsg=CMSG_FIRSTHDR(&msg);
	cmsg->cmsg_level=SOL_SOCKET;
	cmsg->cmsg_type=SCM_RIGHTS;
	cmsg->cmsg_len=CMSG_LEN(sizeof(int));
	memcpy(CMSG_DATA(cmsg), &fd, sizeof(int));
	if(sendmsg(sfd, &msg, 0)<0){
		warning("sendmsg failed to send fd %d over %d\n", fd, sfd);
		ans=-1;
	} else{
		ans=0;
	}
	return ans;
}
/**
   Read a file descriptor from socket. Use with stwritefd();
*/
int streadfd(int sfd, int* fd){
	int ans=0;
	char buf[1]={0};
	struct iovec iov;
	iov.iov_base=buf;
	iov.iov_len=1;
	struct msghdr msg={0};
	msg.msg_iov=&iov;
	msg.msg_iovlen=1;
	char cms[CMSG_SPACE(sizeof(int))];
	msg.msg_control=(caddr_t)cms;
	msg.msg_controllen=CMSG_LEN(sizeof(int));
	struct cmsghdr* cmsg=CMSG_FIRSTHDR(&msg);
	if(recvmsg(sfd, &msg, 0)<0){
		warning("recvmsg failed to receive fd over %d\n", sfd);
		ans=-1;
	} else{
		memmove(fd, CMSG_DATA(cmsg), sizeof(int));
		ans=0;
	}
	return ans;
}
/*
  Check whether the client has shutdown a connection (TCP)
  Returns 1 if closed;
*/

int stcheck(int sfd){
	int ans=0;
	char buf[4];
	ssize_t len=recv(sfd, buf, sizeof buf, MSG_DONTWAIT|MSG_PEEK);
	if(len==0){
		dbg_time("Client has orderly shutdown\n");
		ans=1;
	}else if(len==-1){
		dbg_time("recv failed with errno %d: %s\n", errno, strerror(errno));
	}
	/*
	//the following Implementation is not correct.
	  //poll does not detect closed socket used for writing.
	struct pollfd data={sfd, POLLIN|POLLPRI|POLLOUT, 0 };
	if(poll(&data, 1, 1)){
	dbg("data.revents=%d\n", data.revents);
	if(data.revents&POLLERR || data.revents & POLLHUP || data.revents & POLLNVAL){
		warning("socket %d is no longer valid\n", sfd);
		ans=1;
	}
	}*/
	return ans;
}
