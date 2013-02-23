/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   \file sockio.c

   Routines handle socket i/o.
*/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include "misc.h"
#include "common.h"
#include "sockio.h"
static __attribute__((constructor))void init(){
    signal(SIGPIPE, SIG_IGN);
}
int stwrite(int sfd, const void *p, size_t len){
    int nwrite;
    const char *start=p;
    do{
	errno=0;
	nwrite=write(sfd, start, len);
	if(nwrite<=0){
	    warning("Write socket %d failed with error %d: %s\n",
		    sfd, errno, strerror(errno));
	    return -1;
	}else if(errno){
	    warning("stwrite: %s\n", strerror(errno));
	}
	len-=nwrite;
	start+=nwrite;
    }while(len>0);
    return 0;
}
int stread(int sfd, void *p, size_t len){
    int nread;
    char *start=p;
    do{
	errno=0;
	nread=read(sfd, start, len);
	if(nread<=0){
	    if(nread<0){
		warning3("Read socket %d failed with error %d: %s\n",sfd, errno, strerror(errno));
	    }else{
		warning3("stread: socket %d closed\n", sfd);
	    }
	    return -1;
	}else if(errno){
	    warning("stread: %s\n", strerror(errno));
	}
	len-=nread;
	start+=nread;
    }while (len>0);
    return 0;
}
/*Prevent calling read/write in this file from now on. use stread/stwrite instead */
#define read READ_IS_PROHIBITED
#define write WRITE_IS_PROHIBITED
/**
   Write a string to socket
*/
int stwritestr(int sfd, const char *str){
    if(str){
	int len=strlen(str)+1;
	return stwriteint(sfd, len) || stwrite(sfd, str, len);
    }else{
	return stwriteint(sfd, 0);
    }
}
/**
   Read a string from socket
*/
int streadstr(int sfd, char **str){
    int len;
    int ans=streadint(sfd, &len);
    if(!ans && len>0){
	*str=calloc(1, sizeof(char)*len);
	ans=stread(sfd, *str, len);
	if(ans){
	    free(*str);
	    *str=NULL;
	}
    }else{
	*str=NULL;
    }
    return ans;
}
/**
   Write a string array to socket
 */
int stwritestrarr(int sfd, const char *const *str, int nstr){
    int ans=stwriteint(sfd, nstr);
    for(int i=0; i<nstr && !ans; i++){
	ans=stwritestr(sfd, str[i]);
    }
    return ans;
}
/**
   Read a string array from socket
*/
int streadstrarr(int sfd, char ***str, int *nstr){
    int ans=streadint(sfd, nstr);
    if(ans) return ans;
    *str=calloc(*nstr, sizeof(char*));
    for(int istr=0; istr<*nstr && !ans; istr++){
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
    struct cmsghdr *cmsg=CMSG_FIRSTHDR(&msg);
    cmsg->cmsg_level = SOL_SOCKET;
    cmsg->cmsg_type = SCM_RIGHTS;
    cmsg->cmsg_len = CMSG_LEN(sizeof(int));
    *(int*)CMSG_DATA(cmsg)=fd;
    if(sendmsg(sfd, &msg, 0)<0){
	warning("sendmsg failed to send fd %d over %d\n", fd, sfd);
	ans=-1;
    }else{
	ans=0;
    }
    return ans;
}
/**
   Read a file descriptor from socket. Use with stwritefd();
*/
int streadfd(int sfd, int *fd){
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
    struct cmsghdr *cmsg=CMSG_FIRSTHDR(&msg);
    if(recvmsg(sfd, &msg, 0)<0){
	warning("recvmsg failed to receive fd over %d\n", sfd);
	ans=-1;
    }else{
	memmove(fd, CMSG_DATA(cmsg), sizeof(int));
	warning("recvmsg received %d. ans=%d\n", *fd, ans);
	ans=0;
    }
    return ans;
}
