/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "misc.h"
#include "common.h"
#include "sockio.h"
static __attribute__((constructor))void init(){
    signal(SIGPIPE, SIG_IGN);
}
int stwrite(int sfd, const void *p, size_t len){
    if(write(sfd, p, len)!=len){
	perror("stwrite");
	warning3("Write socket %d failed with error %d\n", sfd, errno);
	return -1;
    }else{
	return 0;
    }
}
int stread(int sfd, void *p, size_t len){
    int nread;
    char *start=p;
    do{
	nread=read(sfd, start, len);
	if(nread<=0){
	    perror("stread");
	    warning3("Read socket %d failed with error %d\n",sfd, errno);
	    return -1;
	}
	len-=nread;
	start+=nread;
    }while (len>0);
    return 0;
}
/*Prevent calling read/write in this file from now on. use stread/stwrite instead */
#define read READ_IS_PROHIBITED
#define write WRITE_IS_PROHIBITED
int stwriteint(int sfd, int cmd){
    return stwrite(sfd, &cmd, sizeof(int));
}
int stwriteintarr(int sfd, int* cmd, unsigned int len){
    return stwrite(sfd,cmd,len*sizeof(int));
}
int streadint(int sfd, int *cmd){
    return stread(sfd, cmd, sizeof(int));
}
int streadintarr(int sfd, int* cmd, unsigned int len){
    return stread(sfd,cmd,len*sizeof(int));
}
int stwritestr(int sfd, const char *str){
    if(str){
	int len=strlen(str)+1;
	return stwriteint(sfd, len) || stwrite(sfd, str, len);
    }else{
	return stwriteint(sfd, 0);
    }
}

int streadstr(int sfd, char **str){
    int len;
    int err=streadint(sfd, &len);
    if(!err && len>0){
	*str=calloc(1, sizeof(char)*len);
	err=stread(sfd, *str, len);
	if(err){
	    free(*str);
	    *str=NULL;
	}
    }else{
	*str=NULL;
    }
    return err;
}
