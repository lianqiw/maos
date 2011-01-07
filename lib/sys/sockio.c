/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <unistd.h>
#include <string.h>
#include "misc.h"
#include "common.h"
#include "sockio.h"

/**
   Socket IO
*/
int swrite(int *fd, const void *p, size_t len){
    if(*fd==-1){
	warning3("Trying to write to a closed socket\n");
	return -1;
    }
    if(write(*fd, p, len)!=len){
	perror("swrite");
	warning3("Write failed. close socket %d\n", *fd);
	close(*fd);
	*fd=-1;
	return -1;
    }else{
	return 0;
    }
}
int sread(int *fd, void *p, size_t len){
    if(*fd==-1){
	warning3("Trying to read a closed socket\n");
	return -1;
    }
    if(read(*fd, p, len)!=len){
	perror("sread");
	warning3("Read failed. close socket %d\n",*fd);
	close(*fd);
	*fd=-1;
	return -1;
    }else{
	return 0;
    }
}
//Provent calling read/write in this file from now on. use sread/swrite instead
#define read READ_IS_PROHIBITED
#define write WRITE_IS_PROHIBITED
void swriteint(int *fd, int cmd){
    swrite(fd, &cmd, sizeof(int));
}
void swriteintarr(int *fd, int* cmd, unsigned int len){
    swrite(fd,cmd,len*sizeof(int));
}
int sreadint(int *fd){
    int cmd=-1;
    sread(fd, &cmd, sizeof(int));
    return cmd;
}
void sreadintarr(int *fd, int* cmd, unsigned int len){
    sread(fd,cmd,len*sizeof(int));
}
void swritestr(int *fd, const char *str){
    if(str){
	int len=strlen(str)+1;
	swriteint(fd, len);
	swrite(fd, str, len);
    }else{
	swriteint(fd, 0);
    }
}

char *sreadstr(int *fd){
    int len;
    char *str;
    len=sreadint(fd);
    if(len>0){
	str=calloc(1, sizeof(char)*len);
	if(sread(fd, str, len)){
	    free(str);
	    str=NULL;
	}
    }else{
	str=NULL;
    }
    return str;
}
