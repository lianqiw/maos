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
#include <sys/types.h>
#include <signal.h>
#ifdef __linux__
#include <linux/limits.h>
#else
#include <limits.h>
#endif
#include "common.h"
#include "misc.h"
#include "io.h"

void writeintsock(int fd, int cmd){
    if(write(fd, &cmd, sizeof(int))!=sizeof(int))
	warning3("Write failed\n");
}
int readint(int fd){
    int cmd;
    if(read(fd, &cmd, sizeof(int))!=sizeof(int)){
	warning3("Read failed\n");
	cmd=-1;
    }
    return cmd;
}
void writestr(int fd, const char *str){
    if(str){
	int len=strlen(str)+1;
	writeintsock(fd, len);
	if(write(fd, str, len)!=len)
	    warning3("Write failed\n");
    }else{
	writeintsock(fd, 0);
    }
}
char *readstr(int fd){
    int len;
    char *str;
    len=readint(fd);
    if(len){
	str=calloc(1, sizeof(char)*len);
	if(read(fd, str, len)!=len){
	    warning3("Read failed\n");
	    free(str);
	    str=NULL;
	}
    }else{
	str=NULL;
    }
    return str;
}
