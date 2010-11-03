/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <execinfo.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <utime.h>
#include <fcntl.h>           /* For O_* constants */
#include <limits.h>
#include <errno.h>
#include <stdarg.h>
#include <dirent.h>
#include "common.h"
#include "misc.h"
#include "path.h"
#include "common.h"

/**
   Obtain the hostname of current machine. The returnned string
   should not be freed.
*/
__thread char fnglobal[PATH_MAX];//thread local storage
char* FF(const char *format, ...){
    va_list ap;
    va_start(ap, format);
    vsnprintf(fnglobal, sizeof(fnglobal), format, ap);
    va_end(ap);
    if(strlen(fnglobal)==0) 
	return NULL;
    else
	return fnglobal;
}



char *mybasename(const char *fn){
    /**
       Obtain the basename of a file. The returnned string must be freed.
    */
    char fn2[strlen(fn)+1];
    strcpy(fn2,fn);
    if(fn2[strlen(fn2)-1]=='/')
	fn2[strlen(fn2)-1]='\0';
    char* sep=strrchr(fn2,'/');
    if(!sep) sep=fn2; else sep++;
    char *bn=malloc(strlen(sep)+1);
    strcpy(bn,sep);
    return bn;
}
/**
   Copy a file from file stream src to dest.
*/
static void copyfile_fp(FILE *dest, FILE *src){
    char buffer[4096];
    size_t br,bw;
    while(!feof(src)){
	br=fread(buffer, 1, 4096, src);
	if((bw=fwrite(buffer, 1, br, dest))!=br){
	    error("copyfile: Write failed %ld of %ld written.\n", 
		  (long)bw,(long)br);
	}
    }
}
/**
   Copy a file from src to dest
*/
void copyfile(const char *dest, const char *src){
    FILE *psrc=fopen(src,"rb");
    if(!psrc){
	error("Open source file failed\n");
    }
    FILE *pdest=fopen(dest,"wb");
    if(!pdest){
	error("Open destination file failed\n");
    }
    copyfile_fp(pdest, psrc);
}


int check_suffix(const char *fn, const char *suffix){
    /**
       Check the suffix of a file.
    */
    if(!fn) return 0;
    if(strncmp(fn+strlen(fn)-strlen(suffix),suffix,strlen(suffix))){
	return 0;
    }else{
	return 1;
    }
}

char *argv2str(int argc, char **argv){
    char *cwd=mygetcwd();
    int slen=strlen(cwd)+2;
    for(int iarg=0; iarg<argc; iarg++){
	slen+=1+strlen(argv[iarg]);
    }
    char *scmd=calloc(slen, sizeof(char));
    char *home=getenv("HOME");//don't free
    if(!strncmp(cwd,home,strlen(home))){
	strcpy(scmd,"~");
	strcat(scmd,cwd+strlen(home));
    }else{
	strcpy(scmd,cwd);
    }
    strcat(scmd,"/");
    for(int iarg=0; iarg<argc; iarg++){
	strcat(scmd,argv[iarg]);
	strcat(scmd," ");
    }
    if(strlen(scmd)>slen-1) error("Overflow\n");
    free(cwd);
    return scmd;
}
void print_file(const char *fnin){
    char *fn=search_file(fnin);
    if(!fn){
	warning("%s not found\n", fnin);
	return;
    }
    FILE *fp;
    if(!(fp=fopen(fn,"r"))){
	error("Open %s failed\n",fn);
    }
    copyfile_fp(stderr, fp);
    /*const int nmax=128;
    char line[nmax];
    while(fgets(line, nmax, fp)){
	fprintf(stderr,"%s\n",line);
	}*/
    fclose(fp);
    free(fn);
}
