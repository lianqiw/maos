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
#include <sys/stat.h>//added in apple.
#include <utime.h>
#include <fcntl.h>           /* For O_* constants */
#include <limits.h>
#include <errno.h>
#include <dirent.h>
#include <time.h>
#include <signal.h>
#include <stdarg.h>
#include "common.h"
#include "thread.h"
#include "misc.h"
/**
   Get current time in seconds.
*/
int myclocki(){
    time_t a;
    return time(&a);
}
/**
   Get current time in ascii string for easy print. The string
   contains spaces and is not suitable to use in filename. The
   returned string should not be modified.  */
char *myasctime(void){
    time_t a;
    time(&a);
    char *st=ctime(&a);
    st[strlen(st)-1]='\0';//remove final \n
    return st;
}
/**
   Get furrent time in ascii string that doesn't contain
   spaces. Suitable for use in filenames. The returnned string
   should be freed. */
char *strtime(void){
    char str[64];
    time_t t=myclocki();
    struct tm *tmp=localtime(&t);//don't free tmp
    strftime(str,64,"%F-%H%M%S",tmp);
    char *dir=strdup(str);
    return dir;
}
const char *myhostname(void){
    static int inited=0;
    static char host[60];
    PNEW(lock);
    LOCK(lock);
    if(!inited){
	if(gethostname(host,60)){
	    warning("Unable to get hostname\n");
	    sprintf(host,"localhost");
	}
	inited=1;
    }
    UNLOCK(lock);
    return host;
}

/**
   Get current time in milli-second resolution.
*/
double myclockd(void){
    struct timeval tk;
    gettimeofday(&tk,NULL);
    return (double)tk.tv_sec+(double)tk.tv_usec*1e-6;
}/**
   Get current directory. The returnned string must be freed.
*/
char *mygetcwd(void){
    char *cwd0=malloc(sizeof(char)*PATH_MAX);
    if(!getcwd(cwd0,PATH_MAX)) 
	error("Error getting current directory\n");
    cwd0=realloc(cwd0,strlen(cwd0)+1);
    if(!exist(cwd0)) error("Failed to get current directory\n");
    return cwd0;
}
/**
   Translate a path into absolute path.
*/
char *myabspath(const char *path){
#if _BSD_SOURCE || _XOPEN_SOURCE >= 500
    return realpath(path, NULL);
#else
    char *cpath=mygetcwd();
    if(chdir(path)){
	error("path %s doesn't exist\n",path);
    }
    char *abspath=mygetcwd();
    if(chdir(cpath)){
	error("Unable to cd back to %s\n",cpath);
    }
    free(cpath);
    return abspath;
#endif
}

void mysymlink(const char *fn, const char *fnlink){
    if(!exist(fn)) return;
    remove(fnlink);
    if(symlink(fn, fnlink)){
	warning("Unable to make symlink %s->%s\n",fnlink,fn);
    }
}
int exist(const char *fn){
    /**
       Test whether a file exists.
    */
    struct stat buf;
    return !stat(fn, &buf);
}
void touch(const char *fn){
    /**
       Update a file's modtime to current.
    */
    utimes(fn, NULL);
}

/**
   Concatenate many strings. Argument list must end with NULL.
*/
char *stradd(const char* a, ...){
    char *out;
    int n;
    va_list ap;
    va_start(ap,a);
    n=strlen(a)+1;
    while(1){
	const char *arg=va_arg(ap,const char*);
	if(!arg)
	    break;
	n+=strlen(arg);
    }
    va_end(ap);
    out=malloc(n*sizeof(char));
    strcpy(out,a);
    va_start(ap,a);
    while(1){
	const char *arg=va_arg(ap,const char*);
	if(!arg)
	    break;
	strcat(out,arg);
    }
    va_end(ap);
    return out;
}

/**
   translate a filename into absolute file name that starts with /
*/
void expand_filename(char **fnout, const char *fn){
    //if fn contains leading ~, expand it.
    if(!(*fnout)){
	char *out;
	if(fn[0]=='~'){
	    char *HOME=getenv("HOME");
	    out=stradd(HOME,fn+1,NULL);
	}else{
	    out=mystrdup(fn);
	}
	*fnout=out;
    }else{
	char *out=*fnout;
	if(fn[0]=='~'){
	    char *HOME=getenv("HOME");
	    strcpy(out,HOME);
	    strcat(out,fn+1);
	}else{
	    strcpy(out,fn);
	}
    }
}

#undef strdup
#undef strndup

char *mystrndup(const char *A, int len){
    /**
       declare strndup so my memory mangement mem.c is happy when DEBUG=1
    */
    if(!A) return NULL;
    int ncpy=strlen(A);
    if(ncpy>len) ncpy=len;
    char *B=malloc(sizeof(char)*(ncpy+1));
    memcpy(B,A,sizeof(char)*ncpy);
    B[ncpy]='\0';
    return B;
}
char *mystrdup(const char *A){
    /**
       declare strdup so my memory mangement mem.c is happy when DEBUG=1
    */
    if(!A) return NULL;
    return mystrndup(A, strlen(A));
}
/**
 Compute max, min and sum of a double vector*/
void maxmindbl(const double *restrict p, long N, 
	       double *restrict max, double *restrict min, double *restrict sum){
    double a,b,s;
    long i;
    a=p[0]; b=p[0];s=0;
    for(i=1; i<N; i++){
	s+=p[i];
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    if(max) *max=a;
    if(min) *min=b;
    if(sum) *sum=s;
}
/**
 Compute max, min and sum of a long vector*/
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min, long *restrict sum){
    long a,b,s;
    long i;
    a=p[0]; b=p[0]; 
    s=0;
    for(i=0; i<N; i++){
	s+=p[i];
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    if(max)*max=a; 
    if(min)*min=b; 
    if(sum)*sum=s;
}
/**
 Compute max, min and sum of a complex vector*/
void maxmincmp(const dcomplex *restrict p, long N,
	       double *restrict max, double *restrict min, double *restrict sum){
    double a,b,s;
    long i;
    a=cabs(p[0]); 
    b=cabs(p[0]);
    s=0;
    for(i=0; i<N; i++){
	double tmp=cabs(p[i]);
	s+=tmp;
	if(tmp>a) a=tmp;
	if(tmp<b) b=tmp;
    }
    if(max)*max=a; 
    if(min)*min=b; 
    if(sum)*sum=s;
}
void remove_file_older(const char *fndir, long sec){
    /**
       Remove files that are older than sec seconds in folder fndir.
     */
    DIR *dir=opendir(fndir);
    if(!dir){
	error("Unable to open directory %s\n",fndir);
    }
    struct dirent *dp;
    struct stat buf;
    char fnfull[PATH_MAX];
    long sec2=myclocki()-sec;
    while((dp=readdir(dir))){
	snprintf(fnfull,PATH_MAX,"%s/%s",fndir,dp->d_name);
	if(stat(fnfull,&buf)){
	    perror("stat");
	    warning("Unable to stat %s\n",fnfull);
	}else if(S_ISREG(buf.st_mode)){
	    if(buf.st_mtime<sec2){
		remove(fnfull);
		info2("Remove %s. %ld seconds old\n", fnfull, myclocki()-buf.st_mtime);
	    }
	}
    }
    closedir(dir);
}
void mymkdir(const char *format, ...){
    /**
       Make dirs recursively. like mkdir -p in bash
    */
    format2fn;
    if(!fn) return;
    if(fn[strlen(fn)-1]=='/')
	fn[strlen(fn)-1]='/';
    if(mkdir(fn, 0700)==-1){
	if(errno==EEXIST){
	    return;
	}else if(errno==ENOENT){
	    char *tmp=rindex(fn,'/');
	    if(!tmp){
		error("Unable to mkdir '%s'\n",fn);
	    }
	    tmp[0]='\0';
	    mymkdir("%s",fn);
	    tmp[0]='/';
	    if(mkdir(fn,0700)==-1&&errno!=EEXIST){
		error("Unable to mkdir '%s'\n",fn);
	    }
	}
    }
}
