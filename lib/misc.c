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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <utime.h>
#include <fcntl.h>           /* For O_* constants */
#include <limits.h>
#include <errno.h>
#include <stdarg.h>
#include <dirent.h>
#include <time.h>
#include "common.h"
#include "thread.h"
#include "process.h"
#include "misc.h"
#include "path.h"

/**
   Obtain the basename of a file. The returnned string must be freed.
*/
char *mybasename(const char *fn){
    if(!fn || strlen(fn)==0) return NULL;
    char fn2[PATH_MAX];
    strcpy(fn2,fn);
    //If this is a folder, remove the last /
    if(fn2[strlen(fn2)-1]=='/')
	fn2[strlen(fn2)-1]='\0';
    char* sep=strrchr(fn2,'/');
    if(!sep){
	sep=fn2; 
    }else{
	sep++;
    }
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


/**
   Check the suffix of a file.
*/
int check_suffix(const char *fn, const char *suffix){
    if(!fn || !suffix) return 0;
    int lfn=strlen(fn);
    int lsu=strlen(suffix);
    if(lfn < lsu) return 0;
    if(mystrcmp(fn+lfn-lsu,suffix)){
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
    if(!mystrcmp(cwd,HOME)){
	strcpy(scmd,"~");
	strcat(scmd,cwd+strlen(HOME));
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
    fclose(fp);
    free(fn);
}

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
const char *myasctime(void){
    static char st[64];
    time_t a;
    time(&a);
    ctime_r(&a, st);
    st[strlen(st)-1]='\0';//remove final \n
    return st;
}
/**
   Get furrent time in ascii string that doesn't contain
   spaces. Suitable for use in filenames. The returnned string
   must be freed. */
char *strtime(void){
    char str[64];
    time_t t=myclocki();
    struct tm tmp;
    localtime_r(&t,&tmp);//don't free tmp
    strftime(str,64,"%F-%H%M%S",&tmp);
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
    char cwd0[PATH_MAX];
    if(!getcwd(cwd0,PATH_MAX)) 
	error("Error getting current directory\n");
    return strdup(cwd0);
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
/**
   Test whether a file exists.
*/
int exist(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf);
}

/**
   Test whether fn is directory
*/
int isdir(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISDIR(buf.st_mode);
}

/**
   Test whether fn is ordinary file
*/
int isfile(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISREG(buf.st_mode);
}

/**
   Test whether fn is a symbolic link
*/
int islink(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISLNK(buf.st_mode);
}
/**
 * Compute length of file in Bytes
 */
off_t flen(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return stat(fn, &buf)?0:buf.st_size;
}
/**
   Return the modification time of the file
 */
time_t fmtime(const char *fn){
    struct stat buf;
    if(!fn || stat(fn, &buf)) return 0;
    return buf.st_ctime;
}
void touch(const char *fn){
    /**
       Update a file's modtime to current.
    */
    if(utimes(fn, NULL)){
	if(errno==ENOENT){
	    FILE *fp=fopen(fn, "w");
	    fclose(fp);
	}
    }
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
	    out=stradd(HOME,fn+1,NULL);
	}else{
	    out=strdup(fn);
	}
	*fnout=out;
    }else{
	char *out=*fnout;
	if(fn[0]=='~'){
	    strcpy(out,HOME);
	    strcat(out,fn+1);
	}else{
	    strcpy(out,fn);
	}
    }
}
char *mystrndup(const char *A, int len){
    int len2=strlen(A);
    if(len2<len) len=len2;
    char *B=malloc(len+1);
    memcpy(B,A,len);
    B[len]='\0';
    return B;
}

#if USE_MEM == 1
#undef strdup
/**
   declare strdup so my memory mangement mem.c is happy when DEBUG=1. Handles
NULL pointer correctly.  */
char *mystrdup(const char *A){
    if(!A){
	return NULL;
    }else{
	int nlen=strlen(A);
	char *B=malloc(nlen+1);
	memcpy(B,A,nlen+1);
	return B;
    }
}
#endif
char* (*strdup0)(const char *)=strdup;

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
	if(!stat(fnfull,&buf) && S_ISREG(buf.st_mode) && buf.st_mtime<sec2){
	    remove(fnfull);
	    info2("Remove %s. %ld days old\n", fnfull, (myclocki()-buf.st_mtime)/3600/24);
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
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.
 */
int mystrcmp(const char *a, const char *b){
    if(!a || !b) return 1;
    int la=strlen(a);
    int lb=strlen(b);
    if(la==0 || lb==0 || la<lb){
	return 1;
    }else{
	return strncmp(a,b,lb);
    }
}


/**
   Compute the factorial. Overflow LONG if n>20, so we use double as output.*/
double factorial(long n){
    double fact=1;
    while(n>1){
	fact*=n--;
    }
    if(isinf(fact)){
	error("Factorial overflows\n");
    }
    return fact;
}
/**
   Make the fd close on exec.
*/
void cloexec(int fd){
    int oldflag=fcntl(fd, F_GETFD, 0);
    if(oldflag!=-1){
	oldflag |= FD_CLOEXEC;
	fcntl(fd, F_SETFD, oldflag);
    }
}
/**
   wrap of nanosleep
*/
void mysleep(double sec){
    struct timespec ts;
    ts.tv_sec=(time_t)trunc(sec);
    ts.tv_nsec=(long)((sec-ts.tv_sec)*1e9);
    nanosleep(&ts, NULL);
}

/**
   Pause
*/
void mypause(void){
    info2("Press ENTER key to continue.\n");
    int ans;
    while((ans=getchar())!='\n');
}
