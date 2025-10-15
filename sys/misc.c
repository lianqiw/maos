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
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h> //fstat
#include <sys/statvfs.h>
#include <sys/resource.h>
#include <utime.h>
#include <fcntl.h>           /* For O_* constants */
#include <limits.h>
#include <errno.h>
#include <dirent.h>
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include "common.h"
#include "thread.h"
#include "process.h"
#include "misc.h"
#include "path.h"
#include "bin.h"
#include "scheduler_client.h"
/**
   Obtain the dirname of a path. See mybasename().
*/
char* mydirname(const char* fn){
	if(!fn||strlen(fn)==0) return NULL;
	char fn2[PATH_MAX];
	strncpy(fn2, fn, PATH_MAX-1);
	/*If this is a folder, remove the last / */
	if(fn2[strlen(fn2)-1]=='/')
		fn2[strlen(fn2)-1]='\0';
	char* sep=strrchr(fn2, '/');
	if(!sep){
		fn2[0]='.'; fn2[1]='\0';
	} else{
		sep[0]='\0';
	}
	return mystrdup(fn2);
}
/**
   Copy a file from file stream src to dest.
*/
static int copyfile_fp(FILE* src, FILE* dest){
	int ans=0;
	char buffer[4096];
	size_t br, bw;
	while(!feof(src)){
		br=fread(buffer, 1, 4096, src);
		if((bw=fwrite(buffer, 1, br, dest))!=br){
			warning("copyfile: Write failed %ld of %ld written: %s\n",	(long)bw, (long)br, strerror(errno));
			ans=-1;
		}
	}
	return ans;
}
/**
   Copy a file from src to dest
*/
int copyfile(const char* src, const char* dest){
	int ans=0;
	FILE* psrc=fopen(src, "rb");
	if(psrc){
		FILE* pdest=fopen(dest, "wb");
		if(pdest){
			ans=copyfile_fp(psrc, pdest);
			fclose(pdest);
		} else{
			warning("Open destination file failed: %s\n", strerror(errno));
			ans=1;
		}
		fclose(psrc);
	}else{
		warning("Open source file failed: %s\n", strerror(errno));
		ans=1;
	}
	return ans;
}


/**
   Check the suffix of a file.
*/
int check_suffix(const char* fn, const char* suffix){
	if(!fn||!suffix) return 0;
	size_t lfn=strlen(fn);
	size_t lsu=strlen(suffix);
	if(lfn<lsu) return 0;
	if(mystrcmp(fn+lfn-lsu, suffix)){
		return 0;
	} else{
		return 1;
	}
}
/**
   Convert argc, argv to a single string, prefixed by the current directory.
*/
char* argv2str(int argc, const char* argv[], const char* delim){
	if(!argc) return NULL;
	size_t slen=strlen(DIRSTART)+2;
	if(!delim) delim=" ";
	for(int iarg=0; iarg<argc; iarg++){
		slen+=strlen(delim)+strlen(argv[iarg]);
	}
	char* scmd=mycalloc(slen, char);
	
	const char* exename=mybasename(argv[0]);
	snprintf(scmd, slen, "%s/%s%s", DIRSTART, exename, delim);
	for(int iarg=1; iarg<argc; iarg++){
		if(argv[iarg]&&strlen(argv[iarg])>0){
			strcat(scmd, argv[iarg]);
			strcat(scmd, delim);
		}
	}
	if(strlen(scmd)>slen-1) error("Overflow\n");
	char* scmdend=scmd+strlen(scmd);
	if(delim[0]!='\n'){
		for(char* p=scmd; p<scmdend; p++){
			if(p[0]=='\n') p[0]=' ';
		}
	}
	while(scmdend>scmd+1&&scmdend[0]==delim[0]&&scmdend[-1]==delim[0]){
		scmdend[0]='\0';
		scmdend--;
	}
	return scmd;
}
/**
   Print the content of a file.
*/
void print_file(const char* fnin){
	char* fn=search_file(fnin, 1);
	if(!fn){
		warning("%s not found\n", fnin);
		return;
	}
	FILE* fp;
	if(!(fp=fopen(fn, "r"))){
		warning("Open %s failed\n", fn);
	}
	copyfile_fp(fp, stdout);
	fflush(stdout);
	fclose(fp);
	free(fn);
}

/**
   Get current time in seconds as an integer.
*/
time_t myclocki(){
	time_t a;
	return time(&a);
}
/**
   Get current time in ascii string for easy print. The string
   contains spaces and is not suitable to use in filename. The
   returned string should not be modified.  */
const char* myasctime(time_t at){
	static char st[20];
	if(!at) at=time(NULL);
	struct tm am;
	localtime_r(&at, &am);
	snprintf(st, sizeof(st), "%04d/%02d/%02d %02d:%02d:%02d",
		am.tm_year+1900, am.tm_mon+1, am.tm_mday, am.tm_hour, am.tm_min, am.tm_sec);
	return st;
}
/**
   Get furrent time in ascii string that doesn't contain
   spaces. Suitable for use in filenames. The returnned string
   must be freed. */
char* strtime_pid(void){
	char str[64];
	time_t t=myclocki();
	struct tm tmp;
	localtime_r(&t, &tmp);/*don't free tmp */
	strftime(str, 64, "%F-%H%M%S", &tmp);
	char pid[20];
	snprintf(pid, 20, "-%05d", (int)getpid());
	char* dir=stradd(str, pid, NULL);
	return dir;
}

/**
   Get current time in nano-second resolution.
*/
double myclockd(void){
	static time_t t0=0;//to avoid precision error when cast to float
#if defined(_POSIX_TIMERS) && _POSIX_TIMERS > 0
	struct timespec tk;
	clock_gettime(CLOCK_MONOTONIC, &tk);
	if(!t0) t0=tk.tv_sec;
	return (double)(tk.tv_sec-t0)+(double)tk.tv_nsec*1e-9;
#else
	struct timeval tk;
	gettimeofday(&tk, NULL);
	if(!t0) t0=tk.tv_sec;
	return (double)(tk.tv_sec-t0)+(double)tk.tv_usec*1e-6;
#endif
}
/**
   Get current directory. The returnned string must be freed.
*/
char* mygetcwd(void){
	char cwd0[PATH_MAX];
	if(!getcwd(cwd0, PATH_MAX)){
		dbg("Error getting current directory: %s\n", strerror(errno));
		strncpy(cwd0, getenv("PWD"), PATH_MAX); cwd0[PATH_MAX-1]=0;
	}
	return strdup(cwd0);
}
/**
Translate a path into absolute path. The caller shall free the returned string.
*/
char *myabspath(const char *path){
	if(!path) return mygetcwd();
	if(path[0]=='/') return strdup(path);
	char path2[PATH_MAX];
	if(path[0]=='~' && path[1]=='/'){
		snprintf(path2, PATH_MAX, "%s/%s", HOME, path+2);
	}else{
		if(path[0]=='.'){
			if(path[1]=='/'){
				path+=2;
			}else if(path[1]=='\0'){
				path+=1;
			}
		}
		if(!getcwd(path2, PATH_MAX)){
			strncpy(path2, getenv("PWD"), PATH_MAX); path2[PATH_MAX-1]=0;
		}
		while(!mystrcmp(path, "../")){
			char *tmp=strrchr(path2, '/');
			if(tmp){
				tmp[0]='\0';
				path+=3;
				while(path[0]=='/') path++;
			}else{
				error("Relative path `%s` has too many parent levels from `%s`\n", path, path2);
				break;
			}
		}
		if(path[0]!='\0'){
			strncat(path2, "/", PATH_MAX-strlen(path2)-1);
			strncat(path2, path, PATH_MAX-strlen(path2)-1);
		}
	}
	if(!exist(path2)){
		warning("File or directory %s does not exist\n", path2);
	}
	return strdup(path2);
}
/**
 * Create symbolic link
 * */
int mysymlink(const char* source, const char* dest){
	int ans=0;
	if(!exist(source)) ans=-1;
	if(!ans && exist(dest)) ans=remove(dest);
	if(!ans) ans=symlink(source, dest);
	if(ans)	warning("Unable to symlink %s to %s: %s\n", source, dest, strerror(errno));
	return ans;
}
/**
 * Create hard link
 * */
int mylink(const char* source, const char* dest){
	int ans=0;
	if(!exist(source)) ans=-1;
	if(!ans&&exist(dest)) ans=remove(dest);
	if(!ans) ans=link(source, dest);
	if(ans){
		dbg("Unable to link %s to %s, try copyfile instead: %s\n", source, dest, strerror(errno));
		ans=copyfile(source, dest);
	}
	return ans;
}
/**
   Test whether a file exists.
*/
int exist(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return !lstat(fn, &buf);
}
/*
	Update the modification time of a file and create it if not exist.
*/
void touch(const char *format, ...){
	format2fn;
	if(!exist(fn)){
		creat(fn, 0666);
	}else{
		utimes(fn, NULL);
	}
}
/**
   Test whether fn is directory
*/
int isdir(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return !stat(fn, &buf)&&S_ISDIR(buf.st_mode);
}

/**
   Test whether fn is ordinary file
*/
int isfile(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return !stat(fn, &buf)&&S_ISREG(buf.st_mode);
}

/**
   Test whether fn is a symbolic link
*/
int islink(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return !stat(fn, &buf)&&S_ISLNK(buf.st_mode);
}
/**
   Test whether fd is a socket
*/
int issock(int fd){
	if(fd==-1) return 0;
	struct stat buf;
	return !fstat(fd, &buf)&&S_ISSOCK(buf.st_mode);
}
/**
 * Compute length of file in Bytes
 */
size_t flen(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return stat(fn, &buf)?0:(size_t)buf.st_size;
}
/**
   Return the modification time of the file
 */
time_t fmtime(const char* fn){
	struct stat buf;
	if(!fn||stat(fn, &buf)) return 0;
	return buf.st_mtime;
}
/**
   Concatenate many strings. Argument list must end with NULL.
*/
char* stradd(const char* a, ...){
	char* out;
	size_t n=0;
	va_list ap;
	va_start(ap, a);
	if(a){
		n=strlen(a)+1;
	}
	for(const char* arg=va_arg(ap, const char*); arg; arg=va_arg(ap, const char*)){
		n+=strlen(arg);
	}
	va_end(ap);
	out=mycalloc(n, char);
	if(a){
		strcpy(out, a);
	}
	va_start(ap, a);
	for(const char* arg=va_arg(ap, const char*); arg; arg=va_arg(ap, const char*)){
		strcat(out, arg);
	}
	va_end(ap);
	return out;
}
/**
   Concatenate many strings, like stradd, but arguments are an array of char*
*/
char* strnadd(int argc, const char* argv[], const char* delim){
	if(!argc) return NULL;
	size_t slen=1+strlen(delim);
	for(int iarg=0; iarg<argc; iarg++){
		slen+=strlen(delim)+strlen(argv[iarg]);
	}
	char* scmd=mycalloc(slen, char);
	for(int iarg=0; iarg<argc; iarg++){
		if(argv[iarg]&&strlen(argv[iarg])>0){
			strcat(scmd, argv[iarg]);
			strcat(scmd, delim);
		}
	}
	char* scmdend=scmd+strlen(scmd)-1;
	while(scmdend>scmd+1&&scmdend[0]==delim[0]&&scmdend[-1]==delim[0]){
		scmdend[0]='\0';
		scmdend--;
	}
	return scmd;
}
/**
   translate a filename into absolute file name that starts with /
*/
char* expand_filename(const char* fn){
	char* out;
	if(fn[0]=='~'){
		out=stradd(HOME, fn+1, NULL);
	} else{
		out=strdup(fn);
	}
	return out;
}
/**
   Duplicate a string. Check for NULL.  Do not call strndup to avoid recursive deadlock.
 */
char* mystrndup(const char* A, size_t len){
	if(!len) return NULL;
	size_t len2=strlen(A);
	if(!len2) return NULL;
	if(len2<len) len=len2;
	char* B=(char*)malloc(len+1);
	memcpy(B, A, len);
	B[len]='\0';
	return B;
}

/**
   declare strdup so my memory mangement mem.c is happy when DEBUG=1. Handles
NULL pointer correctly.  Do not call strdup to avoid recursive deadlock. */
char* mystrdup(const char* A){
	if(!A){
		return NULL;
	} else{
		size_t nlen=strlen(A);
		char* B=(char*)malloc(nlen+1);
		memcpy(B, A, nlen+1);
		return B;
	}
}
/**
 * @brief replaces snprintf to check for return and print a warning if there is truncation.
 * 
 */
int mysnprintf(char* restrict str, size_t size, const char *restrict format, ...){
	va_list ap;
	va_start(ap, format);
	int n=vsnprintf(str, size, format, ap);
	va_end(ap);	
	if(n<0){
		warning("snprintf failed with error %d\n", n);
	}else if(n>(ssize_t)size){
		warning("snprintf is truncated: %s (need %d, got %zu)\n", str, n, size);
	}
	return n;
} 

/**
   Remove files that are older than sec seconds in folder fndir. If sec==0,
   remove everything.
*/
void remove_file_older(const char* fndir, int level, long sec){
	DIR* dir=opendir(fndir);
	if(!dir){
		warning("Unable to open directory %s\n", fndir);
		return;
	} else{
	//info("Cleaning %s\n", fndir);
	}
	struct dirent* dp;
	struct stat buf;
	char fnfull[PATH_MAX];
	long sec2=myclocki()-sec;
	while((dp=readdir(dir))){
		if(dp->d_name[0]=='.') continue;
		snprintf(fnfull, PATH_MAX, "%s/%s", fndir, dp->d_name);
		if(!stat(fnfull, &buf)){
			if(S_ISDIR(buf.st_mode)){
				if(level!=0) remove_file_older(fnfull, level-1, sec);
			} else if(S_ISREG(buf.st_mode)&&(buf.st_atime<=sec2||sec==0)){
				if(check_suffix(fnfull, ".bin")||check_suffix(fnfull, ".lock")){
					remove(fnfull);
					//info("Remove %s. %ld days old\n", fnfull, (long)(myclocki()-buf.st_mtime)/3600/24);
				}
			}
		}
	}
	closedir(dir);
}

/**
   Make dirs recursively. like mkdir -p in bash
*/
void mymkdir(const char* format, ...){
	format2fn;
	if(!fn) return;
	while(fn[strlen(fn)-1]=='/')
		fn[strlen(fn)-1]='\0';
	if(mkdir(fn, 0777)==-1&&errno!=EEXIST){
	//perror("mkdir");
		char* tmp=strrchr(fn, '/');
		if(!tmp){
			error("Unable to mkdir '%s'\n", fn);
		}
		tmp[0]='\0';
		mymkdir("%s", fn);

		tmp[0]='/';
		if(mkdir(fn, 0777)==-1&&errno!=EEXIST){
			error("Unable to mkdir '%s'\n", fn);
		}
	}
}
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.
 */
int mystrcmp(const char* a, const char* b){
	if(!a||!b) return 1;
	size_t la=strlen(a);
	size_t lb=strlen(b);
	if(la==0||lb==0||la<lb){
		return 1;
	} else{
		return strncmp(a, b, lb);
	}
}

/**
   Make the fd close on exec.
*/
void cloexec(int fd){
	int oldflag=fcntl(fd, F_GETFD, 0);
	if(oldflag!=-1){
		oldflag|=FD_CLOEXEC;
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
* 	Pause execution and listen input from fd1 and/or fd2 for continuation. Returns new pause flag.
 */
int mypause(int fd1,/**<first file no, usually 0 for stdin*/
			int fd2 /**<second file no, usually created from a pipe t*/
){
	int ans=1;
	errno=0;
	info2("Press enter to step, c to resume:");
	fd_set active_fd_set;
	FD_ZERO(&active_fd_set);
	if(fd1>0|| (fd1==0 && !detached) ){
		FD_SET(fd1, &active_fd_set);
	}
	if(fd2>0 && fd2 != fd1){
		FD_SET(fd2, &active_fd_set);
	}
	info_errno("mypause");
	while(1){
		int select_ans=select(FD_SETSIZE, &active_fd_set, NULL, NULL, 0);
		if(select_ans<0){
			warning("Select failed: %s\n", strerror(errno));
			return ans;
		}
		info_errno("select");
		for(int i=0; i<FD_SETSIZE; i++){
			if(FD_ISSET(i, &active_fd_set)){
				char key;
				int nread=read(i, &key, 1);
				info_errno("read");
				if(nread==1){
					if(key=='c'){
						info2("Resuming\n");
						return 0;
					} else{
						if(key>='0' && key<='9'){
							ans=key-'0';
						}
					}
				}
				info2("Continuing...\n");
				return ans;
			}
		}
	}
	return ans;
}
/**
   Return available space of mounted file system in bytes.
*/
long available_space(const char* path){
	struct statvfs buf;
	if(statvfs(path, &buf)){
		perror("statvfs");
		return 0;
	} else{
		return (long)buf.f_bsize*(long)buf.f_bavail;
	}
}
/**
   Extract a string constant from the command line, and output the position
   where the string terminates.*/
static char* cmd_string(char* start, char** end2){
	char* input=start;
	char* end;
	while(isspace((int)input[0])||input[0]=='\n') input++;
	if(input[0]=='\''||input[0]=='"'){
		end=strchr(input+1, input[0]);/*find matching quote. */
		input[0]=' ';
		input++;
		if(!end){
			error("String does not end\n");
		}
	} else{
		end=input;
		while(!isspace((int)end[0])&&end[0]!='\n'&&end[0]!='\0') end++;
	}
	int noteos=0;
	if(end[0]!='\0'){
		end[0]='\0';
		noteos=1;
	}
	char* out=strdup(input);
	memset(input, ' ', strlen(input));
	if(noteos){
		end[0]=' ';
		*end2=end+1;
	} else{
		*end2=end;
	}
	return out;
}
/**
   Parse command line arguments. The remaining string contains whatever is not yet parsed.
   This is more relaxed than the built in getopd
*/
void parse_argopt(char* cmds, argopt_t* options){
	char* cmds_end=cmds+(cmds?strlen(cmds):0);
	char* start=cmds;
	while(start<cmds_end){
		if(isspace((int)start[0])||start[0]=='\n'){
			start[0]=' ';
			start++;
			continue;
		}
		if(options&&start[0]=='-'){
			char* start0=start;
			char key='0';
			char* value;
			int iopt=-1;
			start++;
			if(start[0]=='-'){/*long option, replace with short ones. */
				start++;
				for(int i=0; (options[i].name); i++){
					if(!mystrcmp(start, options[i].name)){
						key=options[i].key;
						start+=strlen(options[i].name);
						while(isspace((int)start[0])||start[0]=='\n'){
							start[0]=' ';
							start++;
						}
						if(start[0]=='='){
							start[0]=' ';
							start++;
						}
						iopt=i;
						break;
					}
				}
			} else{
				key=start[0];
				start++;
				for(int i=0; (options[i].name); i++){
					if(key==options[i].key){
						iopt=i;
						break;
					}
				}
			}
			if(iopt==-1){
				continue;/*we don't want this key. */
			}
			if((options[iopt].valtype)){//expects a value.
				value=start;
				while(value[0]=='\n'||isspace((int)value[0])){
					value[0]=' ';
					value++;
				}
				if(value[0]=='\0'){
					value=NULL;
				}
			} else{
				value=NULL;
			}
			int isfun=(options[iopt].isfun);
			switch(options[iopt].type){
			case 0:/*no result needed */
				break;
			case M_INT:{
				if(options[iopt].valtype==2){//needs an array
					if(isfun) error("Not implemented yet\n");
					int val=(int)strtol(value, &start, 10);
					int** tmp=(int**)options[iopt].val;
					int* nval=(int*)options[iopt].nval;
					int i=*nval;
					//uncomment the following to avoid duplicates
					/*for(i=0; i<*nval; i++){
						if((*tmp)[i]==val) break;
					}*/
					if(i==*nval){
						(*nval)++;
						*tmp=myrealloc(*tmp, (size_t)*nval, int);
						(*tmp)[(*nval)-1]=val;
					}
				} else{
					int val=value?(int)strtol(value, &start, 10):1;
					if(value==start) val=1;//no entry, default to 1.
					if(isfun){/*Is function */
						void (*tmp)(int)=(void (*)(int))options[iopt].val;
						tmp(val);
					} else{
						int* tmp=(int*)options[iopt].val;
						*tmp=val;
					}
				}
			}
					  break;
			case M_DBL:{
				if(options[iopt].valtype==2){//needs an array
					if(isfun) error("Not implemented yet\n");
					double val=strtod(value, &start);
					double** tmp=(double**)options[iopt].val;
					int* nval=(int*)options[iopt].nval;
					(*nval)++;
					*tmp=myrealloc(*tmp, (size_t)*nval, double);
					(*tmp)[(*nval)-1]=(int)val;
				} else{
					double val=value?strtod(value, &start):1;
					if(isfun){/*Is function */
						void (*tmp)(double)=(void (*)(double))options[iopt].val;
						tmp(val);
					} else{
						double* tmp=(double*)options[iopt].val;
						*tmp=val;
					}
				}
			}
					  break;
			case M_STR:{
				char* val=value?cmd_string(value, &start):strdup("Unknown");
				if(isfun){
					void (*tmp)(char*)=(void (*)(char*))options[iopt].val;
					tmp(val);
					free(val);
				} else{
					char** tmp=(char**)options[iopt].val;
					free(*tmp); *tmp=val;
				}
			}
					  break;
			default:
				error("Unknown type");
			}/*switch */
			/*Empty the string that we already parsed. */
			memset(start0, ' ', start-start0);
		} else if(start[0]=='='){/*equal sign found, key=value */
			/*create a \n before the key. */
			int skipspace=1;
			for(char* start2=start-1; start2>=cmds; start2--){
				if(isspace((int)*start2)||*start2=='\n'){
					if(!skipspace){
						*start2='\n';
						break;
					}
				} else{
					skipspace=0;
				}
			}
			start++;
		} else if(!mystrcmp(start, ".conf")){ /*.conf found. */
			/*create a \n before the key. and \n after .conf */
			for(char* start2=start-1; start2>=cmds; start2--){
				if(isspace((int)*start2)||*start2=='\n'){
					//check whether -c is before *.conf
					char *start3=start2;
					//skip continuous space
					while(start3-1>=cmds && isspace((int)start3[-1])){
						start3--;
					}
					if(start3-2>=cmds&&start3[-2]=='-'&&start3[-1]=='c'){
						if(start2-3>=cmds){
							start2[-3]='\n';
						}
					}else{
						*start2='\n';
					}
					break;
				}
			}
			start+=5;
			start[0]='\n';
			start++;
		} else if(start[0]=='['){/*make sure we don't split brackets that are part of value. */
			const char *save_start=start;
			int nopen=1; //opening brackets not closed.
			for(start++; start[0] && nopen; start++){
				if(start[0]=='['){
					nopen++;
				}else if(start[0]==']'){
					nopen--;
				}
			}
			if(!start[0] && nopen){
				error("Bracket is not closed: {%s} from {%s}\n", save_start, cmds);
			}
			/*
			char* bend=strchr(start+1, ']');
			char* bnextstart=strchr(start+1, '[');
			if(bend){
				if(!bnextstart||bend<bnextstart)){
					for(; start<bend+1; start++){
						if(start[0]=='\n') start[0]=' ';
					}
				}else{

				}
			} else{
				
				start++;
			}*/
		} else if(start[0]=='\''||start[0]=='"'){/*make sure we don't split strings that are part of value. */
			char* quoteend=strchr(start, start[0]);
			if(quoteend){
				start=quoteend+1;
			} else{
				warning("Quote is not closed\n");
				start++;
			}
		} else if(start[0]=='-' && start[1]=='o'){//-o dir
			if(start-1>=cmds) start[-1]='\n';
			start+=2;
		} else{
			start++;
		}
	}
}

/**
   Set scheduling priorities for the process to enable real time behavior.
*/
void set_realtime(int icpu, int niceness){
	(void)icpu;
	(void)niceness;
	//Set CPU affinity.
	/*
	  //Deprecated. Use external tools to do so, like openmp env's
#ifdef __linux__
	if(icpu>0){
	cpu_set_t cpuset={{0}};
	CPU_SET(icpu, &cpuset);
	sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
	}
	//lock data in memory, avoid swapping.
	mlockall(MCL_FUTURE | MCL_CURRENT);
#endif
	*/
	//fail stack
	struct rlimit rl;
	if(!getrlimit(RLIMIT_STACK, &rl)){
		const int NSTACK=(int)(rl.rlim_cur/2);
		char tmp[NSTACK];
		memset(tmp, 0, NSTACK);
	}
	//Set only if we are root.
	if(getuid()==0){
		info("Set priority to -20\n");
		setpriority(PRIO_PROCESS, (id_t)getpid(), -20);//is this necessary?
#ifdef __linux__
		struct sched_param param;
		sched_getparam(getpid(), &param);
		param.sched_priority=sched_get_priority_max(SCHED_FIFO)-1;
		sched_setscheduler(getpid(), SCHED_FIFO, &param);
#endif
	} else{
		warning("Please run program as setsid or as root to lift priority\n");
	}
}
void free_strarr(char **str, int n){
	if(str){
		for(int i=0; i<n; i++){
			free(str[i]);
		}
		free(str);
	}
}
#undef strdup
char* (*strdup0)(const char*)=strdup;
const int default_color_table[]={0x0000FF,
			   0xFF0000,
			   0x00FF00,
			   0x009999,
			   0x00FFFF,
			   0x9900CC,
			   0xFFCC00,
			   0xFF00FF,
			   0x000000,
			   0x666666,
};
void print_version(void){
	extern const char *GIT_VERSION;
	info2("SRC: %s v%s %s\n", SRCDIR, PACKAGE_VERSION, GIT_VERSION);
	char exe[PATH_MAX];
	if(!get_job_progname(exe, PATH_MAX, 0)){
		info2("BUILT: %s by %s on %s", BUILDDIR, COMPILER, myasctime(fmtime(exe)));
	} else{
		info2("BUILT: %s by %s on %s %s", BUILDDIR, COMPILER, __DATE__, __TIME__);//__DATE__ and __TIME__ is only applicable to this specific file
	}
#ifdef __OPTIMIZE__
#define OPT_STR "+O3"
#else
#define OPT_STR "+O0"
#endif
#if CPU_SINGLE
#define CPU_FP "F32"
#else
#define CPU_FP "F64"
#endif
	info2(" CPU(" CPU_FP "," OPT_STR ")");
#if USE_CUDA
#if CUDA_DOUBLE
#define GPU_FP "F64"
#else
#define GPU_FP "F32"
#endif
	info2(" with CUDA(v%d," GPU_FP ")\n", USE_CUDA);
#else
	info2(" w/o CUDA\n");
#endif
	info("Launched at %s in %s with PID %ld.\n", myasctime(0), HOST, (long)getpid());
#if !MAOS_DISABLE_SCHEDULER
	extern uint16_t PORT;
	info("The browser based job monitor and drawdaemon can be accessed at http://localhost:%d\n", PORT);
#endif	
}
void mystrrep(char *str, const char *prefix, const char *substitute){
	if(!str || !prefix || !substitute) return;
	size_t nsub=strlen(substitute);
	size_t npre=strlen(prefix);
	size_t nstr=strlen(str);
	if(nstr<npre || nsub>npre){
		return;
	}
	if(!mystrcmp(str, prefix)){
		memcpy(str, substitute, nsub);
		memmove(str+nsub, str+npre, nstr-npre+1);
	}
}
/**
 * @brief Convert time in seconds to string representation with minimal length
 * 
 * @param tmp 	target memory
 * @param stmp 	length of target memory
 * @param sec 	time in seconds
 */
int sec2str(char*tmp, long stmp, double sec){
	if(stmp<0) return 0;
	if(sec<0) sec=0; //ignore negative
	int offset=0;
	long hr=sec/3600; sec-=hr*3600;
	long m=sec/60; sec-=m*60;
	if(hr){
		offset+=snprintf(tmp+offset, stmp-offset, "%ldh", hr);
		sec=0;//ignore remaining seconds
		if(hr>5){
			m=0;
		}
	}
	if(m){
		offset+=snprintf(tmp+offset, stmp-offset, "%ldm", m);
		if(m>5){
			sec=0;//ignore seconds
		}else{
			sec=round(sec);//ignore fraction seconds
		}
	}
	if(sec>10){
		offset+=snprintf(tmp+offset, stmp-offset, "%.0f", sec);
	}else if(sec>0){
		offset+=snprintf(tmp+offset, stmp-offset, "%.2g", sec);
	}else if(offset==0){
		offset+=snprintf(tmp+offset, stmp-offset, "0");
	}
	return offset;
}
/**
 * @brief Rename prefix_host_pid.suffix to prefix_done.suffix. 
 * If already exists and count_in>=0 increment the counter and use prefix_done.counter.suffix as filename
 * 
 * @param prefix 
 * @param suffix 
 * @param count_in if >=0, do not increase
 * @return int the counter
 */
static int rename_done(const char *prefix, const char *suffix, int count_in){
	char fn[PATH_MAX];
	char fn1[PATH_MAX];
	char fn2[PATH_MAX];
	int count=count_in<0?0:count_in;
	snprintf(fn, sizeof(fn), "%s_%s_%ld.%s", prefix, HOST, (long)getpid(), suffix);
	if(!exist(fn)){
		dbg("%s does not exist\n", fn);
		return -1;
	}
	do{
		if(count==0){
			snprintf(fn2, sizeof(fn2), "%s_done.%s", prefix, suffix);
		} else{
			snprintf(fn2, sizeof(fn2), "%s_done.%d.%s", prefix, count, suffix);
		}
	} while(exist(fn2) && count_in<0 && (count=count+1));
	snprintf(fn1, sizeof(fn1), "%s_done.%s", prefix, suffix);
	if(count>0){
		if(rename(fn1, fn2)){
			dbg("Rename %s to %s failed\n", fn1, fn2);
		}
	}
	if(rename(fn, fn1)){
		dbg("Rename %s to %s failed\n", fn, fn2);
	} else{
		snprintf(fn, sizeof(fn), "%s_recent.%s", prefix, suffix);
		remove(fn);
		mysymlink(fn1, fn);
	}
	return count;
}
void rename_log(int sig, const char *exe){
	if(sig==0){//success
		int count=rename_done("run", "log", -1);
		if(count!=-1){
			rename_done(exe, "conf", count);
		}
	}else{//killed or error
		char fn[PATH_MAX];
		snprintf(fn, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
		char fn2[PATH_MAX];
		snprintf(fn2, PATH_MAX, "run_%s_%ld.err", HOST, (long)getpid());
		(void)rename(fn, fn2);
		mysymlink(fn2, "run_recent.log");
	}
}
