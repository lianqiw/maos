/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include <ctype.h>
#include "io.h"
static void write_timestamp(file_t *fp);

static const char *myasctime(void){
    static char st[64];
    time_t a;
    time(&a);
    ctime_r(&a, st);
    st[strlen(st)-1]='\0';
    return st;
}
static int exist(const char *fn){
    /**
       Test whether a file exists.
    */
    struct stat buf;
    return !stat(fn, &buf);
}
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.
 */
static int mystrcmp(const char *a, const char *b){
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
   Check the suffix of a file.
*/
static int check_suffix(const char *fn, const char *suffix){
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
/**
   Test whether fn is a symbolic link
*/
static int islink(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISLNK(buf.st_mode);
}
static char* procfn(const char *fn, const char *mod,const int gzip){
    if(!fn){
	info("fn is empty\n");
	return NULL;
    }
    char *fn2;
    if(fn[0]=='~'){
	char *HOME=getenv("HOME");
	fn2=malloc(strlen(HOME)+strlen(fn)+16);
	strcpy(fn2,HOME);
	strcat(fn2,fn+1);
    }else{
	fn2=malloc(strlen(fn)+16);
	strcpy(fn2,fn);
    } 
    /*If there is no recognized suffix, add .bin in the end.*/
    int nosuffix=0;
    if(!check_suffix(fn2,".bin") && !check_suffix(fn2, ".bin.gz")
       &&!check_suffix(fn2,".fits") && !check_suffix(fn2, ".fits.gz")){
	strncat(fn2, ".bin", 4);
	nosuffix=1;
    }
    if(mod[0]=='r' || mod[0]=='a'){
	if(!exist(fn2)){/*If does not exist.*/
	    if(!check_suffix(fn2, ".gz")){
		strncat(fn2, ".gz", 3);
	    }else if (check_suffix(fn, ".gz")){
		/*ended with bin.gz, change to .bin*/
		fn2[strlen(fn2)-3]='\0';
	    }
	    if(!exist(fn2)){
		if(nosuffix){
		    fn2[strlen(fn2)-7]='\0';/*remove the added .bin, .gz*/
		    strncat(fn2, ".fits", 5);/*replace with .fits*/
		    if(!exist(fn2)){
			return NULL;
		    }
		}else{
		    return NULL;
		}
	    }
	}
    }else if (mod[0]=='w'){
	if(islink(fn2)){
	    /*remove old file to avoid write over a symbolic link.*/
	    if(remove(fn2)){
		info("Failed to remove %s\n", fn2);
		return NULL;
	    }
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
}
/*stripped down version of io.c*/
file_t* zfopen(const char *fn, char *mod){
    char *fn2=procfn(fn, mod, 1);
    if(!fn2){
	info("%s does not exist\n", fn);
	return NULL;
    }
    file_t* fp=calloc(1, sizeof(file_t));
    /*check fn instead of fn2. if end of .bin or .fits, disable compressing.*/
    /*Now open the file to get a fd number that we can use to lock on the
      file.*/
    switch(mod[0]){
    case 'r':/*read only */
	fp->fd=open(fn2, O_RDONLY);
	break;
    case 'w':/*write */
    case 'a':
	fp->fd=open(fn2, O_RDWR | O_CREAT, 0600);
	if(fp->fd!=-1 && flock(fp->fd, LOCK_EX|LOCK_NB)){
	    info("Trying to write to a file that is already opened for writing: %s\n", fn2);
	    close(fp->fd);
	    return NULL;
	}else{
	    if(mod[0]=='w' && ftruncate(fp->fd, 0)){/*Need to manually truncate the file. */
		warning("Truncating %s failed\n", fn2);
	    }
	}
	break;
    default:
	info("Unknown mod=%s\n", mod);
	return NULL;
    }
    if(fp->fd==-1){
	info("Unable to open file %s\n", fn2);
	return NULL;
    }
    if(mod[0]=='w'){
	if(check_suffix(fn2, ".bin") || check_suffix(fn2, ".fits")){
	    fp->isgzip=0;
	}else{
	    fp->isgzip=1;
	}
    }else{ 
	uint16_t magic;
	if(read(fp->fd, &magic, sizeof(uint16_t))!=sizeof(uint16_t)){
	    info("Read magic failed\n");
	    close(fp->fd); return NULL;
	}
	if(magic==0x8b1f){
	    fp->isgzip=1;
	}else{
	    fp->isgzip=0;
	}
	lseek(fp->fd, 0, SEEK_SET);
    }
    if(fp->isgzip){
	if(!(fp->p=gzdopen(fp->fd,mod))){
	    info("Error gzdopen for %s\n",fn2);
	    close(fp->fd); return NULL;
	}
    }else{
	if(!(fp->p=fdopen(fp->fd,mod))){
	    info("Error fdopen for %s\n",fn2);
	    close(fp->fd); return NULL;
	}
    }
    if(check_suffix(fn2, ".fits") || check_suffix(fn2, ".fits.gz")){
	fp->isfits=1;
    }  
    if(mod[0]=='w' && !fp->isfits){
	write_timestamp(fp);
    }
    free(fn2);
    return fp;
}
void zfclose(file_t *fp){
    if(fp->isgzip){
	gzclose((voidp)fp->p);
    }else{
	if(fclose((FILE*)fp->p)){
	    perror("fclose\n");
	}
    }
    free(fp);
}
static inline void zfwrite_do(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    if(fp->isgzip){
	if(gzwrite((voidp)fp->p, ptr, size*nmemb)!=size*nmemb){
	    perror("gzwrite");
	    error("write failed\n");
	}
    }else{
	if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    perror("fwrite");
	    error("writefailed\n");
	}
    }
}
/**
   Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
   fwrite. Follows the interface of fwrite.
 */
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits && size>1){
	/* write a block of 2880 bytes each time, with big-endianness.*/
	const int bs=2880;
	char junk[bs];
	int length=size*nmemb;
	int nb=(length+bs-1)/bs;
	char *in=(char*)ptr;
	for(int ib=0; ib<nb; ib++){
	    int nd=length<bs?length:bs;
	    switch(size){
	    case 2:
		for(int i=0; i<nd; i+=2){
		    junk[i]=in[i+1];
		    junk[i+1]=in[i];
		}
		break;
	    case 4:
		for(int i=0; i<nd; i+=4){
		    junk[i]  =in[i+3];
		    junk[i+1]=in[i+2];
		    junk[i+2]=in[i+1];
		    junk[i+3]=in[i  ];
		}
		break;
	    case 8:
	    case 16:
		for(int i=0; i<nd; i+=8){
		    junk[i  ]=in[i+7];
		    junk[i+1]=in[i+6];
		    junk[i+2]=in[i+5];
		    junk[i+3]=in[i+4];
		    junk[i+4]=in[i+3];
		    junk[i+5]=in[i+2];
		    junk[i+6]=in[i+1];
		    junk[i+7]=in[i  ];
		}
		break;
	    default:
		error("Invalid\n");
	    }
	    /* use bs instead of nd to test tailing blanks*/
	    in+=bs; length-=bs;
	    if(length<0){
		memset(junk+nd, 0, (bs-nd)*sizeof(char));
	    }
	    zfwrite_do(junk, sizeof(char), bs, fp);
	}
    }else{
	zfwrite_do(ptr, size, nmemb, fp);
    }
}
void zfwrite_dcomplex(const double* pr, const double *pi,const size_t nmemb, file_t *fp){
    dcomplex *tmp=malloc(sizeof(dcomplex)*nmemb);
    long i;
    for(i=0; i<nmemb; i++){
	tmp[i].x=pr[i];
	tmp[i].y=pi[i];
    }
    zfwrite(tmp, sizeof(dcomplex), nmemb, fp);
    free(tmp);
}
void zfwrite_fcomplex(const float* pr, const float *pi,const size_t nmemb, file_t *fp){
    fcomplex *tmp=malloc(sizeof(fcomplex)*nmemb);
    long i;
    for(i=0; i<nmemb; i++){
	tmp[i].x=pr[i];
	tmp[i].y=pi[i];
    }
    zfwrite(tmp, sizeof(fcomplex), nmemb, fp);
    free(tmp);
}
static inline int zfread_do(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(fp->isgzip){
	return gzread((voidp)fp->p, ptr, size*nmemb)>0?0:-1;
    }else{
	return fread(ptr, size, nmemb, (FILE*)fp->p)==nmemb?0:-1;
    }
}
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread. It does byte ordering from big endian
   to small endian in case we are reading fits file.  */
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits && size>1){/*need to do byte swapping.*/
	const long bs=2880;
	char junk[bs];
	long length=size*nmemb;
	long nb=(length+bs-1)/bs;
	char *out=(char*)ptr;
	for(int ib=0; ib<nb; ib++){
	    if(zfread_do(junk, sizeof(char), bs, fp)) return -1;
	    int nd=length<bs?length:bs;
	    switch(size){
	    case 2:
		for(int i=0; i<nd; i+=2){
		    out[i]  =junk[i+1];
		    out[i+1]=junk[i];
		}
		break;
	    case 4:
		for(int i=0; i<nd; i+=4){
		    out[i]  =junk[i+3];
		    out[i+1]=junk[i+2];
		    out[i+2]=junk[i+1];
		    out[i+3]=junk[i  ];
		}
		break;
	    case 8:
	    case 16:
		for(int i=0; i<nd; i+=8){
		    out[i  ]=junk[i+7];
		    out[i+1]=junk[i+6];
		    out[i+2]=junk[i+5];
		    out[i+3]=junk[i+4];
		    out[i+4]=junk[i+3];
		    out[i+5]=junk[i+2];
		    out[i+6]=junk[i+1];
		    out[i+7]=junk[i  ];
		}
		break;
	    default:
		error("Invalid\n");
	    }
	    out+=bs;
	    length-=bs;
	}
	return 0;
    }else{
	return zfread_do(ptr, size, nmemb, fp);
    }
}

/**
   Move the current position pointer, like fseek
*/
int zfseek(file_t *fp, long offset, int whence){
    if(fp->isfits){/*in fits, we need to offset integer block of 2880.*/
	const long bs=2880;
	long nb=(offset+bs-1)/bs;
	offset=nb*bs;
    }
    if(fp->isgzip){
	return gzseek((voidp)fp->p,offset,whence)==-1?-1:0;
    }else{
	return fseek((FILE*)fp->p,offset,whence);
    }
}
long zftell(file_t *fp){
    if(fp->isgzip){
	return gztell((voidp)fp->p);
    }else{
	return ftell((FILE*)fp->p);
    }
}
int zfeof(file_t *fp){
    return zfseek(fp, 1, SEEK_CUR);
}
/**
   Write the magic into file. Also write a dummy header to make data alignment to 8 bytes.
*/
static void write_bin_magic(uint32_t magic, file_t *fp){
    if(fp->isfits) error("fits file is not supported\n");
    uint32_t magic2=M_SKIP;
    zfwrite(&magic2, sizeof(uint32_t), 1, fp);
    zfwrite(&magic,  sizeof(uint32_t), 1, fp);
}
/**
   Append to the file the header to the end of the file(or rather, the
   tailer). First write magic number, then the length of the header, then the
   header, then the length of the header again, then the magic number again. The
   two set of strlen and header are used to identify the header from the end of
   the file and also to verify that what we are reading are indeed header. The
   header may be written multiple times. They will be concatenated when
   read. The header should contain key=value entries just like the configuration
   files. The entries should be separated by new line charactor. */
static void write_bin_header(const char *header, file_t *fp){
    if(!header) return;
    if(fp->isfits) error("fits file is not supported\n");
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(header)+1;
    /*make header 8 byte alignment.*/
    char *header2=strdup(header);
    if(nlen % 8 != 0){
	nlen=(nlen/8+1)*8;
	header2=realloc(header2, nlen);
    }
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(header2, 1, nlen, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    free(header2);
}

static void write_timestamp(file_t *fp){
    if(fp->isfits) error("Not supported\n");
    char header[128];
    snprintf(header,128, "Created by write on %s\n",
	     myasctime());
    write_bin_header(header, fp);
}

/**
   Obtain the current magic number. If it is a header, read it out if output of
header is not NULL.  The header will be appended to the output header.*/
static uint32_t read_bin_magic(file_t *fp, char **header){
    uint32_t magic,magic2;
    uint64_t nlen, nlen2;
    if(fp->isfits) error("fits file is not supported\n");
    while(1){
	/*read the magic number.*/
	if(zfread2(&magic, sizeof(uint32_t), 1, fp)) return 0;
	/*If it is header, read or skip it.*/
	if(magic==M_SKIP){
	    continue;
	}else if(magic==M_HEADER){
	    zfread(&nlen, sizeof(uint64_t), 1, fp);
	    if(nlen>0){
		if(header){
		    char header2[nlen];
		    zfread(header2, 1, nlen, fp);
		    header2[nlen-1]='\0'; /*make sure it is NULL terminated.*/
		    if(*header){
			*header=realloc(*header, ((*header)?strlen(*header):0)+strlen(header2)+1);
			strncat(*header, header2, nlen);
		    }else{
			*header=strdup(header2);
		    }
		}else{
		    zfseek(fp, nlen, SEEK_CUR);
		}
	    }
	    zfread(&nlen2, sizeof(uint64_t),1,fp);
	    zfread(&magic2, sizeof(uint32_t),1,fp);
	    if(magic!=magic2 || nlen!=nlen2){
	      error("Header verification failed: %u<>%u, %lu<>%lu\n", magic, magic2, (unsigned long)nlen, (unsigned long)nlen2);
	    }
	}else{ /*otherwise return the magic number*/
	    return magic;
	}
    }/*while*/
}
/**
   Write fits header. extra is extra header that will be put in fits comment
*/
static void
write_fits_header(file_t *fp, const char *str, uint32_t magic, int count, ...){
    uint64_t naxis[count];
    va_list ap;
    va_start (ap, count);              /*Initialize the argument list. */
    int empty=0;
    for (int i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	if((*addr)==0) empty=1;
	naxis[i]=*addr;
    }
    va_end(ap);

    if(empty) count=0;

    int bitpix;
    switch(magic){
    case M_FLT:
	bitpix=-32;
	break;
    case M_DBL:
	bitpix=-64;
	break;
    case M_INT32:
	bitpix=32;
	break;
    case M_INT16:
	bitpix=16;
	break;
    case M_INT8:
	bitpix=8;
	break;
    default:
	bitpix=0;
	error("Data type is not yet supported.\n");
    }
    const int nh=36;
    char header[nh][80];
    memset(header, ' ', sizeof(char)*36*80);
    int hc=0;
    if(fp->isfits==1){
	snprintf(header[hc], 80, "%-8s= %20s", "SIMPLE", "T");    header[hc][30]=' '; hc++;
	fp->isfits++;
    }else{
	snprintf(header[hc], 80, "%-8s= %s", "XTENSION", "'IMAGE   '");    header[hc][20]=' '; hc++;
    }
    snprintf(header[hc], 80, "%-8s= %20d", "BITPIX", bitpix); header[hc][30]=' '; hc++;
    snprintf(header[hc], 80, "%-8s= %20d", "NAXIS", count);   header[hc][30]=' '; hc++;
#define FLUSH_OUT /*write the block and reset */	\
	if(hc==nh){					\
	    zfwrite(header, sizeof(char), 36*80, fp);	\
	    memset(header, ' ', sizeof(char)*36*80);	\
	    hc=0;					\
	}
    for (int i = 0; i < count; i++){
	FLUSH_OUT;
	snprintf(header[hc], 80, "%-5s%-3d= %20lu", "NAXIS", i+1, (unsigned long)(naxis[i])); header[hc][30]=' '; hc++;
    }
    if(str){/*We wrap our header in COMMENT section to avoid dealing with name
	     * length limit*/
	const char *str2=str+strlen(str);
	while(str<str2){
	    char *nl=strchr(str, '\n');
	    int length;
	    if(nl){
		length=nl-str+1;
		if(length<=70) nl[0]=';';
	    }else{
		length=strlen(str);
	    }
	    if(length>70) length=70;
	    FLUSH_OUT;
	    strncpy(header[hc], "COMMENT   ", 10);
	    strncpy(header[hc]+10, str, length);
	    hc++;
	    str+=length;
	}
    }
    FLUSH_OUT;
    snprintf(header[hc], 80, "%-8s", "END"); header[hc][8]=' '; hc++;
    zfwrite(header, sizeof(char), 36*80, fp);
#undef FLUSH_OUT
}

/**
   Read fits header
 */
int read_fits_header(file_t *fp, char **str, uint32_t *magic, uint64_t *nx, uint64_t *ny){
    char line[82];//extra space for \n and \0
    int end=0;
    int page=0;
    int bitpix=0;
    int naxis=0;
    while(!end){
	int start=0;
	if(page==0){
	    /*First read mandatory fits headers.*/
	    if(zfread2(line, 1, 80, fp)) return -1; line[80]='\0';
	    if(strncmp(line, "SIMPLE", 6) && strncmp(line, "XTENSION= 'IMAGE", 16)){
		info("Garbage in fits file at %ld:\n", zftell(fp));
		info("%s\n", line);
		return -1;
	    }
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &bitpix)!=1) error("Unable to determine bitpix\n");
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &naxis)!=1) error("Unable to determine naxis\n");
	    if(naxis>2) error("Data type not supported\n");
	    if(naxis>0){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (long unsigned int*)nx)!=1) error("Unable to determine nx\n");
	    }else{
		*nx=0;
	    }
	    if(naxis>1){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (long unsigned int*)ny)!=1) error("Unable to determine ny\n");
	    }else{
		*ny=1;
	    }
	    start=3+naxis;
	}
	for(int i=start; i<36; i++){
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(!strncmp(line, "END",3)){
		end=1;
	    }else{
		char *hh=line;
		int length=80;
		int newline=1;
		if(!strncmp(line, "COMMENT", 7)){
		    hh=line+10;
		    length-=10;
		    newline=0;
		}
		//Remove trailing space.
		for(int j=length-1; j>=0; j--){
		    if(isspace((int)hh[j])){
			hh[j]='\0';
			length--;
		    }else{
			if(newline){
			    hh[j+1]='\n';
			    hh[j+2]='\0';
			}
			break;
		    }
		}
		if(length>0){
		    if(*str){
			*str=realloc(*str, strlen(*str)+length+1+newline);
		    }else{
			*str=malloc(length+1+newline); (*str)[0]='\0';
		    }
		    strcat(*str, hh);
		}
	    }
	}
	page++;
    }
    switch(bitpix){
    case -32:
	*magic=M_FLT;
	break;
    case -64:
	*magic=M_DBL;
	break;
    case 32:
	*magic=M_INT32;
	break;
    case 16:
	*magic=M_INT16;
	break;
    case 8:
	*magic=M_INT8;
	break;
    default:
	error("bitpix=%d is not yet handled.\n", bitpix);
    }
    return 0;
}

/**
   Write multiple long numbers into the file. To write three numbers, a, b, c,
call with zfwritelarr(fp, 3, &a, &b, &c); */
void zfwritelarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              /*Initialize the argument list. */
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	zfwrite(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       /* Clean up.  */
}
/**
   Read multiple long numbers from the file. To read three numbers, a, b, c,
   call with zfreadlarr(fp, 3, &a, &b, &c);
 */
void zfreadlarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              /*Initialize the argument list. */
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  /*Get the next argument value.   */
	zfread(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       /* Clean up.  */
}
/**
  A unified header writing routine for .bin and .fits files. It write the array
information and string header if any.  */
void write_header(const header_t *header, file_t *fp){
    if(fp->isfits){
	if(!iscell(header->magic)){
	    write_fits_header(fp, header->str, header->magic, 2, &header->nx, &header->ny);
	}
    }else{
	if(header){
	    write_bin_header(header->str, fp);
	}
	write_bin_magic(header->magic, fp);
	zfwritelarr(fp, 2, &header->nx, &header->ny);
    }
}
/**
   A unified header reading routine for .bin and .fits files. It read the array
information and string header if any.  Return error signal.*/
int read_header2(header_t *header, file_t *fp){
    int ans;
    header->str=NULL;
    if(fp->isfits){
	ans=read_fits_header(fp, &header->str, &header->magic, &header->nx, &header->ny);
    }else{
	header->magic=read_bin_magic(fp, &header->str);
	if(header->magic==0){
	    ans=0;
	}else{
	    ans=0;
	    zfreadlarr(fp, 2, &header->nx, &header->ny);
	}
    }
    return ans;
}
/**
   calls read_header2 and abort if error happens.*/
void read_header(header_t *header, file_t *fp){
    if(read_header2(header, fp)){
	error("read_header failed\n");
    }
}
