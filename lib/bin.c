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

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <ctype.h> /*isspace */
#include "common.h"
#include "thread.h"
#include "process.h"
#include "bin.h"
#include "path.h"
#include "misc.h"
#include "readcfg.h"
/**
   \file bin.c Defines our custom file format .bin or zipped .bin.gz and the
   basic IO functions. All file read/write operators are through functions in
   this file.
 */
/*
  2009-10-01: switch from popen of gz to zlib to read/write compressed files.
  with default compression retio, the write time for compressed file is from 3
  to 10 times slower with a mean of 4.  the reduction factor is generally
  between 10 to 100 times in size with a mean of 20 while some is only 5.

  Take RL->M for example, without compression, it takes 3.81 seconds in
  uncompressed mode, the ondisk size of 402,333,476. it takes 14.02 seconds in
  compressed mode with on disk size of 44,791,294. On the other hand, manually
  compress the uncompressed file using gz takes 26.8 seconds, and the on disk
  size is 44,894,765, comparable to written by zlib. The speed is 28 Mib/s for
  compressed.

  The speed to write uncompressed file is near 250MiB/s which is obviously due
  to buffering. The speed to write compressed file is on the order of 20-53 MiB

  I write to give a compression level of 2, the speed for RLM is about doubled.
  The size is not about 56,775,376, not bad.

  Notice that the stored data by 64 machines using type long can not be read by
  32 bit due to different byte length. Modified the routine to first 64 bit
  read/write.  nx,ny, are uint64_t for 64 and 32 bit machines magic is uint32_t
  always to be backward compatible.

*/
/**
   contains information about opened files.
*/
struct file_t{
    void *p;
    int isgzip;/*Is the file zipped.*/
    int isfits;/*Is the file fits.*/
    char *fn;
    int fd;/*For locking. */
#if IO_TIMMING == 1
    struct timeval tv1;
#endif
};
/**
   Describes the information about mmaped data. Don't unmap a segment of mmaped
   memory, which causes the whole page to be unmapped. Instead, reference count
   the mmaped file and unmap the segment when the nref dropes to 1.
 */
struct mmap_t{
    void *p;  /**<points to the beginning of mmaped memory for this type of data.*/
    long n;   /**<length of mmaped memory.*/
    long nref;/**<Number of reference.*/
    int fd;   /**<file descriptor. close it if not -1.*/
};

/**
   Process the input file name and return file names that can be open to
   read/write. If the file name does not end with .bin or .bin.gz it will add to
   the end .bin or .bin.gz depending on the value of defaultgzip. For read only
   access, it will also look into the path for files.
 */
char* procfn(const char *fn, const char *mod, const int defaultgzip){
    char *fn2;
    if(fn[0]=='~'){
	fn2=malloc(strlen(HOME)+strlen(fn)+16);
	strcpy(fn2,HOME);
	strcat(fn2,fn+1);
    }else{
	fn2=malloc(strlen(fn)+16);
	strcpy(fn2,fn);
    }
    /*If there is no recognized suffix, add .bin in the end. */
    if(!check_suffix(fn2,".bin")  && !check_suffix(fn2, ".bin.gz")
       && !check_suffix(fn2,".fits") && !check_suffix(fn2, ".fits.gz")){
	strncat(fn2, ".bin", 4);
    }
    if(mod[0]=='r' || mod[0]=='a'){
	char *fnr=NULL;
	if(!(fnr=search_file(fn2))){/*If does not exist. */
	    if(!check_suffix(fn2, ".gz")){
		/*does not end with .gz, add gz*/
		strncat(fn2, ".gz", 3);
	    }else if (check_suffix(fn, ".gz")){
		/*ended with gz, remove .gz*/
		fn2[strlen(fn2)-3]='\0';
	    }
	    fnr=search_file(fn2);
	}
	free(fn2);
	fn2=fnr;
    }else if (mod[0]=='w'){/*for write, no suffix. we append .bin or .bin.gz */
	if(islink(fn2)){/*remove old file to avoid write over a symbolic link. */
	    if(remove(fn2)){
		error("Failed to remove %s\n", fn2);
	    }
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
}
int zfexist(const char *format, ...){
    format2fn;
    char *fn2=procfn(fn, "rb", 0);
    int ans=0;
    if(fn2){ 
	ans=1;
	free(fn2);
    }
    return ans;
}
void zftouch(const char *format, ...){
    format2fn;
    char *fn2=procfn(fn, "rb", 0);
    if(utimes(fn2, NULL)){
	perror("zftouch failed");
    }
    free(fn2);
}
PNEW(lock);
/**
   Open the file and return a file_t struct that either contains a file
   descriptor (for .bin) or a zlib file pointer.
 */
file_t* zfopen(const char *fn, char *mod){
    LOCK(lock);
    file_t* fp=calloc(1, sizeof(file_t));
    const char* fn2=fp->fn=procfn(fn,mod,1);
    if(!fp->fn){
	error("%s does not exist for read\n", fn);
	_exit(1);
    }
#if IO_TIMMING == 1
    gettimeofday(&(fp->tv1), NULL);
#endif
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
	    error("Trying to write to a file that is already opened for writing: %s\n", fn2);
	}else{
	    if(mod[0]=='w' && ftruncate(fp->fd, 0)){/*Need to manually truncate the file. */
		warning2("Truncating %s failed\n", fn2);
	    }
	}
	break;
    default:
	error("Unknown mod=%s\n", mod);
    }
    if(fp->fd==-1){
	error("Unable to open file %s\n", fn2);
	_exit(1);
    }
    /*check fn instead of fn2. if end of .bin or .fits, disable compressing.*/
    if(mod[0]=='w'){
	if(check_suffix(fn, ".bin") || check_suffix(fn, ".fits")){
	    fp->isgzip=0;
	}else if (check_suffix(fn, ".gz")){
	    fp->isgzip=1;
	}else{
	    fp->isgzip=0;
	}
    }else{ 
	uint16_t magic;
	if(read(fp->fd, &magic, sizeof(uint16_t))!=sizeof(uint16_t)){
	    error("Unable to read %s.\n", fp->fn);
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
	    error("Error gzdopen for %s\n",fn2);
	}
    }else{
	if(!(fp->p=fdopen(fp->fd,mod))){
	    error("Error fdopen for %s\n",fn2);
	}
    }
    if(check_suffix(fn, ".fits") || check_suffix(fn, ".fits.gz")){
	fp->isfits=1;
    }
    if(mod[0]=='w' && !fp->isfits){
	write_timestamp(fp);
    }
    UNLOCK(lock);
    return fp;
}

/**
 * Return the underlining filename
 */
const char *zfname(file_t *fp){
    return fp->fn;
}
/**
   Return whether it's fits file.
*/
int zfisfits(file_t *fp){
    return fp->isfits?1:0;
}
/**
   Close the file.
 */
void zfclose(file_t *fp){
    LOCK(lock);
#if IO_TIMMING == 1
    long fpos;
    if(fp->isgzip){
	fpos=gztell((voidp)fp->p);
    }else{
	fpos=ftell((FILE*)fp->p);
    }
    double size=(double)fpos/1024./1024.;
    struct timeval tv2;
    gettimeofday(&tv2, NULL);
    double tm=(double)tv2.tv_sec+(double)tv2.tv_usec*1e-6
	-((double)fp->tv1.tv_sec+(double)fp->tv1.tv_usec*1e-6);
    const char *type[2]={"\033[33mUncompress", "\033[32mCompressed"};
    info2("%s File %7.2f MiB takes %5.2f s at %6.1f MiB/s: "
	  "%s\033[00m\n", type[fp->isgzip], size, tm, size/tm, fp->fn);
#endif
    if(fp->isgzip){
	gzclose((voidp)fp->p);
    }else{
	fclose((FILE*)fp->p);
    }
    free(fp->fn);
    free(fp);
    UNLOCK(lock);
}
static inline void zfwrite_do(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    if(fp->isgzip){
	if(gzwrite((voidp)fp->p, ptr, size*nmemb)!=size*nmemb){
	    perror("gzwrite");
	    error("write to %s failed\n", fp->fn);
	}
    }else{
	if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
	    perror("fwrite");
	    error("write to %s failed\n", fp->fn);
	}
    }
}
/**
   Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
   fwrite. Follows the interface of fwrite.
 */
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits){
	int length=size*nmemb;
	if(!length) return;
	/* write a block of 2880 bytes each time, with big-endianness.*/
	const int bs=2880;
	char junk[bs];
	int nb=(length+bs-1)/bs;
	char *in=(char*)ptr;
	for(int ib=0; ib<nb; ib++){
	    int nd=length<bs?length:bs;
	    switch(size){
	    case 1:
		memcpy(junk, in, nd);
		break;
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
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread.
 */
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(zfread2(ptr, size, nmemb, fp)){
	error("Error happened while reading %s\n", fp->fn);
    }
}
/**
   Move the current position pointer, like fseek
*/
int zfseek(file_t *fp, long offset, int whence){
    if(fp->isgzip){
	return gzseek((voidp)fp->p,offset,whence)<0?-1:0;
    }else{
	return fseek((FILE*)fp->p,offset,whence);
    }
}
/**
   Move the file position pointer to the beginning
*/
void zfrewind(file_t *fp){
    if(fp->isgzip){
	if(gzrewind((voidp)fp->p)){
	    error("Failed to rewind\n");
	}
    }else{
	rewind((FILE*)fp->p);
    }
}
/**
   Tell position in file.
 */
int zfpos(file_t *fp){
    if(fp->isgzip){
	return gztell((voidp)fp->p);
    }else{
	return ftell((FILE*)fp->p);
    }
}
/**
   Tell whether we end of file is reached.
 */
int zfeof(file_t *fp){
    return zfseek(fp, 1, SEEK_SET)<0?-1:0;
}
/**
   Flush the buffer.
*/
void zflush(file_t *fp){
    if(fp->isgzip){
	gzflush(fp->p,4);
    }else{
	fflush(fp->p);
    }
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
   Obtain the current magic number. If it is a header, read it out if output of
   header is not NULL.  The header will be appended to the output header.

   In file, the header is located before the data magic.
*/
static uint32_t 
read_bin_magic(file_t *fp, char **header){
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
		/*zfseek failed in cygwin (gzseek). so we already readin the header instead.*/
		char header2[nlen];
		zfread(header2, 1, nlen, fp);
		header2[nlen-1]='\0'; /*make sure it is NULL terminated. */
		if(header){
		    if(*header){
			*header=realloc(*header, ((*header)?strlen(*header):0)+strlen(header2)+1);
			strncat(*header, header2, nlen);
		    }else{
			*header=strdup(header2);
		    }
		}
	    }
	    zfread(&nlen2, sizeof(uint64_t),1,fp);
	    zfread(&magic2, sizeof(uint32_t),1,fp);
	    if(magic!=magic2 || nlen!=nlen2){
		info("magic=%u, magic2=%u, nlen=%lu, nlen2=%lu\n", 
		     magic, magic2, (unsigned long)nlen, (unsigned long)nlen2);
		error("Header verification failed\n");
	    }
	}else{ /*otherwise return the magic number*/
	    return magic;
	}
    }/*while*/
}

/**
   Write the magic into file. Also write a dummy header to make data alignment to 8 bytes.
*/
static void 
write_bin_magic(uint32_t magic, file_t *fp){
    if(fp->isfits) error("fits file is not supported\n");
    uint32_t magic2=M_SKIP;
    zfwrite(&magic2, sizeof(uint32_t), 1, fp);
    zfwrite(&magic,  sizeof(uint32_t), 1, fp);
}
/**
   Append the header to current position in the file.  First write magic number, then the
   length of the header, then the header, then the length of the header again,
   then the magic number again. The two set of strlen and header are used to
   identify the header from the end of the file and also to verify that what we
   are reading are indeed header. The header may be written multiple
   times. They will be concatenated when read. The header should contain
   key=value entries just like the configuration files. The entries should be
   separated by new line charactor. 
 */
/*
 * Don't need to write M_SKIP here because we have a pair of magic
 */
static void 
write_bin_header(const char *header, file_t *fp){
    if(!header) return;
    if(fp->isfits) error("fits file is not supported\n");
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(header)+1;
    /*make header 8 byte alignment. */
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

    int bitpix=0;
    switch(magic){
    case M_FLT:
	bitpix=-32;
	break;
    case M_DBL:
	bitpix=-64;
	break;
    default:
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
    if(str){
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
    char line[81];
    int end=0;
    int page=0;
    int bitpix=0;
    int naxis=0;
    while(!end){
	int start=0;
	if(page==0){
	    if(zfread2(line, 1, 80, fp)) return -1; line[80]='\0';
	    if(strncmp(line, "SIMPLE", 6) && strncmp(line, "XTENSION= 'IMAGE", 16)){
		error("Fits header is not recognized\n");
	    }
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &bitpix)!=1) error("Unable to determine bitpix\n");
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &naxis)!=1) error("Unable to determine naxis\n");
	    if(naxis>2) error("Data type not supported\n");
	    if(naxis>0){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (unsigned long *)nx)!=1) error("Unable to determine nx\n");
	    }else{
		*nx=0;
	    }
	    if(naxis>1){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (unsigned long *)ny)!=1) error("Unable to determine ny\n");
	    }else{
		*ny=0;
	    }
	    start=3+naxis;
	}
	for(int i=start; i<36; i++){
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(!strncmp(line, "COMMENT", 7)){
		for(int j=79; j>9; j--){
		    if(isspace((int)line[j])){
			line[j]='\0';
		    }else{
			break;
		    }
		}
		int length=strlen(line+10);
		if(*str){
		    *str=realloc(*str, strlen(*str)+length+1);
		}else{
		    *str=malloc(length+1); (*str)[0]='\0';
		}
		strcat(*str, line+10);
	    }else if(!strncmp(line, "END",3)){
		end=1;
	    }
	}
    }
    switch(bitpix){
    case -32:
	*magic=M_FLT;
	break;
    case -64:
	*magic=M_DBL;
	break;
    default:
	error("Invalid\n");
    }
    return 0;
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
	    ans=-1;
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
	error("read_header failed for %s\n", fp->fn);
    }
}
/**
 * Get the length of string in header, rounded to multiple of 8.
 */
uint64_t bytes_header(const char *header){
    if(header){
	uint64_t len=strlen(header)+1;
	if(len%8 != 0){
	    len=(len/8+1)*8;
	}
	return len+24;
    }else{
	return 0;
    }
}
/**
   Write the time stamp as header into current location in the file.
*/
void write_timestamp(file_t *fp){
    char header[128];
    snprintf(header,128, "Created by MAOS Version %s on %s in %s\n",
	     PACKAGE_VERSION, myasctime(), myhostname());
    write_bin_header(header, fp);
}

/**
   Search and return the value correspond to key. NULL if not found. Do not free the
   returned pointer. The key must be preceeded by space (isspace), and succeeded by = sign. */
const char *search_header(const char *header, const char *key){
    if(!header) return NULL;
    const char *ans=NULL;
    const char *val=header;
    while(val[0]!='\0' && (val=strstr(val, key))){
	if(val>header){
	  char prev=*(val-1);
	  if(!isspace((int)prev) && prev!=';' && prev !=','){
		val=val+strlen(key);
		continue;/*Invalid */
	    }
	}
	val=val+strlen(key);
	while(val[0]==' ') val++;
	if(val[0] == '='){
	    ans=val+1;
	    break;
	}
    }
    return ans;
}
/**
   Read a number from the header with key
*/
double search_header_num(const char *header, const char *key){
    if(!header) return 0;
    const char *val=search_header(header, key);
    if(val){
	return readstr_num(val, NULL);
    }else{
	return NAN;/*not found. */
    }
}


/**
   Write an 1-d or 2-d array into the file. First write a magic number that
   represents the data type. Then write two numbers representing the
   dimension. Finally write the data itself.
 */
void do_write(const void *fpn,     /**<[in] The file pointer*/
	      const int isfn,      /**<[in] Is this a filename or already opened file*/
	      const size_t size,   /**<[in] Size of each element*/
	      const uint32_t magic,/**<[in] The magic number. see bin.h*/ 
	      const char *str,     /**<[in] The header as string*/
	      const void *p,       /**<[in] The data of the array*/
	      const uint64_t nx,   /**<[in] Number of rows. this index changes fastest*/
	      const uint64_t ny    /**<[in] Number of columns. 1 for vector*/
	      ){
    file_t *fp;
    header_t header={magic, nx, ny, (char*)str};
    if(isfn){
	fp=zfopen((char*)fpn, "wb");
    }else{
	fp=(file_t*) fpn;
    }
    write_header(&header, fp);
    if(nx*ny>0) zfwrite(p, size, nx*ny, fp);
    if(isfn) zfclose(fp);
}
/**
   Write a double array of size nx*ny to file.
*/
void writedbl(const double *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(double), M_DBL, NULL, p, nx, ny);
}
/**
   Write a double array of size nx*ny to file.
*/
void writeflt(const float *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(float), M_FLT, NULL, p, nx, ny);
}
/**
   Write a double complex array of size nx*ny to file.
*/
void writecmp(const dcomplex *p, long nx,long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(dcomplex), M_CMP, NULL, p, nx, ny);
}
/**
   Write a float complex array of size nx*ny to file.
*/
void writefcmp(const fcomplex *p, long nx,long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(fcomplex), M_ZMP, NULL, p, nx, ny);
}
/**
   Write a int array of size nx*ny to file.
 */
void writeint(const int *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(int), M_INT32, NULL, p, nx, ny);
}
/**
   Write a long array of size nx*ny to file.
 */
void writelong(const long *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(long), sizeof(long)==8?M_INT64:M_INT32, NULL, p, nx, ny);
}
/**
   Write spint array of size nx*ny to file. 
 */
void writespint(const spint *p, long nx, long ny, const char *format,...){
    format2fn;
    do_write(fn, 1, sizeof(spint), M_SPINT, NULL, p, nx, ny);
}
void readspintdata(file_t *fp, uint32_t magic, spint *out, long len){
    uint32_t size=0;
    switch(magic){
    case M_INT64:
	size=8;
	break;
    case M_INT32:
	size=4;
	break;
    case M_DBL:/*saved by matlab. */
	size=-8;
	break;
    default:
	error("This is not a valid sparse spint file. magic=%x\n", magic);
    }
    if(sizeof(spint)==size){/*Matched int. */
	zfread(out, sizeof(spint), len, fp);
    }else{
	size=abs(size);
	void *p=malloc(size*len);
	zfread(p, size, len, fp);
	switch(magic){
	case M_INT64:{
	    uint64_t *p2=p;
	    for(unsigned long j=0; j<len; j++){
		out[j]=(spint)p2[j];
	    }
	}
	    break;
	case M_INT32:{
	    uint32_t *p2=p;
	    for(unsigned long j=0; j<len; j++){
		out[j]=(spint)p2[j];
	    }
	}
	    break;
	case M_DBL:{
	    double *p2=p;
	    for(unsigned long j=0; j<len; j++){
		out[j]=(spint)p2[j];
	    }
	}
	    break;
	}
    }
}
/**
   Read spint array of size nx*ny from file. Optionally convert from other formats.
 */
spint *readspint(file_t *fp, long* nx, long* ny){
    uint32_t magic=read_bin_magic(fp, NULL);
    uint64_t nx2, ny2;
    zfreadlarr(fp, 2, &nx2, &ny2);
    spint *out=NULL;
    if(nx!=0 && ny!=0){
	*nx=(long)nx2;
	*ny=(long)ny2;
	out=malloc(nx2*ny2*sizeof(spint));
	readspintdata(fp, magic, out, nx2*ny2);
    }else{
	*nx=0;
	*ny=0;
    }
    return out;
}
/**
 * Unreference the mmaped memory. When the reference drops to zero, unmap it.
 */
void mmap_unref(struct mmap_t *in){
    if(in->nref>1){
	in->nref--;
	/*info("%p: nref decreased to %ld\n", in->p, in->nref); */
    }else{
	/*info("%p: nref is 1, unmap\n", in->p); */
	munmap(in->p, in->n);
	if(in->fd!=-1){
	    close(in->fd);
	}
	free(in);
    }
}
/**
   Create a mmap_t object.
*/
struct mmap_t *mmap_new(int fd, void *p, long n){
    struct mmap_t *out=calloc(1, sizeof(struct mmap_t));
    out->p=p;
    out->n=n;
    out->nref=1;
    out->fd=fd;
    return out;
}
/**
   Add a reference to a mmap_t.
 */
mmap_t*mmap_ref(mmap_t *in){
    in->nref++;
    return in;
}
