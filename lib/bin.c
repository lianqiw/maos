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
#include <ctype.h> //isspace
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
    int isgzip;
#if USE_PTHREAD > 0
    pthread_mutex_t lock;
#endif
    char *fn;
    int fd;//For locking.
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
    //If there is no recognized suffix, add .bin in the end.
    if(!check_suffix(fn2,".bin") && !check_suffix(fn2, ".bin.gz")){
	strncat(fn2, ".bin", 4);
    }
    if(mod[0]=='r' || mod[0]=='a'){
	char *fnr=NULL;
	if(!(fnr=search_file(fn2))){//If does not exist.
	    if(check_suffix(fn2, ".bin")){
		//ended with .bin, change to .bin.gz
		strncat(fn2, ".gz", 3);
		if(!(fnr=search_file(fn2))){
		    return NULL;
		}
	    }else{
		//ended with bin.gz, change to .bin
		fn2[strlen(fn2)-3]='\0';
		if(!(fnr=search_file(fn2))){
		    return NULL;
		}
	    }
	}
	free(fn2);
	fn2=fnr;
    }else if (mod[0]=='w'){//for write, no suffix. we append .bin or .bin.gz
	if(islink(fn2)){//remove old file to avoid write over a symbolic link.
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
    PINIT(fp->lock);
    fp->fn=procfn(fn,mod,1);
    if(!fp->fn){
	error("%s does not exist for read\n", fn);
    }
    const char* fn2=fp->fn;
#if IO_TIMMING == 1
    gettimeofday(&(fp->tv1), NULL);
#endif
    /*Now open the file to get a fd number that we can use to lock on the
      file.*/
    switch(mod[0]){
    case 'r'://read only
	fp->fd=open(fn2, O_RDONLY);
	break;
    case 'w'://write
    case 'a':
	fp->fd=open(fn2, O_RDWR | O_CREAT, 0600);
	if(fp->fd!=-1 && flock(fp->fd, LOCK_EX|LOCK_NB)){
	    error("Trying to write to a file that is already opened for writing: %s\n", fn2);
	}else{
	    if(mod[0]=='w' && ftruncate(fp->fd, 0)){//Need to manually truncate the file.
		warning2("Truncating %s failed\n", fn2);
	    }
	}
	break;
    default:
	error("Unknown mod=%s\n", mod);
    }
    if(fp->fd==-1){
	error("Unable to open file %s\n", fn2);
    }
    if(check_suffix(fn2, ".bin") && mod[0]=='w'){
	fp->isgzip=0;
	if(!(fp->p=fdopen(fp->fd,mod))){
	    error("Error fdopen for %s\n",fn2);
	}
    }else{ 
	fp->isgzip=1;
	if(!(fp->p=gzdopen(fp->fd,mod))){
	    error("Error gzdopen for %s\n",fn2);
	}
    }
    if(mod[0]=='w'){
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
   Close the file.
 */
void zfclose(file_t *fp){
    LOCK(lock);
    PDEINIT(fp->lock);
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
    fprintf(stderr, "%s File %7.2f MiB takes %5.2f s at %6.1f MiB/s: "
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
/**
   Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
   fwrite. Follows the interface of fwrite.
 */
void zfwrite(const void* ptr, const size_t size, 
	     const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    LOCK(fp->lock);
    if(fp->isgzip){
	gzwrite((voidp)fp->p, ptr, size*nmemb);
    }else{
	size_t ct=fwrite(ptr, size, nmemb, (FILE*)fp->p);
	if(ct!=nmemb)
	    error_write;
    }
    UNLOCK(fp->lock);
}
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread.
 */
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    LOCK(fp->lock);
    if(fp->isgzip){
	int status;
	if((status=gzread((voidp)fp->p, ptr, size*nmemb))<1){
	    if(status==-1){
		warning("Error happened in reading\n");
	    }else{
		warning("End of File encoutnered\n");
	    }
	}
    }else{
	size_t ct=fread(ptr, size, nmemb, (FILE*)fp->p);
	if(ct!=nmemb){
	    error("Reading from %s failed\n", fp->fn);
	}
	if(feof((FILE*)fp->p)){
	    warning("End of File encountered in %s\n",fp->fn);
	}
    }
    UNLOCK(fp->lock);
}
/**
   Move the current position pointer, like fseek
*/
int zfseek(file_t *fp, long offset, int whence){
    if(fp->isgzip){
	int res=gzseek((voidp)fp->p,offset,whence);
	if(res<0)
	    return 1;
	else
	    return 0;
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
   Tell whether we end of file is reached.
 */
int zfeof(file_t *fp){
    LOCK(fp->lock);
    int ans;
    long cpos, fpos;
    if(fp->isgzip){
	cpos=gztell((voidp)fp->p);
	gzseek((voidp)fp->p,0,SEEK_END);
	fpos=gztell((voidp)fp->p);
    }else{
	cpos=ftell((FILE*)fp->p);
	fseek((FILE*)fp->p, 0, SEEK_END);    
	fpos=ftell((FILE*)fp->p);
    }
    if(cpos!=fpos){
	fprintf(stderr,"There are unread bytes near the end of the file\n"
		"Current position is %ld, File END is %ld\n",cpos, fpos);
	ans=1;
    }else{
	ans=0;
    }
    UNLOCK(fp->lock);
    return ans;
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
   Obtain the current magic number. If it is a header, read it out if output of
   header is not NULL.  The header will be appended to the output header.

   In file, the header is located before the data magic.
*/
uint32_t read_magic(file_t *fp, char **header){
    uint32_t magic,magic2;
    uint64_t nlen, nlen2;
    while(1){
	/*read the magic number.*/
	zfread(&magic, sizeof(uint32_t), 1, fp);
	/*If it is header, read or skip it.*/
	if(magic==M_SKIP){
	    continue;
	}else if(magic==M_HEADER){
	    zfread(&nlen, sizeof(uint64_t), 1, fp);
	    if(nlen>0){
		if(header){
		    char header2[nlen];
		    zfread(header2, 1, nlen, fp);
		    header2[nlen-1]='\0'; //make sure it is NULL terminated.
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
void write_magic(uint32_t magic, file_t *fp){
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
void write_header(const char *header, file_t *fp){
    if(!header) return;
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(header)+1;
    //make header 8 byte alignment.
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
    write_header(header, fp);
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
	    if(!isspace(val[-1]) && val[-1]!=';' && val[-1] !=','){
		val=val+strlen(key);
		continue;//Invalid
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
	return NAN;//not found.
    }
}
/**
   Write multiple long numbers into the file. To write three numbers, a, b, c,
call with zfwritelarr(fp, 3, &a, &b, &c); */
void zfwritelarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              //Initialize the argument list.
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  //Get the next argument value.  
	zfwrite(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       // Clean up. 
}
/**
   Read multiple long numbers from the file. To read three numbers, a, b, c,
   call with zfreadlarr(fp, 3, &a, &b, &c);
 */
void zfreadlarr(file_t *fp, int count, ...){
    va_list ap;
    int i;
    va_start (ap, count);              //Initialize the argument list.
    for (i = 0; i < count; i++){
	uint64_t *addr=va_arg (ap, uint64_t*);  //Get the next argument value.  
	zfread(addr, sizeof(uint64_t), 1, fp);
    }
    va_end (ap);                       // Clean up. 
}
/**
   Write an 1-d or 2-d array into the file. First write a magic number that
   represents the data type. Then write two numbers representing the
   dimension. Finally write the data itself.
 */
void do_write(const void *fpn,     /**<[in] The file pointer*/
	      const int isfn,    /**<[in] Is this a filename or already opened file*/
	      const size_t size,   /**<[in] Size of each element*/
	      const uint32_t magic,/**<[in] The magic number. see bin.h*/ 
	      const void *p,       /**<[in] The data of the array*/
	      const uint64_t nx,   /**<[in] Number of rows. this index changes fastest*/
	      const uint64_t ny    /**<[in] Number of columns. 1 for vector*/
	      ){
    file_t *fp;
    if(isfn){
	fp=zfopen((char*)fpn, "wb");
    }else{
	fp=(file_t*) fpn;
    }
    //Write a dummy magic so that our header is multiple of 8 bytes long.
    write_magic(magic, fp);
    if(p && nx>0 && ny>0){
	zfwritelarr(fp, 2, &nx, &ny);
	zfwrite(p, size, nx*ny,fp);
    }else{
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }
    if(isfn) zfclose(fp);
}

/**
   Write a double array of size nx*ny to file.
*/
void writedbl(const double *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(double), M_DBL, p, nx, ny);
}
/**
   Write a double complex array of size nx*ny to file.
*/
void writecmp(const dcomplex *p, long nx,long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(dcomplex), M_CMP, p, nx, ny);
}
/**
   Write a int array of size nx*ny to file.
 */
void writeint(const int *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(int), M_INT32, p, nx, ny);
}
/**
   Write a long array of size nx*ny to file.
 */
void writelong(const long *p, long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(long), sizeof(long)==8?M_INT64:M_INT32, p, nx, ny);
}
/**
   Write spint array of size nx*ny to file. 
 */
void writespint(const spint *p, long nx, long ny, const char *format,...){
    format2fn;
    do_write(fn, 1, sizeof(spint), M_SPINT, p, nx, ny);
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
    case M_DBL://saved by matlab.
	size=-8;
	break;
    default:
	error("This is not a valid sparse spint file. magic=%x\n", magic);
    }
    if(sizeof(spint)==size){//Matched int.
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
    uint32_t magic=read_magic(fp, NULL);
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
	//info("%p: nref decreased to %ld\n", in->p, in->nref);
    }else{
	//info("%p: nref is 1, unmap\n", in->p);
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
