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

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "common.h"
#include "common.h"
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

  Notice that the stored data by 64 machines using type long
  can not be read by 32 bit due to different byte
  length. Modified the routine to first 64 bit read/write.
  nx,ny, are uint64_t for 64 and 32 bit machines
  magic is uint32_t always to be backward compatible.
*/

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
    if(mod[0]=='r'){
	char *fnr=NULL;
	if(!(fnr=search_file(fn2))){
	    if(check_suffix(fn2, ".bin.gz")){
		//ended with bin.gz
		fn2[strlen(fn2)-3]='\0';
		if(!(fnr=search_file(fn2))){
		    error("Neither %s nor %s exist\n", fn,fn2);
		}
	    }else if(check_suffix(fn2, ".bin")){
		//ended with .bin
		strncat(fn2, ".gz", 3);
		if(!(fnr=search_file(fn2))){
		    error("Neither %s nor %s exist\n", fn,fn2);
		}
	    }else{//no recognized suffix
		strncat(fn2, ".bin", 4);
		if(!(fnr=search_file(fn2))){
		    strncat(fn2,".gz",3);
		    if(!(fnr=search_file(fn2))){
			error("Neither %s, %s.bin, nor %s.bin.gz exist\n",
			      fn,fn,fn);
		    }
		}
	    }
	}
	free(fn2);
	fn2=fnr;
    }else if (mod[0]=='w'){//for write, no suffix. we append .bin or .bin.gz
	if(!(check_suffix(fn2, ".bin.gz") || check_suffix(fn2, ".bin"))){
	    if(defaultgzip){
		strncat(fn2,".bin.gz",7);
	    }else{
		strncat(fn2,".bin",4);
	    }
	}
	if(exist(fn2)){//remove old file to avoid write over a symbolic link.
	    remove(fn2);
	}
    }else if(mod[0]=='a'){//open for appending
	if(!(check_suffix(fn2, ".bin.gz") || check_suffix(fn2, ".bin"))){
	    if(defaultgzip){
		strncat(fn2,".bin.gz",7);
	    }else{
		strncat(fn2,".bin",4);
	    }
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
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
    const char* fn2=fp->fn;
#if IO_TIMMING == 1
    gettimeofday(&(fp->tv1), NULL);
#endif
    if(check_suffix(fn2, ".gz")){
	fp->isgzip=1;
	if(!(fp->p=gzopen(fn2,mod))){
	    error("Error gzopen for %s\n",fn2);
	}
    }else{ 
	fp->isgzip=0;
	if(!(fp->p=fopen(fn2,mod))){
	    error("Error fopen for %s\n",fn2);
	}
    }
    if(mod[0]=='w'){
	write_timestamp(fp);
    }
    UNLOCK(lock);
    return fp;
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
header is not NULL.  The header will be appended to the output header.*/
uint32_t read_magic(file_t *fp, char **header){
    uint32_t magic,magic2;
    uint64_t nlentot=0;
    uint64_t nlen, nlen2;
    if(header && *header){
	nlentot=strlen(*header)+1;
    }
    while(1){
	/*read the magic number.*/
	zfread(&magic, sizeof(uint32_t), 1, fp);
	/*If it is header, read or skip it.*/
	if(magic==M_HEADER){
	    zfread(&nlen, sizeof(uint64_t), 1, fp);
	    if(nlen>0){
		if(header){
		    *header=realloc(*header, nlentot+nlen);
		    if(nlentot>0){
			(*header)[nlentot-1]='\n';
		    }
		    zfread(*header+nlentot, 1, nlen, fp);
		    nlentot+=nlen;
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
   Append the header to current position in the file.  First write magic number, then the
   length of the header, then the header, then the length of the header again,
   then the magic number again. The two set of strlen and header are used to
   identify the header from the end of the file and also to verify that what we
   are reading are indeed header. The header may be written multiple
   times. They will be concatenated when read. The header should contain
   key=value entries just like the configuration files. The entries should be
   separated by new line charactor. */
void write_header(const char *header, file_t *fp){
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(header)+1;
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(header, 1, nlen, fp);
    zfwrite(&nlen, sizeof(uint64_t), 1, fp);
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
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
   Append to the file the header to end of an already saved file.
*/
void write_header_end(const char *header,const char *format,...){
    format2fn;
    file_t *fp=zfopen((char*)fn, "ab");//open for append.
    if(zfseek(fp, 0, SEEK_END)){
	error("Failed to seek to the end of file");
    }else{
	write_header(header, fp);
    }
}

/**
   Search and return the value correspond to key. NULL if not found. Do not free the
   returned pointer.  */
const char *search_header(const char *header, const char *key){
    char *val=strstr(header, key);
    if(val){
	char *equal=strchr(val,'=');
	if(equal){
	    return equal+1;
	}else{
	    error("= not found after key in %s\n", header);
	    return NULL;
	}
    }else{
	return NULL;
    }
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
	return 0;//not found.
    }
}
/**
   Read header from the end of the file. Difference headers will be combined. Do
not use. Seek back is very slow in big file.  */
char *read_header_end(const char *format,...){
    format2fn;
    file_t *fp=zfopen((char*)fn, "rb");
    zfseek(fp, 0, SEEK_END);//move to the end;
    char *strout=NULL;
    int nlenout=0;
    uint32_t magic, magic2;

    while(1){
        if(zfseek(fp, -sizeof(uint32_t), SEEK_CUR)){//move the the end
	    warning("Error seek to (another) header\n");
	    break;
	}
	zfread(&magic, sizeof(uint32_t), 1, fp);
	if(magic!=M_HEADER){
	    warning("There is no more header in file %s\n",fp->fn);
	    break;
	}
	//move 8 byte back to get the string length
	zfseek(fp, -(sizeof(uint64_t)+sizeof(uint32_t)), SEEK_CUR);
	uint64_t nlen, nlen2;
	//after read, pointer is after the nlen.
	zfread(&nlen, sizeof(uint64_t), 1, fp);
	//seek to the beginning of the first header
	zfseek(fp, -(nlen+sizeof(uint64_t)*2+sizeof(uint32_t)), SEEK_CUR);
	zfread(&magic2, sizeof(uint32_t), 1, fp);
	zfread(&nlen2, sizeof(uint64_t), 1, fp);
	if(magic!=magic2 || nlen!=nlen2){
	    warning("This is not a header\n");
	    break;
	}
	if(nlen>0){
	    //enlarge the string.
	    info("Reading header with length %lu\n",nlen);
	    if(nlenout>0){
		strout[nlenout-1]='\n';//turn \0 to \n.
	    }
	    strout=realloc(strout, nlenout+nlen);
	    zfread(strout+nlenout, 1, nlen, fp);
	    strout[nlenout-1]='\0';
	    nlenout+=nlen;
	}
	//now seek to before our first header, in case there is another header.
	if(zfseek(fp, -(nlen+sizeof(uint64_t)+sizeof(uint32_t)), SEEK_CUR)){
	    warning("This file contains only headers\n");
	    break;
	}
    }
    zfclose(fp);
    return strout;
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
	      const int isfile,    /**<[in] Is this a filename or already opened file*/
	      const size_t size,   /**<[in] Size of each element*/
	      const uint32_t magic,/**<[in] The magic number. see bin.h*/ 
	      const void *p,       /**<[in] The data of the array*/
	      const uint64_t nx,   /**<[in] Number of rows. this index changes fastest*/
	      const uint64_t ny    /**<[in] Number of columns. 1 for vector*/
	      ){
    file_t *fp;
    if(isfile){
	fp=zfopen((char*)fpn, "wb");
    }else{
	fp=(file_t*) fpn;
    }
    zfwrite(&magic, sizeof(uint32_t),1,fp);
    if(p && nx>0 && ny>0){
	zfwritelarr(fp, 2, &nx, &ny);
	zfwrite(p, size, nx*ny,fp);
    }else{
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }
    if(isfile) zfclose(fp);
}
/**
   Read an 1-d or 2-d array from the file. First a magic number is read from the
   file and compared with magicin. If mismatch, operation will fail. Then the
   dimension and actual data are read and output to pnx, pny, p
 */
void do_read(const void *fpn,      /**<[in]  The file pointer*/
	     const int isfile,     /**<[in]  Is this a filename or already opened file*/
	     const size_t size,    /**<[in]  Size of each element*/
	     const uint32_t magic, /**<[in] The magic number wanted*/
	     void **p,             /**<[out] The address of the pointer of the array*/
	     uint64_t *pnx,        /**<[out] Return number of rows*/
	     uint64_t *pny         /**<[out] Return number of columns*/
	     ){
    file_t *fp;
    if(isfile){
	fp=zfopen((char*)fpn, "rb");
    }else{
	fp=(file_t*) fpn;
    }
    uint32_t magic2;
    zfread(&magic2, sizeof(uint32_t),1,fp);
    if(magic2!=magic){
	error("File is wrong format. magic is %x, wanted %x\n",magic2, magic);
    }
    zfreadlarr(fp, 2, pnx, pny);
    if(*pnx==0 || *pny==0){
	*p=NULL;
    }else{
	if(*p) warning("*p should be null\n");
	*p=malloc(size*(*pnx)*(*pny));
	zfread(*p, size, (*pnx)*(*pny),fp);
    }
    if(isfile) zfclose(fp);
}
/**
   Write a double array of size nx*ny to file.
*/
void writedbl(const double *p,long nx, long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(double), M_DBL, p, (uint64_t)nx, (uint64_t)ny);
}
/**
   Write a double complex array of size nx*ny to file.
*/
void writecmp(const dcomplex *p, long nx,long ny, const char*format,...){
    format2fn;
    do_write(fn, 1, sizeof(dcomplex), M_CMP, p, (uint64_t)nx, (uint64_t)ny);
}
/**
   Write a 64 bit integer array of size nx*ny to file.
*/
void writeint64(const int64_t *p, long nx, long ny, const char *format,...){
    format2fn;
    do_write(fn, 1, sizeof(int64_t), M_INT64, p, (uint64_t)nx, (uint64_t)ny);
}
/**
   Write a 32 bit integer array of size nx*ny to file.
*/
void writeint32(const int32_t *p, long nx, long ny, const char *format,...){
    format2fn;
    do_write(fn, 1, sizeof(int32_t), M_INT32, p, (uint64_t)nx, (uint64_t)ny);
}
/**
   Read a double array of size nx*ny from file
*/
void readdbl(double **p, long *nx, long *ny, const char *format,...){
    format2fn;
    uint64_t nx2,ny2;
    do_read(fn, 1, sizeof(double), M_DBL, (void**)p, &nx2, &ny2);
    *nx=nx2;
    *ny=ny2;
    if(*nx!=nx2 || *ny!=ny2){
	warning("Data overflow\n");
    }
}

/**
   Read a 64 bit integer array of size nx*ny from file.
*/
void readint64(int64_t **p, long *nx, long *ny, 
		 const char *format,...){
    format2fn;
    uint64_t nx2,ny2;
    do_read(fn, 1, sizeof(uint64_t), M_INT64, (void**)p, &nx2, &ny2);
    *nx=nx2;
    *ny=ny2;
    if(*nx!=nx2 || *ny!=ny2){
	warning("Data overflow\n");
    }
}

/**
   Read a 32 bit integer array of size nx*ny from file.
*/
void readint32(int32_t **p, long *nx, long *ny, 
		 const char *format,...){
    format2fn;
    uint64_t nx2,ny2;
    do_read(fn, 1,  sizeof(uint32_t), M_INT32, (void**)p, &nx2, &ny2);
    *nx=nx2;
    *ny=ny2;
    if(*nx!=nx2 || *ny!=ny2){
	warning("Data overflow\n");
    }
}
