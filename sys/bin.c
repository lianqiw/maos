/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdint.h>
#include <limits.h>
#include <ctype.h> /*isspace */
#include "process.h"
#include "misc.h"
#include "path.h"
#include "thread.h"
#include "bin.h"

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
/*
  contains information about opened files.
*/
struct file_t{
    void *p;   /**<FILE* or voidp when gzipped*/
    int isgzip;/**<Is the file zipped.*/
    int isfits;/**<Is the file fits.*/
    char *fn;  /**<The disk file name*/
    int fd;    /**<The underlying file descriptor, used for locking. */
};
/*
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
   Disables file saving when set to 1.
*/
int disable_save=0;

/*
  Process the input file name and return file names that can be open to
  read/write. If the file name does not end with .bin or .bin.gz it will add to
  the end .bin or .bin.gz depending on the value of defaultgzip. For read only
  access, it will also look into the path for files.
*/
static char* procfn(const char *fn, const char *mod, const int defaultgzip){
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
    if(mod[0]=='r'){
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
	if(!fnr){
	  
	}
    }else if (mod[0]=='w' || mod[0]=='a'){
	if(disable_save){//When saving is disabled, allow writing to cache folder.
	    char fncache[PATH_MAX];
	    snprintf(fncache, PATH_MAX, "%s/.aos/cache", HOME);
	    if(mystrcmp(fn2, fncache)){
		warning("Saving is disabled for %s.\n", fn2);
		free(fn2);
		fn2=0;
	    }
	}
	if(fn2 && islink(fn2)){/*remove old file to avoid write over a symbolic link. */
	    if(remove(fn2)){
		error("Failed to remove %s\n", fn2);
	    }
	}
    }else{
	error("Invalid mode\n");
    }
    return fn2;
}
/**
   Test whether a bin file exist.
*/
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
/**
   Update the modification time of a bin file.
*/
void zftouch(const char *format, ...){
    format2fn;
    char *fn2=procfn(fn, "rb", 0);
    if(fn2 && utimes(fn2, NULL)){
	perror("zftouch failed");
    }
    free(fn2);
}
PNEW(lock);
/*
  Open a bin file from a fd that may be a socket.
*/
static file_t* zfdopen(int sock, const char *mod){
    file_t* fp=calloc(1, sizeof(file_t));
    fp->isgzip=0;
    fp->fd=sock;
    if(fp->isgzip){
	if(!(fp->p=gzdopen(fp->fd,mod))){
	    error("Error gzdopen for %d\n",sock);
	}
    }else{
	if(!(fp->p=fdopen(fp->fd,mod))){
	    error("Error fdopen for %d\n",sock);
	}
    }
    return fp;
}

/**
   Open a bin file for read/write access. 

   Whether the file is gzipped or not is automatically determined when open for
   reading. When open for writing, if the file name ends with .bin or .fits,
   will not be gzipped. If the file has no suffix, the file will be gzipped and
   .bin is appened to file name.
*/
file_t* zfopen_try(const char *fn, const char *mod){
    LOCK(lock);
    file_t* fp=calloc(1, sizeof(file_t));
    const char* fn2=fp->fn=procfn(fn,mod,1);
    if(!fn2){
	if(mod[0]=='r'){
	    error("%s does not exist for read\n", fn);
	}
	free(fp); fp=0;
	goto end;
    }
    /*Now open the file to get a fd number that we can use to lock on the
      file.*/
    switch(mod[0]){
    case 'r':/*read only */
	if((fp->fd=open(fn2, O_RDONLY))==-1){
	    perror("open for read");
	}else{
	    futimes(fp->fd, NULL);
	}
	break;
    case 'w':/*write */
    case 'a':
	if((fp->fd=open(fn2, O_RDWR | O_CREAT, 0666))==-1){
	    perror("open for write");
	}else{
	    if(flock(fp->fd, LOCK_EX|LOCK_NB)){
		perror("flock");
		close(fp->fd);
		fp->fd=-1;
	    }else if(mod[0]=='w' && ftruncate(fp->fd, 0)){/*Need to manually truncate the file. */
		perror("ftruncate");
	    }
	}
	break;
    default:
	error("Unknown mod=%s\n", mod);
    }
    if(fp->fd==-1){
	error("Unable to open file %s for %s\n", fn2, mod[0]=='r'?"reading":"writing");
	free(fp->fn);
	free(fp);
	fp=0;
	goto end;
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
    /*if(mod[0]=='w' && !fp->isfits){
	write_timestamp(fp);
	}*/
  end:
    UNLOCK(lock);
    return fp;
}
file_t *zfopen(const char *fn, const char *mod){
    file_t *fp=zfopen_try(fn, mod);
    if(!fp){
	error("Open file %s for %s failed\n", fn, mod[0]=='r'?"read":"write");
    }
    return fp;
}
/**
   Return the underlining filename
*/
const char *zfname(file_t *fp){
    return fp->fn;
}
/**
   Return 1 when it is a fits file.
*/
int zfisfits(file_t *fp){
    return fp->isfits?1:0;
}
/**
   Close the file.
*/
void zfclose(file_t *fp){
    LOCK(lock);
    if(fp->isgzip){
	gzclose((voidp)fp->p);
    }else{
	fclose((FILE*)fp->p);
    }
    free(fp->fn);
    free(fp);
    UNLOCK(lock);
}
/*
  Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
  fwrite. Follows the interface of fwrite.
*/
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
   Handles byteswapping in fits file format then call zfwrite_do to do the actual writing.
*/
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp){
    /*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
    if(fp->isfits && BIGENDIAN==0){
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
/**
   Read from the file. If in gzip mode, calls gzread, otherwise, calls
   fread. Follows the interface of fread.
*/
static inline int zfread_do(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(fp->isgzip){
	return gzread((voidp)fp->p, ptr, size*nmemb)>0?0:-1;
    }else{
	return fread(ptr, size, nmemb, (FILE*)fp->p)==nmemb?0:-1;
    }
}
/**
   Handles byteswapping in fits file format then call zfread_do to do the actual writing.
*/
int zfread_try(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
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
   Wraps zfread_try and do error checking.
*/
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
    if(zfread_try(ptr, size, nmemb, fp)){
	perror("zfread");
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
   Tell position pointer in file.
*/
int zfpos(file_t *fp){
    if(fp->isgzip){
	return gztell((voidp)fp->p);
    }else{
	return ftell((FILE*)fp->p);
    }
}
/**
   Return 1 if end of file is reached. 
*/
int zfeof(file_t *fp){
    return zfseek(fp, 1, SEEK_CUR)<0?1:0;
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

/*
  Obtain the current magic number, string header, and array size from a bin
  format file. In file, the string header is located before the data magic.
*/
static int
read_bin_header(header_t *header, file_t *fp){
    uint32_t magic,magic2;
    uint64_t nlen, nlen2;
    if(fp->isfits) error("fits file is not supported\n");
    while(1){
	/*read the magic number.*/
	if(zfread_try(&magic, sizeof(uint32_t), 1, fp)) return -1;
	/*If it is hstr, read or skip it.*/
	if((magic&M_SKIP)==M_SKIP){
	    continue;
	}else if(magic==M_HEADER){
	    zfread(&nlen, sizeof(uint64_t), 1, fp);
	    if(nlen>0){
		/*zfseek failed in cygwin (gzseek). so we always read in the hstr instead.*/
		char hstr2[nlen];
		zfread(hstr2, 1, nlen, fp);
		hstr2[nlen-1]='\0'; /*make sure it is NULL terminated. */
		if(header->str){
		    header->str=realloc(header->str,((header->str)?strlen(header->str):0)+strlen(hstr2)+1);
		    strncat(header->str, hstr2, nlen);
		}else{
		    header->str=strdup(hstr2);
		}
	    }
	    zfread(&nlen2, sizeof(uint64_t),1,fp);
	    zfread(&magic2, sizeof(uint32_t),1,fp);
	    if(magic!=magic2 || nlen!=nlen2){
		info("magic=%u, magic2=%u, nlen=%lu, nlen2=%lu\n", 
		     magic, magic2, (unsigned long)nlen, (unsigned long)nlen2);
		error("Header string verification failed\n");
	    }
	}else{ //Finish
	    header->magic=magic;
	    zfread(&header->nx, sizeof(uint64_t), 1, fp);
	    zfread(&header->ny, sizeof(uint64_t), 1, fp);
	    return 0;
	}
    }/*while*/
    return -1;
}

/*
  Append the header to current position in the file.  
   
  First write magic number, then the length of the header, then the header,
  then the length of the header again, then the magic number again. The two set
  of strlen and header are used to identify the header from the end of the file
  and also to verify that what we are reading are indeed header. The header may
  be written multiple times. They will be concatenated when read. The header
  should contain key=value entries just like the configuration files. The
  entries should be separated by new line charactor.
  
  Don't need to write M_SKIP here because we have a pair of magic. The length
  of header is rounded to multiple of 8 bytes.
*/
static void 
write_bin_headerstr(const char *str, file_t *fp){
    if(!str || !strlen(str)) return;
    assert(!fp->isfits);
    uint32_t magic=M_HEADER;
    uint64_t nlen=strlen(str)+1;
    uint64_t nlen2=(nlen%8)?((nlen/8+1)*8):(nlen);
    char zero[8]={0};    
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    zfwrite(&nlen2, sizeof(uint64_t), 1, fp);
    /*make str 8 byte alignment. */
    zfwrite(str, 1, nlen, fp);
    if(nlen2>nlen)
	zfwrite(zero, 1, nlen2-nlen, fp);
    zfwrite(&nlen2, sizeof(uint64_t), 1, fp);
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
}

/*
  Write fits header. str is extra header that will be put in fits comment
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
    switch(magic & 0xFFFF){
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
	error("Data type is not yet supported. magic=%x\n", magic);
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
	snprintf(header[hc], 80, "%-5s%-3d= %20lu", "NAXIS", i+1, 
		 (unsigned long)(naxis[i])); header[hc][30]=' '; hc++;
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
/*
  Read fits header.
*/
static int 
read_fits_header(header_t *header, file_t *fp){
    //, char **str, uint32_t *magic, uint64_t *nx, uint64_t *ny){
    char line[82];//extra space for \n \0
    int end=0;
    int page=0;
    int bitpix=0;
    int naxis=0;
    while(!end){
	int start=0;
	if(page==0){
	    if(zfread_try(line, 1, 80, fp)) return -1; line[80]='\0';
	    if(strncmp(line, "SIMPLE", 6) && strncmp(line, "XTENSION= 'IMAGE", 16)){
		warning("Garbage in fits file %s\n", fp->fn);
		return -1;
	    }
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &bitpix)!=1) error("Unable to determine bitpix\n");
	    zfread(line, 1, 80, fp); line[80]='\0';
	    if(sscanf(line+10, "%20d", &naxis)!=1) error("Unable to determine naxis\n");
	    if(naxis>2) error("Data type not supported\n");
	    if(naxis>0){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (unsigned long *)&header->nx)!=1) error("Unable to determine nx\n");
	    }else{
		header->nx=0;
	    }
	    if(naxis>1){
		zfread(line, 1, 80, fp); line[80]='\0';
		if(sscanf(line+10, "%20lu", (unsigned long *)&header->ny)!=1) error("Unable to determine ny\n");
	    }else{
		header->ny=0;
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
		    if(header->str){
			header->str=realloc(header->str, strlen(header->str)+length+1+newline);
		    }else{
			header->str=malloc(length+1+newline); (header->str)[0]='\0';
		    }
		    strcat(header->str, hh);
		}
	    }
	}
	page++;
    }
    switch(bitpix){
    case -32:
	header->magic=M_FLT;
	break;
    case -64:
	header->magic=M_DBL;
	break;
    case 32:
	header->magic=M_INT32;
	break;
    case 16:
	header->magic=M_INT16;
	break;
    case 8:
	header->magic=M_INT8;
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
	if(header->magic!=MCC_ANY){
	    write_fits_header(fp, header->str, header->magic, 2, &header->nx, &header->ny);
	}
    }else{
	if(header->str){
	    write_bin_headerstr(header->str, fp);
	}
	//2015/01/27: We converted dummy header to a subtype.
	uint32_t magic2=M_SKIP;
	zfwrite(&magic2, sizeof(uint32_t), 1, fp);
	zfwrite(&header->magic,  sizeof(uint32_t), 1, fp);
	zfwrite(&header->nx, sizeof(uint64_t), 1, fp);
	zfwrite(&header->ny, sizeof(uint64_t), 1, fp);
    }
}
/**
   A unified header reading routine for .bin and .fits files. It read the array
   information and string header if any.  Return non zero value if reading failed*/
int read_header2(header_t *header, file_t *fp){
    int ans;
    header->str=NULL;
    if(fp->isfits){
	ans=read_fits_header(header, fp);
    }else{
	ans=read_bin_header(header, fp);
    }
    return ans;
}
/**
   calls read_header2 and abort if error happens.*/
void read_header(header_t *header, file_t *fp){
    if(read_header2(header, fp)){
	error("read_header failed for %s. Empty file?\n", fp->fn);
    }
}
/**
 * Get the length of mem to storge the header and its dimension, rounded to multiple of 8.
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
/*void write_timestamp(file_t *fp){
    char header[128];
    snprintf(header,128, "Created by MAOS Version %s on %s in %s\n",
	     PACKAGE_VERSION, myasctime(), HOST);
    write_bin_headerstr(header, fp);
    }*/


/**
   Write an 1-d or 2-d array into the file. First write a magic number that
   represents the data type. Then write two numbers representing the
   dimension. Finally write the data itself.
*/
void writearr(const void *fpn,     /**<[in] The file pointer*/
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
    if(!fp){
	warning("fp is empty\n");
	return;
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
    writearr(fn, 1, sizeof(double), M_DBL, NULL, p, nx, ny);
}

/**
   Unreference the mmaped memory. When the reference drops to zero, unmap it.
*/
void mmap_unref(struct mmap_t *in){
    if(in->nref>1){
	in->nref--;
    }else{
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

/**
   Open a file for write with mmmap. We don't provide a access control here for
   generic usage of the function. Lock on a special dummy file for access
   control.
*/
int mmap_open(char *fn, int rw){
    if(rw && disable_save){
	warning("Saving is disabled for %s\n", fn);
    }
    char *fn2=procfn(fn,rw?"w":"r",0);
    if(!fn2) return -1;
    if(fn2 && strlen(fn2)>=7&&!strncmp(fn2+strlen(fn2)-7,".bin.gz",7)){
	error("new_mmap does not support gzip\n");
    }
    int fd;
    if(rw){
	fd=open(fn2, O_RDWR|O_CREAT, 0600);
	if(fd!=-1 && ftruncate(fd, 0)){/*truncate the file. */
	    error("Unable to ftruncate file to 0 size\n");
	}
    }else{
	fd=open(fn2, O_RDONLY);
    }
    /*in read only mode, allow -1 to indicate failed. In write mode, fail.*/
    if(fd==-1 && rw){
	perror("open");
	error("Unable to create file %s\n", fn2);
    }
 
    free(fn2);
    return fd;
}

/**
   Initialize the header in the mmaped file. If header is not null, header0 points to its location in mmaped file. Upon exit, p0 points to the location of data p.
*/
void mmap_header_rw(char **p0, char **header0, uint32_t magic, long nx, long ny, const char *header){
    char *p=*p0;
    /*Always have a header to align the data. */
    if(header){
	uint64_t nlen=bytes_header(header)-24;
	((uint32_t*)p)[0]=(uint32_t)M_HEADER; p+=4;
	((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
	*header0=p;
	memcpy(p, header, strlen(header)+1); p+=nlen;
	((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
	((uint32_t*)p)[0]=(uint32_t)M_HEADER;p+=4;
    }else{
	*header0=NULL;
    }
    ((uint32_t*)p)[0]=(uint32_t)M_SKIP; p+=4;
    ((uint32_t*)p)[0]=(uint32_t)magic; p+=4;
    ((uint64_t*)p)[0]=(uint64_t)nx; p+=8;
    ((uint64_t*)p)[0]=(uint64_t)ny; p+=8;
    *p0=p;
}

/**
   Initialize the dimension from the header in the mmaped file.
*/
void mmap_header_ro(char **p0, uint32_t *magic, long *nx, long *ny, char **header0){
    char *p=*p0;
    char *header=NULL;
    while(((uint32_t*)p)[0]==M_HEADER){
	p+=4;
	long nlen=((uint64_t*)p)[0];p+=8;
	header=p;
	p+=nlen;
	if(nlen == ((uint64_t*)p)[0]){
	    p+=8;
	    if(((uint32_t*)p)[0]==M_HEADER){
		p+=4;
	    }else{
		header=NULL;
	    }
	}else{
	    header=NULL;
	}
	if(!header){
	    error("Parse head failed\n");
	}
    }
    if(!header){
	p=*p0;
    }
    if(((uint32_t*)p)[0]==M_SKIP) p+=4;
    *magic=((uint32_t*)p)[0]; p+=4;
    *nx=((uint64_t*)p)[0]; p+=8;
    *ny=((uint64_t*)p)[0]; p+=8;
    *p0=p;
    if(header0) *header0=header;
}

