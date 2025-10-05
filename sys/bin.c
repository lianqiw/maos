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
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdint.h>
#include <limits.h>
#include <ctype.h> /*isspace */
#include <errno.h>
#include "process.h"
#include "misc.h"
#include "path.h"
#include "thread.h"
#include "bin.h"
#include "scheduler_client.h"
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

#define USE_ZLIB_H 0
#if USE_ZLIB_H
#include <zlib.h> /*zlib.h in ubuntu sucks */
#else
#ifdef __cplusplus
extern "C"{
#endif
	typedef void* voidp;
	voidp gzopen(const char* path, const char* mod);
	voidp gzdopen(int fd, const char* mod);
	long  gztell(voidp gzfile);
	int   gzclose(voidp gzfile);
	int   gzwrite(voidp gzfile, const void* buf, unsigned len);
	int   gzread(voidp gzfile, voidp buf, unsigned len);
	int   gzseek(voidp file, long offset, int whence);
	int   gzrewind(voidp file);
	int   gzflush(voidp gzfile, int flush);
	const char* gzerror(voidp gzfile, int* error);
#ifdef __cplusplus
}
#endif
#endif
/*
  contains information about opened files.
*/
struct file_t{
	char* fn;  /**<The disk file name*/
	const char *msg; /**<Error message if any*/
	char *gstr;/**<Top level keyword string for fits file. */
	voidp gp;   /**<only used when gzipped*/
	int isgzip;/**<Is the file zipped.*/
	int isfits;/**<Is the file fits.*/
	int fd;    /**<The underlying file descriptor, used for locking. */
	int err;   /**<error happened. 1: eof, 2: other error*/
	int mode;  /**<0: read, 1: write*/
};
/*
  Describes the information about mmaped data. Don't unmap a segment of mmaped
  memory, which causes the whole page to be unmapped. Instead, reference count
  the mmaped file and unmap the segment when the nref dropes to 1.
*/
struct mem_t{
	void* mem;  /**<points to the beginning of memory.*/
	char* shm;/**<unlink shared memory when delete if set*/
	long n;   /**<length of mmaped memory.*/
	int kind; /**<Kind of memory. 0: heap, 1: mmap*/
	unsigned int nref; /**<Number of reference.*/
};

void swap_array(void* out, void* in, long bytes, int size){
	switch(size){
	case 1:
		if(out!=in){
			memcpy(out, in, bytes);
		}
		break;
	case 2:
	{
		uint16_t* pout=out;
		uint16_t* pin=in;
		long len=bytes>>1;
		for(long i=0; i<len; i++){
			pout[i]=swap2bytes(pin[i]);
		}
	}
	break;
	case 4:
	{
		uint32_t* pout=out;
		uint32_t* pin=in;
		long len=bytes>>2;
		for(long i=0; i<len; i++){
			pout[i]=swap4bytes(pin[i]);
		}
	}
	break;
	default://this is the longest data unit we need to handle.
	{
		uint64_t* pout=out;
		uint64_t* pin=in;
		long len=bytes>>3;
		for(long i=0; i<len; i++){
			pout[i]=swap8bytes(pin[i]);
		}
	}
	break;
	}
}
/**
   Disables file saving when set to 1.
*/
int disable_save=0;

/*
  Process the input file name and return file names that can be open to
  read/write. If the file name does not end with .bin or .bin.gz it will add to
  the end .bin. For read only
  access, it will also look into the path for files.
*/
static char* procfn(const char* fn, const char* mod){
	char* fn2;
	if(!fn) return NULL;
	if(fn[0]=='~'){
		fn2=(char*)malloc(strlen(HOME)+strlen(fn)+16);
		strcpy(fn2, HOME);
		strcat(fn2, fn+1);
	} else{
		fn2=(char*)malloc(strlen(fn)+16);
		strcpy(fn2, fn);
	}
	if(fn2[strlen(fn2)-1]=='\n'){//avoid accidental newline in filename
		fn2[strlen(fn2)-1]='\0';
	}
	/*If there is no recognized suffix, add .bin in the end. */
	if(mod[0]=='r'){
		char* fnr=NULL;
		if(!(fnr=search_file(fn2, 1))){/*If does not exist. */
			if(!check_suffix(fn2, ".bin")&&!check_suffix(fn2, ".bin.gz")
				&&!check_suffix(fn2, ".fits")&&!check_suffix(fn2, ".fits.gz")){
			 //Try adding suffix.
				const char* suf[]={".bin",".fits",".bin.gz",".fits.gz"};
				for(unsigned int is=0; is<4; is++){
					strcat(fn2, suf[is]);
					if(!(fnr=search_file(fn2, 1))){
						fn2[strlen(fn2)-strlen(suf[is])]='\0';
					} else{
						break;
					}
				}
			}
		}
		//It is ok to have NUL fnr.
		free(fn2);
		fn2=fnr;
	} else if(mod[0]=='w'||mod[0]=='a'){
		if(!check_suffix(fn2, ".bin")&&!check_suffix(fn2, ".bin.gz")
			&&!check_suffix(fn2, ".fits")&&!check_suffix(fn2, ".fits.gz")){
			strcat(fn2, ".bin");
		}
		/*if (check_suffix(fn2, ".gz")){
			fn2[strlen(fn2)-3]='\0';
			}*/
		if(disable_save&&mystrcmp(fn2, DIRCACHE)){
			//When saving is disabled, allow writing to cache folder.
			{
				static int printed=0;
				if(!printed){
					warning("Saving is disabled for %s.\n", fn2);
					print_backtrace();
					printed=1;
				}
			}
			free(fn2);
			fn2=NULL;
		}
		if(fn2&&islink(fn2)){/*remove old file to avoid write over a symbolic link. */
			if(remove(fn2)){
				warning("Failed to remove %s\n", fn2);
			}
		}
	} else{
		warning("Invalid mode\n");
	}
	return fn2;
}
/**
   Test whether a bin file exist.
*/
int zfexist(const char* format, ...){
	format2fn;
	char* fn2=procfn(fn, "rb");
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
void zftouch(const char* format, ...){
	format2fn;
	char* fn2=procfn(fn, "rb");
	if(fn2&&utimes(fn2, NULL)){
		perror("zftouch failed");
	}
	free(fn2);
}
//PNEW(lock);
/*
  convert a socket or fd into file_t.
*/

file_t* zfdopen(int fd){
	file_t* fp=mycalloc(1, file_t);
	if(fp){
		fp->isgzip=0;
		fp->fd=fd;
	}
	return fp;
}
/**
   Adaptively call open or shm_open depending on the name.
*/
static inline int myopen(const char* name, int oflag, mode_t mode){
	if(IS_SHM(name)){
	/*if(!mystrcmp(name, "/shm")){
		name+=4;
	}
	if(name[0]=='/'){
		name++;
	}
	char name2[NAME_MAX];
	snprintf(name2, NAME_MAX, "/maos_%d_%s", PID, name);
	name=name2;*/
		return shm_open(name, oflag, mode);
	} else{
		return open(name, oflag, mode);
	}
}
//#define OPEN(name, flags...) (IS_SHM(name)?(info("shm_open: %s\n", name),shm_open(name, flags)):open(name, flags))

/**
   Open a bin file for read/write access.

   Whether the file is gzipped or not is automatically determined when open for
   reading. When open for writing, if the file name ends with .bin or .fits,
   will not be gzipped. If the file has no suffix, the file will be gzipped and
   .bin is appened to file name.
*/
file_t* zfopen(const char* fni, const char* mod){
	//LOCK(lock);//nothing to protect
	file_t* fp=mycalloc(1, file_t);
	if(!fp) return NULL;
	fp->fd=-1;
	const char* fn2=fp->fn=procfn(fni, mod);
	if(!fn2){
		if(mod[0]=='r'){
			warning("%s does not exist for read\n", fni);
			printpath();
		}
		goto fail;
	}
	if(mod[0]=='w'||mod[0]=='a'){
		if(check_suffix(fn2, ".gz")){
			fp->isgzip=1;
			if(mod[0]=='w'){
				fp->fn[strlen(fn2)-3]='\0';
			}
		} else{
			fp->isgzip=0;
		}
	}
	/*Now open the file to get a fd number that we can use to lock on the file.*/
	switch(mod[0]){
	case 'r':/*read only */
		if((fp->fd=myopen(fn2, O_RDONLY, 0666))==-1){
			perror("open for read");
		} else{//file exist
			fp->mode=0;
			if(!mystrcmp(fn2, DIRCACHE)){
				utimes(fn2, NULL);
			}
		}
		break;
	case 'w':/*write */
	case 'a':
		if((fp->fd=myopen(fn2, O_RDWR|O_CREAT, 0666))==-1){
			perror("open for write");
		} else{
			if(lockf(fp->fd, F_TLOCK, 0)){//2018-03-07: Removed LOCK_NB. 2021-01-11 replaced with lockf
				perror("lockf");
				close(fp->fd);
				fp->fd=-1;
			} else if(mod[0]=='w'&&ftruncate(fp->fd, 0)){/*Need to manually truncate the file. */
				perror("ftruncate");
			}
			fp->mode=1;
		}
		break;
	default:
		warning("Unknown mod=%s. Supported are r, w, a\n", mod);
		goto fail;
	}
	if(fp->fd==-1){
		warning("Unable to open file %s for %s (%s)\n", fn2, mod[0]=='r'?"reading":"writing", strerror(errno));
		goto fail;
	}
	fcntl(fp->fd, F_SETFD, FD_CLOEXEC);
	/*check fn instead of fn2. if end of .bin or .fits, disable compressing.*/
	if(mod[0]=='r'){
		uint16_t magic;
		if(read(fp->fd, &magic, sizeof(uint16_t))!=sizeof(uint16_t)){
			warning("Unable to read magic from %s.\n", fn2);
			goto fail;
		} else{
			if(magic==0x8b1f){
				fp->isgzip=1;
			} else{
				fp->isgzip=0;
			}
		}
		lseek(fp->fd, 0, SEEK_SET);
	}
	if(fp->isgzip){
		if(!(fp->gp=gzdopen(fp->fd, mod))){
			warning("Error gzdopen for %s\n", fn2);
			goto fail;
		}
	}
	if(check_suffix(fn2, ".fits")||check_suffix(fn2, ".fits.gz")){
		fp->isfits=1;
	}
	/*if(mod[0]=='w' && !fp->isfits){
	write_timestamp(fp);
	}*/
	//UNLOCK(lock);
	return fp;
fail:
  //UNLOCK(lock);
	if(fp->fd!=-1) close(fp->fd);
	free(fp->fn); free(fp); fp=NULL;
	return fp;
}
/**
 * print error message
 * */
static void zferrprint(file_t *fp){
	switch(fp->err){
	case 0:
		break;//dbg("%s: Clear error\n", fp->fn);
	case 1:
		dbg3("%s: EOF encountered\n", fp->fn);
		break;
	case 2:
		warning("%s: Unknown error\n", fp->fn);
		break;
	default:
		warning("%s: Invalid entry: %d\n", fp->fn, fp->err);
	}
}
/**
 * set error field
 * */
static void zferr(file_t *fp, int err){
	if(!fp) return;
	fp->err=err;
	zferrprint(fp);
}
/**
 * Return the error number
 * */
int zferrno(file_t *fp){
	return fp?fp->err:-1;
}
/**
   Open a bin file or socket for read/write access.

   Whether the file is gzipped or not is automatically determined when open for
   reading. When open for writing, if the file name ends with .bin or .fits,
   will not be gzipped. If the file has no suffix, the file will be gzipped and
   .bin is appened to file name.
*/
/*file_t* zfopen(const char* fn, const char* mode){
	file_t* fp;
	fp=zfopen_try(fn, mode);
	if(!fp&&!disable_save){
		error("Open file %s for %s failed\n", fn, mode[0]=='r'?"read":"write");
	}
	return fp;
}*/
/**
   Return the underlining filename
*/
const char* zfname(file_t* fp){
	return fp?fp->fn:NULL;
}

/**
   Return 1 when it is a fits file.
*/
int zfisfits(file_t* fp){
	return (fp&&fp->isfits)?1:0;
}
/**
   Close the file.
*/
void zfclose(file_t* fp){
	if(!fp) return;
	if(fp->fn && fp->mode==0){
		long pos=zfpos(fp);
		long len=zflen(fp);
		if(pos<len){
			warning("zfclose %s: there is extra data not read. pos=%ld, len=%ld\n", fp->fn, pos, len);
		}
	}
	//LOCK(lock);//causes race condition with zfopen_try
	if(fp->isgzip){
		gzclose(fp->gp);
	} else{
		close(fp->fd);
	}
	//close(fp->fd);//keep fd open
	free(fp->gstr);
	free(fp->fn);
	free(fp);
	//UNLOCK(lock);
}
/**
 * Wrap normal and gzip file write operation.
 * */
int zfwrite_wrap(const void* ptr, const size_t tot, file_t* fp){
	if(!fp||fp->fd<0||fp->err||!ptr){
		warning("zfwrite_wrap: error encountered, cancelled.\n");
		return (fp&&fp->err)?fp->err:-1;
	} else if(fp->isgzip){
		return gzwrite(fp->gp, ptr, tot);
	} else{
		return write(fp->fd, ptr, tot);
	}
}
/*
  Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
  fwrite. Follows the interface of fwrite.
*/
int zfwrite_do(const void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	size_t tot=size*nmemb;
	while(!fp->err){
		int count=zfwrite_wrap(ptr, tot, fp);
		if(count>=0){//=0 is not an error
			if((size_t)count<tot){
				ptr+=count;
				tot-=count;
			} else{
				break;
			}
		} else{
			zferr(fp, 2);
		}
	}
	if(fp->err){
		if(errno) perror("zfwrite_do");
		if(errno==ENOSPC){
			error("not space left in device\n");
		}else{
			warning("write %lu bytes to %s failed\n", tot, fp->fn);
		}
	}
	return fp->err;
}

/**
   Handles byteswapping in fits file format then call zfwrite_do to do the actual writing.
*/
int zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	/*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
	if(!ptr || !size || !nmemb || !fp || fp->err){
		warning("zfwrite: invalid input or error in file.\n");
		return -1;
	}
	int ans=0;
	if(fp->isfits&&BIGENDIAN==0){
		int length=size*nmemb;
		/* write a block of 2880 bytes each time, with big-endianness.*/
		const int bs=2880;
		char junk[bs];
		int nb=(length+bs-1)/bs;
		char* in=(char*)ptr;
		for(int ib=0; ib<nb&&!ans; ib++){
			int nd=length<bs?length:bs;
			swap_array(junk, in, nd, size);
			/* use bs instead of nd to test tailing blanks*/
			in+=bs; length-=bs;
			if(length<0){
				memset(junk+nd, 0, (bs-nd)*sizeof(char));
			}
			ans=zfwrite_do(junk, sizeof(char), bs, fp);
		}
	} else{
		ans=zfwrite_do(ptr, size, nmemb, fp);
	}
	return ans;
}
/**
 * Wrap normal and gzip file operation.
 * */
int zfread_wrap(void* ptr, const size_t tot, file_t* fp){
	if(!fp||fp->fd<0||fp->err||!ptr){
		warning("zfread_wrap: error encountered\n");
		return -1;
	} else if(fp->isgzip){
		return gzread(fp->gp, ptr, tot);
	} else{
		return read(fp->fd, ptr, tot);
	}
}
/**
   Read from the file. Handle partial read.
*/
int zfread_do(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	ssize_t tot=size*nmemb;
	while(!fp->err){
		long count=zfread_wrap(ptr, tot, fp);
		if(count>0){
			if(count<tot){
				ptr+=count;
				tot-=count;
			} else{
				break;
			}
		} else if(count==0){//eof
			zferr(fp, 1);
		} else{//unknown error
			warning("%s: unknown error happend during reading.\n", fp->fn);
			zferr(fp, 2);
		}
	}
	return fp->err;
}
/**
   Handles byteswapping in fits file format then call zfread_do to do the actual writing.
   return 0 when succeed and -1 when error.
*/
int zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	int ans=0;
	if(fp->err){
		ans=-1;
	} else if(fp->isfits&&size>1){/*need to do byte swapping.*/
		const long bs=2880;
		long length=size*nmemb;
		long nb=(length+bs-1)/bs;
		char* out=(char*)ptr;
		for(int ib=0; ib<nb; ib++){
			int nd=length<bs?length:bs;
			ans=zfread_do(out, sizeof(char), nd, fp);
			if(nd<bs){//read the padding
				zfseek(fp, bs-nd, SEEK_CUR);
			}
			if(ans) break;
			swap_array(out, out, nd, size);
			out+=bs;
			length-=bs;
		}
	} else{
		ans=zfread_do(ptr, size, nmemb, fp);
	}
	return ans;
}
/**
   Move the current position pointer, like fseek
*/
long zfseek(file_t* fp, long offset, int whence){
	if(fp->isgzip){
		return gzseek(fp->gp, offset, whence);
	} else{
		return lseek(fp->fd, offset, whence);
	}
}
/**
   Move the file position pointer to the beginning
*/
void zfrewind(file_t* fp){
	if(fp->isgzip){
		if(!gzrewind(fp->gp)){
			zferr(fp, 0);
		}
	} else{
		if(lseek(fp->fd, 0, SEEK_SET)!=-1){
			zferr(fp, 0);
		}
	}
}
/**
   Tell position pointer in file.
*/
long zfpos(file_t* fp){
	/*if(fp->isgzip){
		return gztell(fp->gp);
	} else{
		return lseek(fp->fd, 0, SEEK_CUR);
	}*/
	return zfseek(fp, 0, SEEK_CUR);
}
/**
   Return 1 if end of file is reached.
*/
int zfeof(file_t* fp){
	if(fp->err) return fp->err==1;//already EOF error
	long pos=zfpos(fp);
	long end=zfseek(fp, 0, SEEK_END);
	zfseek(fp, pos, SEEK_SET);
	return pos!=end;
}
/**
 * Return file length
 * */
long zflen(file_t *fp){
	long pos=zfpos(fp);
	long end=zfseek(fp, 0, SEEK_END);
	zfseek(fp, pos, SEEK_SET);
	return end;
}

/**
   Flush the buffer.
*/
/*
void zflush(file_t* fp){
	if(fp->isgzip){
		gzflush(fp->gp, 4);
	}
}*/
#define CHECK_EOF(A, errmsg) if(A){if(fp->err==1) goto read_error_eof; else  {fp->msg=errmsg;goto read_error;}}//eof is not error
#define CHECK_ERR(A, errmsg) if(A){if(!fp->err) zferr(fp, 2); fp->msg=errmsg; goto read_error;} //eof is error
/*
  Obtain the current magic number, string header, and array size from a bin
  format file. In file, the string header is located before the data magic.
*/
static int
read_bin_header(header_t* header, file_t* fp){
	uint32_t magic=0, magic2=0;
	uint64_t nlen=0, nlen2=0;
	memset(header, 0, sizeof(header_t));

	if(fp->isfits){
		warning("Please use read_fits_header for fits file.\n");
		return -1;
	}
	while(!fp->err){
		/*read the magic number.*/
		CHECK_EOF(zfread(&magic, sizeof(uint32_t), 1, fp),"Read frame start failed.");
		/*If it is hstr, read or skip it.*/
		if(magic==M_EOD){//end of data. useful in socket i/o to indicate end of current dataset
			return 1;
		}else if((magic&M_SKIP)==M_SKIP){//dummy
			continue;
		} else if(magic==M_COMMENT){//comment
			CHECK_ERR(zfread(&nlen, sizeof(uint64_t), 1, fp),"Read length of comment failed");
			if(nlen>0){
				/*zfseek failed in cygwin (gzseek). so we always read in the hstr instead.*/
				char hstr2[nlen];
				CHECK_ERR(zfread(hstr2, 1, nlen, fp), "Read comment failed");
				hstr2[nlen-1]='\0'; /*make sure it is NULL terminated. */
				if(header->str){
					header->str=(char*)realloc(header->str, ((header->str)?strlen(header->str):0)+strlen(hstr2)+1);
					strncat(header->str, hstr2, nlen);
				} else{
					header->str=strdup(hstr2);
				}
			}
			CHECK_ERR(zfread(&nlen2, sizeof(uint64_t), 1, fp)||nlen2!=nlen, "Read length of comment 2 failed");
			CHECK_ERR(zfread(&magic2, sizeof(uint32_t), 1, fp)||magic2!=magic, "Read data type 2 failed");
		} else{ //Finish
			header->magic=magic;
			CHECK_ERR(zfread(&header->nx, sizeof(uint64_t), 1, fp), "Read nx failed");
			CHECK_ERR(zfread(&header->ny, sizeof(uint64_t), 1, fp), "Real ny failed");
			if((magic&0x6400)!=0x6400){
				warning("wrong magic number=%x at %ld. cancelled.\n", magic, zfpos(fp));
				zferr(fp, 2);
				goto read_error;
			} else{
				return 0;
			}
		}
	}/*while*/
read_error:
	if(fp->err){
		warning("read_bin_header: error=%d.\n", fp->err);
	}
read_error_eof:
	header->magic=0;
	header->nx=0;
	header->ny=0;
	free(header->str); header->str=NULL;
	return fp->err?fp->err:-1;
}

//emulate mempcpy in linux
#define MEMPCPY(dest, src, type)\
    dest=(char *)memcpy(dest, &src, sizeof(type)) + sizeof(type);

/**
	Write bin file header for keyword and data dimension.
	Updated from previous version to reduce syscalls.

	First write magic number, then the length of the header, then the header,
	then the length of the header again, then the magic number again. The two set
	of strlen and header are used to identify the header from the end of the file
	and also to verify that what we are reading are indeed header. The header may
	be written multiple times. They will be concatenated when read. The header
	should contain key=value entries just like the configuration files. The
	entries should be separated by new line charactor.
*/
static void write_bin_header(const header_t* header, file_t* fp){
	const uint64_t nlen8=header->str?(((strlen(header->str)+8)>>3)<<3):0;
	uint64_t nbyte=(nlen8?(nlen8+4*2+8*2):0)+3*8;
	char* buf=mycalloc(nbyte, char);
	char*tmp=buf;
	if(nlen8){
		uint32_t mc=M_COMMENT;
		//Don't need to write M_SKIP here because we have a pair of magic.
		//The length of header->str is rounded to multiple of 8 bytes.
		MEMPCPY(tmp, mc, uint32_t);
		MEMPCPY(tmp, nlen8, uint64_t);
		strncpy(tmp, header->str, nlen8); tmp[nlen8-1]=0;tmp+=nlen8;
		MEMPCPY(tmp, nlen8, uint64_t);
		MEMPCPY(tmp, mc, uint32_t);
	}
	uint32_t mp=M_SKIP;
	MEMPCPY(tmp, mp, uint32_t);
	MEMPCPY(tmp, header->magic, uint32_t);
	MEMPCPY(tmp, header->nx, uint64_t);
	MEMPCPY(tmp, header->ny, uint64_t);
	zfwrite(buf, nbyte, 1, fp);
	free(buf);
}
/*
  Write fits header. str is extra header that will be put in fits comment
*/
static void
write_fits_header(file_t* fp, const char* str, uint32_t magic, int count, ...){
	uint64_t naxis[count+1];
	va_list ap;
	va_start(ap, count);              /*Initialize the argument list. */
	int empty=0;
	for(int i=0; i<count; i++){
		uint64_t* addr=va_arg(ap, uint64_t*);  /*Get the next argument value.   */
		if((*addr)==0) empty=1;
		naxis[i]=*addr;
	}
	va_end(ap);

	if(empty) count=0;
	int iscomplex=0;
	int bitpix=0;
	switch(magic&0xFFFF){
	case M_FLT:
		bitpix=-32;
		break;		
	case M_ZMP:
		bitpix=-32;
		iscomplex=1;
		break;
	case M_DBL:
		bitpix=-64;
		break;
	case M_CMP:
		bitpix=-64;
		iscomplex=1;
		break;
	case M_INT64:
		bitpix=64;
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
		warning("Data type is not yet supported. magic=%x\n", magic);
		return;
	}
	if(iscomplex){//insert 2 to beginning of naxis
		for(int i=count; i>0; i--){
			naxis[i]=naxis[i-1];
		}
		naxis[0]=2;
		count++;
	}
	const int nh=36;//each fits page can only contain 36 headers.
	char header[nh][80];
	memset(header, ' ', sizeof(char)*36*80);
	int hc=0;
	if(fp->isfits==1){
		snprintf(header[hc], 80, "%-8s= %20s", "SIMPLE", "T");    header[hc][30]=' '; hc++;
		fp->isfits++;
	} else{
		snprintf(header[hc], 80, "%-8s= %s", "XTENSION", "'IMAGE   '");    header[hc][20]=' '; hc++;
	}
	snprintf(header[hc], 80, "%-8s= %20d", "BITPIX", bitpix); header[hc][30]=' '; hc++;
	snprintf(header[hc], 80, "%-8s= %20d", "NAXIS", count);   header[hc][30]=' '; hc++;

#define FLUSH_OUT /*write the page if ready */		\
	if(hc==nh){					\
	    zfwrite(header, sizeof(char), 36*80, fp);	\
	    memset(header, ' ', sizeof(char)*36*80);	\
	    hc=0;					\
	}
	for(int i=0; i<count; i++){
		FLUSH_OUT;
		snprintf(header[hc], 80, "%-5s%-3d= %20lu", "NAXIS", i+1,
			(unsigned long)(naxis[i])); header[hc][30]=' '; hc++;
	}
	if(fp->isfits==1){//Write the extend keyword which does not mendate extension to be present.
		snprintf(header[hc], 80, "%-8s= %20s", "EXTEND", "T");    header[hc][30]=' '; hc++;
	} else{
		snprintf(header[hc], 80, "%-8s= %20s", "PCOUNT", "0");    header[hc][30]=' '; hc++;
		snprintf(header[hc], 80, "%-8s= %20s", "GCOUNT", "1");    header[hc][30]=' '; hc++;
	}
	if(str){
		const char* str_end=str+strlen(str);
		while(isspace(str[0])&&str<str_end) str++;
		while(str&&str<str_end){
			FLUSH_OUT;
			const char* nl=strchr(str, '\n');//separation of keys
			const char* nc=strchr(str, ';');//separation of keys
			const char* eq=strchr(str, '=');
			if(!nl) nl=str_end;
			if(nc&&nc<nl) nl=nc;
			if(eq>nl) eq=0;
			int length;
			if(eq){
				length=nl-eq;
			} else{
				length=nl-str+1;
			}
			if(length>70) length=70;//each line can contain maximum 70 values

			if(eq){//there is an equal sign.
				//strncpy(header[hc], str, MIN(8, (eq-str)));//replaced by below
				for(int ic=0; ic<MIN(8, (eq-str)); ic++){//fits standard requires upper case keywords
					if(!str[ic]) break;
					header[hc][ic]=toupper((unsigned char)str[ic]);
				}
				header[hc][8]='=';
				header[hc][9]=' ';
				//strncpy(header[hc]+10, eq+1, length);//replaced by below
				int slash_found=0;//convert char before / to upper case for the FITS standard
				for(int ic=0; ic<length; ic++){
					if(eq[1+ic]=='/'){
						slash_found=1;
					}
					header[hc][10+ic]=slash_found?eq[1+ic]:toupper((unsigned char)eq[1+ic]);
					if(!eq[1+ic]) break;
				}
			} else{
				strncpy(header[hc], "COMMENT   ", 11);
				strncpy(header[hc]+10, str, length);
			}
			if(nl){//Replace \n by space
				header[hc][10+length-1]=' ';
			}
			hc++;
			//update str to after \n
			if(eq){
				str=eq+length+1;
			} else{
				str+=length;
			}
			while(str<str_end && isspace(str[0])) str++;
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
read_fits_header(header_t* header, file_t* fp){
	char line[82];//extra space for \n \0
	int end=0;
	int page=0;
	int bitpix=0;
	int naxis=0;
	int nc=1;//1: default. 2: is complex number
	while(!end){
		int start=0;
		if(page==0){
			CHECK_EOF(zfread(line, 1, 80, fp), "Read line 1 failed");
			line[80]='\0';
			CHECK_ERR(strncmp(line, "SIMPLE", 6)&&strncmp(line, "XTENSION= 'IMAGE", 16), "Invalid frame data type");
			CHECK_ERR(zfread(line, 1, 80, fp), "Read line 2 failed"); line[80]='\0';
			CHECK_ERR(sscanf(line+10, "%20d", &bitpix)!=1, "Read bitpix failed");
			CHECK_ERR(zfread(line, 1, 80, fp), "Read line 3 failed"); line[80]='\0';
			CHECK_ERR(sscanf(line+10, "%20d", &naxis)!=1, "Read naxis failed");
			CHECK_ERR(naxis>3, "Unable to handle naxis>3");
			if(naxis>2){
				CHECK_ERR(zfread(line, 1, 80, fp), "Read line 4 failed"); line[80]='\0';
				CHECK_ERR(sscanf(line+10, "%20d", &nc)!=1, "Read nc failed");
			}
			if(naxis>0){
				CHECK_ERR(zfread(line, 1, 80, fp), "Read line 4 failed"); line[80]='\0';
				CHECK_ERR(sscanf(line+10, "%20lu", (unsigned long*)&header->nx)!=1, "Read nx failed");
			} else{
				header->nx=0;
			}
			if(naxis>1){
				CHECK_ERR(zfread(line, 1, 80, fp), "Read line 5 failed"); line[80]='\0';
				CHECK_ERR(sscanf(line+10, "%20lu", (unsigned long*)&header->ny)!=1, "Read ny failed");
			} else{
				header->ny=0;
			}
			start=3+naxis;
		}
		int was_comment=0;
		for(int i=start; i<36; i++){
			CHECK_ERR(zfread(line, 1, 80, fp), "Read data line failed"); 
			line[80]='\0';
			if(!strncmp(line, "END", 3)){
				end=1;
			} else{
				char* hh=line;
				int length=80;
				int is_comment=0;
				//Ignore a few keys
				if(!strncmp(line, "EXTEND", 6)){
					continue;
				} else if(!strncmp(line, "PCOUNT", 6)){
					continue;
				} else if(!strncmp(line, "GCOUNT", 6)){
					continue;
				} else if(!strncmp(line, "COMMENT", 7)){
					hh=line+10;
					length-=10;
					is_comment=1;
				}
				//Remove trailing space.
				while(length>0&&isspace((int)hh[length-1])){
					hh[length-1]='\0';
					length--;
				}

				if(length>0){
					if(header->str){
						header->str=myrealloc(header->str, strlen(header->str)+length+3, char);
						if(!(is_comment&&was_comment)){
							strcat(header->str, "\n\0");
						}
					} else{
						header->str=(char*)malloc(length+3); (header->str)[0]='\0';
					}
					strcat(header->str, hh);
				}
				was_comment=is_comment;
			}
		}
		page++;
	}
	switch(bitpix){
	case -32:
		header->magic=nc==2?M_ZMP:M_FLT;
		break;
	case -64:
		header->magic=nc==2?M_CMP:M_DBL;
		break;
	case 64:
		header->magic=M_INT64;
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
		zferr(fp, 2);
		warning("Unsupported bitpix=%d.\n", bitpix);
		goto read_error;
	}
	return 0;
read_error:
	if(fp->err){
		warning("read_fits_header: error=%d\n", fp->err);
	}
read_error_eof:
	header->magic=0;
	header->nx=0;
	header->ny=0;
	return fp->err?fp->err:-1;
}
/**
   A unified header writing routine for .bin and .fits files. It write the array
   information and string header if any.  */
void write_header(const header_t* header, file_t* fp){
	if(fp->isfits){
		if(header->magic!=MCC_ANY){
			write_fits_header(fp, header->str?header->str:fp->gstr, header->magic, 2, &header->nx, &header->ny);
		}else if(header->str){
			fp->gstr=strdup(header->str);//save global keyword to be used for each page
		}
	} else{
		write_bin_header(header, fp);
	}
}
/**
   A unified header reading routine for .bin and .fits files. It read the array
   information and string header if any.  Return non zero value if reading failed*/
int read_header(header_t* header, file_t* fp){
	int ans;
	header->str=NULL;
	if(fp->err){
		dbg("fp->err=%d, read_header canceled.\n", fp->err);
		ans=-1;
	} else if(fp->isfits){
		ans=read_fits_header(header, fp);
	} else{
		ans=read_bin_header(header, fp);
	}
	return ans;
}
/**
 * Get the length of mem to storge the header and its dimension, rounded to multiple of 8.
 */
uint64_t bytes_header(const char* keywords){
	if(keywords){
		uint64_t len=strlen(keywords)+1;
		if(len%8!=0){
			len=(len/8+1)*8;
		}
		return len+24;
	} else{
		return 0;
	}
}
/**
   Write the time stamp as header into current location in the file.
*/
/*void write_timestamp(file_t *fp){
	char header[128];
	snprintf(header,128, "Created by MAOS Version %s on %s in %s\n",
		 PACKAGE_VERSION, myasctime(0), HOST);
	write_bin_headerstr(header, fp);
	}*/


/**
   Write an 1-d or 2-d array into the file. First write a magic number that
   represents the data type. Then write two numbers representing the
   dimension. Finally write the data itself.

   returns the file offset of the array start
*/
long writearr(const void* fpn,     /**<[in] The file pointer*/
	const int isfn,      /**<[in] Is this a filename or already opened file*/
	const size_t size,   /**<[in] Size of each element*/
	const uint32_t magic,/**<[in] The magic number. see bin.h*/
	const char* str,     /**<[in] The header as string*/
	const void* p,       /**<[in] The data of the array*/
	const uint64_t nx,   /**<[in] Number of rows. this index changes fastest*/
	const uint64_t ny    /**<[in] Number of columns. 1 for vector*/
){
	file_t* fp;
	if(isfn){
		fp=zfopen((char*)fpn, "wb");
	} else{
		fp=(file_t*)fpn;
	}
	if(!fp){
		if(!disable_save){
			warning_once("writearr failed: fp is empty.\n");
		}else{
			warning_once("writearr failed: empty fp.\n");
		}
		return -1;
	}
	header_t header={magic, nx, ny, (char *)str};
	write_header(&header, fp);
	//long pos=zfpos(fp);
	int ans=0;
	if(nx*ny>0){
		if((ans=zfwrite(p, size, nx*ny, fp))){
			error("writearr failed\n");
		}
	}
	//else dbg("writearr: size is zero: %s %ldx%ld\n", zfname(fp), nx, ny);
	if(isfn) zfclose(fp);
	return ans;
}


/**
   Unreference the mmaped memory. When the reference drops to zero free or unmap it.
*/
unsigned int mem_unref(mem_t** pin){
	if(!pin) return 0;
	mem_t* in=*pin;
	if(!in) return 0;
	*pin=0;//this is important to avoid dangling pointer.
	unsigned int nref=0;
	if(!(nref=atomic_sub_fetch(&(in->nref), 1))){//deallocate
		switch(in->kind){
		case 0:
			free(in->mem);
			break;
		case 1:
			msync(in->mem, in->n, MS_SYNC);
			munmap(in->mem, in->n);
			break;
		default:
			warning("invalid in->kind=%d\n", in->kind);
		}
		if(in->shm){
			shm_unlink(in->shm);
		}
		free(in);
	}
	return nref;
}
/**
   Create a mem_t object.

   Notice that nref is set to 0. Default is set to reference heap memory.
*/
mem_t* mem_new(void* p){
	mem_t* out=mycalloc(1, mem_t);
	if(out){
		out->mem=p;
	}
	return out;
}
/**
   Add a reference to a mem_t.
*/
mem_t* mem_ref(mem_t* in){
	if(in){
		atomic_add_fetch(&in->nref, 1);
	}
	return in;
}
/**
   Check whether it is a referenced data.
*/
int mem_nref(const mem_t* in){
	return in?(in->kind+in->nref):0;
}
/**
   Replace internal vector. If p is null, checking only. Use with care.
*/
void mem_replace(mem_t* in, void* p){
	if(!in||in->nref>1||in->kind){
		error("Replacing referenced or mmaped memory\n");
	} else{
		in->mem=p;
	}
}
/**
   Return internal pointer
*/
void* mem_p(const mem_t* in){
	return in?in->mem:NULL;
}
/**
   Open a file for write with mmmap. We don't provide an access control here for
   generic usage of the function. Lock on a special dummy file for access
   control.
*/
mem_t* mmap_open(const char* fn, size_t msize, int rw){
	int fd=-1;
	char* fn2=NULL;
	if(!fn){
		if(!rw){
			error("fn cannot be null for reading.\n");return NULL;
		}
		if(!msize){
			error("size cannot be zero when fn is NULL.\n"); return NULL;
		}
	} else{
		if(rw&&disable_save&&!IS_SHM(fn)&&mystrcmp(fn, DIRCACHE)){
			warning("Saving is disabled for %s\n", fn);
			print_backtrace();
			fn=NULL;
		}
	}
	if(!rw||fn){
		fn2=procfn(fn, rw?"w":"r");
		if(!fn2) return NULL;
		if(!check_suffix(fn2, ".bin")){
			error("mmap only support .bin file\n");
			return NULL;
		}

		if(rw){
			fd=myopen(fn2, O_RDWR|O_CREAT, 0600);
			/*First truncate the file to 0 to delete old data. */
			if(fd==-1||ftruncate(fd, msize)==-1){
				error("Unable to open or ftruncate file %s to %zu size\n", fn2, msize);
			}
		} else{
			fd=myopen(fn2, O_RDONLY, 0600);
		}
		/*in read only mode, allow -1 to indicate failed. In write mode, fail.*/
		if(fd==-1){
			perror("open");
			if(rw){
				error("Unable to create file %s\n", fn2);
			}
			return NULL;
		}

		if(!rw&&!msize){
			msize=flen(fn2);
		}

	}
	//Notice that changed made in MAP_PRIVATE mode are not saved to file. So cannot be used
	//mmap should not be used for writing as it causes a lot of disk activity during small updates.
	void* p=mmap(NULL, msize, (rw?PROT_WRITE:0)|PROT_READ, /*(fd==-1?MAP_ANONYMOUS:0)|*/MAP_SHARED, fd, 0);
	if(fd!=-1) close(fd);//it is ok to close fd after mmap.
	mem_t* mem=0;
	if(p==MAP_FAILED){
		perror("mmap");
		error("mmap failed\n");
	} else{
		mem=mem_new(p);
		mem->kind=1;
		mem->n=msize;
		if(IS_SHM(fn2)){
			mem->shm=fn2;
			fn2=NULL;
		}
	}
	free(fn2);
	return mem;
}

/**
   Initialize the header in the mmaped file. If header is not null, header0 points to its location in mmaped file. Upon exit, p0 points to the location of data p.
*/
void mmap_write_header(char** p0, uint32_t magic, long nx, long ny, const char* keywords){
	char* p=*p0;
	/*Always have a header to align the data. */
	if(keywords){
		uint64_t nlen=bytes_header(keywords)-24;
		((uint32_t*)p)[0]=(uint32_t)M_COMMENT; p+=4;
		((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
		memcpy(p, keywords, strlen(keywords)+1); p+=nlen;
		((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
		((uint32_t*)p)[0]=(uint32_t)M_COMMENT;p+=4;
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
void mmap_read_header(char** p0, uint32_t* magic, long* nx, long* ny, const char** pkeywords){
	char* p=*p0;
	char* keywords=NULL;
	while(((uint32_t*)p)[0]==M_COMMENT){
		p+=4;
		uint64_t nlen=((uint64_t*)p)[0];p+=8;
		keywords=p;
		p+=nlen;
		if(nlen==((uint64_t*)p)[0]){
			p+=8;
			if(((uint32_t*)p)[0]==M_COMMENT){
				p+=4;
			} else{
				keywords=NULL;
			}
		} else{
			keywords=NULL;
		}
		if(!keywords){
			error("Parse head failed\n");
		}
	}
	if(!keywords){
		p=*p0;
	}
	if(((uint32_t*)p)[0]==M_SKIP) p+=4;
	*magic=((uint32_t*)p)[0]; p+=4;
	*nx=((uint64_t*)p)[0]; p+=8;
	*ny=((uint64_t*)p)[0]; p+=8;
	*p0=p;
	if(pkeywords) *pkeywords=keywords;
}
#define USE_ASYNC 0
#if USE_ASYNC
#include <aio.h>
#endif
struct async_t{
	file_t* fp;
	const char* p;  //original memory
	long pos; //start position of array in file
	long prev;//previous saved index 
#if USE_ASYNC
	struct aiocb aio;
#endif	
};
/**
 * Allocates space in fp for our block of data by writing header first and skip to the end of data offset.
 * It has the same interface as writearr except only fp is allowed
 * Use aio_* api if USE_ASYNC is set, otherwise use pwrite()
 * aio_* api launches threads for async writing.
 * pwrite is much simpler for our usecase. There is no impact on simulation step timing.
 * we don't use mmap to write to file as it causes excessive traffic to NFS file system.
 * */
async_t* async_init(file_t* fp, const size_t size, const uint32_t magic,
			  const char* str, const void* p, const uint64_t nx, const uint64_t ny){
	//dbg("%s: initializing async info\n", fp->fn);
	if(!fp) return NULL;
#if 1 //use hole to reserve space for the current block.
	header_t header={magic, nx, ny, (char*)str};
	write_header(&header, fp);
	long pos=zfpos(fp);
	if(nx&&ny){//need to write after seek to create a hole.
		char zero=0;
		zfseek(fp, size*nx*ny-sizeof(char), SEEK_CUR);
		zfwrite(&zero, 1, sizeof(char), fp);
	}
#else
	long pos=writearr(fp, 0, size, magic, str, p, nx, ny);
#endif
	if(pos<0){
		warning("pos shall not be negative\n");
		return NULL;
	}
	async_t* async=mycalloc(1, struct async_t);
	if(async){
		async->fp=fp;
		async->p=p;
		async->pos=pos;//start of data block. do not modify
		async->prev=0; //length of data block written by previous async_write() call.
#if USE_ASYNC	
		async->aio.aio_fildes=fp->fd;
#endif		
	}
	return async;
}
/**
 * Wait until previous async write are finished.
 * */
#if USE_ASYNC	
static void async_sync(async_t* async, int wait){
	long ans;
	if(async->aio.aio_nbytes){//check previous result
		/*avoid using aio_fsync(). Calling aio_fsync breaks the code because if
		aio_fsync is called after aio_write, aio_return returns status of
		aio_fsync, not aio_write. aio_fsync is not available in macos.*/
		int retry=0;
		while((ans=aio_error(&async->aio))==EINPROGRESS&&wait){
			const struct aiocb* paio=&async->aio;
			//aio_suspend waits for async io to finish. it may also timeout or be intterupted
			if((ans=aio_suspend(&paio, 1, NULL))){
				dbg("%s: aio_suspend returns %ld\n", zfname(async->fp), ans);
			}
			retry++;
		}
		if(ans!=EINPROGRESS){
			if(retry>1) dbg("%s: aio synced after %d retries\n", zfname(async->fp), retry);
			if((ans=aio_return(&async->aio))!=(long)async->aio.aio_nbytes){
				dbg("%s: wrote %ld bytes out of %ld, mark eof.\n",
					zfname(async->fp), ans, async->aio.aio_nbytes);
				zferr(async->fp, 2);
			}
			async->aio.aio_nbytes=0;
		}
	}
}
#endif	

/**
 * Write to file asynchronously. offset is the end of block that is ready to write.
 * */
int async_write(async_t* async, long offset, int wait){
	if(!async||async->fp->err){
		if(!async){
			dbg("async is empty.\n");
			return 0;
		}else{
			dbg("%s: previous aio_write returned err=%d.\n", zfname(async->fp), async->fp->err);
			return -1;
		}
	}
	long ans=0;
#if USE_ASYNC
	//synchronize previous operation
	if(async->aio.aio_nbytes){
		async_sync(async, wait);
	}
	if(!wait&&async->aio.aio_nbytes){
		//dbg("%s: skip this write\n", zfname(async->fp));
		return 0;
	}
	if(offset>async->prev){//queue next.
		memset(&async->aio, 0, sizeof(struct aiocb));//macos recommends zeroing the buffer.
		async->aio.aio_fildes=async->fp->fd;
		async->aio.aio_offset=async->pos+async->prev;
		async->aio.aio_buf=(void*)(async->p+async->prev);
		async->aio.aio_nbytes=offset-async->prev;
		//dbg("%s: Async writing from %ld to %ld at offset %ld\n", zfname(async->fp), async->prev, offset, async->aio.aio_offset);
		if((ans=aio_write(&async->aio))){
			if(errno==EAGAIN){//mark as temporary failure.
				async->aio.aio_nbytes=0;
			}else{//mark as permanent failure
				zferr(async->fp, 2);
				ans=-1;
				warning("%s: aio_write failed with errno %d:%s.\n", zfname(async->fp), errno, strerror(errno));
			}
		}else{//update previous marker only if success
			async->prev=offset;
		}
	}
#else
	(void) wait;
	int nwrite=pwrite(async->fp->fd, async->p+async->prev, offset-async->prev, async->prev+async->pos);
	if(nwrite<0){
		perror("pwrite");
		warning("pwrite returns %d\n", nwrite);
		ans=-1;
	}else if(nwrite<offset-async->prev){
		async->prev+=nwrite;
	}else{
		async->prev=offset;
	}
#endif	
	return ans;
}
void async_free(async_t* async){
	if(!async) return;
#if USE_ASYNC
	async_sync(async, 1);
#endif	
	free(async);
}
