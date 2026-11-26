/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <inttypes.h>
#include "io.h"
//static void write_timestamp(file_t *fp);
/**
   Format the current local time as a string like ctime() without the trailing newline.

   @return Pointer to a static buffer; not thread safe and must not be freed.
*/
const char* myasctime(void){
	static char st[64];
	time_t a;
	time(&a);
	ctime_r(&a, st);
	st[strlen(st)-1]='\0';
	return st;
}
/**
   Test whether a file exists.

   @param[in] fn File name.
   @return 1 if the file exists, 0 otherwise.
*/
static int exist(const char* fn){
	struct stat buf;
	return !stat(fn, &buf);
}
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.

   @param[in] a String to compare.
   @param[in] b Reference string; its length sets the number of characters compared.
   @return 0 if the first strlen(b) characters of a equal b; nonzero otherwise (also when either is NULL or empty or a is shorter).
 */
static int mystrcmp(const char* a, const char* b){
	if(!a||!b) return 1;
	int la=strlen(a);
	int lb=strlen(b);
	if(la==0||lb==0||la<lb){
		return 1;
	} else{
		return strncmp(a, b, lb);
	}
}
/**
   Check the suffix of a file.

   @param[in] fn File name.
   @param[in] suffix Suffix to check.
   @return 1 if fn ends with suffix, 0 otherwise.
*/
static int check_suffix(const char* fn, const char* suffix){
	if(!fn||!suffix) return 0;
	int lfn=strlen(fn);
	int lsu=strlen(suffix);
	if(lfn<lsu) return 0;
	if(mystrcmp(fn+lfn-lsu, suffix)){
		return 0;
	} else{
		return 1;
	}
}

/**
   Make dirs recursively. like mkdir -p in bash

   @param[in,out] fn Directory path. Trailing slashes are removed in place; the string is temporarily modified during recursion.
*/
void mymkdir(char* fn){
	if(!fn) return;
	while(fn[strlen(fn)-1]=='/')
		fn[strlen(fn)-1]='\0';
	if(mkdir(fn, 0777)==-1&&errno!=EEXIST){
		perror("mkdir");
		char* tmp=strrchr(fn, '/');
		if(!tmp){
			error("Unable to mkdir '%s'\n", fn);
		}
		tmp[0]='\0';
		mymkdir(fn);

		tmp[0]='/';
		if(mkdir(fn, 0777)==-1&&errno!=EEXIST){
			error("Unable to mkdir '%s'\n", fn);
		}
	}
}
/**
   Obtain the dirname of a path.

   @param[in] fn Path.
   @return Newly allocated directory part ("." if fn has no '/'), or NULL if fn is NULL or empty. The caller frees it.
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
	return strdup(fn2);
}
/**
   Test whether fn is a symbolic link

   @param[in] fn File name.
   @return 1 if fn is a symbolic link, 0 otherwise.
*/
static int islink(const char* fn){
	if(!fn) return 0;
	struct stat buf;
	return !lstat(fn, &buf)&&S_ISLNK(buf.st_mode);
}
/**
   Process a file name for opening: expand leading ~, append .bin if no known
   suffix. For reading, try alternatives with or without .gz and .fits. For
   writing, remove an existing symbolic link or create the parent directory.

   @param[in] fn File name.
   @param[in] mod Open mode: "r", "a", or "w".
   @return Newly allocated processed file name, or NULL if the file is not found (read) or cannot be prepared (write). The caller frees it.
*/
static char* procfn(const char* fn, const char* mod){
	if(!fn){
		dbg("fn is empty\n");
		return NULL;
	}
	char* fn2;
	if(fn[0]=='~'){
		char* HOME=getenv("HOME");
		fn2=(char*)malloc(strlen(HOME)+strlen(fn)+16);
		strcpy(fn2, HOME);
		strcat(fn2, fn+1);
	} else{
		fn2=(char*)malloc(strlen(fn)+16);
		strcpy(fn2, fn);
	}
	/*If there is no recognized suffix, add .bin in the end.*/
	int nosuffix=0;
	if(!check_suffix(fn2, ".bin")&&!check_suffix(fn2, ".bin.gz")
		&&!check_suffix(fn2, ".fits")&&!check_suffix(fn2, ".fits.gz")){
		strcat(fn2, ".bin");
		nosuffix=1;
	}
	if(mod[0]=='r'||mod[0]=='a'){
		if(!exist(fn2)){/*If does not exist.*/
			if(!check_suffix(fn2, ".gz")){
				strcat(fn2, ".gz");
			} else if(check_suffix(fn, ".gz")){
			/*ended with bin.gz, change to .bin*/
				fn2[strlen(fn2)-3]='\0';
			}
			if(!exist(fn2)){
				if(nosuffix){
					fn2[strlen(fn2)-7]='\0';/*remove the added .bin, .gz*/
					strcat(fn2, ".fits");/*replace with .fits*/
					if(!exist(fn2)){
						free(fn2);
						return NULL;
					}
				} else{
					free(fn2);
					return NULL;
				}
			}
		}
	} else if(mod[0]=='w'){
		if(islink(fn2)){
			/*remove old file to avoid write over a symbolic link.*/
			if(remove(fn2)){
				dbg("Failed to remove %s\n", fn2);
				free(fn2);
				return NULL;
			}
		} else{
			char* fd=mydirname(fn2);
			if(!exist(fd)){
				mymkdir(fd);
			}
			free(fd);
		}
	} else{
		error("Invalid mode\n");
	}
	return fn2;
}
/**
   Open a file for reading or writing. Stripped down version of io.c. Files are
   gzip compressed when writing unless the name ends with .fits. When reading,
   gzip is detected from the magic number.

   @param[in] fn File name.
   @param[in] mod Open mode: "r", "w", or "a".
   @return Newly allocated file handle, or NULL on failure. Close it with zfclose().
*/
file_t* zfopen(const char* fn, const char* mod){
	char* fn2=procfn(fn, mod);
	if(!fn2){
		dbg("%s does not exist\n", fn);
		return NULL;
	}
	file_t* fp=(file_t*)calloc(1, sizeof(file_t));
	/*check fn instead of fn2. if end of .bin or .fits, disable compressing.*/
	/*Now open the file to get a fd number that we can use to lock on the
	  file.*/
	switch(mod[0]){
	case 'r':/*read only */
		fp->fd=open(fn2, O_RDONLY);
		break;
	case 'w':/*write */
	case 'a':
		fp->fd=open(fn2, O_RDWR|O_CREAT, 0600);
		if(fp->fd!=-1&&flock(fp->fd, LOCK_EX|LOCK_NB)){
			dbg("Trying to write to a file that is already opened for writing: %s\n", fn2);
			close(fp->fd);
			free(fp); free(fn2);
			return NULL;
		} else{
			if(mod[0]=='w'&&ftruncate(fp->fd, 0)){/*Need to manually truncate the file. */
				warning("Truncating %s failed\n", fn2);
			}
		}
		break;
	default:
		dbg("Unknown mod=%s\n", mod);
		free(fp); free(fn2);
		return NULL;
	}
	if(fp->fd==-1){
		dbg("Unable to open file %s\n", fn2);
		free(fp); free(fn2);
		return NULL;
	}
	if(mod[0]=='w'){
		if(check_suffix(fn2, ".fits")){
			fp->isgzip=0;
		} else{
			fp->isgzip=1;
		}
	} else{
		uint16_t magic;
		if(read(fp->fd, &magic, sizeof(uint16_t))!=sizeof(uint16_t)){
			dbg("Read magic failed\n");
			close(fp->fd); free(fp); free(fn2); return NULL;
		}
		if(magic==0x8b1f){
			fp->isgzip=1;
		} else{
			fp->isgzip=0;
		}
		lseek(fp->fd, 0, SEEK_SET);
	}
	if(fp->isgzip){
		if(!(fp->p=gzdopen(fp->fd, mod))){
			dbg("Error gzdopen for %s\n", fn2);
			close(fp->fd); free(fp); free(fn2); return NULL;
		}
	} else{
		if(!(fp->p=fdopen(fp->fd, mod))){
			dbg("Error fdopen for %s\n", fn2);
			close(fp->fd); free(fp); free(fn2); return NULL;
		}
	}
	if(check_suffix(fn2, ".fits")||check_suffix(fn2, ".fits.gz")){
		fp->isfits=1;
	}
	/*if(mod[0]=='w' && !fp->isfits){
	write_timestamp(fp);
	}*/
	free(fn2);
	return fp;
}
/**
   Close the file and free fp.

   @param[in] fp File handle; freed on return. NULL is ignored.
*/
void zfclose(file_t* fp){
	if(!fp) return;
	if(fp->isgzip){
		gzclose((gzFile)fp->p);
	} else{
		if(fclose((FILE*)fp->p)){
			perror("fclose\n");
		}
	}
	free(fp);
}
/**
   Write raw data with gzwrite or fwrite. Closes the file and aborts on failure.

   @param[in] ptr Data to write.
   @param[in] size Size of each element in bytes.
   @param[in] nmemb Number of elements.
   @param[in,out] fp File handle.
*/
static inline void zfwrite_do(const void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	if(fp->isgzip){
		if(gzwrite((gzFile)fp->p, ptr, size*nmemb)!=(long)(size*nmemb)){
			zfclose(fp);
			perror("gzwrite");
			error("write failed\n");
		}
	} else{
		if(fwrite(ptr, size, nmemb, (FILE*)fp->p)!=nmemb){
			zfclose(fp);
			perror("fwrite");
			error("writefailed\n");
		}
	}
}
/**
   Write to the file. If in gzip mode, calls gzwrite, otherwise, calls
   fwrite. Follows the interface of fwrite. For fits files, data is converted to
   big endian and written in blocks of 2880 bytes padded with zeros. Complex
   numbers must be written with size of each real component and twice the
   number of elements so that both components are byte swapped separately.

   @param[in] ptr Data to write.
   @param[in] size Size of each element in bytes.
   @param[in] nmemb Number of elements.
   @param[in,out] fp File handle.
 */
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	/*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
	if(fp->isfits&&size>1){
	/* write a block of 2880 bytes each time, with big-endianness.*/
		const int bs=2880;
		char junk[bs];
		long length=(long)(size*nmemb);
		long nb=(length+bs-1)/bs;
		char* in=(char*)ptr;
		for(long ib=0; ib<nb; ib++){
			int nd=length<bs?(int)length:bs;
			switch(size){
			case 2:
				for(int i=0; i<nd; i+=2){
					junk[i]=in[i+1];
					junk[i+1]=in[i];
				}
				break;
			case 4:
				for(int i=0; i<nd; i+=4){
					junk[i]=in[i+3];
					junk[i+1]=in[i+2];
					junk[i+2]=in[i+1];
					junk[i+3]=in[i];
				}
				break;
			case 8:
			case 16:
				for(int i=0; i<nd; i+=8){
					junk[i]=in[i+7];
					junk[i+1]=in[i+6];
					junk[i+2]=in[i+5];
					junk[i+3]=in[i+4];
					junk[i+4]=in[i+3];
					junk[i+5]=in[i+2];
					junk[i+6]=in[i+1];
					junk[i+7]=in[i];
				}
				break;
			default:
				zfclose(fp);
				error("Invalid\n");
			}
			/* use bs instead of nd to test tailing blanks*/
			in+=bs; length-=bs;
			if(length<0){
				memset(junk+nd, 0, (bs-nd)*sizeof(char));
			}
			zfwrite_do(junk, sizeof(char), bs, fp);
		}
	} else{
		zfwrite_do(ptr, size, nmemb, fp);
	}
}
/**
   Read raw data with gzread or fread.

   @param[out] ptr Buffer to receive the data.
   @param[in] size Size of each element in bytes.
   @param[in] nmemb Number of elements.
   @param[in,out] fp File handle.
   @return 0 on success, -1 on failure or end of file.
*/
static inline int zfread_do(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	if(fp->isgzip){
		return gzread((gzFile)fp->p, ptr, size*nmemb)>0?0:-1;
	} else{
		return fread(ptr, size, nmemb, (FILE*)fp->p)==nmemb?0:-1;
	}
}
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread. It does byte ordering from big endian
   to small endian in case we are reading fits file.

   @param[out] ptr Buffer to receive the data.
   @param[in] size Size of each element in bytes.
   @param[in] nmemb Number of elements.
   @param[in,out] fp File handle.
   @return 0 on success, -1 on failure or end of file.
*/
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp){
	/*a wrapper to call either fwrite or gzwrite based on flag of isgzip*/
	if(fp->isfits&&size>1){/*need to do byte swapping.*/
		const long bs=2880;
		char junk[bs];
		long length=size*nmemb;
		long nb=(length+bs-1)/bs;
		char* out=(char*)ptr;
		for(long ib=0; ib<nb; ib++){
			if(zfread_do(junk, sizeof(char), bs, fp)) return -1;
			int nd=length<bs?(int)length:(int)bs;
			switch(size){
			case 2:
				for(int i=0; i<nd; i+=2){
					out[i]=junk[i+1];
					out[i+1]=junk[i];
				}
				break;
			case 4:
				for(int i=0; i<nd; i+=4){
					out[i]=junk[i+3];
					out[i+1]=junk[i+2];
					out[i+2]=junk[i+1];
					out[i+3]=junk[i];
				}
				break;
			case 8:
			case 16:
				for(int i=0; i<nd; i+=8){
					out[i]=junk[i+7];
					out[i+1]=junk[i+6];
					out[i+2]=junk[i+5];
					out[i+3]=junk[i+4];
					out[i+4]=junk[i+3];
					out[i+5]=junk[i+2];
					out[i+6]=junk[i+1];
					out[i+7]=junk[i];
				}
				break;
			default:
				zfclose(fp);
				error("Invalid\n");
			}
			out+=bs;
			length-=bs;
		}
		return 0;
	} else{
		return zfread_do(ptr, size, nmemb, fp);
	}
}

/**
   Move the current position pointer, like fseek. For fits files, offset is
   rounded up to a multiple of 2880 bytes.

   @param[in,out] fp File handle.
   @param[in] offset Offset in bytes.
   @param[in] whence SEEK_SET, SEEK_CUR, or SEEK_END.
   @return 0 on success, -1 on failure.
*/
int zfseek(file_t* fp, long offset, int whence){
	if(fp->isfits){/*in fits, we need to offset integer block of 2880.*/
		const long bs=2880;
		long nb=(offset+bs-1)/bs;
		offset=nb*bs;
	}
	if(fp->isgzip){
		return gzseek((gzFile)fp->p, offset, whence)==-1?-1:0;
	} else{
		return fseek((FILE*)fp->p, offset, whence);
	}
}
/**
   Return the current position in the file, like ftell.

   @param[in] fp File handle.
   @return Current position in bytes (uncompressed position for gzip files), or -1 on error.
*/
long zftell(file_t* fp){
	if(fp->isgzip){
		return gztell((gzFile)fp->p);
	} else{
		return ftell((FILE*)fp->p);
	}
}
/**
   Test end of file by seeking forward by one byte.

   @param[in,out] fp File handle. The position is advanced on success.
   @return Nonzero if the seek fails, 0 otherwise.
*/
int zfeof(file_t* fp){
	return zfseek(fp, 1, SEEK_CUR);
}
/**
   Write the magic into file. Also write a dummy header to make data alignment to 8 bytes.

   @param[in] magic Magic number identifying the data type.
   @param[in,out] fp File handle. Must not be a fits file.
*/
static void write_bin_magic(uint32_t magic, file_t* fp){
	if(fp->isfits){
		zfclose(fp);
		error("fits file is not supported\n");
	}
	uint32_t magic2=M_SKIP;
	zfwrite(&magic2, sizeof(uint32_t), 1, fp);
	zfwrite(&magic, sizeof(uint32_t), 1, fp);
}
/**
   Append to the file the header to the end of the file(or rather, the
   tailer). First write magic number, then the length of the header, then the
   header, then the length of the header again, then the magic number again. The
   two set of strlen and header are used to identify the header from the end of
   the file and also to verify that what we are reading are indeed header. The
   header may be written multiple times. They will be concatenated when
   read. The header should contain key=value entries just like the configuration
   files. The entries should be separated by new line charactor.

   @param[in] keywords Header string. Nothing is written if NULL.
   @param[in,out] fp File handle. Must not be a fits file.
*/
static void write_bin_header(const char* keywords, file_t* fp){
	if(!keywords) return;
	if(fp->isfits){
		zfclose(fp);
		error("fits file is not supported\n");
	}
	uint32_t magic=M_COMMENT;
	uint64_t nlen=strlen(keywords)+1;
	/*make header 8 byte alignment.*/
	uint64_t nlen2=(nlen/8+1)*8;
	char* keywords2=(char*)calloc(nlen2, sizeof(char));
	memcpy(keywords2, keywords, nlen);
	zfwrite(&magic, sizeof(uint32_t), 1, fp);
	zfwrite(&nlen2, sizeof(uint64_t), 1, fp);
	zfwrite(keywords2, 1, nlen2, fp);
	zfwrite(&nlen2, sizeof(uint64_t), 1, fp);
	zfwrite(&magic, sizeof(uint32_t), 1, fp);
	free(keywords2);
}

/*static void write_timestamp(file_t *fp){
	if(fp->isfits) error("Not supported\n");
	char header[128];
	snprintf(header,128, "Created by write on %s\n",
		 myasctime());
	write_bin_header(header, fp);
	}*/

/**
   Obtain the current magic number. If it is a header, read it out if output of
   header is not NULL.  The header will be appended to the output header.

   @param[in,out] fp File handle. Must not be a fits file.
   @param[in,out] pkeywords If not NULL, headers read are appended to *pkeywords (allocated if NULL). The caller frees it.
   @return The magic number of the next data block, or 0 at end of file.
*/
uint32_t read_bin_magic(file_t* fp, char** pkeywords){
	uint32_t magic, magic2;
	uint64_t nlen, nlen2;
	if(fp->isfits){
		zfclose(fp);
		error("fits file is not supported\n");
	}
	while(1){
	/*read the magic number.*/
		if(zfread2(&magic, sizeof(uint32_t), 1, fp)) return 0;
		/*If it is header, read or skip it.*/
		if(magic==M_SKIP){
			continue;
		} else if(magic==M_COMMENT){
			zfread(&nlen, sizeof(uint64_t), 1, fp);
			if(nlen>0){
				if(pkeywords){
					char keywords[nlen];
					zfread(keywords, 1, nlen, fp);
					keywords[nlen-1]='\0'; /*make sure it is NULL terminated.*/
					if(*pkeywords){
						*pkeywords=(char*)realloc(*pkeywords, sizeof(char)*(((*pkeywords)?strlen(*pkeywords):0)+strlen(keywords)+1));
						strncat(*pkeywords, keywords, nlen);
					} else{
						*pkeywords=strdup(keywords);
					}
				} else{
					zfseek(fp, nlen, SEEK_CUR);
				}
			}
			zfread(&nlen2, sizeof(uint64_t), 1, fp);
			zfread(&magic2, sizeof(uint32_t), 1, fp);
			if(magic!=magic2||nlen!=nlen2){
				zfclose(fp);
				error("Header verification failed: %u<>%u, %lu<>%lu\n", magic, magic2, (unsigned long)nlen, (unsigned long)nlen2);
			}
		} else{ /*otherwise return the magic number*/
			return magic;
		}
	}/*while*/
}
/**
   Write fits header.

   @param[in,out] fp File handle. fp->isfits is incremented after the primary HDU.
   @param[in] str Extra header that will be put in fits comment. May be NULL.
   @param[in] magic Magic number identifying the data type.
   @param[in] ndim Number of dimensions.
   @param[in] dims Dimensions of length ndim.
*/
static void
write_fits_header(file_t* fp, const char* str, uint32_t magic, uint64_t ndim, mwSize* dims){
	int bitpix;
	switch(magic){
	case M_FLT:
		bitpix=-32;
		break;
	case M_DBL:
	case M_DSP64:
		bitpix=-64;
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
	case M_UINT8:
		bitpix=8;
		break;
	default:
		bitpix=0;
		zfclose(fp);
		error("Data type magic=%x is not yet supported.\n", magic);
	}
	const int nh=36;
	char header[nh][80];
	memset(header, ' ', sizeof(char)*36*80);
	int hc=0;
	const int primary=(fp->isfits==1);
	if(primary){
		snprintf(header[hc], 80, "%-8s= %20s", "SIMPLE", "T");    header[hc][30]=' '; hc++;
		fp->isfits++;
	} else{
		snprintf(header[hc], 80, "%-8s= %s", "XTENSION", "'IMAGE   '");    header[hc][20]=' '; hc++;
	}
	snprintf(header[hc], 80, "%-8s= %20d", "BITPIX", bitpix); header[hc][30]=' '; hc++;
	snprintf(header[hc], 80, "%-8s= %20ld", "NAXIS", (long)ndim);   header[hc][30]=' '; hc++;

#define FLUSH_OUT /*write the block and reset */	\
	if(hc==nh){					\
	    zfwrite(header, sizeof(char), 36*80, fp);	\
	    memset(header, ' ', sizeof(char)*36*80);	\
	    hc=0;					\
	}
	for(unsigned long i=0; i<ndim; i++){
		FLUSH_OUT;
		snprintf(header[hc], 80, "%-5s%-3lu= %20ld", "NAXIS", i+1, (long)(dims[i])); header[hc][30]=' '; hc++;
	}
	if(primary){//Write the extend keyword which does not mendate extension to be present.
		snprintf(header[hc], 80, "%-8s= %20s", "EXTEND", "T");    header[hc][30]=' '; hc++;
	} else{
		snprintf(header[hc], 80, "%-8s= %20s", "PCOUNT", "0");    header[hc][30]=' '; hc++;
		snprintf(header[hc], 80, "%-8s= %20s", "GCOUNT", "1");    header[hc][30]=' '; hc++;
	}
	if(str){/*We wrap our header in COMMENT section to avoid dealing with name
		 * length limit*/
		const char* str2=str+strlen(str);
		while(str<str2){
			const char* nl=strchr(str, '\n');
			int length;
			if(nl){
				length=nl-str+1;
			} else{
				length=strlen(str);
			}
			if(length>70) length=70;
			FLUSH_OUT;
			strncpy(header[hc], "COMMENT   ", 11);
			strncpy(header[hc]+10, str, length);
			if(nl){
				header[hc][10+length-1]=';';
			}
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
   Read fits header of the next HDU.

   @param[in,out] fp File handle.
   @param[in,out] str Non-standard keywords and comments are appended to *str (allocated if NULL). The caller frees it.
   @param[out] magic Data type magic number.
   @param[out] ndim Number of dimensions.
   @param[out] dims Dimensions; newly allocated when ndim>0. The caller frees it.
   @return 0 on success, -1 at end of file or on error.
 */
int read_fits_header(file_t* fp, char** str, uint32_t* magic, uint64_t* ndim, mwSize** dims){
	char line[82];//extra space for \n and \0
	int end=0;
	int page=0;
	int bitpix=0;
	while(!end){
		int start=0;
		if(page==0){
			/*First read mandatory fits headers.*/
			if(zfread2(line, 1, 80, fp)) return -1;
			start++;
			line[80]='\0';
			if(strncmp(line, "SIMPLE", 6)&&strncmp(line, "XTENSION= 'IMAGE", 16)){
				dbg("Garbage in fits file at %ld:\n", zftell(fp));
				dbg("%s\n", line);
				return -1;
			}
			zfread(line, 1, 80, fp); line[80]='\0'; start++;
			if(sscanf(line+10, "%20d", &bitpix)!=1){
				zfclose(fp);
				error("Unable to determine bitpix\n");
			}
			zfread(line, 1, 80, fp); line[80]='\0'; start++;
			if(sscanf(line+10, "%" SCNu64, ndim)!=1){
				zfclose(fp);
				error("Unable to determine naxis\n");
			}
			if(*ndim>0){
				*dims=(mwSize*)calloc(MAX(2, *ndim), sizeof(mwSize));
				for(uint64_t idim=0; idim<*ndim; idim++){
					do{
					//skip illegal lines.
						zfread(line, 1, 80, fp); line[80]='\0'; start++;
					} while(strncmp(line, "NAXIS", 5));
					if(sscanf(line+10, "%20lu", &((*dims)[idim]))!=1){
						zfclose(fp);
						error("Unable to determine nx\n");
					}
				}
			}
		}
		for(int i=start; i<36; i++){
			zfread(line, 1, 80, fp); line[80]='\0';
			if(!strncmp(line, "END", 3)){
				end=1;
			} else{
				char* hh=line;
				int length=80;
				int newline=1;
				if(!strncmp(line, "EXTEND", 6)){
					continue;
				} else if(!strncmp(line, "PCOUNT", 6)){
					continue;
				} else if(!strncmp(line, "GCOUNT", 6)){
					continue;
				} else if(!strncmp(line, "COMMENT", 7)){
					hh=line+10;
					length-=10;
					newline=0;
				}
				//Remove trailing space.
				for(int j=length-1; j>=0; j--){
					if(isspace((int)hh[j])){
						hh[j]='\0';
						length--;
					} else{
						if(newline){
							hh[j+1]='\n';
							hh[j+2]='\0';
						}
						break;
					}
				}
				if(length>0){
					if(*str){
						*str=(char*)realloc(*str, (strlen(*str)+length+1+newline)*sizeof(char));
					} else{
						*str=(char*)malloc((length+1+newline)*sizeof(char)); (*str)[0]='\0';
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
	case 64:
		*magic=M_INT64;
		break;
	case 32:
		*magic=M_INT32;
		break;
	case 16:
		*magic=M_INT16;
		break;
	case 8:
		*magic=M_UINT8;
		break;
	default:
		zfclose(fp);
		error("bitpix=%d is not yet handled.\n", bitpix);
	}
	return 0;
}

/**
   Write multiple long numbers into the file. To write three numbers, a, b, c,
   call with zfwritelarr(fp, 3, &a, &b, &c);

   @param[in,out] fp File handle.
   @param[in] count Number of uint64_t pointers that follow.
*/
void zfwritelarr(file_t* fp, int count, ...){
	va_list ap;
	int i;
	va_start(ap, count);              /*Initialize the argument list. */
	for(i=0; i<count; i++){
		uint64_t* addr=va_arg(ap, uint64_t*);  /*Get the next argument value.   */
		zfwrite(addr, sizeof(uint64_t), 1, fp);
	}
	va_end(ap);                       /* Clean up.  */
}
/**
   Read multiple long numbers from the file. To read three numbers, a, b, c,
   call with zfreadlarr(fp, 3, &a, &b, &c);

   @param[in,out] fp File handle.
   @param[in] count Number of uint64_t pointers that follow.
 */
void zfreadlarr(file_t* fp, int count, ...){
	va_list ap;
	int i;
	va_start(ap, count);              /*Initialize the argument list. */
	for(i=0; i<count; i++){
		uint64_t* addr=va_arg(ap, uint64_t*);  /*Get the next argument value.   */
		zfread(addr, sizeof(uint64_t), 1, fp);
	}
	va_end(ap);                       /* Clean up.  */
}
/**
   A unified header writing routine for .bin and .fits files. It write the array
   information and string header if any.

   @param[in] header Header with magic, dimensions, and string header.
   @param[in,out] fp File handle.
*/
void write_header(const header_t* header, file_t* fp){
	if(fp->isfits){
		if(!iscell(header->magic)){
			write_fits_header(fp, header->str, header->magic, header->ndim, header->dims);
		}
	} else{
		if(header){
			write_bin_header(header->str, fp);
		}
		write_bin_magic(header->magic, fp);
		if(header->ndim>2){
			warning("Third dimension and more will not be saved.\n");
		}
		if(header->dims){
			zfwrite(header->dims, sizeof(uint64_t), 2, fp);
		} else{
			uint64_t temp[2]={0,0};
			zfwrite(temp, sizeof(uint64_t), 2, fp);
		}
	}
}
/**
   A unified header reading routine for .bin and .fits files. It read the array
   information and string header if any.

   @param[out] header Receives magic, ndim, dims, ntot, and str. header->str is reset and newly allocated; dims is allocated if NULL. The caller frees them.
   @param[in,out] fp File handle.
   @return 0 on success (magic is 0 at end of a bin file), -1 on fits read error.
*/
int read_header2(header_t* header, file_t* fp){
	int ans;
	header->str=NULL;
	if(fp->isfits){
		ans=read_fits_header(fp, &header->str, &header->magic, &header->ndim, &header->dims);
	} else{
		header->magic=read_bin_magic(fp, &header->str);
		if(header->magic==0){
			ans=0;
			header->ndim=0;//EOF. do not reuse ndim from a previous call
		} else{
			ans=0;
			header->ndim=2;
			if(!header->dims) header->dims=(mwSize*)calloc(2, sizeof(mwSize));
			zfread(header->dims, sizeof(mwSize), 2, fp);
		}
	}
	header->ntot=header->ndim>0?1:0;
	for(uint64_t idim=0; idim<header->ndim; idim++){
		header->ntot*=header->dims[idim];
	}
	return ans;
}
/**
   calls read_header2 and abort if error happens.

   @param[out] header Receives the header. See read_header2().
   @param[in,out] fp File handle.
*/
void read_header(header_t* header, file_t* fp){
	if(read_header2(header, fp)){
		zfclose(fp);
		error("read_header failed\n");
	}
}
/**
   Parse an integer from header str

   @param[in] str Header string.
   @param[in] name Keyword including the separator, e.g. "Sparse=".
   @return The integer following name, or 0 if not found or str or name is NULL.
 */
int search_keyword_int(const char* str, const char* name){
	if(!str||!name) return 0;
	const char* tmp=strstr(str, name);
	if(tmp){
		return strtol(tmp+strlen(name), NULL, 10);
	} else{
		return 0;
	}
}
/**
   Parse a floating point number from header str

   @param[in] str Header string.
   @param[in] name Keyword including the separator.
   @return The number following name, or 0 if not found.
*/
double search_keyword_dbl(const char* str, const char* name){
	const char* tmp=strstr(str, name);
	if(tmp){
		return strtod(tmp+strlen(name), NULL);
	} else{
		return 0;
	}
}
