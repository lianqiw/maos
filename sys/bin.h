/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_BIN_H_
#define AOS_BIN_H_
#include <sys/time.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <pthread.h>
#include "common.h"
/**
   \file bin.h

   Defines our custom file format .bin that may be gzipped and the basic IO
   functions. All file read/write operators are through functions in this
   file. The routines can also operate on .fits files.
 */
#define IO_TIMMING 0
/*The definitions here should not be changed once set for backward/foreward compatibility. */
#define M_CSP64  0x6400  /*sparse complex */
#define M_SP64   0x6401  /*sparse */
#define M_DBL    0x6402  /*double */
#define M_INT64  0x6403  /*int 64 array */
#define M_CMP    0x6404  /*complex array */
#define M_INT32  0x6405  /*int 32 array */
#define M_CSP32  0x6406  /*sparse complex 32 bit */
#define M_SP32   0x6407  /*sparse 32 bit */
#define M_FLT    0x6408  /*single precision float. */
#define M_ZMP    0x6409  /*single precision complex */
#define M_INT8   0x640A  /*int 8  array */
#define M_INT16  0x640B  /*int 16 array */

#define M_SSP64  0x6430  /*single precision float + int64 */
#define M_SSP32  0x6431  /*single precision float + int 32 */
#define M_ZSP64  0x6432  /*single precision complex + int64 */
#define M_ZSP32  0x6433  /*single precision complex + int32 */

/*The individual MC_* and MCC_* have been deprecated. Use MCC_ANY for cell arrays of any type */
#define MCC_ANY  0x6421  /*cell of any thing */
#define M_HEADER 0x6500  /*header. */
#define M_SKIP   0x6600  /*the padding of magic number. */
#define iscell(magic) (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420)

#define USE_ZLIB_H 0
#if USE_ZLIB_H
#include <zlib.h> /*zlib.h in ubuntu sucks */
#else
typedef void* voidp;
char* procfn(const char *fn, const char *mod,const int defaultgzip);
voidp gzopen(const char *path, const char *mod);
voidp gzdopen(int fd, const char *mod);
long  gztell(voidp gzfile);
int   gzclose(voidp gzfile);
int   gzwrite(voidp gzfile, const void* buf, unsigned len);
long  gztell(voidp gzfile);
int   gzread(voidp gzfile, voidp buf, unsigned len);
int   gzseek(voidp file, long offset, int whence);
int   gzrewind(voidp file);
int   gzflush(voidp gzfile, int flush);
#endif
typedef struct file_t file_t;
typedef struct mmap_t mmap_t;
typedef struct {
    uint32_t magic;
    uint64_t nx;
    uint64_t ny;
    char *str;
}header_t;
/*
  wrapping for standard and zlib functions that 
  open, close, write and read files, test end of file.
  Make the following function private so that other 
  routine won't call them directly to avoid conflict.
*/
extern int disable_save; //saving to disk will be disabled when set to nonzero.
/*
  The following functions takes long type integers.
*/
int zfexist(const char *format, ...); 
void zftouch(const char *format, ...);
int zfeof(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void zfrewind(file_t *fp);
file_t* zfdopen(int sock, char *mod);
file_t *zfopen(const char *fn, char *mod);
file_t *zfdopen(int fd, char *mod);
const char *zfname(file_t *fp);
int zfisfits(file_t *fp);
void zfclose(file_t *fp);
void zflush(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
const char *search_header(const char *header, const char *key);
double search_header_num(const char *header, const char *key);
double search_header_num_valid(const char *header, const char *key);
uint64_t bytes_header(const char *header);
void write_timestamp(file_t *fp);
//uint32_t read_magic(file_t *fp, char **header);
//void write_magic(uint32_t magic, file_t *fp);
void write_header(const header_t *header, file_t *fp);
int read_header2(header_t *header, file_t *fp);
void read_header(header_t *header, file_t *fp);
/**Check whether the header refers to a cell. If yes, return NULL. nx, ny are assigned to the dimension.*/
INLINE header_t *check_cell(header_t *header, long *nx, long *ny){
    header_t *headerc=NULL;
    if(iscell(header->magic)){
	*nx=header->nx;
	*ny=header->ny;
    }else{
	*nx=1;
	*ny=1;
	headerc=header;
    }
    return headerc;
}
/*
  convenient function to write multiple long numbers
*/
void zfwritelarr(file_t* fp,int count, ...);
void zfreadlarr(file_t* fp,int count, ...);

void do_write(const void *fpn, const int isfile, const size_t size, const uint32_t magic,
	      const char *header, const void *p, const uint64_t nx, const uint64_t ny);
void writedbl(const double *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writeflt(const float *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writecmp(const dcomplex *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writefcmp(const fcomplex *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writeint(const int *p, long nx, long ny, const char*format,...) CHECK_ARG(4);
void writelong(const long *p, long nx, long ny, const char*format,...) CHECK_ARG(4);
void writespint(const spint *p, long nx, long ny, const char *format,...) CHECK_ARG(4);
void readspintdata(file_t *fp, uint32_t magic, spint *out, long len);
spint *readspint(file_t *fp, long* nx, long* ny);
void mmap_unref(struct mmap_t *in);
struct mmap_t *mmap_new(int fd, void *p, long n);
mmap_t*mmap_ref(mmap_t *in);
#endif
