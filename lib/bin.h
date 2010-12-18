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

#ifndef AOS_BIN_H_
#define AOS_BIN_H_
#include <sys/time.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <pthread.h>
#include "common.h"
#define IO_TIMMING 0
//The definitions here should not be changed once set for backward/foreward compatibility.
#define M_CSP64  0x6400  //sparse complex
#define M_SP64   0x6401  //sparse
#define M_DBL    0x6402  //double
#define M_INT64  0x6403  //int 64 array
#define M_CMP    0x6404  //complex array
#define M_INT32  0x6405  //int 32 array
#define M_CSP32  0x6406  //sparse complex 32 bit
#define M_SP32   0x6407  //sparse 32 bit
//The individual MC_* and MCC_* have been deprecated. Use MCC_ANY for cell arrays of any type
/*
#define MC_CSP   0x6410  //complex sparse cell
#define MC_SP    0x6411  //sparse cell.
#define MC_DBL   0x6412  //double cell
#define MC_INT64 0x6413  //int64 cell
#define MC_CMP   0x6414  //complex cell
#define MC_INT32 0x6415  //int32 cell

#define MCC_DBL  0x6422  //cell of dcell
#define MCC_CMP  0x6424  //cell of ccell
*/
#define MCC_ANY  0x6421  //cell of any thing
#define M_HEADER 0x6500  //header.

#define iscell(magic) (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420)

#define USE_ZLIB_H 0
#if USE_ZLIB_H
#include <zlib.h> //zlib.h in ubuntu sucks
#else
typedef void* voidp;
char* procfn(const char *fn, const char *mod,const int defaultgzip);
voidp gzopen(const char *path, const char *mod);
long  gztell(voidp gzfile);
int   gzclose(voidp gzfile);
int   gzwrite(voidp gzfile, const void* buf, unsigned len);
long  gztell(voidp gzfile);
int   gzread(voidp gzfile, voidp buf, unsigned len);
int   gzseek(voidp file, long offset, int whence);
int   gzrewind(voidp file);
int   gzflush(voidp gzfile, int flush);
#endif
/**
   contains information about opened files.
*/
typedef struct file_t{
    void *p;
    int isgzip;
#if USE_PTHREAD > 0
    pthread_mutex_t lock;
#endif
    char *fn;
#if IO_TIMMING == 1
    struct timeval tv1;
#endif
}file_t;
/*
  wrapping for standard and zlib functions that 
  open, close, write and read files, test end of file.
  Make the following function private so that other 
  routine won't call them directly to avoid conflict.
*/

/*
  The following functions takes long type integers.
*/
int zfeof(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void zfrewind(file_t *fp);
file_t* zfopen(const char *fn, char *mod);
void zfclose(file_t *fp);
void zflush(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
void write_header_end(const char *header,const char *format,...);
char *read_header_end(const char *format,...);
const char *search_header(const char *header, const char *key);
double search_header_num(const char *header, const char *key);
void write_header(const char *header, file_t *fp);
void write_timestamp(file_t *fp);
uint32_t read_magic(file_t *fp, char **header);
/*
  convenient function to write multiple long numbers
*/
void zfwritelarr(file_t* fp,int count, ...);
void zfreadlarr(file_t* fp,int count, ...);

void do_write(const void *fpn, const int isfile, const size_t size, const uint32_t magic, 
	      const void *p, const uint64_t nx, const uint64_t ny);
void do_read(const void *fpn, const int isfile, const size_t size, const uint32_t magicin, 
	     void **p, uint64_t *nx, uint64_t *ny);

void writedbl(const double *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writecmp(const dcomplex *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void readdbl(double **p, long* nx, long* ny, const char* format,...) CHECK_ARG(4);
void writeint64(const int64_t *p, long nx, long ny, const char* format,...) CHECK_ARG(4);
void writeint32(const int32_t *p, long nx, long ny, const char *format,...) CHECK_ARG(4);
void readint64(int64_t **p, long *nx, long *ny, const char* format,...) CHECK_ARG(4);
void readint32(int32_t **p, long *nx, long *ny, const char* format,...) CHECK_ARG(4);
#endif
