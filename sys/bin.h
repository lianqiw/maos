/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <limits.h>
#include <stdint.h>
#include <pthread.h>
#include "common.h"
/**
   \file bin.h

   Defines our custom file format .bin that may be gzipped and the basic IO
   functions. All file read/write operators are through functions in this
   file. The routines can also operate on .fits files.
 */
/*The definitions here should not be changed once set for backward/foreward compatibility. */
typedef enum M_ID{
  M_0      =0x0,     /*undefined*/
	M_CSP64  =0x6400,  /*sparse complex */
	M_DSP64  =0x6401,  /*sparse double*/
	M_DBL    =0x6402,  /*double */
	M_INT64  =0x6403,  /*int 64 array */
	M_CMP    =0x6404,  /*complex array */
	M_INT32  =0x6405,  /*int 32 array */
	M_CSP32  =0x6406,  /*sparse complex with 32 bit integer */
	M_DSP32  =0x6407,  /*sparse double  with 32 bit integer*/
	M_FLT    =0x6408,  /*single precision float. */
	M_ZMP    =0x6409,  /*single precision complex */
	M_INT8   =0x640A,  /*int 8  array */
	M_INT16  =0x640B,  /*int 16 array */

  //From 0x6410 to 0x6424 are cell arrays that are not longer used except 0x6421
	MCC_ANY  =0x6421,  /*cell of any thing */

	M_SSP64  =0x6430,  /*single precision float + int64 */
	M_SSP32  =0x6431,  /*single precision float + int32 */
	M_ZSP64  =0x6432,  /*single precision complex + int64 */
	M_ZSP32  =0x6433,  /*single precision complex + int32 */

	M_COMMENT=0x6500,  /*data comments. */
	M_SKIP   =0x6600,   /*the padding of magic number. */

	M_MAP32     =0x016408, /*map_t, compatible with M_FLT*/
  M_RECTMAP32 =0x026408, /*map_t, compatible with M_FLT*/
  M_LOC32     =0x036408, /*loc_t with float data*/

  M_MAP64     =0x016402, /*map_t, compatible with M_DBL*/
  M_RECTMAP64 =0x026402, /*map_t, compatible with M_DBL*/
  M_LOC64     =0x036402, /*loc_t with double data*/
  M_EOD       =0x64FF    /*Indicate end of valid data */
}M_ID;
#if LONG_MAX==2147483647L //long is 32 bit
#define M_LONG M_INT32
#elif LONG_MAX==9223372036854775807L
#define M_LONG M_INT64 //long is 64 bit
#else
#error "Unknown long size"
#endif
#define M_INT M_INT32
#define M_STR M_INT8

#if defined(DLONG)
typedef long spint;
#define M_SPINT M_LONG
#define M_DSP M_DSP64
#define M_CSP M_CSP64
#define M_SSP M_SSP64
#define M_ZSP M_ZSP64
#else
typedef int spint;
#define M_SPINT M_INT
#define M_DSP M_DSP32
#define M_CSP M_CSP32
#define M_SSP M_SSP32
#define M_ZSP M_ZSP32
#endif


typedef struct file_t file_t;
typedef struct {
  union{
    M_ID id;//this must be the first element because we cast header_t to uint32_t.
    uint32_t magic;
  };
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
extern int disable_save; ///saving to disk will be disabled when set to nonzero.
/*
  The following functions takes long type integers.
*/
int zfexist(const char* format,...) CHECK_ARG(1);
void zftouch(const char *format,...) CHECK_ARG(1);
int zfeof(file_t *fp);
long zfpos(file_t *fp);
long zfseek(file_t *fp, long offset, int whence);
void zfrewind(file_t *fp);
long zflen(file_t *fp);
file_t *zfopen(const char *fn, const char *mod);
file_t *zfdopen(int fd);
const char *zfname(file_t *fp);
int zfisfits(file_t *fp);
void zfclose(file_t *fp);
int zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
int zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
int zferrno(file_t* fp);
uint64_t bytes_header(const char *header);
void write_timestamp(file_t *fp);
void write_header(const header_t *header, file_t *fp);
int read_header(header_t *header, file_t *fp);
/**Check whether the header refers to a cell. If yes, return NULL. nx, ny are assigned to the dimension.*/
long writearr(const void *fpn, const int isfile, const size_t size, const uint32_t magic,
	      const char *header, const void *p, const uint64_t nx, const uint64_t ny);
typedef struct mem_t mem_t;
struct mem_t *mem_new(void *p)__attribute__((warn_unused_result));
void mem_unref(mem_t **in);
mem_t*mem_ref(mem_t *in)__attribute__((warn_unused_result));
void mem_replace(mem_t *in, void *p);
int mem_nref(const mem_t *in);
void* mem_p(const mem_t *in);
mem_t* mmap_open(const char *fn, size_t msize, int rw);
void mmap_write_header(char **p0, uint32_t magic, long nx, long ny, const char *header);
void mmap_read_header(char **p0, uint32_t *magic, long *nx, long *ny, const char **header0);
#define IS_SHM(name) (name && ((name[0]=='/' && !strchr(name+1, '/')) || !mystrcmp(name, "/shm")))

//For parallel/async IO.
typedef struct async_t async_t;
async_t* async_init(file_t* fp, const size_t size, const uint32_t magic,
      const char* str, const void* p, const uint64_t nx, const uint64_t ny);
void async_write(async_t* async, long offset, int wait);
void async_free(async_t* async);
#endif
