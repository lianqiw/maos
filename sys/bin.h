/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define M_CSP64  0x6400  /*sparse complex */
#define M_DSP64  0x6401  /*sparse double*/
#define M_DBL    0x6402  /*double */
#define M_INT64  0x6403  /*int 64 array */
#define M_CMP    0x6404  /*complex array */
#define M_INT32  0x6405  /*int 32 array */
#define M_CSP32  0x6406  /*sparse complex with 32 bit integer */
#define M_DSP32  0x6407  /*sparse double  with 32 bit integer*/
#define M_FLT    0x6408  /*single precision float. */
#define M_ZMP    0x6409  /*single precision complex */
#define M_INT8   0x640A  /*int 8  array */
#define M_INT16  0x640B  /*int 16 array */

//From 0x6410 to 0x6424 are cell arrays that are not longer used except 0x6421
#define MCC_ANY  0x6421  /*cell of any thing */

#define M_SSP64  0x6430  /*single precision float + int64 */
#define M_SSP32  0x6431  /*single precision float + int32 */
#define M_ZSP64  0x6432  /*single precision complex + int64 */
#define M_ZSP32  0x6433  /*single precision complex + int32 */

#define M_HEADER 0x6500  /*header. */
#define M_SKIP   0x6600  /*the padding of magic number. */

#if LONG_MAX==2147483647L //long is 32 bit
#define M_LONG M_INT32
#elif LONG_MAX==9223372036854775807L
#define M_LONG M_INT64 //long is 64 bit
#else
#error "Unknown long size"
#endif
#define M_INT M_INT32
#define M_STR M_INT8

#ifdef DLONG
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

#define USE_ZLIB_H 0
#if USE_ZLIB_H
#include <zlib.h> /*zlib.h in ubuntu sucks */
#else
#ifdef __cplusplus
extern "C"{
#endif
typedef void* voidp;
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
const char*gzerror(voidp gzfile, int *error);
#ifdef __cplusplus
}
#endif
#endif
typedef struct file_t file_t;
typedef struct {
    uint32_t magic;//this must be the first element because we cast header_t to uint32_t.
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
int zfexist(const char *format, ...); 
void zftouch(const char *format, ...);
int zfeof(file_t *fp);
int zfpos(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void zfrewind(file_t *fp);
file_t *zfopen_try(const char *fn, const char *mod);
file_t *zfopen(const char *fn, const char *mod);
file_t *zfdopen(int fd, const char *mod);
const char *zfname(file_t *fp);
int zfisfits(file_t *fp);
void zfclose(file_t *fp);
void zflush(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, const file_t *fp);
void zfread(void* ptr, const size_t size, const size_t nmemb, const file_t* fp);
uint64_t bytes_header(const char *header);
void write_timestamp(file_t *fp);
void write_header(const header_t *header, file_t *fp);
int read_header2(header_t *header, file_t *fp);
void read_header(header_t *header, file_t *fp);
/**Check whether the header refers to a cell. If yes, return NULL. nx, ny are assigned to the dimension.*/
void writearr(const void *fpn, const int isfile, const size_t size, const uint32_t magic,
	      const char *header, const void *p, const uint64_t nx, const uint64_t ny);
typedef struct mem_t mem_t;
struct mem_t *mem_new(void *p);
void mem_unref(struct mem_t *in);
mem_t*mem_ref(mem_t *in);
void mem_replace(mem_t *in, void *p);
int mem_isref(const mem_t *in);
void* mem_p(const mem_t *in);
mem_t* mmap_open(char *fn, size_t msize, int rw);
void mmap_write_header(char **p0, uint32_t magic, long nx, long ny, const char *header);
void mmap_read_header(char **p0, uint32_t *magic, long *nx, long *ny, const char **header0);
#define IS_SHM(name) ((name[0]=='/' && !strchr(name+1, '/')) || !mystrcmp(name, "/shm"))
#endif
