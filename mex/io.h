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
#ifndef _MEX_IO_H
#define _MEX_IO_H
#include <zlib.h>
#include <signal.h>
#if __GNUC__ && defined(__STDC_UTF_16__) && !defined(__cplusplus)
typedef int16_t char16_t;
#endif
#include <string.h>
#include <mex.h>
#include <stdint.h>
#include <errno.h>
#ifndef MWSIZE_MAX
//#if !defined(MX_API_VER) || MX_API_VER < 0x07030000 //MATLAB2016 dropped this def.
typedef unsigned int mwIndex;
#endif
//#if !MX_HAS_INTERLEAVED_COMPLEX //Matlab R2018a uses interleaved memory in 2018a and after
double *mxGetPiIsDeprecated(mxArray*);
//#endif
#define MAX(A,B) (A)>(B)?(A):(B)
typedef struct {float x; float y;} fcomplex;
typedef struct {double x; double y;} dcomplex;
#define dbg(A...) fprintf(stderr,A)
#define warning(A...)							\
    do{									\
	fprintf(stderr,"\033[01;33m%-15s:%-3d\t", __FILE__, __LINE__);	\
	fprintf(stderr,A);						\
	fprintf(stderr,"\033[00;00m");					\
    }while(0)
#define error(A...) \
    do{									\
	fprintf(stderr, "\033[01;31m%-15s:%-3d Fatal error\t",__FILE__, __LINE__); \
	fprintf(stderr, A);						\
	fprintf(stderr,"\033[00;00m");					\
	mexErrMsgTxt("Error happened\n");				\
    }while(0)

#define error2(A...) \
    do{									\
	fprintf(stderr, "\033[01;31m");					\
	fprintf(stderr, A);						\
	fprintf(stderr,"\033[00;00m");					\
	mexErrMsgTxt("Error happened\n");				\
    }while(0)

typedef struct file_t{
    void *p;
    int fd;
    int isgzip;
    int isfits;
    int eof;
}file_t;
#define M_CSP64 0x6400
#define M_DSP64  0x6401
#define M_DBL   0x6402
#define M_INT64 0x6403
#define M_CMP   0x6404 //double complex
#define M_INT32 0x6405
#define M_CSP32 0x6406
#define M_DSP32  0x6407
#define M_FLT   0x6408
#define M_ZMP   0x6409 //float complex
#define M_INT8  0x640A
#define M_INT16 0x640B

#define M_SSP64  0x6430  /*single precision float + int64 */
#define M_SSP32  0x6431  /*single precision float + int32 */
#define M_ZSP64  0x6432  /*single precision complex + int64 */
#define M_ZSP32  0x6433  /*single precision complex + int32 */

/*
#define MC_CSP  0x6410
#define MC_DSP  0x6411
#define MC_DBL  0x6412
#define MC_INT64 0x6413
#define MC_CMP 0x6414
#define MC_INT32 0x6415
*/

#define MCC_ANY  0x6421
/*#define MCC_DBL 0x6422
#define MCC_CMP 0x6424
*/
#define M_HEADER 0x6500
#define M_SKIP   0x6600

#define MAT_SP 0xFF01
#define MAT_CSP 0xFF02
#define iscell(magic) (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420)

typedef struct {
    uint32_t magic;
    uint64_t ntot;
    uint64_t ndim;
    mwSize *dims;
    char *str;
}header_t;
file_t* zfopen(const char *fn, const char *mod);
void zfclose(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
void zfwrite_dcomplex(const double* pr, const double *pi,const size_t nmemb, file_t *fp);
void zfwrite_fcomplex(const float* pr, const float *pi,const size_t nmemb, file_t *fp);
int zfread2(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
/**
   Read from the file. if in gzip mode, calls gzread, otherwise calls
   fread. follows the interface of fread.
 */
#define zfread(ptr,size,nmemb,fp)					\
    if(zfread2(ptr,size,nmemb,fp)) {					\
	fp->eof=1;							\
	warning("Error happened while reading: %s\n", strerror(errno));	\
    }
int zfeof(file_t *fp);
long zftell(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void zfwritelarr(file_t *fp, int count, ...);
void zfreadlarr(file_t* fp,int count, ...);
void write_header(const header_t *header, file_t *fp);
int read_header2(header_t *header, file_t *fp);
void read_header(header_t *header, file_t *fp);
int search_header_int(const char *header, const char *name);
double search_header_dbl(const char *header, const char *name);
#endif
