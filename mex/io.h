#ifndef _MEX_IO_H
#define _MEX_IO_H
#include <complex.h>
#include <zlib.h>
#include <signal.h>
#include <mex.h>
#include <matrix.h>
#include <stdint.h>
#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef unsigned int mwIndex;
#endif
typedef double __complex__ dcomplex;
typedef float __complex__ fcomplex;
#define info(A...) fprintf(stderr,A)
#define warning(A...) {							\
	fprintf(stderr,"\033[01;33m%-15s:%-3d\t", __FILE__, __LINE__);	\
	fprintf(stderr,A);fprintf(stderr,"\033[00;00m");}
#define error(A...) {							\
	fprintf(stderr, "\033[01;31m%-15s:%-3d Fatal error\t",__FILE__, __LINE__); \
	fprintf(stderr, A); fprintf(stderr,"\033[00;00m");		\
	mexErrMsgTxt("Error happend\n");}
typedef struct file_t{
    void *p;
    int isgzip;
    int eof;
}file_t;
#define M_CSP64 0x6400
#define M_SP64  0x6401
#define M_DBL   0x6402
#define M_INT64 0x6403
#define M_CMP   0x6404
#define M_INT32 0x6405
#define M_CSP32 0x6406
#define M_SP32  0x6407
#define M_FLT   0x6408
#define M_FCMP  0x6409
#define MC_CSP  0x6410
#define MC_SP   0x6411
#define MC_DBL  0x6412
#define MC_INT64 0x6413
#define MC_CMP 0x6414
#define MC_INT32 0x6415
#define MCC_ANY  0x6421
#define MCC_DBL 0x6422
#define MCC_CMP 0x6424
#define M_HEADER 0x6500
#define M_SKIP   0x6600

#define MAT_SP 0xFF01
#define MAT_CSP 0xFF02
const char *myasctime(void);
file_t* zfopen(const char *fn, char *mod);
void zfclose(file_t *fp);
void zfwrite(const void* ptr, const size_t size, const size_t nmemb, file_t *fp);
void zfwrite_complex(const double* pr, const double *pi,const size_t nmemb, file_t *fp);
void zfread(void* ptr, const size_t size, const size_t nmemb, file_t* fp);
int test_eof(file_t *fp);
int zfseek(file_t *fp, long offset, int whence);
void write_header(const char *header, file_t *fp);
void write_timestamp(file_t *fp);
uint32_t read_magic(file_t *fp, char **header);
void write_magic(uint32_t magic, file_t *fp);
#endif
