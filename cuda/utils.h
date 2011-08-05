#ifndef AOS_CUDA_UTILS_H
#define AOS_CUDA_UTILS_H

#define DO(A) if((A)!=cudaSuccess) error("(cuda) %s\n", cudaGetErrorString(cudaGetLastError()));

#define CONCURRENT 1
#if CONCURRENT
#define CUDA_SYNC_STREAM				\
    while(cudaStreamQuery(stream)!=cudaSuccess){	\
	if(THREAD_RUN_ONCE){/* no jobs to do*/		\
	    cudaStreamSynchronize(stream);		\
	    break;					\
	}						\
    }
#else
#define CUDA_SYNC_STREAM cudaStreamSynchronize(stream)
#endif
#define cudaCallocHost(P,N) ({DO(cudaMallocHost(&(P),N)); DO(cudaMemset(P,0,N));})
#define cudaCalloc(P,N) ({DO(cudaMalloc(&(P),N));DO(cudaMemset(P,0,N));})
#define TO_IMPLEMENT error("Please implement")

#define TIMING 0
#if TIMING == 1
extern int nstream;
#define STREAM_NEW(stream) ({DO(cudaStreamCreate(&stream));info2("nstream=%d\n",lockadd(&nstream,1)+1);})
#define STREAM_DONE(stream) ({DO(cudaStreamDestroy(stream));info2("nstream=%d\n",lockadd(&nstream,-1)-1);})
#else
#define STREAM_NEW(stream) DO(cudaStreamCreate(&stream))
#define STREAM_DONE(stream) DO(cudaStreamDestroy(stream))
#endif
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])
#define MYSPARSE 0
/* Private to accphi.cu. Do not include in maos*/
typedef struct{
    float (*loc)[2];//in device.
    float dx;
    int nloc;
}culoc_t;

typedef struct{
    cudaArray *ca;//3D array. for layered texture
    float **p;//float array.
    float *ht;
    float *vx;
    float *vy;
    float *ox;
    float *oy;
    float *dx;
    float *iac;
    int* cubic;
    int* nx;
    int* ny;
    int nlayer;
}cumap_t;


typedef struct{
    int *p;
    int *i;
    float *x;
    int nx;
    int ny;
    int nzmax;
}cusp_t;
void gpu_info(void);
void map2gpu(map_t **source, int nps, cumap_t *dest, int type);
void sp2gpu(cusp_t *dest, dsp *src);
void gpu_calc_ptt(double *rmsout, double *coeffout, 
		  const double ipcc, const dmat *imcc,
		  const float (*restrict loc)[2], 
		  const int nloc,
		  const float *restrict phi,
		  const float *restrict amp,
		  cudaStream_t stream
		  );
void gpu_calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out, int nmod,
		     double MCC_fcp, double ht, double scale,
		     double thetax, double thetay,
		     const double ipcc, const dmat *imcc,
		     const float (*restrict loc)[2], 
		     const int nloc,
		     const float *restrict phi,
		     const float *restrict amp,
		     cudaStream_t stream);
void gpu_loc2dev(float (* restrict *dest)[2], loc_t *src);
void gpu_dbl2dev(float * restrict *dest, double *src, int n);
void gpu_dbl2flt(float * restrict *dest, double *src, int n);
void gpu_long2dev(int * restrict *dest, long *src, int n);
void gpu_spint2dev(int * restrict *dest, spint *src, int n);
void gpu_spint2int(int * restrict *dest, spint *src, int n);
void gpu_dev2dbl(double * restrict *dest, float *src, int n);
void cuspmul(float *y, cusp_t *A, float *x, float alpha,cudaStream_t stream);
void cusptmul(float *y, cusp_t *A, float *x, float alpha,cudaStream_t stream);
__global__ void fscale_do(float *v, int n, float alpha);
void gpu_writeflt(float *p, int nx, int ny, const char *format, ...);
void gpu_writeint(int *p, int nx, int ny, const char *format, ...);
#endif
