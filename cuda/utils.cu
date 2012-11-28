/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "utils.h"
#include "curmat.h"
#include "cucmat.h"
#include <errno.h>
#include <pthread.h>
#if defined(HAS_NVML) && HAS_NVML==1
extern "C"{
    /*taken from nvml.h*/
    typedef struct nvmlDevice_st* nvmlDevice_t;
    typedef struct nvmlMemory_st 
    {
	unsigned long long total; //!< Total installed FB memory (in bytes)
	unsigned long long free; //!< Unallocated FB memory (in bytes)
	unsigned long long used; //!< Allocated FB memory (in bytes). Note that the driver/GPU always sets aside a small amount of memory for bookkeeping
    } nvmlMemory_t;
    int nvmlDeviceGetHandleByIndex(unsigned int index, nvmlDevice_t *device);
    int nvmlDeviceGetMemoryInfo(nvmlDevice_t device, nvmlMemory_t *memory);
    int nvmlInit();
    int nvmlShutdown();
}
#endif
const char *cufft_str[]={
    "success", 
    "invalid plan",
    "allocation failed",
    "",
    "invalid value",
    "internal errlr",
    "exec failed (error elsewhere caused cufft to fail)",
    "setup failed"
    "invalid size"
};

#if CUDA_VERSION < 4010
pthread_mutex_t cufft_mutex=PTHREAD_MUTEX_INITIALIZER;
#endif
int gpu_recon;/**<GPU for reconstruction*/
int NGPU=0;
int* GPUS=NULL;
int nstream=0;
cudata_t *cudata_all=NULL;/*for all GPU. */
static cusparseMatDescr_t spdesc=NULL;
#ifdef __APPLE__
pthread_key_t cudata_key;
#else
__thread cudata_t *cudata=NULL;/*for current thread and current GPU */
#endif

static __attribute((constructor)) void init(){
#ifdef __APPLE__
    pthread_key_create(&cudata_key, NULL);
#endif
    DO(cusparseCreateMatDescr(&spdesc));
    cusparseSetMatType(spdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(spdesc, CUSPARSE_INDEX_BASE_ZERO);
}
/**
   Get GPU info.
*/
void gpu_info(){
    struct cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    info("name=%s\n"
	 "TotalGlobalMem=%d\n"
	 "SharedMemPerBlock=%d\n"
	 "regsPerBlock=%d\n"
	 "warpSize=%d",
	 prop.name,
	 (int)prop.totalGlobalMem,
	 (int)prop.sharedMemPerBlock,
	 prop.regsPerBlock,
	 prop.warpSize);
}
/**
   Print memory consumption.
*/
void gpu_print_mem(const char *msg){
    size_t fr, tot;
    cudaDeviceSynchronize();
    DO(cudaMemGetInfo(&fr, &tot));
    info2("GPU (%d) mem used %'lu MB (%s)\n",(int)(cudata-cudata_all),(tot-fr)/1024/1024, msg);
}
/**
   Get available memory.
*/
long gpu_get_mem(void){
    size_t fr, tot;
    DO(cudaMemGetInfo(&fr, &tot));
    return (long)fr;
}
static int cmp_gpu_info(const long *a, const long *b){
    return (int)(b[1]-a[1]);
}
/**
   Initialize GPU. Return 1 if success.
   if gpus is not null, it is of length ngpu. gpus specifies gpu index to use.
   if gpus is null, ngpu specifies number of gpus to use. all if 0.
*/
int gpu_init(int *gpus, int ngpu){
    int ans, ngpu_tot=0;//total number of GPUs.
    if((ans=cudaGetDeviceCount(&ngpu_tot)) || ngpu_tot==0){//no GPUs available.
	info2("No GPUs available. ans=%d\n", ans);
	return 0;
    }
    NGPU=0;
    /*
      User specified exact GPUs to use. We check every entry. 
      If <0 is found, do not use any GPU.
      If >=ngpu_tot is found, skip the GPU.
      If duplicates are found, use only once.
     */
    if(gpus && ngpu>0){
	if(!GPUS) GPUS=(int*)malloc(ngpu*sizeof(int));
	for(int ig=0; ig<ngpu; ig++){
	    if(gpus[ig]<0){
		info2("CUDA is disabled by user.\n");
		free(GPUS); GPUS=NULL; 
		return 0;
	    }else{
		if(gpus[ig]>=ngpu_tot){
		    warning2("Skip GPU %d: not exist\n", gpus[ig]);
		}else{
		    GPUS[NGPU++]=gpus[ig];
		    /* Enable the following to disallow use GPUs in multiple threads
		      int j;
		    for(j=0; j<NGPU; j++){
			if(GPUS[j]==gpus[ig]){
			    warning2("Skip GPU %d: duplicated\n", gpus[ig]);
			    break;
			}
		    }
		    if(j==NGPU){
			GPUS[NGPU++]=gpus[ig];
			}*/
		}
	    }
	}
    }else{
	int repeat=1;
	if(ngpu<=0){
	    repeat=0;
	    ngpu=ngpu_tot;
	}
	GPUS=(int*)calloc(ngpu, sizeof(int));
	register_deinit(NULL, GPUS);
	/*For each GPU, query the available memory.*/
	long (*gpu_info)[2]=(long(*)[2])calloc(2*ngpu_tot, sizeof(long));
#if defined(HAS_NVML) && HAS_NVML==1
	nvmlDevice_t dev;
	nvmlMemory_t mem;
	if(nvmlInit()){
	    warning("nvml init failed\n");
	}
#endif
	for(int ig=0; ig<ngpu_tot; ig++){
	    gpu_info[ig][0]=ig;
#if defined(HAS_NVML) && HAS_NVML==1
	    if(nvmlDeviceGetHandleByIndex(ig, &dev) == 0 
	       && nvmlDeviceGetMemoryInfo(dev, &mem) == 0){
		gpu_info[ig][1]=mem.free;
	    }else
#endif
		{
		    cudaSetDevice(ig);//this allocates context.
		    gpu_info[ig][1]=gpu_get_mem();
		    //cudaDeviceReset(); We already started simulation. Do not reset.
		}
	}
#if defined(HAS_NVML) && HAS_NVML==1
	nvmlShutdown();
#endif
	/*sort so that gpus with higest memory is in the front.*/
	qsort(gpu_info, ngpu_tot, sizeof(long)*2, (int(*)(const void*, const void *))cmp_gpu_info);
	for(int igpu=0; igpu<ngpu_tot; igpu++){
	    info2("GPU %d has mem %.1f GB\n", igpu, gpu_info[igpu][1]/1024/1024/1024.);
	}
	for(int i=0, igpu=0; i<ngpu; i++, igpu++){
	    if(igpu==ngpu_tot || gpu_info[igpu][1]<500000000){
		if(repeat){
		    igpu=0; //reset to beginning.
		}else{
		    break; //stop
		}
	    }
	    GPUS[NGPU++]=(int)gpu_info[igpu][0];
	}
	free(gpu_info);
    }
    if(NGPU) {
	gpu_recon=0;/*first gpu in GPUS*/
	cudata_all=(cudata_t*)calloc(NGPU, sizeof(cudata_t));
	register_deinit(NULL, cudata_all);
	info2("Using GPU");
	for(int i=0; GPUS && i<NGPU; i++){
	    info2(" %d", GPUS[i]);
	    gpu_set(i);
	    //Reserve memory in GPU so the next maos will not pick this GPU.
	    DO(cudaMalloc(&cudata->reserve, 500000000));
	}
	info2("\n");
    }
    return NGPU;
}

/**
   Clean up device.
*/
void gpu_cleanup(void){
    for(int ig=0; ig<NGPU; ig++){
	cudaSetDevice(GPUS[ig]);
	cudaDeviceReset();
    }
}
/**
   Convert double array to device memory (float)
*/

void cp2gpu(float * restrict *dest, double *src, int n){
    if(!src) return;
    float *tmp=(float*)malloc(n*sizeof(float));
    for(int i=0; i<n; i++){
	tmp[i]=(float)src[i];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, n*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(float),cudaMemcpyHostToDevice));
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/

void cp2gpu(fcomplex * restrict *dest, dcomplex *restrict src, int n){
    if(!src) return;
    fcomplex *tmp=(fcomplex*)malloc(n*sizeof(fcomplex));
    for(int i=0; i<n; i++){
	tmp[i]=(make_cuFloatComplex)(cuCreal(src[i]), cuCimag(src[i]));
    }
    if(!*dest){
	DO(cudaMalloc((fcomplex**)dest, n*sizeof(fcomplex)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(fcomplex),cudaMemcpyHostToDevice));
    free(tmp);
}

/**
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use float
   array. Allow multiple calling to override the data.  */
void cp2gpu(cumap_t ***dest0, map_t **source, int nps){
    if(nps==0) return;
    if(!*dest0){
	*dest0=new cumap_t*[nps];
	for(int ips=0; ips<nps; ips++){
	    (*dest0)[ips]=new cumap_t(source[ips]->nx, source[ips]->ny);
	}	
    }
    cumap_t **dest=*dest0;
    for(int ips=0; ips<nps; ips++){
	dest[ips]->vx=source[ips]->vx;
	dest[ips]->vy=source[ips]->vy;
	dest[ips]->ht=source[ips]->h;
	dest[ips]->ox=source[ips]->ox;
	dest[ips]->oy=source[ips]->oy;
	dest[ips]->dx=source[ips]->dx;
	int nx=source[ips]->nx;
	int ny=source[ips]->ny;
	cp2gpu(&dest[ips]->p, source[ips]->p, nx*ny);
    }
    CUDA_SYNC_DEVICE;
}

/*
  Convert a host dsp array to GPU sprase array. Both are in CSC format. 
*/
void cp2gpu(cusp **dest0, dsp *src){
    if(!src) return;
    if(!*dest0) *dest0=(cusp*)calloc(1, sizeof(cusp));
    cusp *dest=*dest0;
    dest->nx=src->m;
    dest->ny=src->n;
    dest->nzmax=src->nzmax;
    dest->p=NULL; dest->i=NULL; dest->x=NULL;
    cp2gpu(&dest->p, src->p, src->n+1);
    cp2gpu(&dest->i, src->i, src->nzmax);
    cp2gpu(&dest->x, src->x, src->nzmax);
}
void cp2gpu(cuspcell **dest0, spcell *src){
    if(!src) return;
    if(!*dest0){
	*dest0=cuspcellnew(src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest0)->p[i], src->p[i]);
    }
}
__global__ void cuspmul_do(float *y, cusp *A, float *x, float alpha){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<A->ny; i+=step){
	for(int j=A->p[i]; j<A->p[i+1]; j++){
	    atomicAdd(&y[A->i[j]], A->x[j]*x[i]*alpha);
	}
    }
}
static const char *scsrmv_err[]={
    "Success",
    "Not initialized",
    "Allocation failed",
    "Invalid value",
    "Archtecture mismatch",
    "Mapping error",
    "Execution failed",
    "Internal error",
    "Matrix type not supported"
};
/*
  y=A*x where A is sparse. x, y are vectors. Slow for GS0.
*/
void cuspmul(float *y, cusp *A, float *x, float alpha, 
#if MYSPARSE ==1
	     cudaStream_t stream
#else
	     cusparseHandle_t handle
#endif
	     ){
#if MYSPARSE ==1
    cuspmul_do<<<DIM(A->nx, 256), 0, stream>>>(y,A,x,alpha);
#else
    int status=cusparseScsrmv(handle, CUSPARSE_OPERATION_TRANSPOSE, 
			      A->ny, A->nx, alpha, spdesc,
			      A->x, A->p, A->i, x, 1.f, y);
    if(status!=0){
	error("cusparseScsrmv failed with status '%s'\n", scsrmv_err[status]);
    }
#endif
}

/*
  y=A'*x where A is sparse. x, y are vectors
*/
__global__ void cusptmul_do(float *y, int icol, cusp *A, float *x, float alpha){
    __shared__ float val;
    if(threadIdx.x==0) val=0;
    int i=blockIdx.x * blockDim.x + threadIdx.x;
    int j=i+A->p[icol];
    atomicAdd(&val, A->x[j]*x[A->i[j]]);
    if(threadIdx.x==0) y[icol]+=val*alpha;
}
/*
  Does not work right yet. Try to launch a block for each column and n items in each block.
*/
void cusptmul(float *y, cusp *A, float *x, float alpha, 
#if MYSPARSE ==1
	      cudaStream_t stream
#else
	      cusparseHandle_t handle
#endif
	      ){
#if MYSPARSE == 1
    for(int i=0; i<A->ny; i++){
	cusptmul_do<<<1, A->p[i+1]-A->p[i], 0, stream>>>(y,i,A,x,alpha);
    }
    warning("Not working correctly yet\n");
#else
    int status=cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, 
			      A->ny, A->nx, alpha, spdesc,
			      A->x, A->p, A->i, x, 1.f, y);
    if(status!=0){
	error("cusparseScsrmv failed with status %d\n", status);
    }
#endif
}

/**
   Convert a source loc_t to device memory.
*/
void cp2gpu(float (* restrict *dest)[2], loc_t *src){
    float (*tmp)[2]=(float(*)[2])malloc(src->nloc*2*sizeof(float));
    for(int iloc=0; iloc<src->nloc; iloc++){
	tmp[iloc][0]=(float)src->locx[iloc];
	tmp[iloc][1]=(float)src->locy[iloc];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, src->nloc*2*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, src->nloc*2*sizeof(float),cudaMemcpyHostToDevice));
    free(tmp);
}

/**
   Convert dmat array to device memory.
*/
void cp2gpu(float * restrict *dest, dmat *src){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx*src->ny);
}
/**
   Convert dmat array to curmat
*/
void cp2gpu(curmat *restrict *dest, dmat *src){
    if(!src){
	curzero(*dest);
	return;
    }
    float *pdest=NULL;
    if(*dest){
	pdest=(*dest)->p;
	assert(src->nx*src->ny==(*dest)->nx*(*dest)->ny);
    }
    cp2gpu(&pdest, src->p, src->nx*src->ny);
    if(!*dest){
	*dest=curnew(src->nx, src->ny, pdest);
    }
}
void cp2gpu(curmat *restrict *dest, float *src, int nx, int ny, cudaStream_t stream){
    if(!src){
	curzero(*dest);
	return;
    }
    float *pdest=NULL;
    if(*dest){
	pdest=(*dest)->p;
	assert(nx*ny==(*dest)->nx*(*dest)->ny);
    }else{
	cudaMalloc(&pdest, nx*ny*sizeof(float));
    }
    if(stream){
	DO(cudaMemcpyAsync(pdest, src, nx*ny*sizeof(float),cudaMemcpyHostToDevice, stream));
    }else{
	DO(cudaMemcpy(pdest, src, nx*ny*sizeof(float),cudaMemcpyHostToDevice));
    }
    if(!*dest){
	*dest=curnew(nx, ny, pdest);
    }
}
/*
  convert cmat to cucmat
*/
void cp2gpu(cucmat *restrict *dest, cmat *src){
    if(!src){
	czero(*dest);
	return;
    }
    if(!*dest){
	*dest=cucnew(src->nx, src->ny);
    }else{
	assert(src->nx*src->ny==(*dest)->nx*(*dest)->ny);
    }
    cp2gpu(&(*dest)->p, (dcomplex*)src->p, (int)(src->nx*src->ny));
}
/**
   Convert dcell to curcell
*/
void cp2gpu(curcell *restrict *dest, dcell *src){
    if(!src) {
	dzero(*dest);
	return;
    }
    if(!*dest) {
	long nc=src->nx*src->ny;
	long nx[nc];
	long ny[nc];
	for(long i=0; i<nc; i++){
	    if(src->p[i]){
		nx[i]=src->p[i]->nx;
		ny[i]=src->p[i]->ny;
	    }else{
		nx[i]=0;
		ny[i]=0;
	    }
	}
	*dest=curcellnew(src->nx, src->ny, nx, ny);
    }else if((*dest)->nx!=src->nx || (*dest)->ny!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      (*dest)->nx, (*dest)->ny, src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest)->p[i], src->p[i]);
    }
}
/**
   Convert dcell to curcell
*/
void cp2gpu(cuccell *restrict *dest, ccell *src){
    if(!src) {
	dzero(*dest);
	return;
    }
    if(!*dest) {
	long nc=src->nx*src->ny;
	long nx[nc];
	long ny[nc];
	for(long i=0; i<nc; i++){
	    if(src->p[i]){
		nx[i]=src->p[i]->nx;
		ny[i]=src->p[i]->ny;
	    }else{
		nx[i]=0;
		ny[i]=0;
	    }
	}
	*dest=cuccellnew(src->nx, src->ny, nx, ny);
    }else if((*dest)->nx!=src->nx || (*dest)->ny!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      (*dest)->nx, (*dest)->ny, src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest)->p[i], src->p[i]);
    }
}
/**
   Convert dmat array to device memory.
*/
void cp2gpu(fcomplex * restrict *dest, cmat *src){
    if(src){
	cp2gpu(dest, (dcomplex*)src->p, src->nx*src->ny);
    }
}
/**
   Convert double array to device memory (float)
*/
void dbl2flt(float * restrict *dest, double *src, int n){
    if(!src) return;
    if(!*dest){
	cudaMallocHost((float**)dest, n*sizeof(float));
    }
    for(int i=0; i<n; i++){
	(*dest)[i]=(float)src[i];
    }
}
/**
   Convert long array to device int
*/
void cp2gpu(int * restrict *dest, long *src, int n){
    if(!src) return;
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(long)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((long)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
void cp2gpu(int *restrict *dest, int *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    DO(cudaMemcpy(*dest, src, sizeof(int)*n, cudaMemcpyHostToDevice));
}
/**
   Convert long array to device int
*/
void cp2gpu(int * restrict *dest, spint *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(spint)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((spint)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
/**
   Convert device (float) array and add to host double.
   dest = alpha * dest + beta *src;
*/
void cp2cpu(double * restrict *dest, double alpha, float *src, double beta, int n, 
	    cudaStream_t stream, pthread_mutex_t *mutex){
    float *tmp=(float*)malloc4async(n*sizeof(float));
    DO(cudaMemcpyAsync(tmp, src, n*sizeof(float), cudaMemcpyDeviceToHost, stream));
    if(!*dest){
	*dest=(double*)malloc(sizeof(double)*n);
    }
    double *restrict p=*dest;
    CUDA_SYNC_STREAM;
    if(mutex) LOCK(*mutex);
    for(int i=0; i<n; i++){
	p[i]=p[i]*alpha+beta*tmp[i];
    }
    if(mutex) UNLOCK(*mutex);
    free4async(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_write(float *p, int nx, int ny, const char *format, ...){
    format2fn;
    float *tmp=(float*)malloc(nx*ny*sizeof(float));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
    writeflt(tmp,nx,ny,"%s",fn);
    free(tmp);
}

/*
  Write float on gpu to file
*/
void gpu_write(fcomplex *p, int nx, int ny, const char *format, ...){
    format2fn;
    fcomplex *tmp=(fcomplex*)malloc(nx*ny*sizeof(fcomplex));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(fcomplex), cudaMemcpyDeviceToHost);
    writefcmp((float complex*)tmp,nx,ny,"%s",fn);
    free(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_write(int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    writeint(tmp,nx,ny,"%s",fn);
    free(tmp);
}

void cp2gpu(cumuv_t *out, MUV_T *in){
    if(!in->M) error("in->M should not be NULL\n");
    spcell *Mt=spcelltrans(in->M);
    cp2gpu(&(out)->Mt, Mt);
    cp2gpu(&(out)->U, in->U);
    cp2gpu(&(out)->V, in->V);
    spcellfree(Mt);
}
void cp2cpu(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream, pthread_mutex_t *mutex){
    if(!in){
	if(*out) dzero(*out);
	return;
    }
    if(!*out) {
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx*(*out)->ny==in->nx*in->ny);
    }
    cp2cpu(&(*out)->p, alpha, in->p, beta, in->nx*in->ny, stream, mutex);
}
void cp2cpu(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream, pthread_mutex_t *mutex){
    if(!in){
	if(*out) dcellzero(*out);
	return;
    }
    if(!*out) *out=dcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], alpha, in->p[i], beta, stream, mutex);
    }
}
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream){
    if(!in) {
	if(*out) szero(*out);
	return;
    }
    if(!*out) *out=snew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(float), cudaMemcpyDeviceToHost, stream));
    if(in->header) (*out)->header=strdup(in->header);
}


void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream){
    if(!in){
	if(*out) zzero(*out);
	return;
    }
    if(!*out) *out=znew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(fcomplex), cudaMemcpyDeviceToHost, stream));
    if(in->header) (*out)->header=strdup(in->header);
}

void cp2cpu(scell **out, const curcell *in, cudaStream_t stream){
    if(!in){
	if(*out) scellzero(*out);
	return;
    }
    if(!*out) *out=scellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], in->p[i], stream);
    }
}

void cp2cpu(zcell **out, const cuccell *in, cudaStream_t stream){
    if(!in){
	if(*out) zcellzero(*out);
	return;
    }
    if(!*out) *out=zcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], in->p[i], stream);
    }
}
void cellarr_cur(struct cellarr *ca, const curmat *A, cudaStream_t stream){
    smat *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_smat(ca, tmp);
    sfree(tmp);
}

void cellarr_cuc(struct cellarr *ca, const cucmat *A, cudaStream_t stream){
    zmat *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_zmat(ca, tmp);
    zfree(tmp);
}

void cellarr_curcell(struct cellarr *ca, const curcell *A, cudaStream_t stream){
    scell *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_scell(ca, tmp);
    scellfree(tmp);
}

void cellarr_cuccell(struct cellarr *ca, const cuccell *A, cudaStream_t stream){
    zcell *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_zcell(ca, tmp);
    zcellfree(tmp);
}
