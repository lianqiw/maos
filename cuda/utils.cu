#include "utils.h"
#include "curmat.h"
#include <pthread.h>
static cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
pthread_mutex_t cufft_mutex=PTHREAD_MUTEX_INITIALIZER;
extern cusparseMatDescr_t cuspdesc;

int nstream=0;
int NG1D=64; /**<Optimum number of blocks. Twice of multi processors.*/
int NG2D=8; /**<Optimum number of blocks. Twice of multi processors.*/

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
    DO(cudaMemGetInfo(&fr, &tot));
    info2("GPU mem used %5lu MiB (%s)\n",(tot-fr)/1024/1024, msg);
}
/**
   Get available memory.
*/
size_t gpu_get_mem(void){
    size_t fr, tot;
    DO(cudaMemGetInfo(&fr, &tot));
    return fr;
}
/**
   Initialize GPU. Return 1 if success.
 */
int gpu_init(int igpu){
    if(igpu<0) return 0;
    if(igpu==INT_MAX){//automatic
	int ngpu=0;
	DO(cudaGetDeviceCount(&ngpu));
	switch(ngpu){
	case 0:
	    warning2("No GPU is available for computing\n");
	    igpu=-1;
	    break;
	case 1:{//1 device. check whethere it is emulation.
	    cudaDeviceProp prop;
	    cudaGetDeviceProperties(&prop, 0);
	    if(prop.major==9999){
		warning("The only device is in emulation mode. Won't use it\n");
		igpu=-1;
	    }else{
		igpu=0;
	    }
	}
	    break;
	default:{//There are multiple devices.
	    cudaDeviceProp prop;
	    prop.canMapHostMemory=1;
	    prop.major=2;
	    int ans=cudaChooseDevice(&igpu, &prop);
	    if(ans==cudaErrorInvalidValue){
		warning("cudaChooseDevice Failed. Will not use GPU computing.\n");
		igpu=-1;
	    }
	}
	    break;
	}
    }
    if(igpu>-1){
	int ans=cudaSetDevice(igpu);
	switch(ans){
	case cudaSuccess:{
	    cudaDeviceProp prop;
	    cudaGetDeviceProperties(&prop, igpu);
	    NG1D=prop.multiProcessorCount*2;
	    NG2D=(int)round(sqrt((double)NG1D));
	    if(0){
		int nn=16777216;
		//int nn=1024;
		curmat *temp=curnew(nn,1);
		curset(temp,1,0);
		curmat *res=curnew(1,1);
		TIC;tic;
		cursum2(res->p, temp, 0);
		cudaStreamSynchronize(0);
		toc("inn %d", nn);
		curshow(res, 0);
		CUDA_SYNC_DEVICE;
		cudaStreamSynchronize(0);
		exit(0);
	    }
	}return 1;
	case cudaErrorInvalidDevice:
	    error2("Invalid GPU device %d\n", igpu);
	    _exit(1);
	    break;
	case cudaErrorSetOnActiveProcess:
	    warning2("Error set on active process\n");
	    break;
	default:
	    error2("Unknown error\n");
	    _exit(1);
	}
    }
    gpu_print_mem("gpu init");
    return 0;
}
/**
   Clean up device.
*/
void gpu_cleanup(void){
    cudaDeviceReset();
}
/**
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use float
array. Allow multiple calling to override the data.  */
void gpu_map2dev(cumap_t **dest0, map_t **source, int nps, int type){
    if(nps==0) return;
    if(!*dest0){
	*dest0=(cumap_t*)calloc(1, sizeof(cumap_t));
    }else if((*dest0)->nlayer!=nps){
	error("Mismatch. nlayer=%d, nps=%d\n", (*dest0)->nlayer, nps);
    }
    cumap_t *dest=*dest0;
    dest->nlayer=nps;
    int nx0=source[0]->nx;
    int ny0=source[0]->ny;
    if(!dest->vx){//data is not initialized.
	if(type==1){//all layers must be same size.
	    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
	    DO(cudaMalloc3DArray(&dest->ca, &channelDesc, make_cudaExtent(nx0, ny0, nps), 
				 cudaArrayLayered));
	}else{
	    DO(cudaMallocHost(&(dest->p), nps*sizeof(float*)));
	    //memory in device.
	    for(int ips=0; ips<nps; ips++){
		DO(cudaMalloc(&(dest->p[ips]), source[ips]->nx*source[ips]->ny*sizeof(float)));
	    }
	}
        dest->vx=new float[nps];
	dest->vy=new float[nps]; 
	dest->ht=new float[nps];
	dest->ox=new float[nps];
	dest->oy=new float[nps];
	dest->dx=new float[nps];
	dest->nx=new int[nps];
	dest->ny=new int[nps];
	for(int ips=0; ips<nps; ips++){
	    if(type==1 && source[ips]->nx!=nx0 && source[ips]->ny!=ny0){
		error("Only support map_t arrays of the same size if type==1\n");
	    }
	    dest->nx[ips]=source[ips]->nx;
	    dest->ny[ips]=source[ips]->ny;
	}
    }
    
    float *tmp=NULL;
    if(type==1) tmp=new float[nx0*ny0*nps];
    
    for(int ips=0; ips<nps; ips++){
	int nx=source[ips]->nx;
	int ny=source[ips]->ny;
	dest->vx[ips]=source[ips]->vx;
	dest->vy[ips]=source[ips]->vy;
	dest->ht[ips]=source[ips]->h;
	dest->ox[ips]=source[ips]->ox;
	dest->oy[ips]=source[ips]->oy;
	dest->dx[ips]=source[ips]->dx;

	if(type==1){//cudaArray
	    for(long ix=0; ix<(long)nx0*(long)ny0; ix++){
		tmp[ips*nx0*ny0+ix]=(float)source[ips]->p[ix];
	    }
	}else{//Flat memory
	    gpu_dbl2dev(&dest->p[ips], source[ips]->p, nx*ny);
	}
    }
    if(type==1){
	struct cudaMemcpy3DParms par={0};
	par.srcPos = make_cudaPos(0,0,0);
	par.dstPos = make_cudaPos(0,0,0);
	par.dstArray=dest->ca;
	par.srcPtr = make_cudaPitchedPtr(tmp, nx0*sizeof(float), nx0, ny0);
	par.extent = make_cudaExtent(nx0, ny0, nps);
	par.kind   = cudaMemcpyHostToDevice;
	DO(cudaMemcpy3D(&par));
	delete [] tmp;
    }
    CUDA_SYNC_DEVICE;
}

/*
  Convert a host dsp array to GPU sprase array. Both are in CSC format. 
*/
void gpu_sp2dev(cusp **dest0, dsp *src){
    if(!*dest0) *dest0=(cusp*)calloc(1, sizeof(cusp));
    cusp *dest=*dest0;
    dest->nx=src->m;
    dest->ny=src->n;
    dest->nzmax=src->nzmax;
    dest->p=NULL; dest->i=NULL; dest->x=NULL;
#if MYSPARSE == 1
    gpu_spint2int(&dest->p, src->p, src->n+1);
    gpu_spint2dev(&dest->i, src->i, src->nzmax);
    gpu_dbl2dev(&dest->x, src->x, src->nzmax);
#else
    gpu_spint2dev(&dest->p, src->p, src->n+1);
    gpu_spint2dev(&dest->i, src->i, src->nzmax);
    gpu_dbl2dev(&dest->x, src->x, src->nzmax);
#endif
}
void gpu_spcell2dev(cuspcell **dest0, spcell *src){
    if(!*dest0){
	*dest0=cuspcellnew(src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	gpu_sp2dev(&(*dest0)->p[i], src->p[i]);
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
    cuspmul_do<<<A->nx/256, 256, 0, stream>>>(y,A,x,alpha);
#else
    int status=cusparseScsrmv(handle, CUSPARSE_OPERATION_TRANSPOSE, 
			      A->ny, A->nx, alpha, cuspdesc,
			      A->x, A->p, A->i, x, 1.f, y);
    if(status!=0){
	error("cusparseScsrmv failed with status %d\n", status);
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
			      A->ny, A->nx, alpha, cuspdesc,
			      A->x, A->p, A->i, x, 1.f, y);
    if(status!=0){
	error("cusparseScsrmv failed with status %d\n", status);
    }
#endif
}
__global__ static void calc_ptt_do( float *cc,
				    const float (*restrict loc)[2], 
				    const int nloc,
				    const float *restrict phi,
				    const float *restrict amp){
    __shared__ float ccb[4];//for each block.
    if(threadIdx.x<4){
	ccb[threadIdx.x]=0.f;
    }
    __syncthreads();
    float cci[4]={0.f,0.f,0.f,0.f};//for each thread
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float tmp=phi[i]*amp[i];
	cci[0]+=tmp;
	cci[1]+=tmp*loc[i][0];
	cci[2]+=tmp*loc[i][1];
	cci[3]+=tmp*phi[i];
    }
    //Add results to shared value in each block.
    atomicAdd(&ccb[0], cci[0]);
    atomicAdd(&ccb[1], cci[1]);
    atomicAdd(&ccb[2], cci[2]);
    atomicAdd(&ccb[3], cci[3]);
    __syncthreads();//Wait until all threads in this block is done.
    if(threadIdx.x<4){//This is the first thread of a block. add block result to global.
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x]);
    }
}
__global__ static void calc_ngsmod_do( float *cc,
				       const float (*restrict loc)[2], 
				       const int nloc,
				       const float *restrict phi,
				       const float *restrict amp){
    int step=blockDim.x * gridDim.x; 
    float cci[7]={0,0,0,0,0,0,0};//for each thread
    __shared__ float ccb[7];//for each block.
    if(threadIdx.x<7){
	ccb[threadIdx.x]=0.f;
    }
    __syncthreads();
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float tmp=phi[i]*amp[i];
	cci[0]+=tmp;
	cci[1]+=tmp*loc[i][0];
	cci[2]+=tmp*loc[i][1];
	cci[3]+=tmp*loc[i][0]*loc[i][0];
	cci[4]+=tmp*loc[i][1]*loc[i][1];
	cci[5]+=tmp*loc[i][0]*loc[i][1];
	cci[6]+=tmp*phi[i];
    }
    //Add results to shared value in each block.
    atomicAdd(&ccb[0], cci[0]);
    atomicAdd(&ccb[1], cci[1]);
    atomicAdd(&ccb[2], cci[2]);
    atomicAdd(&ccb[3], cci[3]);
    atomicAdd(&ccb[4], cci[4]);
    atomicAdd(&ccb[5], cci[5]);
    atomicAdd(&ccb[6], cci[6]);
    __syncthreads();//Wait until all threads in this block is done.
    if(threadIdx.x<7){//This is the first thread of a block. add block result to global.
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x]);
    }
}
/*
  Let M be the modal matrix of pistion/tip/tilt. Calculate M'*diag(amp)*phi
  where amp is the amptliude weighting.  */
void gpu_calc_ptt(double *rmsout, double *coeffout, 
		  const double ipcc, const dmat *imcc,
		  const float (*restrict loc)[2], 
		  const int nloc,
		  const float *restrict phi,
		  const float *restrict amp,
		  cudaStream_t stream
		  ){
        //sum with 16 blocks, each with 256 threads.
    float *cc;
    cudaCalloc(cc, 4*sizeof(float), stream);
    calc_ptt_do<<<16, 256, 0, stream>>>(cc, loc, nloc, phi, amp);
    CUDA_SYNC_STREAM;
    float ccb[4];
    cudaMemcpy(ccb, cc, 4*sizeof(float), cudaMemcpyDefault);
    cudaFree(cc); cc=NULL;
    double coeff[3], tot;
    coeff[0]=ccb[0]; coeff[1]=ccb[1]; coeff[2]=ccb[2]; tot=ccb[3];
    if(coeffout){
	dmulvec3(coeffout, imcc, coeff);
    }
    if(rmsout){
	double pis=ipcc*coeff[0]*coeff[0];//piston mode variance
	double ptt=dwdot3(coeff, imcc, coeff);//p/t/t mode variance.
	rmsout[0]=tot-pis;//PR
	rmsout[1]=ptt-pis;//TT
	rmsout[2]=tot-ptt;//PTTR	
    }
}
void gpu_calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out, int nmod,
		     double MCC_fcp, double ht, double scale,
		     double thetax, double thetay,
		     const double ipcc, const dmat *imcc,
		     const float (*restrict loc)[2], 
		     const int nloc,
		     const float *restrict phi,
		     const float *restrict amp,
		     cudaStream_t stream){
    float *cc;
    double tot=0;
    cudaCalloc(cc, 7*sizeof(float), stream);
    if(nmod==2){//single DM.
	calc_ptt_do<<<16,256,0,stream>>>(cc, loc, nloc, phi, amp);
    }else if(nmod==5){//AHST mode
	calc_ngsmod_do<<<16,256,0,stream>>>(cc, loc, nloc, phi, amp);
    }else{
	TO_IMPLEMENT;
    }
    CUDA_SYNC_STREAM;
    float ccb[7];
    cudaMemcpy(ccb, cc, 7*sizeof(float), cudaMemcpyDefault);
    cudaFree(cc); 
    tot=ccb[nmod==2?3:6];
    
    double coeff[6];
    coeff[0]=ccb[0]; coeff[1]=ccb[1]; 
    coeff[2]=ccb[2]; coeff[3]=ccb[3];
    coeff[4]=ccb[4]; coeff[5]=ccb[5];
    
    if(pttrcoeff_out){
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, imcc, coeff, 1);
    }
    if(pttr_out){
	//compute TT removed wavefront variance as a side product
	double pis=ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, imcc, coeff);
	pttr_out[0]=tot-pis;//PR
	pttr_out[1]=ptt-pis;//TT
	pttr_out[2]=tot-ptt;//PTTR
    }
    //don't use +=. need locking
    ngsmod_out[0]=coeff[1];
    ngsmod_out[1]=coeff[2];
    const double scale1=1.-scale;
    if(nmod==5){
	ngsmod_out[2]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
		       -2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	ngsmod_out[3]=(scale1*(coeff[3]-coeff[4])
		       -2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
	ngsmod_out[4]=(scale1*(coeff[5])
		       -scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
    }
}
/**
   Convert a source loc_t to device memory.
*/
void gpu_loc2dev(float (* restrict *dest)[2], loc_t *src){
    float (*tmp)[2]=(float(*)[2])malloc(src->nloc*2*sizeof(float));
    for(int iloc=0; iloc<src->nloc; iloc++){
	tmp[iloc][0]=(float)src->locx[iloc];
	tmp[iloc][1]=(float)src->locy[iloc];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, src->nloc*2*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, src->nloc*2*sizeof(float),cudaMemcpyDefault));
    cudaDeviceSynchronize();
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/
void gpu_dbl2dev(float * restrict *dest, double *src, int n){
    if(!src) return;
    float *tmp=(float*)malloc(n*sizeof(float));
    for(int i=0; i<n; i++){
	tmp[i]=(float)src[i];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, n*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(float),cudaMemcpyDefault));
    cudaDeviceSynchronize();
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/
void gpu_cmp2dev(fcomplex * restrict *dest, dcomplex *src, int n){
    if(!src) return;
    fcomplex *tmp=(fcomplex*)malloc(n*sizeof(fcomplex));
    for(int i=0; i<n; i++){
	tmp[i]=(make_cuFloatComplex)(cuCreal(src[i]), cuCimag(src[i]));
    }
    if(!*dest){
	DO(cudaMalloc((fcomplex**)dest, n*sizeof(fcomplex)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(fcomplex),cudaMemcpyDefault));
    cudaDeviceSynchronize();
    free(tmp);
}
/**
   Convert dmat array to device memory.
*/
void gpu_dmat2dev(float * restrict *dest, dmat *src){
    if(!src) return;
    gpu_dbl2dev(dest, src->p, src->nx*src->ny);
}
/**
   Convert dmat array to curmat
*/
void gpu_dmat2cu(curmat *restrict *dest, dmat *src){
    if(!src){
	dzero(*dest);
	return;
    }
    if(!*dest){
	*dest=curnew(src->nx, src->ny);
    }
    gpu_dbl2dev(&(*dest)->p, src->p, src->nx*src->ny);
}
/**
   Convert dcell to curcell
*/
void gpu_dcell2cu(curcell *restrict *dest, dcell *src){
    if(!src) {
	dzero(*dest);
	return;
    }
    if(!*dest) {
	*dest=curcellnew(src->nx, src->ny);
    }else if((*dest)->nx!=src->nx || (*dest)->ny!=src->ny){
	error("Mismatch: %dx%d vs %ldx%ld\n", 
	      (*dest)->nx, (*dest)->ny, src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	gpu_dmat2cu(&(*dest)->p[i], src->p[i]);
    }
}
/**
   Convert dmat array to device memory.
*/
void gpu_cmat2dev(fcomplex * restrict *dest, cmat *src){
    if(src){
	gpu_cmp2dev(dest, (dcomplex*)src->p, src->nx*src->ny);
    }
}
/**
   Convert double array to device memory (float)
*/
void gpu_dbl2flt(float * restrict *dest, double *src, int n){
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
void gpu_long2dev(int * restrict *dest, long *src, int n){
    if(!src) return;
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(long)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyDefault));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((long)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyDefault));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
void gpu_int2dev(int *restrict *dest, int *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    DO(cudaMemcpy(*dest, src, sizeof(int)*n, cudaMemcpyDefault));
}
/**
   Convert long array to device int
*/
void gpu_spint2dev(int * restrict *dest, spint *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(spint)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyDefault));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((spint)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyDefault));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
/**
   Convert long array to device int
*/
void gpu_spint2int(int * restrict *dest, spint *src, int n){
    info("sizeof(spint)=%ld\n", sizeof(spint));

    if(!*dest){
	DO(cudaMallocHost((int**)dest, n*sizeof(int)));
    }

    for(int i=0; i<n; i++){
	(*dest)[i]=(int)src[i];
    }
}
/**
   Convert device (float) array to host double.
*/
void gpu_dev2dbl(double * restrict *dest, float *src, int n, cudaStream_t stream){
    TIC;tic;
    float *tmp=(float*)malloc(n*sizeof(float));
    DO(cudaMemcpyAsync(tmp, src, n*sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_SYNC_STREAM;
    if(!*dest){
	*dest=(double*)malloc(sizeof(double)*n);
    }
    double *restrict p=*dest;
    for(int i=0; i<n; i++){
	p[i]=tmp[i];
    }
    free(tmp);
}
/**
   scale vector by alpha.
*/
__global__ void fscale_do(float *v, int n, float alpha){
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	v[i]*=alpha;
    }
}

/*
  Write float on gpu to file
*/
void gpu_writeflt(float *p, int nx, int ny, const char *format, ...){
    format2fn;
    float *tmp=(float*)malloc(nx*ny*sizeof(float));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(float), cudaMemcpyDefault);
    cudaDeviceSynchronize();
    writeflt(tmp,nx,ny,"%s",fn);
    free(tmp);
}

/*
  Write float on gpu to file
*/
void gpu_writefcmp(fcomplex *p, int nx, int ny, const char *format, ...){
    format2fn;
    fcomplex *tmp=(fcomplex*)malloc(nx*ny*sizeof(fcomplex));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(fcomplex), cudaMemcpyDefault);
    cudaDeviceSynchronize();
    writefcmp((float complex*)tmp,nx,ny,"%s",fn);
    free(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_writeint(int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDefault);
    cudaDeviceSynchronize();
    writeint(tmp,nx,ny,"%s",fn);
    free(tmp);
}

/**
   Compute the dot product of two vectors
*/
/*
__global__ static void dot_do(float *res, const float *restrict a, const float *restrict b, const int n){
    float sumt=0;
    __shared__ float sumb;
    if(threadIdx.x == 0) sumb=0;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sumt+=a[i]*b[i];
	}
    atomicAdd(&sumb, sumt);
    __syncthreads();
    if(threadIdx.x==0){
	atomicAdd(res, sumb);
    }
}
float gpu_dot(const float *restrict a, const float *restrict b, const int n, cudaStream_t stream){
    float *res=NULL;
    cudaCallocHost(res, 1*sizeof(float), stream);
    dot_do<<<MAX(MIN(n/256, 32), 1), MIN(n, 256), 0, stream>>>(res, a, b, n);
    CUDA_SYNC_STREAM;
    float result=res[0];
    DO(cudaFreeHost(res));
    return result;
    }*/

void gpu_muv2dev(cumuv_t *out, MUV_T *in){
    if(!in->M) error("in->M should not be NULL\n");
    spcell *Mt=spcelltrans(in->M);
    gpu_spcell2dev(&(out)->Mt, Mt);
    gpu_dcell2cu(&(out)->U, in->U);
    gpu_dcell2cu(&(out)->V, in->V);
    spcellfree(Mt);
}
void gpu_cur2d(dmat **out, const curmat *in, cudaStream_t stream){
    if(!*out) *out=dnew(in->nx, in->ny);
    gpu_dev2dbl(&(*out)->p, in->p, in->nx*in->ny, stream);
}
void gpu_curcell2d(dcell **out, const curcell *in, cudaStream_t stream){
    if(!*out) *out=dcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	gpu_cur2d(&(*out)->p[i], in->p[i], stream);
    }
}

void gpu_cur2s(smat **out, const curmat *in, cudaStream_t stream){
    if(!*out) *out=snew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(float), cudaMemcpyDeviceToHost, stream));
}


void gpu_cuc2z(zmat **out, const cucmat *in, cudaStream_t stream){
    if(!*out) *out=znew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(fcomplex), cudaMemcpyDeviceToHost, stream));
}

void gpu_curcell2s(scell **out, const curcell *in, cudaStream_t stream){
    if(!*out) *out=scellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	gpu_cur2s(&(*out)->p[i], in->p[i], stream);
    }
}

void gpu_cuccell2z(zcell **out, const cuccell *in, cudaStream_t stream){
    if(!*out) *out=zcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	gpu_cuc2z(&(*out)->p[i], in->p[i], stream);
    }
}
void cellarr_cur(struct cellarr *ca, const curmat *A, cudaStream_t stream){
    smat *tmp=NULL;
    gpu_cur2s(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_smat(ca, tmp);
    sfree(tmp);
}

void cellarr_cuc(struct cellarr *ca, const cucmat *A, cudaStream_t stream){
    zmat *tmp=NULL;
    gpu_cuc2z(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_zmat(ca, tmp);
    zfree(tmp);
}
