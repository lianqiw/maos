extern "C"
{
#include <cuda.h>
#include "gpu.h"
#include "utils.h"
}
static cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
int nstream=0;
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
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use float array.
*/
void map2gpu(map_t **source, int nps, cumap_t *dest, int type){
    if(nps==0) return;
    if(dest->nlayer!=0 && dest->nlayer!=nps){
	error("Mismatch. nlayer=%d, nps=%d\n", dest->nlayer, nps);
    }
    dest->nlayer=nps;
    int nx0=source[0]->nx;
    int ny0=source[0]->ny;
    if(!dest->vx){
	if(type==1){
	    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
	    DO(cudaMalloc3DArray(&dest->ca, &channelDesc, make_cudaExtent(nx0, ny0, nps), 
				 cudaArrayLayered));
	}else{
	    DO(cudaMallocHost(&(dest->p), nps*sizeof(float*)));
	    //memory in device.
	    for(int ips=0; ips<nps; ips++){
		DO(cudaMalloc(&(dest->p[ips]), nx0*ny0*sizeof(float)));
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
    }
    
    float *tmp=NULL;
    if(type==1) tmp=new float[nx0*ny0*nps];
    
    for(int ips=0; ips<nps; ips++){
	int nx=source[ips]->nx;
	int ny=source[ips]->ny;
	if(type==1 && (nx!=nx0 || ny!=ny0)){
	    error("Only support map_t arrays of the same size if type==1\n");
	}
	dest->vx[ips]=source[ips]->vx;
	dest->vy[ips]=source[ips]->vy;
	dest->ht[ips]=source[ips]->h;
	dest->ox[ips]=source[ips]->ox;
	dest->oy[ips]=source[ips]->oy;
	dest->dx[ips]=source[ips]->dx;
	dest->nx[ips]=source[ips]->nx;
	dest->ny[ips]=source[ips]->ny;
	if(type==1){//cudaArray
	    for(long ix=0; ix<(long)nx0*(long)ny0; ix++){
		tmp[ips*nx0*ny0+ix]=(float)source[ips]->p[ix];
	    }
	}else{//Flat memory
	    gpu_dbl2dev(&dest->p[ips], source[ips]->p, nx*ny);
	    /*
	      float *tmp2=new float[nx*ny];
	      for(int iloc=0; iloc<nx*ny; iloc++){
	      tmp2[iloc]=(float)source[ips]->p[iloc];
	      }
	      cudaMemcpy(dest->p[ips], tmp2, nx*ny*sizeof(float), cudaMemcpyHostToDevice);
	      free(tmp2);*/
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
}

/*
  Convert a host dsp array to GPU sprase array. Both are in CSC format. 
*/
void sp2gpu(cusp_t *dest, dsp *src){
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
__global__ void cuspmul_do(float *y, cusp_t *A, float *x, float alpha){
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
void cuspmul(float *y, cusp_t *A, float *x, float alpha, cudaStream_t stream){
    cuspmul_do<<<A->nx/256, 256, 0, stream>>>(y,A,x,alpha);
}

/*
  y=A'*x where A is sparse. x, y are vectors
*/
__global__ void cusptmul_do(float *y, int icol, cusp_t *A, float *x, float alpha){
    __shared__ float val;
    if(threadIdx.x==0) val=0;
    int i=blockIdx.x * blockDim.x + threadIdx.x;
    int j=i+A->p[icol];
    atomicAdd(&val, A->x[j]*x[A->i[j]]);
    if(threadIdx.x==0) y[icol]+=val*alpha;
}
/*
  Does not work yet. Try to launch a block for each column and n items in each block.
*/
void cusptmul(float *y, cusp_t *A, float *x, float alpha, cudaStream_t stream){
    for(int i=0; i<A->ny; i++){
	cusptmul_do<<<1, A->p[i+1]-A->p[i], 0, stream>>>(y,i,A,x,alpha);
    }
}
__global__ static void calc_ptt_do( float *cc,
				    const float (*restrict loc)[2], 
				    const int nloc,
				    const float *restrict phi,
				    const float *restrict amp){
    __shared__ float ccb[4];//for each block.
    if(threadIdx.x==0){
	ccb[0]=ccb[1]=ccb[2]=ccb[3]=0.f;
    }
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
    if(threadIdx.x==0){//This is the first thread of a block. add block result to global.
	atomicAdd(&cc[0], ccb[0]);
	atomicAdd(&cc[1], ccb[1]);
	atomicAdd(&cc[2], ccb[2]);
	atomicAdd(&cc[3], ccb[3]);
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
    if(threadIdx.x==0){
	ccb[0]=ccb[1]=ccb[2]=ccb[3]=ccb[4]=ccb[5]=ccb[6]=0;
    }
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
    if(threadIdx.x==0){//This is the first thread of a block. add block result to global.
	atomicAdd(&cc[0], ccb[0]);
	atomicAdd(&cc[1], ccb[1]);
	atomicAdd(&cc[2], ccb[2]);
	atomicAdd(&cc[3], ccb[3]);
	atomicAdd(&cc[4], ccb[4]);
	atomicAdd(&cc[5], ccb[5]);
	atomicAdd(&cc[6], ccb[6]);
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
    cudaCallocHost(cc, 4*sizeof(float));
    calc_ptt_do<<<16, 256, 0, stream>>>(cc, loc, nloc, phi, amp);
    CUDA_SYNC_STREAM;
    double coeff[3], tot;
    coeff[0]=cc[0]; coeff[1]=cc[1]; coeff[2]=cc[2]; tot=cc[3];
    cudaFree(cc);
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
    cudaCallocHost(cc, 7*sizeof(float));
    cc[0]=cc[1]=cc[2]=cc[3]=cc[4]=cc[5]=cc[6]=0;
    if(nmod==2){//single DM.
	calc_ptt_do<<<16,256,0,stream>>>(cc, loc, nloc, phi, amp);
    }else if(nmod==5){//AHST mode
	calc_ngsmod_do<<<16,256,0,stream>>>(cc, loc, nloc, phi, amp);
    }else{
	TO_IMPLEMENT;
    }
    CUDA_SYNC_STREAM;
    tot=cc[nmod==2?3:6];

    double coeff[6];
    coeff[0]=cc[0]; coeff[1]=cc[1]; 
    coeff[2]=cc[2]; coeff[3]=cc[3];
    coeff[4]=cc[4]; coeff[5]=cc[5];
    
    cudaFree(cc); 
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
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/
void gpu_dbl2dev(float * restrict *dest, double *src, int n){
    float *tmp=(float*)malloc(n*sizeof(float));
    for(int i=0; i<n; i++){
	tmp[i]=(float)src[i];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, n*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(float),cudaMemcpyDefault));
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/
void gpu_dbl2flt(float * restrict *dest, double *src, int n){
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
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(long)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyDefault));
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((long)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyDefault));
	free(tmp);
    }
}
/**
   Convert long array to device int
*/
void gpu_spint2dev(int * restrict *dest, spint *src, int n){
    info("sizeof(spint)=%ld\n", sizeof(spint));
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(spint)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyDefault));
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((spint)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyDefault));
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
void gpu_dev2dbl(double * restrict *dest, float *src, int n){
    TIC;tic;
    float *tmp=(float*)malloc(n*sizeof(float));
    DO(cudaMemcpy(tmp, src, n*sizeof(float), cudaMemcpyDefault));
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
    cudaMemcpy(tmp, p, nx*ny*sizeof(float), cudaMemcpyDefault);
    writeflt(tmp,nx,ny,"%s",fn);
    free(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_writeint(int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDefault);
    writeint(tmp,nx,ny,"%s",fn);
    free(tmp);
}
