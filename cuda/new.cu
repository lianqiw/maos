#include <cuda.h>
#include <stdlib.h>
const int NG1D=64;
const int NG2D=8;
const int WRAP_SIZE=32; /*The wrap size is currently always 32 */
const int REDUCE_WRAP=4;
const int REDUCE_WRAP_LOG2=2;
const int DIM_REDUCE=WRAP_SIZE*REDUCE_WRAP; /*dimension to use in reduction. */
const int REDUCE_STRIDE=WRAP_SIZE+WRAP_SIZE/2+1;
extern "C" __global__ void sum_do(float * res, const float *a, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(res, sb[0]);
    }
}
/*
  In each block, we first do the reduction in each warp. This avoid syncthreads and if test. Then we copy results from each wrap to the first wrap and do the reduction again.
*/
extern "C" __global__ void sum2_do(float * res, const float *a, const int n){
    __shared__ float sb[REDUCE_WRAP*REDUCE_STRIDE];
    const int idx=threadIdx.x;
    const int wrap=idx/WRAP_SIZE; //which wrap
    const int jdx=(WRAP_SIZE-1) & idx;//index within this wrap
    volatile float *s=sb+REDUCE_STRIDE*wrap+jdx+WRAP_SIZE/2;
    s[-16]=0;
    //Read in vector from global mem
    register float sum=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + idx; i<n; i+=step){
	sum+=a[i];
	}
    s[0]=sum;
    //Handle each wrap without sync
#pragma unroll
    for(int i=0; i<5; i++){
	int offset=1<<i;
	sum += s[-offset];//every thread retrives current value
	s[0] = sum;//every thread write new value.
    }
    __syncthreads();//synchronize different wraps*/
    if(idx<REDUCE_WRAP){//use a few threads for reduce
	float sum2=sb[REDUCE_STRIDE * idx + WRAP_SIZE/2 + WRAP_SIZE - 1];
	//reuse sb for size of REDUCE_WRAP+REDUCE_WRAP/2;
	sb[idx]=0;
	volatile float *s2 = sb + REDUCE_WRAP/2 + idx;
	s2[0]=sum2;
#pragma unroll	
	for(int i=0; i<REDUCE_WRAP_LOG2; i++){
	    int offset=1<<i;
	    sum2+=s2[-offset];
	    s2[0]=sum2;
	}
	if(idx+1==REDUCE_WRAP){
	    atomicAdd(res, sum2);
	}
    }
}
