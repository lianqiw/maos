extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#define AOS_CUDA_ACCPHI_CU
#undef EPS
#define EPS 1e-5
/**
   First written on 2011-07.
   Port wfsgrad.c perfevl, etc to GPU.

   Change log:

   1) memcpy failed for cc when declared with float cc[5] with invalid argument
   when using cudaMemcpyDeviceToHost. Changed to cudaMemcpyDefault works. But
   the consequency is that cc in kernel have only zero value. When use float*
   cc, the memcpy succeed, but cc is null in kernel. It seems we can only pass
   cc into the kernel.

   2) copying DM information to cuda messes up atm because gpu_dm2gpu used texRefatm.
*/

#define ATM_TEXTURE 1 //Use texture for ATM. Same speed as not after make p in device memory.
#define DM_TEXTURE 0  //Use texture for DM. not critical.

#if ATM_TEXTURE
texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefatm;
#endif
#if DM_TEXTURE
texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefdm;
#endif


static cumap_t cuatm={0};//array of cumap_t;
static cumap_t cudm={0};
static float *cc=NULL;



/*
  Question: when I declare CC as a float* here and initialize once, I get
  segmentation error due to cc appear as NULL using cuda-gdb.  

  Serious problem found: At time step 1, ray tracing through atm works well with
  stored memory. But at step 2, the ray tracing gives zero value. Implies memory
  has been cleared. Same problem with cc. it is cleared to all zeros.

*/

/*
  Transfer atmospheric data to GPU.
*/
void gpu_atm2gpu(map_t **atm, int nps){
    TIC;tic;
#if ATM_TEXTURE
    gpu_map2dev(&cuatm, atm, nps, 1);
    texRefatm.addressMode[0] = cudaAddressModeWrap;
    texRefatm.addressMode[1] = cudaAddressModeWrap;
    texRefatm.filterMode     = cudaFilterModeLinear;
    texRefatm.normalized     = true; 
    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
    DO(cudaBindTextureToArray(texRefatm, cuatm.ca, channelDesc));
#else
    gpu_map2dev(&cuatm, atm, nps, 2);
#endif
    toc2("atm to gpu");//0.4 second.
}
/*
  Copy DM commands to GPU.
*/
void gpu_dm2gpu(map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
#if DM_TEXTURE
    gpu_map2dev(&cudm, dmreal, ndm, 1);
    texRefdm.addressMode[0] = cudaAddressModeClamp;
    texRefdm.addressMode[1] = cudaAddressModeClamp;
    texRefdm.filterMode     = cudaFilterModePoint;
    texRefdm.normalized     = false; 
    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
    DO(cudaBindTextureToArray(texRefdm, cudm.ca, channelDesc));
#else
    gpu_map2dev(&cudm, dmreal, ndm, 2);
#endif
    if(dmcfg && !cudm.cubic){
	cudm.cubic=new int[ndm];
	cudm.iac=new float[ndm];
	for(int idm=0; idm<ndm; idm++){
	    cudm.cubic[idm]=dmcfg[idm].cubic;
	    cudm.iac[idm]=dmcfg[idm].iac;
	}
    }
    CUDA_SYNC_DEVICE;
}

#define KARG_COMMON const float (*restrict loc)[2], const int nloc, const float dx, const float dy, const float dispx, const float dispy, const float alpha
#if ATM_TEXTURE
/*
  Ray tracing from texture to atm.
*/
__global__ static void prop_atm(float *restrict out, const int ilayer, 
				KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	out[i]+=tex2DLayered(texRefatm, x, y, ilayer)*alpha;
    }
}
#endif
#if DM_TEXTURE 
/*
  Ray tracing from texture to dm.
*/
__global__ static void prop_dm_linear(float *restrict out,  const int ilayer, 
				KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x;i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	out[i]+=tex2DLayered(texRefdm, x, y, ilayer)*alpha;
    }
}
/*
  Ray tracing from texture to dm with cubic influence functions..
*/
__global__ static void prop_dm_cubic(float *restrict out, const int ilayer,
				     KARG_COMMON, const float *cc){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x); x=x-ix;
	int iy=floorf(y); y=y-iy;
	float fx[4],fy, sum=0;
   
	fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
	fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
	fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
	fx[3]=x*x*(cc[3]+cc[4]*x);		
	fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
	fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
	fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
	fx[3]=x*x*(cc[3]+cc[4]*x);		
		
	fy=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y)); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*tex2DLayered(texRefdm, kx+ix, iy-1, ilayer);
	}

	fy=cc[0]+y*y*(cc[1]+cc[2]*y); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*tex2DLayered(texRefdm, kx+ix, iy, ilayer);
	}

	fy=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y)); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*tex2DLayered(texRefdm, kx+ix, iy+1, ilayer);
	}

	fy=y*y*(cc[3]+cc[4]*y); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*tex2DLayered(texRefdm, kx+ix, iy+2, ilayer);
	}

	out[i]+=sum*alpha;
    }
}
#endif
#if !ATM_TEXTURE || !DM_TEXTURE 
//This is memory bound. So increasing # of points processed does not help.
__global__ static void prop_linear(float *restrict out, const float *restrict in, const int nx, const int ny,
				      KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x);
	int iy=floorf(y);
	x=x-ix; y=y-iy;
	if(ix>=0 && ix<nx-1 && iy>=0 && iy<ny-1){
	    out[i]+=(in[iy*nx+ix]*(1-x)+in[iy*nx+ix+1]*x)*(1-y)
		+(in[(iy+1)*nx+ix]*(1-x)+in[(iy+1)*nx+ix+1]*x)*y;
	}
    }
}
//This is memory bound. So increasing # of points processed does not help.
__global__ static void prop_linear_wrap(float *restrict out, const float *restrict in, const int nx, const int ny,
				      KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x);
	int iy=floorf(y);
	x=x-ix; y=y-iy;
	while(ix<0) ix=ix+nx; 
	while(iy<0) iy=iy+ny;
	while(ix>nx-1) ix=ix-nx; 
	while(iy>ny-1) iy=iy-ny;
	int ix1=(ix==nx-1)?0:ix+1;
	int iy1=(iy==ny-1)?0:iy+1;
	out[i]+=(in[iy*nx+ix]*(1-x)+in[iy*nx+ix1]*x)*(1-y)
	    +(in[(iy1)*nx+ix]*(1-x)+in[(iy1)*nx+ix1]*x)*y;
    }
}
#endif
#if !DM_TEXTURE 
//This is memory bound. So increasing # of points processed does not help.
__global__ static void prop_cubic(float *restrict out, const float *restrict in, const int nx, const int ny,
				     KARG_COMMON, const float *cc){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x); x=x-ix;
	int iy=floorf(y); y=y-iy;
	float fx[4],fy;
	float sum=0;
	if(ix<1 || ix>nx-3 || iy<1 || iy>ny-3){
	    return;//out of range.
	}

	fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
	fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
	fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
	fx[3]=x*x*(cc[3]+cc[4]*x);		

	fy=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y)); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*in[(iy-1)*nx+kx+ix];
	}

	fy=cc[0]+y*y*(cc[1]+cc[2]*y); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*in[iy*nx+kx+ix];
	}

	fy=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y)); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*in[(iy+1)*nx+kx+ix];
	}

	fy=y*y*(cc[3]+cc[4]*y); 
#pragma unroll
	for(int kx=-1; kx<3; kx++){
	    sum+=fx[kx+1]*fy*in[(iy+2)*nx+kx+ix];
	}
	out[i]+=sum*alpha;
    }
}
#endif

/**
   Ray tracing of atm.
 */
void gpu_atm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax,const float thetay,
		 const float mispx, const float mispy, const float dtisim, const float atmalpha, cudaStream_t stream){

    if(fabs(atmalpha)<EPS) return;
    for(int ips=0; ips<cuatm.nlayer; ips++){
	const float dx=cuatm.dx[ips];
	const float du=1.f/dx;
#if ATM_TEXTURE
	const float scx=du/(float)cuatm.nx[ips];
	const float scy=du/(float)cuatm.ny[ips];
#define offset 0.5	    
#else
#define scx du
#define scy du
#define offset 0
#endif
	const float ht=cuatm.ht[ips];
	const float vx=cuatm.vx[ips];
	const float vy=cuatm.vy[ips];
	const float dispx=(ht*thetax+mispx-vx*dtisim-cuatm.ox[ips]+offset*dx)*scx;
	const float dispy=(ht*thetay+mispy-vy*dtisim-cuatm.oy[ips]+offset*dx)*scy;
	const float scale=1.f-ht/hs;

#define COMM loc,nloc,scale*scx,scale*scy, dispx, dispy, atmalpha
#if ATM_TEXTURE
#define FUN prop_atm
#define KARG ips, COMM
#else
#define FUN prop_linear_wrap
#define KARG cuatm.p[ips], cuatm.nx[ips], cuatm.ny[ips], COMM
#endif
	FUN <<<nloc/256, 256, 0, stream>>> (phiout, KARG);
#undef KARG
#undef FUN
#undef COMM
    }    
}

/**
  Ray tracing of dm
*/
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax, const float thetay,
		   const float mispx, const float mispy, const float dmalpha, cudaStream_t stream){
    if(fabs(dmalpha)<EPS && !cudm.nx) return;
    for(int idm=0; idm<cudm.nlayer; idm++){
	const int cubic=cudm.cubic[idm];
	const float iac=cudm.iac[idm];
	if(cubic){
	    if(!cc) {
		DO(cudaMallocHost(&cc, 5*sizeof(float)));
	    }
	    float cubicn=1.f/(1.f+2.f*iac);
	    cc[0]=1.f*cubicn;
	    cc[1]=(4.f*iac-2.5f)*cubicn; 
	    cc[2]=(1.5f-3.f*iac)*cubicn;		       
	    cc[3]=(2.f*iac-0.5f)*cubicn;			
	    cc[4]=(0.5f-iac)*cubicn; 
	}
	const float dx=cudm.dx[idm];
	const float du=1.f/dx;
	//Bind automatically unbinds previous.
	    
	//In cubic, cudaFilterModePoint does not need 0.5 offset.
	const float ht=cudm.ht[idm];
	const float dispx=(ht*thetax+mispx-cudm.ox[idm])*du;
	const float dispy=(ht*thetay+mispy-cudm.oy[idm])*du;
	const float scale=1.f-ht/hs;

#define COMM loc,nloc,scale*du,scale*du, dispx, dispy, dmalpha
#if DM_TEXTURE
#define KARG idm, COMM
#define FUNC prop_dm_cubic
#define FUNL prop_dm_linear
#else
#define KARG cudm.p[idm],cudm.nx[idm],cudm.ny[idm], COMM
#define FUNC prop_cubic
#define FUNL prop_linear
#endif
	if (cubic){//128 is a good number for cubic.
	    FUNC <<<128, 256, 0, stream>>>(phiout, KARG, cc);
	}else{
	    FUNL <<<nloc/256, 256, 0, stream>>>(phiout, KARG);
	}
#undef KARG
#undef COMM
#undef FUNC
#undef FUNL
    }//idm
}
