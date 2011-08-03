extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
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
   
   Lesson: 

   1) Do not declare file scope variables of __device__. Do not work.  

   2) When call cudaBindTextureToArray in multi-threaded algorithms, the binding
   will conflict between different threads.

   3) For layered 2D texture., must use cudaMalloc3DArray with cudaArrayLayered
*/

#define ATM_TEXTURE 0 //Use texture for ATM. Same speed as not after make p in device memory.
#define DM_TEXTURE 0  //Use texture for DM. not critical.
#define TO_IMPLEMENT error("Please implement")
#if ATM_TEXTURE
texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefatm;
#endif
#if DM_TEXTURE
texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefdm;
#endif
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
#define cudaCallocHost(P,N) ({cudaMallocHost(&P,N); cudaMemset(P,0,N);})
#define cudaCalloc(P,N) ({DO(cudaMalloc(&P,N));DO(cudaMemset(P,0,N));})

static cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
typedef struct{
    float (*loc)[2];//in device.
    float **phi;//in device.
    float dx;
    int nloc;
}culoc_t;
static int cunpowfs;
static culoc_t *cusaloc=NULL;//on host
static float (*cuplocs)[2]=NULL;
static float *cupamp=NULL;
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

static cumap_t cuatm={0};//array of cumap_t;
static cumap_t cudm={0};
static float *cc=NULL;

#define TIMING 0
#if TIMING == 1
static int nstream=0;
#define STREAM_NEW(stream) ({DO(cudaStreamCreate(&stream));info2("nstream=%d\n",lockadd(&nstream,1)+1);})
#define STREAM_DONE(stream) ({DO(cudaStreamDestroy(stream));info2("nstream=%d\n",lockadd(&nstream,-1)-1);})
#else
#define STREAM_NEW(stream) DO(cudaStreamCreate(&stream))
#define STREAM_DONE(stream) DO(cudaStreamDestroy(stream))
#endif
/*
  Question: when I declare CC as a float* here and initialize once, I get
  segmentation error due to cc appear as NULL using cuda-gdb.  

  Serious problem found: At time step 1, ray tracing through atm works well with
  stored memory. But at step 2, the ray tracing gives zero value. Implies memory
  has been cleared. Same problem with cc. it is cleared to all zeros.

*/
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
static void map2gpu(map_t **source, int nps, cumap_t *dest, int type){
    if(nps==0) return;
    if(dest->nlayer!=0 && dest->nlayer!=nps){
	error("Mismatch. nlayer=%d, nps=%d\n", dest->nlayer, nps);
    }
    dest->nlayer=nps;
    int nx0=source[0]->nx;
    int ny0=source[0]->ny;
    if(!dest->vx){
	if(type==1){
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
	if(type==1){
	    for(long ix=0; ix<(long)nx0*(long)ny0; ix++){
		tmp[ips*nx0*ny0+ix]=(float)source[ips]->p[ix];
	    }
	}else{
	    float *tmp2=new float[nx*ny];
	    for(int iloc=0; iloc<nx*ny; iloc++){
		tmp2[iloc]=(float)source[ips]->p[iloc];
	    }
	    cudaMemcpy(dest->p[ips], tmp2, nx*ny*sizeof(float), cudaMemcpyHostToDevice);
	    free(tmp2);
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
  Transfer atmospheric data to GPU.
*/
void gpu_atm2gpu(map_t **atm, int nps){
    TIC;tic
    cudaThreadSynchronize();
#if ATM_TEXTURE
    map2gpu(atm, nps, &cuatm, 1);
    texRefatm.addressMode[0] = cudaAddressModeWrap;
    texRefatm.addressMode[1] = cudaAddressModeWrap;
    texRefatm.filterMode     = cudaFilterModeLinear;
    texRefatm.normalized     = true; 
    DO(cudaBindTextureToArray(texRefatm, cuatm.ca, channelDesc));
#else
    map2gpu(atm, nps, &cuatm, 2);
#endif
    toc2("atm to gpu");//0.4 second.
}
/*
  Copy DM commands to GPU.
*/
void gpu_dm2gpu(map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    DO(cudaThreadSynchronize());
#if DM_TEXTURE
    map2gpu(dmreal, ndm, &cudm, 1);
    texRefdm.addressMode[0] = cudaAddressModeClamp;
    texRefdm.addressMode[1] = cudaAddressModeClamp;
    texRefdm.filterMode     = cudaFilterModePoint;
    texRefdm.normalized     = false; 
    DO(cudaBindTextureToArray(texRefdm, cudm.ca, channelDesc));
#else
    map2gpu(dmreal, ndm, &cudm, 2);
#endif
    if(dmcfg && !cudm.cubic){
	cudm.cubic=new int[ndm];
	cudm.iac=new float[ndm];
	for(int idm=0; idm<ndm; idm++){
	    cudm.cubic[idm]=dmcfg[idm].cubic;
	    cudm.iac[idm]=dmcfg[idm].iac;
	}
    }
}
/*
  Save POWFS loc to GPU.
*/
void gpu_saloc2gpu(int npowfs, /**<Total number of powfs*/
		  int ipowfs, /**<Index of this ipowfs*/
		  int nwfs,   /**<Total number of wfs for this powfs*/
		  loc_t *loc  /**<The loc grid*/
		  ){
    DO(cudaThreadSynchronize());
    if(!cusaloc){
	cusaloc=(culoc_t*)calloc(npowfs, sizeof(culoc_t));
	cunpowfs=npowfs;
    }else if(cunpowfs!=npowfs){
	error("npowfs mismatch\n");
    }
    if(cusaloc[ipowfs].loc){
	error("Already initialized\n");
    }
    DO(cudaMalloc(&cusaloc[ipowfs].loc, loc->nloc*sizeof(float)*2));
    cusaloc[ipowfs].phi=new float*[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	DO(cudaMalloc(&(cusaloc[ipowfs].phi[iwfs]), loc->nloc*sizeof(float)));//on device.
    }
    cusaloc[ipowfs].nloc=loc->nloc;
    cusaloc[ipowfs].dx=loc->dx;
    float (*loctmp)[2]=new float[loc->nloc][2];
    for(int i=0; i<loc->nloc; i++){
	loctmp[i][0]=loc->locx[i];
	loctmp[i][1]=loc->locy[i];
    }
    DO(cudaMemcpy(cusaloc[ipowfs].loc, loctmp, loc->nloc*sizeof(float)*2, cudaMemcpyDefault));
    delete [] loctmp;
}

/*
  save aper_locs, aper_amp to GPU.
*/
void gpu_plocs2gpu(loc_t *loc, dmat *amp){
    if(cuplocs) error("Already Copied.\n");
    int nloc=loc->nloc;
    float (*tmp)[2]=(float(*)[2])malloc(sizeof(float)*2*nloc);
    float *amp0=(float*)malloc(sizeof(float)*nloc);
    for(int iloc=0; iloc<nloc; iloc++){
	tmp[iloc][0]=(float)loc->locx[iloc];
	tmp[iloc][1]=(float)loc->locy[iloc];
	amp0[iloc]=amp->p[iloc];
    }
    DO(cudaMalloc(&cuplocs, sizeof(float)*2*nloc));
    DO(cudaMemcpy(cuplocs, tmp, sizeof(float)*2*nloc, cudaMemcpyHostToDevice));
    DO(cudaMalloc(&cupamp, sizeof(float)*nloc));
    DO(cudaMemcpy(cupamp, amp0, sizeof(float)*nloc, cudaMemcpyHostToDevice));
    free(tmp);
    free(amp0);
}

#define KARG_COMMON float (*restrict loc)[2], int nloc, float dx, float dy, float dispx, float dispy, float alpha
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
	if(x>=0 && x<nx-1 && y>=0 && y<ny-1){
	    out[i]+=(in[iy*nx+ix]*(1-x)+in[iy*nx+ix+1]*x)*(1-y)
		+(in[(iy+1)*nx+ix]*(1-x)+in[(iy+1)*nx+ix+1]*x)*y;
	}
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
static void atm2loc(float *phiout, float (*loc)[2], int nloc, float hs, float thetax, float thetay,
		    float mispx, float mispy, float dtisim, float atmalpha, cudaStream_t stream){

    if(fabs(atmalpha)<EPS) return;
    for(int ips=0; ips<cuatm.nlayer; ips++){
	float dx=cuatm.dx[ips];
	float du=1.f/dx;
#if ATM_TEXTURE
	float scx=du/(float)cuatm.nx[ips];
	float scy=du/(float)cuatm.ny[ips];
#define offset 0.5	    
#else
#define scx du
#define scy du
#define offset 0
#endif
	float ht=cuatm.ht[ips];
	float vx=cuatm.vx[ips];
	float vy=cuatm.vy[ips];
	float dispx=(ht*thetax+mispx-vx*dtisim-cuatm.ox[ips]+offset*dx)*scx;
	float dispy=(ht*thetay+mispy-vy*dtisim-cuatm.oy[ips]+offset*dx)*scy;
	float scale=1.f-ht/hs;

#define COMM loc,nloc,scale*scx,scale*scy, dispx, dispy, atmalpha
#if ATM_TEXTURE
#define FUN prop_atm
#define KARG ips, COMM
#else
#define FUN prop_linear
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
static void dm2loc(float *phiout, float (*loc)[2], int nloc, float hs, float thetax, float thetay,
		   float mispx, float mispy, float dmalpha, cudaStream_t stream){
    if(fabs(dmalpha)<EPS && !cudm.nx) return;
    for(int idm=0; idm<cudm.nlayer; idm++){
	int cubic=cudm.cubic[idm];
	float iac=cudm.iac[idm];
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
	float dx=cudm.dx[idm];
	float du=1.f/dx;
	//Bind automatically unbinds previous.
	    
	//In cubic, cudaFilterModePoint does not need 0.5 offset.
	float ht=cudm.ht[idm];
	float dispx=(ht*thetax+mispx-cudm.ox[idm])*du;
	float dispy=(ht*thetay+mispy-cudm.oy[idm])*du;
	float scale=1.f-ht/hs;

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
__global__ static void calc_ptt_do( float *cc,
				    const float (*restrict loc)[2], 
				    const int nloc,
				    const float *restrict phi,
				    const float *restrict amp){
    int step=blockDim.x * gridDim.x; 
    float cci[4]={0.f,0.f,0.f,0.f};//for each thread
    __shared__ float ccb[4];//for each block.
    if(threadIdx.x==0){
	ccb[0]=ccb[1]=ccb[2]=ccb[3]=0.f;
    }
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
static void gpu_calc_ptt(double *rmsout, double *coeffout, 
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
static void gpu_calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
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
  Ray tracing for WFS. \todo Expand to do gradients in GPU without transfering
  data back to CPU.
 */
void gpu_wfs(gpu_wfs_t *info){
    //TIC;tic;
    const PARMS_T *parms=info->parms;
    const POWFS_T *powfs=info->powfs;
    int iwfs=info->iwfs;
    int isim=info->isim;
    int ipowfs=parms->wfs[iwfs].powfs;
    int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    float hs=parms->powfs[ipowfs].hs;
    float thetax=parms->wfs[iwfs].thetax;
    float thetay=parms->wfs[iwfs].thetay;
    float mispx=powfs[ipowfs].misreg[wfsind][0];
    float mispy=powfs[ipowfs].misreg[wfsind][1];
    float dtisim=parms->sim.dt*isim;
    float atmalpha=info->atmalpha;
    float dmalpha=info->dmalpha;
    float (*loc)[2]=cusaloc[ipowfs].loc;
    int nloc=cusaloc[ipowfs].nloc;
    float *phiout=cusaloc[ipowfs].phi[wfsind];
    DO(cudaMemset(phiout, 0, nloc*sizeof(float)));
    cudaStream_t stream;
    STREAM_NEW(stream);
    atm2loc(phiout, loc, nloc, hs, thetax, thetay, mispx, mispy, dtisim, atmalpha, stream);
    dm2loc(phiout, loc, nloc, hs, thetax, thetay, mispx, mispy, dmalpha, stream);
    CUDA_SYNC_STREAM;
    //DO(cudaStreamSynchronize(stream));
    //toc2("gpu part");
    float *phiout2=(float*)malloc(nloc*sizeof(float));
    DO(cudaMemcpy(phiout2, phiout, nloc*sizeof(float), cudaMemcpyDefault));
    for(int i=0; i<nloc; i++){
	info->phi->p[i]+=phiout2[i];//this assignment is bottleneck.
    }
    free(phiout2);
    STREAM_DONE(stream);
    //toc2("final");
}
/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
 */
void gpu_perfevl(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int ievl=info->start;
    assert(info->end==info->start+1);//only one evl.
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int do_psf=(parms->evl.psfmean || parms->evl.psfhist) && isim>=parms->evl.psfisim;
    const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
    const int nloc=aper->locs->nloc;
    const double thetax=parms->evl.thetax[ievl];
    const double thetay=parms->evl.thetay[ievl];
    //Setup pointers for easy usage
    PDMAT(simu->olmp->p[ievl],polmp);//OL mode for each dir
    PDMAT(simu->olep->p[ievl],polep);//OL error for each dir
    PDMAT(simu->clmp->p[ievl],pclmp);
    PDMAT(simu->clep->p[ievl],pclep);

    float *iopdevl; 
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    cudaCalloc(iopdevl, aper->locs->nloc*sizeof(float));
    cudaStream_t stream;
    STREAM_NEW(stream);

    if(parms->sim.idealevl){
	error("Please finished by: \n"
	      "1) make aloc square, \n"
	      "2) send dmproj to this file by calling gpu_dm2gpu\n");
    }else if(simu->atm && !parms->sim.wfsalias){
	atm2loc(iopdevl, cuplocs, nloc, parms->evl.hs[ievl], thetax, thetay, 
		parms->evl.misreg[0], parms->evl.misreg[1], isim*dt, 1, stream);
    }
    CUDA_SYNC_STREAM;
    if(simu->telws){//Wind shake
	TO_IMPLEMENT;
    }
    if(simu->surfevl && simu->surfevl->p[ievl]){
	TO_IMPLEMENT;
    }
    if(save_evlopd){
	TO_IMPLEMENT;
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(nmod==3){
	gpu_calc_ptt(polep[isim], polmp[isim], aper->ipcc, aper->imcc,
		 cuplocs, nloc, iopdevl, cupamp, stream);
    }else{
	TO_IMPLEMENT;
    }
  
    if(parms->evl.psfmean &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
			     ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	TO_IMPLEMENT;
    }
    
    if(parms->sim.evlol) goto end;
    
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	dm2loc(iopdevl, cuplocs, nloc, parms->evl.hs[ievl], thetax, thetay,
	       parms->evl.misreg[0], parms->evl.misreg[1], -1, stream);
	CUDA_SYNC_STREAM;
	if(imoao>-1){
	    TO_IMPLEMENT;
	}
    }
    if(save_evlopd){
	TO_IMPLEMENT;
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(parms->tomo.split){
	if(parms->ndm<=2){
	    PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
	    if(nmod==3){
		gpu_calc_ngsmod(pclep[isim], pclmp[isim], pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl, cupamp, stream);
	    }else{
		gpu_calc_ngsmod(NULL, NULL, pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl, cupamp, stream);	
		TO_IMPLEMENT;
	    }
	}
    }else{
	if(nmod==3){
	    gpu_calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,
			 cuplocs, nloc, iopdevl, cupamp, stream);
	}else{
	    TO_IMPLEMENT;
	}
    }
    if(parms->evl.psf[ievl] && isim>=parms->evl.psfisim && do_psf){
	TO_IMPLEMENT;
    }
 end:
    STREAM_DONE(stream);
    cudaFree(iopdevl);
}
