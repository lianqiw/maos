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

/*
  Ray tracing with matched spacing. Reverse, from out to in. out is xloc, in is ploc.
*/
__global__ static void prop_grid_match_do(float *restrict out, int nxout,
					  const float *restrict in, int nxin, 
					  float fracx, float fracy,
					  float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    float fracx1=1.f-fracx;
    float fracy1=1.f-fracy;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    out[ix+iy*nxout]+=
		alpha*(+(in[ix+    iy*nxin]*fracx1+in[ix+1+    iy*nxin]*fracx)*fracy1
		       +(in[ix+(iy+1)*nxin]*fracx1+in[ix+1+(iy+1)*nxin]*fracx)*fracy);
	}
    }
}

__global__ static void prop_grid_nomatch_do(float *restrict out, int nxo, const float *restrict in, 
					    int nxi, float dispx, float dispy, float ratio, float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	float jy;
	float fracy=modff(dispy+iy*ratio, &jy);
	int ky=(int)jy;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    float jx;
	    float fracx=modff(dispx+ix*ratio, &jx);
	    int kx=(int)jx;
	    out[ix+iy*nxo]+=
		alpha*(+(in[kx+      ky*nxi]*(1.f-fracx)+
			 in[kx+1+    ky*nxi]*fracx)*(1.f-fracy)
		       +(in[kx  +(ky+1)*nxi]*(1.f-fracx)+
			 in[kx+1+(ky+1)*nxi]*fracx)*fracy);
	}
    }
}
__global__ static void prop_grid_nomatch_trans_do(const float *restrict out, int nxo, float *restrict in, 
						  int nxi, float dispx, float dispy, float ratio, float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	float jy;
	float fracy=modff(dispy+iy*ratio, &jy);
	int ky=(int)jy;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    float jx;
	    float fracx=modff(dispx+ix*ratio, &jx);
	    int kx=(int)jx;
	    float temp=out[ix+iy*nxo]*alpha;
	    atomicAdd(&in[kx+      ky*nxi], temp*(1.f-fracx)*(1.f-fracy));
	    atomicAdd(&in[kx+1    +ky*nxi], temp*fracx*(1.f-fracy));
	    atomicAdd(&in[kx+  (ky+1)*nxi], temp*(1.f-fracx)*fracy);
	    atomicAdd(&in[kx+1+(ky+1)*nxi], temp*fracx*fracy);
	}
    }
}

/*
  Ray tracing with over sampling. Reverse, from out to in. out is xloc, in is
ploc. confirmed to agree with HXW'.  */
__global__ static void prop_grid_os2_trans_do(const float *restrict out, int nxout,
					      float *restrict in, int nxin, 
					      float fracx, float fracy,
					      float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int ax=fracx<0.5f?0:1;
    int ay=fracy<0.5f?0:1;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
	    //odd and even points are different.
#pragma unroll 
	    for(int by=0; by<2; by++){
		int iy2=iy*2+by;
		int iy3=iy+ay*by;
		float fracy2=fracy+(0.5f-ay)*by;
		float fracy21=1.f-fracy2;
		for(int bx=0; bx<2; bx++){
		    int ix2=ix*2+bx;
		    int ix3=ix+ax*bx;
		    float fracx2=fracx+(0.5f-ax)*bx;
		    float fracx21=1.f-fracx2;
		    if(ix2<nx && iy2<ny){
			float a=out[ix2+(iy2)*nxout]*alpha;
			atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
			atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
			atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
			atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
		    }
		}
	    }
	}
    }
}
/*
  Ray tracing with over sampling. Forward, from out to in. out is xloc, in is
ploc. confirmed to agree with HXW'.  */
__global__ static void prop_grid_os2_do(float *restrict out, int nxout,
					const float *restrict in, int nxin, 
					float fracx, float fracy,
					float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int ax=fracx<0.5f?0:1;
    int ay=fracy<0.5f?0:1;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
	    //odd and even points are different.
#pragma unroll 
	    for(int by=0; by<2; by++){
		int iy2=iy*2+by;
		int iy3=iy+ay*by;
		float fracy2=fracy+(0.5f-ay)*by;
		float fracy21=1.f-fracy2;
		for(int bx=0; bx<2; bx++){
		    int ix2=ix*2+bx;
		    int ix3=ix+ax*bx;
		    float fracx2=fracx+(0.5f-ax)*bx;
		    float fracx21=1.f-fracx2;
		    if(ix2<nx && iy2<ny){
			out[ix2+(iy2)*nxout]+=
			    alpha*(+in[ix3+    (iy3)*nxin]*fracx21*fracy21
				   +in[ix3+1+  (iy3)*nxin]*fracx2*fracy21
				   +in[ix3+  (iy3+1)*nxin]*fracx21*fracy2
				   +in[ix3+1+(iy3+1)*nxin]*fracx2*fracy2);
		    }
		}
	    }
	}
    }
}
/* Do the ray tracing when in/out grid matches in sampling. Handle offsets correctly.*/
/*
void prop_grid_match(curmat *out, float oxo, float oyo,
		     const curmat *in, float oxi, float oyi, float dxi,
		     float dispx, float dispy,
		     float alpha, cudaStream_t stream){
    const float dx1=1.f/dxi;
    dispx=(dispx-oxi+oxo)*dx1;
    dispy=(dispy-oyi+oyo)*dx1;
    const int offx=(int)floorf(dispx); dispx-=offx;
    const int offy=(int)floorf(dispy); dispy-=offy;
    const int nxo=out->nx;
    const int nyo=out->ny;
    const int nxi=in->nx;
    const int nyi=in->ny;
    int offx1, offx2;
    int offy1, offy2;
    int nx, ny;
    if(offx>0){
	offx1=0;
	offx2=offx;
	nx=nxi-offx-1; if(nx>nxo) nx=nxo;
    }else{
	offx1=-offx;
	offx2=0;
	nx=nxo+offx; if(nx>nxi-1) nx=nxi-1;
    }
    if(offy>0){
	offy1=0;
	offy2=offy;
	ny=nyi-offy-1; if(ny>nyo) ny=nyo;
    }else{
	offy1=-offy;
	offy2=0;
	ny=nyo+offy; if(ny>nyi-1) ny=nyi-1;
    }
    prop_grid_match_do<<<DIM2(nx,ny,16,8),0,stream>>>
	(out->p+offx1+offy1*nxo, nxo, in->p+offx2+offy2*nxi, nxi, dispx, dispy, alpha, nx, ny);
}
*/

/**
   Do the ray tracing
   from in to out if trans=='n'
   from out to in if trans=='t'
 */
void gpu_prop_grid(curmat *out, float oxo, float oyo, float dxo,
		   curmat *in, float oxi, float oyi, float dxi,
		   float dispx, float dispy,
		   float alpha, char trans, cudaStream_t stream){
    const float dxi1=1.f/dxi;
    const float ratio=dxo*dxi1;
    if(fabs(ratio-1.f)<1.e-4 && trans=='t'){
	gpu_prop_grid(in, oxi, oyi, dxi, out, oxo, oyo, dxo, -dispx, -dispy, alpha,'n', stream);
	return;
    }
    const int nxo=out->nx;
    const int nyo=out->ny;
    const int nxi=in->nx;
    const int nyi=in->ny;
    const float ratio1=1.f/ratio;
    //offset of origin in input grid spacing.
    dispx=(dispx-oxi+oxo)*dxi1;
    dispy=(dispy-oyi+oyo)*dxi1;
    int offx1=0, offy1=0;
    //if output is bigger than input.
    if(dispx<0){
	offx1=(int)ceilf(-dispx*ratio1);
	dispx+=offx1*ratio;
    }
    if(dispy<0){
	offy1=(int)ceilf(-dispy*ratio1);
	dispy+=offy1*ratio;
    }
    //convert offset into input grid coordinate. -1e-4 to avoid laying on the last point.
    int nx=(int)floorf((nxi-1-dispx-1e-4)*ratio1)+1;
    int ny=(int)floorf((nyi-1-dispy-1e-4)*ratio1)+1;

    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
    int offx2=(int)floorf(dispx); dispx-=offx2;
    int offy2=(int)floorf(dispy); dispy-=offy2;
    if(trans=='n'){
	if(fabs(ratio-1.f)<1.e-4){
	    prop_grid_match_do<<<DIM2(nx, ny, 8), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx, ny);
	}else if(fabs(ratio-0.5f)<1.e-4){
	    prop_grid_os2_do<<<DIM2(nx, ny, 8), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx, ny);
	}else{
	    prop_grid_nomatch_do<<<DIM2(nx, ny, 8), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, ratio, alpha, nx, ny);
	}
    }else if(trans=='t'){
	if(fabs(ratio-1.f)<1.e-4){
	    error("Please revert the input/output and call with trans='n'\n");
	}else if(fabs(ratio-0.5f)<1.e-4){
	    prop_grid_os2_trans_do<<<DIM2(nx, ny, 8), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx, ny);
	}else{
	    prop_grid_nomatch_trans_do<<<DIM2(nx, ny, 8), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, ratio, alpha, nx, ny);
	}
    }else{
	error("Invalid trans=%c\n", trans);
    }
}
