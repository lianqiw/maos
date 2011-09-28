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
//Removed option of making DM in texture
#if ATM_TEXTURE
texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefatm;
#else
static int WRAP_ATM;
#endif

/*
  Question: when I declare CC as a float* here and initialize once, I get
  segmentation error due to cc appear as NULL using cuda-gdb.  

  Serious problem found: At time step 1, ray tracing through atm works well with
  stored memory. But at step 2, the ray tracing gives zero value. Implies memory
  has been cleared. Same problem with cc. it is cleared to all zeros.

*/
typedef struct{
    int ips;
    float *next_atm;
    int nx0;
    int ny0;
    int offx;
    int offy;
    map_t *atm;
}atm_prep_t;
void atm_prep(atm_prep_t *data){
    PNEW(lock);
    const int ips=data->ips;
    LOCK(lock);//make sure we only read one layer at a time.
    TIC;tic;
    const int nx0=data->nx0;
    const int ny0=data->ny0;
    const map_t *atm=data->atm;
    const int nxi=atm->nx;
    const int nyi=atm->ny;
    int offx=data->offx;
    int offy=data->offy;
    offx=offx%nxi; if(offx<0) offx+=nxi;
    offy=offy%nyi; if(offy<0) offy+=nyi;
    
    int mx, my;
    if(offx+nx0>nxi){
	mx=nxi-offx;
    }else{
	mx=nx0;
    }
    if(offy+ny0>nyi){
	my=nyi-offy;
    }else{
	my=ny0;
    }
    PDMAT(atm, pin);
    float (*pout)[nx0]=(float(*)[nx0])(data->next_atm);
    for(int iy=0; iy<my; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(float)pin[iy+offy][ix+offx];
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(float)pin[iy+offy][ix+offx-nxi];
	}
    }
    for(int iy=my; iy<ny0; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(float)pin[iy+offy-nyi][ix+offx];
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(float)pin[iy+offy-nyi][ix+offx-nxi];
	}
    }
    toc("Layer %d: Preparing atm", ips);
    UNLOCK(lock);
    free(data);//allocated in parent thread. we free it here.
}
/**
   Transfer atmosphere or update atmosphere in GPU.
*/
void gpu_atm2gpu_new(map_t **atm, const PARMS_T *parms, int iseed, int isim){
    if(parms->atm.evolve){
	TO_IMPLEMENT;
	//test whether screen changed. transfer if changed.
    }
    const int nps=parms->atm.nps;
    static int nx0=0,ny0=0;
    static int iseed0=-1;
    if(!nx0){
	if(parms->atm.nxm!=parms->atm.nx || parms->atm.nym!=parms->atm.ny){//user specified atmosphere size.
	    long avail=gpu_get_mem();
	    info("Available memory is %d\n", avail);
	    long nxa=(avail-600000000)/nps/sizeof(float);//we are able to host this amount.
	    if(nxa<0){
		error("GPU does not have enough memory\n");
	    }
	    info2("GPU can host %d %dx%d atmosphere\n", nps, (int)round(sqrt(nxa)), (int)round(sqrt(nxa)));
	    if(nxa>parms->atm.nx*parms->atm.ny){//we can host all atmosphere.
		info2("We can host all layers in full.\n");
		nx0=parms->atm.nx;
		ny0=parms->atm.ny;
	    }else{
		info2("We will host %dx%d in GPU\n", nx0, ny0);
		nx0=parms->atm.nxm;
		ny0=parms->atm.nym;
	    }
	}else{
	    nx0=ny0=parms->atm.nxm;
	}
    }
    //The atm in GPU is the same as in CPU.
    if(nx0==parms->atm.nx && ny0==parms->atm.ny){
#if !ATM_TEXTURE
	WRAP_ATM=1;
#endif
	if(iseed0!=iseed){
	    gpu_atm2gpu(atm, nps);
	    iseed0=iseed;
	}
	return;
    }
#if !ATM_TEXTURE
    WRAP_ATM=0;
#endif
    static int need_init=1;
    static int *next_isim=NULL;//next time step to update this atmosphere.
    static float *next_ox=NULL;
    static float *next_oy=NULL;
    static float **next_atm=NULL;
    static pthread_t *next_threads=NULL;
    if(need_init){
	need_init=0;
	//The atm in GPU is smaller than in CPU.
	next_isim=(int*)calloc(nps, sizeof(int));
	next_ox=(float*)calloc(nps, sizeof(float));
	next_oy=(float*)calloc(nps, sizeof(float));
	next_threads=(pthread_t*)calloc(nps, sizeof(pthread_t));
	next_atm=(float **)calloc(nps, sizeof(void*));

	for(int im=0; im<NGPU; im++){//Loop over all GPUs.
	    gpu_set(im);
	    cumap_t *cuatm=cudata->atm=(cumap_t*)calloc(1, sizeof(cumap_t));
	    cuatm->nlayer=nps;
#if ATM_TEXTURE
	    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
	    DO(cudaMalloc3DArray(&cuatm->ca, &channelDesc, make_cudaExtent(nx0, ny0, nps), 
				 cudaArrayLayered));
	    texRefatm.addressMode[0] = cudaAddressModeClamp;//Wrap; don't use wrap since we are cropped
	    texRefatm.addressMode[1] = cudaAddressModeClamp;//Wrap;
	    texRefatm.filterMode     = cudaFilterModeLinear;
	    texRefatm.normalized     = true; 
	    DO(cudaBindTextureToArray(texRefatm, cuatm->ca, channelDesc));
#else
	    DO(cudaMallocHost(&(cuatm->p), nps*sizeof(float*)));
	    for(int ips=0; ips<nps; ips++){
		DO(cudaMalloc(&(cuatm->p[ips]), nx0*ny0*sizeof(float)));
	    }
#endif
	    cuatm->vx=new float[nps];
	    cuatm->vy=new float[nps]; 
	    cuatm->ht=new float[nps];
	    cuatm->ox=new float[nps];
	    cuatm->oy=new float[nps];
	    cuatm->dx=new float[nps];
	    cuatm->nx=new int[nps];
	    cuatm->ny=new int[nps];
	    for(int ips=0; ips<nps; ips++){
		cuatm->nx[ips]=nx0;
		cuatm->ny[ips]=ny0;
	    }
	}//for im
    }//if need_init;
    const double dt=parms->sim.dt;
    const double dx=parms->atm.dx;
    if(iseed0!=iseed){//A new seed update vx, vy, ht, etc.
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    cumap_t *cuatm=cudata->atm;
	    for(int ips=0; ips<nps; ips++){
		cuatm->vx[ips]=atm[ips]->vx;
		cuatm->vy[ips]=atm[ips]->vy;
		cuatm->ht[ips]=atm[ips]->h;
		cuatm->dx[ips]=atm[ips]->dx;
		cuatm->ox[ips]=INFINITY;//place holder
		cuatm->oy[ips]=INFINITY;
		next_isim[ips]=isim;//right now.
		if(next_atm[ips]){
		    warning("next_atm[%d] is not empty!!!\n", ips);
		    cudaFreeHost(next_atm[ips]);
		    next_atm[ips]=NULL;
		}
		//copy from below.
		if(atm[ips]->vx>0){//align right
		    next_ox[ips]=(parms->atm.nxn/2+1-parms->atm.nxm)*dx-cuatm->vx[ips]*dt*next_isim[ips];
		}else{//align left
		    next_ox[ips]=(-parms->atm.nxn/2)*dx-cuatm->vx[ips]*dt*next_isim[ips];
		}
		if(atm[ips]->vy>0){//align right
		    next_oy[ips]=(parms->atm.nyn/2+1-parms->atm.nym)*dx-cuatm->vy[ips]*dt*next_isim[ips];
		}else{//align left
		    next_oy[ips]=(-parms->atm.nyn/2)*dx-cuatm->vy[ips]*dt*next_isim[ips];
		}
	    }
	}
    }
    for(int ips=0; ips<nps; ips++){
	/*Load atmosphere to float memory in advance. This takes time if atm is
	  stored in file.*/
	if(isim>next_isim[ips]-100 && !next_atm[ips]){
	    //pinned memory is faster for copying to GPU.
	    cudaMallocHost(&(next_atm[ips]), sizeof(float)*nx0*ny0);
	    const int nxi=atm[ips]->nx;
	    const int nyi=atm[ips]->ny;
	    int offx=(int)round((next_ox[ips]-atm[ips]->ox)/dx);
	    int offy=(int)round((next_oy[ips]-atm[ips]->oy)/dx);
	    next_ox[ips]=atm[ips]->ox+offx*dx;
	    next_oy[ips]=atm[ips]->oy+offy*dx;
	    offx=offx%nxi; if(offx<0) offx+=nxi;
	    offy=offy%nyi; if(offy<0) offy+=nyi;
	    atm_prep_t *data=(atm_prep_t*)calloc(1, sizeof(atm_prep_t));//cannot use local variable.
	    data->ips=ips;
	    data->next_atm=next_atm[ips];
	    data->nx0=nx0;
	    data->ny0=ny0;
	    data->offx=offx;
	    data->offy=offy;
	    data->atm=atm[ips];
	    //launch an independent thread to pull data in. thread will exit when it is done.
	    pthread_create(&next_threads[ips], NULL, (void *(*)(void *))atm_prep, data);
	}//need to cpy atm to next_atm.
	if(isim==next_isim[ips]){
	    //need to copy atm to gpu. and update next_isim
	    TIC;tic;
	    pthread_join(next_threads[ips], NULL);
	    tic("Layer %d wait ",ips);
	    for(int im=0; im<NGPU; im++){
		tic;
		gpu_set(im);
		cumap_t *cuatm=cudata->atm;
		cuatm->ox[ips]=next_ox[ips];
		cuatm->oy[ips]=next_oy[ips];
#if ATM_TEXTURE
		/*The extent field defines the dimensions of the transferred area in
		  elements. If a CUDA array is participating in the copy, the extent
		  is defined in terms of that array's elements. If no CUDA array is
		  participating in the copy then the extents are defined in elements
		  of unsigned char.*/
		struct cudaMemcpy3DParms par={0};
		par.srcPos = make_cudaPos(0,0,0);
		par.srcPtr = make_cudaPitchedPtr((float*)next_atm[ips], nx0*sizeof(float), nx0, ny0);
		par.dstPos = make_cudaPos(0,0,ips);
		par.dstArray=cuatm->ca;
		par.extent = make_cudaExtent(nx0, ny0, 1);
		par.kind   = cudaMemcpyHostToDevice;
		DO(cudaMemcpy3D(&par));
#else
		DO(cudaMemcpy(cuatm->p[ips], (float*)next_atm[ips], 
			      nx0*ny0*sizeof(float), cudaMemcpyHostToDevice));
#endif
		CUDA_SYNC_DEVICE;
		int offx=(int)round((next_ox[ips]-atm[ips]->ox)/dx);
		int offy=(int)round((next_oy[ips]-atm[ips]->oy)/dx);
		toc2("Step %d: Copying layer %d to GPU %d: offx=%d, offy=%d", 
		     isim, ips, GPUS[im], offx, offy);tic;
	    }//for im
	    cudaFreeHost(next_atm[ips]);
	    next_atm[ips]=NULL;
	    //Update next_isim.
	    long isim1, isim2;
	    if(atm[ips]->vx>0){//align right.
		isim1=(long)floor(-(next_ox[ips]+(parms->atm.nxn/2)*dx)/(atm[ips]->vx*dt));
	    }else{//align left
		isim1=(long)floor(-(next_ox[ips]+(nx0-(parms->atm.nxn/2+1))*dx)/(atm[ips]->vx*dt));
	    }
	    if(atm[ips]->vy>0){//align right.
		isim2=(long)floor(-(next_oy[ips]+(parms->atm.nyn/2)*dx)/(atm[ips]->vy*dt));
	    }else{//align left
		isim2=(long)floor(-(next_oy[ips]+(ny0-(parms->atm.nyn/2+1))*dx)/(atm[ips]->vy*dt));
	    }
	    next_isim[ips]=isim1<isim2?isim1:isim2;
	    if(next_isim[ips]>parms->sim.end){//no need to do
		next_isim[ips]=INT_MAX;
	    }
	    if(atm[ips]->vx>0){//align right
		next_ox[ips]=(parms->atm.nxn/2+1-parms->atm.nxm)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }else{//align left
		next_ox[ips]=(-parms->atm.nxn/2)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }
	    if(atm[ips]->vy>0){//align right
		next_oy[ips]=(parms->atm.nyn/2+1-parms->atm.nym)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }else{//align left
		next_oy[ips]=(-parms->atm.nyn/2)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }
	}
    }
    iseed0=iseed;
}
/**
   Transfer atmospheric data to GPU.
*/
void gpu_atm2gpu(map_t **atm, int nps){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	gpu_print_mem("atm in");
	TIC;tic;
#if ATM_TEXTURE
	gpu_map2dev(&cudata_all[im]->atm, atm, nps, 1);
	texRefatm.addressMode[0] = cudaAddressModeWrap;
	texRefatm.addressMode[1] = cudaAddressModeWrap;
	texRefatm.filterMode     = cudaFilterModeLinear;
	texRefatm.normalized     = true; 
	cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
	DO(cudaBindTextureToArray(texRefatm, cudata_all[im]->atm->ca, channelDesc));
#else
	gpu_map2dev(&cudata_all[im]->atm, atm, nps, 2);
#endif
	toc2("atm to gpu");//0.4 second.
	gpu_print_mem("atm out");
    }
}
/**
   Copy DM commands to GPU.
*/
static void gpu_dm2gpu(cumap_t **cudm, map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    gpu_map2dev(cudm, dmreal, ndm, 2);
    if(dmcfg && !(*cudm)->cubic){
	(*cudm)->cubic=new int[ndm];
	(*cudm)->iac=new float[ndm];
	(*cudm)->cc=new float*[ndm];
	for(int idm=0; idm<ndm; idm++){
	    (*cudm)->cubic[idm]=dmcfg[idm].cubic;
	    (*cudm)->iac[idm]=dmcfg[idm].iac;
	    DO(cudaMallocHost(&((*cudm)->cc[idm]), 5*sizeof(float)));
	    float *cc=(*cudm)->cc[idm];
	    const float iac=(*cudm)->iac[idm];
	    float cubicn=1.f/(1.f+2.f*iac);
	    cc[0]=1.f*cubicn;
	    cc[1]=(4.f*iac-2.5f)*cubicn; 
	    cc[2]=(1.5f-3.f*iac)*cubicn;		       
	    cc[3]=(2.f*iac-0.5f)*cubicn;			
	    cc[4]=(0.5f-iac)*cubicn; 
	}
    }
    CUDA_SYNC_DEVICE;
}
void gpu_dmreal2gpu(map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	gpu_dm2gpu(&cudata_all[im]->dmreal, dmreal, ndm, dmcfg);
    }
}
void gpu_dmproj2gpu(map_t **dmproj, int ndm, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	gpu_dm2gpu(&cudata_all[im]->dmproj, dmproj, ndm, dmcfg);
    }
}

#define KARG_COMMON const float (*restrict loc)[2], const int nloc, const float dx, const float dy, const float dispx, const float dispy, const float alpha
#if ATM_TEXTURE
/*
  Ray tracing from texture to atm.
*/
__global__ void prop_atm(float *restrict out, const int ilayer, 
			 KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	out[i]+=tex2DLayered(texRefatm, x, y, ilayer)*alpha;
    }
}
#endif

//This is memory bound. So increasing # of points processed does not help.
__global__ void prop_linear(float *restrict out, const float *restrict in, const int nx, const int ny,
				   KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x);
	int iy=floorf(y);
	x=x-ix; y=y-iy;
	if(ix>=0 && ix<nx-1 && iy>=0 && iy<ny-1){
	    out[i]+= alpha*((+in[iy*nx+ix]*(1.f-x) +in[iy*nx+ix+1]*x)*(1.f-y)
			    +(+in[(iy+1)*nx+ix]*(1.f-x) +in[(iy+1)*nx+ix+1]*x)*y);
	}
    }
}
//This is memory bound. So increasing # of points processed does not help.
__global__ void prop_linear_nocheck(float *restrict out, const float *restrict in, 
				    const int nx, const int ny, KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float x=loc[i][0]*dx+dispx;
	float y=loc[i][1]*dy+dispy;
	int ix=floorf(x);
	int iy=floorf(y);
	x=x-ix; y=y-iy;
	out[i]+=alpha*((in[iy*nx+ix]*(1-x)+in[iy*nx+ix+1]*x)*(1-y)
		       +(in[(iy+1)*nx+ix]*(1-x)+in[(iy+1)*nx+ix+1]*x)*y);
    }
}
//This is memory bound. So increasing # of points processed does not help.
__global__ void prop_linear_wrap(float *restrict out, const float *restrict in, const int nx, const int ny,
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
	out[i]+=alpha*((in[iy*nx+ix]*(1-x)+in[iy*nx+ix1]*x)*(1-y)
		       +(in[(iy1)*nx+ix]*(1-x)+in[(iy1)*nx+ix1]*x)*y);
    }
}

//This is memory bound. So increasing # of points processed does not help.
__global__ void prop_cubic(float *restrict out, const float *restrict in, const int nx, const int ny,
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

/**
   Ray tracing of atm.
*/
void gpu_atm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, 
		 const float thetax,const float thetay,
		 const float mispx, const float mispy, const float dtisim, const float atmalpha, cudaStream_t stream){
    cumap_t *cuatm=cudata->atm;
    if(fabs(atmalpha)<EPS) return;
    for(int ips=0; ips<cuatm->nlayer; ips++){
	const float dx=cuatm->dx[ips];
	const float du=1.f/dx;
#if ATM_TEXTURE
	const float scx=du/(float)cuatm->nx[ips];
	const float scy=du/(float)cuatm->ny[ips];
#define offset 0.5	    
#else
#define scx du
#define scy du
#define offset 0
#endif
	const float ht=cuatm->ht[ips];
	const float vx=cuatm->vx[ips];
	const float vy=cuatm->vy[ips];
	const float dispx=(ht*thetax+mispx-vx*dtisim-cuatm->ox[ips]+offset*dx)*scx;
	const float dispy=(ht*thetay+mispy-vy*dtisim-cuatm->oy[ips]+offset*dx)*scy;
	const float scale=1.f-ht/hs;

#define COMM loc,nloc,scale*scx,scale*scy, dispx, dispy, atmalpha
#if ATM_TEXTURE
	prop_atm<<<DIM(nloc,256), 0, stream>>> (phiout, ips, COMM);
#else
	if(WRAP_ATM){
	    prop_linear_wrap<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm->p[ips], cuatm->nx[ips], cuatm->ny[ips], COMM);
	}else{//we are gauranteed.
	    prop_linear_nocheck<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm->p[ips], cuatm->nx[ips], cuatm->ny[ips], COMM);
	}
#endif
#undef COMM
    }    
}

/**
   Ray tracing of dm
*/
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, cumap_t *cudm,
		const float hs, const float thetax, const float thetay,
		const float mispx, const float mispy, const float dmalpha, cudaStream_t stream){
    if(fabs(dmalpha)<EPS && !cudm->nx) return;
    for(int idm=0; idm<cudm->nlayer; idm++){
	const int cubic=cudm->cubic[idm];
	const float dx=cudm->dx[idm];
	const float du=1.f/dx;
	//Bind automatically unbinds previous.
	    
	//In cubic, cudaFilterModePoint does not need 0.5 offset.
	const float ht=cudm->ht[idm];
	const float dispx=(ht*thetax+mispx-cudm->ox[idm])*du;
	const float dispy=(ht*thetay+mispy-cudm->oy[idm])*du;
	const float scale=1.f-ht/hs;

#define COMM loc,nloc,scale*du,scale*du, dispx, dispy, dmalpha
#define KARG cudm->p[idm],cudm->nx[idm],cudm->ny[idm], COMM
	if (cubic){//128 is a good number for cubic.
	    prop_cubic<<<DIM(nloc,128), 0, stream>>>(phiout, KARG, cudm->cc[idm]);
	}else{
	    prop_linear<<<DIM(nloc,256), 0, stream>>>(phiout, KARG);
	}
#undef KARG
#undef COMM

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
/*
  __global__ static void prop_grid_os2_trans_old_do(const float *restrict out, int nxout,
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
  #pragma unroll 
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
  }*/
/**
   Do a single column (along x).
*/

__global__ static void prop_grid_os2_trans_col_do(float *restrict out, int nxout,
						  float *restrict in, int nxin, 
						  float fracx, float fracy2,
						  float alpha, int nx){
    int stepx=blockDim.x*gridDim.x;
    int ax=fracx<0.5f?0:1;
    const int iy3=0;
    float fracy21=1.f-fracy2;
    for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
#pragma unroll
	for(int bx=0; bx<2; bx++){
	    int ix2=ix*2+bx;
	    if(ix2<nx){
		int ix3=ix+ax*bx;
		float fracx2=fracx+(0.5f-ax)*bx;
		float fracx21=1.f-fracx2;
		float a = out[ix2]*alpha;
		atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
		atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
		atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
		atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
	    }
	}
    }
}
/**
   Do a single row (along y).
*/
__global__ static void prop_grid_os2_trans_row_do(float *restrict out, int nxout,
						  float *restrict in, int nxin, 
						  float fracx2, float fracy,
						  float alpha, int ny){
    int stepy=blockDim.y*gridDim.y;
    int ay=fracy<0.5f?0:1;
    const int ix3=0;
    float fracx21=1.f-fracx2;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
#pragma unroll
	for(int by=0; by<2; by++){
	    int iy2=iy*2+by;
	    if(iy2<ny){
		int iy3=iy+ay*by;
		float fracy2=fracy+(0.5f-ay)*by;
		float fracy21=1.f-fracy2;
		float a = out[iy2*nxout]*alpha;
		atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
		atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
		atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
		atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
	    }
	}
    }
}
/**
   Try to improve prop_grid_os2_trans_do using shared memory.
*/
__global__ static void prop_grid_os2_trans_share_do(float *restrict out, int nxout,
						    float *restrict in, int nxin, 
						    float fracx, float fracy,
						    float alpha, int nx, int ny){
    extern __shared__ float cachein[];//caching.
    const int ind=threadIdx.x+threadIdx.y*blockDim.x;
    cachein[ind]=0;
    __syncthreads();
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<ny; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<nx; ix+=stepx){
	    int ix2=ix*2;
	    int iy2=iy*2;

	    float v_0_0=alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
				 +out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(1.f-fracy)
			       +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
				 +out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f-fracy));
	    float v_1_0=alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
				 +out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(1.f-fracy)
			       +(+out[ix2+  (iy2+1)*nxout]*(fracx)
				 +out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f-fracy));
	    float v_0_1=alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
				 +out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(fracy)
			       +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
				 +out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f+fracy));
	    float v_1_1=alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
				 +out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(fracy)
			       +(+out[ix2+  (iy2+1)*nxout]*(fracx)
				 +out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f+fracy));
	    atomicAdd(&cachein[ind], v_0_0);
	    if(threadIdx.x+1<blockDim.x){
		atomicAdd(&cachein[ind+1], v_1_0);
		if(threadIdx.y+1<blockDim.y){
		    atomicAdd(&cachein[ind+1+blockDim.x], v_1_1);
		}else{
		    atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1);
		}
	    }else{
		atomicAdd(&in[ix+1      +iy*nxin], v_1_0);
		atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1);
	    }
	    if(threadIdx.y+1<blockDim.y){
		atomicAdd(&cachein[ind+blockDim.x], v_0_1);
	    }else{
		atomicAdd(&in[ix    +(iy+1)*nxin], v_0_1);
	    }
	    
	    //atomicAdd(&in[ix        +iy*nxin], v_0_0);
	    //atomicAdd(&in[ix+1      +iy*nxin], v_1_0);
	    //atomicAdd(&in[ix    +(iy+1)*nxin], v_0_1);
	    //atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1);
	    
	    __syncthreads();
	    atomicAdd(&in[ix+iy*nxin], cachein[ind]);
	    cachein[ind]=0;
	    __syncthreads();
	}
    }
}	    
/*
  Ray tracing with over sampling. Forward, from out to in. out is xloc, in is
  ploc. confirmed to agree with HXW'.  */
__global__ static void prop_grid_os2_trans_do(float *restrict out, int nxout,
					      float *restrict in, int nxin, 
					      float fracx, float fracy,
					      float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<ny; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<nx; ix+=stepx){
	    int ix2=ix*2;
	    int iy2=iy*2;
	    atomicAdd(&in[ix        +iy*nxin], alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
							+out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(1.f-fracy)
						      +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
							+out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f-fracy)));
	    
	    atomicAdd(&in[ix+1      +iy*nxin], alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
							+out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(1.f-fracy)
						      +(+out[ix2+  (iy2+1)*nxout]*(fracx)
							+out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f-fracy)));
	    
	    atomicAdd(&in[ix    +(iy+1)*nxin], alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
							+out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(fracy)
						      +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
							+out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f+fracy)));
	    
	    atomicAdd(&in[(ix+1)+(iy+1)*nxin], alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
							+out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(fracy)
						      +(+out[ix2+  (iy2+1)*nxout]*(fracx)
							+out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f+fracy)));
	    
	    /*float in00=in[ix  +iy*nxin]*alpha;
	      float in01=in[ix+1+iy*nxin]*alpha;
	      float in10=in[ix    +(iy+1)*nxin]*alpha;
	      float in11=in[(ix+1)+(iy+1)*nxin]*alpha;
	      out[ix2+iy2*nxout]+=(+(in00*(1.f-fracx)+in01*fracx)*(1.f-fracy)
	      +(in10*(1.f-fracx)+in11*fracx)*fracy);
	      out[ix2+1+iy2*nxout]+=(+(in00*(0.5f-fracx)+in01*(0.5f+fracx))*(1.f-fracy)
	      +(in10*(0.5f-fracx)+in11*(0.5f+fracx))*fracy);
	      out[ix2+(iy2+1)*nxout]+=(+(in00*(1.f-fracx)+in01*fracx)*(0.5f-fracy)
	      +(in10*(1.f-fracx)+in11*fracx)*(0.5f+fracy));
	      out[ix2+1+(iy2+1)*nxout]+=(+(in00*(0.5f-fracx)+in01*(0.5f+fracx))*(0.5f-fracy)
	      +(in10*(0.5f-fracx)+in11*(0.5f+fracx))*(0.5f+fracy));*/
	}
    }
}
/*
  __global__ static void prop_grid_os2_old_do(float *restrict out, int nxout,
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
  #pragma unroll 
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
    if(fabs(ratio-1.f)<1.e-5 && trans=='t'){
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
    int offx1=0, offy1=0;//for output. fine sampling.
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
    int nx=(int)floorf((nxi-1-dispx-1e-5)*ratio1)+1;
    int ny=(int)floorf((nyi-1-dispy-1e-5)*ratio1)+1;

    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
    int offx2=(int)floorf(dispx); dispx-=offx2;//for input. coarse sampling.
    int offy2=(int)floorf(dispy); dispy-=offy2;
    if(trans=='n'){
	if(fabs(ratio-1.f)<1.e-5){
	    prop_grid_match_do<<<DIM2(nx, ny, 32), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx, ny);
	}else{
	    prop_grid_nomatch_do<<<DIM2(nx, ny, 32), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, ratio, alpha, nx, ny);
	}
    }else if(trans=='t'){
	if(fabs(ratio-1.f)<1.e-5){
	    error("Please revert the input/output and call with trans='n'\n");
	}else if(fabs(ratio-0.5f)<1.e-5){
	    if(dispy>0.5f){//do a single col first.
		prop_grid_os2_trans_col_do<<<DIM2(nx,1,256),0,stream>>>
		    (out->p+offy1*nxo+offx1, nxo,
		     in->p+offy2*nxi+offx2, nxi, 
		     dispx, dispy, alpha, nx);
		dispy-=0.5f;
		offy1+=1; ny-=1;
		offy2+=1;
	    }
	 
	    if(dispx>0.5f){//do a single row first
		prop_grid_os2_trans_row_do<<<DIM2(1,ny,256),0,stream>>>
		    (out->p+offy1*nxo+offx1, nxo,
		     in->p+offy2*nxi+offx2, nxi, 
		     dispx, dispy, alpha, ny);
		dispx-=0.5f;
		offx1+=1; nx-=1;
		offx2+=1;
	    }
	    int nx2=nx>>1; 
	    int ny2=ny>>1;
	    
#define BS 16 // cannot be 32.
	    prop_grid_os2_trans_share_do<<<DIM2(nx2, ny2, BS), (BS)*(BS)*sizeof(float), stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx2, ny2);
	    if(ny & 1 == 1){//do the last col that is left over
		prop_grid_os2_trans_col_do<<<DIM2(nx,1,256),0,stream>>>
		    (out->p+(offy1+ny2*2)*nxo+offx1, nxo,
		     in->p+(offy2+ny2)*nxi+offx2, nxi, 
		     dispx, dispy, alpha, nx);
	    }
	    if(nx & 1 == 1){//do the last row that is left over
		prop_grid_os2_trans_row_do<<<DIM2(1,ny,256),0,stream>>>
		    (out->p+offy1*nxo+(offx1+nx2*2), nxo,
		     in->p+offy2*nxi+(offx2+nx2), nxi, 
		     dispx, dispy, alpha, ny);
	    }
	}else{
	    prop_grid_nomatch_trans_do<<<DIM2(nx, ny, 32), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, ratio, alpha, nx, ny);
	}
    }else{
	error("Invalid trans=%c\n", trans);
    }
}
