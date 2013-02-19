/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <cuda.h>
extern "C"
{
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"

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


   2012-02-24
   Retired ATM_TEXTURE since it does not improve performance. map_t is not a child class of curmat
*/
static int WRAP_ATM;

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
    int isim;
    int isim_next;
    map_t *atm;
}atm_prep_t;
static void atm_prep(atm_prep_t *data){
    PNEW(lock);
    const int ips=data->ips;
    LOCK(lock);/*make sure we only read one layer at a time. */
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
    toc2("Step %d: Layer %d: Preparing atm for step %d", data->isim, ips, data->isim_next);
    UNLOCK(lock);
    free(data);/*allocated in parent thread. we free it here. */
}
/**
   Transfer atmospheric data to GPU.
*/
static void gpu_atm2gpu_full(map_t **atm, int nps){
    if(!atm) return;
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	gpu_print_mem("atm in full");
	TIC;tic;
	cudata->nps=nps;
	cp2gpu(&cudata->atm, atm, nps);
	toc2("atm to gpu full");/*0.4 second. */
	gpu_print_mem("atm out");
    }
}
/**
   Transfer atmosphere or update atmosphere in GPU.
*/
void gpu_atm2gpu(map_t **atm, const PARMS_T *parms, int iseed, int isim){
    if(!atm) return;
    if(parms->atm.evolve){
	TO_IMPLEMENT;
	/*test whether screen changed. transfer if changed. */
    }
    const int nps=parms->atm.nps;
    static int nx0=0,ny0=0;
    static int iseed0=-1;
    if(!nx0){
	/*The minimum size to cover the meta-pupil*/
	long nxn=parms->atm.nxn;
	long nyn=parms->atm.nyn;
	
	long avail_min, avail_max;
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    long availi=gpu_get_mem();
	    if(igpu==0){
		avail_min=availi;
		avail_max=availi;
	    }else{
		if(avail_min>availi){
		    avail_min=availi;
		}
		if(avail_max<availi){
		    avail_max=availi;
		}
	    }
	}
	long spare=300*1024*1024;
	long need=spare+nps*sizeof(float)*nxn*nyn;
	info2("Available memory is %ld (min) %ld (max). Min atm is %ldx%ld, need %ld\n", avail_min, avail_max, nxn, nyn, need);
	if(avail_min<need){
	    if(avail_max<need){
		error("ALL GPUs does not have enough memory\n");
	    }else{
		char *gcmd=NULL;
		for(int igpu=0; igpu<NGPU; igpu++){
		    extern int *GPUS;
		    gpu_set(igpu);
		    if(gpu_get_mem()>need){
			char tmp[8];
			snprintf(tmp,8," -g%d", GPUS[igpu]);
			char *tmp2=gcmd;
			gcmd=stradd(tmp, tmp2, (void*)NULL);
			free(tmp2);
		    }
		}
		error("Please rerun maos with %s\n", gcmd);
	    }
	    _Exit(0);
	}else{
	    /*we are able to host this amount. */
	    long nxa=(long)roundf(sqrt((avail_min-spare)/nps/sizeof(float)));
	    info2("GPU can host %d %ldx%ld atmosphere\n", nps, nxa, nxa);
	    if(nxa*nxa>parms->atm.nx*parms->atm.ny){/*we can host all atmosphere. */
		nx0=parms->atm.nx;
		ny0=parms->atm.ny;
	    }else{
		nx0=nxa;
		ny0=nxa;
	    }
	    info2("We will host %dx%d in GPU\n", nx0, ny0);
	}
    }
    /*The atm in GPU is the same as in CPU. */
    if(nx0==parms->atm.nx && ny0==parms->atm.ny){
	WRAP_ATM=1;
	if(iseed0!=iseed){
	    gpu_atm2gpu_full(atm, nps);
	    iseed0=iseed;
	}
	return;
    }
    WRAP_ATM=0;
    static int need_init=1;
    static int *next_isim=NULL;/*next time step to update this atmosphere. */
    static float *next_ox=NULL;
    static float *next_oy=NULL;
    static float **next_atm=NULL;
    static pthread_t *next_threads=NULL;
    if(need_init){
	need_init=0;
	/*The atm in GPU is smaller than in CPU. */
	next_isim=(int*)calloc(nps, sizeof(int));
	next_ox=(float*)calloc(nps, sizeof(float));
	next_oy=(float*)calloc(nps, sizeof(float));
	next_threads=(pthread_t*)calloc(nps, sizeof(pthread_t));
	next_atm=(float **)calloc(nps, sizeof(void*));

	for(int im=0; im<NGPU; im++){/*Loop over all GPUs. */
	    gpu_set(im);
	    cudata->atm=(cumap_t**)calloc(nps, sizeof(cumap_t*));
	    cudata->nps=nps;
	    for(int ips=0; ips<nps; ips++){
		cudata->atm[ips]=new cumap_t(nx0, ny0);
	    }
	}/*for im */
    }/*if need_init; */
    const double dt=parms->sim.dt;
    const double dx=parms->atm.dx;
    if(iseed0!=iseed){/*A new seed or initialization update vx, vy, ht, etc. */
	iseed0=iseed;
    	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    cumap_t **cuatm=cudata->atm;
	    for(int ips=0; ips<nps; ips++){
		cuatm[ips]->vx=atm[ips]->vx;
		cuatm[ips]->vy=atm[ips]->vy;
		cuatm[ips]->ht=atm[ips]->h;
		cuatm[ips]->dx=atm[ips]->dx;
		cuatm[ips]->ox=INFINITY;/*place holder */
		cuatm[ips]->oy=INFINITY;
	    }
	}
	for(int ips=0; ips<nps; ips++){
	    if(next_atm[ips]){
		cudaFreeHost(next_atm[ips]);
		next_atm[ips]=NULL;
	    }
	    next_isim[ips]=isim;/*right now. */
	    /*copy from below. */
	    if(atm[ips]->vx>0){/*align right */
		next_ox[ips]=(parms->atm.nxn/2+1-nx0)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }else{/*align left */
		next_ox[ips]=(-parms->atm.nxn/2)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }
	    if(atm[ips]->vy>0){/*align right */
		next_oy[ips]=(parms->atm.nyn/2+1-ny0)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }else{/*align left */
		next_oy[ips]=(-parms->atm.nyn/2)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }
	}
    }
    for(int ips=0; ips<nps; ips++){
	/*Load atmosphere to float memory in advance. This takes time if atm is
	  stored in file.*/
	if(isim>next_isim[ips]-100 && !next_atm[ips]){
	    /*pinned memory is faster for copying to GPU. */
	    cudaMallocHost(&(next_atm[ips]), sizeof(float)*nx0*ny0);
	    const int nxi=atm[ips]->nx;
	    const int nyi=atm[ips]->ny;
	    int offx=(int)round((next_ox[ips]-atm[ips]->ox)/dx);
	    int offy=(int)round((next_oy[ips]-atm[ips]->oy)/dx);
	    next_ox[ips]=atm[ips]->ox+offx*dx;
	    next_oy[ips]=atm[ips]->oy+offy*dx;
	    offx=offx%nxi; if(offx<0) offx+=nxi;
	    offy=offy%nyi; if(offy<0) offy+=nyi;
	    atm_prep_t *data=(atm_prep_t*)calloc(1, sizeof(atm_prep_t));/*cannot use local variable. */
	    data->ips=ips;
	    data->next_atm=next_atm[ips];
	    data->nx0=nx0;
	    data->ny0=ny0;
	    data->offx=offx;
	    data->offy=offy;
	    data->isim=isim;
	    data->isim_next=next_isim[ips];
	    data->atm=atm[ips];
	    /*launch an independent thread to pull data in. thread will exit when it is done. */
	    pthread_create(&next_threads[ips], NULL, (void *(*)(void *))atm_prep, data);
	}/*need to cpy atm to next_atm. */
	if(isim==next_isim[ips]){
	    /*need to copy atm to gpu. and update next_isim */
	    TIC;tic;
	    pthread_join(next_threads[ips], NULL);
	    toc2("Step %d: Layer %d wait for transfering",isim, ips);
	    for(int im=0; im<NGPU; im++){
		tic;
		gpu_set(im);
		cumap_t **cuatm=cudata->atm;
		cuatm[ips]->ox=next_ox[ips];
		cuatm[ips]->oy=next_oy[ips];
		DO(cudaMemcpy(cuatm[ips]->p, (float*)next_atm[ips], nx0*ny0*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_SYNC_DEVICE;
		int offx=(int)round((next_ox[ips]-atm[ips]->ox)/dx);
		int offy=(int)round((next_oy[ips]-atm[ips]->oy)/dx);
		toc2("Step %d: Copying layer %d size %dx%d to GPU %d: offx=%d, offy=%d", 
		     isim, ips, nx0, ny0, im, offx, offy);tic;
	    }/*for im */
	    cudaFreeHost(next_atm[ips]);
	    next_atm[ips]=NULL;
	    /*Update next_isim. */
	    long isim1, isim2;
	    if(fabs(atm[ips]->vx)<EPS){
		isim1=INT_MAX;
	    }else if(atm[ips]->vx>0){/*align right. */
		isim1=(long)floor(-(next_ox[ips]+(parms->atm.nxn/2)*dx)/(atm[ips]->vx*dt));
	    }else{/*align left */
		isim1=(long)floor(-(next_ox[ips]+(nx0-(parms->atm.nxn/2+1))*dx)/(atm[ips]->vx*dt));
	    }
	    if(fabs(atm[ips]->vy)<EPS){
		isim2=INT_MAX;
	    }else if(atm[ips]->vy>0){/*align right. */
		isim2=(long)floor(-(next_oy[ips]+(parms->atm.nyn/2)*dx)/(atm[ips]->vy*dt));
	    }else{/*align left */
		isim2=(long)floor(-(next_oy[ips]+(ny0-(parms->atm.nyn/2+1))*dx)/(atm[ips]->vy*dt));
	    }
	    next_isim[ips]=isim1<isim2?isim1:isim2;
	    if(next_isim[ips]>parms->sim.end){/*no need to do */
		next_isim[ips]=INT_MAX;
	    }
	    if(atm[ips]->vx>0){/*align right */
		next_ox[ips]=(parms->atm.nxn/2+1-nx0)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }else{/*align left */
		next_ox[ips]=(-parms->atm.nxn/2)*dx-atm[ips]->vx*dt*next_isim[ips];
	    }
	    if(atm[ips]->vy>0){/*align right */
		next_oy[ips]=(parms->atm.nyn/2+1-ny0)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }else{/*align left */
		next_oy[ips]=(-parms->atm.nyn/2)*dx-atm[ips]->vy*dt*next_isim[ips];
	    }
	    info2("Step %d: next update layer %d in step %d\n", isim, ips, next_isim[ips]);
	}
    }
}

/**
   Copy DM configurations to GPU.
*/
float* gpu_dmcubic_cc(float iac){
    float cc[5];
    float cubicn=1.f/(1.f+2.f*iac);
    cc[0]=1.f*cubicn;
    cc[1]=(4.f*iac-2.5f)*cubicn; 
    cc[2]=(1.5f-3.f*iac)*cubicn;		       
    cc[3]=(2.f*iac-0.5f)*cubicn;			
    cc[4]=(0.5f-iac)*cubicn; 
    float *ret=NULL;
    DO(cudaMalloc(&ret, 5*sizeof(float)));
    cudaMemcpy(ret, cc, 5*sizeof(float), cudaMemcpyHostToDevice);
    return ret;
}
/**
   Copy DM commands to GPU.
*/
static void gpu_dm2gpu(cumap_t ***cudm, map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    cp2gpu(cudm, dmreal, ndm);
    if(dmcfg){
	for(int idm=0; idm<ndm; idm++){
	    if(dmcfg[idm].cubic && !(*cudm)[idm]->cubic_cc){
		(*cudm)[idm]->cubic_cc=gpu_dmcubic_cc(dmcfg[idm].iac);
	    }
	}
    }
    CUDA_SYNC_DEVICE;
}
void gpu_dmreal2gpu(map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->ndm=ndm;
	gpu_dm2gpu(&cudata->dmreal, dmreal, ndm, dmcfg);
    }
}
void gpu_dmproj2gpu(map_t **dmproj, int ndm, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->ndm=ndm;
	gpu_dm2gpu(&cudata->dmproj, dmproj, ndm, dmcfg);
    }
}

#define KARG_COMMON const float (*restrict loc)[2], const int nloc, const float dx, const float dy, const float dispx, const float dispy, const float alpha

/*This is memory bound. So increasing # of points processed does not help. */
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
/*This is memory bound. So increasing # of points processed does not help. */
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
/*This is memory bound. So increasing # of points processed does not help. */
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

/*This is memory bound. So increasing # of points processed does not help. */
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
	    continue;/*out of range. */
	}
	/*cc need to be in device memory for sm_13 to work.*/
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
    cumap_t **cuatm=cudata->atm;
    if(fabs(atmalpha)<EPS) return;
    for(int ips=0; ips<cudata->nps; ips++){
	const float dx=cuatm[ips]->dx;
	const float du=1.f/dx;
	const float ht=cuatm[ips]->ht;
	const float vx=cuatm[ips]->vx;
	const float vy=cuatm[ips]->vy;
	const float dispx=(ht*thetax+mispx-vx*dtisim-cuatm[ips]->ox)*du;
	const float dispy=(ht*thetay+mispy-vy*dtisim-cuatm[ips]->oy)*du;
	const float scale=1.f-ht/hs;
#define COMM loc,nloc,scale*du,scale*du, dispx, dispy, atmalpha
	if(WRAP_ATM){
	    prop_linear_wrap<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips]->p, cuatm[ips]->nx, cuatm[ips]->ny, COMM);
	}else{/*we are gauranteed. */
	    prop_linear_nocheck<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips]->p, cuatm[ips]->nx, cuatm[ips]->ny, COMM);
	}
#undef COMM
    }    
}

/**
   Ray tracing of dm
*/
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, cumap_t **cudm, int ndm,
		const float hs, const float thetax, const float thetay,
		const float mispx, const float mispy, const float dmalpha, cudaStream_t stream){
    if(fabs(dmalpha)<EPS) return;
    for(int idm=0; idm<ndm; idm++){
	assert(cudm[idm]->ny>1);//prevent accidentally pass in a vector
	const float dx=cudm[idm]->dx;
	const float du=1.f/dx;
	const float ht=cudm[idm]->ht;
	const float dispx=(ht*thetax+mispx-cudm[idm]->ox)*du;
	const float dispy=(ht*thetay+mispy-cudm[idm]->oy)*du;
	const float scale=1.f-ht/hs;
#define COMM loc,nloc,scale*du,scale*du, dispx, dispy, dmalpha
#define KARG cudm[idm]->p,cudm[idm]->nx,cudm[idm]->ny, COMM
	if (cudm[idm]->cubic_cc){//128 is a good number for cubic. 
	    prop_cubic<<<DIM(nloc,128), 0, stream>>>(phiout, KARG, cudm[idm]->cubic_cc);
	}else{
	    prop_linear<<<DIM(nloc,256), 0, stream>>>(phiout, KARG);
	}
#undef KARG
#undef COMM
    }/*idm */
}
