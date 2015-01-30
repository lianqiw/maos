/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
extern "C"
{
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "cudata.h"
/**
   First written on 2011-07.
   Port wfsgrad.c perfevl, etc to GPU.

   Change log:

   1) memcpy failed for cc when declared with Real cc[5] with invalid argument
   when using cudaMemcpyDeviceToHost. Changed to cudaMemcpyDefault works. But
   the consequency is that cc in kernel have only zero value. When use Real*
   cc, the memcpy succeed, but cc is null in kernel. It seems we can only pass
   cc into the kernel.

   2) copying DM information to cuda messes up atm because gpu_dm2gpu used texRefatm.


   2012-02-24
   Retired ATM_TEXTURE since it does not improve performance. map_t is not a child class of curmat
*/
static int WRAP_ATM;

/*
  Question: when I declare CC as a Real* here and initialize once, I get
  segmentation error due to cc appear as NULL using cuda-gdb.  

  Serious problem found: At time step 1, ray tracing through atm works well with
  stored memory. But at step 2, the ray tracing gives zero value. Implies memory
  has been cleared. Same problem with cc. it is cleared to all zeros.

*/


typedef struct{
    int ips;
    Real *next_atm;
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
    typedef Real pout_t[nx0];
    pout_t *pout=(pout_t*)data->next_atm;
    for(int iy=0; iy<my; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(Real)pin[iy+offy][ix+offx];
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(Real)pin[iy+offy][ix+offx-nxi];
	}
    }
    for(int iy=my; iy<ny0; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(Real)pin[iy+offy-nyi][ix+offx];
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(Real)pin[iy+offy-nyi][ix+offx-nxi];
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
    TIC;tic;
    for(int im=0; im<NGPU; im++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
    {
	gpu_set(im);
	gpu_print_mem("atm in full");
	cudata->nps=nps;
	cp2gpu(&cudata->atm, atm, nps);
	gpu_print_mem("atm out");
    }
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
    toc2("atm to gpu");
}
/**
   Transfer atmosphere or update atmosphere in GPU.
*/
void gpu_atm2gpu(const mapcell *atmc, const dmat *atmscale, const PARMS_T *parms, int iseed, int isim){
    if(!atmc) return;
    map_t **atm=atmc->p;
    const int nps=parms->atm.nps;
    static int nx0=0,ny0=0;
    static int iseed0=-1;
    if(iseed0!=iseed && atmscale){
	dfree(cudata_t::atmscale);
	if(atmscale) cudata_t::atmscale=ddup(atmscale);
    }
    if(!nx0){
	/*The minimum size to cover the meta-pupil*/
	long nxn=parms->atm.nxn;
	long nyn=parms->atm.nyn;
	
	long avail_min=LONG_MAX, avail_max=0;
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    long availi=gpu_get_mem();
	    if(avail_min>availi){
		avail_min=availi;
	    }
	    if(avail_max<availi){
		avail_max=availi;
	    }
	}
	avail_min-=256<<20; /*reserve 256 MiB*/
	long need=nps*sizeof(Real)*nxn*nyn;
	info2("Min atm is %ldx%ld, available memory is %ld~%ld MB, need at least %ldMB\n", 
	      nxn, nyn, avail_min>>20, avail_max>>20, need>>20);
	if(avail_min<need){
	    if(avail_max<need){
		error2("No GPU has enough memory\n");
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
	    long nxa=(long)floor(sqrt((avail_min)/nps/sizeof(Real)));
	    info2("GPU can host %d %ldx%ld atmosphere\n", nps, nxa, nxa);
	    if(nxa*nxa>parms->atm.nx*parms->atm.ny){/*we can host all atmosphere. */
		nx0=parms->atm.nx;
		ny0=parms->atm.ny;
	    }else{
		nx0=MIN(nxa, nxn*2);
		ny0=MIN(nxa, nyn*2);
	    }
	    info2("We will host %dx%d in GPU, taking %zd MiB\n", 
		  nx0, ny0, (nx0*ny0*nps*sizeof(Real))>>20);
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
    static Real *next_ox=NULL;
    static Real *next_oy=NULL;
    static Real **next_atm=NULL;
    static pthread_t *next_threads=NULL;
    if(need_init){
	need_init=0;
	/*The atm in GPU is smaller than in CPU. */
	next_isim=(int*)calloc(nps, sizeof(int));
	next_ox=(Real*)calloc(nps, sizeof(Real));
	next_oy=(Real*)calloc(nps, sizeof(Real));
	next_threads=(pthread_t*)calloc(nps, sizeof(pthread_t));
	next_atm=(Real **)calloc(nps, sizeof(void*));

	for(int im=0; im<NGPU; im++){/*Loop over all GPUs. */
	    gpu_set(im);
	    cudata->atm=new cumap_t[nps];
	    cudata->nps=nps;
	    for(int ips=0; ips<nps; ips++){
		cudata->atm[ips].p=new curmat(nx0, ny0);
		cudata->atm[ips].nx=nx0;
		cudata->atm[ips].ny=ny0;

	    }
	}/*for im */
    }/*if need_init; */
    const double dt=parms->sim.dt;
    const double dx=parms->atm.dx;
    if(iseed0!=iseed){/*A new seed or initialization update vx, vy, ht, etc. */
	iseed0=iseed;

    	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    cumap_t *cuatm=cudata->atm;
	    for(int ips=0; ips<nps; ips++){
		/*Do not copy over nx, ny from atm as cuatm is smaller*/
		cuatm[ips].vx=atm[ips]->vx;
		cuatm[ips].vy=atm[ips]->vy;
		cuatm[ips].ht=atm[ips]->h;
		cuatm[ips].dx=atm[ips]->dx;
		cuatm[ips].dy=atm[ips]->dx;
		cuatm[ips].ox=INFINITY;/*place holder */
		cuatm[ips].oy=INFINITY;
	    }
	}
	for(int ips=0; ips<nps; ips++){
	    if(next_atm[ips]){
		free(next_atm[ips]); next_atm[ips]=NULL;
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
    for(int ips=0; ips<nps; ips++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
    {
	/*Load atmosphere to Real memory in advance. This takes time if atm is
	  stored in file.*/
	if(isim>next_isim[ips]-100 && !next_atm[ips]){
	    /*pinned memory is faster for copying to GPU. */
	    next_atm[ips]=(Real*)malloc(sizeof(Real)*nx0*ny0);
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
	    for(int im=0; im<NGPU; im++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
	    {
		gpu_set(im);
		cumap_t *cuatm=cudata->atm;
		cuatm[ips].ox=next_ox[ips];
		cuatm[ips].oy=next_oy[ips];
		DO(cudaMemcpy(cuatm[ips].p->p, (Real*)next_atm[ips],
			      nx0*ny0*sizeof(Real), cudaMemcpyHostToDevice));
	    }/*for im */
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
	    free(next_atm[ips]);next_atm[ips]=0;
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
	    toc2("Step %d: Layer %d transfered. next in step %d. ",isim, ips, next_isim[ips]); 
}//if isim
}//for ips
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
}

/**
   Copy DM configurations to GPU.
*/
curmat* gpu_dmcubic_cc(Real iac){
    Real cc[5];
    Real cubicn=1.f/(1.f+2.f*iac);
    cc[0]=1.f*cubicn;
    cc[1]=(4.f*iac-2.5f)*cubicn; 
    cc[2]=(1.5f-3.f*iac)*cubicn;		       
    cc[3]=(2.f*iac-0.5f)*cubicn;			
    cc[4]=(0.5f-iac)*cubicn; 
    curmat *res=curnew(5,1);
    cudaMemcpy(res->p, cc, 5*sizeof(Real), cudaMemcpyHostToDevice);
    return res;
}
/**
   Copy DM commands to GPU.
*/
static void gpu_dm2gpu(cumap_t **cudm, map_t **dmreal, int ndm, DM_CFG_T *dmcfg){
    cp2gpu(cudm, dmreal, ndm);
    if(dmcfg){
	for(int idm=0; idm<ndm; idm++){
	    if(dmcfg[idm].cubic && !(*cudm)[idm].cubic_cc){
		(*cudm)[idm].cubic_cc=gpu_dmcubic_cc(dmcfg[idm].iac);
	    }
	    (*cudm)[idm].ht+=dmcfg[idm].vmisreg;
	}
    }
}
void gpu_dmreal2gpu(mapcell *dmreal, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++)
#if _OPENMP >= 200805
#pragma omp task
#endif
    {
	gpu_set(im);
	cudata->ndm=dmreal->nx;
	gpu_dm2gpu(&cudata->dmreal, dmreal->p, dmreal->nx, dmcfg);
    }
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
}
void gpu_dmproj2gpu(mapcell *dmproj, DM_CFG_T *dmcfg){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->ndm=dmproj->nx;
	gpu_dm2gpu(&cudata->dmproj, dmproj->p, dmproj->nx, dmcfg);
    }
}

#define KARG_COMMON const Real (*restrict loc)[2], const int nloc, const Real dxi, const Real dyi, const Real dispx, const Real dispy, const Real alpha

/*This is memory bound. So increasing # of points processed does not help. */
__global__ void prop_linear(Real *restrict out, const Real *restrict in, const int nx, const int ny,
			    KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real x=loc[i][0]*dxi+dispx;
	Real y=loc[i][1]*dyi+dispy;
	int ix=Z(floor)(x);
	int iy=Z(floor)(y);
	x=x-ix; y=y-iy;
	if(ix>=0 && ix<nx-1 && iy>=0 && iy<ny-1){
	    Real tmp=((+in[iy*nx+ix]*(1.f-x) +in[iy*nx+ix+1]*x)*(1.f-y)
		      +(+in[(iy+1)*nx+ix]*(1.f-x) +in[(iy+1)*nx+ix+1]*x)*y);
	    add_valid(out[i], alpha, tmp);
	}
    }
}
/*This is memory bound. So increasing # of points processed does not help. */
__global__ void prop_linear_nocheck(Real *restrict out, const Real *restrict in, 
				    const int nx, const int ny, KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real x=loc[i][0]*dxi+dispx;
	Real y=loc[i][1]*dyi+dispy;
	int ix=Z(floor)(x);
	int iy=Z(floor)(y);
	x=x-ix; y=y-iy;
	out[i]+=alpha*((in[iy*nx+ix]*(1-x)+in[iy*nx+ix+1]*x)*(1-y)
		       +(in[(iy+1)*nx+ix]*(1-x)+in[(iy+1)*nx+ix+1]*x)*y);
    }
}
/*This is memory bound. So increasing # of points processed does not help. */
__global__ void prop_linear_wrap(Real *restrict out, const Real *restrict in, const int nx, const int ny,
				 KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real x=loc[i][0]*dxi+dispx;
	Real y=loc[i][1]*dyi+dispy;
	int ix=Z(floor)(x);
	int iy=Z(floor)(y);
	x=x-ix; y=y-iy;
	while(ix<0) ix+=nx; 
	while(iy<0) iy+=ny;
	while(ix>nx-1) ix-=nx; 
	while(iy>ny-1) iy-=ny;
	int ix1=(ix==nx-1)?0:(ix+1);
	int iy1=(iy==ny-1)?0:(iy+1);
	out[i]+=alpha*((in[iy*nx+ix]*(1.f-x)+in[iy*nx+ix1]*x)*(1.f-y)
		       +(in[(iy1)*nx+ix]*(1.f-x)+in[(iy1)*nx+ix1]*x)*y);
    }
}

/*This is memory bound. So increasing # of points processed does not help. */
__global__ void prop_cubic(Real *restrict out, const Real *restrict in, const int nx, const int ny,
			   KARG_COMMON, const Real *cc){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real x=loc[i][0]*dxi+dispx;
	Real y=loc[i][1]*dyi+dispy;
	int ix=Z(floor)(x); x=x-ix;
	int iy=Z(floor)(y); y=y-iy;
	Real fx[4],fy;
	Real sum=0;
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
	add_valid(out[i], alpha, sum);
    }
}

/**
   Ray tracing of atm.
*/
void gpu_atm2loc(Real *phiout, culoc_t *loc, const Real hs, const Real thetax,const Real thetay,
		 const Real mispx, const Real mispy, const Real dt, const int isim, Real atmalpha, cudaStream_t stream){
    cumap_t *cuatm=cudata->atm;
    if(Z(fabs)(atmalpha)<EPS) return;
    if(cudata_t::atmscale){
	atmalpha*=cudata_t::atmscale->p[isim];
    }
    for(int ips=0; ips<cudata->nps; ips++){
	const Real dx=cuatm[ips].dx;
	const Real dy=cuatm[ips].dy;
	const Real ht=cuatm[ips].ht;
	const Real vx=cuatm[ips].vx;
	const Real vy=cuatm[ips].vy;
	const Real dispx=(ht*thetax+mispx-vx*dt*isim-cuatm[ips].ox)/dx;
	const Real dispy=(ht*thetay+mispy-vy*dt*isim-cuatm[ips].oy)/dy;
	const Real scale=1.f-ht/hs;
	const int nloc=loc->nloc;
#define COMM loc->p,loc->nloc,scale/dx,scale/dy, dispx, dispy, atmalpha
	if(WRAP_ATM){
	    prop_linear_wrap<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips].p->p, cuatm[ips].nx, cuatm[ips].ny, COMM);
	}else{/*we are gauranteed. */
	    prop_linear_nocheck<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips].p->p, cuatm[ips].nx, cuatm[ips].ny, COMM);
	}
#undef COMM
    }    
}

/**
   Ray tracing of dm
*/
void gpu_dm2loc(Real *phiout, culoc_t **locarr, cumap_t *cudm, int ndm,
		const Real hs, const Real thetax, const Real thetay,
		const Real mispx, const Real mispy, const Real dmalpha, cudaStream_t stream){
    if(Z(fabs)(dmalpha)<EPS) return;
    for(int idm=0; idm<ndm; idm++){
	assert(cudm[idm].ny>1);//prevent accidentally pass in a vector
	const Real dx=cudm[idm].dx;
	const Real dy=cudm[idm].dy;
	const Real ht=cudm[idm].ht;
	const Real dispx=(ht*thetax+mispx-cudm[idm].ox)/dx;
	const Real dispy=(ht*thetay+mispy-cudm[idm].oy)/dy;
	const Real scale=1.f-ht/hs;
	const Real (*loc)[2]=locarr[idm]->p;
	const int nloc=locarr[idm]->nloc;
#define COMM loc, nloc,scale/dx,scale/dy, dispx, dispy, dmalpha
#define KARG cudm[idm].p->p,cudm[idm].nx,cudm[idm].ny, COMM
	if (cudm[idm].cubic_cc){//128 is a good number for cubic. 
	    prop_cubic<<<DIM(nloc,128), 0, stream>>>(phiout, KARG, cudm[idm].cubic_cc->p);
	}else{
	    prop_linear<<<DIM(nloc,256), 0, stream>>>(phiout, KARG);
	}
#undef KARG
#undef COMM
    }/*idm */
}

/**
   Convert NGS mode vector to aperture grid for science directions.  */
void gpu_ngsmod2science(curmat *opd, Real (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream){
    if(ngsmod->nmod==2){
	curaddptt(opd, loc, 0, mod[0]*alpha, mod[1]*alpha, stream);
    }else{
	const Real ht=ngsmod->ht;
	const Real scale=ngsmod->scale;
	Real focus;
	if(ngsmod->nmod>5){
	    focus=mod[5];
	    if(!ngsmod->ahstfocus){
		focus+=mod[2]*(1.-scale);
	    }
	}else{
	    focus=mod[2]*(1.-scale);
	}
	add_ngsmod_do<<<DIM(opd->nx*opd->ny, 256), 0, stream>>>
	    (opd->p, loc, opd->nx*opd->ny, 
	     mod[0], mod[1], mod[2], mod[3], mod[4], focus,
	     thetax, thetay, scale, ht, alpha);
    }
}
