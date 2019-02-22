/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
    int nx0;
    int ny0;
    int offx;
    int offy;
    int isim;
    int isim_next;
    map_t *atm;
    Real ox;
    Real oy;
    Real *next_atm;
    pthread_t threads;
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
    typedef Real pout_t[nx0];
    pout_t *pout=(pout_t*)data->next_atm;
    for(int iy=0; iy<my; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(Real)P(atm, ix+offx, iy+offy);
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(Real)P(atm, ix+offx-nxi, iy+offy);
	}
    }
    for(int iy=my; iy<ny0; iy++){
	for(int ix=0; ix<mx; ix++){
	    pout[iy][ix]=(Real)P(atm, ix+offx, iy+offy-nyi);
	}
	for(int ix=mx; ix<nx0; ix++){
	    pout[iy][ix]=(Real)P(atm, ix+offx-nxi, iy+offy-nyi);
	}
    }
    toc2("Step %d: Layer %d: Preparing atm for step %d", data->isim, ips, data->isim_next);
    UNLOCK(lock);
}
/**
   Transfer atmospheric data to GPU.
*/
static void gpu_atm2gpu_full(const mapcell *atm){
    if(!atm) return;
    TIC;tic;
    for(int im=0; im<NGPU; im++)
#if _OPENMP >= 200805 
#pragma omp task
#endif
    {
	gpu_set(im);
	gpu_print_mem("atm in full");
	cudata->nps=atm->nx;
	cp2gpu(cudata->atm, atm);
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
    static int iseed0=-1;
    if(iseed0!=iseed){
	dfree(cuglobal->atmscale);
	if(atmscale) cuglobal->atmscale=ddup(atmscale);
    }
    /*The minimum size to cover the meta-pupil*/
    const long nxn=parms->atm.nxnmax;
    static int nx0=0,ny0=0;
    if(!nx0){
	if(parms->dbg.fullatm){
	    info("Always host full atmosphere in GPU. Set dbg.fullatm=0 to disable\n");
	    nx0=parms->atm.nx;
	    ny0=parms->atm.ny;
	}else{
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
	    avail_min-=256<<19;//reserve 128MiB 
	    avail_max-=256<<19;
	    long need=nps*sizeof(Real)*nxn*nxn;
	    info("Min atm is %ldx%ld, available memory is %ld~%ld MB, need at least %ldMB\n", 
		 nxn, nxn, avail_min>>20, avail_max>>20, need>>20);
	    if(avail_min<need){
		if(avail_max<need){
		    error("No GPU has enough memory\n");
		}else{
		    char *gcmd=NULL;
		    for(int igpu=0; igpu<NGPU; igpu++){
			//extern Array<int> GPUS;
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
		info("GPU can host %d %ldx%ld atmosphere\n", nps, nxa, nxa);
		if(nxa*nxa>parms->atm.nx*parms->atm.ny){/*we can host all atmosphere. */
		    nx0=parms->atm.nx;
		    ny0=parms->atm.ny;
		}else{
		    nx0=MIN(nxa, nxn*2);
		    ny0=MIN(nxa, nxn*2);
		}
	    }
	    info("We will host %dx%d in GPU, taking %zd MiB\n", 
		 nx0, ny0, (nx0*ny0*nps*sizeof(Real))>>20);
	}
    }
    /*The atm in GPU is the same as in CPU. */
    if(nx0==parms->atm.nx && ny0==parms->atm.ny){
	WRAP_ATM=1;
	if(iseed0!=iseed){
	    gpu_atm2gpu_full(atmc);
	    iseed0=iseed;
	}
	return;
    }
    WRAP_ATM=0;
    const double dt=parms->sim.dt;
    const double dx=parms->atm.dx;
    static atm_prep_t *prep_data=NULL;
    /*static int *next_isim=NULL;
    static Real *next_ox=NULL;
    static Real *next_oy=NULL;
    static Real **next_atm=NULL;
    static pthread_t *next_threads=NULL;*/
    if(iseed0==-1){//Initializing data
	/*The atm in GPU is smaller than in CPU. */
	prep_data=(atm_prep_t*)calloc(nps, sizeof(atm_prep_t));
	for(int ips=0; ips<nps; ips++){
	    prep_data[ips].nx0=nx0;
	    prep_data[ips].ny0=ny0;
	}
	for(int im=0; im<NGPU; im++){/*Loop over all GPUs. */
	    gpu_set(im);
	    cudata->atm=cumapcell(nps, 1);
	    cudata->nps=nps;
	    for(int ips=0; ips<nps; ips++){
		cudata->atm[ips].p=curmat(nx0, ny0);
		cudata->atm[ips].nx=nx0;
		cudata->atm[ips].ny=ny0;
	    }
	}/*for im */
    }/*if need_init; */
    if(iseed0!=iseed){/*A new seed or initialization update vx, vy, ht, etc. */
	iseed0=iseed;
    	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    cumapcell &cuatm=cudata->atm;
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
	    if(prep_data[ips].next_atm){
		warning("next_atm is not empty\n");
		free(prep_data[ips].next_atm);
		prep_data[ips].next_atm=NULL;
	    }
	    prep_data[ips].isim_next=isim;/*right now. */
	}
    }
    for(int ips=0; ips<nps; ips++) {
	/*Load atmosphere to Real memory in advance. This takes time if atm is
	  stored in file.*/
	if(isim+100 > prep_data[ips].isim_next && !prep_data[ips].next_atm){
	    const int nxni=parms->atm.nxn->p[ips];
	    prep_data[ips].next_atm=(Real*)malloc(sizeof(Real)*nx0*ny0);
	    const int margin=1;//to avoid going out of range due to interpolation.
	    //Compute the origin of the subset phase screen to be copied to GPU.
	    if(atm[ips]->vx>0){/*align right */
		prep_data[ips].ox=(nxni/2+1-nx0+margin)*dx-atm[ips]->vx*dt*prep_data[ips].isim_next;
	    }else{/*align left */
		prep_data[ips].ox=(-nxni/2-margin)*dx-atm[ips]->vx*dt*prep_data[ips].isim_next;
	    }
	    if(atm[ips]->vy>0){/*align top */
		prep_data[ips].oy=(nxni/2+1-ny0+margin)*dx-atm[ips]->vy*dt*prep_data[ips].isim_next;
	    }else{/*align bottom */
		prep_data[ips].oy=(-nxni/2-margin)*dx-atm[ips]->vy*dt*prep_data[ips].isim_next;
	    }
	    //Compute the offset of the subset phase screen
	    int offx=(int)round((prep_data[ips].ox-atm[ips]->ox)/dx);
	    int offy=(int)round((prep_data[ips].oy-atm[ips]->oy)/dx);
	    prep_data[ips].ox=atm[ips]->ox+offx*dx;
	    prep_data[ips].oy=atm[ips]->oy+offy*dx;

	    //Wrapping
	    const int nxi=atm[ips]->nx;
	    const int nyi=atm[ips]->ny;
	    offx=offx%nxi; if(offx<0) offx+=nxi;
	    offy=offy%nyi; if(offy<0) offy+=nyi;

	    prep_data[ips].ips=ips;
	    prep_data[ips].isim=isim;
	    prep_data[ips].offx=offx;
	    prep_data[ips].offy=offy;
	    prep_data[ips].atm=atm[ips];
	    
	    /*launch an independent thread to pull data in. thread will exit when it is done. */
	    pthread_create(&prep_data[ips].threads, NULL, (void *(*)(void *))atm_prep, prep_data+ips);
	}/*need to cpy atm to next_atm. */
    }
    for(int ips=0; ips<nps; ips++) {
	if(isim==prep_data[ips].isim_next){
	    /*need to copy atm to gpu. and update next_isim */
	    TIC;tic;
	    pthread_join(prep_data[ips].threads, NULL);
	    for(int im=0; im<NGPU; im++)
	    {
		gpu_set(im);
		cumapcell &cuatm=cudata->atm;
		cuatm[ips].ox=prep_data[ips].ox;
		cuatm[ips].oy=prep_data[ips].oy;
		DO(cudaMemcpy(cuatm[ips](), prep_data[ips].next_atm,
			      nx0*ny0*sizeof(Real), cudaMemcpyHostToDevice));
	    }/*for im */
	    free(prep_data[ips].next_atm);prep_data[ips].next_atm=0;
	    /*Update next_isim. */
	    long isim1, isim2;
	    const int nxni=parms->atm.nxn->p[ips];
	    if(fabs(atm[ips]->vx)<EPS){
		isim1=INT_MAX;
	    }else if(atm[ips]->vx>0){/*align right. */
		isim1=(long)floor(-(prep_data[ips].ox+(nxni/2)*dx)/(atm[ips]->vx*dt));
	    }else{/*align left */
		isim1=(long)floor(-(prep_data[ips].ox+(nx0-(nxni/2+1))*dx)/(atm[ips]->vx*dt));
	    }
	    if(fabs(atm[ips]->vy)<EPS){
		isim2=INT_MAX;
	    }else if(atm[ips]->vy>0){/*align right. */
		isim2=(long)floor(-(prep_data[ips].oy+(nxni/2)*dx)/(atm[ips]->vy*dt));
	    }else{/*align left */
		isim2=(long)floor(-(prep_data[ips].oy+(ny0-(nxni/2+1))*dx)/(atm[ips]->vy*dt));
	    }
	    prep_data[ips].isim_next=isim1<isim2?isim1:isim2;
	    if(prep_data[ips].isim_next>parms->sim.end){/*no need to do */
		prep_data[ips].isim_next=INT_MAX;
	    }
	 
	    toc2("Step %d: Layer %d transfered. next in step %d. ",isim, ips, prep_data[ips].isim_next); 
	}//if isim
    }//for ips
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
}

/**
   Copy DM configurations to GPU.
*/
curmat gpu_dmcubic_cc(Real iac){
    if(iac){
	Real cc[5];
	Real cubicn=1.f/(1.f+2.f*iac);
	cc[0]=1.f*cubicn;
	cc[1]=(4.f*iac-2.5f)*cubicn; 
	cc[2]=(1.5f-3.f*iac)*cubicn;		       
	cc[3]=(2.f*iac-0.5f)*cubicn;			
	cc[4]=(0.5f-iac)*cubicn; 
	curmat res(5,1);
	cudaMemcpy(res, cc, 5*sizeof(Real), cudaMemcpyHostToDevice);
	return res;
    }else{
	return curmat();
    }
}

void gpu_dmreal2gpu(mapcell *dmreal){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cp2gpu(cudata->dmreal, dmreal);
    }
}
void gpu_dmproj2gpu(mapcell *dmproj){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cp2gpu(cudata->dmproj, dmproj);
    }
}


/*This is memory bound. So increasing # of points processed does not help. */
__global__ void prop_linear(Real *restrict out, const Real *restrict in,
			    const int nx, const int ny, KARG_COMMON){
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
__global__ void prop_linear(Real *restrict out, const Comp *restrict in,
			    const int nx, const int ny, KARG_COMMON){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real x=loc[i][0]*dxi+dispx;
	Real y=loc[i][1]*dyi+dispy;
	int ix=Z(floor)(x);
	int iy=Z(floor)(y);
	x=x-ix; y=y-iy;
	if(ix>=0 && ix<nx-1 && iy>=0 && iy<ny-1){
	    Real tmp=((+in[iy*nx+ix].x*(1.f-x) +in[iy*nx+ix+1].x*x)*(1.f-y)
		      +(+in[(iy+1)*nx+ix].x*(1.f-x) +in[(iy+1)*nx+ix+1].x*x)*y);
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
__global__ void prop_linear_wrap(Real *restrict out, const Real *restrict in,
				 const int nx, const int ny, KARG_COMMON){
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
__global__ void prop_cubic(Real *restrict out, const Real *restrict in,
			   const int nx, const int ny, KARG_COMMON, const Real *cc){
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
void gpu_atm2loc(Real *phiout, const culoc_t &loc, Real hs, Real hc, Real thetax, Real thetay,
		  Real mispx, Real mispy, Real dt, int isim, Real atmalpha, cudaStream_t stream){
    cumapcell &cuatm=cudata->atm;
    if(Z(fabs)(atmalpha)<EPS) return;
    if(cuglobal->atmscale){
	atmalpha*=cuglobal->atmscale->p[isim];
    }
    for(int ips=0; ips<cudata->nps; ips++){
	const Real dx=cuatm[ips].dx;
	const Real dy=cuatm[ips].dy;
	const Real ht=cuatm[ips].ht-hc;
	const Real vx=cuatm[ips].vx;
	const Real vy=cuatm[ips].vy;
	const Real dispx=(ht*thetax+mispx-vx*dt*isim-cuatm[ips].ox)/dx;
	const Real dispy=(ht*thetay+mispy-vy*dt*isim-cuatm[ips].oy)/dy;
	const Real scale=1.f-ht/hs;
	if(scale<0) continue;
	const int nloc=loc.Nloc();
	
#define COMM loc(),loc.Nloc(),scale/dx,scale/dy, dispx, dispy, atmalpha
	if(WRAP_ATM){
	    prop_linear_wrap<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips](), cuatm[ips].nx, cuatm[ips].ny, COMM);
	}else{/*we are gauranteed. */
	    prop_linear_nocheck<<<DIM(nloc,256), 0, stream>>>
		(phiout, cuatm[ips](), cuatm[ips].nx, cuatm[ips].ny, COMM);
	}
#undef COMM
    }
}
void gpu_map2loc(const cumap_t &map, const culoc_t &loc, Real *phiout,
		 Real alpha, Real dispx, Real dispy, Real scale, int wrap, cudaStream_t stream){
    (void)wrap;
    if(scale<0) return;
    dispx=(dispx-map.ox)/map.dx;
    dispy=(dispy-map.oy)/map.dy;
    const int nloc=loc.Nloc();
    if (map.cubic_cc){//128 is a good number for cubic. 
	prop_cubic<<<DIM(nloc,128), 0, stream>>>
	    (phiout, map(),map.nx,map.ny, loc(), loc.Nloc(),scale/map.dx,scale/map.dy, dispx, dispy, alpha, 
	     map.cubic_cc());
    }else{
	prop_linear<<<DIM(nloc,256), 0, stream>>>
	    (phiout, map(),map.nx,map.ny, loc(), loc.Nloc(),scale/map.dx,scale/map.dy, dispx, dispy, alpha);
    }
}
    
/**
   Ray tracing of dm. use locondm for each dm. so that distortion can be properly accounted for. Use the other version if no distortion.
*/
void gpu_dm2loc(Real *phiout, const Array<culoc_t> &locondm, const cumapcell &cudm, int ndm,
		Real hs,  Real hc, Real thetax, Real thetay, Real mispx, Real mispy, Real alpha, cudaStream_t stream){
    for(int idm=0; idm<ndm; idm++){
	assert(cudm[idm].ny>1);//prevent accidentally pass in a vector
	const Real ht=cudm[idm].ht-hc;
	gpu_map2loc(cudm[idm], locondm[idm], phiout, 
		    alpha, ht*thetax+mispx, ht*thetay+mispy,
		    1-ht/hs, 0, stream);
    }/*idm */
}
/**
   Ray tracing of dm. 
*/
void gpu_dm2loc(Real *phiout, const culoc_t &locout, const cumapcell &cudm, int ndm,
		Real hs, Real hc, Real thetax, Real thetay, Real mispx, Real mispy, Real alpha, cudaStream_t stream){
    for(int idm=0; idm<ndm; idm++){
	assert(cudm[idm].ny>1);//prevent accidentally pass in a vector
	const Real ht=cudm[idm].ht-hc;
	gpu_map2loc(cudm[idm], locout, phiout, 
		    alpha, ht*thetax+mispx, ht*thetay+mispy,
		    1-ht/hs, 0, stream);
    }/*idm */
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void gpu_ngsmod2science(curmat &opd, Real (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream){
    if(ngsmod->nmod==2){
	curaddptt(opd, loc, 0, mod[0]*alpha, mod[1]*alpha, stream);
    }else{
	const Real ht=ngsmod->ht;
	const Real scale=ngsmod->scale;

	Real focus=0, ps1=0, ps2=0, ps3=0,astigx=0,astigy=0;
	if(ngsmod->indfocus){
	    focus+=mod[ngsmod->indfocus];
	}
	if(ngsmod->indps){
	    if(!ngsmod->ahstfocus){
		focus+=mod[ngsmod->indps]*(1.f-scale);
	    }
	    ps1=mod[ngsmod->indps];
	    ps2=mod[ngsmod->indps+1];
	    ps3=mod[ngsmod->indps+2];
	}
	if(ngsmod->indastig){
	    astigx=mod[ngsmod->indastig];
	    astigy=mod[ngsmod->indastig+1];
	}

	add_ngsmod_do<<<DIM(opd.N(), 256), 0, stream>>>
	    (opd(), loc, opd.N(),
	     mod[0], mod[1], ps1, ps2, ps3, astigx, astigy, focus,
	     thetax, thetay, scale, ht, alpha);
    }
}
