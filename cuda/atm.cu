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
#include "cudata.h"

/*
  2019-09-03

  Separated from accphi.cu.

  Include routines that copy atm and dm from CPU to GPU memory
*/



static void atm_prep(atm_prep_t *data){
    PNEW(lock);
    LOCK(lock);/*make sure we only read one layer at a time. */
    //TIC;tic;
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
    //toc2("Step %d: Layer %d: Preparing atm for step %d", data->isim, ips, data->isim_next);
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
		{
		    nx0=MIN(parms->atm.nx,MIN(nxa, nxn*2));
		    ny0=MIN(parms->atm.ny,MIN(nxa, nxn*2));
		}
	    }
	    info("We will host %dx%d in GPU, taking %zd MiB\n", 
		 nx0, ny0, (nx0*ny0*nps*sizeof(Real))>>20);
	}
    }
    /*The atm in GPU is the same as in CPU. */
    if(nx0==parms->atm.nx && ny0==parms->atm.ny){
	cuglobal->atm_full=1;
	if(iseed0!=iseed){
	    gpu_atm2gpu_full(atmc);
	    iseed0=iseed;
	}
	return;
    }
    cuglobal->atm_full=0;
    const real dt=parms->sim.dt;
    const real dx=parms->atm.dx;

    if(iseed0==-1){//Initializing data
	/*The atm in GPU is smaller than in CPU. */
	cuglobal->atm_prep_data.init(nps);
	for(int ips=0; ips<nps; ips++){
	    cuglobal->atm_prep_data[ips].nx0=nx0;
	    cuglobal->atm_prep_data[ips].ny0=ny0;
	}
	for(int im=0; im<NGPU; im++){/*Loop over all GPUs. */
	    gpu_set(im);
	    cudata->atm=cumapcell(nps, 1);
	    for(int ips=0; ips<nps; ips++){
		cudata->atm[ips].p=curmat(nx0, ny0);
		cudata->atm[ips].nx=nx0;
		cudata->atm[ips].ny=ny0;
	    }
	}/*for im */
    }/*if need_init; */
    Array<atm_prep_t>& prep_data=cuglobal->atm_prep_data;
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
	/*Load atmosphere to Real memory in advance. This takes time if atm is stored in file.*/
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
	    
	    /*launch an independent thread to pull data in. thread will be joined when needed.*/
	    pthread_create(&prep_data[ips].threads, NULL, (void *(*)(void *))atm_prep, prep_data+ips);
	}/*need to cpy atm to next_atm. */
    }
    for(int ips=0; ips<nps; ips++) {
	if(isim==prep_data[ips].isim_next){
	    /*need to copy atm to gpu. and update next_isim */
	    //TIC;tic;
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
	    if(prep_data[ips].isim_next>=parms->sim.end){/*no need to do */
		prep_data[ips].isim_next=INT_MAX;
	    }
	 
	    //toc2("Step %d: Layer %d transfered. next in step %d. ",isim, ips, prep_data[ips].isim_next); 
	}//if isim
    }//for ips
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
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
