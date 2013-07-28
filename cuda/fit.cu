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
#include "utils.h"
#include "accphi.h"
#include "fit.h"
#define TIMING 0
#if TIMING <1
#undef EVENT_INIT
#undef EVENT_TIC
#undef EVENT_TOC
#undef EVENT_PRINT
#define EVENT_INIT(A)
#define EVENT_TIC(A)
#define EVENT_TOC
#define EVENT_PRINT(A...)
#else
#define EVENT_PRINT(A...) fprintf(stderr, A);EVENT_DEINIT
#endif

#if TIMING <2
#define EVENT2_INIT(A)
#define EVENT2_TIC(A)
#define EVENT2_TOC
#define EVENT2_PRINT(A...)
#else
#define EVENT2_INIT EVENT_INIT
#define EVENT2_TIC EVENT_TIC
#define EVENT2_TOC EVENT_TOC
#define EVENT2_PRINT(A...) fprintf(stderr, A);EVENT_DEINIT
#endif

void cufit_grid::init(const PARMS_T *parms, const RECON_T *recon){
    /*Initialize*/
    const int ndm=parms->ndm;
    const int npsr=recon->npsr;
    nfit=parms->fit.nfit;
    /*Setup various grid*/
 
    if(parms->sim.idealfit){
	floc=new culoc_t(recon->floc);
	fit_thetax=new float[nfit];
	fit_thetay=new float[nfit];
	for(int ifit=0; ifit<nfit; ifit++){
	    fit_thetax[ifit]=parms->fit.thetax[ifit];
	    fit_thetay[ifit]=parms->fit.thetay[ifit];
	}
    }
    /*Various*/
    if(recon->fitNW){
	dmat *_fitNW=dcell2m(recon->fitNW);
	cp2gpu(&fitNW, _fitNW);
	dfree(_fitNW);
	dotNW=curnew(fitNW->ny, 1);
    }
    if(recon->actslave){
	cp2gpu(&actslave, recon->actslave, 1);
    }
    if(parms->fit.cachedm){
	long acnx[ndm], acny[ndm];
	for(int idm=0; idm<ndm; idm++){
	    acnx[idm]=grid->acmap[idm].nx;
	    acny[idm]=grid->acmap[idm].ny;
	}
	dmcache=curcellnew(ndm, 1, acnx, acny);
    }
    if(parms->fit.cachex){
	long xcnx[npsr], xcny[npsr];
	for(int ips=0; ips<npsr; ips++){
	    xcnx[ips]=grid->xcmap[ips].nx;
	    xcny[ips]=grid->xcmap[ips].ny;
	}
	xcache=curcellnew(npsr, 1, xcnx, xcny);
    } 
    cp2gpu(&fitwt, recon->fitwt);
 
    long fnx[nfit],fny[nfit];
    for(int ifit=0; ifit<nfit; ifit++){
	fnx[ifit]=grid->fmap.nx;
	fny[ifit]=grid->fmap.ny;
    }
    opdfit=curcellnew(nfit, 1, fnx, fny);
    opdfit2=curcellnew(nfit, 1, fnx, fny);
    opdfitv=curnew(grid->fmap.nx*grid->fmap.ny, nfit, opdfit->m->p, 0);
    opdfit2v=curnew(grid->fmap.nx*grid->fmap.ny, nfit, opdfit2->m->p, 0);
    pis=curnew(1, parms->fit.nfit);

    /*Data for ray tracing*/
    //dm -> floc
    PROP_WRAP_T *hadata=new PROP_WRAP_T[ndm*nfit];
    DO(cudaMalloc(&this->hadata, sizeof(PROP_WRAP_T)*ndm*nfit));
    //dm: amap->acmap
    PROP_WRAP_T *ha0data=NULL, *ha1data=NULL;
    if(parms->fit.cachedm){
	ha0data=new PROP_WRAP_T[ndm];
	ha1data=new PROP_WRAP_T[ndm*nfit];
	DO(cudaMalloc(&this->ha0data, sizeof(PROP_WRAP_T)*ndm));
	DO(cudaMalloc(&this->ha1data, sizeof(PROP_WRAP_T)*ndm*nfit));
    }
    if(!parms->sim.idealfit){
    //xmap -> fmap
	PROP_WRAP_T *hxpdata=new PROP_WRAP_T[npsr*nfit];
	PROP_WRAP_T *hxp0data=NULL, *hxp1data=NULL;
	DO(cudaMalloc(&this->hxpdata, sizeof(PROP_WRAP_T)*npsr*nfit));
	if(parms->fit.cachex){
	    hxp0data=new PROP_WRAP_T[npsr];
	    hxp1data=new PROP_WRAP_T[npsr*nfit];
	    DO(cudaMalloc(&this->hxp0data, sizeof(PROP_WRAP_T)*npsr));
	    DO(cudaMalloc(&this->hxp1data, sizeof(PROP_WRAP_T)*npsr*nfit));
	}
    
	for(int ipsr=0; ipsr<npsr; ipsr++){
	    const float ht=recon->ht->p[ipsr];
	    if(parms->fit.cachex){
		gpu_prop_grid_prep(hxp0data+ipsr, grid->xcmap[ipsr], grid->xmap[ipsr],
				   0,0, NULL);
		hxp0data[ipsr].togpu(&this->hxp0data[ipsr]);
	    }
	    for(int ifit=0;  ifit<nfit; ifit++){
		const float hs=parms->fit.hs[ifit];
		const float thetax=(float)parms->fit.thetax[ifit];
		const float thetay=(float)parms->fit.thetay[ifit];
		const float scale=1.f-ht/hs;
		cugrid_t fmapscale=grid->fmap*scale;
		gpu_prop_grid_prep(hxpdata+ifit+ipsr*nfit, fmapscale, grid->xmap[ipsr],
				   thetax*ht, thetay*ht, NULL);
		hxpdata[ifit+ipsr*nfit].togpu(&this->hxpdata[ifit+ipsr*nfit]);
		if(parms->fit.cachex){
		    gpu_prop_grid_prep(hxp1data+ifit+ipsr*nfit, fmapscale, grid->xcmap[ipsr],
				       thetax*ht, thetay*ht, NULL);
		    hxp1data[ifit+ipsr*nfit].togpu(&this->hxp1data[ifit+ipsr*nfit]);
		}
	    }
	}
	delete [] hxpdata;
	delete [] hxp0data;
	delete [] hxp1data;
    }
    for(int idm=0; idm<ndm; idm++){
	const float ht=parms->dm[idm].ht;
	if(parms->fit.cachedm){
	    gpu_prop_grid_prep(ha0data+idm, grid->acmap[idm], grid->amap[idm],
			       0, 0, grid->cubic_cc->p[idm]);
	    ha0data[idm].togpu(&this->ha0data[idm]);
	}
	for(int ifit=0; ifit<nfit; ifit++){
	    const float hs=parms->fit.hs[ifit];
	    const float thetax=(float)parms->fit.thetax[ifit];
	    const float thetay=(float)parms->fit.thetay[ifit];
	    const float scale=1.f-ht/hs;
	    cugrid_t fmapscale=grid->fmap*scale;
	    if(parms->fit.cachedm){
		gpu_prop_grid_prep(ha1data+ifit+idm*nfit, fmapscale, grid->acmap[idm],
				   thetax*ht, thetay*ht, NULL);
		ha1data[ifit+idm*nfit].togpu(&this->ha1data[ifit+idm*nfit]);
	    }
	    gpu_prop_grid_prep(hadata+ifit+idm*nfit, fmapscale, grid->amap[idm],
			       thetax*ht, thetay*ht, grid->cubic_cc->p[idm]);	
	    hadata[ifit+idm*nfit].togpu(&this->hadata[ifit+idm*nfit]);
	}
    }
    delete [] hadata;
    delete [] ha0data;
    delete [] ha1data;
}

/*
  Todo: share the ground layer which is both matched and same.
*/
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ void 
apply_W_do(float *outs, const float *ins, const int *W0f, float alpha, int nx, int n){
    float *out=outs+blockIdx.y*nx*nx;
    const float *in=ins+blockIdx.y*nx*nx;
    const int step=blockDim.x * gridDim.x;
    for(int k=blockIdx.x * blockDim.x + threadIdx.x; k<n; k+=step){
	int i=W0f[k];
	out[i]+=alpha*(in[i]
		       +0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
		       +0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
    }
}
/**
do HXp operation, opdfit+=Hxp*xin*/
void cufit_grid::do_hxp(const curcell *xin, stream_t &stream){
    EVENT2_INIT(4);
    EVENT2_TIC(0);
    opdfit->m->zero(stream);
    EVENT2_TIC(1);
    if(!hxpdata){//idealfit
	for(int ifit=0; ifit<nfit; ifit++){
	    gpu_atm2loc(opdfit->p[ifit]->p, floc->p, floc->nloc, INFINITY,
			fit_thetax[ifit], fit_thetay[ifit], 
			0, 0, grid->dt*grid->isim, 1, stream);
	}
    }else{
	const int npsr=xin->nx;
	if(xcache){
	    xcache->m->zero(stream);
	    gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
		(hxp0data, xcache->pm, xin->pm, 1, npsr, 1, NULL, 'n');
	    EVENT2_TIC(2);
	    gpu_prop_grid_do<<<dim3(3,3,nfit),dim3(16,16),0,stream>>>
		(hxp1data, opdfit->pm, xcache->pm, nfit, npsr, 1, NULL, 'n');
	    EVENT2_TIC(3);
	}else{
	    EVENT2_TIC(2);
	    gpu_prop_grid_do<<<dim3(3,3,nfit),dim3(16,16),0,stream>>>
		(hxpdata, opdfit->pm, xin->pm, nfit, npsr, 1, NULL, 'n');
	    EVENT2_TIC(3);
	}
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_hxp: zero: %.3f cache: %.3f prop: %.3f tot: %.3f\n", times[1], times[2], times[3], times[0]);
}
/**
do HXp' operation, xout+=alpha*Hxp'*opdfit2*/
void cufit_grid::do_hxpt(const curcell *xout, float alpha, stream_t &stream){
    const int npsr=xout->nx;
    EVENT2_INIT(4);
    EVENT2_TIC(0);
    if(xcache){
	xcache->m->zero(stream);
	gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
	    (hxp1data, opdfit2->pm, xcache->pm, nfit, npsr, alpha, fitwt->p, 't');
	EVENT2_TIC(2);
	gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
	    (hxp0data, xcache->pm, xout->pm, 1, npsr, alpha, NULL, 't');
	EVENT2_TIC(3);
    }else{
	gpu_prop_grid_do<<<dim3(3,3,npsr), dim3(16,16),0,stream>>>
	    (hxpdata, opdfit2->pm, xout->pm, nfit, npsr, alpha, fitwt->p, 't');
	EVENT2_TIC(3);
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_hxpt: zero: %.3f prop: %.3f cache: %.3f tot: %.3f\n", times[1], times[2], times[3], times[0]);
}
/**
   res_vec[i]+=sum_j(as[i][j]*b[j]);
 */
__global__ void 
inn_multi_do(float *res_vec, const float *as, const float *b, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    const float *a=as+n*blockIdx.y;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i]*b[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(&res_vec[blockIdx.y], sb[0]);
    }
}
/**
   a[i][j]=b[j]*beta1[i]*beta2;
 */
__global__ void 
assign_multi_do(float *as, const float *b, float *beta1, float beta2, int n){
    float *a=as+blockIdx.y*n;
    float beta=beta1[blockIdx.y]*beta2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=b[i]*beta;
    }
}
/**
   opdfit2 = (W0-W1*W1')*opdfit */
void cufit_grid::do_w(stream_t &stream){
    EVENT2_INIT(6);
    EVENT2_TIC(0);
    EVENT2_TIC(1);
    pis->zero(stream);
    const int nxf=opdfit->p[0]->nx;
    const int nf=opdfit->p[0]->nx*opdfit->p[0]->ny;
    const W01_T *W01=grid->W01;
    //inn_multi_do takes only 28 us while curmv takes 121 us.
    inn_multi_do<<<dim3(32, nfit), dim3(DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>
	(pis->p, opdfit->m->p, W01->W1->p, nf);
    EVENT2_TIC(2);
    //assign_multi_do takes 8 us while curmm takes 28.6 us
    assign_multi_do<<<dim3(32, nfit), dim3(256), 0, stream>>>
	(opdfit2->m->p, W01->W1->p, pis->p, -1, nf);
    EVENT2_TIC(3);
    if(W01->nW0f){ //19us
	apply_W_do<<<dim3(16, nfit), dim3(256,1), 0, stream>>> 		
	    (opdfit2->m->p, opdfit->m->p, W01->W0f, 
	     W01->W0v, nxf, W01->nW0f); 
    }
    EVENT2_TIC(4);
    if(W01->W0p){ //53 us *********Can we use B/W map to avoid this matrix?***********
	cuspmul(opdfit2->m->p, W01->W0p, 
		opdfit->m->p, nfit, 'n', 1.f, stream); //80us
    }
    EVENT2_TIC(5);
    EVENT2_TOC;
    EVENT2_PRINT("do_w: zero: %.3f W1_dot: %.3f W1_add: %.3f W0f: %.3f W0p: %.3f tot: %.3f\n",
		times[1], times[2], times[3], times[4], times[5], times[0]);
}

/**
   opdfit+=Ha*xin;
*/
void cufit_grid::do_ha(const curcell *xin, stream_t &stream){
    opdfit->m->zero(stream); 
    const int ndm=xin->nx;
    EVENT2_INIT(3);
    EVENT2_TIC(0);
    if(dmcache){
	/*xout->dmcache*/ 
	dmcache->m->zero(stream); 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (ha0data, dmcache->pm, xin->pm, 
	     1, ndm, 1.f,NULL,'n');//25 us
	EVENT2_TIC(1);
	/*dmcache->opdfit*/ 
	gpu_prop_grid_do<<<dim3(8,8,nfit),dim3(16,16),0,stream>>> 
	    (ha1data, opdfit->pm, dmcache->pm, 
	     nfit, ndm, 1.f,NULL, 'n'); // 87 us
	EVENT2_TIC(2);
    }else{ 
	EVENT2_TIC(1);
	/*xout->opfit*/ 
	gpu_prop_grid_do<<<dim3(8,8,nfit),dim3(16,16),0,stream>>> 
	    (hadata, opdfit->pm, xin->pm, 
	     nfit, ndm, 1.f,NULL, 'n'); 
	EVENT2_TIC(2);
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_ha: cache: %.3f prop: %.3f tot: %.3f\n", times[1], times[2], times[0]);
}

/**
   xout+=alpha*HA'*opdfit2*/
void cufit_grid::do_hat(curcell *xout,  float alpha, stream_t &stream){
    const int ndm=xout->nx;
    EVENT2_INIT(3);
    EVENT2_TIC(0);
    if(dmcache){ 
	/*opdfit2->dmcache*/ 
	dmcache->m->zero(stream); 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (ha1data, opdfit2->pm, dmcache->pm, 
	     nfit, ndm, 1.,fitwt->p, 't'); //105 us
	EVENT2_TIC(1);
	/*dmcache->xout*/ 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (ha0data, dmcache->pm, xout->pm, 
	     1, ndm, alpha, NULL,'t'); //278 us. ******NEED Improvement****
	EVENT2_TIC(2);
    }else{ 
	/*opfit2->xout	*/ 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (hadata, opdfit2->pm, xout->pm, 
	     nfit, ndm, alpha, fitwt->p, 't'); 
	EVENT2_TIC(1);
	EVENT2_TIC(2);
    } 
    EVENT2_TOC;
    EVENT2_PRINT("do_hat: prop: %.3f cache: %.3f tot: %.3f\n", times[1], times[2], times[0]);
}

/*
  Right hand size operator. 
*/
void cufit_grid::R(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
    if(!*xout){
	*xout=curcellnew(grid->ndm, 1, grid->anx, grid->any);
    }else{
	curscale((*xout)->m, beta, stream);
    }
    do_hxp(xin, stream);//153 us
    do_w(stream);//123 us
    do_hat(*xout, alpha, stream);//390 us
}
void cufit_grid::Rt(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
    if(!*xout){
	*xout=curcellnew(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curscale((*xout)->m, beta, stream);
    }
    do_ha(xin, stream);
    do_w(stream);
    do_hxpt(*xout, alpha, stream);
}
void cufit_grid::L(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
    const int ndm=grid->ndm;
    EVENT_INIT(6);
    EVENT_TIC(0);
    if(!*xout){
	*xout=curcellnew(ndm, 1, grid->anx, grid->any);
    }else{
	curscale((*xout)->m, beta, stream);
    }   
    do_ha(xin, stream);//112 us
    EVENT_TIC(1);
    do_w(stream);//108 us
    EVENT_TIC(2);
    do_hat(*xout, alpha, stream);//390 us
    EVENT_TIC(3);
    if(fitNW){
	curmv(dotNW->p, 0, fitNW, xin->m->p, 't', 1, stream);
	curmv((*xout)->m->p, 1, fitNW, dotNW->p, 'n', alpha, stream);
    }
    EVENT_TIC(4);
    if(actslave){
	cuspmul((*xout)->m->p, actslave, xin->m->p, 1,'n', alpha, stream);
    }
    EVENT_TIC(5);
    EVENT_TOC;
    EVENT_PRINT("FitL HA: %.3f W: %.3f HAT: %.3f NW: %.3f SL: %.3f tot: %.3f\n",
	  times[1],times[2],times[3],times[4],times[5],times[0]);
}
