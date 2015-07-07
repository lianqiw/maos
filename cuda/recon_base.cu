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
#include "utils.h"
#include "recon_base.h"
#include "curmat.h"
#include "accphi.h"
namespace cuda_recon{


W01_T::W01_T(const dsp *R_W0, const dmat *R_W1, int R_nxx)
    :nxx(R_nxx),pis(0){
    if(!R_W0 || !R_W1){
	error("R0, R1 must not be empty\n");
    }
    cp2gpu(W1, R_W1);
    {
	/*W0 of partially illuminates subaps are stored as sparse matrix in W0f in
	  GPU. W0 of fully illuminated subaps are stored in W0p.*/
	spint *pp=R_W0->p;
	spint *pi=R_W0->i;
	double *px=R_W0->x;
	dsp *W0new=dspnew(R_W0->m, R_W0->n, R_W0->nzmax);
	spint *pp2=W0new->p;
	spint *pi2=W0new->i;
	double *px2=W0new->x;
	int *full;
	DO(cudaMallocHost(&full, R_W0->n*sizeof(int)));
	//#define W0_BW 1
	double W1max=dmax(R_W1);
	double thres=W1max*(1.f-1e-6);
	W0v=(Real)(W1max*4./9.);//max of W0 is 4/9 of max of W1. 
	int count=0;
	int count2=0;
	for(int ic=0; ic<R_W0->n; ic++){
	    pp2[ic]=count;
	    if(R_W1->p[ic]>thres){
		full[count2]=ic;
		count2++;
	    }else{
		int nv=pp[ic+1]-pp[ic];
		memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
		memcpy(px2+count, px+pp[ic], sizeof(double)*nv);
		count+=nv;
	    }
	}
	pp2[R_W0->n]=count;
	W0new->nzmax=count;
	W0p=cusp(W0new, 1);
	cp2gpu(W0f, full, count2, 1);
	dspfree(W0new);
	cudaFreeHost(full);
    }
}
/**
   res_vec[i]+=sum_j(as[i][j]*b[j]);
*/
__global__ void 
inn_multi_do(Real *res_vec, const Real *as, const Real *b, const int n){
    extern __shared__ Real sb[];
    sb[threadIdx.x]=0;
    const Real *a=as+n*blockIdx.y;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=blockDim.x * gridDim.x){
	sb[threadIdx.x]+=a[i]*b[i];
    }
    for(int step=(blockDim.x>>1);step>0;step>>=1){
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
assign_multi_do(Real *as, const Real *b, Real *beta1, Real beta2, int n){
    Real *a=as+blockIdx.y*n;
    Real beta=beta1[blockIdx.y]*beta2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=b[i]*beta;
    }
}
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ void 
apply_W0_do(Real *outs, const Real *ins, const int *W0f, Real W0v, int nx, int ntot, int nW0f){
    Real *out=outs+blockIdx.y*ntot;
    const Real *in=ins+blockIdx.y*ntot;
    const int step=blockDim.x * gridDim.x;
    for(int k=blockIdx.x * blockDim.x + threadIdx.x; k<nW0f; k+=step){
	int i=W0f[k];
	out[i]+=W0v*(in[i]
		       +0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
		       +0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
    }
}
void W01_T::apply(Real *restrict xout, const Real *xin, int ndir, stream_t &stream){
    if(pis.Nx() < ndir){
	pis=curmat(ndir,1);
    }else{
	pis.zero(stream);
    }
    //Apply W1: Piston removal
    inn_multi_do<<<dim3(32, ndir), dim3(DIM_REDUCE), DIM_REDUCE*sizeof(Real), stream>>>
	(pis.P(), xin, W1.P(), W1.Nx());
    assign_multi_do<<<dim3(32, ndir), dim3(256), 0, stream>>>
	(xout, W1.P(), pis.P(), -1, W1.Nx());
    //Apply W0: bilinar-weighting
    if(W0f){
	apply_W0_do<<<dim3(16, ndir), dim3(256,1), 0, stream>>> 		
	    (xout, xin, W0f.P(), W0v, nxx, W1.Nx(), W0f.Nx());
    }
    if(W0p){
	cuspmul(xout, W0p, xin, ndir, 'n', 1.f, stream);
    }
}

curecon_geom::curecon_geom(const PARMS_T *parms, const RECON_T *recon)
    :npsr(0),ndm(0),delay(0),isim(0),reconisim(0),
     xnx(0),xny(0),anx(0),any(0), anloc(0),ngrad(0), dt(0){
    if(!parms) return;
    ndm=parms->ndm;
    npsr=parms->sim.idealfit?parms->atm.nps:recon->npsr;
    pmap=(recon->pmap);
    fmap=(recon->fmap);
    /*Setup various grid*/
    amap=cugridcell(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	amap[idm]=(recon->amap->p[idm]);
    }
    if(recon->xloc){
	xmap=cugridcell(npsr,1);
	for(int ipsr=0; ipsr<npsr; ipsr++){
	    xmap[ipsr]=recon->xmap->p[ipsr];
	}
	xnx=recon->xnx->p;
	xny=recon->xny->p;
    }
    if(recon->xcmap){
	xcmap=cugridcell(npsr,1);
	for(int ipsr=0; ipsr<npsr; ipsr++){
	    xcmap[ipsr]=(recon->xcmap->p[ipsr]);
	}
    }
    W01=W01_T(recon->W0, recon->W1, fmap.nx);
    anx=recon->anx->p;
    any=recon->any->p;
    anloc=recon->anloc->p;
    ngrad=recon->ngrad->p;
    dt=parms->sim.dt;
    delay=2;//2 frame delay
}

map_l2d::map_l2d(const cugrid_t &out, dir_t *dir, int _ndir, //output.
		 const cugridcell &in, int _nlayer,//input. layers.
		 Real dt){//directions and star height.
    nlayer=_nlayer;
    ndir=_ndir;
    PROP_WRAP_T *hdata_cpu=new PROP_WRAP_T[nlayer*ndir];
    DO(cudaMalloc(&hdata, sizeof(PROP_WRAP_T)*nlayer*ndir));
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	const Real ht=in[ilayer].ht;
	for(int idir=0; idir<ndir; idir++){
	    if(!dir[idir].skip){
		const Real dispx=dir[idir].thetax*ht+in[ilayer].vx*dt;
		const Real dispy=dir[idir].thetay*ht+in[ilayer].vy*dt;
		const Real hs=dir[idir].hs;
		const Real scale=1.f-ht/hs;
		cugrid_t outscale=out.Scale(scale);
		gpu_map2map_prep(hdata_cpu+idir+ilayer*ndir, outscale, in[ilayer],
				   dispx, dispy, in[ilayer].cubic_cc);
	    }
	    hdata_cpu[idir+ilayer*ndir].togpu(hdata+idir+ilayer*ndir);
	}
    }
    delete [] hdata_cpu;
}

map_l2l::map_l2l(const cugridcell &out, const cugridcell &in, int _nlayer){//input. layers.
    nlayer=_nlayer;
    ndir=0;//this is laye to layer.
    PROP_WRAP_T *hdata_cpu=new PROP_WRAP_T[nlayer];
    DO(cudaMalloc(&hdata, sizeof(PROP_WRAP_T)*nlayer));
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	if(Z(fabs)(out[ilayer].ht-in[ilayer].ht)>EPS){
	    error("Layer height miX(mat)ch\n");
	}
	gpu_map2map_prep(hdata_cpu+ilayer, out[ilayer], in[ilayer],
			   0, 0, in[ilayer].cubic_cc);
	hdata_cpu[ilayer].togpu(hdata+ilayer);
    }
    delete [] hdata_cpu;
}

}//namespace
