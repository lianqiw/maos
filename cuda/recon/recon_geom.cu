/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "../sim/accphi.h"
#include "recon_geom.h"



w01_t::w01_t(const dsp* R_W0, const dmat* R_W1, int R_nxx){
	nxx=R_nxx;
	if(!R_W0||!R_W1){
		error("R0, R1 must not be empty\n");
	}
	cp2gpu(W1, R_W1);
	if(1){
	/*W0 of partially illuminates points are stored as sparse matrix in W0f in
	  GPU. W0 of fully illuminated points are stored in W0p.*/
		spint* pp=R_W0->pp;
		spint* pi=R_W0->pi;
		real* px=R_W0->px;
		dsp* W0new=dspnew(R_W0->nx, R_W0->ny, R_W0->nzmax);
		spint* pp2=W0new->pp;
		spint* pi2=W0new->pi;
		real* px2=W0new->px;
		//int *full;
		Array<int> full(R_W0->ny, 1);
		//#define W0_BW 1
		real W1max=dmax(R_W1);
		real thres=W1max*(1.f-1e-6);
		W0v=(Real)(W1max*4./9.);//max of W0 is 4/9 of max of W1. 
		int count=0;
		int count2=0;
		for(int ic=0; ic<R_W0->ny; ic++){
			pp2[ic]=count;
			if(R_W1->p[ic]>thres){
				full[count2]=ic;
				count2++;
			} else{
				int nv=pp[ic+1]-pp[ic];
				memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
				memcpy(px2+count, px+pp[ic], sizeof(real)*nv);
				count+=nv;
			}
		}
		pp2[R_W0->ny]=count;
		W0new->nzmax=count;
		//W0new is the transpose of W0p.
		dsp* W0new2=dsptrans(W0new); dspfree(W0new);
		W0p=cusp(W0new2, 1);
		cp2gpu(W0f, full(), count2, 1);
		dspfree(W0new2);
		//cudaFreeHost(full);
	} else{
		W0p=cusp(R_W0, 1);
	}
}
/**
   res_vec[i]+=sum_j(as[i][j]*b[j]);
*/
__global__ void
inn_multi_do(Real* res_vec, const Real* as, const Real* b, const int n){
	extern __shared__ Real sb[];
	sb[threadIdx.x]=0;
	const Real* a=as+n*blockIdx.y;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
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
assign_multi_do(Real* as, const Real* b, Real* beta1, Real beta2, int n){
	Real* a=as+blockIdx.y*n;
	Real beta=beta1[blockIdx.y]*beta2;
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]=b[i]*beta;
	}
}
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ void
apply_W0_do(Real* outs, const Real* ins, const int* W0f, Real W0v, int nx, int ntot, int nW0f){
	Real* out=outs+blockIdx.y*ntot;
	const Real* in=ins+blockIdx.y*ntot;
	const int step=blockDim.x*gridDim.x;
	for(int k=blockIdx.x*blockDim.x+threadIdx.x; k<nW0f; k+=step){
		int i=W0f[k];
		out[i]+=W0v*(in[i]
			+0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
			+0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
	}
}
/**
 * xout=w0*xin - w1*(w1'*xin).
 * Notice that the weighting per direction is applied at the next step (HA^T)
 * */
void w01_t::apply(curcell& xout, const curcell& xin, stream_t& stream) const{
	int ndir=xin.Nx();
	if(pis.Nx()<ndir){
		pis=curmat(ndir, 1);
	} else{
		pis.Zero(stream);
	}
	//Apply W1: Piston removal
	inn_multi_do<<<dim3(32, ndir), dim3(DIM_REDUCE), DIM_REDUCE*sizeof(Real), stream>>>
		(pis(), xin.M()(), W1(), W1.Nx());
	assign_multi_do<<<dim3(32, ndir), dim3(256), 0, stream>>>
		(xout.M()(), W1(), pis(), -1, W1.Nx());
	//Apply W0: bilinar-weighting
	if(W0f){
		apply_W0_do<<<dim3(16, ndir), dim3(256, 1), 0, stream>>>
			(xout.M()(), xin.M()(), W0f(), W0v, nxx, W1.Nx(), W0f.Nx());
	}
	if(W0p){
		cuspmul(xout.M()(), W0p, xin.M()(), ndir, 'n', 1.f, stream);
	}
}

curecon_geom::curecon_geom(const parms_t* parms, const recon_t* recon)
	:npsr(0), ndm(0), delay(0), reconisim(0),
	W01(recon->W0, recon->W1, recon->fmap->nx), xnx(0), xny(0), anx(0), any(0), anloc(0), ngrad(0), dt(0){
	dbg("update reconstructor geometry.\n");
	ndm=parms->ndm;
	npsr=parms->sim.idealtomo?parms->atm.nps:recon->npsr;
	pmap=recon->pmap;
	fmap=recon->fmap;
	/*Setup various grid*/
	amap=cugridcell(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		amap[idm]=(recon->amap->p[idm]);
	}
	if(recon->xloc){
		xmap=cugridcell(npsr, 1);
		for(int ipsr=0; ipsr<npsr; ipsr++){
			xmap[ipsr]=recon->xmap->p[ipsr];
		}
		xnx=P(recon->xnx);
		xny=P(recon->xny);
	}
	if(recon->xcmap){
		xcmap=cugridcell(npsr, 1);
		for(int ipsr=0; ipsr<npsr; ipsr++){
			xcmap[ipsr]=(recon->xcmap->p[ipsr]);
		}
	}
	anx=P(recon->anx);
	any=P(recon->any);
	anloc=P(recon->anloc);
	ngrad=P(recon->ngrad);
	dt=parms->sim.dt;
	delay=2;//2 frame delay
}


