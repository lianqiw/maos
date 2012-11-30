/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "accphi.h"
#include "cucmat.h"
#include "pcg.h"

#define TIMING 3

/*
  If merge the operation in to gpu_prop_grid_adaptive_do, need to do atomic
operation because previous retracing starts at offo, not 0.  */
__global__ void gpu_laplacian_do(GPU_PROP_GRID_T *datai, float **inall, float **in2all, int nwfs, float alpha){
    float *restrict in, *restrict out;
    {
	const int ips=blockIdx.z;
	in=in2all[ips];
	out=inall[ips];
	datai+=nwfs*ips;
	//GPU_PROP_GRID_T *datai=data+nwfs*ips;
    }
    const int nx=datai->nxi;
    const int ny=datai->nyi;
    const float alpha2=datai->l2c*alpha;
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;
    const int nx1=nx-1;
    const int ny1=ny-1;
    const int ix0=blockIdx.x*blockDim.x+threadIdx.x;
    const int iy0=blockIdx.y*blockDim.y+threadIdx.y;
    for(int iy=iy0; iy<ny; iy+=stepy){
	int iy1 =iy+1; if(iy1>ny1) iy1-=ny;
	int iy2 =iy+2; if(iy2>ny1) iy2-=ny;
	int iy1p=iy-1; if(iy1p<0) iy1p+=ny;
	int iy2p=iy-2; if(iy2p<0) iy2p+=ny;
	for(int ix=ix0; ix<nx; ix+=stepx){
	    int ix1 =ix+1; if(ix1>nx1) ix1-=nx;
	    int ix2 =ix+2; if(ix2>nx1) ix2-=nx;
	    int ix1p=ix-1; if(ix1p<0) ix1p+=nx;
	    int ix2p=ix-2; if(ix2p<0) ix2p+=nx;
	    out[ix+iy*nx]+=alpha2*(20.f*in[ix+iy*nx]
				   -8.f*(in[ix1p+iy*nx]+in[ix1+iy*nx]+in[ix+iy1p*nx]+in[ix+iy1*nx])
				   +2.f*(in[ix1p+iy1p*nx]+in[ix1+iy1p*nx]+in[ix1p+iy1*nx]+in[ix1+iy1*nx])
				   +(in[ix+iy2p*nx]+in[ix2p+iy*nx]+in[ix2+iy*nx]+in[ix+iy2*nx]));
	}
    }
    if(datai->zzi>-1){/*piston constaint*/
	if(threadIdx.x==0 && threadIdx.y==0 && blockIdx.x==0 && blockIdx.y==0){
	    out[datai->zzi]+=in[datai->zzi]*datai->zzv*alpha;
	}
    }
}
__global__ void gpu_prop_grid_adaptive_do(GPU_PROP_GRID_T *data, float **outall, float **inall, int nwfs, int nps, 
					  float alpha, char trans){
    /*
      Each block handles a specific part of some wfs
    */
    int nn;
    if(trans=='t'){
	nn=nwfs;
    }else{
	nn=nps;
    }
    const int ix0=blockIdx.x*blockDim.x+threadIdx.x;
    const int iy0=blockIdx.y*blockDim.y+threadIdx.y;
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;

    for(int ii=0; ii<nn; ii++){
	float *restrict out;
	float *restrict in;

	GPU_PROP_GRID_T *datai;
	{
	    int iwfs, ips;
	    if(trans=='t'){
		ips=blockIdx.z;
		iwfs=ii;
	    }else{
		iwfs=blockIdx.z;
		ips=ii;
	    }
	    datai=data+iwfs+nwfs*ips;
	    if(datai->trans=='r'){/*reverse input/output*/
		out=inall[ips]+datai->offo;
		in=outall[iwfs]+datai->offi;
	    }else{
		out=outall[iwfs]+datai->offo;
		in=inall[ips]+datai->offi;
	    }
	}
	const int nx=datai->nx;
	const int ny=datai->ny;
	if(nx==0) continue;//skip empty wfs
	const int nxo=datai->nxo;
	const int nxi=datai->nxi;
	
	
	if(fabsf(datai->ratio-1.f)<EPS){//Matched. always forward prop.
	    const float fracx=datai->dispx;
	    const float fracy=datai->dispy;
	    const float fracx1=1.f-fracx;//this reduces register usage.
	    const float fracy1=1.f-fracy;
	    /*During reverse operation, for different iwfs, the offo is
	      different causing same thread to handle different memory in out
	      for different iwfs. This causes problem without synchronization/atomic operation*/
	    if(datai->trans=='r'){
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			atomicAdd(&out[ix+iy*nxo], 
				  alpha*(+(in[ix+    iy*nxi]*fracx1+in[ix+1+    iy*nxi]*fracx)*fracy1
					 +(in[ix+(iy+1)*nxi]*fracx1+in[ix+1+(iy+1)*nxi]*fracx)*fracy));
		    }
		}
	    }else{
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			out[ix+iy*nxo]+=
			    alpha*(+(in[ix+    iy*nxi]*fracx1+in[ix+1+    iy*nxi]*fracx)*fracy1
				   +(in[ix+(iy+1)*nxi]*fracx1+in[ix+1+(iy+1)*nxi]*fracx)*fracy);
		    }
		}
	    }
	}else{//Generic
	    const float ratio=datai->ratio;
	    const float dispx=datai->dispx;
	    const float dispy=datai->dispy;
	    if(datai->trans=='t'){
		for(int iy=iy0; iy<ny; iy+=stepy){
		    float fracy;
		    int ky;
		    {
			float temp;
			fracy=modff(dispy+iy*ratio, &temp);
			ky=(int)temp;
		    }
		    for(int ix=ix0; ix<nx; ix+=stepx){
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
	    }else{
		for(int iy=iy0; iy<ny; iy+=stepy){
		    float fracy; int ky;
		    {
			float jy;
			fracy=modff(dispy+iy*ratio, &jy);
			ky=(int)jy;
		    }
		    for(int ix=ix0; ix<nx; ix+=stepx){
			float fracx; int kx;
			{
			    float jx;
			    fracx=modff(dispx+ix*ratio, &jx);
			    kx=(int)jx;
			}
			out[ix+iy*nxo]+=
			    alpha*(+(in[kx+      ky*nxi]*(1.f-fracx)+
				     in[kx+1+    ky*nxi]*fracx)*(1.f-fracy)
				   +(in[kx  +(ky+1)*nxi]*(1.f-fracx)+
				     in[kx+1+(ky+1)*nxi]*fracx)*fracy);
		    }
		}
	    }
	}
    }
}

/*
  The third grid dimension tells the wfs to handle. 
*/
#define DIM_GP 128
__global__ static void gpu_gp_do(GPU_GP_T *data, float **gout, float *ttout, float *dfout, float **wfsopd, int ptt){
    __shared__ float gx[DIM_GP];
    __shared__ float gy[DIM_GP];
    __shared__ float gdf[DIM_GP];
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int nsa=datai->nsa;
    const int step=blockDim.x * gridDim.x;
    int (*restrict saptr)[2]=datai->saptr;
    float *restrict g=gout[iwfs];
    float GPscale=datai->GPscale;
    if(wfsopd){
	const float *restrict map=wfsopd[iwfs];
	const short2 *restrict pxy=datai->GPp;
	int nx=datai->nxp;
	/*GP operation.*/
	if(pos==1){
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
		int ix=saptr[isa][0];
		int iy=saptr[isa][1];
	
		const short2 *restrict pxy2=pxy+isa*4;
		g[isa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].x
		    +map[iy*nx+ix+1]*pxy2[1].x
		    +map[(iy+1)*nx+ix ] *pxy2[2].x
		    +map[(iy+1)*nx+ix+1]*pxy2[3].x);

		g[isa+nsa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].y
		    +map[iy*nx+ix+1]*pxy2[1].y
		    +map[(iy+1)*nx+ix ] *pxy2[2].y
		    +map[(iy+1)*nx+ix+1]*pxy2[3].y);
	    }/*for isa */
	}else{
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
		int ix=saptr[isa][0];
		int iy=saptr[isa][1];
		const short2 *restrict pxy2=pxy+isa*9;
		g[isa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].x
		    +map[iy*nx+ix+1]*pxy2[1].x
		    +map[iy*nx+ix+2]*pxy2[2].x
		    +map[(iy+1)*nx+ix ] *pxy2[3].x
		    +map[(iy+1)*nx+ix+1]*pxy2[4].x
		    +map[(iy+1)*nx+ix+2]*pxy2[5].x
		    +map[(iy+2)*nx+ix ] *pxy2[6].x
		    +map[(iy+2)*nx+ix+1]*pxy2[7].x
		    +map[(iy+2)*nx+ix+2]*pxy2[8].x);
		g[isa+nsa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].y
		    +map[iy*nx+ix+1]*pxy2[1].y
		    +map[iy*nx+ix+2]*pxy2[2].y
		    +map[(iy+1)*nx+ix ] *pxy2[3].y
		    +map[(iy+1)*nx+ix+1]*pxy2[4].y
		    +map[(iy+1)*nx+ix+2]*pxy2[5].y
		    +map[(iy+2)*nx+ix ] *pxy2[6].y
		    +map[(iy+2)*nx+ix+1]*pxy2[7].y
		    +map[(iy+2)*nx+ix+2]*pxy2[8].y);
	    }/*for isa */
	}
    }
    /* Global TT, Diff-Focus projection. Modifed from previous kernel so that
       each thread handle the same subaperture as previous gradient operation to
       avoid synchronization */
    if(datai->PTT && ptt){ //temp
	float (*restrict PTT)[2]=(float(*)[2])datai->PTT;
	gx[threadIdx.x]=0;
	gy[threadIdx.x]=0;
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
	    float tmpx=PTT[isa][0]*g[isa];
	    float tmpy=PTT[isa][1]*g[isa];
	    gx[threadIdx.x]+=tmpx+PTT[isa+nsa][0]*g[isa+nsa];
	    gy[threadIdx.x]+=tmpy+PTT[isa+nsa][1]*g[isa+nsa];
	}
	for(int step=(DIM_GP>>1); step>0; step>>=1){
	    __syncthreads();
	    if(threadIdx.x<step){
		gx[threadIdx.x]+=gx[threadIdx.x+step];
		gy[threadIdx.x]+=gy[threadIdx.x+step];
	    }
	}
	if(threadIdx.x==0){
	    atomicAdd(&ttout[iwfs*2], -gx[0]);
	    atomicAdd(&ttout[iwfs*2+1], -gy[0]);
	}
    }
    if(datai->PDF && ptt){
	float *restrict PDF=datai->PDF;//the diagonal block
	const int ipowfs=datai->ipowfs;
	float scale=-1./(datai->nwfs-1);//PDF*scale gives the off diagonal block
	for(int irow=0; irow<nwfs; irow++){//skip first row
	    if(data[irow].ipowfs!=ipowfs || data[irow].jwfs==0) continue;//different group or first wfs.
	    gdf[threadIdx.x]=0;
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
		gdf[threadIdx.x]+=PDF[isa]*g[isa]+PDF[isa+nsa]*g[isa+nsa];
	    }
	    for(int step=(DIM_GP>>1); step>0; step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
		    gdf[threadIdx.x]+=gdf[threadIdx.x+step];
		}
	    }
	    if(threadIdx.x==0){
		if(datai->PTT){//adjust by TT coupling
		    gdf[0]-=datai->PDFTT[0]*gx[0]+datai->PDFTT[1]*gy[0];
		}
		if(irow!=iwfs){
		    gdf[0]*=scale;
		}
		atomicAdd(&dfout[irow], -gdf[0]);
	    }
	}
    }
}
__global__ static void gpu_gpt_do(GPU_GP_T *data, float **wfsopd, float *ttin, float *dfin, float **gin, int ptt){
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int step=blockDim.x * gridDim.x;
    const int nsa=datai->nsa;
    int (*saptr)[2]=datai->saptr;
    const float (*restrict neai)[3]=datai->neai;
    float dxp=datai->dxp;
    float oxp=datai->oxp;
    float oyp=datai->oyp;
    float focus=0;
    if(datai->PDF && ptt){
	if(iwfs==0){
	    for(int id=1; id<nwfs; id++){
		focus+=dfin[id];
	    }
	}else{
	    focus=-dfin[iwfs];
	}
    }
    const float *restrict g=gin[iwfs];
    float *restrict map=wfsopd[iwfs];
    float GPscale=datai->GPscale;
    const short2 *restrict pxy=datai->GPp;
    float ttx=0, tty=0;
    if(datai->PTT && ptt){
	ttx=ttin[iwfs*2+0];
	tty=ttin[iwfs*2+1];
    }
    const int nx=datai->nxp;
    if(pos==1){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    float cx=neai[isa][0];
	    float cy=neai[isa][1];
	    float cxy=neai[isa][2];
	    float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    float tmp=cxy*gx;//save data
	    gx=GPscale*(cx*gx+cxy*gy);
	    gy=GPscale*(tmp+cy*gy);
	    const short2 *restrict pxy2=pxy+isa*4;
	    atomicAdd(&map[iy    *nx+ix],   gx*pxy2[0].x + gy*pxy2[0].y);
	    atomicAdd(&map[iy    *nx+ix+1], gx*pxy2[1].x + gy*pxy2[1].y);
	    atomicAdd(&map[(iy+1)*nx+ix],   gx*pxy2[2].x + gy*pxy2[2].y);
	    atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[3].x + gy*pxy2[3].y);
	}
    }else if(pos==2){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    float cx=neai[isa][0];
	    float cy=neai[isa][1];
	    float cxy=neai[isa][2];
	    float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    float tmp=cxy*gx;
	    gx=GPscale*(cx*gx+cxy*gy);
	    gy=GPscale*(tmp+cy*gy);
	    const short2 *restrict pxy2=pxy+isa*9;
	    atomicAdd(&map[iy    *nx+ix],   gx*pxy2[0].x + gy*pxy2[0].y);
	    atomicAdd(&map[iy    *nx+ix+1], gx*pxy2[1].x + gy*pxy2[1].y);
	    atomicAdd(&map[iy    *nx+ix+2], gx*pxy2[2].x + gy*pxy2[2].y);
	    atomicAdd(&map[(iy+1)*nx+ix],   gx*pxy2[3].x + gy*pxy2[3].y);
	    atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[4].x + gy*pxy2[4].y);
	    atomicAdd(&map[(iy+1)*nx+ix+2], gx*pxy2[5].x + gy*pxy2[5].y);
	    atomicAdd(&map[(iy+2)*nx+ix],   gx*pxy2[6].x + gy*pxy2[6].y);
	    atomicAdd(&map[(iy+2)*nx+ix+1], gx*pxy2[7].x + gy*pxy2[7].y);
	    atomicAdd(&map[(iy+2)*nx+ix+2], gx*pxy2[8].x + gy*pxy2[8].y);
	}
    }
}
/*Only do tt, NEA. do not to GP'*/
__global__ static void gpu_nea_do(GPU_GP_T *data, float *ttin, float *dfin, float **gin){
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int step=blockDim.x * gridDim.x;
    const int nsa=datai->nsa;
    int (*saptr)[2]=datai->saptr;
    const float (*restrict neai)[3]=datai->neai;
    float dxp=datai->dxp;
    float oxp=datai->oxp;
    float oyp=datai->oyp;
    float focus=0;
    if(datai->PDF){
	if(iwfs==0){
	    for(int id=1; id<nwfs; id++){
		focus+=dfin[id];
	    }
	}else{
	    focus=-dfin[iwfs];
	}
    }
    float *restrict g=gin[iwfs];
    float ttx=ttin[iwfs*2+0];
    float tty=ttin[iwfs*2+1];
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	float cx=neai[isa][0];
	float cy=neai[isa][1];
	float cxy=neai[isa][2];
	float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	g[isa]=cx*gx+cxy*gy;
	g[isa+nsa]=cxy*gx+cy*gy;
    }
}

/*
  Tomography right hand size matrix. Computes xout = xout *beta + alpha * Hx' G' C * xin.
  xout is zeroed out before accumulation.
*/
void gpu_TomoR(curcell **xout, float beta, const void *A, const curcell *grad, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }
    curcell *opdx=*xout;
    const int nwfs=grad->nx;
    curcell *opdwfs=curecon->opdwfs;
    curzero(opdwfs->m, stream);
    if(fabsf(beta)<EPS){
	curzero(opdx->m, stream);
    }else if(fabsf(beta-1)>EPS){
	curscale(opdx->m, beta, stream);
    }
    curmat *ttf=curecon->ttf;
    curzero(ttf, stream);
    //just does low rank terms
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, NULL, 1);
    gpu_gpt_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, opdwfs->pm, ttf->p, ttf->p+nwfs*2, grad->pm, 1);
    gpu_prop_grid_adaptive_do<<<dim3(3,3,recon->npsr), dim3(16,16), 0, stream>>>
	(curecon->hxtdata, opdwfs->pm, opdx->pm, nwfs, recon->npsr, alpha, 't');
}

void gpu_TomoRt(curcell **gout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int nwfs=parms->nwfsr;
    if(!*gout){
	*gout=curcellnew(parms->nwfs, 1, recon->ngrad, (long*)NULL);
    }
    curcell *grad=*gout;
    curcell *opdwfs=curecon->opdwfs;
    curmat *ttf=curecon->ttf;
    curzero(opdwfs->m, stream);
    gpu_prop_grid_adaptive_do<<<dim3(3,3, nwfs), dim3(16,16), 0, stream>>>
      (curecon->hxdata, opdwfs->pm, xin->pm, nwfs, recon->npsr, alpha, 'n');
    curzero(ttf, stream);
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, opdwfs->pm, 1);
    gpu_nea_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, ttf->p, ttf->p+nwfs*2, grad->pm);
}

/*
  Tomography left hand size matrix. Computes xout = beta*xout + alpha * Hx' G' C Gp Hx * xin.
  xout is zeroed out before accumulation.
*/
void gpu_TomoL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int nwfs=parms->nwfsr;
    curcell *grad=curecon->grad;
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }
#if TIMING==2
    EVENT_INIT(6);
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    curcell *opdx=*xout;
    curcell *opdwfs=curecon->opdwfs;
    int ptt=!parms->recon.split || parms->dbg.splitlrt; 
    curmat *ttf=curecon->ttf;
    curzero(opdwfs->m, stream);
    gpu_prop_grid_adaptive_do<<<dim3(3,3, nwfs), dim3(16,16), 0, stream>>>
	(curecon->hxdata, opdwfs->pm, xin->pm, nwfs, recon->npsr, 1, 'n');
    RECORD(1);
    curzero(ttf, stream);
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, opdwfs->pm, ptt);
    RECORD(2);
    curzero(opdwfs->m, stream);
    gpu_gpt_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(curecon->gpdata, opdwfs->pm, ttf->p, ttf->p+nwfs*2, grad->pm, ptt);
    RECORD(3);
    if(fabsf(beta)<EPS){
	curzero(opdx->m, stream);
    }else if(fabsf(beta-1.f)>EPS){
	curscale(opdx->m, beta, stream);
    }
    gpu_prop_grid_adaptive_do<<<dim3(3,3,recon->npsr), dim3(16,16), 0, stream>>>
	(curecon->hxtdata, opdwfs->pm, opdx->pm, nwfs, recon->npsr, alpha, 't');
    RECORD(4);
    /*This could be in parallel to above ones. */
    gpu_laplacian_do<<<dim3(3,3,recon->npsr),dim3(16,16), 0, stream>>>
	(curecon->hxdata, opdx->pm, xin->pm, nwfs, alpha);
    RECORD(5);
#if TIMING==2
    EVENT_TOC;
    info2("TomoL: Hx %.0f, Gp %.0f, Gpt %.0f, Hxt %.0f, L2 %.0f\n", 
	  times[1], times[2], times[3], times[4], times[5]);
#endif
    //overhead of TomoL 27 micro-seconds (timing without synchornization).
}
/**
   Wrap of the tomography operation

   grad is the gradient input.
   opdx is the right hand side vector computed from grad. Allow NULL.
   opdr is the tomography result.
*/
double gpu_tomo_do(const PARMS_T *parms,const RECON_T *recon, curcell *opdr, curcell *opdx, curcell *grad, stream_t &stream){
    curcell *rhs=NULL;
    if(opdx){
	rhs=opdx;
    }
    double res=0;
    curecon_t *curecon=cudata->recon;
    curmat *tmp=NULL;
#if TIMING>=2
    EVENT_INIT(3);
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    gpu_TomoR(&rhs, 0, recon, grad, 1, stream);
    RECORD(1);
    switch(parms->tomo.alg){
    case 0:
	if(!opdr->m){
	    error("opdr must be continuous\n");
	}
	if(!rhs->m){
	    error("rhs must be continuous\n");
	}
	cuchol_solve(opdr->m->p, curecon->RCl, curecon->RCp, rhs->m->p, stream);
	if(curecon->RUp){
	    tmp=curnew(curecon->RVp->ny, 1);
	    curmv(tmp->p, 0, curecon->RVp, rhs->m->p, 't', -1, stream);
	    curmv(opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, stream);
	}
	break;
    case 1:{
	G_PREFUN prefun=NULL;
	void *predata=NULL;
	if(parms->tomo.precond==1){
	    prefun=gpu_Tomo_fdprecond;
	    predata=(void*)recon;
	}
	if((res=gpu_pcg(&opdr, gpu_TomoL, recon, prefun, predata, rhs, &curecon->cgtmp_tomo, 
			parms->recon.warm_restart, parms->tomo.maxit, stream))>1){
	    warning("Tomo CG not converge.\n");
	}
    }break;
    case 2:
	curmv(opdr->m->p, 0, curecon->RMI, rhs->m->p, 'n', 1, stream);
	break;
    default:
	error("Invalid");
    }
#if TIMING>=2
    RECORD(2);
    EVENT_TOC;
    info2("Tomo RHS: %6.0f, LHS: %6.0f, Tot: %6.0f\n", times[1], times[2], times[0]);
#endif
    if(!opdx){
	curcellfree(rhs);
    }
    if(tmp) curfree(tmp);
    return res;
}
void gpu_tomo_test(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    const PARMS_T *parms=simu->parms;
    stream_t &stream=curecon->cgstream[0];//had to be cgstream because fdpcg uses it intrisically
    RECON_T *recon=simu->recon;
    /*Debugging. */
    dcell *rhsc=NULL;
    dcell *lc=NULL;
    dcell *rtc=NULL;
    curcell *rhsg=NULL;
    curcell *lg=NULL;
    curcell *rtg=NULL;
    if(1){
	muv(&rhsc, &recon->RR, simu->gradlastol, 1);
    }else{
	rhsc=dcellread("../load_rhs");
	cp2gpu(&rhsg, rhsc);
	dcellscale(rhsc, 1.e12);
	for(int i=0; i<rhsc->nx; i++){
	    if(rhsc->p[i]->nx != recon->xnx[i]*recon->xny[i]){
		error("Loaded RHS has wrong dimension\n");
	    }
	    rhsc->p[i]->nx=rhsc->p[i]->nx*rhsc->p[i]->ny;
	    rhsc->p[i]->ny=1;
	    rhsg->p[i]->nx=recon->xnx[i];
	    rhsg->p[i]->ny=recon->xny[i];
	}
    }
    dcellwrite(rhsc, "CPU_TomoR");
    muv_trans(&rtc, &recon->RR, rhsc, 1);
    dcellwrite(rtc, "CPU_TomoRt");
    muv(&lc, &recon->RL, rhsc, 1);
    dcellwrite(lc, "CPU_TomoL");
    muv(&lc, &recon->RL, rhsc, -1);
    dcellwrite(lc, "CPU_TomoL2");
    dcellzero(lc);
    muv(&lc, &recon->RL, rhsc, 1);
    dcellwrite(lc, "CPU_TomoL3");
    if(parms->tomo.alg==1 && parms->tomo.precond==1){
	dcell *lp=NULL;
	fdpcg_precond(&lp, recon, lc);
	dcellwrite(lp, "CPU_TomoP");
	fdpcg_precond(&lp, recon, lc);
	dcellwrite(lp, "CPU_TomoP2");
    }
    dcellzero(lc);
    recon->desplitlrt=1;//temporary.
    for(int i=0; i<10; i++){
	muv_solve(&lc, &recon->RL, NULL, rhsc);
	dcellwrite(lc, "CPU_Tomo_%d", i);
    }
	
    if(!rhsg){
	cp2gpu(&curecon->gradin, simu->gradlastol);
	gpu_TomoR(&rhsg, 0, recon, curecon->gradin, 1, stream);
    }
    curcellwrite(rhsg, "GPU_TomoR");
    gpu_TomoR(&rhsg, 0, recon, curecon->gradin, 1, stream);
    curcellwrite(rhsg, "GPU_TomoR1");
    gpu_TomoRt(&rtg, 0, recon, rhsg, 1, stream);
    curcellwrite(rtg, "GPU_TomoRt");
    gpu_TomoL(&lg, 0, recon, rhsg, 1,stream);
    curcellwrite(lg, "GPU_TomoL");
    gpu_TomoL(&lg, 0, recon, rhsg, 1,stream);
    curcellwrite(lg, "GPU_TomoL1"); 
    gpu_TomoL(&lg, 1, recon, rhsg, -1,stream);
    curcellwrite(lg, "GPU_TomoL2");
    gpu_TomoL(&lg, 0, recon, rhsg, 1,stream);
    curcellwrite(lg, "GPU_TomoL3");
    if(parms->tomo.alg==1 && parms->tomo.precond==1){
	curcell *lp=NULL;
	gpu_Tomo_fdprecond(&lp, recon, lg, stream);
	curcellwrite(lp, "GPU_TomoP");
	gpu_Tomo_fdprecond(&lp, recon, lg, stream);
	curcellwrite(lp, "GPU_TomoP2");
	//exit(0);
    }
    G_PREFUN prefun=NULL;
    void *predata=NULL;
    if(parms->tomo.precond==1){
	prefun=gpu_Tomo_fdprecond;
	predata=(void*)recon;
    }
    curcellzero(lg, 0);
    for(int i=0; i<10; i++){
	gpu_pcg(&lg, (G_CGFUN)gpu_TomoL, (void*)recon, prefun, predata, rhsg, &curecon->cgtmp_tomo,
		simu->parms->recon.warm_restart, parms->tomo.maxit, stream);
	CUDA_SYNC_STREAM;
	curcellwrite(lg, "GPU_Tomo_%d", i);
    }
    CUDA_SYNC_DEVICE;
    exit(0);
}
