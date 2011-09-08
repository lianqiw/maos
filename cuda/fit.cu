extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "recon.h"
#include "accphi.h"
#define SYNC_FIT  for(int ifit=0; ifit<parms->fit.nfit; ifit++){ cudaStreamSynchronize(curecon->fitstream[ifit]); }

/**



   fit.cu is not used. uses cumuv currently. Implement in the future if necessary.






 */

/*
  Todo: share the ground layer which is both matched and same.
*/
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ static void apply_W_do(float *restrict out, const float *restrict in, const int *W0f, 
				  float alpha, int nx, int n){
    const int step=blockDim.x * gridDim.x;
    for(int k=blockIdx.x * blockDim.x + threadIdx.x; k<n; k+=step){
	int i=W0f[k];
	out[i]+=alpha*(in[i]
		       +0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
		       +0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
    }
}
/*
__global__ prop_grid_cubic_do_2(float *restrict out, int nxo,
				const float *in, int nxi,
				float dispx, float dispy, float dispx2, float dispy2,
				int nx, int ny,
				float iac, float alpha){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<ny; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<nx; ix+=stepx){

	}
    }

    }*/
void gpu_prop_grid_cubic(curmat *out, float oxo, float oyo, float dxo,
			 const curmat *in, float oxi, float oyi, float dxi,
			 float dispx, float dispy,float iac,
			 float alpha,  cudaStream_t stream){
    //offset of origin.
    dispx=(dispx-oxi+oxo);
    dispy=(dispy-oyi+oyo);
    int offx1=0, offy1=0;
    const double dxo1=1./dxo;
    //if output is bigger than input.
    if(dispx<0){
	offx1=(int)ceilf(-dispx*dxo1);
	dispx+=offx1*dxo;
    }
    if(dispy<0){
	offy1=(int)ceilf(-dispy*dxo1);
	dispy+=offy1*dxo;
    }
    int nxi=in->nx;
    int nyi=in->ny;
    int nxo=out->nx;
    int nyo=out->ny;
    int nx=(int)floorf(((nxi-1)*dxi-dispx)*dxo1)+1;
    int ny=(int)floorf(((nyi-1)*dxi-dispy)*dxo1)+1;
	
    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
  
    if(fabsf(2*dxo-dxi)<EPS){
	dispx=dispx/dxi;
	dispy=dispy/dxi;
	int offx2=floor(dispx); dispx-=offx2;
	int offy2=floor(dispy); dispy-=offy2;
	float dispx2=dispx+(dispx>0.5)?-0.5:0.5;
	float dispy2=dispy+(dispy>0.5)?-0.5:0.5;
	TO_IMPLEMENT;
	/*
	prop_grid_cubic_do_2<<<DIM2(nx,ny,16,8),0,stream>>>
	    (out->p+offy1*nxo+offx1, nxo,
	     in->p+offy2*nxi+offx2, nxi,
	     dispx, dispy, dispx2, dispy2,
	     nx, ny, iac, alpha);*/
    }else{
	TO_IMPLEMENT;
    }
}
void gpu_FitR(curcell **xout, const void *A, const curcell *xin, const float alpha){
    SIM_T *simu=(SIM_T*)A;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int nfit=parms->fit.nfit;
    const int npsr=recon->npsr;
    curcell *opdfit=curecon->opdfit;
    curcell *opdfit2=curecon->opdfit2;
    float oxp=recon->pmap->ox;
    float oyp=recon->pmap->oy;
    float dxp=recon->pmap->dx;
    const int nxp=recon->pmap->nx;
    const int nyp=recon->pmap->ny;
    const int np=nxp*nyp;
    float *pis;
    cudaMalloc(&pis, nfit*sizeof(float));
    cudaMemset(pis, 0,  nfit*sizeof(float));
    for(int ifit=0; ifit<nfit; ifit++){
	double hs=parms->fit.ht[ifit];
	float thetax=parms->fit.thetax[ifit];
	float thetay=parms->fit.thetay[ifit];
	curzero(opdfit->p[ifit], curecon->fitstream[ifit]);
	curzero(opdfit2->p[ifit], curecon->fitstream[ifit]);
	for(int ips=0; ips<npsr; ips++){
	    const double ht = recon->ht->p[ips];
	    float scale=1.f-ht/hs;
	    gpu_prop_grid(opdfit->p[ifit], oxp*scale, oyp*scale, dxp*scale, 
			  xin->p[ips], recon->xmap[ips]->ox, recon->xmap[ips]->oy, recon->xmap[ips]->dx,
			  thetax*ht, thetay*ht,
			  1.f,'n', curecon->fitstream[ifit]);
	}
	inn_do<<<DIM(np, 256), 0, curecon->fitstream[ifit]>>>
	    (&pis[ifit], opdfit->p[ifit]->p, curecon->W1->p, np);
	add2_do<<<DIM(np, 256), 0, curecon->fitstream[ifit]>>>
	    (opdfit2->p[ifit]->p, curecon->W1->p, &pis[ifit], -1.f, np);
	cuspmul(opdfit2->p[ifit]->p, curecon->W0p, opdfit->p[ifit]->p, 1.f, curecon->fitsphandle[ifit]);
	apply_W_do<<<DIM(np, 256)>>>
	    (opdfit2->p[ifit]->p, opdfit->p[ifit]->p, curecon->W0f, curecon->W0v, nxp, curecon->nW0f);
    }
    SYNC_FIT;
    cudaFree(pis);
}

