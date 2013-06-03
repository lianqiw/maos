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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include <cusparse.h>
#include <cufft.h>
#include "cucmat.h"
#include "wfs.h"
/*
  Timing results for TMT NFIRAOS case per LGS WFS:
  Embedding takes about 1 ms.
  FFT takes about 2 ms
  CCWM takes about 1ms
  realpart takes about 1ms

  Total takes about 12 ms.
*/

#undef TIMING
#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define ctoc(A)
#else
#define ctoc(A) CUDA_SYNC_STREAM; toc2(A)
#endif

/**
   Embed amp*exp(2*pi*i*opd). input is nxin*nxin, output is nxout*nxout;
*/
__global__ static void sa_embed_wvf_do(fcomplex *restrict wvf, 
				    const float *restrict opd, const float *restrict amp, 
				    const float wvl, const int nxin, const int nxout){
    const int isa=blockIdx.x;
    const int pad=(nxout-nxin)>>1;
    const int skipin=isa*nxin*nxin;
    const int skipout=isa*nxout*nxout+pad;
    const float pi2l=2.f*M_PI/wvl;
    for(int iy=threadIdx.y; iy<nxin; iy+=blockDim.y){
	const int skipin2=skipin+iy*nxin;
	const int skipout2=skipout+(iy+pad)*nxout;
	for(int ix=threadIdx.x; ix<nxin; ix+=blockDim.x){
	    /*test sinpi later */
	    float s,c;
	    sincosf(pi2l*opd[skipin2+ix], &s, &c);
	    wvf[skipout2+ix]=make_cuComplex(amp[skipin2+ix]*c, amp[skipin2+ix]*s);
	}
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void sa_cpcorner_do(fcomplex *restrict out, int noutx,  int nouty,
				   const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    const int isa=blockIdx.x;
    out+=isa*noutx*nouty;
    in+=isa*ninx*niny;
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    out[iy*noutx+ix]                   = in[iy*ninx+ix];
	    out[iy*noutx+(noutx-1-ix)]         = in[iy*ninx+(ninx-1-ix)];
	    out[(nouty-1-iy)*noutx+(noutx-1-ix)] = in[(niny-1-iy)*ninx+(ninx-1-ix)];
	    out[(nouty-1-iy)*noutx+(ix)]       = in[(niny-1-iy)*ninx+(ix)];
	}
    }
}

/**
   Embed or crop an array to another array. Preserve center. 
*/
__global__ void sa_cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
			    const fcomplex *restrict in, int ninx, int niny, float scale){
    int nx, ny, nskipoutx, nskipouty, nskipinx, nskipiny;
    if(noutx<ninx){
	nx=noutx;
	nskipoutx=0;
	nskipinx=(ninx-noutx)>>1;
    }else{
	nx=ninx;
	nskipoutx=(noutx-ninx)>>1;
	nskipinx=0;
    }
    if(nouty<niny){
	ny=nouty;
	nskipouty=0;
	nskipiny=(niny-nouty)>>1;
    }else{
	ny=niny;
	nskipouty=(nouty-niny)>>1;
	nskipiny=0;
    }
    const int isa=blockIdx.x;
    out+=isa*noutx*nouty+nskipouty*(noutx)+nskipoutx;
    in+=isa*ninx*niny+nskipiny*(ninx)+nskipinx;
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    out[iy*noutx+ix]=make_cuFloatComplex(scale*cuCrealf(in[iy*ninx+ix]),scale*cuCimagf(in[iy*ninx+ix]));
	}
    }
}
/**
   abs2 to real.
*/
__global__ static void sa_abs2real_do(fcomplex *wvf, const int nx, float alpha){
    const int isa=blockIdx.x;
    wvf+=nx*nx*isa;
    for(int iy=threadIdx.y; iy<nx; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    float r=cuCrealf(wvf[iy*nx+ix]);
	    float i=cuCimagf(wvf[iy*nx+ix]);
	    wvf[iy*nx+ix]=make_cuFloatComplex((r*r+i*i)*alpha, 0);
	}
    }
}
/**
   FFT Shift.
*/
__global__ static void sa_fftshift_do(fcomplex *wvf, const int nx, const int ny){
    const int isa=blockIdx.x;
    wvf+=nx*ny*isa;
    int nx2=nx>>1;
    int ny2=ny>>1;
    for(int iy=threadIdx.y; iy<ny2; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx2; ix+=blockDim.x){
	    fcomplex tmp;
	    tmp=wvf[ix+iy*nx];
	    wvf[ix+iy*nx]=wvf[(ix+nx2)+(iy+ny2)*nx];
	    wvf[(ix+nx2)+(iy+ny2)*nx]=tmp;
	    tmp=wvf[ix+(iy+ny2)*nx];
	    wvf[ix+(iy+ny2)*nx]=wvf[(ix+nx2)+iy*nx];
	    wvf[(ix+nx2)+iy*nx]=tmp;
	}
    }
}
/**
   FFT Shift from complex to real.
*/
__global__ static void sa_acc_real_fftshift_do(float *restrict out, const fcomplex *restrict wvf, 
					    int nx, int ny, float alpha){
    const int isa=blockIdx.x;
    wvf+=nx*ny*isa;
    out+=nx*ny*isa;
    int nx2=nx>>1;
    int ny2=ny>>1;
    for(int iy=threadIdx.y; iy<ny2; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx2; ix+=blockDim.x){
	    out[ix+iy*nx]+=alpha*cuCrealf(wvf[(ix+nx2)+(iy+ny2)*nx]);
	    out[(ix+nx2)+(iy+ny2)*nx]+=alpha*cuCrealf(wvf[ix+iy*nx]);
	    out[ix+(iy+ny2)*nx]+=alpha*cuCrealf(wvf[(ix+nx2)+iy*nx]);
	    out[(ix+nx2)+iy*nx]+=alpha*cuCrealf(wvf[ix+(iy+ny2)*nx]);
	}
    }
}
/**
   Rotate and embed.
 */
__global__ static void sa_embed_rot_do(fcomplex *restrict out, const int noutx, const int nouty,
				    const fcomplex *restrict in, const int ninx, const int niny, const float* srot){
    const int isa=blockIdx.x;
    out+=isa*noutx*nouty;
    in+=isa*ninx*niny;
    float theta=srot[isa];
    float sx, cx;
    sincosf(theta, &sx, &cx);
    const int niny2=niny/2;
    const int ninx2=ninx/2;
    const int nouty2=nouty/2;
    const int noutx2=noutx/2;

    for(int iy=threadIdx.y; iy<nouty; iy+=blockDim.y){
	int y0=iy-nouty2;
	for(int ix=threadIdx.x; ix<noutx; ix+=blockDim.x){
	    int x0=ix-noutx2;
	    float x=(cx*x0-sx*y0)+ninx2;
	    float y=(sx*x0+cx*y0)+niny2;
	    int jx=floorf(x); x=x-jx;
	    int jy=floorf(y); y=y-jy;
	    if(jx>=0 && jx<ninx-1 && jy>=0 && jy<niny-1){
		out[iy*noutx+ix]=
		    make_cuFloatComplex((cuCrealf(in[jy*ninx+jx])*(1-x)
					 +cuCrealf(in[jy*ninx+jx+1])*x)*(1-y)
					+(cuCrealf(in[(jy+1)*ninx+jx])*(1-x)
					  +cuCrealf(in[(jy+1)*ninx+jx+1])*x)*y, 0);
	    }
	}
    }
} 
/**
   Multiple each OTF with another. 
*/
__global__ static void sa_ccwm_do(fcomplex *otf, const int notfx, const int notfy, 
			       fcomplex **lotfcs, int each){
    const int isa=blockIdx.x;
    otf+=notfx*notfy*isa;
    const fcomplex *restrict lotfc=each?(fcomplex*)lotfcs:lotfcs[isa];
    for(int iy=threadIdx.y; iy<notfy; iy+=blockDim.y){
	const int skip=iy*notfx;
	fcomplex *restrict otf2=otf+skip;
	const fcomplex *restrict lotfc2=lotfc+skip; 
	for(int ix=threadIdx.x; ix<notfx; ix+=blockDim.x){
	    otf2[ix]=cuCmulf(otf2[ix], lotfc2[ix]);
	}
    }
}
/**
   Multiple an otf with another 1-d otf along each column
*/
__global__ static void sa_ccwmcol_do(fcomplex *otf, const int notfx, const int notfy,
				  fcomplex *const *etfs, int each){
    const int isa=blockIdx.x;
    otf+=notfy*notfx*isa;
    const fcomplex *restrict etf=etfs[each?0:isa];
    for(int iy=threadIdx.y; iy<notfy; iy+=blockDim.y){
	fcomplex *restrict otf2=otf+iy*notfx;
	for(int ix=threadIdx.x; ix<notfx; ix+=blockDim.x){
	    otf2[ix]=cuCmulf(otf2[ix], etf[ix]);
	}
    }
}
/**
   Take the real part. Notice we are using = instead of +=
*/
__global__ static void sa_realpart_do(float *out, const fcomplex*restrict in, int ninx, int niny){
    const int isa=blockIdx.x;
    in+=isa*ninx*niny;
    out+=isa*ninx*niny;
    for(int iy=threadIdx.y; iy<niny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<ninx; ix+=blockDim.x){
	    out[ix+iy*ninx]=cuCrealf(in[ix+iy*ninx]);
	}
    }
}
/**
   Take the real part and accumulate to output
*/
__global__ static void sa_acc_real_do(float *out, const fcomplex*restrict in, int ninx, int niny, float alpha){
    const int isa=blockIdx.x;
    in+=isa*ninx*niny;
    out+=isa*ninx*niny;
    for(int iy=threadIdx.y; iy<niny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<ninx; ix+=blockDim.x){
	    out[ix+iy*ninx]+=cuCrealf(in[ix+iy*ninx])*alpha;
	}
    }
}
/**
   Do the role of si. input psfr is sampled with notfx*notfy, with sampling dtheta.
   Output ints is sampled with pixpsax*pixpsay, at pixtheta.
*/
__global__ static void sa_si_rot_do(float *restrict ints, int pixpsax, int pixpsay, 
				 int pixoffx, int pixoffy, float pixthetax, float pixthetay,
				 const float *restrict psfr, float dtheta, int notfx, int notfy,
				 const float *restrict srot, float alpha){
    float pxo=-(pixpsax*0.5-0.5+pixoffx)*pixthetax;
    float pyo=-(pixpsay*0.5-0.5+pixoffy)*pixthetay;
    int isa=blockIdx.x;
    float dispx=notfx/2;
    float dispy=notfy/2;
    float dtheta1=1.f/dtheta;
    float sx=0, cx=1.f;
    if(srot){
	sincosf(srot[isa], &sx, &cx);
    }
    ints+=isa*pixpsax*pixpsay;
    psfr+=isa*notfx*notfy;
    for(int iy=threadIdx.y; iy<pixpsay; iy+=blockDim.y){
	float y0=iy*pixthetay+pyo;
	for(int ix=threadIdx.x; ix<pixpsax; ix+=blockDim.x){
	    float x0=ix*pixthetax+pxo;
	    float x=(cx*x0-sx*y0)*dtheta1+dispx;
	    float y=(sx*x0+cx*y0)*dtheta1+dispy;
	    int jx=floorf(x); x=x-jx;
	    int jy=floorf(y); y=y-jy;
	    if(jx>=0 && jx<notfx-1 && jy>=0 && jy<notfy-1){
		ints[iy*pixpsax+ix]+=alpha*((psfr[jy*notfx+jx]*(1.-x)+psfr[jy*notfx+jx+1]*x)*(1.-y)
					    +(psfr[(jy+1)*notfx+jx]*(1.-x)+psfr[(jy+1)*notfx+jx+1]*x)*y);
	    }
	}
    }
}

/**
   Add tip/tilt to the OTF for each subaps. exp(-2*pi*sx/nx)*exp(-2*pi*sy/ny).
   peak of otf is in corner.
 */
__global__ static void sa_add_otf_tilt_corner_do(fcomplex *restrict otf, int nx, int ny, 
					      float *restrict gx, float *restrict gy, float gscale){
    int isa=blockIdx.x;
    float sx=gx[isa]*gscale;
    float sy=gy[isa]*gscale;
    fcomplex *restrict otfi=otf+isa*nx*ny;
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    float phase=-2.f*M_PI*(+(float)(ix-(ix*2>=nx?nx:0))/(float)nx*sx
				   +(float)(iy-(iy*2>=ny?ny:0))/(float)ny*sy);
	    float s,c;
	    sincosf(phase, &s, &c);
	    fcomplex otfii=otfi[ix+iy*nx];
	    otfi[ix+iy*nx]=make_cuFloatComplex(cuCrealf(otfii)*c - cuCimagf(otfii)*s,
					       cuCrealf(otfii)*s + cuCimagf(otfii)*c);
	}
    }
}
/**
   Do physical wfs images in GPU. please check wfsints() in CPU code for comments.
*/
void gpu_wfsints(SIM_T *simu, float *phiout, curmat *gradref, int iwfs, int isim, cudaStream_t stream){
    TIC;tic;
    cuwloc_t *cupowfs=cudata->powfs;
    cuwfs_t *cuwfs=cudata->wfs;
    const PARMS_T *parms=simu->parms;
    const POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    const float hs=parms->powfs[ipowfs].hs;
    const float dtisim=parms->sim.dt*isim;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int ncompx=powfs[ipowfs].ncompx;/*necessary size to build detector image. */
    const int ncompy=powfs[ipowfs].ncompy;
    const int notf=MAX(ncompx,ncompy);
    const int nx=powfs[ipowfs].pts->nx;
    const int nwvf=nx*parms->powfs[ipowfs].embfac;
    const int pixpsax=powfs[ipowfs].pixpsax;
    const int pixpsay=powfs[ipowfs].pixpsay;
    const float pixthetax=parms->powfs[ipowfs].radpixtheta;
    const float pixthetay=parms->powfs[ipowfs].pixtheta;
    const float siglev=parms->wfs[iwfs].siglevsim;
    const float *restrict const srot1=parms->powfs[ipowfs].radrot?cuwfs[iwfs].srot:NULL;
    const int multi_dtf=(parms->powfs[ipowfs].llt&&!parms->powfs[ipowfs].radrot 
			 && parms->powfs[ipowfs].radpix);
    const float *restrict const srot2=multi_dtf?cuwfs[iwfs].srot:NULL;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    curcell *restrict const ints=cuwfs[iwfs].ints;
    curmat *lltopd=NULL;
    curcell *pistatout=NULL;
    if(parms->powfs[ipowfs].pistatout && isim>=parms->powfs[ipowfs].pistatstart){
	pistatout=cuwfs[iwfs].pistatout;
    }
    cuccell *wvfout=NULL;
    const int wvf_n=notf/2+2;//was notf/2
    if(parms->powfs[ipowfs].psfout){
	wvfout=cuccellnew(nsa, nwvl, wvf_n, wvf_n);
    }
    float norm_psf=sqrt(powfs[ipowfs].areascale)/((float)powfs[ipowfs].pts->nx*nwvf);
    float norm_pistat=norm_psf*norm_psf/((float)notf*notf);
    float norm_ints=siglev*norm_psf*norm_psf/((float)ncompx*ncompy);
    /* Do msa subapertures in a batch to avoid using too much memory.*/
    fcomplex *psf, *wvf, *otf, *psfout=NULL, *psfstat=NULL;

    fcomplex *lotfc=NULL;
    fcomplex *lwvf=NULL;

    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
	int nlx=powfs[ipowfs].llt->pts->nx;
	lltopd=curnew(nlx, nlx);
	if(cuwfs[iwfs].lltncpa){
	    cudaMemcpyAsync(lltopd->p, cuwfs[iwfs].lltncpa, 
			    sizeof(float)*nlx*nlx, MEMCPY_D2D, stream);
	}else{
	    cudaMemsetAsync(lltopd->p, 0, sizeof(float)*nlx*nlx, stream);
	}
	const int illt=parms->powfs[ipowfs].llt->i[wfsind];
	const double thetaxl=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox[illt]/hs;
	const double thetayl=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy[illt]/hs;
	gpu_atm2loc(lltopd->p, cupowfs[ipowfs].llt->loc, cupowfs[ipowfs].llt->nloc, 
		    hs, thetaxl, thetayl, 
		    parms->powfs[ipowfs].llt->misreg[0], 
		    parms->powfs[ipowfs].llt->misreg[1], 
		    dtisim, 1, stream);
	if((simu->uptreal && simu->uptreal->p[iwfs]) ||pistatout||parms->sim.uptideal){
	    float ttx,tty;
	    if(pistatout||parms->sim.uptideal){
		//warning("Remove tip/tilt in uplink ideally\n");
		float *lltg;
		cudaCallocHost(lltg, 2*sizeof(float), stream);
		cuztilt<<<1,dim3(16,16),0,stream>>>(lltg, lltopd->p, 1, cupowfs[ipowfs].llt->dx, 
						    cupowfs[ipowfs].llt->nxsa, cuwfs[iwfs].lltimcc,
						    cupowfs[ipowfs].llt->pts, cuwfs[iwfs].lltamp, 1.f);
		CUDA_SYNC_STREAM;
		ttx=-lltg[0];
		tty=-lltg[1];
		cudaFreeHost(lltg);
	    }else{
		ttx=simu->uptreal->p[iwfs]->p[0];
		tty=simu->uptreal->p[iwfs]->p[1];
	    }
	    /* copy uptreal to output  */
	    PDMAT(simu->uptcmds->p[iwfs], puptcmds);
	    puptcmds[isim][0]=ttx;
	    puptcmds[isim][1]=tty;
	    /* add tip/tilt to opd  */
	    const double dx=powfs[ipowfs].llt->pts->dx;
	    const double ox=powfs[ipowfs].llt->pts->origx[0];
	    const double oy=powfs[ipowfs].llt->pts->origy[0];
	    add_tilt_do<<<1, dim3(16,16), 0, stream>>>(lltopd->p, nlx, nlx, ox, oy, dx, ttx, tty);
	}/*if upt */
	ctoc("llt opd");
	int nlwvf=nlx*parms->powfs[ipowfs].embfac;
	cudaMalloc(&lwvf, nlwvf*nlwvf*sizeof(fcomplex));
	if(nlwvf != notf){
	    cudaMalloc(&lotfc, notf*notf*sizeof(fcomplex));
	}else{
	    lotfc=lwvf;
	}
	if(parms->save.wfsopd[iwfs]){
	    cellarr_cur(simu->save->wfslltopd[iwfs], isim, lltopd, stream);
	}
    }/*if has llt */

    /* Now begin physical optics preparation*/
    int isotf=(lltopd || pistatout);
    int msa=cuwfs[iwfs].msa;/* number of subaps to process at each time.*/
    cudaMalloc(&wvf, sizeof(fcomplex)*nwvf*nwvf*msa);
    if(nwvf==notf){
	psf=wvf;
    }else{
	cudaMalloc(&psf, sizeof(fcomplex)*notf*notf*msa);
    }
    if(srot1 || ncompx!=notf || ncompy!=notf){
	cudaMalloc(&otf, sizeof(fcomplex)*ncompx*ncompy*msa);
	if(isotf){/*There is an additional pair of FFT.*/
	    norm_ints/=((float)notf*notf);
	}
    }else{
	otf=psf;
    }
    if(wvfout){
	cudaMalloc(&psfout, sizeof(fcomplex)*notf*notf*msa);
    }
    if(pistatout){
	cudaMalloc(&psfstat, sizeof(fcomplex)*notf*notf*msa);
    }
    /* Now begin physical optics  */
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	float wvl=parms->powfs[ipowfs].wvl[iwvl];
	float dtheta=wvl/(nwvf*powfs[ipowfs].pts->dx);
	if(lltopd){ /*First calculate LOTF */
	    int nlx=powfs[ipowfs].llt->pts->nx;
	    int nlwvf=nlx*parms->powfs[ipowfs].embfac;
	    cudaMemsetAsync(lwvf, 0, sizeof(fcomplex)*nlwvf*nlwvf, stream);
	    if(nlwvf<notf){
		cudaMemsetAsync(lotfc, 0, sizeof(fcomplex)*notf*notf, stream);
	    }
	    sa_embed_wvf_do<<<1,dim3(16,16),0,stream>>>
		(lwvf, lltopd->p, cuwfs[iwfs].lltamp, wvl, nlx, nlwvf);
	    /*Turn to PSF. peak in corner */
	    CUFFT(cuwfs[iwfs].lltplan_wvf, lwvf, CUFFT_FORWARD);
	    sa_abs2real_do<<<1,dim3(16,16),0,stream>>>(lwvf, nlwvf, 1./(float)(nlwvf*nlwvf));
	    /*Turn to OTF. peak in corner*/
	    /*Use backward to make lotfc the conjugate of otf. peak is in corner. */
	    CUFFT(cuwfs[iwfs].lltplan_wvf, lwvf, CUFFT_INVERSE);
	    if(lwvf!=lotfc){
		sa_cpcorner_do<<<1, dim3(16,16),0,stream>>>(lotfc, notf, notf, lwvf, nlwvf, nlwvf);
	    }
	}
	ctoc("llt otf");

	for(int isa=0; isa<nsa; isa+=msa){
	    int ksa=MIN(msa, nsa-isa);/*total number of subapertures left to do. */
	    /*embed amp/opd to wvf */
	    cudaMemsetAsync(wvf, 0, sizeof(fcomplex)*ksa*nwvf*nwvf, stream);
	    if(notf>nwvf){
		cudaMemsetAsync(psf, 0, sizeof(fcomplex)*ksa*notf*notf, stream);
	    }
	    sa_embed_wvf_do<<<ksa, dim3(16,16),0,stream>>>
		(wvf, phiout+isa*nx*nx, cuwfs[iwfs].amp+isa*nx*nx, wvl, nx, nwvf);
	    ctoc("embed");
	    /* turn to complex psf, peak in corner */
	    CUFFT(cuwfs[iwfs].plan1, wvf, CUFFT_FORWARD);
	    /* copy big psf to smaller psf covering detector focal plane. */
	    if(psf!=wvf){
		sa_cpcorner_do<<<ksa, dim3(16,16),0,stream>>>
		    (psf, notf, notf, wvf, nwvf, nwvf);
	    }
	    //gpu_write(psf, notf, notf*ksa, "psf_out_1");
	    ctoc("psf");
	    if(wvfout){
		cudaMemsetAsync(psfout, 0, sizeof(fcomplex)*ksa*notf*notf, stream);
		CUFFT2(cuwfs[iwfs].plan2, psf, psfout, CUFFT_INVERSE);
		sa_cpcenter_do<<<ksa,dim3(16,16),0,stream>>>
		    (wvfout->p[isa+nsa*iwvl]->p, wvf_n, wvf_n, 
		     psfout, notf, notf, norm_psf/(notf*notf));
	    }
	    /* abs2 part to real, peak in corner */
	    sa_abs2real_do<<<ksa,dim3(16,16),0,stream>>>(psf, notf, 1);
	    ctoc("abs2real");
	    //gpu_write(psf, notf, notf*ksa, "psf_out_2");
	    if(isotf){
		/* turn to otf. peak in corner */
		CUFFT(cuwfs[iwfs].plan2, psf, CUFFT_FORWARD);
		ctoc("fft to otf");
		if(pistatout){
		    cudaMemcpyAsync(psfstat, psf, sizeof(fcomplex)*notf*notf*ksa, 
				    MEMCPY_D2D, stream);
		    sa_add_otf_tilt_corner_do<<<ksa,dim3(16,16),0,stream>>>
			(psfstat, notf,notf, gradref->p+isa, gradref->p+nsa+isa, -1.f/dtheta);
		    CUFFT(cuwfs[iwfs].plan2, psfstat, CUFFT_INVERSE);/*back to PSF. peak in corner*/
		    if(parms->sim.skysim){/*want peak in corner*/
			sa_acc_real_do<<<ksa,dim3(16,16),0,stream>>>
			    (pistatout->p[isa+nsa*iwvl]->p, psfstat, notf, notf, norm_pistat);
		    }else{/*want peak in center*/
			sa_acc_real_fftshift_do<<<ksa,dim3(16,16),0,stream>>>
			    (pistatout->p[isa+nsa*iwvl]->p, psfstat, notf, notf, norm_pistat);
		    }
		}
		if(lltopd){/*multiply with uplink otf. */
		    sa_ccwm_do<<<ksa,dim3(16,16),0,stream>>>(psf, notf, notf, (fcomplex**)lotfc, 1);
		    ctoc("ccwm with lotfc");
		}
		/* is OTF now. */
	    }
	    //gpu_write(psf, notf, notf*ksa, "psf_out_3");
	    ctoc("before ints");
	    if(ints){
		if(!isotf || otf!=psf){/*rotate PSF, or turn to OTF first time. */
		    if(isotf){/*srot1 is true. turn otf back to psf for rotation. */
			CUFFT(cuwfs[iwfs].plan2, psf, CUFFT_INVERSE);
			ctoc("fft to psf");
		    }
		    if(otf!=psf){
			cudaMemsetAsync(otf, 0, sizeof(fcomplex)*ksa*ncompx*ncompy, stream);
		    }
		    if(srot1){/*rotate and embed psf*/
			sa_fftshift_do<<<ksa, dim3(16,16),0,stream>>>
			    (psf, notf, notf);/*shift to center */
			sa_embed_rot_do<<<ksa, dim3(16,16), 0, stream>>>
			    (otf, ncompx, ncompy, psf, notf, notf, srot1?srot1+isa:NULL);
			sa_fftshift_do<<<ksa, dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy);/*shift back to corner */
		    }else if(otf!=psf){/*copy the psf corner*/
			sa_cpcorner_do<<<ksa, dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy, psf, notf, notf);
			ctoc("cpcorner");
		    }
		    /*Turn PSF to OTF. */
		    CUFFT(cuwfs[iwfs].plan3, otf,CUFFT_FORWARD);
		    ctoc("fft to otf");
		}
		/*now we have otf. multiple with etf, dtf. */
		if(cuwfs[iwfs].dtf[iwvl].etf){
		    if(cuwfs[iwfs].dtf[iwvl].etfis1d){
			sa_ccwmcol_do<<<ksa,dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].etf+isa, 0);
		    }else{
			ctoc("before ccwm");
			sa_ccwm_do<<<ksa,dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].etf+isa, 0);
			ctoc("ccwm");
		    }
		}
		/*multiple with nominal */
		if(cuwfs[iwfs].dtf[iwvl].nominal){
		    sa_ccwm_do<<<ksa,dim3(16,16),0,stream>>>
			(otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].nominal+isa, 0);
		    ctoc("nominal");
		}
		/*back to spatial domain. */
		CUFFT(cuwfs[iwfs].plan3, otf, CUFFT_INVERSE);
		ctoc("fft");
		float *psfr;
		DO(cudaMalloc(&psfr, sizeof(float)*ncompx*ncompy*ksa));
		sa_realpart_do<<<ksa, dim3(16,16),0,stream>>>(psfr, otf, ncompx, ncompy);
		ctoc("realpart");
		sa_si_rot_do<<<ksa, dim3(16,16),0,stream>>>
		    (ints->p[isa]->p, pixpsax, pixpsay, 
		     parms->powfs[ipowfs].pixoffx, parms->powfs[ipowfs].pixoffy,
		     pixthetax, pixthetay, psfr, dtheta, ncompx, ncompy, srot2?srot2+isa:NULL, 
		     norm_ints*parms->wfs[iwfs].wvlwts[iwvl]);
		ctoc("final");
		DO(cudaFree(psfr));
	    }/*if ints. */
	}/*for isa block loop */
    }/*for iwvl */
    if(lltopd){
	curfree(lltopd);
	if(lwvf!=lotfc) cudaFree(lotfc);
	cudaFree(lwvf);
    }
    if(otf!=psf) cudaFree(otf);
    if(psf!=wvf) cudaFree(psf);
    cudaFree(wvf);
    if(psfstat)  cudaFree(psfstat);
    if(parms->powfs[ipowfs].psfout){
	cellarr_cuccell(simu->save->wfspsfout[iwfs], isim, wvfout, stream);
	cuccellfree(wvfout);
	cudaFree(psfout);
    }
}
