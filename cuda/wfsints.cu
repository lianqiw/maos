/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "wfs.h"
/*
  Timing results for TMT NFIRAOS case per LGS WFS:
  Embedding takes about 1 ms.
  FFT takes about 2 ms
  CCWM takes about 1ms
  realpart takes about 1ms

  Total takes about 12 ms.
*/
#define ctoc(A) //CUDA_SYNC_STREAM; toc(A)
/**
   Embed amp*exp(2*pi*i*opd). input is nxin*nxin, output is nxout*nxout;
*/
__global__ static void embed_wvf_do(fcomplex *restrict wvf, 
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
	    //test sinpi later
	    float s,c;
	    sincosf(pi2l*opd[skipin2+ix], &s, &c);
	    wvf[skipout2+ix]=make_cuComplex(amp[skipin2+ix]*c, amp[skipin2+ix]*s);
	}
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void cpcorner_do(fcomplex *restrict out, int noutx,  int nouty,
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
__global__ void cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
				   const fcomplex *restrict in, int ninx, int niny){
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
	    out[iy*noutx+ix]=in[iy*ninx+ix];
	}
    }
}
/**
   abs2 to real.
*/
__global__ static void abs2real_do(fcomplex *wvf, const int nx, float alpha){
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
__global__ static void fftshift_do(fcomplex *wvf, const int nx, const int ny){
    const int isa=blockIdx.x;
    wvf+=nx*nx*isa;
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
   Rotate and embed.
 */
__global__ static void embed_rot_do(fcomplex *restrict out, const int noutx, const int nouty,
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
__global__ static void ccwm_do(fcomplex *otf, const int notfx, const int notfy, 
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
__global__ static void ccwmcol_do(fcomplex *otf, const int notfx, const int notfy,
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
   Take the real part
*/
__global__ static void realpart_do(float *out, const fcomplex*restrict in, int ninx, int niny){
    const int isa=blockIdx.x;
    in+=isa*ninx*niny;
    out+=isa*ninx*niny;
    for(int iy=threadIdx.y; iy<niny; iy+=blockDim.y){
    	const int skip=iy*ninx;
	float *restrict out2=out+skip;
	const fcomplex *restrict in2=in+skip;
	for(int ix=threadIdx.x; ix<ninx; ix+=blockDim.x){
	    out2[ix]=cuCrealf(in2[ix]);
	}
    }
}
/**
   Do the role of si. input psfr is sampled with ncompx*ncompy, with sampling dtheta.
   Output ints is sampled with pixpsax*pixpsay, at pixtheta.
*/
__global__ static void si_rot_do(float *restrict ints, int pixpsax, int pixpsay, int pixoffx, int pixoffy, float pixtheta, 
				 const float *restrict psfr, float dtheta, int ncompx, int ncompy,
				 const float *restrict srot, float alpha){
    float pxo=-(pixpsax*0.5-0.5+pixoffx)*pixtheta;
    float pyo=-(pixpsay*0.5-0.5+pixoffy)*pixtheta;
    int isa=blockIdx.x;
    float dispx=ncompx/2;
    float dispy=ncompy/2;
    float dtheta1=1.f/dtheta;
    float sx=0, cx=1.f;
    if(srot){
	sincosf(srot[isa], &sx, &cx);
    }
    ints+=isa*pixpsax*pixpsay;
    psfr+=isa*ncompx*ncompy;
    for(int iy=threadIdx.y; iy<pixpsay; iy+=blockDim.y){
	float y0=iy*pixtheta+pyo;
	for(int ix=threadIdx.x; ix<pixpsax; ix+=blockDim.x){
	    float x0=ix*pixtheta+pxo;
	    float x=(cx*x0-sx*y0)*dtheta1+dispx;
	    float y=(sx*x0+cx*y0)*dtheta1+dispy;
	    int jx=floorf(x); x=x-jx;
	    int jy=floorf(y); y=y-jy;
	    if(jx>=0 && jx<ncompx-1 && jy>=0 && jy<ncompy-1){
		ints[iy*pixpsax+ix]+=alpha*((psfr[jy*ncompx+jx]*(1-x)+psfr[jy*ncompx+jx+1]*x)*(1-y)
					    +(psfr[(jy+1)*ncompx+jx]*(1-x)+psfr[(jy+1)*ncompx+jx+1]*x)*y);
	    }
	}
    }
}
/**
   Add tip/tilt to the array. OPD=OPD+x*ttx+y*tty, where x=ix*dx+ox, y=iy*dy+oy; 
*/
__global__ static void add_tilt_do(float *opd, int nx, int ny, float ox, float oy, float dx, float ttx, float tty){
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	float vty=(oy+iy*dx)*tty;
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    opd[ix+iy*nx]+=vty+(ox+ix*dx)*ttx;
	}
    }
}
void wfsints(SIM_T *simu, float *phiout, int iwfs, int isim, cudaStream_t stream){
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
    const int hasllt=(parms->powfs[ipowfs].llt!=NULL);
    const int ncompx=powfs[ipowfs].ncompx;//necessary size to build detector image.
    const int ncompy=powfs[ipowfs].ncompy;
    const int nx=powfs[ipowfs].pts->nx;
    const int npsf=nx*parms->powfs[ipowfs].embfac;
    const int pixpsax=powfs[ipowfs].pixpsax;
    const int pixpsay=powfs[ipowfs].pixpsay;
    const float pixtheta=parms->powfs[ipowfs].pixtheta;
    const float siglev=parms->wfs[iwfs].siglevsim;
    const float *restrict const srot1=parms->powfs[ipowfs].radrot?cuwfs[iwfs].srot:NULL;
    const int multi_dtf=(parms->powfs[ipowfs].llt&&!parms->powfs[ipowfs].radrot &&parms->powfs[ipowfs].radpix);
    const float *restrict const srot2=multi_dtf?cuwfs[iwfs].srot:NULL;
    float *restrict const ints=cuwfs[iwfs].ints->p;
    float *lltopd=NULL;
    dcell *pistatout=NULL;
    if(parms->powfs[ipowfs].pistatout
       &&isim>=parms->powfs[ipowfs].pistatstart){
	if(!simu->pistatout[iwfs]){
	    simu->pistatout[iwfs]
		=dcellnew(nsa,parms->powfs[ipowfs].nwvl);
	}
	pistatout=simu->pistatout[iwfs];
    }
    ccell *psfout=NULL;
    cellarr *psfoutcellarr=NULL;
    cellarr *ztiltoutcellarr=NULL;
    if(parms->powfs[ipowfs].psfout){
	psfout=simu->wfspsfout[iwfs];
	psfoutcellarr=simu->save->wfspsfout[iwfs];
	ztiltoutcellarr=simu->save->ztiltout[iwfs];
    }
    if(powfs[ipowfs].llt && parms->powfs[ipowfs].trs){
	int nlx=powfs[ipowfs].llt->pts->nx;
	cudaMalloc(&lltopd, sizeof(float)*nlx*nlx);
	if(cuwfs[iwfs].lltncpa){
	    cudaMemcpyAsync(lltopd, cuwfs[iwfs].lltncpa, sizeof(float)*nlx*nlx, cudaMemcpyDefault, stream);
	}else{
	    cudaMemsetAsync(lltopd, 0, sizeof(float)*nlx*nlx, stream);
	}
	const int illt=parms->powfs[ipowfs].llt->i[wfsind];
	const double thetaxl=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox[illt]/hs;
	const double thetayl=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy[illt]/hs;
	gpu_atm2loc(lltopd, cupowfs[ipowfs].llt->loc, cupowfs[ipowfs].llt->nloc, 
		    hs, thetaxl, thetayl, 
		    parms->powfs[ipowfs].llt->misreg[0], 
		    parms->powfs[ipowfs].llt->misreg[1], 
		    dtisim, 1, stream);
	if((simu->uptreal && simu->uptreal->p[iwfs]) ||pistatout||parms->sim.uptideal){
	    float ttx,tty;
	    if(pistatout||parms->sim.uptideal){
		warning("Remove tip/tilt in uplink ideally\n");
		float *lltg;
		cudaCallocHost(lltg, 2*sizeof(float), stream);
		cuztilt<<<1,dim3(16,16),0,stream>>>(lltg, lltopd, 1, cupowfs[ipowfs].llt->dx, 
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
	    // copy uptreal to output 
	    PDMAT(simu->uptcmds->p[iwfs], puptcmds);
	    puptcmds[isim][0]=ttx;
	    puptcmds[isim][1]=tty;
	    // add tip/tilt to opd 
	    const double dx=powfs[ipowfs].llt->pts->dx;
	    const double ox=powfs[ipowfs].llt->pts->origx[0];
	    const double oy=powfs[ipowfs].llt->pts->origy[0];
	    add_tilt_do<<<1, dim3(16,16), 0, stream>>>(lltopd, nlx, nlx, ox, oy, dx, ttx, tty);
	}//if upt
	ctoc("llt opd");
    }//if has llt
    // Now begin physical optics 
    int rotpsfotf=(hasllt && parms->powfs[ipowfs].radrot);
    //CUDA_SYNC_STREAM;
    for(int iwvl=0; iwvl<parms->powfs[ipowfs].nwvl; iwvl++){
	float wvl=parms->powfs[ipowfs].wvl[iwvl];
	float dtheta=wvl/(npsf*powfs[ipowfs].pts->dx);
	fcomplex *lotfc=NULL;
	int nlx, nlpsf;
	if(lltopd){
	    //First calculate LOTF
	    nlx=powfs[ipowfs].llt->pts->nx;
	    nlpsf=nlx*parms->powfs[ipowfs].embfac;
	    cudaCalloc(lotfc, nlpsf*nlpsf*sizeof(fcomplex), stream);
	    embed_wvf_do<<<1,dim3(16,16),0,stream>>>(lotfc, lltopd, cuwfs[iwfs].lltamp, wvl, nlx, nlpsf);
	    //Turn to PSF
	    //CUDA_SYNC_DEVICE;
	    CUFFT(cuwfs[iwfs].lltplan1, lotfc, CUFFT_FORWARD);
	    if(nlpsf != npsf){
		//crop or embed the array.
		warning("Crop array from %d to %d\n", nlpsf, npsf);
		fcomplex *lotfc2=NULL;
		cudaCalloc(lotfc2, npsf*npsf*sizeof(fcomplex), stream);
		cpcorner_do<<<1, dim3(16,16),0,stream>>>(lotfc2, npsf, npsf, lotfc, nlpsf, nlpsf);
		//CUDA_SYNC_STREAM;
		cudaFree(lotfc);
		lotfc=lotfc2;
	    }
	    abs2real_do<<<1,dim3(16,16),0,stream>>>(lotfc, npsf, 1./(float)(nlpsf*nlpsf));
	    //Use backward to make lotfc the conjugate of otf. peak is in corner.
	    CUFFT(cuwfs[iwfs].lltplan2, lotfc, CUFFT_INVERSE);
	    //CUDA_SYNC_STREAM;
	}
	ctoc("llt otf");
	    
	float norm_psf=sqrt(powfs[ipowfs].areascale)/((float)powfs[ipowfs].pts->nx*npsf);
	float norm=siglev*norm_psf*norm_psf/((float)ncompx*ncompy);
	float norm_pistat=norm_psf*norm_psf/((float)npsf*npsf);

	/* Do msa subapertures in a batch to avoid using too much memory.*/
	fcomplex *psf, *otf;
	int msa=cuwfs[iwfs].msa;
	cudaCalloc(psf, sizeof(fcomplex)*npsf*npsf*msa, stream);
	if(srot1 || ncompx!=npsf || ncompy!=npsf){
	    cudaCalloc(otf, sizeof(fcomplex)*ncompx*ncompy*msa, stream);
	}else{
	    otf=psf;
	}
	for(int isa=0; isa<nsa; isa+=msa){
	    int ksa=MIN(msa, nsa-isa);//total number of subapertures left to do.
	    //embed amp/opd to wvf
	    cudaMemsetAsync(psf, 0, sizeof(fcomplex)*msa*npsf*npsf, stream);
	    if(otf!=psf) cudaMemsetAsync(otf, 0, sizeof(fcomplex)*msa*ncompx*ncompy, stream);
	    embed_wvf_do<<<ksa, dim3(16,16),0,stream>>>
		(psf, phiout+isa*nx*nx, cuwfs[iwfs].amp+isa*nx*nx, wvl, nx, npsf);
	    ctoc("embed");
	    //turn to complex psf, peak in corner
	    CUFFT(cuwfs[iwfs].plan1, psf, CUFFT_FORWARD);
	    ctoc("psf");
	    if(psfout) TO_IMPLEMENT;
	    //abs2 part to real, peak in corner
	    abs2real_do<<<ksa,dim3(16,16),0,stream>>>(psf, npsf, 1);
	    ctoc("abs2real");

	    if(lltopd || pistatout){
		//turn to otf.
		CUFFT(cuwfs[iwfs].plan1, psf, CUFFT_FORWARD);
		ctoc("fft to otf");
		if(lltopd){//multiply with uplink otf.
		    ccwm_do<<<ksa,dim3(16,16),0,stream>>>(psf, npsf, npsf, (fcomplex**)lotfc, 1);
		}
		ctoc("ccwm with lotfc");
		if(pistatout){
		    TO_IMPLEMENT;
		}
		//is OTF now.
	    }
	    ctoc("before ints");
	    //CUDA_SYNC_STREAM;
	    if(ints){
		if(srot1 || !(lltopd || pistatout) || otf!=psf){
		    //rotate PSF, or turn to OTF first time.
		    if(lltopd || pistatout){//srot1 is true. turn otf back to psf for rotation.
			//Was OTF. Turn to PSF
			norm/=((float)npsf*npsf);
			CUFFT(cuwfs[iwfs].plan1, psf, CUFFT_INVERSE);
			ctoc("fft to psf");
		    }
		    if(srot1){
			fftshift_do<<<ksa, dim3(16,16),0,stream>>>(psf, npsf, npsf);//shift to center
			embed_rot_do<<<ksa, dim3(16,16), 0, stream>>>
			    (otf, ncompx, ncompy, psf, npsf, npsf, srot1?srot1+isa:NULL);
			fftshift_do<<<ksa, dim3(16,16),0,stream>>>(otf, ncompx, ncompy);//shift back to corner
		    }else if(otf!=psf){
			cpcorner_do<<<ksa, dim3(16,16),0,stream>>>(otf, ncompx, ncompy, psf, npsf, npsf);
			ctoc("cpcorner");
		    }
		    //Turn PSF to OTF.
		    CUFFT(cuwfs[iwfs].plan2, otf,CUFFT_FORWARD);
		    ctoc("fft to otf");
		}
		//now we have otf. multiple with etf, dtf.
		if(cuwfs[iwfs].dtf[iwvl].etf){
		    if(cuwfs[iwfs].dtf[iwvl].etfis1d){
			ccwmcol_do<<<ksa,dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].etf+isa, 0);
		    }else{
			ctoc("before ccwm");
			ccwm_do<<<ksa,dim3(16,16),0,stream>>>
			    (otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].etf+isa, 0);
			ctoc("ccwm");
		    }
		}
		//multiple with nominal
		if(cuwfs[iwfs].dtf[iwvl].nominal){
		    ccwm_do<<<ksa,dim3(16,16),0,stream>>>
			(otf, ncompx, ncompy, cuwfs[iwfs].dtf[iwvl].nominal+isa, 0);
		    ctoc("nominal");
		}
		//back to spatial domain.
		CUFFT(cuwfs[iwfs].plan2, otf, CUFFT_INVERSE);
		ctoc("fft");
		float *psfr;
		DO(cudaMalloc(&psfr, sizeof(float)*ncompx*ncompy*ksa));
		realpart_do<<<ksa, dim3(16,16),0,stream>>>(psfr, otf, ncompx, ncompy);
		ctoc("realpart");
		si_rot_do<<<ksa, dim3(16,16),0,stream>>>
		    (ints+isa*pixpsax*pixpsay, pixpsax, pixpsay, 
		     parms->powfs[ipowfs].pixoffx, parms->powfs[ipowfs].pixoffy,
		     pixtheta, psfr, dtheta, ncompx, ncompy, srot2?srot2+isa:NULL, norm);
		ctoc("final");
		//CUDA_SYNC_STREAM;
		cudaFree(psfr);
	    }//if ints.
	}//for isa
	if(otf!=psf) cudaFree(otf);
	cudaFree(psf);
	if(lotfc) cudaFree(lotfc);
    }//for iwvl
    //CUDA_SYNC_STREAM;
    cudaFree(lltopd);
}
