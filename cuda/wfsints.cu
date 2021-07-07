/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#define TIMING 0

#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include "cucmat.h"
#include "wfs.h"
#include "cudata.h"
/*
  Timing results for TMT NFIRAOS case per LGS WFS:
  Embedding takes about 1 ms.
  FFT takes about 2 ms
  CCWM takes about 1ms
  realpart takes about 1ms

  Total takes about 12 ms.
*/


/**
   Embed amp*exp(2*pi*i*opd). input is nxin*nxin, output is nxout*nxout;
*/
__global__ static void sa_embed_wvf_do(Comp* restrict wvf,
	const Real* restrict opd, const Real* restrict amp,
	const Real wvl, const int nxin, const int nxout){
	const int isa=blockIdx.x;
	const int pad=(nxout-nxin)>>1;
	const int skipin=isa*nxin*nxin;
	const int skipout=isa*nxout*nxout+pad;
	const Real pi2l=2.f*M_PI/wvl;
	for(int iy=threadIdx.y; iy<nxin; iy+=blockDim.y){
		const int skipin2=skipin+iy*nxin;
		const int skipout2=skipout+(iy+pad)*nxout;
		for(int ix=threadIdx.x; ix<nxin; ix+=blockDim.x){
			/*test sinpi later */
			Real s, c;
			Z(sincos)(pi2l*opd[skipin2+ix], &s, &c);
			wvf[skipout2+ix].x=amp[skipin2+ix]*c;
			wvf[skipout2+ix].y=amp[skipin2+ix]*s;
		}
	}
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void sa_cpcorner_do(Comp* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny){
	int nx, ny;
	ny=MIN(niny, nouty)>>1;
	nx=MIN(ninx, noutx)>>1;
	const int isa=blockIdx.x;
	out+=isa*noutx*nouty;
	in+=isa*ninx*niny;
	for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
			out[iy*noutx+ix]=in[iy*ninx+ix];
			out[iy*noutx+(noutx-1-ix)]=in[iy*ninx+(ninx-1-ix)];
			out[(nouty-1-iy)*noutx+(noutx-1-ix)]=in[(niny-1-iy)*ninx+(ninx-1-ix)];
			out[(nouty-1-iy)*noutx+(ix)]=in[(niny-1-iy)*ninx+(ix)];
		}
	}
}

/**
   Embed or crop an array to another array. Preserve center.
*/
__global__ void sa_cpcenter_do(Comp* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny, Real scale){
	int nx, ny, nskipoutx, nskipouty, nskipinx, nskipiny;
	if(noutx<ninx){
		nx=noutx;
		nskipoutx=0;
		nskipinx=(ninx-noutx)>>1;
	} else{
		nx=ninx;
		nskipoutx=(noutx-ninx)>>1;
		nskipinx=0;
	}
	if(nouty<niny){
		ny=nouty;
		nskipouty=0;
		nskipiny=(niny-nouty)>>1;
	} else{
		ny=niny;
		nskipouty=(nouty-niny)>>1;
		nskipiny=0;
	}
	const int isa=blockIdx.x;
	out+=isa*noutx*nouty+nskipouty*(noutx)+nskipoutx;
	in+=isa*ninx*niny+nskipiny*(ninx)+nskipinx;
	for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
			out[iy*noutx+ix].x=scale*Z(cuCreal)(in[iy*ninx+ix]);
			out[iy*noutx+ix].y=scale*Z(cuCimag)(in[iy*ninx+ix]);
		}
	}
}
/**
   abs2 to real.
*/
__global__ static void sa_abs2real_do(Comp* wvf, const int nx, Real alpha){
	const int isa=blockIdx.x;
	wvf+=nx*nx*isa;
	for(int iy=threadIdx.y; iy<nx; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
			Real r=Z(cuCreal)(wvf[iy*nx+ix]);
			Real i=Z(cuCimag)(wvf[iy*nx+ix]);
			wvf[iy*nx+ix].x=(r*r+i*i)*alpha;
			wvf[iy*nx+ix].y=0;
		}
	}
}
/**
   FFT Shift from complex to real.
*/
__global__ static void sa_acc_real_fftshift_do(Real* restrict out, const Comp* restrict wvf,
	int nx, int ny, Real alpha){
	const int isa=blockIdx.x;
	wvf+=nx*ny*isa;
	out+=nx*ny*isa;
	int nx2=nx>>1;
	int ny2=ny>>1;
	for(int iy=threadIdx.y; iy<ny2; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<nx2; ix+=blockDim.x){
			out[ix+iy*nx]+=alpha*Z(cuCreal)(wvf[(ix+nx2)+(iy+ny2)*nx]);
			out[(ix+nx2)+(iy+ny2)*nx]+=alpha*Z(cuCreal)(wvf[ix+iy*nx]);
			out[ix+(iy+ny2)*nx]+=alpha*Z(cuCreal)(wvf[(ix+nx2)+iy*nx]);
			out[(ix+nx2)+iy*nx]+=alpha*Z(cuCreal)(wvf[ix+(iy+ny2)*nx]);
		}
	}
}
/**
   Multiple one array with another. OTFs contains repeated arrays each with notf
   numbers. if repeat is 1, dtfs is repeatd for each array.
*/
template <typename T, typename P>
__global__ static void sa_ccwm_do(T* otfs, const int notf, const P* dtfs, int repeat){
	const int isa=blockIdx.x;
	int offset=notf*isa;
	T* restrict otf=otfs+offset;
	const P* dtf=repeat?dtfs:(dtfs+offset);
	for(int ix=threadIdx.x; ix<notf; ix+=blockDim.x){
		otf[ix]*=dtf[ix];
	}
}
/**
   Multiple each OTF with another.
*/
__global__ static void sa_ccwm2_do(Comp* otfs, const int notfx, const int notfy,
	const Comp* dtfs1, Real wt1, const Comp* dtfs2, Real wt2, int repeat){
	const int isa=blockIdx.x;
	int offset=notfx*notfy*isa;
	Comp* restrict otf=otfs+offset;
	const Comp* dtf1=repeat?dtfs1:(dtfs1+offset);
	const Comp* dtf2=repeat?dtfs2:(dtfs2+offset);
	Comp temp;
	for(int iy=threadIdx.y; iy<notfy; iy+=blockDim.y){
		const int skip=iy*notfx;
		for(int ix=threadIdx.x; ix<notfx; ix+=blockDim.x){
			temp.x=Z(cuCreal)(dtf1[ix+skip])*wt1+Z(cuCreal)(dtf2[ix+skip])*wt2;
			temp.y=Z(cuCimag)(dtf1[ix+skip])*wt1+Z(cuCimag)(dtf2[ix+skip])*wt2;
			otf[ix+skip]*=temp;
		}
	}
}
/**
   Take the real part and accumulate to output
*/
__global__ static void sa_acc_real_do(Real* out, const Comp* restrict in, int ninx, int niny, Real alpha){
	const int isa=blockIdx.x;
	in+=isa*ninx*niny;
	out+=isa*ninx*niny;
	for(int iy=threadIdx.y; iy<niny; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<ninx; ix+=blockDim.x){
			out[ix+iy*ninx]+=Z(cuCreal)(in[ix+iy*ninx])*alpha;
		}
	}
}
/**
   Do the role of si. input psfr is sampled with notfx*notfy, with sampling dtheta.
   Output ints is sampled with pixpsax*pixpsay, at pixtheta.
*/
__global__ static void sa_si_rot_do(Real* restrict ints, int pixpsax, int pixpsay,
	const Real* pixoffx, const Real* pixoffy,
	Real pixthetax, Real pixthetay,
	const Comp* restrict otf, Real dtheta, int notfx, int notfy,
	const Real* restrict srot, Real alpha){
	int isa=blockIdx.x;
	Real pxo=-(pixpsax*0.5-0.5+(pixoffx?pixoffx[isa]:0))*pixthetax;
	Real pyo=-(pixpsay*0.5-0.5+(pixoffy?pixoffy[isa]:0))*pixthetay;

	Real dispx=notfx/2;
	Real dispy=notfy/2;
	Real dtheta1=1.f/dtheta;
	Real sx=0, cx=1.f;
	if(srot){
		Z(sincos)(srot[isa], &sx, &cx);
	}
	ints+=isa*pixpsax*pixpsay;
	otf+=isa*notfx*notfy;
	for(int iy=threadIdx.y; iy<pixpsay; iy+=blockDim.y){
		Real y0=iy*pixthetay+pyo;
		for(int ix=threadIdx.x; ix<pixpsax; ix+=blockDim.x){
			Real x0=ix*pixthetax+pxo;
			Real x=(cx*x0-sx*y0)*dtheta1+dispx;
			Real y=(sx*x0+cx*y0)*dtheta1+dispy;
			int jx=Z(floor)(x); x=x-jx;
			int jy=Z(floor)(y); y=y-jy;
			if(jx>=0&&jx<notfx-1&&jy>=0&&jy<notfy-1){
				ints[iy*pixpsax+ix]+=alpha*((Z(cuCreal)(otf[jy*notfx+jx])*(1.-x)
					+Z(cuCreal)(otf[jy*notfx+jx+1])*x)*(1.-y)
					+(Z(cuCreal)(otf[(jy+1)*notfx+jx])*(1.-x)
						+Z(cuCreal)(otf[(jy+1)*notfx+jx+1])*x)*y);
			}
		}
	}
}

/**
   Add tip/tilt to the OTF for each subaps. exp(-2*pi*sx/nx)*exp(-2*pi*sy/ny).
   peak of otf is in corner.
 */
__global__ static void sa_add_otf_tilt_corner_do(Comp* restrict otf, int nx, int ny,
	Real* restrict gx, Real* restrict gy, Real gscale){
	int isa=blockIdx.x;
	Real sx=gx[isa]*gscale;
	Real sy=gy[isa]*gscale;
	Comp* restrict otfi=otf+isa*nx*ny;
	for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
		for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
			Real phase=-2.f*M_PI*(+(Real)(ix-(ix*2>=nx?nx:0))/(Real)nx*sx
				+(Real)(iy-(iy*2>=ny?ny:0))/(Real)ny*sy);
			Real s, c;
			Z(sincos)(phase, &s, &c);
			Comp otfii=otfi[ix+iy*nx];
			otfi[ix+iy*nx].x=Z(cuCreal)(otfii)*c-Z(cuCimag)(otfii)*s;
			otfi[ix+iy*nx].y=Z(cuCreal)(otfii)*s+Z(cuCimag)(otfii)*c;
		}
	}
}
/**
   Do physical wfs images in GPU. please check wfsints() in CPU code for comments.

   stream is no longer an input parameter, as the FFT plan depends on it.
*/
void wfsints(sim_t* simu, Real* phiout, curmat& gradref, int iwfs, int isim){
	//TIC;tic;
	Array<cupowfs_t>& cupowfs=cudata->powfs;
	Array<cuwfs_t>& cuwfs=cuglobal->wfs;
	stream_t& stream=cuwfs[iwfs].stream;
	const parms_t* parms=simu->parms;
	const powfs_t* powfs=simu->powfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	const Real hs=parms->wfs[iwfs].hs;
	const Real hc=parms->wfs[iwfs].hc;
	const int nsa=powfs[ipowfs].saloc->nloc;
	const int notfx=powfs[ipowfs].notfx;/*necessary size to build detector image. */
	const int notfy=powfs[ipowfs].notfy;
	const int notf=MAX(notfx, notfy);
	const int nx=powfs[ipowfs].pts->nx;
	const int nwvf=nx*parms->powfs[ipowfs].embfac;
	const int pixpsax=powfs[ipowfs].pixpsax;
	const int pixpsay=powfs[ipowfs].pixpsay;
	const Real pixthetax=parms->powfs[ipowfs].radpixtheta;
	const Real pixthetay=parms->powfs[ipowfs].pixtheta;
	const Real siglev=parms->wfs[iwfs].sigsim;
	const Real* srot=(parms->powfs[ipowfs].llt&&parms->powfs[ipowfs].radpix)?cuwfs[iwfs].srot():NULL;
	const Real* pixoffx=0, * pixoffy=0;
	if(cupowfs[ipowfs].pixoffx){
		int icol=cupowfs[ipowfs].pixoffx>1?wfsind:0;
		pixoffx=cupowfs[ipowfs].pixoffx()+cupowfs[ipowfs].pixoffx.Nx()*icol;
		pixoffy=cupowfs[ipowfs].pixoffy()+cupowfs[ipowfs].pixoffy.Nx()*icol;
	}
	const int nwvl=parms->powfs[ipowfs].nwvl;
	curcell& ints=cuwfs[iwfs].ints;
	curcell pistatout;
	if(parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart){
		pistatout=cuwfs[iwfs].pistatout;
	}
	cuccell wvfout;
	const int wvf_n=notf/2+2;//was notf/2
	if(parms->powfs[ipowfs].psfout){
		wvfout=cuwfs[iwfs].wvfout;
	}
	Real norm_psf=sqrt(powfs[ipowfs].areascale)/((Real)powfs[ipowfs].pts->nx*nwvf);
	Real norm_pistat=norm_psf*norm_psf/((Real)notf*notf);
	Real norm_ints=siglev*norm_psf*norm_psf/((Real)notfx*notfy);
	/* Do msa subapertures in a batch to avoid using too much memory.*/
	Comp* psf, * wvf, * otf, * psfstat=NULL;

	Comp* lotfc=NULL;
	Comp* lwvf=NULL;
	curmat lltopd;

	if(powfs[ipowfs].llt&&parms->powfs[ipowfs].trs){
		int nlx=powfs[ipowfs].llt->pts->nx;
		lltopd=cuwfs[iwfs].lltopd;
		if(cuwfs[iwfs].lltncpa){
			cucp(lltopd, cuwfs[iwfs].lltncpa, stream);
		} else{
			cuzero(lltopd, stream);
		}
		const int illt=parms->powfs[ipowfs].llt->i->p[wfsind];
		const real thetaxl=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox->p[illt]/hs;
		const real thetayl=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy->p[illt]/hs;
		atm2loc(lltopd, cupowfs[ipowfs].llt.loc,
			hs, hc, thetaxl, thetayl,
			parms->powfs[ipowfs].llt->misreg->p[0],
			parms->powfs[ipowfs].llt->misreg->p[1],
			parms->sim.dt, isim, 1, stream);
		Real ttx=0, tty=0;
		if((simu->fsmreal&&simu->fsmreal->p[iwfs])||pistatout||parms->sim.idealfsm){
			if(pistatout||parms->sim.idealfsm){
			//warning("Remove tip/tilt in uplink ideally\n");
				Real* lltg=cuwfs[iwfs].lltg;
				lltg[0]=lltg[1]=0;
				cuztilt(lltg, lltopd, 1,
					cupowfs[ipowfs].llt.pts.Dxsa(),
					cupowfs[ipowfs].llt.pts.Nxsa(), cuwfs[iwfs].lltimcc,
					cupowfs[ipowfs].llt.pts(), cuwfs[iwfs].lltamp, 1.f, stream);
				CUDA_SYNC_STREAM;
				simu->fsmreal->p[iwfs]->p[0]=-lltg[0];
				simu->fsmreal->p[iwfs]->p[1]=-lltg[1];
			}
			ttx=simu->fsmreal->p[iwfs]->p[0];
			tty=simu->fsmreal->p[iwfs]->p[1];
		}/*if fsmreal */
		if(simu->telws){
			Real tt=simu->telws->p[isim]*parms->powfs[ipowfs].llt->ttrat;
			Real angle=simu->winddir?simu->winddir->p[0]:0;
			ttx+=tt*cosf(angle)*parms->powfs[ipowfs].llt->ttrat;
			tty+=tt*sinf(angle)*parms->powfs[ipowfs].llt->ttrat;
		}
		if(simu->llt_tt&&simu->llt_tt->p[iwfs]){
			ttx+=simu->llt_tt->p[iwfs]->p[isim];//put all to x direction.
		}
		if(ttx!=0&&tty!=0){
			/* add tip/tilt to opd  */
			const real dx=powfs[ipowfs].llt->pts->dx;
			const real ox=powfs[ipowfs].llt->pts->origx[0];
			const real oy=powfs[ipowfs].llt->pts->origy[0];
			add_tilt_do<<<1, dim3(16, 16), 0, stream>>>(lltopd, nlx, nlx, ox, oy, dx, ttx, tty);
		}
		int nlwvf=nlx*parms->powfs[ipowfs].embfac;
		lwvf=cuwfs[iwfs].lltwvf;
		if(nlwvf!=notf){
			lotfc=cuwfs[iwfs].lltotfc;
		} else{
			lotfc=lwvf;
		}
		if(parms->save.wfsopd->p[iwfs]){
			zfarr_push(simu->save->wfslltopd[iwfs], isim, lltopd, stream);
		}
	}/*if has llt */

	/* Now begin physical optics preparation*/
	int isotf=(lltopd||pistatout);
	int msa=cuwfs[iwfs].msa;/* number of subaps to process at each time.*/
	wvf=cuwfs[iwfs].wvf;
	if(nwvf==notf){
		psf=wvf;
	} else{
		psf=cuwfs[iwfs].psf;
	}
	if(notfx!=notf||notfy!=notf){
		otf=cuwfs[iwfs].otf;
		if(isotf){/*There is an additional pair of FFT.*/
			norm_ints/=((Real)notf*notf);
		}
	} else{
		otf=psf;
	}

	if(pistatout){
		psfstat=cuwfs[iwfs].psfstat();
	}
	/* Now begin physical optics  */
	for(int iwvl=0; iwvl<nwvl; iwvl++){

		Real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
		Real dtheta=wvl/(nwvf*powfs[ipowfs].pts->dx);
		if(lltopd){ /*First calculate LOTF */
			int nlx=powfs[ipowfs].llt->pts->nx;
			int nlwvf=nlx*parms->powfs[ipowfs].embfac;
			cudaMemsetAsync(lwvf, 0, sizeof(Comp)*nlwvf*nlwvf, stream);
			if(lotfc!=lwvf){
				cudaMemsetAsync(lotfc, 0, sizeof(Comp)*notf*notf, stream);
			}
			sa_embed_wvf_do<<<1, dim3(16, 16), 0, stream>>>
				(lwvf, lltopd, cuwfs[iwfs].lltamp, wvl, nlx, nlwvf);
				/*Turn to PSF. peak in corner */
			CUFFT(cuwfs[iwfs].lltplan_wvf, lwvf, CUFFT_FORWARD);
			sa_abs2real_do<<<1, dim3(16, 16), 0, stream>>>(lwvf, nlwvf, 1./(Real)(nlwvf*nlwvf));
			/*Turn to OTF. peak in corner*/
			/*Use backward to make lotfc the conjugate of otf. peak is in corner. */
			CUFFT(cuwfs[iwfs].lltplan_wvf, lwvf, CUFFT_INVERSE);
			if(lwvf!=lotfc){
				sa_cpcorner_do<<<1, dim3(16, 16), 0, stream>>>(lotfc, notf, notf, lwvf, nlwvf, nlwvf);
			}
		}

		for(int isa=0; isa<nsa; isa+=msa){
			ctoc_init(30);
			int ksa=MIN(msa, nsa-isa);/*total number of subapertures left to do. */
			/*embed amp/opd to wvf */
			cudaMemsetAsync(wvf, 0, sizeof(Comp)*ksa*nwvf*nwvf, stream);
			if(psf!=wvf){
				cudaMemsetAsync(psf, 0, sizeof(Comp)*ksa*notf*notf, stream);
			}
			sa_embed_wvf_do<<<ksa, dim3(16, 16), 0, stream>>>
				(wvf, phiout+isa*nx*nx, cuwfs[iwfs].amp()+isa*nx*nx, wvl, nx, nwvf);
			ctoc("embed");//0.7 ms
			/* turn to complex psf, peak in corner */
			CUFFT(cuwfs[iwfs].plan1, wvf, CUFFT_FORWARD);
			/* copy big psf to smaller psf covering detector focal plane. */
			if(psf!=wvf){
				sa_cpcorner_do<<<ksa, dim3(16, 16), 0, stream>>>
					(psf, notf, notf, wvf, nwvf, nwvf);
			}
			//gpu_write(psf, notf, notf*ksa, "psf_out_1");
			ctoc("psf");//1.1ms
			if(wvfout){
				cuzero(cuwfs[iwfs].psfout, stream);
				CUFFT2(cuwfs[iwfs].plan2, psf, cuwfs[iwfs].psfout, CUFFT_INVERSE);
				sa_cpcenter_do<<<ksa, dim3(16, 16), 0, stream>>>
					(wvfout[isa+nsa*iwvl], wvf_n, wvf_n,
						cuwfs[iwfs].psfout, notf, notf, norm_psf/(notf*notf));
			}
			/* abs2 part to real, peak in corner */
			sa_abs2real_do<<<ksa, dim3(16, 16), 0, stream>>>(psf, notf, 1);
			ctoc("abs2real");//0.8ms
			//gpu_write(psf, notf, notf*ksa, "psf_out_2");
			if(isotf){
			/* turn to otf. peak in corner */
				CUFFT(cuwfs[iwfs].plan2, psf, CUFFT_FORWARD);
				ctoc("fft to otf");//3.7 ms for notf=84. 2.6 ms for 70. 0.7ms for 64.
				if(pistatout){
					cudaMemcpyAsync(psfstat, psf, sizeof(Comp)*notf*notf*ksa,
						MEMCPY_D2D, stream);
					if(parms->powfs[ipowfs].pistatout){
						sa_add_otf_tilt_corner_do<<<ksa, dim3(16, 16), 0, stream>>>
							(psfstat, notf, notf, gradref()+isa, gradref()+nsa+isa, -1.f/dtheta);
					}
					CUFFT(cuwfs[iwfs].plan2, psfstat, CUFFT_INVERSE);/*back to PSF. peak in corner*/
					if(parms->sim.skysim){/*want peak in corner*/
						sa_acc_real_do<<<ksa, dim3(16, 16), 0, stream>>>
							(pistatout[isa+nsa*iwvl], psfstat, notf, notf, norm_pistat);
					} else{/*want peak in center*/
						sa_acc_real_fftshift_do<<<ksa, dim3(16, 16), 0, stream>>>
							(pistatout[isa+nsa*iwvl], psfstat, notf, notf, norm_pistat);
					}
				}
				if(lltopd){/*multiply with uplink otf. */
					sa_ccwm_do<<<ksa, 256, 0, stream>>>(psf, notf*notf, lotfc, 1);
					ctoc("ccwm with lotfc");//0.66 ms
				}
				/* is OTF now. */
			}
			//gpu_write(psf, notf, notf*ksa, "psf_out_3");
			if(ints){
				if(!isotf||otf!=psf){/*rotate PSF, or turn to OTF first time. */
					if(isotf){/*srot1 is true. turn otf back to psf for rotation. */
						CUFFT(cuwfs[iwfs].plan2, psf, CUFFT_INVERSE);
						ctoc("fft to psf");
					}
					if(otf!=psf){
						cudaMemsetAsync(otf, 0, sizeof(Comp)*ksa*notfx*notfy, stream);
					}
					if(otf!=psf){/*copy the psf corner*/
						sa_cpcorner_do<<<ksa, dim3(16, 16), 0, stream>>>
							(otf, notfx, notfy, psf, notf, notf);
						ctoc("cpcorner");
					}
					/*Turn PSF to OTF. */
					CUFFT(cuwfs[iwfs].plan3, otf, CUFFT_FORWARD);
					ctoc("fft to otf");
				}
				/*now we have otf. multiply with etf, dtf. */
				if(cuwfs[iwfs].dtf[iwvl].etf[1].etf){
					const int dtrat=parms->powfs[ipowfs].llt->coldtrat;
					Real wt2=(Real)(isim%dtrat)/(real)dtrat;
					Real wt1=1.-wt2;
					sa_ccwm2_do<<<ksa, dim3(16, 16), 0, stream>>>
						(otf, notfx, notfy, cuwfs[iwfs].dtf[iwvl].etf[0].etf.Col(isa), wt1, cuwfs[iwfs].dtf[iwvl].etf[1].etf.Col(isa), wt2, 0);
					ctoc("ccwm2");
				} else if(cuwfs[iwfs].dtf[iwvl].etf[0].etf){
					sa_ccwm_do<<<ksa, 256, 0, stream>>>
						(otf, notfx*notfy, cuwfs[iwfs].dtf[iwvl].etf[0].etf.Col(isa), 0);
					ctoc("ccwm");//0.97 ms
				}
				/*multiply with nominal */
				if(cuwfs[iwfs].dtf[iwvl].nominal){
					int repeat=0;
					Comp* pnominal=0;
					if(cuwfs[iwfs].dtf[iwvl].nominal.Ny()==1){
						repeat=1;
						pnominal=cuwfs[iwfs].dtf[iwvl].nominal.Col(0);
					} else{
						repeat=0;
						pnominal=cuwfs[iwfs].dtf[iwvl].nominal.Col(isa);
					}
					sa_ccwm_do<<<ksa, 256, 0, stream>>>
						(otf, notfx*notfy, pnominal, repeat);
					ctoc("nominal");
				}
				/*back to spatial domain. */
				CUFFT(cuwfs[iwfs].plan3, otf, CUFFT_INVERSE);
				ctoc("fft");//3.7 ms
				sa_si_rot_do<<<ksa, dim3(16, 16), 0, stream>>>
					(ints[isa], pixpsax, pixpsay, pixoffx?pixoffx+isa:NULL, pixoffy?pixoffy+isa:NULL,
						pixthetax, pixthetay, otf, dtheta, notfx, notfy, srot?srot+isa:NULL,
						norm_ints*parms->wfs[iwfs].wvlwts->p[iwvl]);
				ctoc("si");//0.1 ms
				ctoc_final("wfsints %d", iwfs);
			}/*if ints. */
		}/*for isa block loop */
	}/*for iwvl */
	if(parms->powfs[ipowfs].psfout){
		zfarr_push(simu->save->wfspsfout[iwfs], isim, wvfout, stream);
	}
}
