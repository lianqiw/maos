/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "../math/cumath.h"
#include "accphi.h"
#include <curand_kernel.h>
#include "wfs.h"
#include "cudata.h"

#if !USE_CPP
extern "C"{
#endif
#include "../../maos/sim.h"
#if !USE_CPP
}
#endif
#define triscalex 0.866025403784439 //sqrt(3.)/2;
#define PYWFS4_GRAD_CALC\
    grad[i]=(ints[i]-ints[i+nsa]+ints[nsa*2+i]-ints[nsa*3+i])*alpha-goff[i]; \
    grad[i+nsa]=(ints[i]+ints[i+nsa]-ints[nsa*2+i]-ints[nsa*3+i])*alpha-goff[i+nsa]

#define PYWFS3_GRAD_CALC\
    grad[i]=(ints[i+nsa]-ints[nsa*2+i])*alpha*triscalex-goff[i]; \
    grad[i+nsa]=(ints[i]-0.5*(ints[i+nsa]+ints[nsa*2+i]))*alpha-goff[i+nsa]

#define PYWFSR_GRAD_CALC\
	for(int ig=0; ig<ng; ig++){\
		grad[nsa*ig+i]=ints[nsa*ig+i]*alpha-goff[nsa*ig+i];\
	}
//4 side pyramid
__global__ static void
pywfs4_grad_0_do(Real* grad, const Real* ints, const Real* saa, const Real siglev, const Real* goff, Real gain, int nsa){
	const Real alpha0=gain/siglev;
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=alpha0/saa[i];
		PYWFS4_GRAD_CALC;
	}
}
__global__ static void
pywfs4_grad_1_do(Real* grad, const Real* ints, const Real* saa, const Real* goff, Real gain, int nsa){
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=gain/(ints[i]+ints[i+nsa]+ints[nsa*2+i]+ints[nsa*3+i]);
		PYWFS4_GRAD_CALC;
	}
}
__global__ static void
pywfs4_grad_2_do(Real* grad, const Real* ints, const Real* saa, const Real* isum, const Real* goff, Real gain, int nsa){
	const Real alpha0=gain*nsa/(*isum);
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=alpha0/saa[i];
		PYWFS4_GRAD_CALC;
	}
}
//3 side pyramid
__global__ static void
pywfs3_grad_0_do(Real* grad, const Real* ints, const Real* saa, const Real siglev, const Real* goff, Real gain, int nsa){
	const Real alpha0=gain/siglev;
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=alpha0/saa[i];
		PYWFS3_GRAD_CALC;
	}
}
__global__ static void
pywfs3_grad_1_do(Real* grad, const Real* ints, const Real* saa, const Real* goff, Real gain, int nsa){
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=gain/(ints[i]+ints[i+nsa]+ints[nsa*2+i]);
		PYWFS3_GRAD_CALC;
	}
}
__global__ static void
pywfs3_grad_2_do(Real* grad, const Real* ints, const Real* saa, const Real* isum, const Real* goff, Real gain, int nsa){
	const Real alpha0=gain*nsa/(*isum);
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		const Real alpha=alpha0/saa[i];
		PYWFS3_GRAD_CALC;
	}
}
//RAW: use normalized ints as grads
__global__ static void
pywfsr_grad_do(Real *grad, const Real *ints, const Real *saa, const Real siglev, const Real *goff, Real gain, int nsa,int ng){
	const Real alpha=gain/siglev;
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<nsa; i+=blockDim.x*gridDim.x){
		PYWFSR_GRAD_CALC;
	}
}

void pywfs_grad(curmat& grad, /**<[out] gradients*/
		const curmat& ints, /**<[in] Intensity*/
		const curmat& saa,  /**<[in] Subaperture normalized area*/
		curmat& isum, /**<[out] Sum intensity*/
		const curmat& goff, /**<[in] Gradient of flat wavefront*/
		const pywfs_t* pywfs,
		cudaStream_t stream){
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	const int ng=pywfs_ng(pycfg);
	if(pycfg->raw || pycfg->nside<3){
		pywfsr_grad_do<<<DIM(ints.Nx(), 256), 0, stream>>>
			(grad, ints, saa, pycfg->siglev, goff, pywfs->gain, ints.Nx(),ng);
	}else switch(pycfg->sigmatch){
	case 0://No siglev correction
		info_once("No siglev correction\n");
		if(pycfg->nside==3){
			pywfs3_grad_0_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, pycfg->siglev, goff, pywfs->gain, ints.Nx());
		}else if(pycfg->nside==4){
			pywfs4_grad_0_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, pycfg->siglev, goff, pywfs->gain, ints.Nx());
		}
		break;
	case 1:
		info_once("Individual correction\n");
		if(pycfg->nside==3){
			pywfs3_grad_1_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, goff, pywfs->gain, ints.Nx());
		}else if(pycfg->nside==4){
			pywfs4_grad_1_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, goff, pywfs->gain, ints.Nx());
		}
		break;
	case 2:
		info_once("Global correction (preferred);\n");
		cursum2(isum, ints, stream);//sum of ints
		if(pycfg->nside==3){
			pywfs3_grad_2_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, isum, goff, pywfs->gain, ints.Nx());
		}else if(pycfg->nside==4){
			pywfs4_grad_2_do<<<DIM(ints.Nx(), 256), 0, stream>>>
				(grad, ints, saa, isum, goff, pywfs->gain, ints.Nx());
		}
		break;
	default:
		error("Invalid sigmatch.\n");
	}
}
/**
   FFT for PYWFS.

   stream is no longer an input parameter, as the FFT plan depends on it.

   2017-06-22: Tried to parallelize the modulation step to multiple streams
   within a GPU. Does not help since each FFT occupies the full GPU, preventing
   stream concurrency.
 */
void pywfs_ints(curmat& ints, curmat& phiout, cuwfs_t& cuwfs, Real siglev){
	//Pyramid WFS
	cupowfs_t* cupowfs=cuwfs.powfs;
	stream_t& stream=cuwfs.stream;
	pywfs_t* pywfs=cupowfs->pywfs;
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	cuzero(cuwfs.pypsf, stream);
	locfft_t* locfft=pywfs->locfft;
	const int nwvl=locfft->wvl->nx;
	const Real dx=locfft->loc->dx;
	const long nembed=locfft->nembed->p[0];
	const long nembed2=nembed/2;
	const long ncomp=pywfs->nominal->nx;
	const long ncomp2=ncomp/2;
	const Real pos_r=pycfg->modulate;
	const int pos_n=pycfg->modulpos;
	const int pos_nr=pycfg->modulring;
	cucmat& otf=cuwfs.pyotf;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		ctoc_init(20);
		cucmat& wvf=cuwfs.pywvf[iwvl];
		const Real wvl=locfft->wvl->p[iwvl];
		cuzero(wvf, stream);
		embed_wvf_do<<<DIM(phiout.Nx(), 256), 0, stream>>>
			(wvf, phiout, cuwfs.amp, cupowfs->embed[iwvl], phiout.Nx(), wvl);
		ctoc("embed");
		CUFFT(cuwfs.plan_fs, wvf, CUFFT_FORWARD);
		ctoc("fft");
		fftshift_do<<<DIM2(wvf.Nx(), wvf.Ny(), 16), 0, stream>>>
			(wvf, wvf.Nx(), wvf.Ny());
		ctoc("shift");
		const Real otfnorm=1./(sqrt(locfft->ampnorm)*locfft->nembed->p[iwvl]);
		Scale(wvf, otfnorm, stream);
		ctoc("scale");
		if(global->setupdone && global->parms->plot.run){
			cucdraw_gpu("Ints", wvf, 1, stream, 1, "PWFS PSF", "x", "y", "wfs %d focus", pywfs->iwfs0);
		}
		//cuwrite(wvf, stream, "gpu_wvf0");
		const Real dtheta=locfft->wvl->p[iwvl]/(dx*nembed);
		for(int ir=0; ir<pos_nr; ir++){
			real pos_ri=pos_r*(ir+1)/pos_nr;
			//Scale number of points by ring size to have even surface brightness
			int pos_ni=pos_n*(ir+1)/pos_nr;
			const Real alpha=pywfs->cfg->wvlwts->p[iwvl]/(ncomp*ncomp*pos_ni*pos_nr);
			for(int ipos=0; ipos<pos_ni; ipos++){
				const Real theta=2*M_PI*ipos/pos_ni;
				const Real posx=cos(theta)*pos_ri;
				const Real posy=sin(theta)*pos_ri;
				const long offy=(long)round(posy/dtheta);//offset of center
				const long offy2=nembed2+offy-ncomp2;//offset of corner
				const long iy0=MAX(-offy2, 0);
				const long ny2=MIN(ncomp, nembed-offy2)-iy0;

				const long offx=(long)round(posx/dtheta);
				const long offx2=nembed/2+offx-ncomp2;
				const long ix0=MAX(-offx2, 0);
				const long nx2=MIN(ncomp, nembed-offx2)-ix0;

				cuzero(otf, stream);
				cwm_do<<<DIM2(nx2, ny2, 16), 0, stream>>>
					(otf()+ix0+iy0*ncomp,
						cupowfs->pyramid[iwvl]()+ix0+iy0*ncomp,
						wvf()+ix0+offx2+(iy0+offy2)*nembed,
						ncomp, nembed, nx2, ny2);
				   //cuwrite(otf, stream, "gpu_wvf1_%d", ipos);
				CUFFT(cuwfs.plan_py, otf, CUFFT_INVERSE);
				//cuwrite(otf, stream, "gpu_wvf2_%d", ipos);
				curaddcabs2(cuwfs.pypsf, otf, alpha, stream);
				//cuwrite(cuwfs.pypsf, stream, "gpu_wvf3_%d", ipos);
			}//iposr: points along a ring.
		}//ir: ring of dithering
		if(global->setupdone&&global->parms->plot.run){
			curdraw_gpu("Ints", cuwfs.pypsf, 1, stream, 1, "PWFS Pupil", "x", "y", "wfs %d pupil", pywfs->iwfs0);
		}
		ctoc("modul");
		embed_do<<<DIM(ncomp*ncomp, 256), 0, stream>>>
			(otf, cuwfs.pypsf, ncomp*ncomp);
		ctoc("embed");
		//cuwrite(otf, stream, "gpu_wvf4");
		CUFFT(cuwfs.plan_py, otf, CUFFT_FORWARD);//otf contains the OTF
		ctoc("fft");
		//cuwrite(otf, stream, "gpu_wvf5");
		cwm_do<<<DIM(ncomp*ncomp, 256), 0, stream>>>
			(otf(), cupowfs->pynominal(), ncomp*ncomp);//applying DTF and then FFT
		ctoc("cwm");
		CUFFT(cuwfs.plan_py, otf, CUFFT_INVERSE);
		ctoc("ifft");
		//cuwrite(otf, stream, "gpu_wvf6");
		//Use ray tracing for si
		const Real dx2=dx*nembed/ncomp;
		const int nsa=cupowfs->saloc.Nloc();
		const Real alpha=(Real)nsa*siglev/(Real)(ncomp*ncomp);
		for(int ind=0; ind<pycfg->nside; ind++){
			Real shx=0, shy=0;
			culoc_t saloc;
			if(cupowfs->msaloc){
				saloc=cupowfs->msaloc[ind];
			} else{
				saloc=cupowfs->saloc;
			}
			if(pywfs->pupilshift){
				shx=P(pywfs->pupilshift, ind, 0)*saloc.Dx();
				shy=P(pywfs->pupilshift, ind, 1)*saloc.Dy();
			}
			Real* pout=ints()+ind*nsa;
			Real offx=P(pywfs->sioff, ind, 0);
			Real offy=P(pywfs->sioff, ind, 1);
			map2loc_linear<<<DIM(nsa, 256), 0, stream>>>
				(pout, otf, otf.Nx(), otf.Ny(),
					saloc(), saloc.Nloc(),
					1./dx2, 1./dx2,
					(offx*ncomp2)-(-ncomp2+0.5)+shx,
					(offy*ncomp2)-(-ncomp2+0.5)+shy, alpha);
		}
		ctoc("sample");
		ctoc_final("pywfs");
		//cuwrite(ints, stream, "gpu_ints"); exit(0);
	}
}
dmat* gpu_pywfs_mkg(const pywfs_t* pywfs, const loc_t* locin, const loc_t* locfft, const dmat* mod, real displacex, real displacey){
	gpu_set(cuglobal->wfsgpu[pywfs->iwfs0]);
	const pywfs_cfg_t *pycfg=pywfs->cfg;
	cuwfs_t& cuwfs=cuglobal->wfs[pywfs->iwfs0];
	cupowfs_t* cupowfs=cuwfs.powfs;
	stream_t& stream=cuwfs.stream;
	mapcell* mapinsq=(mapcell*)cellnew(1, 1);
	dcell* opdin=dcellnew(1, 1);
	const Real siglev=100;//irrelevant in noise free case.
	for(int i=0; i<1; i++){
		mapinsq->p[i]=mapnew2(locin->map);
		opdin->p[i]=dnew(locin->nloc, 1);
	}
	cumapcell cumapin(1, 1);
	cp2gpu(cumapin, mapinsq);
	curmat phiout(locfft->nloc, 1);
	curmat phiout0(locfft->nloc, 1);
	if(pywfs->opdadd){
		cp2gpu(phiout0, pywfs->opdadd);
	}
	culoc_t culocout(locfft);
	const int nsa=cupowfs->saloc.Nloc();
	curmat ints(nsa, pycfg->nside);
	const int ng=pywfs_ng(pycfg);
	curmat grad(nsa*ng, 1);
	curmat grad0(nsa*ng, 1);
	cuzero(ints, stream);
	Add(phiout, (Real)1, phiout0, (Real)1, stream);
	pywfs_ints(ints, phiout, cuwfs, siglev);
	pywfs_grad(grad0, ints, cupowfs->saa(0), cuwfs.isum, cupowfs->pyoff, pywfs, stream);
	TIC;tic;
	//cuwrite(grad0, stream, "grad0_gpu");
	//cuwrite(ints, stream, "ints0_gpu");
	const int nmod=mod?mod->ny:locin->nloc;
	dmat* ggd=dnew(nsa*ng, nmod);

	for(int imod=0; imod<nmod; imod++){
		Real poke=pycfg->poke;
		if(mod){
			dmat* tmp=drefcols(mod, imod, 1);
			//real radial=ceil((sqrt(8.*(imod+1)+1)-3)*0.5)+1;
			real tmax, tmin;
			dmaxmin(tmp, &tmax, &tmin);
			poke/=(tmax-tmin);//sqrt(radial);
			dadd(&P(opdin, 0), 0, tmp, poke);
			dfree(tmp);
		} else{
			opdin->p[0]->p[imod]=poke;
			if(imod>0){
				opdin->p[0]->p[imod-1]=0;
			}
		}
		loc_embed(P(mapinsq,0), locin, P(opdin,0));
		CUDA_SYNC_STREAM;
		cp2gpu(cumapin, mapinsq);
		//cuzero(phiout, stream);
		Copy(phiout, phiout0, stream);
		mapcell2loc(phiout, culocout, cumapin, pywfs->cfg->hs, pywfs->cfg->hc, displacex, displacey, 0, 0, 1, stream);
		//cuwrite(cumapin[0].p, stream, "gpu_cumapin_%d", imod);
		//cuwrite(phiout, stream, "gpu_phiout_%d", imod);
		cuzero(ints, stream);
		pywfs_ints(ints, phiout, cuwfs, siglev);
		//cuwrite(ints, stream, "gpu_ints_%d", imod);
		pywfs_grad(grad, ints, cupowfs->saa(0), cuwfs.isum, cupowfs->pyoff, pywfs, stream);
		Add(grad, (Real)1, grad0, (Real)-1, stream);
		Scale(grad, (Real)1./poke, stream);
		dmat* gradc=drefcols(ggd, imod, 1);
		cp2cpu(&gradc, grad, stream);
		if(imod%((nmod+9)/10)==0){
			Real ts=myclockd()-tk;
			info("%d of %ld. %.2f of %.2f seconds. std(grad)=%g.\n", imod, locin->nloc, ts, ts/(imod+1)*locin->nloc, dstd(gradc));
		}
		dfree(gradc);
	}
	cellfree(opdin);
	cellfree(mapinsq);
	return ggd;
}
