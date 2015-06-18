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
#include "accphi.h"
#include <curand_kernel.h>
#include "wfs.h"
#include "cudata.h"
#include "cucmat.h"
#ifdef __cplusplus
extern "C"{
#endif
#include "../maos/pywfs.h"
#include "../maos/sim.h"
#ifdef __cplusplus
}
#endif

__global__ static void
pywfs_grad_do(Real *grad, Real *ints, Real *saa, Real *isum, Real *goff, Real gain, int nsa){
    Real alpha0=gain*nsa/(*isum);
    for(int i=threadIdx.x + blockIdx.x * blockDim.x; i<nsa; i+=blockDim.x * gridDim.x){
	Real alpha=alpha0/saa[i];
	grad[i]=(-ints[i]+ints[i+nsa]-ints[nsa*2+i]+ints[nsa*3+i])*alpha-goff[i];
	grad[i+nsa]=(-ints[i]-ints[i+nsa]+ints[nsa*2+i]+ints[nsa*3+i])*alpha-goff[i+nsa];
    }
}
void pywfs_grad(curmat *grad, /**<[out] gradients*/
		const curmat *ints, /**<[in] Intensity*/
		const curmat *saa,  /**<[in] Subaperture normalized area*/
		curmat *isum, /**<[out] Sum intensity*/
		const curmat *goff, /**<[in] Gradient of flat wavefront*/
		Real gain,   /**<[in] Gain*/
		cudaStream_t stream){
    cursum2(isum->p, ints, stream);//sum of ints
    pywfs_grad_do<<<DIM(ints->nx, 256), 0, stream>>>
	(grad->p, ints->p, saa->p, isum->p, goff->p, gain, ints->nx);
}
void pywfs_ints(curmat *ints, curmat *phiout, cuwfs_t *cuwfs, Real siglev, cudaStream_t stream){
    //Pyramid WFS
    cupowfs_t *cupowfs=cuwfs->powfs;
    PYWFS_T *pywfs=cupowfs->pywfs;
    cuzero(cuwfs->pypsf, stream);
    locfft_t *locfft=pywfs->locfft;
    const int nwvl=locfft->wvl->nx;
    Real pos_r=pywfs->modulate; 
    Real dx=locfft->loc->dx;
    long nembed=locfft->nembed->p[0];
    long nembed2=nembed/2;
    long ncomp=pywfs->nominal->nx;
    long ncomp2=ncomp/2;
    int pos_n=pywfs->modulpos;
    if(pos_r<=0){
	pos_n=1;
    }
    cucmat *otf=cuwfs->pyotf;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cuwfs->pywvf->p[iwvl];
	Real alpha=pywfs->wvlwts->p[iwvl]/(ncomp*ncomp*pos_n);
	Real wvl=locfft->wvl->p[iwvl];
	cuzero(wvf, stream);
	embed_wvf_do<<<DIM(phiout->nx,256),0,stream>>>
	    (wvf->p, phiout->p, cuwfs->amp, cupowfs->embed[iwvl], phiout->nx, wvl);
	CUFFT(cuwfs->plan_fs, wvf->p, CUFFT_FORWARD);
	fftshift_do<<<DIM2(wvf->nx, wvf->nx, 16),0,stream>>>
	    (wvf->p, wvf->nx, wvf->ny);
	Real otfnorm=1./(sqrt(locfft->ampnorm)*locfft->nembed->p[iwvl]);
	cucscale(wvf, otfnorm, stream);
	Real dtheta=locfft->wvl->p[iwvl]/(dx*nembed);
	for(int ipos=0; ipos<pos_n; ipos++){
	    Real theta=2*M_PI*ipos/pos_n;
	    Real posx=cos(theta)*pos_r;
	    Real posy=sin(theta)*pos_r;
	    long offy=(long)round(posy/dtheta);
	    long offy2=nembed2+offy-ncomp2;
	    long iy0=MAX(-offy2, 0);
	    long ny2=MIN(ncomp+offy2, nembed)-offy2-iy0;

	    long offx=(long)round(posx/dtheta);
	    long offx2=nembed/2+offx-ncomp2;
	    long ix0=MAX(-offx2, 0);
	    long nx2=MIN(ncomp+offx2, nembed)-offx2-ix0;
	    cuzero(otf, stream);
	    cwm_do<<<DIM2(nx2, ny2,16),0,stream>>>
		(otf->p+ix0+iy0*ncomp, 
		 cupowfs->pyramid->p[iwvl]->p+ix0+iy0*ncomp, 
		 wvf->p+ix0+offx2+(iy0+offy2)*nembed,
		 ncomp, nembed, nx2, ny2);
	    CUFFT(cuwfs->plan_py, otf->p, CUFFT_INVERSE);
	    curaddcabs2(&cuwfs->pypsf, 1, otf, alpha, stream);
	}
	embed_do<<<DIM(ncomp*ncomp, 256),0,stream>>>
	    (otf->p, cuwfs->pypsf->p, ncomp*ncomp);
	CUFFT(cuwfs->plan_py, otf->p, CUFFT_FORWARD);
	cwm_do<<<DIM(ncomp*ncomp, 256), 0, stream>>>
	    (otf->p, cupowfs->pynominal->p, ncomp*ncomp);
	CUFFT(cuwfs->plan_py, otf->p, CUFFT_INVERSE);
	//Use ray tracing for si
	Real dx2=dx*nembed/ncomp;
	const int nsa=cupowfs->saloc->nloc;
	for(int iy=0; iy<2; iy++){
	    for(int ix=0; ix<2; ix++){
		Real* pout=ints->p+(ix+iy*2)*nsa;
		prop_linear<<<DIM(nsa, 256), 0, stream>>>
		    (pout, otf->p, otf->nx, otf->ny, 
		     cupowfs->saloc->p, cupowfs->saloc->nloc,
		     1./dx2, 1./dx2, 
		     (((ix-0.5)*ncomp2)-(-ncomp2+0.5)),
		     (((iy-0.5)*ncomp2)-(-ncomp2+0.5)), 
		     (Real)nsa/(ncomp*ncomp)*siglev);
	    }
	}
    }
}
dsp *gpu_pywfs_mkg(const PARMS_T *parms, const POWFS_T *powfs, loc_t *aloc, int iwfs, int idm){
    gpu_set(cudata_t::wfsgpu[iwfs]);
    cuwfs_t *cuwfs=cudata_t::wfs+iwfs;
    cupowfs_t *cupowfs=cuwfs->powfs;
    PYWFS_T *pywfs=cupowfs->pywfs;
    stream_t &stream=*cuwfs->stream;
    mapcell *dmrealsq=(mapcell*)cellnew(parms->ndm,1);
    dcell *dmreal=dcellnew(parms->ndm, 1);
    Real siglev=100;//irrelevant in noise free case.
    for(int i=0; i<parms->ndm; i++){
	dmrealsq->p[i]=mapnew2(aloc->map);
	dmreal->p[i]=dnew(aloc->nloc,1);
    }
    gpu_dmreal2gpu(dmrealsq, &parms->dm[idm]);
    Real poke=1e-6;
    Real poke1=1./poke;
    curmat *phiout=curnew(pywfs->locfft->loc->nloc,1);
    const int nsa=cupowfs->saloc->nloc;
    curmat *ints=curnew(nsa,4);
    curmat *grad=curnew(nsa*2,1);
    curmat *grad0=curnew(nsa*2,1);
    dmat *gradc=dnew(nsa*2,1);

    curmat *opd0=0;
    cp2gpu(&opd0, pywfs->atm);
    if(opd0) curadd(&phiout, 1, opd0, 1, stream);
    cuzero(ints, stream);
    pywfs_ints(ints, phiout, cuwfs, siglev, stream);
    pywfs_grad(grad0, ints, cupowfs->saa, cuwfs->isum, cupowfs->pyoff, pywfs->gain, stream);
    dsp *gg=dspnew(nsa*2, aloc->nloc, nsa*2*aloc->nloc);
    int count=0;
    TIC;tic;
    for(int iloc=0; iloc<aloc->nloc; iloc++){
	dmreal->p[idm]->p[iloc]=poke;
	if(iloc>0){
	    dmreal->p[idm]->p[iloc-1]=0;
	}
	loc_embed(dmrealsq->p[idm], aloc, dmreal->p[idm]->p);
	gpu_dmreal2gpu(dmrealsq, 0);
	if(opd0){
	    curcp(&phiout, opd0, stream);
	}else{
	    cuzero(phiout, stream);
	}
	int ipowfs=parms->wfs[iwfs].powfs;
	gpu_dm2loc(phiout->p, cuwfs->loc_dm, cudata->dmreal, cudata->ndm,
		   parms->powfs[ipowfs].hs, parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, 
		   0, 0, 1, stream);
	cuzero(ints, stream);
	pywfs_ints(ints, phiout, cuwfs, siglev, stream);
	pywfs_grad(grad, ints, cupowfs->saa, cuwfs->isum, cupowfs->pyoff, pywfs->gain,stream);
	curadd(&grad, 1, grad0, -1, stream);
	dzero(gradc);
	cp2cpu(&gradc, grad, stream);
	gg->p[iloc]=count;
	const Real thres=dmaxabs(gradc)*EPS;
	for(int ig=0; ig<gradc->nx; ig++){
	    if(fabs(gradc->p[ig])>thres){
		gg->x[count]=gradc->p[ig]*poke1;
		gg->i[count]=ig;
		count++;
	    }
	}
	if(iloc%10==0){
	    Real ts=myclockd()-tk;
	    info2("%d of %ld. %.2f of %.2f seconds\n", iloc, aloc->nloc, ts, ts/(iloc+1)*aloc->nloc);
	}
    }
    gg->p[aloc->nloc]=count;
    dspsetnzmax(gg, count);
    cufree(grad0);
    cufree(opd0);
    cufree(grad);
    cufree(ints);
    dfree(gradc);
    cufree(phiout);
    cellfree(dmreal);
    cellfree(dmrealsq);
    return gg;
}
