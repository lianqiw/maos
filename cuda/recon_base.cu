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
#include "recon_base.h"
#include "curmat.h"
#include "accphi.h"
namespace cuda_recon{
curecon_geom::curecon_geom(const PARMS_T *parms, const RECON_T *recon)
    :npsr(0),ndm(0),isim(0),isimr(0),cubic_cc(0),
     amap(0),xmap(0),xcmap(0),W01(0),
     xnx(0),xny(0),anx(0),any(0),
     anloc(0),ngrad(0),
     dt(0),delay(0){
    if(!parms) return;
    ndm=parms->ndm;
    npsr=parms->sim.idealfit?parms->atm.nps:recon->npsr;
    pmap.init(recon->pmap);
    fmap.init(recon->fmap);
    /*Setup various grid*/
    amap=new cugrid_t[ndm];
    for(int idm=0; idm<ndm; idm++){
	amap[idm].init(recon->amap[idm]);
	amap[idm].cubic_cc=gpu_dmcubic_cc(parms->dm[idm].iac);
    }
    if(!parms->sim.idealfit){
	xmap=new cugrid_t[npsr];
	for(int ipsr=0; ipsr<npsr; ipsr++){
	    xmap[ipsr].init(recon->xmap[ipsr]);
	}
	xnx=recon->xnx;
	xny=recon->xny;
    }
    if(parms->fit.cachex){
	xcmap=new cugrid_t[npsr];
	for(int ipsr=0; ipsr<npsr; ipsr++){
	    xcmap[ipsr].init(recon->xcmap[ipsr]);
	}
    }
    cubic_cc=curcellnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	if(parms->dm[idm].cubic){
	    cubic_cc->p[idm]=gpu_dmcubic_cc(parms->dm[idm].iac);
	}
    }
    W01=gpu_get_W01(recon->W0, recon->W1);
    anx=recon->anx;
    any=recon->any;
    anloc=recon->anloc;
    ngrad=recon->ngrad;
    dt=parms->sim.dt;
    delay=2;//2 frame delay
}

map_l2d::map_l2d(const cugrid_t &out, int _ndir, //output.
		 const cugrid_t *in, int _nlayer,//input. layers.
		 float *thetaxv, float *thetayv, float *hsv, int *skip, float dt){//directions and star height.
    nlayer=_nlayer;
    ndir=_ndir;
    PROP_WRAP_T *hdata_cpu=new PROP_WRAP_T[nlayer*ndir];
    cudaMalloc(&hdata, sizeof(PROP_WRAP_T)*nlayer*ndir);
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	const float ht=in[ilayer].ht;
	for(int idir=0; idir<ndir; idir++){
	    if(!skip || !skip[idir]){
		const float hs=hsv?hsv[idir]:INFINITY;
		const float dispx=thetaxv[idir]*ht+in[ilayer].vx*dt;
		const float dispy=thetayv[idir]*ht+in[ilayer].vy*dt;
		const float scale=1.f-ht/hs;
		cugrid_t outscale=out*scale;
		gpu_prop_grid_prep(hdata_cpu+idir+ilayer*ndir, outscale, in[ilayer],
				   dispx, dispy, in[ilayer].cubic_cc);
	    }
	    hdata_cpu[idir+ilayer*ndir].togpu(hdata+idir+ilayer*ndir);
	}
    }
    delete [] hdata_cpu;
}

map_l2l::map_l2l(const cugrid_t *out, const cugrid_t *in, int _nlayer){//input. layers.
    nlayer=_nlayer;
    ndir=1;
    PROP_WRAP_T *hdata_cpu=new PROP_WRAP_T[nlayer];
    cudaMalloc(&hdata, sizeof(PROP_WRAP_T)*nlayer);
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	if(fabs(out[ilayer].ht-in[ilayer].ht)>EPS){
	    error("Layer height mismatch\n");
	}
	gpu_prop_grid_prep(hdata_cpu+ilayer, out[ilayer], in[ilayer],
			   0, 0, in[ilayer].cubic_cc);
	hdata_cpu[ilayer].togpu(hdata+ilayer);
    }
    delete [] hdata_cpu;
}
}//namespace
