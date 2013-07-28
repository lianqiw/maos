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
#include "recon_geom.h"
#include "curmat.h"
#include "accphi.h"
void curecon_geom::init(const PARMS_T *parms, const RECON_T *recon){
    ndm=parms->ndm;
    npsr=parms->sim.idealfit?parms->atm.nps:recon->npsr;
    pmap.init(recon->pmap);
    fmap.init(recon->fmap);
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
    amap=new cugrid_t[ndm];
    for(int idm=0; idm<ndm; idm++){
	amap[idm].init(recon->amap[idm]);
    }
    if(parms->fit.cachedm){
	acmap=new cumap_t[ndm];
	for(int idm=0; idm<ndm; idm++){
	    acmap[idm].init(recon->acmap[idm]);
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
    vx=vy=NULL;//no information yet.
}
