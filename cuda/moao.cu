/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "moao.h"
#include "utils.h"
#include "recon.h"
#include "accphi.h"
#include "pcg.h"
#include "cudata.h"
#include "perf.h"

namespace cuda_recon{
cumoao_t::cumoao_t(const PARMS_T *parms, MOAO_T *moao, dir_t *dir, int _ndir, curecon_geom *_grid)
    :cusolve_cg(parms?parms->fit.maxit:0, parms?parms->recon.warm_restart:0),grid(_grid),ndir(_ndir){
    amap=cugridcell(1,1);
    amap[0]=(moao->amap->p[0]);
    if(moao->NW){
	cp2gpu(NW, moao->NW->p[0]);
	dotNW=curmat(NW.Ny(), 1);
    }
    if(moao->actslave){
	actslave=cusp(moao->actslave->p[0], 1);
    }

    dir_t dir0;
    ha.init_l2d(grid->fmap, &dir0, 1, amap);
    opdfit=curcell(1,1,grid->fmap.nx,grid->fmap.ny);
    opdfit2=curcell(1,1,grid->fmap.nx,grid->fmap.ny);

    hxp=Array<map2map>(ndir, 1);
    hap=Array<map2map>(ndir, 1);
    for(int idir=0; idir<ndir; idir++){
	hxp[idir].init_l2d(grid->fmap, dir+idir, 1, grid->xmap);
	hap[idir].init_l2d(grid->fmap, dir+idir, 1, grid->amap);
    }
    rhs=curcell(1,1,amap[0].nx,amap[0].ny);
}
Real cumoao_t::moao_solve(curccell &xout, const curcell &xin, const curcell &ain, stream_t &stream){
    for(int idir=0; idir<ndir; idir++){
	cuzero(opdfit.M(), stream);
	//cuwrite(xin, "xin_%d", idir);
	//dbg("hxp[%d]\n", idir);
	hxp[idir].forward(opdfit.pm, xin.pm, 1.f, NULL, stream);//tomography
	//cuwrite(opdfit, "opdfit0_%d", idir);
	//cuwrite(ain, "ain_%d", idir);
	//dbg("hap[%d]\n", idir);
	hap[idir].forward(opdfit.pm, ain.pm, -1.f, NULL, stream);//minus common DM.
	//cuwrite(opdfit, "opdfit1_%d", idir);
	grid->W01.apply(opdfit2.M()(), opdfit.M()(), opdfit.Nx(), stream);
	//cuwrite(opdfit2, "opdfit2_%d", idir);
	cuzero(rhs.M(), stream);
	ha.backward(opdfit2.pm, rhs.pm, 1, NULL, stream);
	//cuwrite(rhs, "rhs_%d", idir);
	solve(xout[idir], rhs, stream);
	//cuwrite(xout[idir], "xout_%d", idir);
	/*{
	    static int ic=-1; ic++;
	    cuwrite(xout[idir], "xout_%d", ic);
	    cuwrite(rhs, "rhs_%d", ic);
	    cuwrite(xin, "xin_%d", ic);
	    cuwrite(ain, "ain_%d", ic);
	    cuwrite(opdfit, "opd_%d", ic);
	    }*/
    }
    return 0;
}
void cumoao_t::L(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream){
    if(!xout){
	xout=curcell(1, 1, amap[0].nx, amap[0].ny);
    }else{
	curscale(xout.M(), beta, stream);
    }
    cuzero(opdfit.M(), stream);
    ha.forward(opdfit.pm, xin.pm, 1, NULL, stream);
    grid->W01.apply(opdfit2.M()(), opdfit.M()(), opdfit.Nx(), stream);
    ha.backward(opdfit2.pm, xout.pm, alpha, NULL, stream);
    if(NW){
	curmv(dotNW(), 0, NW, xin.M()(), 't', 1, stream);
	curmv(xout.M()(), 1, NW, dotNW(), 'n', alpha, stream);
    }
    if(actslave){
	cuspmul(xout.M()(), actslave, xin.M()(), 1,'n', alpha, stream);
    }
}
}//namespace
/*
  Embed and copy DM commands to GPU.
*/
static void gpu_dm2gpu_embed(curmat &dmgpu, dmat *dmcpu, loc_t *loc, int nx, int ny){
    assert(dmcpu->ny==1);
    Real *pout=(Real*)calloc(nx*ny, sizeof(Real));
    map_t *map=loc->map;
    real *pin=dmcpu->p-1;
    for(long i=0; i<map->nx*map->ny; i++){
	long iphi=(long)map->p[i];
	if(iphi){
	    pout[i]=pin[iphi];
	}
    }
    DO(cudaMemcpy(dmgpu(), pout, nx*ny*sizeof(Real), cudaMemcpyHostToDevice));
    free(pout);
}

/**
   Copy MOAO DM commands from CPU to GPU.*/
void gpu_moao_2gpu(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->gpu.moao){
	error("Invalid use\n");
    }
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    if(parms->gpu.wfs && simu->dm_wfs){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0) continue;
	    MOAO_T *moao=recon->moao+imoao;
	    gpu_set(cuglobal->wfsgpu[iwfs]);
	    if(!cudata->dm_wfs){
		cudata->dm_wfs=Array<cumapcell>(nwfs, 1);
	    }
	    if(!cudata->dm_wfs[iwfs]){
		cudata->dm_wfs[iwfs]=cumapcell(1,1);
		cudata->dm_wfs[iwfs][0]=(recon->moao[imoao].amap->p[0]); 
	    }
	    if(parms->fit.square){
		cp2gpu(cudata->dm_wfs[iwfs][0], simu->dm_wfs->p[iwfs]);
	    }else{
		gpu_dm2gpu_embed(cudata->dm_wfs[iwfs][0], simu->dm_wfs->p[iwfs],
				 moao->aloc->p[0], moao->amap->p[0]->nx, moao->amap->p[0]->ny);
	    }
	}
    }
    if(parms->gpu.evl && simu->dm_evl){
	int imoao=parms->evl.moao;
	MOAO_T *moao=recon->moao+imoao;
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(cuglobal->evlgpu[ievl]);
	    if(!cudata->dm_evl){
		cudata->dm_evl=Array<cumapcell>(nevl, 1);
	    }
	    if(!cudata->dm_evl[ievl]){
		cudata->dm_evl[ievl]=cumapcell(1,1);
		cudata->dm_evl[ievl][0]=(moao->amap->p[0]); 
	    }
	    if(parms->fit.square){
		cp2gpu(cudata->dm_evl[ievl][0], simu->dm_evl->p[ievl]);
	    }else{
		gpu_dm2gpu_embed(cudata->dm_evl[ievl][0], simu->dm_evl->p[ievl],
				 moao->aloc->p[0], moao->amap->p[0]->nx, moao->amap->p[0]->ny);
	    }
	}
    }
}

