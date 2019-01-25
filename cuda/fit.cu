/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "fit.h"
#define TIMING 0
#if TIMING <1
#undef EVENT_INIT
#undef EVENT_TIC
#undef EVENT_TOC
#undef EVENT_PRINT
#define EVENT_INIT(A)
#define EVENT_TIC(A)
#define EVENT_TOC
#define EVENT_PRINT(A...)
#else
#define EVENT_PRINT(A...) info(A);EVENT_DEINIT
#endif

#if TIMING <2
#define EVENT2_INIT(A)
#define EVENT2_TIC(A)
#define EVENT2_TOC
#define EVENT2_PRINT(A...)
#else
#define EVENT2_INIT EVENT_INIT
#define EVENT2_TIC EVENT_TIC
#define EVENT2_TOC EVENT_TOC
#define EVENT2_PRINT(A...) info(A);EVENT_DEINIT
#endif
namespace cuda_recon{
cufit_grid::cufit_grid(const PARMS_T *parms, const RECON_T *recon, const curecon_geom *_grid)
    :cucg_t(parms?parms->fit.maxit:0, parms?parms->recon.warm_restart:0),grid(_grid),
     nfit(0),dir(0){
    if(!parms || !recon) return;
    /*Initialize*/
    const int ndm=parms->ndm;
    const int npsr=recon->npsr;
    nfit=parms->fit.nfit;
  
    if(parms->fit.cachedm){
	acmap=cugridcell(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    acmap[idm]=(recon->acmap->p[idm]);
	}
    }
    if(parms->sim.idealfit){
	idealfit=1;
	floc=culoc_t(recon->floc);
    }
    dir=new dir_t[nfit];
    for(int ifit=0; ifit<nfit; ifit++){
	dir[ifit].thetax=parms->fit.thetax->p[ifit];
	dir[ifit].thetay=parms->fit.thetay->p[ifit];
	dir[ifit].hs=parms->fit.hs->p[ifit];
	dir[ifit].skip=0;
    }
    /*Various*/
    if(recon->fit->NW){
	dmat *_fitNW=dcell2m(recon->fit->NW);
	cp2gpu(fitNW, _fitNW);
	dfree(_fitNW);
	dotNW=curmat(fitNW.Ny(), 1);
    }
    if(recon->fit->actslave){
	cp2gpu(actslave, recon->fit->actslave, 1);
    }
    if(parms->fit.cachedm){
	long acnx[ndm], acny[ndm];
	for(int idm=0; idm<ndm; idm++){
	    acnx[idm]=acmap[idm].nx;
	    acny[idm]=acmap[idm].ny;
	}
	dmcache=curcell(ndm, 1, acnx, acny);
    }
    if(parms->fit.cachex){
	long xcnx[npsr], xcny[npsr];
	for(int ips=0; ips<npsr; ips++){
	    xcnx[ips]=grid->xcmap[ips].nx;
	    xcny[ips]=grid->xcmap[ips].ny;
	}
	xcache=curcell(npsr, 1, xcnx, xcny);
    } 
    cp2gpu(fitwt, recon->fit->wt);
 
    opdfit=curcell(nfit, 1, grid->fmap.nx, grid->fmap.ny);
    opdfit2=curcell(nfit, 1, grid->fmap.nx, grid->fmap.ny);
    /*Data for ray tracing*/
    //dm -> floc
    if(parms->fit.cachex){
	hxp0.Init_l2l(grid->xcmap, grid->xmap);
	hxp1.Init_l2d(grid->fmap, dir, nfit, grid->xcmap);
    }else{
	hxp.Init_l2d(grid->fmap, dir, nfit,grid->xmap);
    }
    if(parms->fit.cachedm){
	ha0.Init_l2l(acmap, grid->amap);
	ha1.Init_l2d(grid->fmap, dir, nfit, acmap);
    }else{
	ha.Init_l2d(grid->fmap, dir, nfit, grid->amap);
    }
}

/*
  Todo: share the ground layer which is both matched and same.
*/

/**
   do HXp operation, opdfit+=Hxp*xin*/
void cufit_grid::do_hxp(const curcell &xin, stream_t &stream){
    cuzero(opdfit.M(), stream);
    if(idealfit){//ideal fiting.
	for(int ifit=0; ifit<nfit; ifit++){
	    gpu_atm2loc(opdfit[ifit](), floc, INFINITY, 0,
			dir[ifit].thetax, dir[ifit].thetay,
			0, 0, grid->dt, grid->reconisim, 1, stream);
	}
    }else{
	if(xcache){//caching
	    cuzero(xcache.M(), stream);
	    hxp0.forward(xcache.pm, xin.pm, 1, NULL, stream);
	    hxp1.forward(opdfit.pm, xcache.pm, 1, NULL, stream);
	}else{
	    hxp.forward(opdfit.pm, xin.pm, 1, NULL, stream);
	}
    }
}
/**
   do HXp' operation, xout+=alpha*Hxp'*opdfit2*/
void cufit_grid::do_hxpt(const curcell &xout, Real alpha, stream_t &stream){
    if(xcache){
	cuzero(xcache.M(), stream);
	hxp1.backward(opdfit2.pm, xcache.pm, alpha, fitwt(), stream);
	hxp0.backward(xcache.pm, xout.pm, alpha, NULL, stream);
    }else{
	hxp.backward(opdfit2.pm, xout.pm, alpha, fitwt(), stream);
    }
}

/**
   opdfit+=Ha*xin;
*/
void cufit_grid::do_ha(const curcell &xin, stream_t &stream){
    cuzero(opdfit.M(), stream); 
    if(dmcache){
	/*xout->dmcache*/ 
	cuzero(dmcache.M(), stream); 
	ha0.forward(dmcache.pm, xin.pm, 1.f, NULL, stream);
	/*dmcache->opdfit*/ 
	ha1.forward(opdfit.pm, dmcache.pm, 1.f, NULL, stream);
    }else{ 
	/*xout->opfit*/ 
	ha.forward(opdfit.pm, xin.pm, 1.f, NULL, stream);
    }
}

/**
   xout+=alpha*HA'*opdfit2*/
void cufit_grid::do_hat(curcell &xout,  Real alpha, stream_t &stream){
    if(dmcache){ 
	/*opdfit2->dmcache*/ 
	cuzero(dmcache.M(), stream); 
	ha1.backward(opdfit2.pm, dmcache.pm, alpha, fitwt(), stream);
	/*dmcache->xout*/ 
	ha0.backward(dmcache.pm, xout.pm, 1, NULL, stream);
    }else{ 
	/*opfit2->xout	*/ 
	ha.backward(opdfit2.pm, xout.pm, alpha, fitwt(), stream);
    } 
}

/*
  Right hand side operator. 
*/
void cufit_grid::R(curcell &xout, Real beta,  curcell &xin, Real alpha, stream_t &stream){
    if(!xout){
	xout=curcell(grid->ndm, 1, grid->anx, grid->any);
    }else{
	curscale(xout.M(), beta, stream);
    }
    do_hxp(xin, stream);//xin->opdfit. 153 us
    //cuwrite(opdfit.M(), "GPU_FitR_x1");
    grid->W01.apply(opdfit2.M()(), opdfit.M()(), opdfit.Nx(), stream);//opdfit->opdfit2. 123 us
    //cuwrite(opdfit2.M(), "GPU_FitR_x2");
    do_hat(xout, alpha, stream);//opdfit2->xout. 390 us
    //cuwrite(xout, "GPU_FitR_x3");
}
void cufit_grid::Rt(curcell &xout, Real beta,  curcell &xin, Real alpha, stream_t &stream){
    if(!xout){
	xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curscale(xout.M(), beta, stream);
    }
    do_ha(xin, stream);
    grid->W01.apply(opdfit2.M()(), opdfit.M()(), opdfit.Nx(), stream);
    do_hxpt(xout, alpha, stream);
}
void cufit_grid::L(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream){
    const int ndm=grid->ndm;
    EVENT_INIT(6);
    EVENT_TIC(0);
    if(!xout){
	xout=curcell(ndm, 1, grid->anx, grid->any);
    }else{
	curscale(xout.M(), beta, stream);
    }   
    do_ha(xin, stream);//112 us
    EVENT_TIC(1);
    grid->W01.apply(opdfit2.M()(), opdfit.M()(), opdfit.Nx(), stream);
    EVENT_TIC(2);
    do_hat(xout, alpha, stream);//390 us
    EVENT_TIC(3);
    if(fitNW){
	curmv(dotNW(), 0, fitNW, xin.M()(), 't', 1, stream);
	curmv(xout.M()(), 1, fitNW, dotNW(), 'n', alpha, stream);
    }
    EVENT_TIC(4);
    if(actslave){
	cuspmul(xout.M()(), actslave, xin.M()(), 1,'n', alpha, stream);
    }
    EVENT_TIC(5);
    EVENT_TOC;
    EVENT_PRINT("FitL HA: %.3f W: %.3f HAT: %.3f NW: %.3f SL: %.3f tot: %.3f\n",
		times[1],times[2],times[3],times[4],times[5],times[0]);
}
}//namespace
