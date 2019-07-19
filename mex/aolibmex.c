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
/**
  Wrap of the maos functions for easy calling from mex.
*/
#include "../lib/aos.h"
ccell *genotfmex(loc_t *loc, const dmat *amp, const dmat *opdbias, const dmat *area, double thres, double wvl, double dtheta, const dmat *cov, double r0, double l0, long ncompx, long ncompy, long nsa, long pttr){
    ccell *out=ccellnew(nsa, 1);
    genotf(out->p, loc, amp, opdbias, area, thres, wvl, dtheta, cov, r0, l0, ncompx, ncompy, nsa, pttr);
    return out;
}
dmat *m3projmex(dmat *mapin_0, char *header, loc_t *locout, double thetax, double thetay, double hs){
    free(mapin_0->header);
    mapin_0->header=header;
    rmap_t *mapin=d2rmap(mapin_0);
    dmat *opd=dnew(locout->nloc, 1);
    m3proj(mapin, opd, locout, thetax, thetay, hs);
    cellfree(mapin);
    cellfree(mapin_0);
    return opd;
}
dmat *mkcirmap(long nx, long ny, double cx, double cy, double r){
    dmat *map=dnew(nx, ny);
    dcircle(map, cx, cy, 1,1, r, 1);
    return map;
}
cn2est_t *cn2estmex(const dmat *wfspair, dmat *wfstheta, const loc_t *saloc,
		  const dmat *saa, const double saat, 
		  const dmat* hs, const dmat *htrecon, int keepht, double l0, dcell *grad){
    if(maxabs(wfstheta->p, wfstheta->nx*wfstheta->ny)>1){
	dmat *tmp=wfstheta;
	wfstheta=ddup(tmp);
	dfree(tmp);
	//Don't scale the matlab one.
	dscale(wfstheta, 1./206265);
    }
    if(grad->nx==1 && grad->ny>1){
	grad->nx=grad->ny;
	grad->ny=1;
    }
    struct cn2est_t *cn2est=cn2est_new(wfspair, wfstheta, saloc, saa, saat, hs, htrecon, keepht, l0);
    cn2est_push(cn2est, grad);
    cn2est_est(cn2est, 1, 0);
    return cn2est;
}
dmat *mtch2(dmat **nea, const dmat *i0, const dmat *gx, const dmat *gy, int cr){
    return mtch(nea, i0, gx, gy, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, cr);
}
dmat *sdepsd(const dmat *ff, const dmat *coeff){
    dmat *psd=0;
    sde_psd(&psd, ff, coeff->p, coeff->nx, coeff->ny);
    return psd;
}
ccell *mkdtfmex(dspcell **si, const dmat *wvls, double dxsa, double embfac, long ncompx, long ncompy, long pixpsax, long pixpsay, double pixthetax, double pixthetay, const dmat* pixoffx, const dmat* pixoffy, double pixblur){
    DTF_T *dtf=mkdtf(wvls, dxsa, embfac, ncompx, ncompy, pixpsax, pixpsay, pixthetax, pixthetay, pixoffx, pixoffy, pixblur, 0, 0, 0);
    int nwvl=wvls->nx;
    ccell *nominal=ccellnew(nwvl, 1);
    *si=dspcellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	nominal->p[iwvl]=cref(dtf[iwvl].nominal->p[0]);
	(*si)->p[iwvl]=dspref(dtf[iwvl].si->p[0]);
	
	cfree(dtf[iwvl].Ux);
	cfree(dtf[iwvl].Uy);
    }
    free(dtf);
    return nominal;
}
///Embed arr according to loc into a 2d array
cell *loc_embed2(loc_t *loc, dmat *arr){
    if(!loc->map){
	loc_create_map(loc);
    }
    if(arr->nx==1 && arr->ny!=1){
	arr->nx=arr->ny;
	arr->ny=1;
    }
    int nx=arr->nx/loc->nloc;
    int ny=arr->ny;
    if(nx*loc->nloc!=arr->nx){
	error("arr has wrong dimension: %ldx%ld. loc length is %ld \n", arr->nx, arr->ny, loc->nloc);
    }
    dcell *dest=dcellnew(nx, ny);
    for(int ix=0; ix<nx*ny; ix++){
	P(dest, ix)=dnew(loc->map->nx, loc->map->ny);
	loc_embed((map_t*)P(dest, ix), loc, arr->p+ix*loc->nloc);
    }
    if(nx==1 && ny==1){
	dmat *dest0=dref(P(dest,0));
	cellfree(dest);
	return (cell*)dest0;
    }else{
	return (cell*)dest;
    }
}
