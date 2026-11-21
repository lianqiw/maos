/*
  Copyright 2009-2026 Lianqi Wang
  
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

#include "../math/mathdef.h"
#include "rayleigh.h"
/**
   Compute rayleigh backscatter patterns for all WFS. All input are in SI unit.

   @param[in] cfg The configuration, including subaperture grid, WFS directions, LLT locations and scattering layers.
   @return Newly allocated nsa*nwfs cell of npsf*npsf backscatter PSF for each subaperture, or NULL if cfg is NULL; the caller frees it.
*/
dcell *rayleigh(rayleigh_t *cfg){
	if(!cfg) return NULL;
	const loc_t *saloc = cfg->saloc;
    const dmat *thetax = cfg->thetax;
    const dmat *thetay = cfg->thetay;
    const dmat *llt_ox = cfg->llt_ox;
    const dmat *llt_oy = cfg->llt_oy;
    const dmat *tau    = cfg->tau;
    
    const real dsa    = cfg->dsa?cfg->dsa:saloc->dx;
    const real sbeam  = cfg->sbeam;
    const real dbeam  = cfg->dbeam;
    const real hs     = cfg->hs;
    const real dtheta = cfg->dtheta;
    const int npsf    = cfg->npsf;
	dbg("dsa=%g, nsa=%ld, dx=%g sbeam=%g, dbeam=%g, hs=%g, dtheta=%g\n", dsa, saloc->nloc, saloc->dx, sbeam, dbeam, hs, dtheta);
	const int nwfs=(int)PN(thetax);
	assert(nwfs==PN(thetay));
	assert(PN(llt_ox)==1 || PN(llt_ox)==nwfs);
	assert(PN(llt_oy)==1 || PN(llt_oy)==nwfs);
	const long nsa=saloc->nloc;
	const real thres0=dbeam*dbeam*0.25;//beam radius squared
	const real gscale=-1./(2*sbeam*sbeam);
	const real dsa2=dsa*0.5;
	const real pray=0.73*(1.0418+0.986)/(4*M_PI);//Rayleigh scatter phase function at theta=pi divided by 4pi
	dcell *out=dcellnew_same(nsa, nwfs, npsf, npsf);
	int npsf2=npsf/2;
	for(int ih=0; ih<NX(tau); ih++){
		const real ht=P(tau, ih, 0);//height
		//total fraction of photons collected by a subaperture with infinity focal plane
		const real alpha=P(tau, ih, 1)*P(tau, ih, 1)*P(tau, ih, 2)*pray*abs2(dsa/ht);
		const real scale=1-ht/hs; //cone effect 
		const real thres_fov=abs2(dbeam*0.5+dtheta*npsf*ht);//threshold for sep2 to check beam is outside
		const real p2s=dtheta*ht;//convert pixel index to spatial at layer
		for(int iwfs=0; iwfs<nwfs; iwfs++){//viewing WFS
			for(int ilgs=0; ilgs<nwfs; ilgs++){//LGS beam
				real x_beam=P(thetax, ilgs)*ht+PR(llt_ox, ilgs)*scale; //beam center at this layer
				real y_beam=P(thetay, ilgs)*ht+PR(llt_oy, ilgs)*scale; //beam center at this layer
				for(long isa=0; isa<nsa; isa++){
					real x_sa=(saloc->locx[isa]+dsa2)*scale+P(thetax, iwfs)*ht; //subaperture optical axis location at this layer.
					real y_sa=(saloc->locy[isa]+dsa2)*scale+P(thetay, iwfs)*ht;
					real sep2=abs2(x_sa-x_beam)+abs2(y_sa-y_beam);
					if(sep2<thres_fov){//inside of the FoV
						for(int ipy=0; ipy<npsf; ipy++){//loop over the PSF points
							real y_pix=y_sa+(ipy-npsf2)*p2s;
							for(int ipx=0; ipx<npsf; ipx++){
								real x_pix=x_sa+(ipx-npsf2)*p2s;
								real pix_sep2=abs2(x_pix-x_beam)+abs2(y_pix-y_beam);
								if(pix_sep2<thres0){//within beam
									P(P(out, isa, iwfs), ipx, ipy)+=exp(gscale*pix_sep2)*alpha;
								}
							}
						}
					}
				}//for jwfs
			}//for isa
		}//for iwfs
	}//for ih
	return out;
}//function
/**
   Free rayleigh_t created by rayleigh_setup().

   @param[in] cfg The configuration to free. Can be NULL.
*/
void rayleigh_free(rayleigh_t *cfg){
	if(!cfg) return;
	locfree(cfg->saloc);
	dfree(cfg->thetax);
	dfree(cfg->thetay);
	dfree(cfg->llt_ox);
	dfree(cfg->llt_oy);
	dfree(cfg->tau);
	free(cfg);
}
/**
   Create a test configuration of rayleigh_t for a telescope of diameter D (default 5 m):
   annular subaperture grid, 6 LGS WFS (one on axis and five on a 35 arcsec ring) and a
   single Rayleigh scattering layer. Free with rayleigh_free().

   @param[in] D Telescope diameter (m). 0 for the default of 5 m.
   @return Newly allocated configuration; free with rayleigh_free().
*/
rayleigh_t *rayleigh_setup(real D){
	rayleigh_t *cfg=mycalloc(1, rayleigh_t);
	if(!D) D=5;
	cfg->saloc=mkannloc(D,1.8,D/2.,0.3);
	const int nwfs=6;
	cfg->llt_ox=dnew(nwfs,1);
	cfg->llt_oy=dnew(nwfs,1);
	dmat *thetax=cfg->thetax=dnew(nwfs,1);
	dmat *thetay=cfg->thetay=dnew(nwfs,1);
	
	real rast=35./206265.;
	for(int iwfs=1; iwfs<nwfs; iwfs++){
		P(thetax, iwfs)=rast*cos((real)iwfs/(nwfs-1)*TWOPI);
		P(thetay, iwfs)=rast*sin((real)iwfs/(nwfs-1)*TWOPI);
	}
	const int nh=1;
	const real dh=200;
	const real h0=30e3;
	const dmat *tau=cfg->tau=dnew(nh,3);
	P(tau, 0, 0)=h0+dh;
	P(tau, 0, 1)=1;
	P(tau, 0, 2)=1e-2;
	for(int ih=1; ih<nh; ih++){
		P(tau, ih, 0)=h0+(ih+1)*dh;
		P(tau, ih, 2)=1e-2;//fraction of scattering
		P(tau, ih, 1)=P(tau, ih-1, 1)*(1-P(tau, ih, 2));
	}

	cfg->dsa=cfg->saloc->dx;
	cfg->sbeam=0.12;
	cfg->dbeam=0.35;
	cfg->hs=180e3;
	cfg->dtheta=0.2/206265;//must be small to sampling the laser beam properly
	cfg->npsf=640;
	return cfg;
}
