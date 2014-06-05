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

#include "maos.h"
#include "setup_recon.h"
#include "recon_utils.h"

/**
   \file setup_recon_dm.c 

   Setup grid and ray tracing operators regarding DM. This is independent of
   1) WFS geometry or noise parameters
   2) Tomography
*/


/**
   Like ploc, but for DM fitting
*/
static void
setup_recon_floc(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    double dxr=parms->atmr.dx/parms->fit.pos;/*sampling of floc */
    if(parms->load.floc){
	warning("Loading floc from %s\n", parms->load.floc);
	recon->floc=locread("%s", parms->load.floc);
	if(fabs(recon->floc->dx-dxr)>dxr*0.01){
	    warning("Loaded floc has unexpected sampling of %g, should be %g\n", 
		    recon->floc->dx, dxr);
	}
    }else{
	double guard=parms->tomo.guard*dxr;
	map_t *fmap=0;
	create_metapupil(&fmap,0,0,parms,0,dxr,dxr,0,guard,0,0,0,parms->fit.square);
	info2("FLOC is %ldx%ld, with sampling of %.2fm\n",fmap->nx,fmap->ny,dxr);
	recon->floc=map2loc(fmap);/*convert map_t to loc_t */
	mapfree(fmap);
    }
    loc_create_map_npad(recon->floc, parms->fit.square?0:1);
    recon->fmap=recon->floc->map;
    /*create the weighting W for bilinear influence function. See [Ellerbroek 2002] */
    if(parms->load.W){
	if(!(zfexist("W0")&&zfexist("W1"))){
	    error("W0 or W1 not exist\n");
	}
	warning("Loading W0, W1");
	recon->W0=spread("W0");
	recon->W1=dread("W1");
    }else{
	/*
	  Compute W0,W1 weighting matrix that can be used to compute piston
	  removed wavefront variance of OPD: RMS=OPD'*(W0-W1*W1')*OPD; W0 is
	  sparse. W1 is a vector. These matrices are used for weighting in DM
	  fitting.
	*/
	double rin=0;
	if(parms->dbg.annular_W && parms->aper.din>0){
	    warning("Define the W0/W1 on annular aperture instead of circular.\n");
	    rin=parms->aper.din/2;
	}
	mkw_annular(recon->floc, 0, 0, rin, parms->aper.d/2,
		    &(recon->W0), &(recon->W1));
    }
    if(parms->save.setup){
	spwrite(recon->W0, "%s/W0",dirsetup);
	dwrite(recon->W1, "%s/W1",dirsetup);
    }
    if(parms->save.setup){
	locwrite(recon->floc, "%s/floc",dirsetup);
    }
    loc_create_stat(recon->floc);
}
/**
   Setup the deformable mirrors grid aloc. This is used for DM fitting.
*/
static void
setup_recon_aloc(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    if(ndm==0) return;
    CALL_ONCE;
    if(parms->fit.cachedm){
	recon->acmap=calloc(ndm, sizeof(map_t*));
    }
    if(parms->load.aloc){
	char *fn=parms->load.aloc;
	warning("Loading aloc from %s\n",fn);
	int ndm0;
	recon->aloc=locarrread(&ndm0,"%s",fn);
	if(ndm0!=ndm) error("Invalid saved aloc");
	for(int idm=0; idm<ndm; idm++){
	    if(fabs(parms->dm[idm].dx-recon->aloc[idm]->dx)>1e-7){
		error("DM[%d]: loaded aloc has dx=%g while dm.dx=%g\n", idm, 
		      recon->aloc[idm]->dx, parms->dm[idm].dx);
	    }
	    double max,min;
	    dmaxmin(recon->aloc[idm]->locx,recon->aloc[idm]->nloc, &max, &min);
	    if(max-min<parms->aper.d){
		warning("DM[%d]: loaded aloc is too small: diameter is %g while aper.d is %g\n", 
		      idm, max-min, parms->aper.d); 
	    }
	}
    }else{
	recon->aloc=calloc(ndm, sizeof(loc_t*));
	/*int nxmax=0, nymax=0; */
	for(int idm=0; idm<ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dx=parms->dm[idm].dx;
	    double dy=parms->dm[idm].dy;
	    double offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
	    double guard=parms->dm[idm].guard*MAX(dx,dy);
	    map_t *map;
	    if(parms->dbg.dmfullfov && !parms->fit.square){//DM covers full fov
		double D=(parms->sim.fov*ht+parms->aper.d+guard*2);
		long nx=D/dx+1;
		long ny=D/dy+1;
		map=mapnew(nx, ny, dx, dy, 0);
		map->h=ht;
		map->ox+=offset*dx;
		mapcircle_symbolic(map, D*0.5);
	    }else{
		create_metapupil(&map,0,0,parms,ht,dx,dy,offset,guard,0,0,0,parms->fit.square);
	    }
	    info2("DM %d: grid is %ld x %ld\n", idm, map->nx, map->ny);
	    recon->aloc[idm]=map2loc(map);
	    mapfree(map);
	}
    }
    recon->amap=calloc(parms->ndm, sizeof(map_t*));
    for(int idm=0; idm<parms->ndm; idm++){
	double ht=parms->dm[idm].ht;
	double offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
	double dx=parms->dm[idm].dx;
	loc_create_map_npad(recon->aloc[idm], parms->fit.square?0:1);
	recon->amap[idm]=recon->aloc[idm]->map;
	recon->amap[idm]->h=ht;
	recon->amap[idm]->cubic=parms->dm[idm].cubic;
	recon->amap[idm]->iac=parms->dm[idm].iac;
	if(parms->fit.cachedm){
	    const double dx2=parms->atmr.dx/parms->fit.pos;
	    const double dy2=dx2;
	    create_metapupil(&recon->acmap[idm],0,0,
			     parms,ht,dx2,dy2,offset*dx/dx2,dx2,0,0,0,parms->fit.square);
	    info("amap origin is %g, %g. acmap is %g, %g\n", 
		 recon->aloc[idm]->map->ox, recon->aloc[idm]->map->oy,
		 recon->acmap[idm]->ox, recon->acmap[idm]->oy);
	}
    }

    if(parms->save.setup){
	locarrwrite(recon->aloc,parms->ndm,"%s/aloc",dirsetup);
	maparrwrite(recon->amap, parms->ndm,"%s/amap", dirsetup);
    }
    recon->aimcc=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc[idm], NULL);
	dinvspd_inplace(recon->aimcc->p[idm]);
    }
    recon->anx=calloc(ndm, sizeof(long));
    recon->any=calloc(ndm, sizeof(long));
    recon->anloc=calloc(ndm, sizeof(long));
    for(int idm=0; idm<ndm; idm++){
	recon->anx[idm]=recon->amap[idm]->nx;
	recon->any[idm]=recon->amap[idm]->ny;
	recon->anloc[idm]=recon->aloc[idm]->nloc;
    }
    /*Dealing with stuck/floating actuators. */
    int anyfloat=0, anystuck=0;
    for(int idm=0; idm<ndm; idm++){
	if(parms->dm[idm].actstuck) anystuck=1;
	if(parms->dm[idm].actfloat) anyfloat=1;
    }
    if(anystuck){
	recon->actstuck=icellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actstuck) continue;
	    recon->actstuck->p[idm]=act_coord2ind(recon->aloc[idm], parms->dm[idm].actstuck);
	}
    }
    if(anyfloat){
	recon->actfloat=icellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actfloat) continue;
	    recon->actfloat->p[idm]=act_coord2ind(recon->aloc[idm], parms->dm[idm].actfloat);
	}
    }
}

/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static void
setup_recon_HA(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    if(parms->load.HA && zfexist(parms->load.HA)){
	warning("Loading saved HA\n");
	recon->HA=spcellread("%s",parms->load.HA);
    }else{
	const int nfit=parms->fit.nfit;
	const int ndm=parms->ndm;
	recon->HA=spcellnew(nfit, ndm);
	PDSPCELL(recon->HA,HA);
	info2("Generating HA ");TIC;tic;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.hs[ifit];
	    for(int idm=0; idm<ndm; idm++){
		const double ht=parms->dm[idm].ht;
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		loc_t *loc=recon->floc;
		if(parms->recon.misreg_dm2sci && parms->recon.misreg_dm2sci[ifit+idm*nfit]){
		    loc=loctransform(loc, parms->recon.misreg_dm2sci[ifit+idm*nfit]);
		}
		HA[idm][ifit]=mkh(recon->aloc[idm], loc, NULL,
				  displace[0], displace[1], 
				  scale,parms->dm[idm].cubic,parms->dm[idm].iac);
		if(loc!=recon->floc){
		    locfree(loc);
		}
	    }
	}
	toc2(" ");
    }
    if(parms->save.setup){
	spcellwrite(recon->HA,"%s/HA",dirsetup);
    }
    recon->actcpl=genactcpl(recon->HA, recon->W1);
    if(recon->actstuck){
	warning2("Apply stuck actuators to HA\n");
	act_stuck(recon->aloc, recon->HA, NULL, recon->actstuck);
    	if(parms->save.setup){
	    spcellwrite(recon->HA,"%s/HA_stuck",dirsetup);
	}
    }
    if(recon->actfloat){
	warning2("Apply float actuators to HA\n");
	act_float(recon->aloc, &recon->HA, NULL, recon->actfloat);
	recon->actinterp=act_float_interp(recon->aloc, recon->actfloat);
	if(parms->save.setup){
	    spcellwrite(recon->HA,"%s/HA_float",dirsetup);
	}
    }
    if(parms->fit.actinterp){
	warning2("Apply slaving actuator operation\n");
	spcell *interp2=act_extrap(recon->aloc, recon->actcpl, 0.5);
	spcelladd(&recon->actinterp, interp2);
	spcellfree(interp2);
    }
    if(parms->save.setup && recon->actinterp){
	spcellwrite(recon->actinterp, "%s/actinterp", dirsetup);
    }
}

/**
   Setup fitting low rank terms that are in the NULL space of DM fitting
   operator. typically include piston on each DM and tip/tilt on certain
   DMs. Becareful with tip/tilt contraint when using CBS.  */
static void 
fit_prep_lrt(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    const int ndm=parms->ndm;
    if(ndm>=3) warning("Low rank terms for 3 or more dms are not tested\n");
    recon->fitNW=dcellnew(ndm,1);
    double scl=recon->fitscl=1./recon->floc->nloc;
    if(fabs(scl)<1.e-15){
	error("recon->fitscl is too small\n");
    }
    /*computed outside. */
    int lrt_tt=parms->fit.lrt_tt;
    int nnw=0;
    if(parms->fit.lrt_piston){
	nnw+=ndm;
    }
    if(lrt_tt){
	nnw+=2*(ndm-1);
    }
    if(nnw==0) return;
    for(int idm=0; idm<ndm; idm++){
	int nloc=recon->aloc[idm]->nloc;
	recon->fitNW->p[idm]=dnew(nloc, nnw);
    }
    int inw=0;/*current column */
    if(parms->fit.lrt_piston){
	info2("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc[idm]->nloc;
	    double *p=recon->fitNW->p[idm]->p+(inw+idm)*nloc;
	    const double *cpl=recon->actcpl->p[idm]->p;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(cpl[iloc]>0.5){
		    p[iloc]=scl;
		}
	    }
	}
	inw+=ndm;
    }
    if(lrt_tt){
	double factor=0;
	if(lrt_tt==1){
	    info2("Adding TT cr on upper DMs to fit matrix.\n");
	    factor=scl*2./parms->aper.d;
	    for(int idm=1; idm<ndm; idm++){
		int nloc=recon->aloc[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+(inw+(idm-1)*2)*nloc;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=recon->actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.5){
			p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }
	}else if(lrt_tt==2){/*Canceling TT. only valid for 2 DMs */
	    warning("Adding Canceling TT cr to fit matrix. Deprecated\n");
	    if(ndm!=2){
		error("Only ndm=2 case is implemented\n");
	    }
	    for(int idm=0; idm<ndm; idm++){
		int nloc=recon->aloc[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+inw*nloc;
		if(idm==0) factor=scl*2/parms->aper.d;
		else if(idm==1) factor=-scl*2./parms->aper.d;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=recon->actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.5){
			p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }

	}
	inw+=2*(ndm-1);
    }
    if(parms->fit.actslave){
	/*
	  2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found
	  that slaving is causing mis-measurement of a few edge
	  actuators. First try to remove W1. Or lower the weight. Revert
	  back.
	  1./floc->nloc is on the same order of norm of Ha'*W*Ha. 
	*/
	recon->actslave=slaving(recon->aloc, recon->actcpl,
				recon->fitNW, recon->actstuck,
				recon->actfloat, .5, 1./recon->floc->nloc);
	if(parms->save.setup){
	    dcellwrite(recon->actcpl, "%s/actcpl", dirsetup);
	    spcellwrite(recon->actslave,"%s/actslave",dirsetup);
	}
    }
    if(parms->save.setup){
	dcellwrite(recon->fitNW,"%s/fitNW",dirsetup);
    }
}

void setup_recon_dm(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    /*setup DM actuator grid */
    setup_recon_aloc(recon,parms);
    /*Grid for DM fitting*/
    setup_recon_floc(recon,parms);
    if(parms->recon.alg==0 || parms->sim.ncpa_calib){
	setup_recon_HA(recon,parms);
	fit_prep_lrt(recon, parms);
    }
}
