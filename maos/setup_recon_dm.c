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

#include "common.h"
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
setup_recon_floc(RECON_T *recon, const PARMS_T *parms, const APER_T *aper){
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
	create_metapupil(&fmap,0,0,parms->dirs,parms->aper.d,0,dxr,dxr,0,guard,0,0,0,parms->fit.square);
	info2("FLOC is %ldx%ld, with sampling of %.2fm\n",fmap->nx,fmap->ny,dxr);
	recon->floc=map2loc(fmap);/*convert map_t to loc_t */
	mapfree(fmap);
	/*Do not restrict fmap to within active pupil. */
    }
    loc_create_map_npad(recon->floc, parms->fit.square?0:1, 0, 0);
    recon->fmap=recon->floc->map;
    /*create the weighting W for bilinear influence function. See [Ellerbroek 2002] */
    if(parms->load.W){
	if(!(zfexist("W0")&&zfexist("W1"))){
	    error("W0 or W1 not exist\n");
	}
	warning("Loading W0, W1");
	recon->W0=dspread("W0");
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
	dspwrite(recon->W0, "%s/W0",dirsetup);
	writebin(recon->W1, "%s/W1",dirsetup);
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
setup_recon_aloc(RECON_T *recon, const PARMS_T *parms, const APER_T *aper){
    const int ndm=parms->ndm;
    if(ndm==0) return;
    if(parms->fit.cachedm){
	recon->acmap=cellnew(ndm, 1);
    }
    if(parms->load.aloc){
	char *fn=parms->load.aloc;
	warning("Loading aloc from %s\n",fn);
	recon->aloc=loccellread("%s",fn);
	if(recon->aloc->nx!=ndm) error("Invalid saved aloc");
	for(int idm=0; idm<ndm; idm++){
	    if(fabs(parms->dm[idm].dx-recon->aloc->p[idm]->dx)>1e-7){
		error("DM[%d]: loaded aloc has dx=%g while dm.dx=%g\n", idm, 
		      recon->aloc->p[idm]->dx, parms->dm[idm].dx);
	    }
	    double max,min;
	    dmaxmin(recon->aloc->p[idm]->locx,recon->aloc->p[idm]->nloc, &max, &min);
	    if(max-min<parms->aper.d){
		warning("DM[%d]: loaded aloc is too small: diameter is %g while aper.d is %g\n", 
		      idm, max-min, parms->aper.d); 
	    }
	}
    }else{
	recon->aloc=cellnew(ndm, 1);
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
		create_metapupil(&map,0,0,parms->dirs, parms->aper.d,ht,dx,dy,offset,guard,0,0,0,parms->fit.square);
	    }
	    info2("DM %d: grid is %ld x %ld\n", idm, map->nx, map->ny);
	    recon->aloc->p[idm]=map2loc(map);
	    mapfree(map);
	}
    }
    recon->amap=cellnew(parms->ndm, 1);
    for(int idm=0; idm<parms->ndm; idm++){
	double ht=parms->dm[idm].ht;
	double offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
	double dx=parms->dm[idm].dx;
	loc_create_map_npad(recon->aloc->p[idm], parms->fit.square?0:1, 0, 0);
	recon->amap->p[idm]=recon->aloc->p[idm]->map;
	recon->amap->p[idm]->h=ht;
	recon->amap->p[idm]->cubic=parms->dm[idm].cubic;
	recon->amap->p[idm]->iac=parms->dm[idm].iac;
	if(parms->fit.cachedm){
	    const double dx2=parms->atmr.dx/parms->fit.pos;
	    const double dy2=dx2;
	    create_metapupil(&recon->acmap->p[idm],0,0, parms->dirs, parms->aper.d,
			     ht,dx2,dy2,offset*dx/dx2,dx2,0,0,0,parms->fit.square);
	    info2("amap origin is %g, %g. acmap is %g, %g\n", 
		 recon->aloc->p[idm]->map->ox, recon->aloc->p[idm]->map->oy,
		 recon->acmap->p[idm]->ox, recon->acmap->p[idm]->oy);
	}
    }

    recon->aimcc=cellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc->p[idm], NULL);
	dinvspd_inplace(recon->aimcc->p[idm]);
    }
    recon->anx=lnew(ndm, 1);
    recon->any=lnew(ndm, 1);
    recon->anloc=lnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	recon->anx->p[idm]=recon->amap->p[idm]->nx;
	recon->any->p[idm]=recon->amap->p[idm]->ny;
	recon->anloc->p[idm]=recon->aloc->p[idm]->nloc;
    }
    /*Dealing with stuck/floating actuators. */
    int anyfloat=0, anystuck=0;
    for(int idm=0; idm<ndm; idm++){
	if(parms->dm[idm].actstuck) anystuck=1;
	if(parms->dm[idm].actfloat) anyfloat=1;
    }
    if(anystuck){
	recon->actstuck=cellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actstuck) continue;
	    recon->actstuck->p[idm]=act_coord2ind(recon->aloc->p[idm], parms->dm[idm].actstuck);
	    double stroke=INFINITY;
	    const int nact=recon->aloc->p[idm]->nloc;
	    if(parms->dm[idm].stroke){
		if(parms->dm[idm].stroke->nx==1){
		    stroke=parms->dm[idm].stroke->p[0];
		    dfree(parms->dm[idm].stroke);
		}else if(parms->dm[idm].stroke->nx!=nact){
		    error("dm.stroke is in wrong format\n");
		}
	    }
	    if(!parms->dm[idm].stroke){
		parms->dm[idm].stroke=dnew(nact, 2);
		for(int iact=0; iact<nact; iact++){
		    parms->dm[idm].stroke->p[iact]=-stroke;
		    parms->dm[idm].stroke->p[iact+nact]=stroke;
		}
		((PARMS_T*)parms)->sim.dmclip=1;
	    }
	    for(int iact=0; iact<nact; iact++){
		double val=recon->actstuck->p[idm]->p[iact];
		if(val){
		    parms->dm[idm].stroke->p[iact]=val*1e-9;
		    parms->dm[idm].stroke->p[iact+nact]=val*1e-9;
		}
	    }
	}
    }
    if(anyfloat){
	recon->actfloat=cellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actfloat) continue;
	    recon->actfloat->p[idm]=act_coord2ind(recon->aloc->p[idm], parms->dm[idm].actfloat);
	}
    }
    if(parms->recon.modal){
	recon->amod=cellnew(ndm, 1);
	recon->anmod=lnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    int nmod=parms->recon.nmod;
	    const long nloc=recon->aloc->p[idm]->nloc;
	    switch(parms->recon.modal){
	    case -2: {//dummy modal control, emulating zonal mode with identity modal matrix
		if(nmod && nmod!=nloc){
		    warning("recon.mod should be 0 or %ld when recon.modal=2 \n",nloc);
		    recon->amod->p[idm]=dnew(nloc, nloc);
		    double val=sqrt(nloc);
		    daddI(recon->amod->p[idm], val);
		    dadds(recon->amod->p[idm], -val/nloc);
		}
	    }
		break;
	    case -1://zernike
		{
		    if(!nmod) nmod=nloc;
		    int rmax=floor((sqrt(1+8*nmod)-3)*0.5);
		    recon->amod->p[idm]=zernike(recon->aloc->p[idm], 0, 0, rmax, 0);
		}
		break;
	    case 1://Karhunen loeve
		recon->amod->p[idm]=KL_vonkarman(recon->aloc->p[idm], nmod, parms->atmr.L0);
		break;
	    default:
		error("Invalid recon.modal");
	    }	    
	    recon->anmod->p[idm]=recon->amod->p[idm]->ny;
	}
	
    }
    if(parms->save.setup){
	writebin(recon->aloc,"%s/aloc",dirsetup);
	writebin(recon->amap, "%s/amap", dirsetup);
	if(parms->recon.modal){
	    writebin(recon->amod, "%s/amod", dirsetup);
	}
    }
}

/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static void
setup_recon_HA(RECON_T *recon, const PARMS_T *parms){
    if(parms->load.HA && zfexist(parms->load.HA)){
	warning("Loading saved HA\n");
	recon->HA=dspcellread("%s",parms->load.HA);
    }else{
	const int nfit=parms->fit.nfit;
	const int ndm=parms->ndm;
	recon->HA=cellnew(nfit, ndm);
	PDSPCELL(recon->HA,HA);
	info2("Generating HA ");TIC;tic;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.hs->p[ifit];
	    for(int idm=0; idm<ndm; idm++){
		const double ht=parms->dm[idm].ht;
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax->p[ifit]*ht;
		displace[1]=parms->fit.thetay->p[ifit]*ht;
		loc_t *loc=recon->floc;
		if(parms->recon.misreg_dm2sci && parms->recon.misreg_dm2sci[ifit+idm*nfit]){
		    loc=loctransform(loc, parms->recon.misreg_dm2sci[ifit+idm*nfit]);
		}
		HA[idm][ifit]=mkh(recon->aloc->p[idm], loc, NULL,
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
	writebin(recon->HA,"%s/HA",dirsetup);
    }
    recon->actcpl=genactcpl(recon->HA, recon->W1);
    //if(1){//new 
    //cpl accounts for floating actuators, but not stuck actuators.
    act_stuck(recon->aloc, recon->actcpl, recon->actfloat);
    //Do not modify HA by floating actuators, otherwise, HA*actinterp will not work.
    act_stuck(recon->aloc, recon->HA, recon->actstuck);
    if(parms->save.setup){
	writebin(recon->HA,"%s/HA_float",dirsetup);
    }
    if(parms->fit.actinterp){
	recon->actinterp=act_extrap(recon->aloc, recon->actcpl, parms->fit.actthres);
    }else if(recon->actfloat){
	warning("There are float actuators, but fit.actinterp is off\n");
    }
    if(recon->actinterp){
	/*
	  DM fitting output a is extrapolated to edge actuators by
	  actinterp*a. The corresponding ray tracing from DM would be
	  HA*actinterp*a. We replace HA by HA*actinterp to take this into
	  account during DM fitting.
	 */
	warning2("Replacing HA by HA*actinterp\n");
	
	dspcell *HA2=0;
	dcellmm(&HA2, recon->HA, recon->actinterp, "nn", 1);
	dspcellfree(recon->HA);
	recon->HA=HA2;
	if(parms->save.setup){
	    writebin(recon->HA,"%s/HA_final",dirsetup);
	}
    }
    
    if(parms->save.setup){
	if(recon->actinterp){
	    writebin(recon->actinterp, "%s/actinterp", dirsetup);
	}
	if(recon->actcpl){
	    writebin(recon->actcpl, "actcpl");
	}
    }
}

/**
   Setup fitting low rank terms that are in the NULL space of DM fitting
   operator. typically include piston on each DM and tip/tilt on certain
   DMs. Becareful with tip/tilt contraint when using CBS.  */
static void 
fit_prep_lrt(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    if(ndm>=3) warning("Low rank terms for 3 or more dms are not tested\n");
    recon->fitNW=cellnew(ndm,1);
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
    dcell* actcpl=dcelldup(recon->actcpl);
    //include stuck actuator
    act_stuck(recon->aloc, actcpl, recon->actstuck);
    for(int idm=0; idm<ndm; idm++){
	int nloc=recon->aloc->p[idm]->nloc;
	recon->fitNW->p[idm]=dnew(nloc, nnw);
    }
    int inw=0;/*current column */
    if(parms->fit.lrt_piston){
	info2("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc->p[idm]->nloc;
	    double *p=recon->fitNW->p[idm]->p+(inw+idm)*nloc;
	    const double *cpl=actcpl->p[idm]->p;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(cpl[iloc]>0.1){
		    //don't count floating or stuck actuators
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
		int nloc=recon->aloc->p[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+(inw+(idm-1)*2)*nloc;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.1){
			p2x[iloc]=recon->aloc->p[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc->p[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }
	}else if(lrt_tt==2){/*Canceling TT. only valid for 2 DMs */
	    warning("Adding Canceling TT cr to fit matrix. Deprecated\n");
	    if(ndm!=2){
		error("Only ndm=2 case is implemented\n");
	    }
	    for(int idm=0; idm<ndm; idm++){
		int nloc=recon->aloc->p[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+inw*nloc;
		if(idm==0) factor=scl*2/parms->aper.d;
		else if(idm==1) factor=-scl*2./parms->aper.d;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.1){
			p2x[iloc]=recon->aloc->p[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc->p[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }

	}
	inw+=2*(ndm-1);
    }
    cellfree(actcpl);
    if(parms->fit.actslave){
	/*
	  2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found
	  that slaving is causing mis-measurement of a few edge
	  actuators. First try to remove W1. Or lower the weight. Revert
	  back.
	  1./floc->nloc is on the same order of norm of Ha'*W*Ha. 
	*/
	TIC;tic;
	recon->actslave=slaving(recon->aloc, recon->actcpl,
				recon->fitNW, recon->actstuck,
				recon->actfloat, parms->fit.actthres, 1./recon->floc->nloc);
	toc2("slaving");
	if(parms->save.setup){
	    writebin(recon->actslave,"%s/actslave",dirsetup);
	}
    }
    if(parms->save.setup){
	writebin(recon->fitNW,"%s/fitNW",dirsetup);
    }
}
/**
   Create initial reconstruction parameters so we can do NCPA calibration.
 */
RECON_T *setup_recon_init(const PARMS_T *parms, const APER_T *aper){
    RECON_T * recon = calloc(1, sizeof(RECON_T)); 
    if(parms->recon.warm_restart){
	info2("Using warm restart\n");
    }else{
	warning2("Do not use warm restart\n");
    }
    /*to be used in tomography. */
    recon->nthread=NTHREAD;
    /*for recon->aloc dimension*/
    recon->ndm=parms->ndm;
    /*setup DM actuator grid */
    setup_recon_aloc(recon,parms,aper);
    /*Grid for DM fitting*/
    setup_recon_floc(recon,parms, aper);
    if(parms->recon.alg==0 || parms->sim.ncpa_calib){
	setup_recon_HA(recon,parms);
	fit_prep_lrt(recon, parms);
    }
    return recon;
}
