/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "recon.h"
#include "fdpcg.h"
#include "ahst.h"
#include "cn2est.h"
#include "recon_utils.h"
#include "moao.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file setup_recon.c Contains routines that setup the wavefront reconstructor
   and DM fitting.  Use parms->wfsr instead of parms->wfs for wfs information,
   which hands GLAO mode correctly.x */
/**
   Setting up PLOC grid, which is a coarse sampled (similar to subaperture
   spacing) grid that defines the circular aperture for wavefront reconstruction
   and DM fitting.*/
static void
setup_recon_ploc(RECON_T *recon, const PARMS_T *parms){
    if(recon->ploc){//free the old one in case of repeated call
	locfree(recon->ploc); recon->ploc=NULL;
    }
    if(parms->load.ploc){//optionally load ploc from the file. see dbg.conf
	char *fn=parms->load.ploc;
	warning("Loading ploc from %s\n",fn);
	recon->ploc=locread("%s",fn);
	recon->pmap=loc2map(recon->ploc);
	free(recon->pmap->p);
	recon->pmap->p=NULL;
    }else{ 
	/*
	  Create a circular PLOC with telescope diameter by calling
	  create_metapupil with height of 0. We don't add any guard points.*/
	double dxr=parms->atmr.dx/parms->tomo.pos;//sampling of ploc
	map_t *pmap=create_metapupil_wrap 
	    (parms,0,dxr,0,0,0,0,parms->fit.square);
	info2("PLOC is %ldx%ld, with sampling of %.2fm\n",pmap->nx,pmap->ny,dxr);
	recon->ploc=map2loc(pmap);//convert map_t to loc_t
	recon->pmap = pmap;
	free(pmap->p); pmap->p=NULL;
	if(parms->save.setup){
	    locwrite(recon->ploc, "%s/ploc",dirsetup);
	}
    }
    loc_create_stat(recon->ploc);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].lo){//this is a low order wfs
	    if(parms->powfs[ipowfs].gtype_recon==0 
	       || parms->powfs[ipowfs].gtype_sim==0){
		//we are using gtilt for low order wfs
		recon->lowfs_gtilt=1;
		if(parms->powfs[ipowfs].gtype_sim==0)
		    warning("Low order POWFS %d is using gtilt in sim. "
			    "This is not recommended\n",ipowfs);
	    }
	}
    }
    //create the weighting W for bilinear influence function. See [Ellerbroek 2002]
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
	mkw_annular(recon->ploc, 0, 0, rin, parms->aper.d/2,
		    &(recon->W0), &(recon->W1));
	if(parms->save.setup){
	    spwrite(recon->W0, "%s/W0",dirsetup);
	    dwrite(recon->W1, "%s/W1",dirsetup);
	}
    }
    
    if(parms->plot.setup){//plot the ploc grid.
	plotloc("FoV",parms,recon->ploc,0, "ploc");
    }
}
/**
   Setup the deformable mirrors grid aloc. This is used for DM fitting.
*/
static void
setup_recon_aloc(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    if(ndm==0) return;
    if(recon->aloc){
	locarrfree(recon->aloc, ndm); recon->aloc=NULL;
    }
    recon->amap=calloc(ndm, sizeof(map_t*));
    recon->aembed=calloc(ndm, sizeof(long*));
    if(parms->load.aloc){
	char *fn=parms->load.aloc;
	warning("Loading aloc from %s\n",fn);
	int naloc;
	recon->aloc=locarrread(&naloc,"%s",fn);
	if(naloc!=ndm) error("Invalid saved aloc");
	for(int idm=0; idm<parms->ndm; idm++){
	    recon->amap[idm]=loc2map(recon->aloc[idm]);
	    recon->aembed[idm]=map2embed(recon->amap[idm]);
	    free(recon->amap[idm]->p); recon->amap[idm]->p=NULL;
	}
    }else{
	recon->aloc=calloc(ndm, sizeof(loc_t*));
	//int nxmax=0, nymax=0;
	for(int idm=0; idm<ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dx=parms->dm[idm].dx;
	    double offset=parms->dm[idm].offset+(parms->dm[idm].order%2)*0.5;
	    const double guard=parms->dm[idm].guard*parms->dm[idm].dx;
	    
	    map_t *map=create_metapupil_wrap
		(parms,ht,dx,offset,guard,0,0,parms->fit.square);
	    info("Dm %d: map is %ld x %ld\n", idm, map->nx, map->ny);
	    recon->aloc[idm]=map2loc(map);
	    recon->amap[idm]=map;
	    recon->aembed[idm]=map2embed(map);
	    free(map->p); map->p=NULL; 
	    free(recon->amap[idm]->nref);recon->amap[idm]->nref=NULL;
	}//idm
        if(parms->save.setup){
	    locarrwrite(recon->aloc,parms->ndm,"%s/aloc",dirsetup);
	}
    }
    if(parms->plot.setup){
	for(int idm=0; idm<ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    plotloc("FoV", parms, recon->aloc[idm], ht, "aloc%d", idm);
	}
    }
    recon->alocm=calloc(ndm, sizeof(loc_t*));
    for(int idm=0; idm<ndm; idm++){
	if(parms->dm[idm].misreg){
	    warning("Creating misregistration for DM.\n");
	    dcell *misreg=dcellread("%s",parms->dm[idm].misreg); 
	    PDCELL(misreg, pmisreg);
	    if(misreg->nx!=2 || misreg->ny!=1)
		error("%s is in wrong format\n",parms->dm[idm].misreg);
	    recon->alocm[idm]=loctransform(recon->aloc[idm],NULL,pmisreg[0]);
	    if(parms->plot.setup){
		plotloc("FoV", parms, recon->alocm[idm], parms->dm[idm].ht, "alocm%d", idm);
	    }
	}else{
	    recon->alocm[idm]=recon->aloc[idm];
	}
    }
    if(parms->save.setup){
	locarrwrite(recon->alocm,parms->ndm,"%s/alocm",dirsetup);
    }
    recon->aimcc=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc[idm], NULL);
	dinvspd_inplace(recon->aimcc->p[idm]);
    }
    //Dealing with stuck/floating actuators.
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
	recon->actinterp=act_float_interp(recon->aloc, recon->actfloat);
	if(parms->save.setup){
	    spcellwrite(recon->actinterp, "%s/actinterp", dirsetup);
	}
    }
}
/**
   Setup the tomography grids xloc which is used for Tomography.
*/
static void 
setup_recon_xloc(RECON_T *recon, const PARMS_T *parms){
    const int npsr=recon->npsr;
    if(recon->xloc){
	locarrfree(recon->xloc, npsr); recon->xloc=NULL;
    }
    recon->xmap=calloc(npsr, sizeof(map_t*));
    if(parms->load.xloc){
	char *fn=parms->load.xloc;
	warning("Loading xloc from %s\n",fn);
	int nxloc;
	recon->xloc=locarrread(&nxloc,"%s",fn);
	if(nxloc!=npsr) 
	    error("Invalid saved file. npsr=%d, nxloc=%d\n",npsr,nxloc);
	for(int ips=0; ips<npsr; ips++){
	    recon->xmap[ips]=loc2map(recon->xloc[ips]);
	    free(recon->xmap[ips]->p); recon->xmap[ips]->p=NULL;
	    free(recon->xmap[ips]->nref);recon->xmap[ips]->nref=NULL;
	}
    }else{
	recon->xloc=calloc(npsr, sizeof(loc_t *));
	info2("Tomography grid is %ssquare:\n", parms->tomo.square?"":"not ");
	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    double dxr=(parms->sim.idealfit)?parms->atm.dx:recon->dx->p[ips];
	    const double guard=parms->tomo.guard*dxr;
	    long nin=0;
	    if(!parms->sim.idealfit && (parms->tomo.precond==1 || parms->tomo.square==2)){
		//FFT in FDPCG prefers power of 2 dimensions.
		nin=nextpow2((long)round(parms->aper.d/recon->dx->p[0]*2.))
		    *recon->os->p[ips]/recon->os->p[0];
	    }
	    map_t *map=create_metapupil_wrap
		(parms,ht,dxr,0,guard,nin,0,parms->tomo.square);
	    info2("layer %d: xloc map is %3ld x %3ld, sampling is %.3f m\n",
		  ips, map->nx,map->ny,dxr);
	    recon->xloc[ips]=map2loc(map);
	    recon->xmap[ips]=map;
	    //Free the data and nref so we are assign pointers to it.
	    free(recon->xmap[ips]->p);recon->xmap[ips]->p=NULL;
	    free(recon->xmap[ips]->nref);recon->xmap[ips]->nref=NULL;
	}
	if(parms->save.setup){
	    locarrwrite(recon->xloc, recon->npsr, "%s/xloc",dirsetup);
	}
    }
    if(parms->plot.setup){
	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    plotloc("FoV",parms,recon->xloc[ips],ht, "xloc%d",ips);
	}
    }
    
    recon->xmcc=dcellnew(npsr,1);
    for(int ipsr=0; ipsr<npsr; ipsr++){
	recon->xmcc->p[ipsr]=loc_mcc_ptt(recon->xloc[ipsr],NULL);
	dinvspd_inplace(recon->xmcc->p[ipsr]);
    }
}
/**
   Setup ray tracing operator from xloc to ploc for guide stars
*/
static void 
setup_recon_HXW(RECON_T *recon, const PARMS_T *parms){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    if(parms->load.HXW){
	warning2("Loading saved HXW\n");
	recon->HXW=spcellread("%s",parms->load.HXW);
	if(recon->HXW->nx!=nwfs || recon->HXW->ny!=npsr){
	    error("Wrong saved HXW\n");
	}
	PDSPCELL(recon->HXW,HXW);
	int nploc=ploc->nloc;
	for(int ips=0; ips<npsr; ips++){
	    int nloc=recon->xloc[ips]->nloc;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(!HXW[ips][iwfs] 
		   || HXW[ips][iwfs]->m!=nploc 
		   || HXW[ips][iwfs]->n!=nloc){
		    error("Wrong saved HXW\n");
		}
	    }
	}
    }else{
	info2("Generating HXW");TIC;tic;
	recon->HXW=spcellnew(nwfs,npsr);
	PDSPCELL(recon->HXW,HXW);
    	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->tomo.split!=2 && parms->powfs[ipowfs].skip){
		//don't need HXW for low order wfs that does not participate in tomography.
		continue;
	    }
	    double  hs = parms->powfs[ipowfs].hs;
	    for(int ips=0; ips<npsr; ips++){
		double  ht = recon->ht->p[ips];
		double  scale=1. - ht/hs;
		double  displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht;
		displace[1]=parms->wfsr[iwfs].thetay*ht;
		HXW[ips][iwfs]=mkh(recon->xloc[ips], ploc, NULL, 
				   displace[0],displace[1],scale,
				   parms->tomo.cubic, parms->tomo.iac);
		
	    }
	}
	toc2(" ");
	if(parms->save.setup){
	    spcellwrite(recon->HXW, "%s/HXW",dirsetup);
	}
    }
    spcellfree(recon->HXWtomo);
    recon->HXWtomo=spcellnew(recon->HXW->nx, recon->HXW->ny);
    PDSPCELL(recon->HXWtomo,HXWtomo);
    PDSPCELL(recon->HXW,HXW);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){//for tomography
	    for(int ips=0; ips<npsr; ips++){
		HXWtomo[ips][iwfs]=spref(HXW[ips][iwfs]);
	    }
	} 
    }
}
/**
   Setup gradient operator from powfs.loc to wavefront sensor for reconstruction
if necessary.  */
static void
setup_recon_GWR(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    if(!parms->dbg.usegwr) return;
    recon->GWR=spcellnew(parms->npowfs, 1);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].gtype_recon==0){
	    if(!powfs[ipowfs].locm && powfs[ipowfs].GS0){
		recon->GWR->p[ipowfs] = spref(powfs[ipowfs].GS0->p[0]);
	    }else{
		double displace[2]={0,0};
		recon->GWR->p[ipowfs] = mkg(powfs[ipowfs].loc, powfs[ipowfs].loc,
					    powfs[ipowfs].amp->p, powfs[ipowfs].saloc,
					    1, 1, displace, 1);
	    }
	}else{
	    double displace[2]={0,0};
	    recon->GWR->p[ipowfs] = mkz(powfs[ipowfs].loc,powfs[ipowfs].amp->p,
					(loc_t*)powfs[ipowfs].pts, 1,1,displace);
	}
    }
}
/**
   Setup gradient operator from ploc to wavefront sensors.
*/
static void
setup_recon_GP(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    spcell *GP=NULL;
    if(parms->load.GP){
	warning2("Loading saved GP\n");
	GP=spcellread("%s",parms->load.GP);
	int nploc=ploc->nloc;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    int nsa=powfs[ipowfs].pts->nsa;
	    if(GP->p[ipowfs]->m!=nsa*2 || GP->p[ipowfs]->n!=nploc){
		error("Wrong saved GP\n");
	    }
	}
    }else{
	GP=spcellnew(parms->npowfs, 1);
	info2("Generating GP with ");TIC;tic;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs==0) continue;
	    /*use ploc as an intermediate plane.  Amplitude must use assumed amp
	      (non-misregistered)*/
	    
	    switch(parms->powfs[ipowfs].gtype_recon){
	    case 0:{ /*Create averaging gradient operator (gtilt) from PLOC,
		       using fine sampled powfs.loc as intermediate plane*/
		double displace[2]={0,0};
		info2(" Gploc");
		GP->p[ipowfs]=mkg(ploc,powfs[ipowfs].loc,powfs[ipowfs].amp->p,
				  powfs[ipowfs].saloc,1,1,displace,1);
		if(parms->save.setup){
		    spwrite(GP->p[ipowfs], "%s/powfs%d_GP", dirsetup, ipowfs);
		}
	    }
		break;
	    case 1:{ /*Create ztilt operator from PLOC, using fine sampled
		       powfs.loc as intermediate plane.*/
		double displace[2]={0,0};
		dsp* ZS0=mkz(powfs[ipowfs].loc,powfs[ipowfs].amp->p,
			     (loc_t*)powfs[ipowfs].pts, 1,1,displace);
		info2(" Zploc");
		dsp *H=mkh(ploc,powfs[ipowfs].loc,powfs[ipowfs].amp->p, 0,0,1,0,0);
		GP->p[ipowfs]=spmulsp(ZS0,H);
		spfree(H);
		spfree(ZS0);
	    }
		break;
	    default:
		error("Invalid gtype_recon\n");
	    }
	}
	toc2(" ");
	if(parms->save.setup){
	    spcellwrite(GP,"%s/GP",dirsetup);
	}
    }
    //assign GP for powfs to recon->GP for each wfs
    spcellfree(recon->GP);
    recon->GP=spcellnew(nwfs,1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	recon->GP->p[iwfs]=spref(GP->p[ipowfs]);
    }
    spcellfree(GP);//assigned to recon->GP already;
}
/**
   Setup gradient operator form aloc for wfs by using GP.
*/
static void
setup_recon_GA(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    const int ndm=parms->ndm;
    spcellfree(recon->GA);
    if(parms->load.GA){
	warning2("Loading saved GA\n");
	recon->GA=spcellread("GA");
	if(recon->GA->nx!=nwfs || recon->GA->ny!=ndm)
	    error("Wrong saved GA\n");
	PDSPCELL(recon->GA,GA);
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc[idm]->nloc;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs = parms->wfsr[iwfs].powfs;
		if(parms->sim.skysim && parms->powfs[ipowfs].lo){
		    continue;
		}
		int nsa=powfs[ipowfs].pts->nsa;
		if(GA[idm][iwfs]->m!=nsa*2 || GA[idm][iwfs]->n!=nloc){
		    error("Wrong saved GA\n");
		}
	    }
	}
    }else{
	info2("Generating GA");TIC;tic;
	recon->GA= spcellnew(nwfs, ndm);
	PDSPCELL(recon->GA,GA);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->sim.skysim && parms->powfs[ipowfs].lo){
		continue;
	    }
	    double  hs = parms->powfs[ipowfs].hs;
	    for(int idm=0; idm<ndm; idm++){
		double  ht = parms->dm[idm].ht;
		double  scale=1. - ht/hs;
		double  displace[2]={0,0};
		if(parms->sim.recon!=2){
		    displace[0]=parms->wfsr[iwfs].thetax*ht;
		    displace[1]=parms->wfsr[iwfs].thetay*ht;
		}
		if(parms->dbg.usegwr){
		    warning("todo: Fix and use mkg directly\n");
		    dsp *H=mkh(recon->aloc[idm], powfs[ipowfs].loc, NULL, 
			       displace[0],displace[1],scale,
			       parms->dm[idm].cubic,parms->dm[idm].iac);
		    GA[idm][iwfs]=spmulsp(recon->GWR->p[ipowfs], H);
		    spfree(H);
		}else{
		    dsp *H=mkh(recon->aloc[idm], ploc, NULL, 
			       displace[0],displace[1],scale,
			       parms->dm[idm].cubic,parms->dm[idm].iac);
		    
		    GA[idm][iwfs]=spmulsp(recon->GP->p[iwfs], H);
		    spfree(H);
		}
	    }//idm
	}
	if(parms->save.setup){
	    spcellwrite(recon->GA, "%s/GA",dirsetup);
	}
    	toc2(" ");
    }
    if(recon->actstuck){
	warning2("Apply stuck actuators to GA\n");
	act_stuck(recon->aloc, recon->GA, NULL, recon->actstuck);
    	if(parms->save.setup){
	    spcellwrite(recon->GA,"%s/GA_stuck",dirsetup);
	}
    }
    if(recon->actfloat){
	warning2("Apply float actuators to GA\n");
	act_float(recon->aloc, &recon->GA, NULL, recon->actfloat);
	if(parms->save.setup){
	    spcellwrite(recon->GA,"%s/GA_float",dirsetup);
	}
    }
    //Create GAlo that only contains GA for low order wfs
    spcellfree(recon->GAlo);
    recon->GAlo=spcellnew(recon->GA->nx, recon->GA->ny);
    recon->GAhi=spcellnew(recon->GA->nx, recon->GA->ny);
    PDSPCELL(recon->GAlo,GAlo);
    PDSPCELL(recon->GAhi,GAhi);
    PDSPCELL(recon->GA,GA);
    for(int idm=0; idm<ndm; idm++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].lo){//for low order wfs
		GAlo[idm][iwfs]=spref(GA[idm][iwfs]);
	    }else{
		GAhi[idm][iwfs]=spref(GA[idm][iwfs]);		
	    }
	}
    }
}
/**
   Crate the xloc to wfs gradient operator.
*/
static void 
setup_recon_GX(RECON_T *recon, const PARMS_T *parms){
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    spcellfree(recon->GX);
    recon->GX=spcellnew(nwfs, npsr);
    PDSPCELL(recon->GX,GX);
    PDSPCELL(recon->HXW,HXW);
    info2("Generating GX");TIC;tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	//gradient from xloc
	for(int ips=0; ips<npsr; ips++){
	    GX[ips][iwfs]=spmulsp(recon->GP->p[iwfs], HXW[ips][iwfs]);
	}//ips
    }
    toc2(" ");
    spcellfree(recon->GXhi);
    spcellfree(recon->GXlo);
    spcellfree(recon->GXtomo);
    spcellfree(recon->GXfocus);
    
    recon->GXtomo=spcellnew(recon->GX->nx, recon->GX->ny);
    PSPCELL(recon->GXtomo,GXtomo);

    recon->GXhi=spcellnew(recon->GX->nx, recon->GX->ny);
    PSPCELL(recon->GXhi, GXhi);

    recon->GXlo=spcellnew(recon->GX->nx, recon->GX->ny);
    PSPCELL(recon->GXlo, GXlo);

    recon->GXfocus=spcellnew(recon->GX->nx, recon->GX->ny);
    PSPCELL(recon->GXfocus,GXfocus);
    
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){//for tomography
	    for(int ips=0; ips<npsr; ips++){
		GXtomo[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	}
	if(parms->powfs[ipowfs].lo){//for low order wfs
	    for(int ips=0; ips<npsr; ips++){
		GXlo[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	 
	}else{//for high order wfs
	    for(int ips=0; ips<npsr; ips++){
		GXhi[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	}
	if(!parms->sim.idealfit && parms->sim.mffocus && parms->sim.closeloop){
	    //for focus tracking.
	    for(int ips=0; ips<npsr; ips++){
		GXfocus[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	}
    }//iwfs
}
/**
   Setup the matrix of the inverse of gradient measurement noise equivalent
   angle covariance matrix. For physical optics wfs, the NEA is computed using
   matched filter output. For geometric optics, the NEA is from input.
*/
static void
setup_recon_saneai(RECON_T *recon, const PARMS_T *parms, 
		   const POWFS_T *powfs){
    const int nwfs=parms->nwfsr;
    spcellfree(recon->sanea);
    spcellfree(recon->saneai);
    spcellfree(recon->saneal);
    spcell *sanea=recon->sanea=spcellnew(nwfs,nwfs);
    spcell *saneal=recon->saneal=spcellnew(nwfs,nwfs);
    spcell *saneai=recon->saneai=spcellnew(nwfs,nwfs);
    info2("saneai:");
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	if(parms->powfs[ipowfs].neareconfile){//taks precedance
	    dmat *nea=dread("%s_wfs%d",parms->powfs[ipowfs].neareconfile,iwfs);//rad
	    /* sanity check */
	    if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy){
		//sanity check on provided nea.
		const int nmtch=powfs[ipowfs].intstat->mtche->ny;
		int indsanea=0;
		if(nmtch==1){
		    indsanea=0;
		}else if(nmtch==parms->powfs[ipowfs].nwfs){
		    indsanea=parms->powfs[ipowfs].wfsind[iwfs];
		}else{
		    error("invalid\n");
		}
		PDCELL(powfs[ipowfs].intstat->saneaixy, saneaixy);
		//Do some error checking.
		for(int isa=0; isa<nsa; isa++){
		    double nea1x=pow(saneaixy[indsanea][isa]->p[0],-0.5);
		    double nea1y=pow(saneaixy[indsanea][isa]->p[3],-0.5);
		    double nea2x=nea->p[isa];
		    double nea2y=nea->p[isa+nsa];
		    if(nea1x<0.5*nea2x || nea1x>2*nea2x || nea1y<0.5*nea2y || nea1y>2*nea2y){
			warning2("iwfs %d, isa %d : Phy nea: %g %g mas. Provided nea: %g %g mas.\n",
				 iwfs,isa,nea1x*206265000, nea1y*206265000, 
				 nea2x*206265000, nea2y*206265000);
		    }
		}
	    }
	    //rad
	    saneal->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea, 2);//rad^2
	    sanea->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea,-1);//rad^-2
	    saneai->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dfree(nea);
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    //Physical optics
	    if(parms->powfs[ipowfs].phytype==1){
		const int nmtch=powfs[ipowfs].intstat->mtche->ny;
		if(parms->sim.glao && nmtch!=1){
		    error("Please average intstat->saneaxy for GLAO mode.\n");
		}
		int indsanea=0;
		if(nmtch==1){
		    indsanea=0;
		}else if(nmtch==parms->powfs[ipowfs].nwfs){
		    indsanea=parms->powfs[ipowfs].wfsind[iwfs];
		}else{
		    error("invalid\n");
		}
		long ldx=powfs[ipowfs].intstat->saneaxy->nx;
		dmat **saneaxy=powfs[ipowfs].intstat->saneaxy->p+indsanea*ldx;
		dmat **sanealxy=powfs[ipowfs].intstat->saneaxyl->p+indsanea*ldx;
		dmat **saneaixy=powfs[ipowfs].intstat->saneaixy->p+indsanea*ldx;
		sanea->p[iwfs+iwfs*nwfs]=nea2sp(saneaxy,nsa);
		saneal->p[iwfs+iwfs*nwfs]=nea2sp(sanealxy,nsa);
		saneai->p[iwfs+iwfs*nwfs]=nea2sp(saneaixy,nsa);
	    }else{
		error("Not implemented yet\n");
	    }
	}else{
	    //nea scales as sqrt(1/dtrat)
	    const double neasq=pow(parms->powfs[ipowfs].nearecon/206265000.,2)/parms->powfs[ipowfs].dtrat;
	    if(neasq<1.e-30) error("nea is too small\n");
	    dmat *nea=dnew(nsa,2);
	    //scale neasq by area^-1. (seeing limited)
	    //scale neasq by area^-2 if diffraction limited
	    //only implementing seeing limited here.
	    PDMAT(nea, neap);
	    double *area=powfs[ipowfs].saa->p;
	    for(int i=0; i<nsa; i++){
		double tmp=neasq/area[i];
		neap[0][i]=neap[1][i]=tmp;
	    }
	    sanea->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea, -1);
	    saneai->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea, -0.5);
	    saneal->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dfree(nea);
	}
    }//iwfs
    
    //Compute the averaged SANEA for each WFS
    recon->neam=dnew(parms->nwfsr, 1);
    double neamhi=0; 
    int counthi=0;
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	const int nsa=powfs[ipowfs].pts->nsa;
	dmat *sanea_iwfs=spdiag(recon->sanea->p[iwfs+iwfs*parms->nwfsr]);
	double area_thres;
	if(nsa>4){
	    area_thres=0.9;
	}else{
	    area_thres=0;
	}
	double nea2_sum=0;
	int count=0;
	for(int isa=0; isa<nsa; isa++){
	    if(powfs[ipowfs].saa->p[isa]>area_thres){
		nea2_sum+=(sanea_iwfs->p[isa])+(sanea_iwfs->p[isa+nsa]);
		count++;
	    }
	}
	dfree(sanea_iwfs);
	recon->neam->p[iwfs]=sqrt(nea2_sum/count/2);//average sanea in radian
	char *neatype;
	if(parms->powfs[ipowfs].neareconfile){
	    neatype="FILE";
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    neatype="mtch";
	}else{
	    neatype="geom";
	}
	info2(" %s(%.2f)", neatype, recon->neam->p[iwfs]*206265000);
	if(!parms->powfs[ipowfs].lo){
	    neamhi+=pow(recon->neam->p[iwfs],2);
	    counthi++;
	}
    }
    recon->neamhi=sqrt(neamhi/counthi);
    info2("\n");

    if(parms->save.setup){
	spcellwrite(recon->sanea, "%s/sanea",dirsetup);
	spcellwrite(recon->saneal,"%s/saneal",dirsetup);
	spcellwrite(recon->saneai,"%s/saneai",dirsetup);
    }
   
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(powfs[ipowfs].intstat){
	    dcellfree(powfs[ipowfs].intstat->sanea);
	    dcellfree(powfs[ipowfs].intstat->saneara);
	    dcellfree(powfs[ipowfs].intstat->saneaxy);
	    dcellfree(powfs[ipowfs].intstat->saneaxyl);
	    dcellfree(powfs[ipowfs].intstat->saneaixy);
	}
    }
}
/**
   setting up global tip/tilt remove operator from LGS gradients.
*/

static void
setup_recon_TTR(RECON_T *recon, const PARMS_T *parms, 
		const POWFS_T *powfs){
    if(!recon->has_ttr) return;
    if(recon->TT){
	dcellfree(recon->TT);
	dcellfree(recon->PTT);
    }
    int nwfs=parms->nwfsr;
    recon->TT=dcellnew(nwfs,nwfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].trs){
	    info2("powfs %d has tip/tilt removed in tomography\n", ipowfs);
	    if(parms->powfs[ipowfs].skip){
		error("This POWFS %d should be included in Tomo.\n", ipowfs);
	    }
	    int nsa=powfs[ipowfs].pts->nsa;
	    dmat *TT=dnew(nsa*2,2);
	    double *TTx=TT->p;
	    double *TTy=TT->p+nsa*2;
	    for(int isa=0; isa<nsa; isa++){
		TTx[isa]=1;
		TTy[isa]=0;
	    }
	    for(int isa=nsa; isa<nsa*2; isa++){
		TTx[isa]=0;
		TTy[isa]=1;
	    }
	    if(parms->sim.glao){
		recon->TT->p[ipowfs*(parms->npowfs+1)]=ddup(TT);
	    }else{
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs[jwfs];
		    recon->TT->p[iwfs*(1+parms->nwfsr)]=ddup(TT);
		}
	    }
	    dfree(TT);
	}
    }
    recon->PTT=dcellpinv(recon->TT,NULL,recon->saneai);
    if(parms->save.setup){
	dcellwrite(recon->TT, "%s/TT",dirsetup);
	dcellwrite(recon->PTT, "%s/PTT",dirsetup);
    }
}
/**
   operator to remove diffrential focus modes that might be caused by sodium layer
   horizontal structure.
*/

static void
setup_recon_DFR(RECON_T *recon, const PARMS_T *parms, 
		const POWFS_T *powfs){

    if(!recon->has_dfr) return;
    if(recon->DF){ 
	dcellfree(recon->DF);
	dcellfree(recon->PDF);
    }
    int nwfs=parms->nwfsr;
    recon->DF=dcellnew(nwfs,nwfs);
    //Then differential focus modes.
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].dfrs){
	    warning("powfs %d has differential focus removed in tomography\n", ipowfs);
	    if(parms->powfs[ipowfs].nwfs<2){
		error("This powfs group has only 1 wfs. Could not remove diff focus\n");
	    }
	    int nsa=powfs[ipowfs].pts->nsa;
	    dmat* DF=dnew(nsa*2,1);
	    //saloc is lower left corner of subaperture. don't have to be the center.
	    memcpy(DF->p, powfs[ipowfs].saloc->locx, sizeof(double)*nsa);
	    memcpy(DF->p+nsa, powfs[ipowfs].saloc->locy, sizeof(double)*nsa);
	    /**
	       postive focus on first wfs. negative focus on diagnonal wfs.
	    */
	    for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs[jwfs];
		if(parms->powfs[ipowfs].skip){
		    error("This WFS %d should be included in Tomo.\n", iwfs);
		}
		dcp(&recon->DF->p[iwfs*nwfs], DF);
		dadd(&recon->DF->p[iwfs+iwfs*nwfs], 0, DF, -1);
	    }
	    dfree(DF);
	}
    }
    recon->PDF=dcellpinv(recon->DF, NULL, recon->saneai);
    if(parms->save.setup){
	dcellwrite(recon->DF, "%s/DF",dirsetup);
	dcellwrite(recon->PDF, "%s/PDF",dirsetup);
    }
}
/**
   wrapps setup_recon_TTR() and setup_recon_DFR() to removal global tip/tilt and
   differential focus.
*/
static void
setup_recon_TTFR(RECON_T *recon, const PARMS_T *parms, 
		 const POWFS_T *powfs){
  
    setup_recon_TTR(recon,parms,powfs);
    setup_recon_DFR(recon,parms,powfs);
    if(recon->DF){
	recon->TTF=dcellcat(recon->TT,recon->DF,2);
    }else{
	recon->TTF=dcellref(recon->TT);
    }
    recon->PTTF=dcellpinv(recon->TTF, NULL, recon->saneai);
    if(parms->save.setup){
	dcellwrite(recon->TTF, "%s/TTF",dirsetup);
	dcellwrite(recon->PTTF, "%s/PTTF",dirsetup);
    }
    //dcellfree(recon->DF);//don't free DF to use in PDF.
    //Keep TT, PTT, used in uplink pointing.
}
/**
   Frees recon->invpsd or recon->fractal
*/
static void free_cxx(RECON_T *recon){
    if(recon->invpsd){
	dcellfree(recon->invpsd->invpsd);
	ccellfree(recon->invpsd->fftxopd);
	free(recon->invpsd);
	recon->invpsd=NULL;
    }
    if(recon->fractal){
	dcellfree(recon->fractal->xopd);
	free(recon->fractal);
	recon->fractal=NULL;
    }
}
/**
   Prepares for tomography
*/
void
setup_recon_tomo_prep(RECON_T *recon, const PARMS_T *parms){
    //Free existing struct if already exist. 
    free_cxx(recon);
    if(parms->tomo.assemble){
	//We need the old copy of L2 when we update the turbulence profile.
	spcellfree(recon->L2save);
	recon->L2save=recon->L2;
    }else{
	spcellfree(recon->L2);
    }
    recon->L2=NULL;
    /*When layers get a weight less than 1%, we put it at 1% to avoid
      regularization unstability issues.*/
    dclip(recon->wt, 0.01, 1);
    //normalize the weights to sum to 1.
    normalize(recon->wt->p, recon->npsr, 1);
    const int npsr=recon->npsr;
    recon->cxx=parms->tomo.cxx;
    //test_cxx(recon, parms);
    switch(parms->tomo.cxx){
    case 0:
	if(parms->load.cxx){
	    recon->L2=spcellread("%s",parms->load.cxx);
	    if(recon->L2->nx!=npsr || recon->L2->ny!=npsr){
		error("Wrong format of loaded L2\n");
	    }
	}else{
	    recon->L2=spcellnew(npsr,npsr);
	    for(int ips=0; ips<npsr; ips++){
		if(parms->tomo.square){//periodic bc
		    recon->L2->p[ips+npsr*ips]=mklaplacian_map
			(recon->xmap[ips]->nx, recon->xmap[ips]->nx,
			 recon->xloc[ips]->dx, recon->r0,
			 recon->wt->p[ips]);
		}else{//reflecive bc
		    recon->L2->p[ips+npsr*ips]=mklaplacian_loc
			(recon->xloc[ips], recon->r0, 
			 recon->wt->p[ips]);
		}
	    }
	    if(fabs(parms->tomo.cxxscale-1)>EPS){
		spcellscale(recon->L2, parms->tomo.cxxscale);
	    }
	    if(parms->save.setup){
		spcellwrite(recon->L2, "%s/L2",dirsetup);
	    }
	}
	break;
    case 1:{
	recon->invpsd=calloc(1, sizeof(INVPSD_T));
	if(parms->load.cxx){
	    recon->invpsd->invpsd=dcellread("%s",parms->load.cxx);
	    if(recon->invpsd->invpsd->nx!=npsr || recon->invpsd->invpsd->ny!=1){
		error("Wrong format of loaded invpsd\n");
	    }
	}else{
	    dcell* invpsd=recon->invpsd->invpsd=dcellnew(npsr,1);
	    for(int ips=0; ips<npsr; ips++){
		long nx=recon->xmap[ips]->nx;
		long ny=recon->xmap[ips]->ny;
		double r0i=recon->r0*pow(recon->wt->p[ips],-3./5.);
		invpsd->p[ips]=turbpsd(nx, ny, recon->xloc[ips]->dx, r0i,
				       recon->l0,-1);
		dscale(invpsd->p[ips], pow((double)(nx*ny),-2)*parms->tomo.cxxscale);
	    }
	    if(parms->save.setup){
		dcellwrite(invpsd, "%s/invpsd",dirsetup);
	    }
	}
	ccell* fftxopd=recon->invpsd->fftxopd=ccellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    fftxopd->p[ips]=cnew(recon->xmap[ips]->nx, recon->xmap[ips]->ny);
	    cfft2plan(fftxopd->p[ips],-1);
	    cfft2plan(fftxopd->p[ips],1);
	}
	recon->invpsd->xloc = recon->xloc;
	recon->invpsd->square = parms->tomo.square;
    }
	break;
    case 2:{
	recon->fractal=calloc(1, sizeof(FRACTAL_T));
	recon->fractal->xloc=recon->xloc;
	recon->fractal->r0=parms->atmr.r0;
	recon->fractal->l0=parms->atmr.l0;
	recon->fractal->wt=parms->atmr.wt;
	recon->fractal->scale=parms->tomo.cxxscale;
	recon->fractal->ninit=parms->tomo.ninit;
	dcell *xopd=recon->fractal->xopd=dcellnew(npsr, 1);
	for(int ips=0; ips<npsr; ips++){
	    int nn=nextpow2(MAX(recon->xmap[ips]->nx, recon->xmap[ips]->ny))+1;
	    xopd->p[ips]=dnew(nn,nn);
	}
    }
	break;
    default:
	error("tomo.cxx=%d is invalid\n", parms->tomo.cxx);
    }
   
    if(parms->tomo.piston_cr){
	/*when add constraint, make sure the order of
	  magnitude are at the same range.*/
	spcellfree(recon->ZZT);
	recon->ZZT=spcellnew(npsr,npsr);
	for(int ips=0; ips<npsr; ips++){
	    double r0=recon->r0;
	    double dx=recon->xloc[ips]->dx;
	    double wt=recon->wt->p[ips];
	    double val=pow(laplacian_coef(r0,wt,dx),2);
	    //info("Scaling of ZZT is %g\n",val);
	    //piston mode eq 47 in Brent 2002 paper
	    int icenter=loccenter(recon->xloc[ips]);
	    int nloc=recon->xloc[ips]->nloc;
	    dsp *ZZT=recon->ZZT->p[ips+npsr*ips]
		=spnew(nloc,nloc,1);
	    int icol;
	    int count=0;
	    for(icol=0; icol<nloc; icol++){
		ZZT->p[icol]=count;
		if(icol==icenter){
		    ZZT->i[count]=icenter;
		    ZZT->x[count]=val;
		    count++;
		}
	    }
	    ZZT->p[nloc]=count;
	}
	if(parms->save.setup){
	    spcellwrite(recon->ZZT, "%s/ZZT",dirsetup);
	}
    }
}
/**
   assemble tomography matrix. In CG mode, this function is not executed if
   tomo.assemble=0, Instead, the algorithm is contained in recon.c. When you
   modify anything, make sure you also do it there.  

   For integrated tomograhy: 
   
   \f$\hat{x}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1}+G_{ngs}^{T}C_{ngs}^{-1}
   G_{ngs})^{-1}(G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}+G_{ngs}^{T}C_{ngs}^{-1}s_{ngs}){\equiv}RL^{-1}RR s.\f$

   For split tomography, the terms regarding NGS are dropped:

   \f$\hat{x}_{lgs}=(G_{lgs}^{T}C_{lgs}^{-1}G_{lgs}+C_{x}^{-1})^{-1}
   G_{lgs}^{T}C_{lgs}^{-1}s_{lgs}\equiv R_L^{-1} R_R s.\f$

   The left hand side of linear equation (inside the inverse) is stored in RL.
   The right hand side of the linear equation (outside of the inverse) is stored
   in RR. The terms regarding the NGS are handled using low rank terms. The
   gradients from LGS and NGS and concatenated to \f$s\f$.

   In the LGS part, there is a global tip/tilt removal operator because LGS is
   insensitive to global tip/tilts.

   For details see http://www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803

*/
void setup_recon_tomo_matrix(RECON_T *recon, const PARMS_T *parms, APER_T *aper){
    //if not cg or forced, build explicitly the tomography matrix.
    int npsr=recon->npsr;
    int nwfs=parms->nwfsr;
    //Free OLD matrices if any.
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    info2("Before assembling tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);

    if(parms->load.tomo){
	//right hand side.
	warning("Loading saved recon->RR\n");
	recon->RR.M=spcellread("RRM");
	if(recon->has_ttr){
	    recon->RR.U=dcellread("RRU");
	    recon->RR.V=dcellread("RRV");
	}
	//Left hand side
	warning("Loading saved recon->RL\n");
	recon->RL.M=spcellread("RLM");
	recon->RL.U=dcellread("RLU");
	recon->RL.V=dcellread("RLV");
	if(parms->tomo.alg==0 && zfexist("RLC")){
	    recon->RL.C=chol_read("RLC");
	}
	if(parms->tomo.alg==2 && zfexist("RLMI")){
	    recon->RL.MI=dread("RLMI");
	}
    }else{
	info2("Building recon->RR\n");
	PDSPCELL(recon->GX,GX);
	const spcell *saneai=recon->saneai;
	/*
	  Reconstruction Right hand side matrix. In split tomography mode, low
	  order NGS are skipped. recon->GXtomo contains GXs that only
	  participate in tomography.
	*/
	spcell *GXtomoT=spcelltrans(recon->GXtomo);
	recon->RR.M=spcellmulspcell(GXtomoT, saneai, 1);
	PDSPCELL(recon->RR.M, RRM);
	/*
	  Tip/tilt and diff focus removal low rand terms for LGS WFS.
	*/
	if(recon->TTF){
	    spcellmulmat(&recon->RR.U, recon->RR.M, recon->TTF, 1);
	    recon->RR.V=dcelltrans(recon->PTTF);
	}
 
	info2("Building recon->RL\n"); //left hand side matrix
	recon->RL.M=spcellmulspcell(recon->RR.M,recon->GXhi,1);
	PDSPCELL(recon->RL.M,RLM);
	if(parms->tomo.piston_cr){ 
	    /*single point piston constraint. no need tikholnov.*/
	    info2("Adding ZZT to RLM\n");
	    for(int ips=0; ips<npsr; ips++){
		spadd(&RLM[ips][ips], recon->ZZT->p[ips+ips*npsr]);
	    }
	}
	/*Apply tikholnov regularization.*/
	if(fabs(parms->tomo.tikcr)>1.e-15){
	    //Estimated from the Formula
	    double maxeig=pow(recon->neamhi * recon->xloc[0]->dx, -2);
	    double tikcr=parms->tomo.tikcr;
	    info2("Adding tikhonov constraint of %g to RLM\n",tikcr);
	    info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	    spcelladdI(recon->RL.M, tikcr*maxeig);
	}
	//add L2 and ZZT
	switch(parms->tomo.cxx){
	case 0://Add L2'*L2 to RL.M
	    for(int ips=0; ips<npsr; ips++){
		dsp* tmp=sptmulsp(recon->L2->p[ips+npsr*ips], recon->L2->p[ips+npsr*ips]);
		if(!tmp){
		    error("L2 is empty!!\n");
		}
		spadd(&RLM[ips][ips], tmp);
		spfree(tmp);
	    }
	    break;
	case 1://Need to apply invpsd separately
	    recon->RL.extra = recon->invpsd;
	    recon->RL.exfun = apply_invpsd;
	    break;
	case 2://Need to apply fractal separately
	    recon->RL.extra = recon->fractal;
	    recon->RL.exfun = apply_fractal;
	}

	//Symmetricize, remove values below 1e-15*max and sort RLM (optional).
	//spcellsym(recon->RL.M);

	//Low rank terms for low order wfs. Only in Integrated tomography.
	dcell *ULo=dcellnew(npsr,nwfs);
	PDCELL(ULo, pULo);
	dcell *VLo=dcellnew(npsr,nwfs);
	PDCELL(VLo, pVLo);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip){
		continue;
	    }
	    if(parms->powfs[ipowfs].lo){
		for(int ips=0; ips<npsr; ips++){
		    spfull(&pULo[iwfs][ips], RRM[iwfs][ips],-1);
		    sptfull(&pVLo[iwfs][ips], GX[ips][iwfs],1);
		}
	    }
	}
	if(!parms->tomo.split || parms->dbg.splitlrt){
	    recon->RL.U=dcellcat(recon->RR.U, ULo, 2);
	    dcell *GPTTDF=NULL;
	    sptcellmulmat(&GPTTDF, recon->GX, recon->RR.V, 1);
	    recon->RL.V=dcellcat(GPTTDF, VLo, 2);
	    dcellfree(GPTTDF);
	}else{
	    warning("Skipping RL Low rank terms in split tomography to suppress noise propagation\n");
	    warning("Skipping RL Low rank terms in split tomography to suppress noise propagation\n");
	    warning("Skipping RL Low rank terms in split tomography to suppress noise propagation\n");
	    warning("Skipping RL Low rank terms in split tomography to suppress noise propagation\n");
	}
	dcellfree(ULo);
	dcellfree(VLo);
	//Remove empty cells.
	dcelldropempty(&recon->RR.U,2);
	dcelldropempty(&recon->RR.V,2);
	dcelldropempty(&recon->RL.U,2);
	dcelldropempty(&recon->RL.V,2);

	long nll=0,nlr=0;
	if(recon->RL.U){
	    /* balance UV. may not be necessary. Just to compare well against
	       laos. */
	    double r0=recon->r0;
	    double dx=recon->xloc[0]->dx;
	    double val=laplacian_coef(r0,1,dx);//needs to be a constant
	    dcellscale(recon->RL.U, 1./val);
	    dcellscale(recon->RL.V, val);
	    /*collect statistics.*/
	    PDCELL(recon->RR.U,RRU);
	    PDCELL(recon->RL.U,RLU);
	    for(int i=0; i<recon->RR.U->ny;i++){
		if(RRU[i][0]) nlr+=RRU[i][0]->ny;
	    }
	    for(int i=0; i<recon->RL.U->ny;i++){
		if(RLU[i][0]) nll+=RLU[i][0]->ny;
	    }
	}
	info2("Tomography number of Low rank terms: %ld in RHS, %ld in LHS\n", nlr,nll);
	if(parms->save.recon){
	    spcellwrite(recon->RR.M,"%s/RRM",dirsetup);
	    dcellwrite(recon->RR.U,"%s/RRU",dirsetup);
	    dcellwrite(recon->RR.V,"%s/RRV",dirsetup);

	    spcellwrite(recon->RL.M,"%s/RLM.bin",dirsetup);//disable compression
	    dcellwrite(recon->RL.U,"%s/RLU",dirsetup);
	    dcellwrite(recon->RL.V,"%s/RLV",dirsetup); 
	}
	spcellfree(GXtomoT);
    }
    if((parms->tomo.alg==0 || parms->tomo.alg==2) && parms->tomo.bgs ){
	/* We need cholesky decomposition in CBS or MVST method. */
	muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
    }
    if(((parms->tomo.alg==0 || parms->tomo.alg==2) && !parms->tomo.bgs)
       ||((parms->sim.ecnn || (parms->tomo.split==2 && !parms->load.mvst)) && !recon->RL.C && !recon->RL.MI)){
	if(parms->load.tomo){
	    if(parms->tomo.alg==0 && zfexist("RLC")){
		recon->RL.C=chol_read("RLC");
	    }
	    if(parms->tomo.alg==2 && zfexist("RLMI")){
		recon->RL.MI=dread("RLMI");
	    }
	}
	if(!recon->RL.C && !recon->RL.MI){
	    muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    if(parms->save.recon){
       	if(recon->RL.C)
	    chol_save(recon->RL.C,"%s/RLC.bin",dirsetup);
	if(recon->RL.MI)
	    dwrite(recon->RL.MI,"%s/RLMI", dirsetup);
	if(recon->RL.CB){
	    for(int ib=0; ib<recon->RL.nb; ib++){
		chol_save(recon->RL.CB[ib],"%s/RLCB_%d.bin",dirsetup, ib);
	    }
	}
	if(recon->RL.MIB){
	    dcellwrite(recon->RL.MIB,"%s/RLMIB", dirsetup);
	}
    }
    if(parms->tomo.assemble && !parms->cn2.tomo){
	//Don't free PTT. Used in forming LGS uplink err
	if(parms->tomo.piston_cr)
	    spcellfree(recon->ZZT);
    }
  
    if(parms->tomo.assemble){
	spcellfree(recon->GP);
    }
    if(parms->tomo.assemble || parms->tomo.square){
	spcellfree(recon->HXWtomo);
    }
    if(parms->sim.ecnn){
	/**
	   We compute the wavefront estimation error covariance in science focal
	   plane due to wavefront measurement noise. Basically we compute
	   Hx*E*Cnn*E'*Hx' where E is the tomography operator, and Hx is ray
	   tracing from tomography grid xloc to science focal plane ploc. Since
	   Cnn is symmetrical and sparse, we can decompose it easily into
	   Cnn=Cnl*Cnl'; We first compute L=Hx*E*Cnl, and the result is simply
	   LL'; This is much faster than computing left and right separately,
	   because 1) the number of points in xloc is larger than in Cnn, so
	   after the tomography right hand side vector is applied, the number of
	   rows is larger than number of columns, this causes the right hand
	   side solver to be much slower. 2) Simply double the computation.

	   For HX opeation, build the sparse matrix and do multiply is way
	   slower than doing ray tracing directly. 

	   For ad hoc split tomography, we need to remove the five NGS modes
	   from here, as well as in time averaging of estimated turbulence.

	   recon->saneal contains Cnl.
	*/
	TIC;tic;
	read_self_cpu();
	
	if(0){
	    //Luc's Method. Solve E^T Hx^T
	    {
		warning("Overriding sanea\n");
		warning("Overriding sanea\n");
		warning("Overriding sanea\n");
		spcellfree(recon->sanea);
		recon->sanea=spcellread("nt");
	    }
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psfr[ievl]) continue;
		dcell *hxt=dcellnew(recon->npsr, 1);
		double hs=parms->evl.hs[ievl];
		for(int ips=0; ips<recon->npsr; ips++){
		    const double ht=recon->ht->p[ips];
		    const double scale=1.-ht/hs;
		    const double dispx=parms->evl.thetax[ievl]*ht;
		    const double dispy=parms->evl.thetay[ievl]*ht;
		    dsp *HXT=mkhb(recon->xloc[ips], aper->locs, NULL,
				 dispx, dispy, scale, 0, 0);
		    spfull(&hxt->p[ips], HXT, 1); spfree(HXT);
		}
		dcell *t1=NULL;
		info("CPU Usage: %.1f HXT   ", read_self_cpu()); toc2(" ");tic;
		muv_direct_solve_cell(&t1, &recon->RL, hxt); dcellfree(hxt);
		info("CPU Usage: %.1f Solve ", read_self_cpu()); toc2(" ");tic;
		dcell *p=NULL;
		muv_t(&p, &recon->RR, t1, 1); dcellfree(t1);
		dcellwrite(p, "p_%d.bin", ievl);
		info("CPU Usage: %.1f RHS   ", read_self_cpu()); toc2(" ");tic;
		dcell *t2=NULL;
		
		spcellmulmat(&t2, recon->sanea, p, 1); 
		dcell *t3=NULL;
		dcellmm(&t3, p, t2, "tn", 1);
		dcellfree(p); dcellfree(t2);
		info("CPU Usage: %.1f MUL   ", read_self_cpu()); toc2(" ");tic;
		dwrite(t3->p[0], "ecnn_new_%d.bin", ievl);
		dcellfree(t3);
	    }
	}
	if(1){
	    dcell *rhs=NULL;
	    muv_sp(&rhs, &recon->RR, recon->saneal, 1);
	    info("CPU Usage: %.1f RHS   ", read_self_cpu()); toc2(" ");tic;
	    dmat *rhs2=dcell2m(rhs); dcellfree(rhs);
	    dwrite(rhs2, "rhs_1");
	    dmat *t1=NULL;
	    muv_direct_solve(&t1, &recon->RL, rhs2); dfree(rhs2);
	    dwrite(t1, "solve_1");
	    info("CPU Usage: %.1f Solve ", read_self_cpu()); toc2(" ");tic;
	    recon->ecnn=dcellnew(parms->evl.nevl, 1);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psfr[ievl]) continue;
		tic;
		char strht[24];
		if(!isinf(parms->evl.hs[ievl])){
		    snprintf(strht, 24, "_%g", parms->evl.hs[ievl]);
		}else{
		    strht[0]='\0';
		}
		/*Build HX for science directions that need ecov.*/
	    
		dmat *x1=dnew(aper->locs->nloc, t1->ny);
		PDMAT(t1, pt1);
		PDMAT(x1, px1);
		int ind=0;
		double hs=parms->evl.hs[ievl];
		for(int ips=0; ips<recon->npsr; ips++){
		    const double ht=recon->ht->p[ips];
		    const double scale=1.-ht/hs;
		    const double dispx=parms->evl.thetax[ievl]*ht;
		    const double dispy=parms->evl.thetay[ievl]*ht;
		    for(int icol=0; icol<t1->ny; icol++){
			prop_nongrid(recon->xloc[ips], &pt1[icol][ind], aper->locs, NULL,
				     px1[icol], 1, dispx, dispy, scale, 0, 0);
		    }
		    ind+=recon->xloc[ips]->nloc;
		}
		info("CPU Usage: %.1f accphi", read_self_cpu()); toc2(" ");tic;
	    
		dmm(&recon->ecnn->p[ievl], x1, x1, "nt", 1);
		dfree(x1);
		info("CPU Usage: %.1f Mul   ", read_self_cpu()); toc2(" ");
		dwrite(recon->ecnn->p[ievl], "ecnn_x%g_y%g%s.bin", 
		       parms->evl.thetax[ievl]*206265,
		       parms->evl.thetay[ievl]*206265, strht);
	    }
	    dfree(t1);
	}
    }
    info2("After assemble tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}
/**
   Update assembled tomography matrix with new L2.
*/
void setup_recon_tomo_update(RECON_T *recon, const PARMS_T *parms){
    setup_recon_tomo_prep(recon, parms); //redo L2, invpsd
    if(parms->tomo.alg==1&&!parms->tomo.assemble){//no need to do anything
	return;
    }
    if(parms->tomo.cxx==0){
	//Need to adjust RLM with the new L2.
	if(!recon->L2save){
	    error("We need the L2save to update the tomography matrix\n");
	}
	PDSPCELL(recon->RL.M,RLM);
	const int npsr=recon->npsr;
	for(int ips=0; ips<npsr; ips++){
	    dsp* LL=sptmulsp(recon->L2->p[ips+npsr*ips], 
			     recon->L2->p[ips+npsr*ips]);
	    dsp *LLold=sptmulsp(recon->L2save->p[ips+npsr*ips], 
				recon->L2save->p[ips+npsr*ips]);
	    if(!LL){
		error("L2 is empty!!\n");
	    }
	    dsp *LLdiff=spadd2(LL,LLold,1,-1);//adjustment to RLM
	    spadd(&RLM[ips][ips], LLdiff);
	    spfree(LLdiff);
	    spfree(LL);
	    spfree(LLold);
	}
    }

    if(parms->tomo.alg==0 || parms->tomo.alg==2){
	//We need cholesky decomposition in CBS or MVST method.
	if(!parms->tomo.bgs){//Full Matrix
	    muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}else{//BGS
	    muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }

    if(parms->tomo.split==2){
	setup_recon_mvst(recon,parms);
    }
}
/**
   Setup ray tracing operator HXF from xloc to aperture ploc along DM fiting directions*/
static void
setup_recon_HXF(RECON_T *recon, const PARMS_T *parms){
    if(parms->load.HXF && zfexist(parms->load.HXF)){
	warning("Loading saved HXF\n");
	recon->HXF=spcellread("%s",parms->load.HXF);
    }else{
	info2("Generating HXF");TIC;tic;
	const int nfit=parms->fit.nfit;
	const int npsr=recon->npsr;
	recon->HXF=spcellnew(nfit, npsr);
	PDSPCELL(recon->HXF,HXF);
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.ht[ifit];
	    for(int ips=0; ips<npsr; ips++){
		const double ht = recon->ht->p[ips];
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		HXF[ips][ifit]=mkh(recon->xloc[ips], recon->ploc, NULL,
				   displace[0], displace[1], scale,
				   parms->tomo.cubic, parms->tomo.iac);
	    }
	}
	if(parms->save.setup){
	    spcellwrite(recon->HXF,"%s/HXF",dirsetup);
	}
	toc2(" ");
    }
}
/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static void
setup_recon_HA(RECON_T *recon, const PARMS_T *parms){
    if(parms->load.HA && zfexist(parms->load.HA)){
	warning("Loading saved HA\n");
	recon->HA=spcellread("%s",parms->load.HA);
    }else{
	const int nfit=parms->fit.nfit;
	const int ndm=parms->ndm;
	recon->HA=spcellnew(nfit, ndm);
	PDSPCELL(recon->HA,HA);
	info2("Generating HA");TIC;tic;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.ht[ifit];
	    for(int idm=0; idm<ndm; idm++){
		const double ht=parms->dm[idm].ht;
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		HA[idm][ifit]=mkh(recon->aloc[idm], recon->ploc, NULL,
				  displace[0], displace[1], 
				  scale,parms->dm[idm].cubic,parms->dm[idm].iac);
	    }
	}
	toc2(" ");
	if(parms->save.setup){
	    spcellwrite(recon->HA,"%s/HA",dirsetup);
	}
    }
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
	if(parms->save.setup){
	    spcellwrite(recon->HA,"%s/HA_float",dirsetup);
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
    recon->fitNW=dcellnew(ndm,1);
    double scl=recon->fitscl;
    if(fabs(scl)<1.e-15){
	error("recon->fitscl is too small\n");
    }
    //computed outside.
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
    int inw=0;//current column
    if(parms->fit.lrt_piston){
	info2("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc[idm]->nloc;
	    double *p=recon->fitNW->p[idm]->p+(inw+idm)*nloc;
	    for(int iloc=0; iloc<nloc; iloc++){
		p[iloc]=scl;
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
		for(int iloc=0; iloc<nloc; iloc++){
		    p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;//x tilt
		    p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;//y tilt
		}
	    }
	}else if(lrt_tt==2){//Canceling TT. only valid for 2 DMs
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
		for(int iloc=0; iloc<nloc; iloc++){
		    p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;//x tilt
		    p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;//y tilt
		}
	    }

	}
	inw+=2*(ndm-1);
    }
  
    if(parms->save.setup){
	dcellwrite(recon->fitNW,"%s/fitNW",dirsetup);
    }
}

/**
   Assemble the DM fitting matrix

   The fitting is done by minimizing \f$||H_X x - H_A a||^2_W\f$ where \f$H_X,
   H_A\f$ are ray tracing operator from tomography grid xloc, and deformable
   mirror grid aloc to pupil grid ploc. The norm is weighted using bilinear
   influence functions within the telescope aperture. We have
   
   \f$a=\left[H_A^T(W_0-W_1 W_1^T)H_A\right]^{-1} H_A^T (W_0-W_1) H_X x\f$

   For details see http://www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803
*/
static void
setup_recon_fit_matrix(RECON_T *recon, const PARMS_T *parms){
    const int nfit=parms->fit.nfit;
    const int ndm=parms->ndm;
    if(ndm==0) return;
    spcell *HATc=spcelltrans(recon->HA);
    PDSPCELL(HATc, HAT);
    PDSPCELL(recon->HA,HA);

    info2("Before assembling fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    //Assemble Fit matrix.
    int npsr=recon->npsr;
    if(parms->load.fit){
	if(!(zfexist("FRM") && zfexist("FRU") && zfexist("FRV"))){
	    error("FRM, FRU, FRV (.bin or .bin.gz) not all exist\n");
	}
	warning("Loading saved recon->FR\n");
	recon->FR.M=spcellread("FRM");
	recon->FR.U=dcellread("FRU");
	recon->FR.V=dcellread("FRV");
    }else{
	if(recon->HXF){
	    info2("Building recon->FR\n");
	    recon->FR.M=spcellnew(ndm, npsr);
	    PDSPCELL(recon->FR.M, FRM);
	    PDSPCELL(recon->HXF, HXF);

	    for(int ips=0; ips<npsr; ips++){
		for(int ifit=0; ifit<nfit; ifit++){
		    if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		    dsp *tmp=spmulsp(recon->W0, HXF[ips][ifit]);
		    for(int idm=0; idm<ndm; idm++){
			spmulsp2(&FRM[ips][idm],HAT[ifit][idm], tmp, 
				 recon->fitwt->p[ifit]);
		    }
		    spfree(tmp);
		}
	    }
	    recon->FR.V=dcellnew(npsr, 1);
	    dmat **FRV=recon->FR.V->p;  
	
	    for(int ips=0; ips<npsr; ips++){
		int nloc=recon->xloc[ips]->nloc;
		FRV[ips]=dnew(nloc,nfit);
		for(int ifit=0; ifit<nfit; ifit++){
		    //notice the sart.
		    if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		    sptmulvec(FRV[ips]->p+ifit*nloc, 
			      HXF[ips][ifit], recon->W1->p, 
			      sqrt(recon->fitwt->p[ifit]));
		}
	    }
	    if(parms->save.recon){
		spcellwrite(recon->FR.M,"%s/FRM",dirsetup);
		dcellwrite(recon->FR.V,"%s/FRV",dirsetup);
	    }
	}else{
	    info("Avoid building recon->FR.M\n");
	    recon->FR.M=NULL;
	    recon->FR.V=NULL;
	}
	//Always need FR.U as it is used to do FL.U, FL.V
	recon->FR.U=dcellnew(ndm, 1);
	dmat **FRU=recon->FR.U->p;
	
	for(int idm=0; idm<ndm; idm++){    
	    int nloc=recon->aloc[idm]->nloc;
	    FRU[idm]=dnew(nloc, nfit);
	    for(int ifit=0; ifit<nfit; ifit++){
		//notice the sart.
		if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		sptmulvec(FRU[idm]->p+ifit*nloc, 
			  HA[idm][ifit], recon->W1->p,
			  sqrt(recon->fitwt->p[ifit]));
	    }
	}
	if(parms->save.recon){
	    dcellwrite(recon->FR.U,"%s/FRU",dirsetup);
	}
    }

    //Compute the proper strength of the low rand terms in FLM
    recon->fitscl=0;
    if(recon->FR.U && recon->FR.U->p){
	dmat **FRU=recon->FR.U->p;
	for(int idm=0; idm<ndm; idm++){
	    //follow loas convention.
	    double scl0=maxabs(FRU[idm]->p,FRU[idm]->nx*FRU[idm]->ny);
	    if(scl0>recon->fitscl) recon->fitscl=scl0;
	}
    }else{
	recon->fitscl=1./recon->ploc->nloc;//scale the constraints. Important!!
    }

    if(parms->load.fit){
	if(!(zfexist("FLM") && zfexist("FLU") && zfexist("FLV"))){
	    error("FLM, FLU, FLV (.bin or .bin.gz) not all exist\n");
	}
	warning("Loading saved recon->FL\n");
	recon->FL.M=spcellread("FLM");
	recon->FL.U=dcellread("FLU");
	recon->FL.V=dcellread("FLV");
    }else{
	fit_prep_lrt(recon,parms);
	if(parms->fit.actslave && parms->fit.alg!=1){
	    /*
	      2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found that slaving is causing mis-measurement of a few edge actuators. First try to remove W1. Or lower the weight. Revert back.
	     */
	    recon->actslave=slaving(recon->aloc, recon->HA, recon->W1, recon->fitNW, recon->actstuck, recon->actfloat, 0.1, 1./recon->ploc->nloc);
	    if(parms->save.setup){
		spcellwrite(recon->actslave,"%s/actslave",dirsetup);
		dcellwrite(recon->fitNW,"%s/fitNW2",dirsetup);
	    }
	}
	info2("Building recon->FL\n");
	recon->FL.M=spcellnew(ndm, ndm);
	dsp *(*FLM)[ndm]=(dsp*(*)[ndm])recon->FL.M->p;
	for(int idm=0; idm<ndm; idm++){
	    for(int ifit=0; ifit<nfit; ifit++){
		if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		dsp *tmp=spmulsp(recon->W0, HA[idm][ifit]);
		for(int jdm=0; jdm<ndm; jdm++){
		    spmulsp2(&FLM[idm][jdm],HAT[ifit][jdm], tmp,
			     recon->fitwt->p[ifit]);
		}
		spfree(tmp);
	    }
	}
	spcellfree(HATc);
	if(fabs(parms->fit.tikcr)>1.e-15){
	    double tikcr=parms->fit.tikcr;
	    /*Estimated from the formula.  1/nloc is due to W0, the other
	      scaling is due to ray tracing between different sampling freq.*/
	    int nact=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		nact+=recon->aloc[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    info2("Adding tikhonov constraint of %g to FLM\n", tikcr);
	    info2("The maximum eigen value is estimated to be around %e\n", maxeig);
	    spcelladdI(recon->FL.M,tikcr*maxeig);
	}

	{//Low rank terms.
	    recon->FL.U=dcellcat_each(recon->FR.U, recon->fitNW, 2);
	    dcell *tmp=NULL;//negative NW.
	    dcelladd(&tmp, 1, recon->fitNW, -1);
	    recon->FL.V=dcellcat_each(recon->FR.U, tmp, 2);
	    dcellfree(tmp);
	}
	if(recon->actslave){
	    spcelladd(&recon->FL.M, recon->actslave);
	}
	//spcellsym(recon->FL.M);
	info2("DM Fit number of Low rank terms: %ld in RHS, %ld in LHS\n",
	      recon->FR.U->p[0]->ny, recon->FL.U->p[0]->ny);
	if(parms->save.recon){
	    spcellwrite(recon->FL.M,"%s/FLM.bin",dirsetup);
	    dcellwrite(recon->FL.U,"%s/FLU",dirsetup);
	    dcellwrite(recon->FL.V,"%s/FLV",dirsetup);
	}
    }
    if((parms->fit.alg==0 || parms->fit.alg==2) && parms->fit.bgs){
	muv_direct_diag_prep(&(recon->FL),(parms->fit.alg==2)*parms->fit.svdthres);
    }
    if(((parms->fit.alg==0 || parms->fit.alg==2) && !parms->fit.bgs)
       ||((parms->tomo.split==2 && !parms->load.mvst) && !recon->FL.C && !recon->RL.MI)){
	if(fabs(parms->fit.tikcr)<1.e-14){
	    warning("tickcr=%g is too small, chol may fail.\n", parms->fit.tikcr);
	}
	muv_direct_prep(&(recon->FL),(parms->fit.alg==2)*parms->fit.svdthres);
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    if(parms->save.recon){
       	if(recon->FL.C)
	    chol_save(recon->FL.C,"%s/FLC.bin",dirsetup);
	if(recon->FL.MI)
	    dwrite(recon->FL.MI,"%s/FLMI", dirsetup);
	if(recon->FL.Up)
	    dwrite(recon->FL.Up, "%s/FLUp", dirsetup);
	if(recon->FL.Vp)
	    dwrite(recon->FL.Vp, "%s/FLVp", dirsetup);
	if(recon->FL.CB){
	    for(int ib=0; ib<recon->FL.nb; ib++){
		chol_save(recon->FL.CB[ib],"%s/FLCB_%d.bin",dirsetup, ib);
	    }
	}
	if(recon->FL.MIB){
	    dcellwrite(recon->FL.MIB,"%s/FLMIB", dirsetup);
	}
    }
    info2("After assemble fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}
/**
   Create the reconstructor to reconstruct the residual focus error due to LGS
   sodium tracking error. Need to average out the focus error caused by
   atmosphere when applying (a low pass filter is applied to the output).  */
static void
setup_recon_focus(RECON_T *recon, POWFS_T *powfs,
		  const PARMS_T *parms){
    int ilgs=-1;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].llt){
	    if(ilgs==-1){
		ilgs=ipowfs;
	    }else{
		warning("There are multiple LGS type\n");
	    }
	}
    }
    if(ilgs==-1){
	warning("There are no LGS with llt. \n");
	return;
    }
    //Create Gfocus: Focus mode -> WFS grad
    dcell *Gfocus=dcellnew(parms->nwfsr, 1);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	dmat *opd=dnew(powfs[ipowfs].loc->nloc,1);
	loc_add_focus(opd->p, powfs[ipowfs].loc, 1);
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	if(parms->powfs[ipowfs].gtype_recon==1){
	    dcell *saimcc;
	    if(powfs[ipowfs].locm){//We can not use realamp as in simulation.
		dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].amp->p);
		saimcc=dcellinvspd_each(mcc);
		dcellfree(mcc);
	    }else{
		saimcc=dcellref(powfs[ipowfs].saimcc[0]);
	    }
	    pts_ztilt(&Gfocus->p[iwfs], powfs[ipowfs].pts, saimcc,
		      powfs[ipowfs].realamp[wfsind], opd->p);
	}else{
	    spmulmat(&Gfocus->p[iwfs], adpind(powfs[ipowfs].GS0,wfsind), opd, 1);
	}
	dfree(opd);
    }
    dmat *GMGngs=NULL;
    dcell *GMngs=dcellnew(1, parms->nwfsr);
    /*Compute focus constructor from NGS Grads. fuse grads
      together to construct a single focus measurement*/
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].trs==0 && parms->powfs[ipowfs].order>1){
	    info2("wfs %d will be used to track focus\n", iwfs);
	}else{
	    continue;
	}
	spmulmat(&GMngs->p[iwfs], recon->saneai->p[iwfs+parms->nwfsr*iwfs], 
		 Gfocus->p[iwfs],1);
	dmm(&GMGngs,Gfocus->p[iwfs], GMngs->p[iwfs], "tn",1);
    }
    dinvspd_inplace(GMGngs);
    dcell *RFngs=recon->RFngs=dcellnew(1, parms->nwfsr);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	if(!Gfocus->p[iwfs]) continue;
	dmm(&RFngs->p[iwfs],GMGngs, GMngs->p[iwfs],"nt",1);
    }
    dfree(GMGngs);
    dcellfree(GMngs);
    /*
      Compute focus constructor from LGS grads. A
      constructor for each LGS.
    */
    dcell *RFlgs=recon->RFlgs=dcellnew(parms->nwfsr, 1);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].llt)
	    continue;
	dmat *GMtmp=NULL;
	dmat *GMGtmp=NULL;
	spmulmat(&GMtmp, recon->saneai->p[iwfs+parms->nwfsr*iwfs], 
		 Gfocus->p[iwfs], 1);
	dmm(&GMGtmp, Gfocus->p[iwfs], GMtmp, "tn",1);
	dinvspd_inplace(GMGtmp);
	dmm(&RFlgs->p[iwfs], GMGtmp, GMtmp, "nt", 1);
	dfree(GMtmp);
	dfree(GMGtmp);
    }
  
    if(parms->save.setup){
	dcellwrite(Gfocus,"%s/Gfocus",dirsetup);
	dcellwrite(RFngs,"%s/RFngs",dirsetup);
	dcellwrite(RFlgs,"%s/RFlgs",dirsetup);
    }
    dcellfree(Gfocus);
}


/**
   compute the MVST split tomography NGS mode reconstructor.

   2010-03-16:
   New Implementation.

   Definition:

   - \f$G_{lgs}\f$ is LGS gradient operator from xloc, contained in recon->GXtomo as spcell
   - \f$G_{ngs}\f$ is NGS gradient operator from xloc, contained in recon->GXL as dense matrix
   - \f$C_{lgs}, C_{ngs}\f$ is the LGS, and NGS measurement noise covariance matrix.
   - \f$\hat{x}_{lgs}\f$ is LGS tomography output
   - \f$\hat{x}_{ngs}\f$ is NGS minimum variance split tomography output
   - \f$a_{ngs}\f$ is NGS minimum variance DM output.
   - \f$F\f$ is the fitting operator
   - \f$H_A\f$ is the ray tracing operator from aloc to ploc, contained in recon->HA.
   - \f$MVModes\f$ is the MVST NGS modes
   - \f$MVRngs\f$ is the MVST NGS reconstructor

   We have
 
   \f{eqnarray*}{
   \hat{x}_{lgs}&=&A^{-1}G^{T}_{lgs}C_{lgs}^{-1}s_{lgs}\\
   Uw&=&A^{-1}G_{ngs}^TC_{ngs}\\
   \hat{x}_{ngs}&=&Uw(1+G_{ngs}Uw)^{-1}(s_{ngs}-G_{ngs}\hat{x}_{lgs});\\
   a_{NGS}&=&F\hat{x}_{ngs}\\
   MVModes&=&F\cdot Uw\\
   MVRngs&=&(1+G_{ngs}Uw)^{-1}
   \f}

   We need to orthnormalize the NGS modes. Propagate the NGS modes \f$F\cdot Uw\f$ to
   fitting directions and compute the cross-coupling matrix with weighting on \f$W_0, W_1\f$, we have

   \f{eqnarray*}{
   Q&=&H_A\cdot F \cdot Uw\\
   M_{CC}&=&Q^T(W_0-W_1 W_1^T)Q\\
   \f}

   Do SVD on \f$MCC\f$ we have \f$M_{CC}=u\Sigma v^T \f$ where \f$u{\equiv}v\f$
   because \f$MCC\f$ is symmetric. Redefine the NGS modes and reconstructor as

   \f{eqnarray*}{
   MVModes&=&F\cdot Uw\cdot u \Sigma^{-1/2}\\
   MVRngs&=&\Sigma^{1/2}u^T (1+G_{ngs}A^{-1}G_{ngs}^T C_{ngs}^{-1})^{-1}
   \f}
*/

void
setup_recon_mvst(RECON_T *recon, const PARMS_T *parms){
    /*
      Notice that: Solve Fitting on Uw and using FUw to form Rngs gives
      slightly different answer than solve fitting after assemble the
      reconstructor. 10^-6 relative difference.
      
      2010-03-10: Bug found and fixed: The MVST with CBS-CBS method gives worst
      performance than integrated tomography. The probelm is in PUm
      computing. I mistakenly called chol_solve, while I should have
      called muv_direct_solve. The former doesn ot apply the low rank
      terms.
    */
    if(parms->tomo.split!=2){
	return;
    }
 
    dcellfree(recon->GXL);
    spcellfull(&recon->GXL, recon->GXlo, 1);

    dcell *U=NULL; 
    dcell *FU=NULL;
    if(parms->load.mvst){
	U=dcellread("mvst_U");
	FU=dcellread("mvst_FU");
    }else{
	//Prepare CBS if not already done
	
	if(!recon->RL.C && !recon->RL.MI){
	    muv_direct_prep(&(recon->RL), 0);
	}
	if(!recon->FL.C && !recon->FL.MI){
	    muv_direct_prep(&(recon->FL), 0);
	}
	
	dcell *GXLT=dcelltrans(recon->GXL);
	muv_direct_solve_cell(&U, &recon->RL, GXLT);
	dcellfree(GXLT);
	dcell *rhs=NULL;
	muv(&rhs, &recon->FR, U, 1);
	muv_direct_solve_cell(&FU, &recon->FL, rhs);
	dcellfree(rhs);
    }
    if(parms->save.mvst || parms->save.setup){
	dcellwrite(U, "%s/mvst_U", dirsetup);
	dcellwrite(FU, "%s/mvst_FU", dirsetup);
    }
    
    dcell *Uw=NULL;
    dcell *FUw=NULL;
    
    dcellmulsp(&Uw, U, recon->saneai, 1);
    dcellmulsp(&FUw, FU, recon->saneai, 1);
    dcell *M=NULL;
    dcellmm(&M, recon->GXL, Uw, "nn", 1);
    dcelladdI(M, 1);
    dcell *Minv=dcellinv(M);
    dcellfree(M);
    if(0){//Orthnormalize the Modes.
	/*
	  Change FUw*Minv -> FUw*(U*sigma^-1/2) * (U*sigma^1/2)'*Minv
	  columes of FUw*(U*sigma^-1/2) are the eigen vectors.

	  U, sigma is the eigen value decomposition of <HA FUw W FUw' HA'>
	*/
	dcell *Q=NULL;//the NGS modes in ploc.
	spcellmulmat(&Q, recon->HA, FUw, 1);
	dcell *QwQc=calcWmcc(Q,Q,recon->W0,recon->W1,recon->fitwt);
	dmat *QwQ=dcell2m(QwQc);
	dmat *QSdiag=NULL, *QU=NULL, *QVt=NULL;
	dsvd(&QU, &QSdiag, &QVt, QwQ);
	if(parms->save.setup) dwrite(QSdiag,"%s/mvst_QSdiag",dirsetup);
	dcwpow(QSdiag, -1./2.);
	dmuldiag(QU,QSdiag);//U*sigma^-1/2
	d2cell(&QwQc,QU,NULL);
	dcell *FUw_keep=FUw;FUw=NULL; 
	dcellmm(&FUw, FUw_keep, QwQc, "nn", 1);
	dcellfree(FUw_keep);
	dcwpow(QSdiag,-2);
	dmuldiag(QU,QSdiag);//U*sigma^1/2 (From U*sigma^(-1/2)*sigma)
	d2cell(&QwQc,QU,NULL);
	dcell *Minv_keep=Minv; Minv=NULL;
	dcellmm(&Minv,QwQc,Minv_keep,"tn", 1);
	dcellfree(Minv_keep);
	dcellfree(Q);
	dcellfree(QwQc);
	dfree(QwQ);
	dfree(QSdiag);
	dfree(QVt);
	dfree(QU);
    }
    recon->MVRngs=dcellreduce(Minv,1);//1xnwfs cell
    recon->MVModes=dcellreduce(FUw,2);//ndmx1 cell
    recon->MVGM=NULL;
    spcellmulmat(&recon->MVGM, recon->GAlo, recon->MVModes, 1);

    dcellfree(Minv);
    dcellfree(U);
    dcellfree(FU);
    dcellfree(FUw);
    dcellfree(Uw);
    if(parms->save.setup){
	dcell *Qn=NULL;
	spcellmulmat(&Qn, recon->HA, recon->MVModes, 1);
	dcell *Qntt=dcellnew(Qn->nx,Qn->ny);
	dmat *TTploc=loc2mat(recon->ploc,1);//TT mode. need piston mode too!
	dmat *PTTploc=dpinv(TTploc,NULL,recon->W0);//TT projector. no need w1 since we have piston.
	dfree(TTploc);
	for(int ix=0; ix<Qn->nx*Qn->ny; ix++){
	    if(!Qn->p[ix]) continue;
	    dmm(&Qntt->p[ix],PTTploc, Qn->p[ix],"nn",1);
	}
	dcellwrite(Qntt,"%s/mvst_modptt", dirsetup);
	dcellfree(Qn);
	dcellfree(Qntt);
	dfree(PTTploc);
	dcellwrite(recon->MVRngs, "%s/mvst_Rngs",dirsetup);
	dcellwrite(recon->MVModes,"%s/mvst_Modes",dirsetup);
    }
    if(parms->dbg.mvstlimit>0){//limit number of modes used.
	warning("MVST: Correction is limited to %d modes\n", parms->dbg.mvstlimit);
	dmat *tmp;
	for(int iy=0; iy<recon->MVRngs->ny; iy++){
	    tmp=recon->MVRngs->p[iy];
	    if(tmp){
		recon->MVRngs->p[iy]=dsub(tmp,0,parms->dbg.mvstlimit,0,tmp->ny);
		dfree(tmp);
	    }
	}
	for(int ix=0; ix<recon->MVModes->nx; ix++){
	    tmp=recon->MVModes->p[ix];
	    if(tmp){
		recon->MVModes->p[ix]=dsub(tmp,0,tmp->nx,0,parms->dbg.mvstlimit);
		dfree(tmp);
	    }
	}
	if(parms->save.setup){
	    dcellwrite(recon->MVRngs, "%s/mvst_Rngs_limit",dirsetup);
	    dcellwrite(recon->MVModes,"%s/mvst_Modes_limit",dirsetup);
	}
    }
    /*
    if(parms->save.setup){
	dcell *QQ=NULL;
	spcellmulmat(&QQ, recon->HA, recon->MVModes,1);
	dcell *MCC=calcWmcc(QQ,QQ,recon->W0,recon->W1,recon->fitwt);
	dcellwrite(MCC,"%s/mvst_MCC",dirsetup);
    
	dcellfree(MCC);
	dcellfree(QQ);
	}*/
}

/**
   Sets up the tomogrpahy turbulence reconstruction structs including wavefront
   reconstructor and DM fitting operator \callgraph
   
   Calls setup_recon_tomo_matrix() 
   and setup_recon_fit_matrix() to setup the tomography and DM
   fitting matrix. 
 
   AHST is handled in setup_ngsmod().

   MVST is handled in setup_recon_mvst().
   
   MOAO is handled in setup_recon_moao().
   
*/
    
void setup_recon_mvr(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    TIC;tic;
    if(parms->cn2.tomo){
	/*Use cn2 estimation results for tomography. Use its ht to build
	  reconstructor matrices.*/
	CN2EST_T *cn2est=recon->cn2est;
	recon->ht=dref(cn2est->htrecon);
	recon->os=dref(cn2est->os);
	recon->wt=dref(cn2est->wtrecon->p[0]);
	//the following will be updated later in simulation.
	dset(recon->wt, 1./recon->wt->nx);//evenly distributed. 
	recon->r0=0.15;//random guess
	recon->l0=30;//random guess
    }else{//use input information from atmr
	recon->wt=dnew(parms->atmr.nps,1);
	recon->ht=dnew(parms->atmr.nps,1);
	recon->os=dnew(parms->atmr.nps,1);
	for(int ips=0; ips<recon->wt->nx; ips++){
	    recon->wt->p[ips]=parms->atmr.wt[ips];
	    recon->ht->p[ips]=parms->atmr.ht[ips];
	    recon->os->p[ips]=parms->atmr.os[ips];
	}
	recon->r0=parms->atmr.r0;
	recon->l0=parms->atmr.l0;
    }
    //sampling of xloc
    recon->dx=dnew(recon->ht->nx, 1);
    for(int iht=0; iht<recon->ht->nx; iht++){
	double scale = 1.0 - recon->ht->p[iht]/parms->atmr.hs;
	recon->dx->p[iht]=(parms->atmr.dx/recon->os->p[iht])*scale;	
    }
    //number of reconstruction layers
    recon->npsr= recon->ht->nx;
   
    //copy over fitwt since we need a dmat
    int nfit=parms->fit.nfit;
    recon->fitwt=dnew(nfit,1);
    memcpy(recon->fitwt->p,parms->fit.wt,sizeof(double)*nfit);
  
    //setup atm reconstruction layer grid
    setup_recon_xloc(recon,parms);
    //setup xloc/aloc to WFS grad
    toc2("Generating xloc");
    if(!parms->sim.idealfit &&
       !(parms->sim.recon==1 && (parms->sim.wfsalias || parms->sim.idealwfs))){
	setup_recon_HXW(recon,parms);
	setup_recon_GX(recon,parms);
	spcellfree(recon->HXW);//only keep HXWtomo for tomography
	//setup inverse noise covariance matrix.
	toc2("Generating GX GA");
	//prepare for tomography setup
	setup_recon_tomo_prep(recon,parms);
	toc2("Prepare tomography");
	if(parms->tomo.assemble || parms->tomo.split==2 || parms->sim.psfr){
	    /*assemble the matrix only if not using CG CG apply the
	      individual matrices on fly to speed up and save memory. */
	    setup_recon_tomo_matrix(recon,parms,aper);
	    toc2("Generating tomography matrix");
	}
	//Fall back function method if .M is NULL
	recon->RL.Mfun=TomoL;
	recon->RL.Mdata=recon;
	recon->RR.Mfun=TomoR;
	recon->RR.Mdata=recon;
	if(parms->tomo.alg==1){//CG
	    switch(parms->tomo.precond){
	    case 0://no preconditioner
		recon->RL.pfun=NULL;
		recon->RL.pdata=NULL;
		break;
	    case 1:
		recon->RL.pfun=fdpcg_precond;
		recon->RL.pdata=(void*)recon;
		break;
	    default:
		error("Invalid tomo.precond");
	    }
	}
	recon->RL.alg = parms->tomo.alg;
	recon->RL.bgs = parms->tomo.bgs;
	recon->RL.warm  = recon->warm_restart;
	recon->RL.maxit = parms->tomo.maxit;
    }
    if(!parms->sim.idealfit){//In idealfit, xloc has high sampling. We avoid HXF.
	setup_recon_HXF(recon,parms);
    }
    setup_recon_HA(recon,parms);
    //always assemble fit matrix, faster if many directions
    setup_recon_fit_matrix(recon,parms);
    toc2("Generating fit matrix");
    //Fall back function method if FR.M is NULL (!HXF<-idealfit)
    recon->FR.Mfun  = FitR;
    recon->FR.Mdata = recon;
    //Fall back function method if FL.M is NULL
    recon->FL.Mfun  = FitL;
    recon->FL.Mdata = recon;
    recon->FL.alg = parms->fit.alg;
    recon->FL.bgs = parms->fit.bgs;
    recon->FL.warm  = recon->warm_restart;
    recon->FL.maxit = parms->fit.maxit;
    //moao
    setup_recon_moao(recon,parms);
    toc2("Preparing moao");
    if(parms->sim.mffocus){
	setup_recon_focus(recon, powfs, parms);
    }
    if(parms->tomo.split){
	//split tomography
	if(parms->tomo.split){
	    if(parms->ndm<=2){
		//setup the ngsmode in both ahst and mvst mode 
		setup_ngsmod(parms,recon,aper,powfs);
	    }else if(parms->tomo.split==1){
		error("Not implemented");
	    }
	}
	if(!parms->sim.idealfit && parms->tomo.split==2){//Need to be after fit
	    setup_recon_mvst(recon,parms);
	}
    }
    if(!parms->sim.idealfit && parms->tomo.precond==1){
	recon->fdpcg=fdpcg_prepare(parms, recon, powfs);
	toc2("Preparing fdpcg");
    }
    //The following have been used in fit matrix.
    dcellfree(recon->fitNW);
    spcellfree(recon->actslave);
    if(recon->FR.M && !parms->sim.wfsalias && !parms->sim.idealwfs){
	spfree(recon->W0); 
	dfree(recon->W1); 
	spcellfree(recon->HA); 
	spcellfree(recon->HXF); 
    }
 
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(!parms->powfs[ipowfs].hasGS0 && powfs[ipowfs].GS0){
	    spcellfree(powfs[ipowfs].GS0);
	    powfs[ipowfs].GS0=NULL;
	}
    }

    /*
      The following arrys are not used after preparation is done.
    */
    spcellfree(recon->GX);
    spcellfree(recon->GXhi);
    spcellfree(recon->GXtomo);//we use HXWtomo instead. faster
    if(!(parms->cn2.tomo && parms->tomo.split==2)){//mvst needs GXlo when updating.
	spcellfree(recon->GXlo);
    }
    if(parms->tomo.alg!=1 || parms->tomo.assemble){
	spcellfree(recon->HXWtomo);
    }
    toc2("setup_recon");
}
/**
   Setup the least square reconstruct by directly inverting GA matrix. 
   The reconstructor is simply the pseudo inverse of GA matrix:
   \f[\hat{x}=(G_a^TC_g^{-1}G_a)^{-1}G_a^TC_g^{-1}\f]

   This is very close to RR except replacing GX with GA.

   We use the tomogrpahy parameters for lsr, since lsr is simply "tomography" onto DM directly.
*/
void setup_recon_lsr(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    spcell *GAlsr;
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;

    if(parms->tomo.split){
	//high order wfs only in split mode.
	GAlsr=recon->GAhi;
	setup_ngsmod(parms,recon,aper,powfs);
    }else{
	//all wfs in integrated mode.
	GAlsr=recon->GA;
    }
    spcell *GAlsrT=spcelltrans(GAlsr);
    info2("Building recon->LR\n");
    recon->LR.M=spcellmulspcell(GAlsrT, recon->saneai, 1);
    spcellfree(GAlsrT);
    /*
      Tip/tilt and diff focus removal low rand terms for LGS WFS.
    */
    if(recon->TTF){
	spcellmulmat(&recon->LR.U, recon->LR.M, recon->TTF, 1);
	recon->LR.V=dcelltrans(recon->PTTF);
    }
    info2("Building recon->LL\n");
    recon->LL.M=spcellmulspcell(recon->LR.M, GAlsr, 1);
    double maxeig=pow(recon->neamhi * recon->aloc[0]->dx, -2);
    if(fabs(parms->tomo.tikcr)>EPS){
	info2("Adding tikhonov constraint of %g to LLM\n", parms->tomo.tikcr);
	info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	spcelladdI(recon->LL.M, parms->tomo.tikcr*maxeig);
    }
    dcell *NW=NULL;
    if(parms->tomo.alg!=2){//Not SVD, need low rank terms for piston/checkboard constraint.
	NW=dcellnew(ndm,1);
	int nmod=2;//two modes.
	for(int idm=0; idm<ndm; idm++){
	    loc_create_map(recon->aloc[idm]);
	    const long nloc=recon->aloc[idm]->nloc;
	    NW->p[idm]=dnew(nloc, ndm*nmod);
	    double *p=NW->p[idm]->p+nmod*idm*nloc;
	    for(long iloc=0; iloc<nloc; iloc++){
		p[iloc]=1;//piston mode
	    }
	    //notice offset of 1 because map start count at 1
	    p=NW->p[idm]->p+(1+nmod*idm)*nloc-1;
	    locmap_t *map=recon->aloc[idm]->map;
	    long (*pmap)[map->nx]=(long(*)[map->nx])map->p;
	    for(long iy=0; iy<map->ny; iy++){
		for(long ix=0; ix<map->nx; ix++){
		    if(pmap[iy][ix]){
			p[pmap[iy][ix]]=(double)2*((iy+ix)&1)-1;
		    }
		}
	    }
	}
	//scale it to match the magnitude of LL.M
	dcellscale(NW, sqrt(maxeig));
	if(parms->save.setup){
	    dcellwrite(NW, "%s/lsrNW",dirsetup);
	}
    }
    if(parms->fit.actslave && parms->tomo.alg!=1){
	//actuator slaving. important. change from 0.5 to 0.1 on 2011-07-14.
	spcell *actslave=slaving(recon->aloc, recon->GAhi, NULL, NW, recon->actstuck, recon->actfloat, 0.1, sqrt(maxeig));
	if(parms->save.setup){
	    if(NW){
		dcellwrite(NW, "%s/lsrNW2",dirsetup);
	    }
	    spcellwrite(actslave,"%s/actslave", dirsetup);
	}
	spcelladd(&recon->LL.M, actslave);
	spcellfree(actslave);
    }
    //Low rank terms for low order wfs. Only in Integrated tomography.
    dcell *ULo=dcellnew(ndm,nwfs);
    PDCELL(ULo, pULo);
    dcell *VLo=dcellnew(ndm,nwfs);
    PDCELL(VLo, pVLo);
    PDSPCELL(recon->LR.M, LRM);
    PDSPCELL(recon->GA, GA);

    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip || !parms->powfs[ipowfs].lo){
	    continue;
	}
	for(int idm=0; idm<ndm; idm++){
	    spfull(&pULo[iwfs][idm], LRM[iwfs][idm],-1);
	    sptfull(&pVLo[iwfs][idm], GA[idm][iwfs],1);
	}
    }
    recon->LL.U=dcellcat(recon->LR.U, ULo, 2);
    dcell *GPTTDF=NULL;
    sptcellmulmat(&GPTTDF, recon->GA, recon->LR.V, 1);
    recon->LL.V=dcellcat(GPTTDF, VLo, 2);
    dcellfree(GPTTDF);
    dcellfree(ULo);
    dcellfree(VLo);
    if(NW){
	info2("Create piston and check board modes that are in NULL space of GA.\n");
	//add to low rank terms.
	dcell *tmp=recon->LL.U;
	recon->LL.U=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellscale(NW, -1);
	tmp=recon->LL.V;
	recon->LL.V=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellfree(NW);
    }
    recon->LL.alg = parms->tomo.alg;
    recon->LL.bgs = parms->tomo.bgs;
    recon->LL.warm = recon->warm_restart;
    recon->LL.maxit = parms->tomo.maxit;
    //Remove empty cells.
    dcelldropempty(&recon->LR.U,2);
    dcelldropempty(&recon->LR.V,2);
    dcelldropempty(&recon->LL.U,2);
    dcelldropempty(&recon->LL.V,2);
    if(parms->save.recon){
	spcellwrite(recon->LR.M,"%s/LRM",dirsetup);
	dcellwrite(recon->LR.U,"%s/LRU",dirsetup);
	dcellwrite(recon->LR.V,"%s/LRV",dirsetup);
	spcellwrite(recon->LL.M,"%s/LLM.bin",dirsetup);//disable compression
	dcellwrite(recon->LL.U,"%s/LLU",dirsetup);
	dcellwrite(recon->LL.V,"%s/LLV",dirsetup); 
    }
    if(parms->tomo.alg==0 || parms->tomo.alg==2){
	if(!parms->tomo.bgs){
	    muv_direct_prep(&recon->LL, (parms->tomo.alg==2)*parms->tomo.svdthres);
	    if(parms->save.recon){
		if(recon->LL.C)
		    chol_save(recon->LL.C, "%s/LLC.bin", dirsetup);
		else
		    dwrite(recon->LL.MI, "%s/LLMI.bin", dirsetup);
	    }
	    spcellfree(recon->LL.M);
	    dcellfree(recon->LL.U);
	    dcellfree(recon->LL.V);	
	}else{
	    muv_direct_diag_prep(&(recon->LL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	    if(parms->save.recon){
		for(int ib=0; ib<recon->LL.nb; ib++){
		    if(recon->LL.CB)
			chol_save(recon->LL.CB[ib],"%s/LLCB_%d.bin",dirsetup, ib);
		    else
			dwrite(recon->LL.MI,"%s/LLMIB_%d.bin", dirsetup, ib);
		}
	    }
	    //Don't free M, U, V
	}
    }
}
/**
   Setup either the minimum variance reconstructor by calling setup_recon_mvr()
   or least square reconstructor by calling setup_recon_lsr() */
RECON_T *setup_recon(const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    RECON_T * recon = calloc(1, sizeof(RECON_T));
    recon->parms=parms;//save a pointer.
    if(parms->cn2.npair){
	/*setup CN2 Estimator. It determines the reconstructed layer heigh can
	  be fed to the tomography */
	recon->cn2est=cn2est_prepare(parms,powfs);
    }
    recon->warm_restart = parms->atm.frozenflow && !parms->dbg.ntomo_maxit;
    //number of deformable mirrors
    recon->ndm = parms->ndm;
    if(recon->warm_restart){
	info2("Using warm restart\n");
    }else{
	warning2("Do not use warm restart\n");
    }
    //to be used in tomography.
    recon->nthread=parms->sim.nthread;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].trs){
	    recon->has_ttr=1;
	    break;
	}
    }
    if(parms->sim.recon!=2){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs<=1) continue;
	    if(parms->powfs[ipowfs].dfrs){
		recon->has_dfr=1;
		break;
	    }
	}   
    }
    //setup DM actuator grid
    setup_recon_aloc(recon,parms);
    //setup pupil coarse grid
    setup_recon_ploc(recon,parms);
    setup_recon_GWR(recon, parms, powfs);
    setup_recon_GP(recon,parms,powfs,aper);
    setup_recon_GA(recon,parms,powfs);
    //assemble noise equiva angle inverse from powfs information
    setup_recon_saneai(recon,parms,powfs);
    //setup LGS tip/tilt/diff focus removal
    setup_recon_TTFR(recon,parms,powfs);
    switch(parms->sim.recon){
    case 0:
	setup_recon_mvr(recon, parms, powfs, aper);
	break;
    case 1:
    case 2:
	setup_recon_lsr(recon, parms, powfs, aper);
	/*if(parms->sim.wfsalias || parms->sim.idealwfs){
	    setup_recon_mvr(recon, parms, powfs, aper);
	    }*/
	break;
    default:
	error("sim.recon=%d is not recognized\n", parms->sim.recon);
    }
    spcellfree(recon->GWR);
    if(parms->sim.recon!=0 || (parms->tomo.assemble && !parms->cn2.tomo)){
	//We already assembled tomo matrix. don't need these matric any more.
	dcellfree(recon->TTF);
	dcellfree(recon->PTTF);
    }
    /* Free arrays that will no longer be used after reconstruction setup is done. */
    spcellfree(recon->sanea); 
    spcellfree(recon->saneal);
    dfree(recon->neam); 
    if(!(parms->tomo.assemble && parms->tomo.alg==1) && !parms->cn2.tomo && !parms->tomo.bgs){
	//We no longer need RL.M,U,V
	spcellfree(recon->RL.M);
	dcellfree(recon->RL.U);
	dcellfree(recon->RL.V);
    }
    if(parms->fit.alg!=1 && !parms->fit.bgs){
	spcellfree(recon->FL.M);
	dcellfree(recon->FL.U);
	dcellfree(recon->FL.V);
    }
    if(parms->tomo.alg==1){
	muv_direct_free(&recon->RL);
    }
    if(parms->fit.alg==1){
	muv_direct_free(&recon->FL);
    }
#if USE_CUDA
    if(use_cuda){
	gpu_setup_recon(parms, powfs, recon);
    }
#endif

    /*
      The following arrys are not used after preparation is done.
    */
    mapfree(aper->ampground);
    return recon;
}
/**
   Free the recon struct.
*/
void free_recon(const PARMS_T *parms, RECON_T *recon){
    free_recon_moao(recon, parms);
    dfree(recon->ht);
    dfree(recon->os);
    dfree(recon->wt);
    dfree(recon->dx);
    dcellfree(recon->MVRngs);
    dcellfree(recon->MVGM);
    dcellfree(recon->MVModes);
    dcellfree(recon->xmcc);
    spcellfree(recon->GX);
    spcellfree(recon->GP);
    spcellfree(recon->GXfocus);
    spcellfree(recon->GA); 
    spcellfree(recon->GAlo);
    spcellfree(recon->GAhi);
    dcellfree(recon->GXL);
    spcellfree(recon->L2);
    spcellfree(recon->L2save);
    free_cxx(recon);
    dcellfree(recon->TT);
    dcellfree(recon->PTT);
    dcellfree(recon->DF);
    dcellfree(recon->PDF);
    dcellfree(recon->TTF);
    dcellfree(recon->PTTF);
    dcellfree(recon->RFlgs);
    dcellfree(recon->RFngs);
    spcellfree(recon->ZZT);
    spcellfree(recon->HXF); 
    spcellfree(recon->HXW);
    spcellfree(recon->HXWtomo);
    spcellfree(recon->HA); 
    spfree(recon->W0); 
    dfree(recon->W1); 
    dcellfree(recon->fitNW);
    dfree(recon->fitwt);
    if(recon->ngsmod){
	dcellfree(recon->ngsmod->GM);
	dcellfree(recon->ngsmod->Rngs);
	dcellfree(recon->ngsmod->Pngs);
	dcellfree(recon->ngsmod->Modes);
	dfree(recon->ngsmod->MCC);
	dcellfree(recon->ngsmod->MCCP);
	dfree(recon->ngsmod->IMCC);
	dfree(recon->ngsmod->IMCC_TT);
	dcellfree(recon->ngsmod->Ptt);
	free(recon->ngsmod);
    }

    int npsr = recon->npsr;
    int ndm = parms->ndm;
   
    locarrfree(recon->xloc, npsr); recon->xloc=NULL;
    maparrfree(recon->xmap, npsr); recon->xmap=NULL;
    locfree(recon->ploc); recon->ploc=NULL;
    for(int idm=0; idm<ndm; idm++){
	if(recon->alocm[idm]!=recon->aloc[idm])
	    locfree(recon->alocm[idm]);
	locfree(recon->alocm[idm]);
	free(recon->aembed[idm]);
    }
    maparrfree(recon->amap, parms->ndm);
    mapfree(recon->pmap);
    free(recon->aembed);
    free(recon->alocm);
    free(recon->aloc);
    icellfree(recon->actstuck);
    icellfree(recon->actfloat);
    dcellfree(recon->aimcc);//used in filter.c
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    muv_free(&recon->FR);
    muv_free(&recon->FL);
    muv_free(&recon->LR);
    muv_free(&recon->LL);
    spcellfree(recon->saneai);

    fdpcg_free(recon->fdpcg);
    cn2est_free(recon->cn2est);
    free(recon);
}


/*
  some tips.
  1) UV balance need to be done carefully. The scaling must be a single number
  2) When add regularizations, the number has to be in the same order as the matrix.
  This is supper important!!!
  3) single point piston constraint in Tomography is not implemented correctly.

  2009-11-22
  Found a difference: Ha in LAOS computes the weighting even if not enough points.
  HA in AOS only computes the weighting only if 4 coupled points all exist.
  AOS used laos geometry that are not enough to cover all points. Will modify mkH/accphi.
  H modified. Performance now agree with laos with AOS W0/W1 and HA

*/
