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
#include "recon.h"
#include "fdpcg.h"
#include "ahst.h"
#include "cn2est.h"
#include "recon_utils.h"
#include "mvm_client.h"
#include "moao.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

/**
   \file setup_recon.c Contains routines that setup the wavefront reconstructor
   and DM fitting.  Use parms->wfsr instead of parms->wfs for wfs information,
   which hands GLAO mode correctly.x 

   TOMOSCALE is used for RLM, RRM, in MVST for M, and in CUDA due to
   limited dynamic range of single precision floating point numbers.*/

/**
   Setting up PLOC grid, which is a coarse sampled (usually halves the
   subaperture spacing) grid that defines the circular aperture for wavefront
   reconstruction and DM fitting.*/
static void
setup_recon_ploc(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    if(parms->load.ploc){/*optionally load ploc from the file. see dbg.conf */
	char *fn=parms->load.ploc;
	warning("Loading ploc from %s\n",fn);
	recon->ploc=locread("%s",fn);
	recon->pmap=loc2map(recon->ploc);
	free(recon->pmap->p);
	recon->pmap->p=NULL;
    }else{ 
	/*
	  Create a circular PLOC with telescope diameter by calling
	  create_metapupil with height of 0. We don't add any guard points. PLOC
	  does not need to be follow XLOC in FDPCG.*/
	double dxr=parms->atmr.dx/parms->tomo.pos;/*sampling of ploc */
	double guard=parms->tomo.guard*dxr;
	map_t *pmap=create_metapupil_wrap(parms,0,dxr,0,guard,0,0,0,parms->tomo.square);
	info2("PLOC is %ldx%ld, with sampling of %.2fm\n",pmap->nx,pmap->ny,dxr);
	recon->ploc=map2loc(pmap);/*convert map_t to loc_t */
	recon->pmap = pmap;
	free(pmap->p); pmap->p=NULL;
	if(parms->save.setup){
	    locwrite(recon->ploc, "%s/ploc",dirsetup);
	}
    }
    loc_create_stat(recon->ploc);
}
static void
setup_recon_floc(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    if(parms->load.floc){
	warning("Loading floc from %s\n", parms->load.floc);
	recon->floc=locread("%s", parms->load.floc);
	recon->fmap=loc2map(recon->floc);
    }else{
	double dxr=parms->atmr.dx/parms->fit.pos;/*sampling of floc */
	double guard=parms->tomo.guard*dxr;
	map_t *fmap=create_metapupil_wrap 
	    (parms,0,dxr,0,guard,0,0,0,parms->fit.square);
	info2("FLOC is %ldx%ld, with sampling of %.2fm\n",fmap->nx,fmap->ny,dxr);
	recon->fmap=fmap;
	recon->floc=map2loc(fmap);/*convert map_t to loc_t */
    }
    free(recon->fmap->p); recon->fmap->p=NULL;
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
	if(parms->save.setup){
	    spwrite(recon->W0, "%s/W0",dirsetup);
	    dwrite(recon->W1, "%s/W1",dirsetup);
	}
	if(parms->save.setup){
	    locwrite(recon->floc, "%s/floc",dirsetup);
	}
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
	/*int nxmax=0, nymax=0; */
	for(int idm=0; idm<ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dx=parms->dm[idm].dx;
	    double offset=parms->dm[idm].offset+(parms->dm[idm].order%2)*0.5;
	    const double guard=parms->dm[idm].guard*parms->dm[idm].dx;
	    
	    map_t *map=create_metapupil_wrap
		(parms,ht,dx,offset,guard,0,0,0,parms->fit.square);
	    info2("DM %d: grid is %ld x %ld\n", idm, map->nx, map->ny);
	    recon->aloc[idm]=map2loc(map);
	    recon->amap[idm]=map;
	    recon->aembed[idm]=map2embed(map);
	    free(map->p); map->p=NULL; 
	    free(recon->amap[idm]->nref);recon->amap[idm]->nref=NULL;
	}/*idm */
        if(parms->save.setup){
	    locarrwrite(recon->aloc,parms->ndm,"%s/aloc",dirsetup);
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
    recon->anx=calloc(ndm, sizeof(int));
    recon->any=calloc(ndm, sizeof(int));
    recon->anloc=calloc(ndm, sizeof(int));
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
    CALL_ONCE;
    const int npsr=recon->npsr;
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
	long nin0=0;
	/*FFT in FDPCG prefers power of 2 dimensions. for embeding and fast FFT*/
	if(parms->tomo.nxbase){
	    nin0=parms->tomo.nxbase;
	}else if(!parms->sim.idealfit && (parms->tomo.precond==1 || parms->tomo.square==2)){
	    /*same square grid dimension in meter on all layers.*/
	    long nx=0, ny=0;
	    for(int ips=0; ips<npsr; ips++){
		long nxi, nyi;
		double dxr=recon->dx->p[ips];
		create_metapupil(parms, recon->ht->p[ips], dxr, 0,
				 &nxi, &nyi, NULL, NULL, NULL, dxr*parms->tomo.guard,
				 0, 0, 0, 1);
		nxi/=recon->os->p[ips];
		nyi/=recon->os->p[ips];
		if(nx<nxi) nx=nxi;
		if(ny<nyi) ny=nyi;
	    }
	    nin0=nextfftsize(MAX(nx, ny));
	
	    while (parms->tomo.precond==1 && (nin0 & 1)){
		nin0=nextfftsize(nin0+1);
	    }
	}
	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    double dxr=(parms->sim.idealfit)?parms->atm.dx:recon->dx->p[ips];
	    const double guard=parms->tomo.guard*dxr;
	    long nin=nin0*recon->os->p[ips];
	    map_t *map=create_metapupil_wrap
		(parms,ht,dxr,0,guard,nin,nin,0,parms->tomo.square);
	    recon->xloc[ips]=map2loc(map);
	    recon->xmap[ips]=map;
	    info2("layer %d: xloc grid is %3ld x %3ld, sampling is %.3f m, %5ld points\n",
		  ips, map->nx,map->ny,dxr, recon->xloc[ips]->nloc);
	    /*Free the data and nref so we are assign pointers to it. */
	    free(recon->xmap[ips]->p);recon->xmap[ips]->p=NULL;
	    free(recon->xmap[ips]->nref);recon->xmap[ips]->nref=NULL;
	}
	if(parms->save.setup){
	    locarrwrite(recon->xloc, recon->npsr, "%s/xloc",dirsetup);
	}
    }
  
    recon->xnx=calloc(recon->npsr, sizeof(long));
    recon->xny=calloc(recon->npsr, sizeof(long));
    for(long i=0; i<recon->npsr; i++){
	recon->xnx[i]=recon->xmap[i]->nx;
	recon->xny[i]=recon->xmap[i]->ny;
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
	    
	    if(parms->recon.split!=2 && parms->powfs[ipowfs].skip && !parms->dbg.deltafocus){
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
	if(!parms->powfs[ipowfs].skip){/*for tomography */
	    for(int ips=0; ips<npsr; ips++){
		HXWtomo[ips][iwfs]=spref(HXW[ips][iwfs]);
	    }
	} 
    }
}
/**
   Setup gradient operator from powfs.loc to wavefront sensor for least square reconstruction.
*/
static void
setup_recon_GWR(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    CALL_ONCE;
    if(!parms->dbg.usegwr) return;
    recon->GWR=spcellnew(parms->npowfs, 1);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].gtype_recon==0){
	    double displace[2]={0,0};
	    recon->GWR->p[ipowfs] = mkg(powfs[ipowfs].gloc, powfs[ipowfs].gloc,
					powfs[ipowfs].gamp->p, powfs[ipowfs].saloc,
					1, 1, displace, 1);
	}else{
	    double displace[2]={0,0};
	    recon->GWR->p[ipowfs] = mkz(powfs[ipowfs].gloc,powfs[ipowfs].gamp->p,
					(loc_t*)powfs[ipowfs].pts, 1,1,displace);
	}
    }
}
/**
   Setup gradient operator from ploc to wavefront sensors.
*/
static void
setup_recon_GP(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    CALL_ONCE;
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
		GP->p[ipowfs]=mkg(ploc,powfs[ipowfs].gloc,powfs[ipowfs].gamp->p,
				  powfs[ipowfs].saloc,1,1,displace,1);
	    }
		break;
	    case 1:{ /*Create ztilt operator from PLOC, using fine sampled
		       powfs.loc as intermediate plane.*/
		double displace[2]={0,0};
		dsp* ZS0=mkz(powfs[ipowfs].gloc,powfs[ipowfs].gamp->p,
			     (loc_t*)powfs[ipowfs].pts, 1,1,displace);
		info2(" Zploc");
		dsp *H=mkh(ploc,powfs[ipowfs].gloc,powfs[ipowfs].gamp->p, 0,0,1,0,0);
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
    /*assign GP for powfs to recon->GP for each wfs */
    recon->GP=GP;
    recon->GP2=spcellnew(nwfs,1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	recon->GP2->p[iwfs]=spref(recon->GP->p[ipowfs]);
    }
}
/**
   Setup gradient operator form aloc for wfs by using GP.
*/
static void
setup_recon_GA(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    CALL_ONCE;
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
	info2("Generating GA ");TIC;tic;
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
		if(!parms->recon.glao){
		    displace[0]=parms->wfsr[iwfs].thetax*ht;
		    displace[1]=parms->wfsr[iwfs].thetay*ht;
		}
		if(parms->dbg.usegwr){
		    warning("todo: Fix and use mkg directly\n");
		    dsp *H=mkh(recon->aloc[idm], powfs[ipowfs].gloc, NULL, 
			       displace[0],displace[1],scale,
			       parms->dm[idm].cubic,parms->dm[idm].iac);
		    GA[idm][iwfs]=spmulsp(recon->GWR->p[ipowfs], H);
		    spfree(H);
		}else{
		    dsp *H=mkh(recon->aloc[idm], ploc, NULL, 
			       displace[0],displace[1],scale,
			       parms->dm[idm].cubic,parms->dm[idm].iac);
		    
		    GA[idm][iwfs]=spmulsp(recon->GP->p[ipowfs], H);
		    spfree(H);
		}
	    }/*idm */
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
    /*Create GAlo that only contains GA for low order wfs */
    spcellfree(recon->GAlo);
    recon->GAlo=spcellnew(recon->GA->nx, recon->GA->ny);
    recon->GAhi=spcellnew(recon->GA->nx, recon->GA->ny);
    PDSPCELL(recon->GAlo,GAlo);
    PDSPCELL(recon->GAhi,GAhi);
    PDSPCELL(recon->GA,GA);
    int nlo=parms->nlopowfs;
    for(int idm=0; idm<ndm; idm++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].lo
	       || (parms->recon.split && nlo==0 && !parms->powfs[ipowfs].trs)){/*for low order wfs */
		GAlo[idm][iwfs]=spref(GA[idm][iwfs]);
	    }
	    if(!parms->powfs[ipowfs].skip){
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
    CALL_ONCE;
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    recon->GX=spcellnew(nwfs, npsr);
    PDSPCELL(recon->GX,GX);
    PDSPCELL(recon->HXW,HXW);
    info2("Generating GX ");TIC;tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	/*gradient from xloc */
	for(int ips=0; ips<npsr; ips++){
	    GX[ips][iwfs]=spmulsp(recon->GP2->p[iwfs], HXW[ips][iwfs]);
	}/*ips */
    }
    toc2(" ");
    spcellfree(recon->GXlo);
    spcellfree(recon->GXtomo);
    spcellfree(recon->GXfocus);
    
    recon->GXtomo=spcellnew(recon->GX->nx, recon->GX->ny);
    PDSPCELL(recon->GXtomo,GXtomo);

    recon->GXlo=spcellnew(recon->GX->nx, recon->GX->ny);
    PDSPCELL(recon->GXlo, GXlo);

    recon->GXfocus=spcellnew(recon->GX->nx, recon->GX->ny);
    PDSPCELL(recon->GXfocus,GXfocus);
    int nlo=parms->nlopowfs;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	for(int ips=0; ips<npsr; ips++){
	    if(!parms->powfs[ipowfs].skip){/*for tomography */
		GXtomo[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	    if(parms->powfs[ipowfs].lo 
	       || (parms->recon.split && nlo==0 && !parms->powfs[ipowfs].trs)){
		/*for low order wfs or extracted t/t for high order ngs wfs.*/
		GXlo[ips][iwfs]=spref(GX[ips][iwfs]);
	    }
	    /*for focus tracking. */
	    if(parms->dbg.deltafocus){
		GXfocus[ips][iwfs]=spref(GX[ips][iwfs]);
		if(!GX[ips][iwfs]){
		    warning("GX[%d][%d] should not be empty\n", ips, iwfs);
		}
	    }
	}
    }/*iwfs */
}
/**
   Setup the matrix of the inverse of gradient measurement noise equivalent
   angle covariance matrix. For physical optics wfs, the NEA is computed using
   matched filter output. For geometric optics, the NEA is from input.
*/
static void
setup_recon_saneai(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    CALL_ONCE;
    const int nwfs=parms->nwfsr;
    spcellfree(recon->sanea);
    spcellfree(recon->saneai);
    spcellfree(recon->saneal);
    spcell *sanea=recon->sanea=spcellnew(nwfs,nwfs);
    spcell *saneal=recon->saneal=spcellnew(nwfs,nwfs);
    spcell *saneai=recon->saneai=spcellnew(nwfs,nwfs);
    info2("Recon NEA:\n");
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs[0];
	int nsa=powfs[ipowfs].pts->nsa;
	int do_ref=0;
	if(parms->powfs[ipowfs].neareconfile){/*taks precedance */
	    dmat *nea=dread("%s_wfs%d",parms->powfs[ipowfs].neareconfile,iwfs);/*rad */
	    if(nea && nea->p[0]<1e-11) {
		error("Supplied NEA is too small\n");
	    }
	    /*rad */
	    saneal->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea, 2);/*rad^2 */
	    sanea->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dcwpow(nea,-1);/*rad^-2 */
	    saneai->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dfree(nea);
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    /*Physical optics */
	    INTSTAT_T *intstat=powfs[ipowfs].intstat;
	    if(intstat->saneaxy){
		const int nnea=intstat->saneaxy->ny;
		if(parms->recon.glao && nnea!=1){
		    error("Please average intstat->saneaxy for GLAO mode.\n");
		}
		int indsanea=0;
		if(nnea==1){
		    indsanea=0;
		}else if(nnea==parms->powfs[ipowfs].nwfs){
		    indsanea=parms->powfs[ipowfs].wfsind[iwfs];
		}else{
		    error("invalid\n");
		}
		if(nnea>1 || iwfs==iwfs0){
		    dmat **saneaxy =intstat->saneaxy ->p+indsanea*nsa;
		    dmat **sanealxy=intstat->saneaxyl->p+indsanea*nsa;
		    dmat **saneaixy=intstat->saneaixy->p+indsanea*nsa;
		    sanea->p[iwfs+iwfs*nwfs]=nea2sp(saneaxy,nsa);
		    saneal->p[iwfs+iwfs*nwfs]=nea2sp(sanealxy,nsa);
		    saneai->p[iwfs+iwfs*nwfs]=nea2sp(saneaixy,nsa);
		}else{
		    do_ref=1;
		}
	    }else{
		error("Not implemented yet\n");
	    }
	}else{
	    /*nea scales as sqrt(1/dtrat) */
	    if(iwfs==iwfs0){
		const double neasq=pow(parms->powfs[ipowfs].nearecon/206265000.,2)/parms->powfs[ipowfs].dtrat;
		if(neasq<1.e-30) error("nea is too small\n");
		dmat *nea=dnew(nsa,2);
		/*scale neasq by area^-1. (seeing limited) */
		/*scale neasq by area^-2 if diffraction limited */
		/*only implementing seeing limited here. */
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
	    }else{
		do_ref=1;
	    }
	}
	if(do_ref){
	    sanea->p[iwfs+iwfs*nwfs] =spref( sanea->p[iwfs0+iwfs0*nwfs]);
	    saneal->p[iwfs+iwfs*nwfs]=spref(saneal->p[iwfs0+iwfs0*nwfs]);
	    saneai->p[iwfs+iwfs*nwfs]=spref(saneai->p[iwfs0+iwfs0*nwfs]);
	}
    }/*iwfs */
    
    /*Compute the averaged SANEA for each WFS */
    recon->neam=dnew(parms->nwfsr, 1);
    double neamhi=0; 
    int counthi=0;
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	const int nsa=powfs[ipowfs].pts->nsa;
	dmat *sanea_iwfs=spdiag(recon->sanea->p[iwfs+iwfs*parms->nwfsr]);
	double area_thres;
	if(nsa>4){
	    area_thres=0.9*parms->powfs[ipowfs].sathruput;
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
	recon->neam->p[iwfs]=sqrt(nea2_sum/count/2);/*average sanea in radian */
	/*double pixtheta=parms->powfs[ipowfs].pixtheta;
	if(recon->neam->p[iwfs]>pixtheta*2 && parms->powfs[ipowfs].usephy){
	//Neglecting WFS whos NEA is greater than twice pixel size in
	//physical optics mode.
	    spfree(recon->saneai->p[iwfs+iwfs*parms->nwfsr]);
	    recon->saneai->p[iwfs+iwfs*parms->nwfsr]=spnewdiag(nsa*2,NULL,0);
	    spfree(recon->saneal->p[iwfs+iwfs*parms->nwfsr]);
	    recon->saneal->p[iwfs+iwfs*parms->nwfsr]=spnewdiag(nsa*2,NULL,0);
	    spfree(recon->sanea->p[iwfs+iwfs*parms->nwfsr]);
	    recon->sanea->p[iwfs+iwfs*parms->nwfsr]=spnewdiag(nsa*2,NULL,INFINITY);
	    recon->neam->p[iwfs]=INFINITY;
	}*/
	char *neatype;
	if(parms->powfs[ipowfs].neareconfile){
	    neatype="FILE";
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    neatype="mtch";
	}else{
	    neatype="geom";
	}
	info2("%s(%.2f) ", neatype, recon->neam->p[iwfs]*206265000);
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
	    warning("why free here crash maos?\n");
	    //dcellfree(powfs[ipowfs].intstat->saneaxy);//why 
	    dcellfree(powfs[ipowfs].intstat->saneaxyl);
	    dcellfree(powfs[ipowfs].intstat->saneaixy);
	}
    }
}
/**
   setting up global tip/tilt remove operator from LGS gradients.
*/

static void
setup_recon_TTR(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    if(!recon->has_ttr) return;
    CALL_ONCE;
    int nwfs=parms->nwfsr;
    recon->TT=dcellnew(nwfs,nwfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].trs 
	   || (parms->recon.split && !parms->powfs[ipowfs].lo)){
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
	    if(parms->recon.glao){
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
setup_recon_DFR(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    if(!recon->has_dfr) return;
    CALL_ONCE;
    if(recon->DF){ 
	dcellfree(recon->DF);
	dcellfree(recon->PDF);
    }
    int nwfs=parms->nwfsr;
    recon->DF=dcellnew(nwfs,nwfs);
    /*Then differential focus modes. */
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].dfrs){
	    if(parms->powfs[ipowfs].nwfs<2){
		error("This powfs group has only 1 wfs. Could not remove diff focus\n");
	    }
	    int nsa=powfs[ipowfs].pts->nsa;
	    dmat* DF=dnew(nsa*2,1);
	    /*saloc is lower left corner of subaperture. don't have to be the center. */
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
setup_recon_TTFR(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    CALL_ONCE;
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
    /*dcellfree(recon->DF);//don't free DF to use in PDF. */
    /*Keep TT, PTT, used in uplink pointing. */
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
   Prepares for tomography. ALlow it to be called multiple times for Cn2 update.
*/
void
setup_recon_tomo_prep(RECON_T *recon, const PARMS_T *parms){
    /*Free existing struct if already exist.  */
    free_cxx(recon);
    if(parms->tomo.assemble){
	/*We need the old copy of L2 when we update the turbulence profile. */
	spcellfree(recon->L2save);
	recon->L2save=recon->L2;
    }else{
	spcellfree(recon->L2);
    }
    recon->L2=NULL;
    /*When layers get a weight less than 1%, we put it at 1% to avoid
      regularization unstability issues.*/
    dclip(recon->wt, 0.01, 1);
    /*normalize the weights to sum to 1. */
    normalize_sum(recon->wt->p, recon->npsr, 1);
    const int npsr=recon->npsr;
    recon->cxx=parms->tomo.cxx;
    /*test_cxx(recon, parms); */
    if(parms->tomo.cxx==0){
	if(parms->load.cxx){
	    recon->L2=spcellread("%s",parms->load.cxx);
	    if(recon->L2->nx!=npsr || recon->L2->ny!=npsr){
		error("Wrong format of loaded L2\n");
	    }
	}else{
	    recon->L2=spcellnew(npsr,npsr);
	    for(int ips=0; ips<npsr; ips++){
		if(parms->tomo.square){/*periodic bc */
		    recon->L2->p[ips+npsr*ips]=mklaplacian_map
			(recon->xmap[ips]->nx, recon->xmap[ips]->nx,
			 recon->xloc[ips]->dx, recon->r0,
			 recon->wt->p[ips]);
		}else{/*reflecive bc */
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
    }
    if(parms->tomo.cxx==1 || (parms->tomo.cxx==2 && parms->tomo.precond==1)){
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
    if(parms->tomo.cxx==2){
	recon->fractal=calloc(1, sizeof(FRACTAL_T));
	recon->fractal->xloc=recon->xloc;
	recon->fractal->r0=parms->atmr.r0;
	recon->fractal->l0=parms->atmr.l0;
	recon->fractal->wt=parms->atmr.wt;
	recon->fractal->scale=parms->tomo.cxxscale;
	recon->fractal->ninit=parms->tomo.ninit;
	dcell *xopd=recon->fractal->xopd=dcellnew(npsr, 1);
	for(int ips=0; ips<npsr; ips++){
	    int nn=nextfftsize(MAX(recon->xmap[ips]->nx, recon->xmap[ips]->ny))+1;
	    xopd->p[ips]=dnew(nn,nn);
	}
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
	    double val=pow(laplacian_coef(r0,wt,dx),2)*1e-6;
	    /*info("Scaling of ZZT is %g\n",val); */
	    /*piston mode eq 47 in Brent 2002 paper */
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

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803

*/
void setup_recon_tomo_matrix(RECON_T *recon, const PARMS_T *parms, APER_T *aper){
    CALL_ONCE;
    /*if not cg or forced, build explicitly the tomography matrix. */
    int npsr=recon->npsr;
    int nwfs=parms->nwfsr;
    /*Free OLD matrices if any. */
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    info2("Before assembling tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);

    if(parms->load.tomo){
	/*right hand side. */
	warning("Loading saved recon->RR\n");
	recon->RR.M=spcellread("RRM");
	if(recon->has_ttr){
	    recon->RR.U=dcellread("RRU");
	    recon->RR.V=dcellread("RRV");
	}
	/*Left hand side */
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
	recon->RR.M=spcellmulspcell(GXtomoT, saneai, TOMOSCALE);
	PDSPCELL(recon->RR.M, RRM);
	/*
	  Tip/tilt and diff focus removal low rand terms for LGS WFS.
	*/
	if(recon->TTF){
	    spcellmulmat(&recon->RR.U, recon->RR.M, recon->TTF, 1);
	    recon->RR.V=dcelltrans(recon->PTTF);
	}
 
	info2("Building recon->RL\n"); /*left hand side matrix */
	recon->RL.M=spcellmulspcell(recon->RR.M,recon->GXtomo,1);
	PDSPCELL(recon->RL.M,RLM);
	if(parms->tomo.piston_cr){ 
	    /*single point piston constraint. no need tikholnov.*/
	    info2("Adding ZZT to RLM\n");
	    for(int ips=0; ips<npsr; ips++){
		spscale(recon->ZZT->p[ips+ips*npsr], TOMOSCALE);
		spadd(&RLM[ips][ips], recon->ZZT->p[ips+ips*npsr]);
	    }
	    spcellfree(recon->ZZT);
	}
	/*Apply tikholnov regularization.*/
	if(fabs(parms->tomo.tikcr)>1.e-15){
	    /*Estimated from the Formula */
	    double maxeig=pow(recon->neamhi * recon->xloc[0]->dx, -2);
	    double tikcr=parms->tomo.tikcr;
	    info2("Adding tikhonov constraint of %g to RLM\n",tikcr);
	    info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	    spcelladdI(recon->RL.M, tikcr*maxeig*TOMOSCALE);
	}
	/*add L2 and ZZT */
	switch(parms->tomo.cxx){
	case 0:/*Add L2'*L2 to RL.M */
	    for(int ips=0; ips<npsr; ips++){
		dsp* tmp=sptmulsp(recon->L2->p[ips+npsr*ips], recon->L2->p[ips+npsr*ips]);
		if(!tmp){
		    error("L2 is empty!!\n");
		}
		spscale(tmp, TOMOSCALE);
		spadd(&RLM[ips][ips], tmp);
		spfree(tmp);
	    }
	    break;
	case 1:/*Need to apply invpsd separately */
	    recon->RL.extra = recon->invpsd;
	    recon->RL.exfun = apply_invpsd;
	    break;
	case 2:/*Need to apply fractal separately */
	    recon->RL.extra = recon->fractal;
	    recon->RL.exfun = apply_fractal;
	}

	/*Symmetricize, remove values below 1e-15*max and sort RLM (optional). */
	/*spcellsym(recon->RL.M); */

	/*Low rank terms for low order wfs. Only in Integrated tomography. */
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
	if(!parms->recon.split || parms->tomo.splitlrt || parms->recon.split==2){
	    recon->RL.U=dcellcat(recon->RR.U, ULo, 2);
	    dcell *GPTTDF=NULL;
	    sptcellmulmat(&GPTTDF, recon->GX, recon->RR.V, 1);
	    recon->RL.V=dcellcat(GPTTDF, VLo, 2);
	    dcellfree(GPTTDF);
	}else{
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	    warning2("Skipping RL Low rank terms in split tomography\n");
	}
	dcellfree(ULo);
	dcellfree(VLo);
	/*Remove empty cells. */
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
	    double val=laplacian_coef(r0,1,dx);/*needs to be a constant */
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

	    spcellwrite(recon->RL.M,"%s/RLM.bin",dirsetup);/*disable compression */
	    dcellwrite(recon->RL.U,"%s/RLU",dirsetup);
	    dcellwrite(recon->RL.V,"%s/RLV",dirsetup); 
	}
	spcellfree(GXtomoT);
    }
    if((parms->tomo.alg==0 || parms->tomo.alg==2) && parms->tomo.bgs ){
	/* We need cholesky decomposition in CBS or MVST method. */
	muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
    }
    if(((parms->tomo.alg==0 || parms->tomo.alg==2) && !parms->tomo.bgs) || parms->sim.ecnn){
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
       	if(recon->RL.C){
	    if(parms->save.recon>1) chol_convert(recon->RL.C, 1);
	    chol_save(recon->RL.C,"%s/RLC.bin",dirsetup);
	}
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
    /*Don't free PTT. Used in forming LGS uplink err */
    info2("After assemble tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}

static void setup_recon_tomo_ecnn(RECON_T *recon, const PARMS_T *parms, APER_T *aper){
    if(!parms->sim.ecnn) return;
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
	/*Luc's Method. Solve E^T Hx^T */
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
	    muv_trans(&p, &recon->RR, t1, 1); dcellfree(t1);
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
/**
   Update assembled tomography matrix with new L2. Called from cn2est when new
   profiles are available.
  */
void setup_recon_tomo_update(RECON_T *recon, const PARMS_T *parms){
    setup_recon_tomo_prep(recon, parms); /*redo L2, invpsd */
#if USE_CUDA
    if(parms->gpu.tomo){
	gpu_update_recon(parms, recon);
    }
#endif
    if(parms->tomo.alg==1&&!parms->tomo.assemble){/*no need to do anything */
	return;
    }
    if(parms->tomo.cxx==0){
	/*Need to adjust RLM with the new L2. */
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
	    dsp *LLdiff=spadd2(LL,LLold,1,-1);/*adjustment to RLM */
	    spadd(&RLM[ips][ips], LLdiff);
	    spfree(LLdiff);
	    spfree(LL);
	    spfree(LLold);
	}
    }

    if(parms->tomo.alg==0 || parms->tomo.alg==2){
	/*We need cholesky decomposition in CBS or MVST method. */
	if(!parms->tomo.bgs){/*Full Matrix */
	    muv_direct_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}else{/*BGS */
	    muv_direct_diag_prep(&(recon->RL), (parms->tomo.alg==2)*parms->tomo.svdthres);
	}
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }

    if(parms->recon.split==2){
	setup_recon_mvst(recon,parms);
    }
}
/**
   Setup ray tracing operator HXF from xloc to aperture ploc along DM fiting directions*/
static void
setup_recon_HXF(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
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
	    double hs=parms->fit.hs[ifit];
	    for(int ips=0; ips<npsr; ips++){
		const double ht = recon->ht->p[ips];
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		HXF[ips][ifit]=mkh(recon->xloc[ips], recon->floc, NULL,
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
		HA[idm][ifit]=mkh(recon->aloc[idm], recon->floc, NULL,
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
    CALL_ONCE;
    const int ndm=parms->ndm;
    if(ndm>=3) warning("Low rank terms for 3 or more dms are not tested\n");
    recon->fitNW=dcellnew(ndm,1);
    double scl=recon->fitscl;
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
		    p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;/*x tilt */
		    p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;/*y tilt */
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
		for(int iloc=0; iloc<nloc; iloc++){
		    p2x[iloc]=recon->aloc[idm]->locx[iloc]*factor;/*x tilt */
		    p2y[iloc]=recon->aloc[idm]->locy[iloc]*factor;/*y tilt */
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

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803 
*/
static void
setup_recon_fit_matrix(RECON_T *recon, const PARMS_T *parms){
    CALL_ONCE;
    const int nfit=parms->fit.nfit;
    const int ndm=parms->ndm;
    if(ndm==0) return;
    spcell *HATc=spcelltrans(recon->HA);
    PDSPCELL(HATc, HAT);
    PDSPCELL(recon->HA,HA);

    info2("Before assembling fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    /*Assemble Fit matrix. */
    int npsr=recon->npsr;
    if(parms->load.fit){
	if(!(zfexist("FRM") && zfexist("FRU") && zfexist("FRV"))){
	    error("FRM, FRU, FRV (.bin) not all exist\n");
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
		    /*notice the sqrt. */
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
	/*Always need FR.U as it is used to do FL.U, FL.V */
	recon->FR.U=dcellnew(ndm, 1);
	dmat **FRU=recon->FR.U->p;
	
	for(int idm=0; idm<ndm; idm++){    
	    int nloc=recon->aloc[idm]->nloc;
	    FRU[idm]=dnew(nloc, nfit);
	    for(int ifit=0; ifit<nfit; ifit++){
		/*notice the sart. */
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

    /*Compute the proper strength of the low rand terms in FLM */
    recon->fitscl=0;
    if(recon->FR.U && recon->FR.U->p){
	dmat **FRU=recon->FR.U->p;
	for(int idm=0; idm<ndm; idm++){
	    /*follow loas convention. */
	    double scl0=maxabs(FRU[idm]->p,FRU[idm]->nx*FRU[idm]->ny);
	    if(scl0>recon->fitscl) recon->fitscl=scl0;
	}
    }else{
	recon->fitscl=1./recon->floc->nloc;/*scale the constraints. Important!! */
    }

    if(parms->load.fit){
	if(!(zfexist("FLM") && zfexist("FLU") && zfexist("FLV"))){
	    error("FLM, FLU, FLV (.bin) not all exist\n");
	}
	warning("Loading saved recon->FL\n");
	recon->FL.M=spcellread("FLM");
	recon->FL.U=dcellread("FLU");
	recon->FL.V=dcellread("FLV");
    }else{
	fit_prep_lrt(recon,parms);
	if(parms->fit.actslave){
	    /*
	      2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found
	      that slaving is causing mis-measurement of a few edge
	      actuators. First try to remove W1. Or lower the weight. Revert
	      back.
	     */
	    recon->actslave=slaving(&recon->actcpl, recon->aloc, recon->HA, recon->W1,
				    recon->fitNW, recon->actstuck,
				    recon->actfloat, 0.1, 1./recon->floc->nloc);
	    if(parms->save.setup){
		dcellwrite(recon->actcpl, "%s/actcpl", dirsetup);
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

	{/*Low rank terms. */
	    recon->FL.U=dcellcat_each(recon->FR.U, recon->fitNW, 2);
	    dcell *tmp=NULL;/*negative NW. */
	    dcelladd(&tmp, 1, recon->fitNW, -1);
	    recon->FL.V=dcellcat_each(recon->FR.U, tmp, 2);
	    dcellfree(tmp);
	}
	if(recon->actslave){
	    spcelladd(&recon->FL.M, recon->actslave);
	}
	/*spcellsym(recon->FL.M); */
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
    if((parms->fit.alg==0 || parms->fit.alg==2) && !parms->fit.bgs){
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
setup_recon_focus(RECON_T *recon, POWFS_T *powfs, const PARMS_T *parms){
    CALL_ONCE;
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
    /*Create GFall: Focus mode -> WFS grad */
    recon->GFall=dcellnew(parms->nwfsr, 1);
    recon->GFlgs=dcellnew(parms->nwfsr, 1);
    recon->GFngs=dcellnew(parms->nwfsr, 1);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	dmat *opd=dnew(powfs[ipowfs].loc->nloc,1);
	loc_add_focus(opd->p, powfs[ipowfs].loc, 1);
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	if(parms->powfs[ipowfs].gtype_recon==1){
	    dcell *saimcc;
	    if(powfs[ipowfs].locm){/*We can not use realamp as in simulation. */
		dcell *mcc=pts_mcc_ptt(powfs[ipowfs].pts, powfs[ipowfs].amp->p);
		saimcc=dcellinvspd_each(mcc);
		dcellfree(mcc);
	    }else{
		saimcc=dcellref(powfs[ipowfs].saimcc[0]);
	    }
	    pts_ztilt(&recon->GFall->p[iwfs], powfs[ipowfs].pts, saimcc,
		      powfs[ipowfs].realamp[wfsind], opd->p);
	}else{
	    spmulmat(&recon->GFall->p[iwfs], adpind(powfs[ipowfs].GS0,wfsind), opd, 1);
	}
	if(parms->powfs[ipowfs].llt){
	    recon->GFlgs->p[iwfs]=dref(recon->GFall->p[iwfs]);
	}else{
	    recon->GFngs->p[iwfs]=dref(recon->GFall->p[iwfs]);
	}
	dfree(opd);
    }
    if(parms->save.setup){
	dcellwrite(recon->GFall,"%s/GFall",dirsetup);
    }
    if(!parms->sim.mffocus) return;

    dmat *GMGngs=NULL;
    dcell *GMngs=dcellnew(1, parms->nwfsr);
    /*Compute focus reconstructor from NGS Grads. fuse grads
      together to construct a single focus measurement*/
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].trs==0 && parms->powfs[ipowfs].order>1){
	    info2("wfs %d will be used to track focus\n", iwfs);
	}else{
	    continue;
	}
	spmulmat(&GMngs->p[iwfs], recon->saneai->p[iwfs+parms->nwfsr*iwfs], 
		 recon->GFall->p[iwfs],1);
	dmm(&GMGngs,recon->GFall->p[iwfs], GMngs->p[iwfs], "tn",1);
    }
    dinvspd_inplace(GMGngs);
    /*A focus reconstructor from all NGS measurements.*/
    dcell *RFngsg=recon->RFngsg=dcellnew(1, parms->nwfsr);
  
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	if(!recon->GFall->p[iwfs]) continue;
	//NGS gradient to Focus mode reconstructor.
	dmm(&RFngsg->p[iwfs],GMGngs, GMngs->p[iwfs],"nt",1);
    }
    dfree(GMGngs);
    dcellfree(GMngs);
    if(parms->dbg.deltafocus){
	/*focus reconstruction for science from opdr using opd best fit.*/
	dmat *Fsci=dnew(recon->floc->nloc, 1);
	loc_add_focus(Fsci->p, recon->floc, 1);
	dmat *WF=NULL;/*WF=(W0-W1*W1')*Fsci*/
	spmulmat(&WF, recon->W0, Fsci, 1);
	dmat *W1F=NULL;
	dmm(&W1F, recon->W1, Fsci, "tn", 1);
	dmm(&WF,  recon->W1, W1F, "nn", -1);
	dmat *FtWF=NULL;
	dmm(&FtWF, Fsci, WF, "tn", 1);
	dinvspd_inplace(FtWF);
	dmat *RFsci=NULL;
	dmm(&RFsci, FtWF, WF, "nt", 1);
	dfree(FtWF); dfree(W1F); dfree(WF); dfree(Fsci); 
	int ifit=parms->fit.indoa;//center fit direction
	//Find the TTF OIWFS
	int iwfs_ttf=-1;
	int ipowfs_ttf=-1;
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].trs==0 && parms->powfs[ipowfs].order>1){
		if(iwfs_ttf==-1){
		    iwfs_ttf=iwfs;
		    ipowfs_ttf=ipowfs;
		}else{
		    error("There are multiple TTF. Revise\n");
		}
	    }
	}
	if(iwfs_ttf==-1){
	    error("There is no TTF.\n");
	}
	//Compute ray tracing operator for TTF WFS.
	double hs=parms->powfs[ipowfs_ttf].hs;
	recon->RFdfx=dcellnew(1, recon->npsr);
	recon->RFdfa=dcellnew(1, parms->ndm);
	PDSPCELL(recon->HXF,HXF);
	PDSPCELL(recon->HA, HA);
	for(int ips=0; ips<recon->npsr; ips++){
	    const double ht = recon->ht->p[ips];
	    const double scale=1.-ht/hs;
	    double displace[2];
	    displace[0]=parms->wfsr[iwfs_ttf].thetax*ht;
	    displace[1]=parms->wfsr[iwfs_ttf].thetay*ht;
	    dsp *HX_TTF=mkh(recon->xloc[ips], recon->floc, NULL,
			    displace[0], displace[1], scale,
			    parms->tomo.cubic, parms->tomo.iac);
	    dmulsp(&recon->RFdfx->p[ips], RFsci, HXF[ips][ifit], 1);
	    dmulsp(&recon->RFdfx->p[ips], RFsci, HX_TTF, -1);
	    spfree(HX_TTF);
	}
	for(int idm=0; idm<parms->ndm; idm++){
	    const double ht=parms->dm[idm].ht;
	    const double scale=1.-ht/hs;
	    double displace[2];
	    displace[0]=parms->wfsr[iwfs_ttf].thetax*ht;
	    displace[1]=parms->wfsr[iwfs_ttf].thetay*ht;
	    dsp *HA_TTF=mkh(recon->aloc[idm], recon->floc, NULL,
			    displace[0], displace[1], scale,
			    parms->dm[idm].cubic,parms->dm[idm].iac);
	    dmulsp(&recon->RFdfa->p[idm], RFsci, HA[idm][ifit], 1);
	    dmulsp(&recon->RFdfa->p[idm], RFsci, HA_TTF, -1);
	    spfree(HA_TTF);
	}
	dfree(RFsci);
    }
    /*
      Compute focus constructor from LGS grads. A constructor for each LGS
      because each LGS may have different range error.
    */
    recon->RFlgsg=dcellnew(parms->nwfsr, parms->nwfsr);
    PDCELL(recon->RFlgsg, RFlgsg);
  
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].llt)
	    continue;
	dmat *GMtmp=NULL;
	dmat *GMGtmp=NULL;
	spmulmat(&GMtmp, recon->saneai->p[iwfs+parms->nwfsr*iwfs], 
		 recon->GFall->p[iwfs], 1);
	dmm(&GMGtmp, recon->GFall->p[iwfs], GMtmp, "tn",1);
	dinvspd_inplace(GMGtmp);
	dmm(&RFlgsg[iwfs][iwfs], GMGtmp, GMtmp, "nt", 1);

	dfree(GMtmp);
	dfree(GMGtmp);
    }
  
    if(parms->save.setup){
	dcellwrite(recon->RFngsg,"%s/RFngsg",dirsetup);
	dcellwrite(recon->RFlgsg,"%s/RFlgsg",dirsetup);
	if(parms->dbg.deltafocus){
	    dcellwrite(recon->RFdfx,"%s/RFdfx",dirsetup);
	}
    }
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
   Uw&=&A^{-1}G_{ngs}^TC_{ngs}^{-1}\\
   \hat{x}_{ngs}&=&Uw(1+G_{ngs}Uw)^{-1}(s_{ngs}-G_{ngs}\hat{x}_{lgs});\\
   a_{NGS}&=&F\hat{x}_{ngs}\\
   MVModes&=&F\cdot Uw\\
   MVRngs&=&(1+G_{ngs}Uw)^{-1}
   \f}

   If we want to orthnormalize the NGS modes. Propagate the NGS modes \f$F\cdot Uw\f$ to
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
    CALL_ONCE;
    TIC;tic;
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
    if(parms->recon.split!=2){
	return;
    }
 
    dcellfree(recon->GXL);
    spcellfull(&recon->GXL, recon->GXlo, 1);
    //NEA of low order WFS.
    dcell *neailo=dcellnew(parms->nwfsr, parms->nwfsr);
    dcell *nealo=dcellnew(parms->nwfsr, parms->nwfsr);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].lo){
	    spfull(&neailo->p[iwfs*(1+parms->nwfsr)],
		   recon->saneai->p[iwfs*(1+parms->nwfsr)], TOMOSCALE);
	    spfull(&nealo->p[iwfs*(1+parms->nwfsr)],
		   recon->sanea->p[iwfs*(1+parms->nwfsr)], 1);
	}
    }
    /* 2012-03-21: Remove focus mode from GL and NEA so no focus is
       measured/estimated by mvst. It is then estimate separately*/
    if(parms->sim.mffocus){
	dcell *focus=NULL;
	dcellmm(&focus, recon->RFngsg, recon->GXL, "nn", 1);
	dcellmm(&recon->GXL, recon->GFngs, focus, "nn", -1);
	dcelldropzero(recon->GXL, 1e-12);
	dcellfree(focus);
	dcellmm(&focus, recon->RFngsg, neailo, "nn", 1);
	dcellmm(&neailo, recon->GFngs, focus, "nn", -1);
	dcelldropzero(neailo, 1e-8);
	dcellfree(focus);
    }

    dcell *U=NULL; 
    dcell *FU=NULL;
    if(parms->load.mvst){
	U=dcellread("mvst_U");
	FU=dcellread("mvst_FU");
    }else{
	/*Prepare CBS if not already done */
	if(!recon->RL.C && !recon->RL.MI){
	    muv_direct_prep(&(recon->RL), 0);
	}
	if(!recon->FL.C && !recon->FL.MI){
	    muv_direct_prep(&(recon->FL), 0);
	}
	toc2("MVST: svd prep");
	dcell *GXLT=dcelltrans(recon->GXL);
	muv_direct_solve_cell(&U, &recon->RL, GXLT);
	dcellfree(GXLT);
	dcell *rhs=NULL;
	muv(&rhs, &recon->FR, U, 1);
	muv_direct_solve_cell(&FU, &recon->FL, rhs);
	dcellfree(rhs);
	toc2("MVST: U, FU");
    
	if(parms->save.mvst || parms->save.setup){
	    dcellwrite(U, "%s/mvst_U", dirsetup);
	    dcellwrite(FU, "%s/mvst_FU", dirsetup);
	}
    }
    dcell *Uw=NULL;
    dcell *FUw=NULL;

    dcellmm(&Uw, U, neailo, "nn", 1);
    dcellmm(&FUw, FU, neailo, "nn", 1);

    dcell *M=NULL;
    dcellmm(&M, recon->GXL, Uw, "nn", 1);
    dcelladdI(M, 1);
    dcell *Minv=dcellinv(M);
    dcellfree(M);
    if(parms->sim.mffocus){
	//Embed a focus removal. Necessary!
	dcell *focus=NULL;
	dcellmm(&focus, Minv, recon->GFngs, "nn", 1);
	dcellmm(&Minv, focus, recon->RFngsg, "nn", -1);
    }
    if(parms->save.setup){
	dcellwrite(Minv, "%s/mvst_Rngs_0",dirsetup);
	dcellwrite(FUw,  "%s/mvst_Modes_0",dirsetup);
    }
    /*Compute the cross coupling matrix of the Modes:
      FUw'*Ha'*W*Fuw*Ha. Re-verified on 2013-03-24.*/
    dcell *QwQc=NULL;
    {
	dcell *Q=NULL;/*the NGS modes in ploc. */
	spcellmulmat(&Q, recon->HA, FUw, 1);
	QwQc=calcWmcc(Q,Q,recon->W0,recon->W1,recon->fitwt);
	dcellfree(Q);
    }
    /*Compute the wavefront error due to measurement noise. Verified on
      2013-03-24. The gain optimization is yet a temporary hack because the way
      PSDs are input.*/
    if(1){
	dcell *RC=NULL;
	dcellmm(&RC, Minv, nealo, "nn", 1);
	dcell *RCRt=NULL;
	dcellmm(&RCRt, RC, Minv, "nt", 1);
	dcell *RCRtQwQ=NULL;
	dcellmm(&RCRtQwQ, RCRt, QwQc, "nn", 1);
	dmat *tmp=dcell2m(RCRtQwQ);
	PDMAT(tmp, ptmp);
	double rss=0;
	for(int i=0; i<tmp->nx; i++){
	    rss+=ptmp[i][i];
	}
	dfree(tmp);
	dcellfree(RCRtQwQ);
	dcellfree(RCRt);
	dcellfree(RC);
	recon->sigmanlo=rss;
	info("rms=%g nm\n", sqrt(rss)*1e9);

	warning("Temporary solution for testing\n");
	dmat *psd_ngs=dread("../../psd_ngs.bin");
	if(1){
	    //need to convert from rad to m2.
	    dmat *psd_ws=dread("%s", parms->sim.wspsd);
	    dmat *psd_ws_m=ddup(psd_ws); 
	    dfree(psd_ws);
	    dmat *psd_ws_y=dnew_ref(psd_ws_m->nx,1,psd_ws_m->p+psd_ws_m->nx);
	    dscale(psd_ws_y, 4./parms->aper.d); dfree(psd_ws_y);
	    add_psd2(&psd_ngs, psd_ws_m); dfree(psd_ws_m);
	}
	dwrite(psd_ngs, "psd_ngs_servo");
	dmat *rss2=dnew(1,1); rss2->p[0]=rss;
	int dtrat=parms->powfs[parms->lopowfs[0]].dtrat;
	dcell *res=servo_optim(psd_ngs, parms->sim.dt, 
			       dtrat, M_PI/4, rss2, 2); 
	dfree(rss2);
	info("dtrat=%d\n", dtrat);
	info("g,a,T was %g,%g,%g\n", parms->sim.eplo->p[0], parms->sim.eplo->p[1], parms->sim.eplo->p[2]);
	memcpy(parms->sim.eplo->p, res->p[0]->p, 3*sizeof(double));
	info("g,a,T=%g,%g,%g\n", parms->sim.eplo->p[0], parms->sim.eplo->p[1], parms->sim.eplo->p[2]);
	info("res=%g, resn=%g nm\n", sqrt(res->p[0]->p[3])*1e9, sqrt(res->p[0]->p[4])*1e9);
	dcellfree(res); 
	dfree(psd_ngs);
    }
    if(1){/*Orthnormalize the Modes. */
	/*
	  Change FUw*Minv -> FUw*(U*sigma^-1/2) * (U*sigma^1/2)'*Minv
	  columes of FUw*(U*sigma^-1/2) are the eigen vectors.

	  U, sigma is the eigen value decomposition of <FUw' HA' W HA FUw>
	*/

	dmat *QSdiag=NULL, *QU=NULL, *QVt=NULL;
	{
	    dmat *QwQ=dcell2m(QwQc);
	    dsvd(&QU, &QSdiag, &QVt, QwQ);
	    dfree(QwQ);
	}
	if(parms->save.setup) {
	    dwrite(QSdiag,"%s/mvst_QSdiag",dirsetup);
	}
	dcwpow(QSdiag, -1./2.);
	dmuldiag(QU,QSdiag);/*U*sigma^-1/2 */
	d2cell(&QwQc,QU,NULL);
	dcell *FUw_keep=FUw;FUw=NULL; 
	dcellmm(&FUw, FUw_keep, QwQc, "nn", 1);
	dcellfree(FUw_keep);
	dcwpow(QSdiag,-2);
	dmuldiag(QU,QSdiag);/*U*sigma^1/2 (From U*sigma^(-1/2)*sigma) */
	d2cell(&QwQc,QU,NULL);
	dcell *Minv_keep=Minv; Minv=NULL;
	dcellmm(&Minv,QwQc,Minv_keep,"tn", 1);
	dcellfree(Minv_keep);
	dfree(QSdiag);
	dfree(QVt);
	dfree(QU);
    }
    dcellfree(QwQc);
    
    recon->MVRngs=dcellreduce(Minv,1);/*1xnwfs cell */
    recon->MVModes=dcellreduce(FUw,2);/*ndmx1 cell */
    recon->MVGM=NULL;
    spcellmulmat(&recon->MVGM, recon->GAlo, recon->MVModes, 1);
    dcellfree(neailo);
    dcellfree(nealo);
    dcellfree(Minv);
    dcellfree(U);
    dcellfree(FU);
    dcellfree(FUw);
    dcellfree(Uw);
    if(parms->save.setup){
	dcell *Qn=NULL;
	spcellmulmat(&Qn, recon->HA, recon->MVModes, 1);
	dcell *Qntt=dcellnew(Qn->nx,Qn->ny);
	dmat *TTploc=loc2mat(recon->floc,1);/*TT mode. need piston mode too! */
	dmat *PTTploc=dpinv(TTploc,NULL,recon->W0);/*TT projector. no need w1 since we have piston. */
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
    if(parms->dbg.mvstlimit>0){/*limit number of modes used. */
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
    toc2("MVST");
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
static void setup_recon_mvr_fit(RECON_T *recon, const PARMS_T *parms, APER_T *aper){
    CALL_ONCE;
    TIC;tic;
    if(parms->recon.alg==0 && !parms->sim.idealfit){/*In idealfit, xloc has high sampling. We avoid HXF. */
	setup_recon_HXF(recon,parms);
    }
    setup_recon_HA(recon,parms);
    
    /*copy over fitwt since we need a dmat */
    int nfit=parms->fit.nfit;
    recon->fitwt=dnew(nfit,1);
    memcpy(recon->fitwt->p,parms->fit.wt,sizeof(double)*nfit);
    
    /*always assemble fit matrix, faster if many directions */
    if(parms->fit.assemble){
	setup_recon_fit_matrix(recon,parms);
    }
    toc2("Generating fit matrix ");
    /*Fall back function method if FR.M is NULL (!HXF<-idealfit) */
    recon->FR.Mfun  = FitR;
    recon->FR.Mdata = recon;
    /*Fall back function method if FL.M is NULL */
    recon->FL.Mfun  = FitL;
    recon->FL.Mdata = recon;
    recon->FL.alg = parms->fit.alg;
    recon->FL.bgs = parms->fit.bgs;
    recon->FL.warm  = parms->recon.warm_restart;
    recon->FL.maxit = parms->fit.maxit;
}
void setup_recon_mvr(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    CALL_ONCE;
    TIC;tic;
    if(parms->cn2.tomo){
	/*Use cn2 estimation results for tomography. Use its ht to build
	  reconstructor matrices.*/
	CN2EST_T *cn2est=recon->cn2est;
	recon->ht=dref(cn2est->htrecon);
	recon->os=dref(cn2est->os);
	recon->wt=dref(cn2est->wtrecon->p[0]);
	/*the following will be updated later in simulation. */
	dset(recon->wt, 1./recon->wt->nx);/*evenly distributed.  */
	recon->r0=0.15;/*random guess */
	recon->l0=30;/*random guess */
    }else{/*use input information from atmr */
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
    /*sampling of xloc */
    recon->dx=dnew(recon->ht->nx, 1);
    for(int iht=0; iht<recon->ht->nx; iht++){
	double scale = 1.0 - recon->ht->p[iht]/parms->atmr.hs;
	recon->dx->p[iht]=(parms->atmr.dx/recon->os->p[iht])*scale;	
    }
    /*number of reconstruction layers */
    recon->npsr= recon->ht->nx;
    /*setup atm reconstruction layer grid */
    setup_recon_xloc(recon,parms);
    /*setup xloc/aloc to WFS grad */
    if(!parms->sim.idealfit){
	setup_recon_HXW(recon,parms);
	setup_recon_GX(recon,parms);
	spcellfree(recon->HXW);/*only keep HXWtomo for tomography */
	/*setup inverse noise covariance matrix. */
	/*prepare for tomography setup */
	setup_recon_tomo_prep(recon,parms);
	if(parms->tomo.assemble || parms->recon.split==2 || parms->sim.psfr){
	    /*assemble the matrix only if not using CG CG apply the
	      individual matrices on fly to speed up and save memory. */
	    setup_recon_tomo_matrix(recon,parms,aper);
	}
	if(parms->sim.ecnn){
	    setup_recon_tomo_ecnn(recon, parms, aper);
	}
	/*Fall back function method if .M is NULL */
	recon->RL.Mfun=TomoL;
	recon->RL.Mdata=recon;
	recon->RR.Mfun=TomoR;
	recon->RR.Mtfun=TomoRt;
	recon->RR.Mdata=recon;
	if(parms->tomo.alg==1){/*CG */
	    switch(parms->tomo.precond){
	    case 0:/*no preconditioner */
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
	recon->RL.warm  = parms->recon.warm_restart;
	recon->RL.maxit = parms->tomo.maxit;
    }
    setup_recon_mvr_fit(recon, parms, aper);
    /*moao */
    setup_recon_moao(recon,parms);
 
    if(parms->recon.split){
	/*split tomography */
	if(parms->recon.split){
	    if(parms->ndm<=2){
		/*setup the ngsmode in both ahst and mvst mode  */
		setup_ngsmod(parms,recon,aper,powfs);
	    }else if(parms->recon.split==1){
		error("Not implemented");
	    }
	}
	if(!parms->sim.idealfit && parms->recon.split==2){/*Need to be after fit */
	    setup_recon_mvst(recon, parms);
	}
    }
    /*todo: can we move to near setting up tomography? */
    if(!parms->sim.idealfit && parms->tomo.precond==1){
	recon->fdpcg=fdpcg_prepare(parms, recon, powfs);
    }

    toc2("setup_recon_mvr");
}
/**
   Setup the least square reconstruct by directly inverting GA matrix. 
   The reconstructor is simply the pseudo inverse of GA matrix:
   \f[\hat{x}=(G_a^TC_g^{-1}G_a)^{-1}G_a^TC_g^{-1}\f]

   This is very close to RR except replacing GX with GA.

   We use the tomogrpahy parameters for lsr, since lsr is simply "tomography" onto DM directly.
*/
void setup_recon_lsr(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    CALL_ONCE;
    spcell *GAlsr;
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;

    if(parms->recon.split){
	/*high order wfs only in split mode. */
	GAlsr=recon->GAhi;
	setup_ngsmod(parms,recon,aper,powfs);
    }else{
	/*all wfs in integrated mode. */
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
    if(fabs(parms->lsr.tikcr)>EPS){
	info2("Adding tikhonov constraint of %g to LLM\n", parms->lsr.tikcr);
	info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	spcelladdI(recon->LL.M, parms->lsr.tikcr*maxeig);
    }
    dcell *NW=NULL;
    if(parms->lsr.alg!=2){
	/* Not SVD, need low rank terms for piston/waffle mode constraint. */
	NW=dcellnew(ndm,1);
	int nmod=2;/*two modes. */
	for(int idm=0; idm<ndm; idm++){
	    loc_create_map(recon->aloc[idm]);
	    const long nloc=recon->aloc[idm]->nloc;
	    NW->p[idm]=dnew(nloc, ndm*nmod);
	    double *p=NW->p[idm]->p+nmod*idm*nloc;
	    for(long iloc=0; iloc<nloc; iloc++){
		p[iloc]=1;/*piston mode */
	    }
	    /*notice offset of 1 because map start count at 1 */
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
	/*scale it to match the magnitude of LL.M */
	dcellscale(NW, sqrt(maxeig));
	if(parms->save.setup){
	    dcellwrite(NW, "%s/lsrNW",dirsetup);
	}
    }
    if(parms->lsr.actslave){
	/*actuator slaving. important. change from 0.5 to 0.1 on 2011-07-14. */
	spcell *actslave=slaving(&recon->actcpl, recon->aloc, recon->GAhi, NULL, NW,
				 recon->actstuck, recon->actfloat, 0.1, sqrt(maxeig));
	if(parms->save.setup){
	    if(NW){
		dcellwrite(NW, "%s/lsrNW2",dirsetup);
	    }
	    spcellwrite(actslave,"%s/actslave", dirsetup);
	}
	spcelladd(&recon->LL.M, actslave);
	spcellfree(actslave);
    }
    /*Low rank terms for low order wfs. Only in Integrated tomography. */
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
	/*add to low rank terms. */
	dcell *tmp=recon->LL.U;
	recon->LL.U=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellscale(NW, -1);
	tmp=recon->LL.V;
	recon->LL.V=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellfree(NW);
    }
    if(parms->lsr.fnreg){
	warning("Loading LSR regularization from file %s.\n", parms->lsr.fnreg);
	spcell *tmp=spcellread("%s", parms->lsr.fnreg);
	spcelladd(&recon->LL.M, tmp);
	spcellfree(tmp);
    }
    recon->LL.alg = parms->lsr.alg;
    recon->LL.bgs = parms->lsr.bgs;
    recon->LL.warm = parms->recon.warm_restart;
    recon->LL.maxit = parms->lsr.maxit;
    /*Remove empty cells. */
    dcelldropempty(&recon->LR.U,2);
    dcelldropempty(&recon->LR.V,2);
    dcelldropempty(&recon->LL.U,2);
    dcelldropempty(&recon->LL.V,2);
    if(parms->save.recon){
	spcellwrite(recon->LR.M,"%s/LRM",dirsetup);
	dcellwrite(recon->LR.U,"%s/LRU",dirsetup);
	dcellwrite(recon->LR.V,"%s/LRV",dirsetup);
	spcellwrite(recon->LL.M,"%s/LLM.bin",dirsetup);/*disable compression */
	dcellwrite(recon->LL.U,"%s/LLU",dirsetup);
	dcellwrite(recon->LL.V,"%s/LLV",dirsetup); 
    }
    if(parms->lsr.alg==0 || parms->lsr.alg==2){
	if(!parms->lsr.bgs){
	    muv_direct_prep(&recon->LL, (parms->lsr.alg==2)*parms->lsr.svdthres);
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
	    muv_direct_diag_prep(&(recon->LL), (parms->lsr.alg==2)*parms->lsr.svdthres);
	    if(parms->save.recon){
		for(int ib=0; ib<recon->LL.nb; ib++){
		    if(recon->LL.CB)
			chol_save(recon->LL.CB[ib],"%s/LLCB_%d.bin",dirsetup, ib);
		    else
			dwrite(recon->LL.MI,"%s/LLMIB_%d.bin", dirsetup, ib);
		}
	    }
	    /*Don't free M, U, V */
	}
    }
}
/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
void setup_recon_lsr_mvm(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    info2("Assembling LSR MVM in CPU\n");
    if(recon->LR.Mfun || parms->lsr.alg==1){
	/*
	  First create an identity matrix. then solve each column one by one. 
	*/
	const int ndm=parms->ndm;
	const int nwfs=parms->nwfsr;
	int ntotgrad=0;
	long *ngrad=recon->ngrad;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    ntotgrad+=powfs[ipowfs].saloc->nloc*2;
	}
	recon->MVM=dcellnew(ndm, nwfs);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(!parms->powfs[ipowfs].skip){
		for(int idm=0; idm<ndm; idm++){
		    recon->MVM->p[idm+ndm*iwfs]=dnew(recon->anloc[idm], powfs[ipowfs].saloc->nloc*2);
		}
	    }
	}

	dcell *res=NULL;
	int curg=0, curwfs=0;
	dmat *eye=dnew(ntotgrad, 1);
	dcell *eyec=d2cellref(eye, ngrad, nwfs);
	for(int ig=0; ig<ntotgrad; ig++){
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", ig, ntotgrad);
	    }else if(ig%100==0){
		info2("%6d of %6d\n", ig, ntotgrad);
	    }
	    if(ig) eye->p[ig-1]=0; eye->p[ig]=1;
	    if(!parms->powfs[parms->wfsr[curwfs].powfs].skip){
		dcellzero(res);
		muv_solve(&res, &recon->LL, &recon->LR, eyec);
	    }
	    for(int idm=0; idm<ndm; idm++){
		dmat *to=recon->MVM->p[idm+curwfs*ndm];
		if(to){
		    int nact=to->nx;
		    memcpy(to->p+curg*nact,res->p[idm]->p, nact*sizeof(double));
		}
	    }
	    curg++;
	    if(curg>=ngrad[curwfs]){
		curwfs++;
		curg=0;
	    }
	}
	dcellfree(res);
	dcellfree(eyec);
	dfree(eye);
    }else{
	dcell *LR=NULL;
	spcellfull(&LR, recon->LR.M, 1);
	if(recon->LR.U && recon->LR.V){
	    dcellmm(&LR, recon->LR.U, recon->LR.V, "nt", -1);
	}
	muv_solve(&recon->MVM, &recon->LL, NULL,  LR);
	dcellfree(LR);
    }
}
typedef struct {
    const PARMS_T *parms;
    RECON_T *recon;
    dcell *MVMt;
    long (*curp)[2];
    long ntotact;
}MVR_MVM_T;
void setup_recon_mvr_mvm_iact(thread_t *info){
    MVR_MVM_T *data=info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;
    const long ntotact=data->ntotact;
    dcell *FLI=NULL;
    dcell *FRT=NULL;
    dcell *RLT=NULL;
    dcell *RRT=NULL;
    dmat *eye=dnew(ntotact, 1);
    dcell *eyec=d2cellref(eye, recon->anloc, ndm);
    long (*curp)[2]=data->curp;
    dcell *MVMt=data->MVMt;
    int nthread=recon->nthread;
    for(long iact=info->start; iact<info->end; iact++){
	int curdm=curp[iact][0];
	int curact=curp[iact][1];
	if(recon->actcpl && recon->actcpl->p[curdm]->p[curact]<EPS){
	    continue;
	}
	if(info->ithread==0){
	    if(!detached){
		info2("%6ld of %6ld\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", iact*nthread, ntotact);
	    }else if(iact % 100==0){
		info2("%6ld of %6ld\n", iact*nthread, ntotact);
	    }
	}
	//TIC;tic;
	dcellzero(FRT);
	dcellzero(RRT);
	/*Apply F_L*/
	eye->p[iact]=1;
	muv_solve(&FLI, &recon->FL, NULL, eyec);
	eye->p[iact]=0;
	/*Apply F_R'*/
	muv_trans(&FRT, &recon->FR, FLI, 1);
	//toc2("fit");
	/*Apply R_L*/
	dcellzero(RLT);//warm restart.
	muv_solve(&RLT, &recon->RL, NULL, FRT);
	/*Apply R_R'*/
	muv_trans(&RRT, &recon->RR, RLT, 1);
	//toc2("tomo");
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    dmat *to=MVMt->p[iwfs+curdm*nwfs];
	    if(to){
		int ng=to->nx;
		memcpy(to->p+curact*ng, RRT->p[iwfs]->p, ng*sizeof(double));
	    }
	}
	/*{
	    dcellwrite(FLI, "cpu_dmfit_%ld", iact);
	    dcellwrite(FRT, "cpu_opdx_%ld", iact);
	    dcellwrite(RLT, "cpu_opdr_%ld", iact);
	    dcellwrite(RRT, "cpu_grad_%ld", iact);
	    }*/
    }
    dcellfree(FLI);
    dcellfree(FRT);
    dcellfree(RLT);
    dcellfree(RRT);
    dfree(eye);
    dcellfree(eyec);
}
/**
   Use the various algorithms recon.alg to assemble a final matrix to multiply
   to gradients to get DM commands.
 */
void setup_recon_mvr_mvm(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    info2("Assembling MVR MVM in CPU\n");
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;
    long ntotact=0;
    for(int idm=0; idm<ndm; idm++){
	ntotact+=recon->anloc[idm];
    }
    long (*curp)[2]=malloc(ntotact*2*sizeof(long));
    int nact=0;
    for(int idm=0; idm<ndm; idm++){
	for(int iact=0; iact<recon->anloc[idm]; iact++){
	    curp[nact+iact][0]=idm;
	    curp[nact+iact][1]=iact;
	}
	nact+=recon->anloc[idm];
    }
    dcell *MVMt=dcellnew(nwfs, ndm);
    for(int idm=0; idm<ndm; idm++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(!parms->powfs[ipowfs].skip){
		MVMt->p[iwfs+idm*nwfs]=dnew(powfs[ipowfs].saloc->nloc*2, recon->anloc[idm]);
	    }
	}
    }
    MVR_MVM_T data={parms, recon, MVMt, curp, ntotact};
    int nthread=recon->nthread;
    thread_t info[nthread];
    thread_prep(info, 0, ntotact, nthread, setup_recon_mvr_mvm_iact, &data);
    CALL_THREAD(info, nthread, 1);
    recon->MVM=dcelltrans(MVMt);
    free(curp);
    dcellfree(MVMt);
}
/**
   Setup either the minimum variance reconstructor by calling setup_recon_mvr()
   or least square reconstructor by calling setup_recon_lsr() */
RECON_T *setup_recon(const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    CALL_ONCE;
    TIC;tic;
    RECON_T * recon = calloc(1, sizeof(RECON_T));
    recon->parms=parms;/*save a pointer. */
    recon->aper=aper;
    recon->powfs=powfs;
    if(parms->cn2.npair){
	/*setup CN2 Estimator. It determines the reconstructed layer heigh can
	  be fed to the tomography */
	recon->cn2est=cn2est_prepare(parms,powfs);
    }
    /*number of deformable mirrors */
    recon->ndm = parms->ndm;
    if(parms->recon.warm_restart){
	info2("Using warm restart\n");
    }else{
	warning2("Do not use warm restart\n");
    }
    /*to be used in tomography. */
    recon->nthread=parms->sim.nthread;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	/*we will remove tip/tilt from the high order NGS wfs in split tomo mode.*/
	if(parms->powfs[ipowfs].trs || (parms->recon.split && !parms->powfs[ipowfs].lo)){
	    recon->has_ttr=1;
	    break;
	}
    }
    if(!parms->recon.glao){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs<=1) continue;
	    if(parms->powfs[ipowfs].dfrs){
		recon->has_dfr=1;
		break;
	    }
	}   
    }
    recon->ngrad=calloc(parms->nwfsr, sizeof(long));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){
	    recon->ngrad[iwfs]=powfs[ipowfs].saloc->nloc*2;
	}
    }
    /*setup DM actuator grid */
    setup_recon_aloc(recon,parms);
    /*setup pupil coarse grid */
    setup_recon_ploc(recon,parms);
    setup_recon_floc(recon,parms);
    setup_recon_GWR(recon, parms, powfs);
    setup_recon_GP(recon,parms,powfs);
    setup_recon_GA(recon,parms,powfs);
    /*assemble noise equiva angle inverse from powfs information */
    setup_recon_saneai(recon,parms,powfs);
    /*setup LGS tip/tilt/diff focus removal */
    setup_recon_TTFR(recon,parms,powfs);
    if(parms->sim.mffocus || parms->sim.ahstfocus){
	/*mvst uses information here*/
	setup_recon_focus(recon, powfs, parms);
    }
    switch(parms->recon.alg){
    case 0:
	setup_recon_mvr(recon, parms, powfs, aper);
	break;
    case 1:
	setup_recon_lsr(recon, parms, powfs, aper);
	if(parms->sim.dmproj){
	    setup_recon_mvr_fit(recon, parms, aper);
	}
	break;
    default:
	error("recon.alg=%d is not recognized\n", parms->recon.alg);
    }
 
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	locfree(powfs[ipowfs].gloc);
	dfree(powfs[ipowfs].gamp);
    }

    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(!parms->powfs[ipowfs].hasGS0 && powfs[ipowfs].GS0){
	    spcellfree(powfs[ipowfs].GS0);
	    powfs[ipowfs].GS0=NULL;
	}
    }
    mapfree(aper->ampground);
    toc2("setup_recon");
    return recon;
}

/*assemble matrix to do matrix vector multiply. Split from setup_recon because GPU may be used.*/
void setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    TIC;tic;
    if(parms->recon.mvm){
	if(parms->load.mvm){
	    recon->MVM=dcellread("%s", parms->load.mvm);
	    for(int idm=0; idm<parms->ndm; idm++){
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		    int ipowfs=parms->wfsr[iwfs].powfs;
		    if(!parms->powfs[ipowfs].skip){
			dmat *temp=recon->MVM->p[idm+iwfs*parms->ndm];
			if(temp->nx!=recon->aloc[idm]->nloc 
			   || temp->ny !=powfs[ipowfs].pts->nsa*2){
			    error("MVM[%d, %d] should have dimension (%ld, %ld), but have (%ld, %ld)\n", idm, iwfs, recon->aloc[idm]->nloc , powfs[ipowfs].pts->nsa*2, temp->nx, temp->ny);
			}
		    }
		}
	    }
	}
	if(parms->gpu.tomo && parms->gpu.fit){
#if USE_CUDA
	    gpu_setup_recon_mvm(parms, recon, powfs);
#endif
	}else if(!parms->load.mvm){
	    if(parms->load.mvmi){
		error("Not handled yet.\n");
	    }
	    if(parms->recon.alg==0){
		setup_recon_mvr_mvm(recon, parms, powfs);
	    }else{
		setup_recon_lsr_mvm(recon, parms, powfs);   
	    }
	}
	if(!parms->load.mvm && (parms->save.setup || parms->save.mvm)){
	    dcellwrite(recon->MVM, "MVM.bin");
	}
	if(parms->sim.mvmport){
	    mvm_client_send_m(parms, recon->MVM);
	}
	if(parms->gpu.tomo && parms->gpu.fit){
	    dcellfree(recon->MVM);
	}
    }
    toc2("setup_recon_mvm");
}
/**
   Free unused object in recon struct after preparation is done.
 */
void free_recon_unused(const PARMS_T *parms, RECON_T *recon){
    spcellfree(recon->GWR);
    if(parms->recon.alg!=0 || (parms->tomo.assemble && !parms->cn2.tomo)){
	/*We already assembled tomo matrix. don't need these matric any more. */
	dcellfree(recon->TTF);
	dcellfree(recon->PTTF);
    }
    /* Free arrays that will no longer be used after reconstruction setup is done. */
    spcellfree(recon->sanea); 
    spcellfree(recon->saneal);
    dfree(recon->neam); 
    if(!(parms->tomo.assemble && parms->tomo.alg==1) && !parms->cn2.tomo && !parms->tomo.bgs){
	/*We no longer need RL.M,U,V */
	spcellfree(recon->RL.M);
	dcellfree(recon->RL.U);
	dcellfree(recon->RL.V);
	spcellfree(recon->RR.M);
	dcellfree(recon->RR.U);
	dcellfree(recon->RR.V);
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

    if(parms->tomo.assemble){
	spcellfree(recon->GP);
	spcellfree(recon->GP2);
    }

    /*The following have been used in fit matrix. */
    if(parms->fit.assemble || parms->gpu.fit){
	dcellfree(recon->fitNW);
	spcellfree(recon->actslave);
    }
    /* when sim.dmproj=1, we need these matrices to use in FR.Mfun*/
    if(recon->FR.M && !parms->sim.dmproj && parms->fit.assemble && !parms->gpu.moao && !parms->sim.ncpa_calib){
	spfree(recon->W0); 
	dfree(recon->W1); 
	spcellfree(recon->HA); 
	spcellfree(recon->HXF); 
    }
    /*
      The following arrys are not used after preparation is done.
    */
    spcellfree(recon->GX);
    spcellfree(recon->GXtomo);/*we use HXWtomo instead. faster */
    spcellfree(recon->GXfocus);
    if(!(parms->cn2.tomo && parms->recon.split==2)){/*mvst needs GXlo when updating. */
	spcellfree(recon->GXlo);
    }
    if(parms->tomo.alg!=1 || (parms->tomo.assemble || (parms->tomo.square && !parms->dbg.tomo_hxw))){
	spcellfree(recon->HXWtomo);
    }
    if(parms->recon.mvm){
	muv_free(&recon->RR);
	muv_free(&recon->RL);
	muv_free(&recon->FR);
	muv_free(&recon->FL);
	muv_free(&recon->LR);
	muv_free(&recon->LL);
	fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
    }
}
/**
   Free the recon struct.
*/
void free_recon(const PARMS_T *parms, RECON_T *recon){
    CALL_ONCE;
    if(!recon) return;
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
    spcellfree(recon->GP2);
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
    dcellfree(recon->GFngs);
    dcellfree(recon->GFlgs);
    dcellfree(recon->GFall);
    dcellfree(recon->RFlgsg);
    dcellfree(recon->RFlgsa);
    dcellfree(recon->RFlgsx);
    dcellfree(recon->RFngsg);
    dcellfree(recon->RFngsa);
    dcellfree(recon->RFngsx);
    dcellfree(recon->RFdfx);
    dcellfree(recon->RFdfa);
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
   
    locarrfree(recon->xloc, npsr);
    maparrfree(recon->xmap, npsr);
    free(recon->xnx);
    free(recon->xny);
    free(recon->anx);
    free(recon->any);
    free(recon->anloc);
    free(recon->ngrad);
    locfree(recon->floc); 
    mapfree(recon->fmap);
    locfree(recon->ploc); 
    mapfree(recon->pmap);
    for(int idm=0; idm<ndm; idm++){
	if(recon->alocm[idm]!=recon->aloc[idm])
	    locfree(recon->alocm[idm]);
	locfree(recon->alocm[idm]);
	free(recon->aembed[idm]);
    }
    maparrfree(recon->amap, parms->ndm);
    free(recon->aembed);
    free(recon->alocm);
    free(recon->aloc);
    icellfree(recon->actstuck);
    icellfree(recon->actfloat);
    spcellfree(recon->actslave);
    dcellfree(recon->actcpl);
    dcellfree(recon->aimcc);/*used in filter.c */
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    muv_free(&recon->FR);
    muv_free(&recon->FL);
    muv_free(&recon->LR);
    muv_free(&recon->LL);
    dcellfree(recon->MVM);
    spcellfree(recon->saneai);
    fdpcg_free(recon->fdpcg); recon->fdpcg=NULL;
    cn2est_free(recon->cn2est);
    dcellfree(recon->opdxadd);
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
