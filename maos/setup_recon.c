/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
/**
   \file setup_recon.c
   Contains routines that setup the wavefront reconstructor and DM fitting.
*/
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
	loc_nxny(&recon->ploc_nx, &recon->ploc_ny, recon->ploc);
    }else{ 
	/*
	  Create a circular PLOC with telescope diameter by calling
	  create_metapupil with height of 0. We don't add any guard points.*/
	double dxr=parms->atmr.dx/parms->tomo.pos;//sampling of ploc
	map_t *pmap=create_metapupil_wrap 
	    (parms,0,dxr,0,0,0,T_PLOC,0,parms->fit.square);
	info2("PLOC is %ldx%ld, with sampling of %.2fm\n",pmap->nx,pmap->ny,dxr);
	recon->ploc=sqmap2loc(pmap);//convert map_t to loc_t
	recon->ploc_nx=pmap->nx;
	recon->ploc_ny=pmap->ny;
	sqmapfree(pmap);//free it.
	if(parms->save.setup){
	    locwrite(recon->ploc, "%s/ploc",dirsetup);
	}
    }
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
	if(!(exist("W0.bin.gz")&&exist("W1.bin.gz"))){
	    error("W0 or W1 not exist\n");
	}
	warning("Loading W0, W1");
	recon->W0=spread("W0.bin.gz");
	recon->W1=dread("W1.bin.gz");
    }else{
	/*
	  Compute W0,W1 weighting matrix that can be used to compute piston
	  removed wavefront variance of OPD: RMS=OPD'*(W0-W1*W1')*OPD; W0 is
	  sparse. W1 is a vector. These matrices are used for weighting in DM
	  fitting.
	*/
	if(parms->dbg.annular_W && parms->aper.din>0){
	    warning("Define the W0/W1 on annular aperture instead of circular.\n");
	    mkw_annular(recon->ploc, 0, 0, 
			parms->aper.din/2, parms->aper.d/2,
			&(recon->W0), &(recon->W1));
	    
	}else{
	    mkw_circular(recon->ploc,0,0,parms->aper.d/2.,
			  &(recon->W0), &(recon->W1));
	}
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
    if(parms->load.aloc){
	char *fn=parms->load.aloc;
	warning("Loading aloc from %s\n",fn);
	int naloc;
	recon->aloc=locarrread(&naloc,"%s",fn);
	if(naloc!=ndm) error("Invalid saved aloc");
	recon->aloc_nx=calloc(ndm, sizeof(long));
	recon->aloc_ny=calloc(ndm, sizeof(long));
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_nxny(&recon->aloc_nx[idm], &recon->aloc_ny[idm], recon->aloc[idm]);
	}
    }else{
	recon->aloc=calloc(ndm, sizeof(loc_t*));
	recon->aloc_nx=calloc(ndm, sizeof(long));
	recon->aloc_ny=calloc(ndm, sizeof(long));
	for(int idm=0; idm<parms->ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dx=parms->dm[idm].dx;
	    double offset=parms->dm[idm].offset+(parms->dm[idm].order%2)*0.5;
	    const double guard=parms->dm[idm].guard*parms->dm[idm].dx;

	    map_t *map=create_metapupil_wrap
		(parms,ht,dx,offset,guard,0,T_ALOC,0,parms->fit.square);
	
	    recon->aloc[idm]=sqmap2loc(map);
	    recon->aloc_nx[idm]=map->nx;
	    recon->aloc_ny[idm]=map->ny;
	    free(map->p);
	    free(map);
	    if(parms->plot.setup){
		plotloc("FoV", parms, recon->aloc[idm], ht, "aloc%d", idm);
	    }
	}//idm
        if(parms->save.setup){
	    locarrwrite(recon->aloc,parms->ndm,"%s/aloc",dirsetup);
	}
    }
    recon->aimcc=dcellnew(ndm,1);
    for(int idm=0; idm<parms->ndm; idm++){
	recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc[idm], NULL);
	dinvspd_inplace(recon->aimcc->p[idm]);
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
    if(parms->load.xloc){
	char *fn=parms->load.xloc;
	warning("Loading xloc from %s\n",fn);
	int nxloc;
	recon->xloc=locarrread(&nxloc,"%s",fn);
	if(nxloc!=npsr) 
	    error("Invalid saved file. npsr=%d, nxloc=%d\n",npsr,nxloc);
	recon->xloc_nx=calloc(npsr, sizeof(long));
	recon->xloc_ny=calloc(npsr, sizeof(long));
	for(int ips=0; ips<npsr; ips++){
	    loc_nxny(&recon->xloc_nx[ips], &recon->xloc_ny[ips], recon->xloc[ips]);
	}
    }else{
	recon->xloc=calloc(npsr, sizeof(loc_t *));
	recon->xloc_nx=calloc(npsr, sizeof(long));
	recon->xloc_ny=calloc(npsr, sizeof(long));

	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    const double dxr=recon->dx->p[ips];
	    const double guard=parms->tomo.guard*dxr;
	    long nin=0;
	    if(parms->tomo.precond==1 || parms->tomo.square==2){
		//FDPCG prefers power of 2 dimensions.
		nin=nextpow2((long)round(parms->aper.d/recon->dx->p[0]*2.))
		    *recon->os->p[ips]/recon->os->p[0];
		//nin=(long)pow(2,ceil(log2(parms->aper.d/recon->dx->p[0]*2)))
		//*recon->os->p[ips]/recon->os->p[0];
		warning("layer %d xloc is set to %ld for FDPCG\n",ips,nin);
	    }
	    map_t *map=create_metapupil_wrap
		(parms,ht,dxr,0,guard,nin,T_XLOC,0,parms->tomo.square);
	    info2("layer %d: xloc map is %ldx%ld, "
		 "with sampling of %.2f m\n",ips,
		 map->nx,map->ny,dxr);
	    recon->xloc[ips]=sqmap2loc(map);
	    recon->xloc_nx[ips]=map->nx;
	    recon->xloc_ny[ips]=map->ny;
	    free(map->p);free(map);
	 
	    if(parms->plot.setup){
		plotloc("FoV",parms,recon->xloc[ips],ht, "xloc%d",ips);
	    }
	}
	if(parms->save.setup){
	    locarrwrite(recon->xloc, recon->npsr, "%s/xloc",dirsetup);
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
    const int nwfs=parms->nwfs;
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
	    int ipowfs = parms->wfs[iwfs].powfs;
	    if(parms->tomo.split!=2 && parms->powfs[ipowfs].skip){
		//don't need HXW for low order wfs that does not participate in tomography.
		continue;
	    }
	    double  hs = parms->powfs[ipowfs].hs;
	    for(int ips=0; ips<npsr; ips++){
		double  ht = recon->ht->p[ips];
		double  scale=1. - ht/hs;
		double  displace[2];
		displace[0]=parms->wfs[iwfs].thetax*ht;
		displace[1]=parms->wfs[iwfs].thetay*ht;
		HXW[ips][iwfs]=mkh(recon->xloc[ips], ploc, NULL, 
				   displace[0],displace[1],scale,
				   parms->tomo.cubic, parms->tomo.iac);
		
	    }
	}
	toc2("done");
	if(parms->save.setup){
	    spcellwrite(recon->HXW, "%s/HXW",dirsetup);
	}
    }
    spcellfree(recon->HXWtomo);
    recon->HXWtomo=spcellnew(recon->HXW->nx, recon->HXW->ny);
    PDSPCELL(recon->HXWtomo,HXWtomo);
    PDSPCELL(recon->HXW,HXW);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){//for tomography
	    for(int ips=0; ips<npsr; ips++){
		HXWtomo[ips][iwfs]=spref(HXW[ips][iwfs]);
	    }
	} 
     }
}
/**
   Setup gradient operator from ploc to wavefront sensors.
 */
static void
setup_recon_GP(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfs;
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
	    //use ploc as an intermediate plane.
	    switch(parms->powfs[ipowfs].gtype_recon){
	    case 0:{ /*Create averaging gradient operator (gtilt) from PLOC,
		       using fine sampled powfs.loc as intermediate plane*/
		double displace[2]={0,0};
		info2(" Gploc");
		GP->p[ipowfs]=mkg(ploc,powfs[ipowfs].loc,powfs[ipowfs].amp,
				  powfs[ipowfs].saloc,1,1,displace,1);
		if(parms->save.setup){
		    spwrite(GP->p[ipowfs], "%s/powfs%d_GP", dirsetup, ipowfs);
		}
	    }
		break;
	    case 1:{ /*Create ztilt operator from PLOC, using fine sampled
		       powfs.loc as intermediate plane*/
		info2(" Zploc");
		dsp *H=mkh(ploc,powfs[ipowfs].loc,powfs[ipowfs].amp, 0,0,1,0,0);
		GP->p[ipowfs]=spmulsp(powfs[ipowfs].ZS0,H);
		spfree(H);
	    }
		break;
	    default:
		error("Invalid gtype_recon\n");
	    }
	}
	toc2(" done");
	if(parms->save.setup){
	    spcellwrite(GP,"%s/GP",dirsetup);
	}
    }
    //assign GP for powfs to recon->GP for each wfs
    spcellfree(recon->GP);
    recon->GP=spcellnew(nwfs,1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfs[iwfs].powfs;
	recon->GP->p[iwfs]=spref(GP->p[ipowfs]);
    }
    spcellfree(GP);//assigned to recon->GP already;
}
/**
   Setup gradient operator form aloc for wfs by using GP and HA
*/
static void
setup_recon_GA(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfs;
    const int ndm=parms->ndm;
    spcellfree(recon->GA);
    if(parms->load.GA){
	warning2("Loading saved GA\n");
	recon->GA=spcellread("GA.bin.gz");
	if(recon->GA->nx!=nwfs || recon->GA->ny!=ndm)
	    error("Wrong saved GA\n");
	PDSPCELL(recon->GA,GA);
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc[idm]->nloc;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs = parms->wfs[iwfs].powfs;
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
	    int ipowfs = parms->wfs[iwfs].powfs;
	    if(parms->sim.skysim && parms->powfs[ipowfs].lo){
		continue;
	    }
	    double  hs = parms->powfs[ipowfs].hs;
	    for(int idm=0; idm<ndm; idm++){
		double  ht = parms->dm[idm].ht;
		double  scale=1. - ht/hs;
		double  displace[2];
		displace[0]=parms->wfs[iwfs].thetax*ht;
		displace[1]=parms->wfs[iwfs].thetay*ht;
		dsp *H=mkh(recon->aloc[idm], ploc, NULL, 
			   displace[0],displace[1],scale,
			   parms->dm[idm].cubic,parms->dm[idm].iac);
		
		GA[idm][iwfs]=spmulsp(recon->GP->p[iwfs],H);
		spfree(H);
	    }//idm
	}
	if(parms->save.setup){
	    spcellwrite(recon->GA, "%s/GA",dirsetup);
	}
    	toc2("done");
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
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].lo){//for low order wfs
		GAlo[idm][iwfs]=spref(GA[idm][iwfs]);
	    }else{
		GAhi[idm][iwfs]=spref(GA[idm][iwfs]);		
	    }
	}
    }
}
/**
   Crate the xloc to wfs gradient operator
*/
static void 
setup_recon_GX(RECON_T *recon, const PARMS_T *parms){
    const int nwfs=parms->nwfs;
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
    toc2("done");
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
	int ipowfs=parms->wfs[iwfs].powfs;
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
	if(!parms->dbg.fitonly && parms->sim.mffocus && parms->sim.closeloop){
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
    const int nwfs=parms->nwfs;
    if(recon->saneai){
	spcellfree(recon->saneai);
    }
    spcell *saneai=recon->saneai=spcellnew(nwfs,nwfs);
    info2("saneai:");
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int nsa=powfs[ipowfs].pts->nsa;
	if(parms->powfs[ipowfs].neareconfile){//taks precedance
	    dmat *nea=dread("%s_wfs%d",parms->powfs[ipowfs].neareconfile,iwfs);//rad
	    if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy){
		//sanity check on provided nea.
		const int nmtch=powfs[ipowfs].intstat->mtche->ny;
		int indsanea=0;
		if(nmtch==1){
		    indsanea=0;
		}else if(nmtch==parms->powfs[ipowfs].nwfs){
		    indsanea=parms->powfs[ipowfs].indwfs[iwfs];
		}else{
		    error("invalid\n");
		}
		PDCELL(powfs[ipowfs].intstat->saneaixy, saneaixy);
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
	 
	    dcwpow(nea,-2);//rad^-2
	    saneai->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,nea->p,1.);
	    dfree(nea);
	}else if((parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].neaphy) && 
		 !parms->powfs[ipowfs].phyusenea){
	    //Physical optics
	    if(parms->powfs[ipowfs].phytype==1){
		saneai->p[iwfs+iwfs*nwfs] =spnew(nsa*2,nsa*2,nsa*4);
		spint *pp=saneai->p[iwfs+iwfs*nwfs]->p;
		spint *pi=saneai->p[iwfs+iwfs*nwfs]->i;
		double *px=saneai->p[iwfs+iwfs*nwfs]->x;
	
		long count=0;
		const int nmtch=powfs[ipowfs].intstat->mtche->ny;
		int indsanea=0;
		if(nmtch==1){
		    indsanea=0;
		}else if(nmtch==parms->powfs[ipowfs].nwfs){
		    indsanea=parms->powfs[ipowfs].indwfs[iwfs];
		}else{
		    error("invalid\n");
		}

		PDCELL(powfs[ipowfs].intstat->saneaixy, saneaixy);
		for(int isa=0; isa<nsa; isa++){
		    pp[isa]=count;
		    pi[count]=isa;//xx
		    px[count]=saneaixy[indsanea][isa]->p[0];
		    count++;
		    pi[count]=isa+nsa;//yx
		    px[count]=saneaixy[indsanea][isa]->p[1];
		    count++;
		}
		for(int isa=0; isa<nsa; isa++){
		    pp[isa+nsa]=count;
		    pi[count]=isa;//xy
		    px[count]=saneaixy[indsanea][isa]->p[2];
		    count++;
		    pi[count]=isa+nsa;//yy
		    px[count]=saneaixy[indsanea][isa]->p[3];
		    count++;
		}
		pp[nsa*2]=count;
	
	    }else{
		error("Not implemented yet\n");
	    }
	}else{
	    const double nea=parms->powfs[ipowfs].nearecon/206265000.;
	    if(nea<1.e-15) error("nea is too small\n");
	    //nea scales as sqrt(1/dtrat) so neaisq scales as dtrat.
	    const double neaisq=pow(nea,-2)*parms->powfs[ipowfs].dtrat;
	    dmat *neai=dnew(nsa,2);
	    //scale neaisq by area. (seeing limited)
	    //scale neaisq by area^2 if diffraction limited
	    //only implementing seeing limited here.
	    double (*neaip)[nsa]=(double(*)[nsa])neai->p;
	    double *area=powfs[ipowfs].pts->area;
	    for(int i=0; i<nsa; i++){
		neaip[0][i]=neaip[1][i]=neaisq*(area[i]);
	    }
	    saneai->p[iwfs+iwfs*nwfs]=spnewdiag(nsa*2,neai->p,1.);
	    dfree(neai);
	}
    }//iwfs
    
    //Compute the averaged SANEA for each WFS
    recon->neam=dnew(parms->nwfs, 1);
    double neamhi=0; 
    int counthi=0;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int nsa=powfs[ipowfs].pts->nsa;
	dmat *sanea=spdiag(recon->saneai->p[iwfs+iwfs*parms->nwfs]);
	double area_thres;
	if(nsa>4){
	    area_thres=0.9;
	}else{
	    area_thres=0;
	}
	double nea2_sum=0;
	int count=0;
	for(int isa=0; isa<nsa; isa++){
	    if(powfs[ipowfs].pts->area[isa]>area_thres){
		nea2_sum+=1./(sanea->p[isa])+1./(sanea->p[isa+nsa]);
		count++;
	    }
	}
	dfree(sanea);
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
	spcellwrite(recon->saneai,"%s/saneai",dirsetup);
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
    int nwfs=parms->nwfs;
    recon->TT=dcellnew(nwfs,nwfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].trs){
	    info2("powfs %d has tip/tilt removed in tomography\n", ipowfs);
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
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs[jwfs];
		if(parms->powfs[ipowfs].skip){
		    error("This WFS %d should be included in Tomo.\n", iwfs);
		}
		dcp(&recon->TT->p[iwfs+iwfs*nwfs], TT);
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
    }
    int nwfs=parms->nwfs;
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
    if(parms->save.setup){
	dcellwrite(recon->DF, "%s/DF",dirsetup);
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
    recon->PTTF=dcellpinv(recon->TTF, NULL,recon->saneai);
    dcellfree(recon->DF);
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
		    if(ips==0){
			info2("Laplacian reg. is using perodic bc\n");
		    }
		    recon->L2->p[ips+npsr*ips]=mklaplacian_map
			(recon->xloc_nx[ips], recon->xloc_ny[ips],
			 recon->xloc[ips]->dx, recon->r0,
			 recon->wt->p[ips]);
		}else{//reflecive bc
		    if(ips==0){
			info2("Laplacian reg. is using reflective bc\n");
		    }
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
		long nx=recon->xloc_nx[ips];
		long ny=recon->xloc_ny[ips];
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
	    fftxopd->p[ips]=cnew(recon->xloc_nx[ips], recon->xloc_ny[ips]);
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
	    int nn=nextpow2(MAX(recon->xloc_nx[ips], recon->xloc_ny[ips]))+1;
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
void setup_recon_tomo_matrix(RECON_T *recon, const PARMS_T *parms){
    //if not cg or forced, build explicitly the tomography matrix.
    int npsr=recon->npsr;
    int nwfs=parms->nwfs;
    //Free OLD matrices if any.
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    info2("Before assembling tomo matrix:\t%.2f MiB\n",get_job_mem()/1024.);

    if(parms->load.tomo){
	//right hand side.
	warning("Loading saved recon->RR\n");
	recon->RR.M=spcellread("RRM.bin.gz");
	if(recon->has_ttr){
	    recon->RR.U=dcellread("RRU.bin.gz");
	    recon->RR.V=dcellread("RRV.bin.gz");
	}
	//Left hand side
	warning("Loading saved recon->RL\n");
	if(exist("RLM.bin"))
	    recon->RL.M=spcellread("RLM.bin");
	else
	    recon->RL.M=spcellread("RLM.bin.gz");
	recon->RL.U=dcellread("RLU.bin.gz");
	recon->RL.V=dcellread("RLV.bin.gz");
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
	}else{
	    /*Apply tikholnov regularization.*/
	    if(fabs(parms->tomo.tikcr)>1.e-15){
		//Estimated from the Formula
		double maxeig=pow(recon->neamhi * recon->xloc[0]->dx, -2);
		double tikcr=parms->tomo.tikcr;
		info2("Adding tikhonov constraint of %g to RLM\n",tikcr);
		info2("The maximum eigen value is estimated to be around %g\n", maxeig);
		
		spcelladdI(recon->RL.M, tikcr*maxeig);
	    }
	}
	//add L2 and ZZT
	switch(parms->tomo.cxx){
	case 0://Add L2'*L2 to RL.M
	     for(int ips=0; ips<npsr; ips++){
		dsp* tmp=sptmulsp(recon->L2->p[ips+npsr*ips], 
				  recon->L2->p[ips+npsr*ips]);
		if(!tmp){
		    error("L2 is empty!!\n");
		}
		spadd(&RLM[ips][ips], tmp);
		spfree(tmp);
	    }
	    break;
	case 1://Need to apply invpsd separately
	    recon->RL.extra = recon->invpsd;
	    recon->RL.exfun = (CGFUN)apply_invpsd;
	    break;
	case 2://Need to apply fractal separately
	    recon->RL.extra = recon->fractal;
	    recon->RL.exfun = (CGFUN)apply_fractal;
	}

	//Symmetricize, remove values below 1e-15*max and sort RLM (optional).
	//spcellsym(recon->RL.M);

	//Low rank terms for low order wfs. Only in Integrated tomography.
	dcell *ULo=dcellnew(npsr,nwfs);
	PDCELL(ULo, pULo);
	dcell *VLo=dcellnew(npsr,nwfs);
	PDCELL(VLo, pVLo);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
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
	recon->RL.U=dcellcat(recon->RR.U, ULo, 2);
	dcell *GPTTDF=NULL;
	sptcellmulmat(&GPTTDF, recon->GX, recon->RR.V, 1);
	recon->RL.V=dcellcat(GPTTDF, VLo, 2);
	dcellfree(GPTTDF);
	dcellfree(ULo);
	dcellfree(VLo);
	//Remove empty cells.
	dcelldropempty(&recon->RR.U,2);
	dcelldropempty(&recon->RR.V,2);
	dcelldropempty(&recon->RL.U,2);
	dcelldropempty(&recon->RL.V,2);

	long nll=0,nlr=0;
	if(recon->RL.U){
	    /*
	      balance UV. may not be necessary. Just to compare well against
	      laos.
	    */
	    double r0=recon->r0;
	    double dx=recon->xloc[0]->dx;
	    double val=laplacian_coef(r0,1,dx);//needs to be a constant
	    dcellscale(recon->RL.U, 1./val);
	    dcellscale(recon->RL.V, val);
	    //collect statistics.
	    PDCELL(recon->RR.U,RRU);
	    PDCELL(recon->RL.U,RLU);
	    for(int i=0; i<recon->RR.U->ny;i++){
		if(RRU[i][0]) nlr+=RRU[i][0]->ny;
	    }
	    for(int i=0; i<recon->RL.U->ny;i++){
		if(RLU[i][0]) nll+=RLU[i][0]->ny;
	    }
	}
	info2("Tomography # of Low rank terms: %ld in RHS, %ld in LHS\n", 
	     nlr,nll);
	if(parms->save.setup && parms->save.recon){
	    spcellwrite(recon->RR.M,"%s/RRM",dirsetup);
	    dcellwrite(recon->RR.U,"%s/RRU",dirsetup);
	    dcellwrite(recon->RR.V,"%s/RRV",dirsetup);

	    spcellwrite(recon->RL.M,"%s/RLM.bin",dirsetup);//disable compression
	    dcellwrite(recon->RL.U,"%s/RLU",dirsetup);
	    dcellwrite(recon->RL.V,"%s/RLV",dirsetup); 
	}
	spcellfree(GXtomoT);
    }
    info2("After assemble matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    if(parms->tomo.alg==0 || parms->tomo.split==2){
	//We need cholesky decomposition in CBS or MVST method.
	muv_chol_prep(&(recon->RL));
	if(parms->save.setup && parms->save.recon){
	    chol_save(recon->RL.C,"%s/RLC",dirsetup);
	}
	info2("After cholesky on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    if(parms->tomo.assemble && !parms->cn2.tomo){
	//Don't free PTT. Used in forming LGS uplink err
	if(parms->tomo.piston_cr)
	    spcellfree(recon->ZZT);
    }
    if((!parms->tomo.assemble || parms->tomo.alg==0) && !parms->cn2.tomo){
	//We just need cholesky factors.
	info2("Freeing RL.M,U,V\n");
	spcellfree(recon->RL.M);
	dcellfree(recon->RL.U);
	dcellfree(recon->RL.V);
    }
    if(parms->tomo.assemble){
	spcellfree(recon->GP);
	spcellfree(recon->HXWtomo);
    }
    info2("At exit :\t%.2f MiB\n",get_job_mem()/1024.);
}
/**
   Update assembled tomography matrix with new L2.
 */
void setup_recon_tomo_matrix_update(RECON_T *recon, const PARMS_T *parms){
    if(parms->tomo.alg==1&&!parms->tomo.assemble){//no need to do anything
	return;
    }
    switch(parms->tomo.cxx){
    case 0:{ //Need to adjust RLM with the new L2.
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
	break;
    case 1://Need to apply invpsd separately
	recon->RL.extra = recon->invpsd;
	recon->RL.exfun = (CGFUN)apply_invpsd;
	break;
    case 2://Need to apply fractal separately
	recon->RL.extra = recon->fractal;
	recon->RL.exfun = (CGFUN)apply_fractal;
    }
}
/**
   Setup ray tracing operator HXF,HA from xloc, aloc to aperture ploc */
static void
setup_recon_HXFHA(RECON_T *recon, const PARMS_T *parms){
    //Greate HXF and HA matrix;
    if(parms->load.HXF && exist(parms->load.HXF)){
	warning("Loading saved HXF\n");
	recon->HXF=spcellread("%s",parms->load.HXF);
    }else{
	info2("Generating HXF");TIC;tic;
	const int nfit=parms->fit.nfit;
	const int npsr=recon->npsr;
	recon->HXF=spcellnew(nfit, npsr);
	PDSPCELL(recon->HXF,HXF);
	for(int ifit=0; ifit<nfit; ifit++){
	    for(int ips=0; ips<npsr; ips++){
		const double scale=1.;
		const double ht = recon->ht->p[ips];
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		HXF[ips][ifit]=mkh(recon->xloc[ips], recon->ploc, NULL,
				   displace[0], displace[1], scale,
				   parms->tomo.cubic, parms->tomo.iac);
	    }
	}
	if(parms->save.setup){
	    spcellwrite(recon->HXF,"%s/HXF.bin.gz",dirsetup);
	}
	toc2("done");
    }
    if(parms->load.HA && exist(parms->load.HA)){
	warning("Loading saved HA\n");
	recon->HA=spcellread("%s",parms->load.HA);
    }else{
	const int nfit=parms->fit.nfit;
	const int ndm=parms->ndm;
	recon->HA=spcellnew(nfit, ndm);
	PDSPCELL(recon->HA,HA);
	info2("Generating HA");TIC;tic;
	for(int ifit=0; ifit<nfit; ifit++){
	    for(int idm=0; idm<ndm; idm++){
		const double scale=1.;
		const double ht=parms->dm[idm].ht;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		HA[idm][ifit]=mkh(recon->aloc[idm], recon->ploc, NULL,
				  displace[0], displace[1], 
				  scale,parms->dm[idm].cubic,parms->dm[idm].iac);
	    }
	}
	toc2("done");
	if(parms->save.setup){
	    spcellwrite(recon->HA,"%s/HA",dirsetup);
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
    recon->NW=dcellnew(ndm,1);
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
	recon->NW->p[idm]=dnew(nloc, nnw);
    }
    int inw=0;//current column
    if(parms->fit.lrt_piston){
	info2("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc[idm]->nloc;
	    double *p=recon->NW->p[idm]->p+(inw+idm)*nloc;
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
		double *p=recon->NW->p[idm]->p+(inw+(idm-1)*2)*nloc;
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
		double *p=recon->NW->p[idm]->p+inw*nloc;
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
	dcellwrite(recon->NW,"%s/NW.bin.gz",dirsetup);
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
    dsp *(*HAT)[ndm]=(dsp*(*)[ndm])HATc->p;
    dsp*(*HXF)[nfit]=(dsp*(*)[nfit])recon->HXF->p;
    dsp*(*HA)[nfit]=(dsp*(*)[nfit])recon->HA->p;
    info2("Before assembling fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    //Assemble Fit matrix.
    int npsr=recon->npsr;
    if(parms->load.fit){
	if(!(exist("FRM.bin.gz") && 
	     exist("FRU.bin.gz") && exist("FRV.bin.gz"))){
	    error("FRM, FRU, FRV (.bin.gz) not all exist\n");
	}
	warning("Loading saved recon->FR\n");
	recon->FR.M=spcellread("FRM.bin.gz");
	recon->FR.U=dcellread("FRU.bin.gz");
	recon->FR.V=dcellread("FRV.bin.gz");
    }else{
	info2("Building recon->FR\n");
	recon->FR.M=spcellnew(ndm, npsr);
	dsp*(*FRM)[ndm]=(dsp *(*)[ndm])recon->FR.M->p;

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
	recon->FR.U=dcellnew(ndm, 1);
	dmat **FRU=recon->FR.U->p;

	recon->FR.V=dcellnew(npsr, 1);
	dmat **FRV=recon->FR.V->p;  
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

	
	if(parms->save.setup && parms->save.recon){
	    spcellwrite(recon->FR.M,"%s/FRM.bin.gz",dirsetup);
	    dcellwrite(recon->FR.U,"%s/FRU.bin.gz",dirsetup);
	    dcellwrite(recon->FR.V,"%s/FRV.bin.gz",dirsetup);
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
    //info2("recon->fitscl=%g\n", recon->fitscl);
    if(parms->load.fit){
	if(!(exist("FLM.bin.gz") && 
	     exist("FLU.bin.gz") && exist("FLV.bin.gz"))){
	    error("FLM, FLU, FLV (.bin.gz) not all exist\n");
	}
	warning("Loading saved recon->FL\n");
	recon->FL.M=spcellread("FLM.bin.gz");
	recon->FL.U=dcellread("FLU.bin.gz");
	recon->FL.V=dcellread("FLV.bin.gz");
    }else{
	//Depends on FRU; Don't move forward
	fit_prep_lrt(recon,parms);
	if(parms->fit.actslave){
	    recon->actslave=act_slaving(recon->aloc, recon->HA, recon->W1, recon->NW,0.1, 1./recon->ploc->nloc);
	    if(parms->save.setup){
		spcellwrite(recon->actslave,"%s/actslave.bin.gz",dirsetup);
		dcellwrite(recon->NW,"%s/NW2.bin.gz",dirsetup);
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
	    recon->FL.U=dcellcat_each(recon->FR.U, recon->NW, 2);
	    dcell *tmp=NULL;//negative NW.
	    dcelladd(&tmp, 1, recon->NW, -1);
	    recon->FL.V=dcellcat_each(recon->FR.U, tmp, 2);
	    dcellfree(tmp);
	}
	if(recon->actslave){
	    spcelladd(&recon->FL.M, recon->actslave);
	}
	//spcellsym(recon->FL.M);
	info2("DM Fit # of Low rank terms: %ld in RHS, %ld in LHS\n",
	     recon->FR.U->p[0]->ny, recon->FL.U->p[0]->ny);
	if(parms->save.setup && parms->save.recon){
	    spcellwrite(recon->FL.M,"%s/FLM.bin.gz",dirsetup);
	    dcellwrite(recon->FL.U,"%s/FLU.bin.gz",dirsetup);
	    dcellwrite(recon->FL.V,"%s/FLV.bin.gz",dirsetup);
	}
    }
    info2("After assemble matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    if(parms->fit.alg==0  || parms->tomo.split==2){
	if(fabs(parms->fit.tikcr)<1.e-14){
	    warning("tickcr=%g is too small or not applied\n", 
		    parms->fit.tikcr);
	}
	muv_chol_prep(&(recon->FL));
	if(parms->save.setup && parms->save.recon){
	    chol_save(recon->FL.C,"%s/FLC",dirsetup);
	}
	info2("After cholesky on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    if(parms->fit.alg==0){
	info2("Freeing FL.M,U,V\n");
	spcellfree(recon->FL.M);
	dcellfree(recon->FL.U);
	dcellfree(recon->FL.V);
    }

    info2("At exit :\t%.2f MiB\n",get_job_mem()/1024.);
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
	if(parms->powfs[ipowfs].hasllt){
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
    dcell *Gfocus=dcellnew(parms->nwfs, 1);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	dmat *opd=dnew(powfs[ipowfs].loc->nloc,1);
	loc_add_focus(opd->p, powfs[ipowfs].loc, 1);
	const int nsa=powfs[ipowfs].pts->nsa;
	Gfocus->p[iwfs]=dnew(nsa*2,1);
	if(parms->powfs[ipowfs].gtype_recon==1){
	    pts_ztilt(Gfocus->p[iwfs]->p, powfs[ipowfs].pts,
		      powfs[ipowfs].imcc, powfs[ipowfs].amp, opd->p);
	}else{
	    spmulmat(&Gfocus->p[iwfs], powfs[ipowfs].GS0, opd, 1);
	}
	dfree(opd);
    }
    dmat *GMGngs=NULL;
    dcell *GMngs=dcellnew(1, parms->nwfs);
    /*Compute focus constructor from NGS Grads. fuse grads
      together to construct a single focus measurement*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].trs==0 && parms->powfs[ipowfs].order>1){
	    info2("wfs %d will be used to track focus\n", iwfs);
	}else{
	    continue;
	}
	spmulmat(&GMngs->p[iwfs], recon->saneai->p[iwfs+parms->nwfs*iwfs], 
		 Gfocus->p[iwfs],1);
	dmm(&GMGngs,Gfocus->p[iwfs], GMngs->p[iwfs], "tn",1);
    }
    dinvspd_inplace(GMGngs);
    dcell *RFngs=recon->RFngs=dcellnew(1, parms->nwfs);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	if(!Gfocus->p[iwfs]) continue;
	dmm(&RFngs->p[iwfs],GMGngs, GMngs->p[iwfs],"nt",1);
    }
    dfree(GMGngs);
    dcellfree(GMngs);
    /*
      Compute focus constructor from LGS grads. A
      constructor for each LGS.
     */
    dcell *RFlgs=recon->RFlgs=dcellnew(parms->nwfs, 1);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].hasllt)
	    continue;
	dmat *GMtmp=NULL;
	dmat *GMGtmp=NULL;
	spmulmat(&GMtmp, recon->saneai->p[iwfs+parms->nwfs*iwfs], 
		 Gfocus->p[iwfs], 1);
	dmm(&GMGtmp, Gfocus->p[iwfs], GMtmp, "tn",1);
	dinvspd_inplace(GMGtmp);
	dmm(&RFlgs->p[iwfs], GMGtmp, GMtmp, "nt", 1);
	dfree(GMtmp);
	dfree(GMGtmp);
    }
  
    if(parms->save.setup){
	dcellwrite(Gfocus,"%s/Gfocus.bin.gz",dirsetup);
	dcellwrite(RFngs,"%s/RFngs.bin.gz",dirsetup);
	dcellwrite(RFlgs,"%s/RFlgs.bin.gz",dirsetup);
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
      called muv_chol_solve. The former doesn ot apply the low rank
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
	dcell *GXLT=dcelltrans(recon->GXL);
	muv_chol_solve_cell(&U, &recon->RL, GXLT);
	dcellfree(GXLT);
	dcell *rhs=NULL;
	muv(&rhs, &recon->FR, U, 1);
	muv_chol_solve_cell(&FU, &recon->FL, rhs);
	dcellfree(rhs);
	if(parms->save.mvst || parms->save.setup){
	    dcellwrite(U, "%s/mvst_U", dirsetup);
	    dcellwrite(FU, "%s/mvst_FU", dirsetup);
	}
	if(parms->tomo.alg!=0){
	    muv_chol_free(&recon->RL);
	}
	if(parms->fit.alg!=0){
	    muv_chol_free(&recon->FL);
	}
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
   
    {
	/*
	  Change FUw*Minv -> FUw*(U*sigma^-1/2) * (U*sigma^1/2)'*Minv
	  columes of FUw*(U*sigma^-1/2) are the eigen vectors.
	 */
	dcell *Q=NULL;//the NGS modes in ploc.
	spcellmulmat(&Q, recon->HA, FUw, 1);
	dcell *QwQc=calcWmcc(Q,Q,recon->W0,recon->W1,recon->fitwt);
	dmat *QwQ=dcell2m(QwQc);
	dmat *QSdiag=NULL, *QU=NULL, *QVt=NULL;
	dsvd(&QSdiag, &QU, &QVt, QwQ);
	if(parms->save.setup) dwrite(QSdiag,"%s/mvst_QSdiag",dirsetup);
	dcwpow(QSdiag, -1./2.);
	dmuldiag(QU,QSdiag);//U*sigma^-1/2
	d2cell(&QwQc,QU,NULL);
	dcell *FUw_keep=FUw;FUw=NULL; 
	dcellmm(&FUw, FUw_keep, QwQc, "nn", 1);
	dcellfree(FUw_keep);
	dcwpow(QSdiag,-2);
	dmuldiag(QU,QSdiag);//U*sigma^1/2
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
    //Make MVRngs 1xn cell instead of nxn cell.
    recon->MVRngs=dcellreduce(Minv,1);
    recon->MVModes=dcellreduce(FUw,2);

  
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
    {
	dcell *Q=NULL;
	spcellmulmat(&Q, recon->HA, recon->MVModes,1);
	dcell *MCC=calcWmcc(Q,Q,recon->W0,recon->W1,recon->fitwt);
	if(parms->save.setup){
	    dcellwrite(MCC,"%s/mvst_MCC",dirsetup);
	}
	dcellfree(MCC);
	dcellfree(Q);
    }
}
/**
   Setup either the minimum variance reconstructor by calling setup_recon_mvr()
or least square reconstructor by calling setup_recon_lsr() */
RECON_T *setup_recon(const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    RECON_T * recon = calloc(1, sizeof(RECON_T));
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
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs<=1) continue;
	if(parms->powfs[ipowfs].dfrs){
	    recon->has_dfr=1;
	    break;
	}
    }   
    //setup DM actuator grid
    setup_recon_aloc(recon,parms);
    //setup pupil coarse grid
    setup_recon_ploc(recon,parms);
    setup_recon_GP(recon,parms,powfs);
    setup_recon_GA(recon,parms,powfs);
    //assemble noise equiva angle inverse from powfs information
    setup_recon_saneai(recon,parms,powfs);
    //setup LGS tip/tilt/diff focus removal
    setup_recon_TTFR(recon,parms,powfs);

    if(parms->sim.recon==0){
	setup_recon_mvr(recon, parms, powfs, aper);
    }else{
	setup_recon_lsr(recon, parms, powfs, aper);
    }
    if(parms->sim.recon!=0 || (parms->tomo.assemble && !parms->cn2.tomo)){
	//We already assembled tomo matrix. don't need these matric any more.
	dcellfree(recon->TTF);
	dcellfree(recon->PTTF);
    }
    return recon;
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
    toc("loc done");
    setup_recon_HXW(recon,parms);
    setup_recon_GX(recon,parms);
    spcellfree(recon->HXW);//only keep HXWtomo for tomography
    //setup inverse noise covariance matrix.
    toc("GX GA");
    //prepare for tomography setup
    setup_recon_tomo_prep(recon,parms);
    toc("tomo_prep");
    setup_recon_HXFHA(recon,parms);
    toc("HXFHA");
    if(parms->tomo.assemble || parms->tomo.split==2){
	/*assemble the matrix only if not using CG CG apply the
	  individual matrices on fly to speed up and save memory. */
	setup_recon_tomo_matrix(recon,parms);
    }
    //always assemble fit matrix, faster if many directions
    setup_recon_fit_matrix(recon,parms);
    toc("fit_matrix");
    //moao
    setup_recon_moao(recon,parms);
    toc("moao");
    if(parms->sim.mffocus){
	setup_recon_focus(recon, powfs, parms);
    }
    if(parms->tomo.split){
	//split tomography
	if(parms->tomo.split && parms->ndm<=2){
	    //setup the ngsmode in both ahst and mvst mode 
	    setup_ngsmod(parms,recon,aper,powfs);
	}else{
	    error("Not implemented");
	}
	if(parms->tomo.split==2){
	    setup_recon_mvst(recon,parms);
	}
    }
    if(parms->tomo.precond==1){
	toc2("before fdpcg_prepare");
	recon->fdpcg=fdpcg_prepare(parms, recon, powfs);
	toc2("fdpcg_prepare");
    }
    //The following have been used in fit matrix.
    dcellfree(recon->NW);
    spcellfree(recon->actslave);
    spfree(recon->W0); 
    dfree(recon->W1); 
    spcellfree(recon->HA); 
 
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(powfs[ipowfs].ZS0){
	    //ZS0 is used in setup_recon. not in simulation, too large a matrix.
	    spfree(powfs[ipowfs].ZS0);
	    powfs[ipowfs].ZS0=NULL;
	}
	if(!parms->powfs[ipowfs].hasGS0 && powfs[ipowfs].GS0){
	    spfree(powfs[ipowfs].GS0);
	    powfs[ipowfs].GS0=NULL;
	}
    }
    //Free memory that will not be used anymore in simulation
    if(recon->ngsmod && recon->ngsmod->Mdm) {
	dcellfree(recon->ngsmod->Mdm);
    }
    /*
       The following arrys are not used after preparation is done.
     */
    sqmapfree(aper->ampground);
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

   This is very close to RR except replacing GX with GA
 */
void setup_recon_lsr(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs, APER_T *aper){
    spcell *GAlsr;
    if(parms->lsr.split){
	//high order wfs only in split mode.
	GAlsr=recon->GAhi;
	setup_ngsmod(parms,recon,aper,powfs);
    }else{
	//all wfs in integrated mode.
	GAlsr=recon->GA;
    }
    const spcell *saneai=recon->saneai;
    spcell *GAlsrT=spcelltrans(GAlsr);
    info2("Building recon->LR\n");
    recon->LR.M=spcellmulspcell(GAlsrT, saneai, 1);
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
    //actuator slaving. important
    spcell *actslave=act_slaving(recon->aloc, recon->GAhi, NULL, NULL,0.5,
				 sqrt(maxeig));
    spcellwrite(actslave,"actslave");
    spcelladd(&recon->LL.M, actslave);
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfs;
    //Low rank terms for low order wfs. Only in Integrated tomography.
    dcell *ULo=dcellnew(ndm,nwfs);
    PDCELL(ULo, pULo);
    dcell *VLo=dcellnew(ndm,nwfs);
    PDCELL(VLo, pVLo);
    PDSPCELL(recon->LR.M, LRM);
    PDSPCELL(recon->GA, GA);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].skip){
	    continue;
	}
	if(parms->powfs[ipowfs].lo){
	    for(int idm=0; idm<ndm; idm++){
		spfull(&pULo[iwfs][idm], LRM[iwfs][idm],-1);
		sptfull(&pVLo[iwfs][idm], GA[idm][iwfs],1);
	    }
	}
    }
    recon->LL.U=dcellcat(recon->LR.U, ULo, 2);
    dcell *GPTTDF=NULL;
    sptcellmulmat(&GPTTDF, recon->GA, recon->LR.V, 1);
    recon->LL.V=dcellcat(GPTTDF, VLo, 2);
    dcellfree(GPTTDF);
    dcellfree(ULo);
    dcellfree(VLo);
    {
	//Create piston and check board modes that are in NULL space of GA.
	dcell *NW=dcellnew(ndm,1);
	int nmod=2;//two modes.
	for(int idm=0; idm<ndm; idm++){
	    loc_create_map(recon->aloc[idm]);
	    const long nloc=recon->aloc[idm]->nloc;
	    NW->p[idm]=dnew(nloc, ndm*nmod);
	    //notice offset of 1 because map start count at 1
	    double *p=NW->p[idm]->p+nmod*idm*nloc-1;
	    for(long iloc=0; iloc<nloc; iloc++){
		p[iloc]=1;//piston mode
	    }
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
	dcellwrite(NW, "NW");
	//add to low rank terms.
	dcell *tmp=recon->LL.U;
	recon->LL.U=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellscale(NW, -1);
	tmp=recon->LL.V;
	recon->LL.V=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
    }
    //Remove empty cells.
    dcelldropempty(&recon->LR.U,2);
    dcelldropempty(&recon->LR.V,2);
    dcelldropempty(&recon->LL.U,2);
    dcelldropempty(&recon->LL.V,2);
    if(parms->save.setup && parms->save.recon){
	spcellwrite(recon->LR.M,"%s/LRM",dirsetup);
	dcellwrite(recon->LR.U,"%s/LRU",dirsetup);
	dcellwrite(recon->LR.V,"%s/LRV",dirsetup);
	spcellwrite(recon->LL.M,"%s/LLM.bin",dirsetup);//disable compression
	dcellwrite(recon->LL.U,"%s/LLU",dirsetup);
	dcellwrite(recon->LL.V,"%s/LLV",dirsetup); 
    }
    if(parms->lsr.alg==0){
	muv_chol_prep(&recon->LL);
	spcellfree(recon->LL.M);
	dcellfree(recon->LL.U);
	dcellfree(recon->LL.V);
    }
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
    dcellfree(recon->TTF);
    dcellfree(recon->PTTF);
    dcellfree(recon->RFlgs);
    dcellfree(recon->RFngs);
    spcellfree(recon->ZZT);
    spcellfree(recon->HXF); 
    spcellfree(recon->HXW);
    spcellfree(recon->HXWtomo);
    spfree(recon->W0); 
    dfree(recon->W1); 
    dcellfree(recon->NW);
    dfree(recon->fitwt);
    if(recon->ngsmod){
	dcellfree(recon->ngsmod->GM);
	dcellfree(recon->ngsmod->Rngs);
	dcellfree(recon->ngsmod->Pngs);
	dcellfree(recon->ngsmod->Mdm);
	dfree(recon->ngsmod->MCC);
	dfree(recon->ngsmod->MCC_OA);
	dfree(recon->ngsmod->IMCC);
	dfree(recon->ngsmod->IMCC_TT);
	dcellfree(recon->ngsmod->Ptt);
	free(recon->ngsmod);
    }

    int npsr = recon->npsr;
    int ndm = parms->ndm;
   
    locarrfree(recon->xloc, npsr); recon->xloc=NULL;
    locarrfree(recon->aloc, ndm); recon->aloc=NULL;
    locfree(recon->ploc); recon->ploc=NULL;
    free(recon->xloc_nx); 
    free(recon->xloc_ny);
    free(recon->aloc_nx);
    free(recon->aloc_ny);
    dcellfree(recon->aimcc);
    muv_free(&recon->RR);
    muv_free(&recon->RL);
    muv_free(&recon->FR);
    muv_free(&recon->FL);
    muv_free(&recon->LR);
    muv_free(&recon->LL);
    dfree(recon->neam);
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
