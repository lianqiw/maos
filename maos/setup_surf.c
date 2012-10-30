/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   \file setup_surf.c
   Setup additional NCPA surface.
*/
#include "maos.h"
#include "setup_surf.h"
#include "recon_utils.h"
#include "setup_powfs.h"
/**
   Setup tilted surface (M3) by ray tracing from the tilted surface to WFS and
   Science grid.

   \todo Merge setup_tsurf, setup_surf, and setup_powfs_ncpa. Be careful about
   which NCPA is compensated by WFS offsets. setup_tsurf and setup_surf act on
   surfaces intersect with most WFS and science, while setup_powfs_ncpa act on
   individual surfaces for each WFS. Maybe it is good idea to keep them
   separate.


   2011-09-07: 
   Relocate simu->surfwfs to powfs.opdadd In order to perserve surface for different seeds
*/
/**
   Propagate tilt surface from star at hs, along direction (thetax, thetay), to loc.
*/
static void tsurf2loc(rectmap_t **tsurf, int ntsurf, dmat *opd, loc_t *locin, double thetax, double thetay, double hs, double rot){
    loc_t *locuse=NULL;
    if(fabs(rot)>1e-10){
	locuse=locdup(locin);
	locrot(locuse, rot);
    }else{
	locuse=locin;
    }
    for(int itsurf=0; itsurf<ntsurf; itsurf++){
	const double alx=tsurf[itsurf]->txdeg/180*M_PI;
	const double aly=tsurf[itsurf]->tydeg/180*M_PI;
	const double ftel=tsurf[itsurf]->ftel;
	const double fexit=tsurf[itsurf]->fexit;
	const double fsurf=tsurf[itsurf]->fsurf;
	const double mag=fexit/ftel;
	const double scalex=-mag;
	const double scaley=mag;
	const double scaleopd=-2;
	const double het=fexit-fsurf;/*distance between exit pupil and M3. */
	rectmap_t *mapsurf=tsurf[itsurf];

	double d_img_focus=1./(1./ftel-1./hs)-ftel;
	/*info2("iwfs%d: d_img_focus=%g\n",iwfs,d_img_focus); */
	double d_img_exit=fexit+d_img_focus;
		
	/*2010-04-02: do not put - sign */
	double bx=thetax*(d_img_focus+ftel)/d_img_exit;
	double by=thetay*(d_img_focus+ftel)/d_img_exit;
	proj_rect_grid(mapsurf,alx,aly,locuse,scalex,scaley, NULL,opd->p,scaleopd, d_img_exit, het, bx, by);
    }
    if(locuse!=locin){
	locfree(locuse);
    }
}

static void 
setup_surf_tilt(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    info("Setting up tilt surface (M3)\n");
    rectmap_t **tsurf=calloc(parms->ntsurf, sizeof(rectmap_t*));
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	char *fn=parms->tsurf[itsurf];
	info("Loading tilt surface from %s\n", fn);
	tsurf[itsurf]=rectmapread("%s",fn); 
    }
    const double rot=-parms->aper.rotdeg/180.*M_PI;

    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	tsurf2loc(tsurf, parms->ntsurf, aper->opdadd->p[ievl], aper->locs, 
		  parms->evl.thetax[ievl], parms->evl.thetay[ievl], parms->evl.hs[ievl], rot);
    }
    if(parms->sim.ncpa_calib){
	for(int ifit=0; ifit<parms->sim.ncpa_ndir; ifit++){
	    tsurf2loc(tsurf, parms->ntsurf, aper->opdfloc->p[ifit], recon->floc, 
		      parms->sim.ncpa_thetax[ifit], parms->sim.ncpa_thetay[ifit], parms->sim.ncpa_hs[ifit], rot);
	}
    }

    for(int iwfs=0; iwfs<parms->nwfs && powfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	loc_t *locwfsin;
	    
	if(powfs[ipowfs].nlocm){
	    error("We don't handle this case yet. Think carefully when to apply shift.\n");
	    int ilocm=powfs[ipowfs].nlocm>1?parms->powfs[ipowfs].wfsind[iwfs]:0;
	    locwfsin=powfs[ipowfs].locm[ilocm];
	}else{
	    locwfsin=powfs[ipowfs].loc;
	}
	tsurf2loc(tsurf, parms->ntsurf, powfs[ipowfs].opdadd->p[wfsind], locwfsin, 
		  parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, parms->powfs[ipowfs].hs, rot);
    }
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	rectmapfree(tsurf[itsurf]);
    }
    free(tsurf);
}

typedef struct{
    const PARMS_T *parms;
    APER_T *aper;
    POWFS_T *powfs;
    RECON_T *recon;
    loc_t *locevl;
    double rot;
    map_t *surf;
    int isurf;
    int nevl;
    int nwfs;
    int nncpa;
    int *evlcover;
    int *wfscover;
    int *ncpacover;
    int opdxcover;
    int do_rot;
}SURF_DATA;

static void prop_surf_evl(thread_t *info){
    SURF_DATA *data=info->data;
    const PARMS_T *parms=data->parms;
    const APER_T *aper=data->aper;
    const map_t *surf=data->surf;
    const double hl=surf->h;
    const int *evlcover=data->evlcover;
    const int do_rot=data->do_rot;
    const int isurf=data->isurf;
    loc_t *locevl=data->locevl;
    for(int ievl=info->start; ievl<info->end; ievl++){
	if(!evlcover[ievl]){
	    warning2("Skip evl %d for surface %s\n", ievl, parms->surf[isurf]);
	    continue;
	}
	const double displacex=parms->evl.thetax[ievl]*hl;
	const double displacey=parms->evl.thetay[ievl]*hl;
	const double scale=1-hl/parms->evl.hs[ievl];
	if(do_rot){
	    prop_grid(surf, locevl, NULL, aper->opdadd->p[ievl]->p, 
		      1, displacex, displacey, scale, 0, 0, 0);
	}else{
	    prop_grid_stat(surf, aper->locs->stat, 
			   aper->opdadd->p[ievl]->p, 
			   1, displacex, displacey, scale, 0, 0, 0);
	}
    }
}

static void prop_surf_ncpa(thread_t *info){
    SURF_DATA *data=info->data;
    const PARMS_T *parms=data->parms;
    const APER_T *aper=data->aper;
    const RECON_T *recon=data->recon;
    const map_t *surf=data->surf;
    const double hl=surf->h;
    const int *ncpacover=data->ncpacover;
    for(int idir=info->start; idir<info->end; idir++){
	if(!ncpacover[idir]) continue;
	const double displacex=parms->sim.ncpa_thetax[idir]*hl;
	const double displacey=parms->sim.ncpa_thetay[idir]*hl;
	const double scale=1.-hl/parms->sim.ncpa_hs[idir];
	prop_grid(surf, recon->floc, NULL,
		  aper->opdfloc->p[idir]->p, 
		  1, displacex, displacey, scale, 0, 0, 0);	
    }
}

static void prop_surf_wfs(thread_t *info){
    SURF_DATA *data=info->data;
    const PARMS_T *parms=data->parms;
    const POWFS_T *powfs=data->powfs;
    const map_t *surf=data->surf;
    const double hl=surf->h;
    const int *wfscover=data->wfscover;
    const int isurf=data->isurf;
    const int do_rot=data->do_rot;
    const double rot=data->rot;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	if(!wfscover[iwfs]){
	    warning2("Skip wfs %d for surface %s\n", iwfs, parms->surf[isurf]);
	    continue;
	}
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	const double hs=parms->powfs[ipowfs].hs;
	const double scale=1.-hl/hs;
	const double displacex=parms->wfs[iwfs].thetax*hl+powfs[ipowfs].misreg[wfsind][0];
	const double displacey=parms->wfs[iwfs].thetay*hl+powfs[ipowfs].misreg[wfsind][1];

	loc_t *locwfs, *locwfsin;
	if(powfs[ipowfs].locm){
	    int ilocm=powfs[ipowfs].nlocm>1?parms->powfs[ipowfs].wfsind[iwfs]:0;
	    locwfsin=powfs[ipowfs].locm[ilocm];
	}else{
	    locwfsin=powfs[ipowfs].loc;
	}
	if(do_rot){
	    locwfs=locdup(locwfsin);
	    locrot(locwfs,rot);
	}else{
	    locwfs=locwfsin;
	}
	prop_grid(surf, locwfs, NULL, powfs[ipowfs].opdadd->p[wfsind]->p, 
		  1, displacex, displacey, scale, 1., 0, 0); 
	if(do_rot){
	    locfree(locwfs);
	}
    }
}

/**
   Setup surface perpendicular to the beam by ray tracing from the surface to
   WFS and Science grid
*/
static void 
setup_surf_perp(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    info2("Setting up surface OPD (M1/M2/M3)\n");
    if(fabs(parms->aper.misreg[0])>EPS || fabs(parms->aper.misreg[1])>EPS){
	warning("Please adjust telescope surface ox, oy to account for misregistration. Not doing "
		"in maos because some surfaces may belong to instrument.\n");
    }
    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    const int nncpa=parms->sim.ncpa_ndir;
    int *evlcover=malloc(nevl*sizeof(int));
    int *wfscover=malloc(nwfs*sizeof(int));
    int *ncpacover=malloc(nncpa*sizeof(int));
    int opdxcover;
    SURF_DATA sdata={parms, aper, powfs, recon, NULL, rot, NULL, 0, nevl, nwfs, nncpa, evlcover, wfscover, ncpacover, 0, 0};
    const int nthread=NTHREAD;
    thread_t  tdata_wfs[nthread], tdata_evl[nthread], tdata_ncpa[nthread];
    thread_prep(tdata_evl, 0, nevl, nthread, prop_surf_evl, &sdata);
    thread_prep(tdata_wfs, 0, nwfs, nthread, prop_surf_wfs, &sdata);
    thread_prep(tdata_ncpa, 0, nncpa, nthread, prop_surf_ncpa, &sdata);
    

    for(int isurf=0; isurf<parms->nsurf; isurf++){
	char *fn=parms->surf[isurf];
	if(!fn) continue;
	info("Loading surface OPD from %s\n", fn);
	map_t *surf=mapread("%s",fn);
	//dwrite((dmat*)surf, "surf_%d", isurf);
	const char *strname=search_header(surf->header, "SURFNAME");
	const char *strevl=search_header(surf->header, "SURFEVL");
	const char *strwfs=search_header(surf->header, "SURFWFS");
	const char *stropdx=search_header(surf->header, "SURFOPDX");
	int do_rot=0;
	if(strname && !strcmp(strname, "M1")){
	    warning("Rotate loc for M1\n");
	    do_rot=(fabs(rot)>1.e-10);
	}
	if(do_rot){
	    locevl=locdup(aper->locs);
	    locrot(locevl,rot);
	}else{
	    locevl=aper->locs;
	}
	if(!strevl){
	    warning2("surf[%d] does not contain SURFEVL\n", isurf);
	    for(int ievl=0; ievl<nevl; ievl++){
		evlcover[ievl]=1;
	    }
	}else{
	    readstr_intarr_nmax(&evlcover, nevl, strevl);
	}
	int evlct=0;
	for(int ievl=0; ievl<nevl; ievl++){
	    evlct+=evlcover[ievl]?1:0;
	}
	if(evlct==0){
	    for(int idir=0; idir<nncpa; idir++){
		ncpacover[idir]=0;
	    }
	}else if(evlct==nevl){
	    for(int idir=0; idir<nncpa; idir++){
		ncpacover[idir]=1;
	    }
	}else{
	    error("Not handled\n");
	}
	if(!strwfs){
	    warning2("surf[%d] does not contain SURFWFS\n", isurf);
	    for(int iwfs=0;iwfs<nwfs; iwfs++){
		wfscover[iwfs]=1;
	    }
	}else{
	    if(nwfs>9 && parms->sim.skysim){
		warning("There are many NGS stars. Replicate surface config from the 8th star\n");
		readstr_intarr_nmax(&wfscover, 9, strwfs);
		wfscover=realloc(wfscover, sizeof(int)*nwfs);
		for(int i=9; i<nwfs; i++){
		    wfscover[i]=wfscover[8];
		}
	    }else{
		readstr_intarr_nmax(&wfscover, nwfs, strwfs);
	    }
	}
	if(!stropdx){
	    opdxcover=1;
	}else{
	    opdxcover=(int)readstr_num(stropdx, NULL);
	}
	sdata.locevl=locevl;
	sdata.surf=surf;
	sdata.opdxcover=opdxcover;
	sdata.do_rot=do_rot;
	sdata.isurf=isurf;
	double hl=surf->h;
	CALL_THREAD(tdata_evl, nthread, 0);
	if(powfs){
	    CALL_THREAD(tdata_wfs, nthread, 0);
	}
	if(parms->sim.ncpa_calib){
	    CALL_THREAD(tdata_ncpa, nthread, 0);
	}
	/*The following cannot be parallelized as done for evl, wfs because the
	  destination opdx is not in sequence.*/
	if(parms->sim.idealfit && recon && opdxcover){
	    if(!recon->opdxadd){
		recon->opdxadd=dcellnew(parms->atmr.nps, 1);
	    }
	    double distmin=INFINITY;
	    int jpsr=-1;
	    /*Select the layer that is closed to the surface. */
	    for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
		double dist=fabs(parms->atmr.ht[ipsr]-hl);
		if(dist < distmin){
		    jpsr=ipsr;
		    distmin=dist;
		}
	    }
	    loc_t *xloc;
	    if(do_rot){
		xloc = locdup(recon->xloc[jpsr]);
		locrot(xloc, rot);
	    }else{
		xloc = recon->xloc[jpsr];
	    }
	    loc_t *surfloc=mksqloc_map(surf);
	    dsp *H=mkhb(xloc, surfloc, NULL, 0, 0, 1, 0, 0);
	    double scale=pow(surf->dx/xloc->dx,2);
	    if(!recon->opdxadd->p[jpsr]){
		recon->opdxadd->p[jpsr]=dnew(xloc->nloc, 1);
	    }
	    spmulvec(recon->opdxadd->p[jpsr]->p, H, surf->p, scale);
	    if(do_rot) locfree(xloc);
	    locfree(surfloc);
	}
	mapfree(surf);
    }
    if(locevl!=aper->locs){
	locfree(locevl);
    }
    free(evlcover);
    free(wfscover);
}

/** We trace rays from Science focal plan OPD to ploc along evaluation
    directions (type=1) or on axis only (type=2).*/
static void FitR_NCPA(dcell **xout, RECON_T *recon, APER_T *aper){
    const PARMS_T *parms=recon->parms;
    dcell *xp=NULL;
    if(aper->opdfloc){
	xp=dcelldup(aper->opdfloc);
    }else{
	xp=dcellnew(parms->sim.ncpa_ndir, 1);
	for(int ievl=0; ievl<parms->sim.ncpa_ndir; ievl++){
	    xp->p[ievl]=dnew(recon->floc->nloc,1);
	    prop_nongrid(aper->locs, aper->opdadd->p[ievl]->p,
			 recon->floc, NULL, xp->p[ievl]->p, 1, 0, 0, 1, 0, 0);
	}
    }
    applyW(xp, recon->W0, recon->W1, parms->sim.ncpa_wt);
    sptcellmulmat_thread(xout, recon->HA_ncpa, xp, 1);
    dcellfree(xp);
}
void FitL_NCPA(dcell **xout, const void *A, 
	       const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    dcell *xp=NULL;
    spcellmulmat_thread(&xp, recon->HA_ncpa, xin, 1.);
    applyW(xp, recon->W0, recon->W1, parms->sim.ncpa_wt);
    sptcellmulmat_thread(xout, recon->HA_ncpa, xp, alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp,recon->fitNW, xin, "tn", 1);
    dcellmm(xout,recon->fitNW, xp, "nn", alpha);
    dcellfree(xp);
    if(recon->actslave){
	spcellmulmat(xout, recon->actslave, xin, 1);
    }
}
static void setup_recon_HAncpa(RECON_T *recon, const PARMS_T *parms){
    const int nevl=parms->sim.ncpa_ndir;
    const int ndm=parms->ndm;
    recon->HA_ncpa=spcellnew(nevl, ndm);
    PDSPCELL(recon->HA_ncpa,HA);
    info2("Generating HA ");TIC;tic;
    for(int ievl=0; ievl<nevl; ievl++){
	double hs=parms->sim.ncpa_hs[ievl];
	for(int idm=0; idm<ndm; idm++){
	    if(parms->sim.ncpa_calib==2 && idm>0){
		continue;
	    }
	    const double ht=parms->dm[idm].ht;
	    const double scale=1.-ht/hs;
	    double displace[2];
	    displace[0]=parms->sim.ncpa_thetax[ievl]*ht;
	    displace[1]=parms->sim.ncpa_thetay[ievl]*ht;
	    HA[idm][ievl]=mkh(recon->aloc[idm], recon->floc, NULL,
			      displace[0], displace[1], 
			      scale,parms->dm[idm].cubic,parms->dm[idm].iac);
	}
    }
    toc2(" ");
    if(parms->save.setup){
	spcellwrite(recon->HA,"%s/HA_ncpa",dirsetup);
    }
}
#include "mtch.h"
#include "genseotf.h"
void setup_surf(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    if(parms->nsurf || parms->ntsurf){
	if(!aper->opdadd){
	    aper->opdadd=dcellnew(parms->evl.nevl,1);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		aper->opdadd->p[ievl]=dnew(aper->locs->nloc, 1);
	    }
	}
	if(!aper->opdfloc && parms->sim.ncpa_calib){
	    aper->opdfloc=dcellnew(parms->sim.ncpa_ndir,1);
	    for(int idir=0; idir<parms->sim.ncpa_ndir; idir++){
		aper->opdfloc->p[idir]=dnew(recon->floc->nloc, 1);
	    }
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(!powfs[ipowfs].opdadd){
		powfs[ipowfs].opdadd=dcellnew(parms->powfs[ipowfs].nwfs, 1);
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    powfs[ipowfs].opdadd->p[jwfs]=dnew(powfs[ipowfs].npts, 1);
		}
	    }
	}

	if(parms->ntsurf>0){
	    setup_surf_tilt(parms, aper, powfs, recon);
	}
	if(parms->nsurf>0){
	    setup_surf_perp(parms, aper, powfs, recon);
	}
    }
    if(parms->sim.ncpa_calib){//calibrate NCPA
	int any_evl=0;
	if(aper->opdadd){
	    for(int i=0; i<parms->evl.nevl; i++){
		if(aper->opdadd->p[i]){
		    any_evl=1;
		}
	    }
	}
	int any_wfs=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(powfs[ipowfs].opdadd){
		any_wfs=1;
	    }
	}
	if(any_evl){
	    if(fabs(parms->aper.rotdeg)>0){
		error("Not handling aper.rotdeg\n");
	    }
	    info("calibrating NCPA\n");
	    setup_recon_HAncpa(recon, parms);
	    dcell *rhs=NULL;
	    FitR_NCPA(&rhs, recon, aper);
	    int maxit=40;
	    pcg(&recon->dm_ncpa, FitL_NCPA, recon, NULL, NULL, rhs, 1, maxit);
	    dcellfree(rhs);
	    dcellwrite(recon->dm_ncpa, "dm_ncpa");
	    spcellfree(recon->HA_ncpa);
	}
	dcell *dm_ncpa=dcellref(recon->dm_ncpa);
	setup_powfs_calib(parms, powfs, recon->aloc, dm_ncpa);

	/* to do 
	   dm flat
	   matched filter
	   genseotf() reentrant. ok
	   genselotf() reentrant. ok
	   gensepsf() reentrant. ok
	   gensei() reentrant. ok.
	   genmtch() reentrant. ok
	*/
	dcellfree(dm_ncpa);
	dcellfree(aper->opdfloc);
    }
    if(parms->save.setup){
	dcellwrite(aper->opdadd, "%s/surfevl.bin",  dirsetup);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    dcellwrite(powfs[ipowfs].opdadd, "%s/surfpowfs_%d.bin", dirsetup,  ipowfs);
	}
	if(recon->opdxadd) dcellwrite(recon->opdxadd, "%s/surfopdx",  dirsetup);
    }
}
