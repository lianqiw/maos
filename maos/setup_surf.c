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
static void 
setup_surf_tilt(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs){
    info("Setting up tilt surface (M3)\n");
    rectmap_t **tsurf=calloc(parms->ntsurf, sizeof(rectmap_t*));
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	char *fn=parms->tsurf[itsurf];
	info("Loading tilt surface from %s\n", fn);
	tsurf[itsurf]=rectmapread("%s",fn); 
    }
    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    int do_rot=(fabs(rot)>1.e-10);

    if(do_rot){
	locevl=locdup(aper->locs);
	locrot(locevl,rot);
    }else{
	locevl=aper->locs;
    }
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
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

	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    double bx=parms->evl.thetax[ievl]/mag;/*2010-04-02: do not put - sign */
	    double by=parms->evl.thetay[ievl]/mag;
	    double d_img_exit=fexit;
	    proj_rect_grid(mapsurf,alx,aly,locevl,scalex,scaley,
			   NULL,aper->opdadd->p[ievl]->p,scaleopd,
			   d_img_exit, het, bx,by);
	}

	for(int iwfs=0; iwfs<parms->nwfs && powfs; iwfs++){
	    const int ipowfs=parms->wfs[iwfs].powfs;
	    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	    const double hs=parms->powfs[ipowfs].hs;
	    loc_t *locwfs, *locwfsin;

	    if(powfs[ipowfs].nlocm){
		error("We don't handle this case yet. Think carefully when to apply shift.\n");
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

	    double d_img_focus=1./(1./ftel-1./hs)-ftel;
	    /*info2("iwfs%d: d_img_focus=%g\n",iwfs,d_img_focus); */
	    double d_img_exit=fexit+d_img_focus;
		
	    /*2010-04-02: do not put - sign */
	    double bx=parms->wfs[iwfs].thetax*(d_img_focus+ftel)/d_img_exit;
	    double by=parms->wfs[iwfs].thetay*(d_img_focus+ftel)/d_img_exit;
	    proj_rect_grid(mapsurf,alx,aly,locwfs,scalex,scaley,
			   NULL,powfs[ipowfs].opdadd->p[wfsind]->p,scaleopd,
			   d_img_exit, het, bx, by);
	    
	    if(do_rot){
		locfree(locwfs);
	    }
	}
    }
    if(do_rot){
	locfree(locevl);
    }
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	rectmapfree(tsurf[itsurf]);
    }
    free(tsurf);
}
/**
   Setup surface perpendicular to the beam by ray tracing from the surface to
   WFS and Science grid
 */
static void 
setup_surf_flat(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    info2("Setting up surface OPD (M1/M2/M3)\n");
    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    int do_rot=(fabs(rot)>1.e-10);
    if(do_rot){
	locevl=locdup(aper->locs);
	locrot(locevl,rot);
    }else{
	locevl=aper->locs;
    }
    const int nevl=parms->evl.nevl;
    const int nwfs=parms->nwfs;
    int *evlcover=malloc(nevl*sizeof(int));
    int *wfscover=malloc(nwfs*sizeof(int));
    int opdxcover;
    for(int isurf=0; isurf<parms->nsurf; isurf++){
	char *fn=parms->surf[isurf];
	if(!fn) continue;
	info("Loading surface OPD from %s\n", fn);
	map_t *surf=mapread("%s",fn);
	surf->ox+=parms->aper.misreg[0];
	surf->oy+=parms->aper.misreg[1];

	const char *strevl=search_header(surf->header, "SURFEVL");
	const char *strwfs=search_header(surf->header, "SURFWFS");
	const char *stropdx=search_header(surf->header, "SURFOPDX");
	if(!strevl){
	    warning2("surf[%d] does not contain SURFEVL\n", isurf);
	    for(int ievl=0; ievl<nevl; ievl++){
		evlcover[ievl]=1;
	    }
	}else{
	    readstr_intarr_nmax(&evlcover, nevl, strevl);
	}
	if(!strwfs){
	    warning2("surf[%d] does not contain SURFWFS\n", isurf);
	    for(int iwfs=0;iwfs<nwfs; iwfs++){
		wfscover[iwfs]=1;
	    }
	}else{
	    readstr_intarr_nmax(&wfscover, nwfs, strwfs);
	}
	if(!stropdx){
	    opdxcover=1;
	}else{
	    opdxcover=(int)readstr_num(stropdx, NULL);
	}
	double hl=surf->h;
	for(int ievl=0; ievl<nevl; ievl++){
	    if(!evlcover[ievl]){
		warning2("Skip evl %d for surface %s\n", ievl, parms->surf[isurf]);
		continue;
	    }
	    if(!aper->opdadd->p[ievl]){
		aper->opdadd->p[ievl]=dnew(aper->locs->nloc, 1);
	    }
	    const double displacex=parms->evl.thetax[ievl]*hl;
	    const double displacey=parms->evl.thetay[ievl]*hl;
	    
	    if(do_rot){
		prop_grid(surf, locevl, NULL, aper->opdadd->p[ievl]->p, 
			  1, displacex, displacey, 1, 0, 0, 0);
	    }else{
		prop_grid_stat(surf, aper->locs->stat, 
			       aper->opdadd->p[ievl]->p, 
			       1, displacex, displacey, 1, 0, 0, 0);
	    }
	}
	for(int iwfs=0; iwfs<parms->nwfs && powfs; iwfs++){
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
    if(do_rot){
	locfree(locevl);
    }
    free(evlcover);
    free(wfscover);
}

void setup_surf(const PARMS_T *parms, APER_T *aper, POWFS_T *powfs, RECON_T *recon){
    if(parms->nsurf<=0 && parms->ntsurf<=0){
	info2("No surfaces to setup\n");
	return;
    }
    if(aper->opdadd){
	error("Who sets aper->opdadd?\n");
    }
    aper->opdadd=dcellnew(parms->evl.nevl,1);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	aper->opdadd->p[ievl]=dnew(aper->locs->nloc, 1);
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	powfs[ipowfs].opdadd=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    powfs[ipowfs].opdadd->p[jwfs]=dnew(powfs[ipowfs].npts, 1);
	}
    }

    if(parms->ntsurf>0){
	setup_surf_tilt(parms, aper, powfs);
    }
    if(parms->surf>0){
	setup_surf_flat(parms, aper, powfs, recon);
    }
    if(parms->save.setup){
	dcellwrite(aper->opdadd, "%s/surfevl.bin", dirsetup);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    dcellwrite(powfs[ipowfs].opdadd, "%s/surfpowfs_%d.bin", dirsetup, ipowfs);
	}
	if(recon->opdxadd) dcellwrite(recon->opdxadd, "%s/surfopdx", dirsetup);
    }
}
