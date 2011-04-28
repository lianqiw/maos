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
 */
void setup_tsurf(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    if(parms->ntsurf<=0) return;
    info("Setting up tilt surface (M3)\n");
    if(!parms->sim.cachesurf){
	error("Not yet finished\n");
    }
    if(!simu->surfwfs){
	simu->surfwfs=dcellnew(parms->nwfs,1);
    }
    if(!simu->surfevl){
	simu->surfevl=dcellnew(parms->evl.nevl,1);
    }
 
    simu->tsurf=calloc(parms->ntsurf, sizeof(rectmap_t*));
    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	char *fn=find_file(parms->tsurf[itsurf]);
	info("Loading tilt surface from %s\n", fn);
	simu->tsurf[itsurf]=rectmapread("%s",fn); 
	free(fn);
    }
    if(!parms->sim.cachesurf){
	return;
    }

    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    int do_rot=(fabs(rot)>1.e-10);

    if(do_rot){
	locevl=locdup(simu->aper->locs);
	locrot(locevl,rot);
    }else{
	locevl=simu->aper->locs;
    }

    for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	const double alx=simu->tsurf[itsurf]->txdeg/180*M_PI;
	const double aly=simu->tsurf[itsurf]->tydeg/180*M_PI;
	const double ftel=simu->tsurf[itsurf]->ftel;
	const double fexit=simu->tsurf[itsurf]->fexit;
	const double fsurf=simu->tsurf[itsurf]->fsurf;
	const double mag=fexit/ftel;
	const double scalex=-mag;
	const double scaley=mag;
	const double scaleopd=-2;
	const double het=fexit-fsurf;//distance between exit pupil and M3.
	rectmap_t *mapsurf=simu->tsurf[itsurf];

	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->surfevl->p[ievl]){
		simu->surfevl->p[ievl]=dnew(simu->aper->locs->nloc, 1);
	    }
	    double bx=parms->evl.thetax[ievl]/mag;//2010-04-02: do not put - sign
	    double by=parms->evl.thetay[ievl]/mag;
	    double d_img_exit=fexit;
	    proj_rect_grid(mapsurf,alx,aly,locevl,scalex,scaley,
			   NULL,simu->surfevl->p[ievl]->p,scaleopd,
			   d_img_exit, het, bx,by);
	}

	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    double hs=parms->powfs[ipowfs].hs;
	    loc_t *locwfs, *locwfsin;
	    if(!simu->surfwfs->p[iwfs]){
		simu->surfwfs->p[iwfs]=dnew(simu->powfs[ipowfs].npts,1);
	    }
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
	    //info2("iwfs%d: d_img_focus=%g\n",iwfs,d_img_focus);
	    double d_img_exit=fexit+d_img_focus;
		
	    //2010-04-02: do not put - sign
	    double bx=parms->wfs[iwfs].thetax*(d_img_focus+ftel)/d_img_exit;
	    double by=parms->wfs[iwfs].thetay*(d_img_focus+ftel)/d_img_exit;
	    proj_rect_grid(mapsurf,alx,aly,locwfs,scalex,scaley,
			   NULL,simu->surfwfs->p[iwfs]->p,scaleopd,
			   d_img_exit, het, bx, by);
	    
	    if(do_rot){
		locfree(locwfs);
	    }
	}
	
        dcellwrite(simu->surfwfs,"surfwfs.bin");
	dcellwrite(simu->surfevl,"surfevl.bin");
	//exit(0);
	
    }
    if(do_rot){
	locfree(locevl);
    }
    if(parms->sim.cachesurf){
	for(int itsurf=0; itsurf<parms->ntsurf; itsurf++){
	    rectmapfree(simu->tsurf[itsurf]);
	}
    }
}
/**
   Setup surface perpendicular to the beam by ray tracing from the surface to
   WFS and Science grid
 */
void setup_surf(SIM_T*simu){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    if(parms->nsurf<=0){
	info2("There are no surface map\n");
	return;
    }else{
	info2("Setting up surface OPD (M1/M2/M3)\n");
    }
    if(!parms->sim.cachesurf){
	error("Not yet finished\n");
    }
    if(!simu->surfwfs){
	simu->surfwfs=dcellnew(parms->nwfs,1);
    }
    if(!simu->surfevl){
	simu->surfevl=dcellnew(parms->evl.nevl,1);
    }
    if(!simu->surfopdx && parms->sim.fitonly){
	simu->surfopdx=dcellnew(parms->atmr.nps, 1);
    }
    simu->surf=calloc(parms->nsurf, sizeof(map_t*));
    for(int isurf=0; isurf<parms->nsurf; isurf++){
	if(!parms->surf[isurf]) continue;
	char *fn=find_file(parms->surf[isurf]);
	info("Loading surface OPD from %s\n", fn);
	simu->surf[isurf]=mapread("%s",fn); free(fn);
	simu->surf[isurf]->ox+=parms->aper.misreg[0];
	simu->surf[isurf]->oy+=parms->aper.misreg[1];
    }
    if(!parms->sim.cachesurf){
	return;
    }

    loc_t *locevl;
    const double rot=-parms->aper.rotdeg/180.*M_PI;
    int do_rot=(fabs(rot)>1.e-10);
    if(do_rot){
	locevl=locdup(simu->aper->locs);
	locrot(locevl,rot);
    }else{
	locevl=simu->aper->locs;
    }
    
  
    for(int isurf=0; isurf<parms->nsurf; isurf++){
	if(!simu->surf[isurf]) continue;
	double hl=simu->surf[isurf]->h;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->surfevl->p[ievl]){
		simu->surfevl->p[ievl]=dnew(simu->aper->locs->nloc, 1);
	    }
	    double displacex=parms->evl.thetax[ievl]*hl;
	    double displacey=parms->evl.thetay[ievl]*hl;
	    
	    if(do_rot){
		prop_grid(simu->surf[isurf], locevl, simu->surfevl->p[ievl]->p, 
			  1, displacex, displacey, 1, 0, 0, 0);
	    }else{
		prop_grid_stat(simu->surf[isurf], simu->aper->locs->stat, 
			       simu->surfevl->p[ievl]->p, 
			       1, displacex, displacey, 1, 0, 0, 0);
	    }
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	    double hs=parms->powfs[ipowfs].hs;
	    const double scale=1.-hl/hs;
	    const double displacex=parms->wfs[iwfs].thetax*hl+powfs[ipowfs].misreg[wfsind][0];
	    const double displacey=parms->wfs[iwfs].thetay*hl+powfs[ipowfs].misreg[wfsind][1];
	    if(!simu->surfwfs->p[iwfs]){
		simu->surfwfs->p[iwfs]=dnew(simu->powfs[ipowfs].npts,1);
	    }
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
	    prop_grid(simu->surf[isurf], locwfs, simu->surfwfs->p[iwfs]->p, 
		      1, displacex, displacey, scale, 1., 0, 0); 
	    if(do_rot){
		locfree(locwfs);
	    }
	}
	if(parms->sim.fitonly){
	    double distmin=INFINITY;
	    int jpsr=-1;
	    //Select the layer that is closed to the surface.
	    for(int ipsr=0; ipsr<parms->atmr.nps; ipsr++){
		double dist=fabs(parms->atmr.ht[ipsr]-hl);
		if(dist < distmin){
		    jpsr=ipsr;
		    distmin=dist;
		}
	    }
	    loc_t *xloc;
	    if(do_rot){
		xloc = locdup(simu->recon->xloc[jpsr]);
		locrot(xloc, rot);
	    }else{
		xloc = simu->recon->xloc[jpsr];
	    }
	    loc_t *surfloc=mksqloc_map(simu->surf[isurf]);
	    dsp *H=mkhb(xloc, surfloc, NULL, 0, 0, 1, 0, 0);
	    double scale=pow(simu->surf[isurf]->dx/xloc->dx,2);
	    if(!simu->surfopdx->p[jpsr]){
		simu->surfopdx->p[jpsr]=dnew(xloc->nloc, 1);
	    }
	    spmulvec(simu->surfopdx->p[jpsr]->p, H, simu->surf[isurf]->p, scale);
	    if(do_rot) locfree(xloc);
	    locfree(surfloc);
	    dwrite(simu->surfopdx->p[jpsr], "surfopdx_%d", jpsr);
	}
    }
    if(do_rot){
	locfree(locevl);
    }
    if(parms->sim.cachesurf){
	maparrfree(simu->surf, parms->nsurf);
    }
}
