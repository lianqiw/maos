/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "../sys/sys.h"
#include "mathdef.h"
#include "loc.h"
#include "map.h"

struct loc_private{
};
/**
   Free pts_t data
*/
void ptsfree_do(pts_t *pts){
    if(pts && pts->nref && !atomicadd(pts->nref, -1)){
	if(pts->origx){
	    free(pts->origx);
	    free(pts->origy);
	}
	if(pts->map){
	    mapfree(pts->map);
	}
	free(pts->nref);
    }
    free(pts);
}

/**
   Free the MAP in loc_t
*/
void loc_free_map(loc_t *loc){
    if(loc->map){
	mapfree(loc->map);
	loc->map=NULL;
    }
}
/**
   Free the stat in loc_t
*/
void loc_free_stat(loc_t *loc){
    if(loc->stat){
	free(loc->stat->cols);
	free(loc->stat);
	loc->stat=NULL;
    }
}
/**
   Free loc_t data
 */
void locfree_do(loc_t *loc){
    if(loc && loc->nref && !atomicadd(loc->nref, -1)){
	loc_free_stat(loc);
	loc_free_map(loc);
	free(loc->locx);
	free(loc->locy);
	free(loc->nref);
    }
    free(loc);
}

/**
   Free loc_t array
*/
void locarrfree_do(loc_t **loc, int nloc){
    if(!loc) return;
    for(int iloc=0; iloc<nloc; iloc++){
	locfree(loc[iloc]);
    }
    free(loc);
}

/**
   Create a loc with nloc elements.
*/
loc_t *locnew(long nloc,real dx, real dy){
    loc_t *loc=mycalloc(1,loc_t);
    loc->id=M_LOC;
    if(nloc>0){
	loc->locx=mycalloc(nloc,real);
	loc->locy=mycalloc(nloc,real);
	loc->nref=mycalloc(1, int); loc->nref[0]=1;
    }
    loc->nloc=nloc;
    if(!isnan(dx)) loc->dx=dx;
    if(!isnan(dy)) loc->dy=dy;

    return loc;
}
/**
   Reference an existing loc
 */
loc_t *locref(const loc_t *in){
    loc_t *out=mycalloc(1, loc_t);
    memcpy(out, in, sizeof(loc_t));
    if(out->nref){
	atomicadd(out->nref, 1);
    }
    return out;
}
/**
   Create a pts with nsa, dsa, nx, dx
*/
pts_t *ptsnew(long nsa, real dsax, real dsay, long nx, real dx, real dy){
    pts_t *pts=mycalloc(1,pts_t);
    pts->id=M_LOC;
    pts->origx=mycalloc(nsa,real);
    pts->origy=mycalloc(nsa,real);
    pts->dsa=dsax;
    pts->dsay=dsay;
    pts->nref=mycalloc(1, int); pts->nref[0]=1;
    pts->nx=nx;
    pts->ny=nx;
    pts->dx=dx;
    pts->dy=dy;
    pts->nsa=nsa;
    return pts;
}
uint32_t lochash(const loc_t *loc, uint32_t key){
    key=hashlittle(loc->locx, loc->nloc*sizeof(real), key);
    key=hashlittle(loc->locy, loc->nloc*sizeof(real), key);
    return key;
}
/**
   Create an vector to embed OPD into square array for FFT purpose. oversize is
2 for fft.  */
lmat *loc_create_embed(long *nembed, const loc_t *loc, real oversize, int fftpad){
    real xmin,xmax,ymin,ymax;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    const real dx_in1=1./loc->dx;
    const real dy_in1=1./loc->dy;
    long nx=(long)round((xmax-xmin)*dx_in1)+1;
    long ny=(long)round((ymax-ymin)*dy_in1)+1;
    long nxy=(long)ceil((nx>ny?nx:ny)*oversize);/*minimum size */
    if(fftpad){
	nxy=nextfftsize(nxy);
    }
    if(*nembed<=0){
	*nembed=nxy;
    }else{
	if(*nembed<(long)(nxy*0.6)){
	    error("Supplied nembed %ld is too small, recommend %ld\n",*nembed, nxy);
	}else if(*nembed<nxy){
	    warning("Supplied nembed %ld is too small, recommend %ld\n",*nembed, nxy);
	}
	nxy=*nembed;
    }
    xmin-=(nxy-nx+1)/2*loc->dx;//odd -> pad on left
    ymin-=(nxy-ny+1)/2*loc->dy;
    lmat *embed=lnew(loc->nloc, 1);
    for(int iloc=0; iloc<loc->nloc; iloc++){
	long ix=(long)round((loc->locx[iloc]-xmin)*dx_in1);
	long iy=(long)round((loc->locy[iloc]-ymin)*dy_in1);
	embed->p[iloc]=ix+iy*nxy;
    }
    return embed;
}
/**
   Create a map for loc with padding of 1. 
*/
void loc_create_map(const loc_t *loc){
    loc_create_map_npad((loc_t*)loc,0,0,0);
}
PNEW(maplock);
DEF_ENV_FLAG(LOC_MAP_EXTEND, 1);

/**
   Create a map for loc so that we can obtain the index in loc by x,y
   coordinate. Useful in ray tracing (accphi.c). If extend is true, the invalid
   region in the map is replaced with closest valid point. This extends the
   coverage of the grid automatically.
*/
void loc_create_map_npad(const loc_t *loc, int npad, int nx, int ny){
    LOCK(maplock);
    if(loc->map){
	if(loc->npad<npad){
	    print_backtrace();
	    warning("loc->map arealdy existed with npad=%d. Cannot redo for %d\n", loc->npad, npad);
	}
	UNLOCK(maplock);
	return;
    }
    if(loc->nloc==0){
	UNLOCK(maplock);
	return;
    }
    ((loc_t*)loc)->npad = npad;/*just record the information. */
    real xmin,xmax,ymin,ymax;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    /*padding the map. normally don't need. */
    if(npad>0){
	xmin-=npad*fabs(loc->dx);
	ymin-=npad*fabs(loc->dy);
	xmax+=npad*fabs(loc->dx);
	ymax+=npad*fabs(loc->dy);
    }
    const real dx_in1=1./fabs(loc->dx);
    const real dy_in1=1./fabs(loc->dy);
    long map_nx=(long) round((xmax-xmin)*dx_in1)+1;
    long map_ny=(long) round((ymax-ymin)*dy_in1)+1;
    if(nx && ny){
	if(map_nx>nx || map_ny>ny){
	    error("Specified size %dx%d is too small, need at least %ldx%ld\n", 
		  nx, ny, map_nx, map_ny);
	}
	xmin-=(nx-map_nx)/2*loc->dx;
	ymin-=(ny-map_ny)/2*loc->dy;
	map_nx=nx;
	map_ny=ny;
    }
    ((loc_t*)loc)->map=mapnew(map_nx, map_ny, loc->dx, loc->dy);
    loc->map->iac=loc->iac;
    loc->map->ox=xmin;
    loc->map->oy=ymin;
    dmat*  pmap=(dmat*)loc->map;
    const real *locx=loc->locx;
    const real *locy=loc->locy;
    for(long iloc=0; iloc<loc->nloc; iloc++){
	int ix=(int)round((locx[iloc]-xmin)*dx_in1);
	int iy=(int)round((locy[iloc]-ymin)*dy_in1);
	P(pmap,ix,iy)=iloc+1;/*start from 1. */
    }
    if(LOC_MAP_EXTEND && loc->nloc<map_nx*map_ny){
	/*
	  For cubic spline interpolation around the edge of the grid, we need
	  "fake" actuators that take the value of its neighbor to avoid edge
	  roll off by using zero for non-existing actuators. This is achieved by
	  filling the vacant locations in the map to loc mapping with negative
	  indices of active neighbor. The neighbors are selected according to 1)
	  whether this is already a ghost actuator, 2) how far away is the
	  neighbor, 3) how far away the neighbor is from the center.

	  The ghost actuator mapping is tweaked to properly handle the boundary
	  region where there are only three actuators. We treat these regions as active.

	  A remain question is which neighbor's value to assume when there are
	  two neighbors that are equally close by.
	 */
	lmat *level=lnew(map_nx, map_ny);
	for(int iy=0; iy<map_ny; iy++){
	    for(int ix=0; ix<map_nx; ix++){
		if(P(pmap,ix,iy)){
		    P(level,ix,iy)=0;//level of existence.
		}
	    }
	}
	/*we do the interpolation with incremental level, for max of npad levels*/
	for(int cur_level=0; cur_level<=npad; cur_level++){
	    int num_found;
	    do{
		num_found=0;
		for(int iy=0; iy<map_ny; iy++){
		    for(int ix=0; ix<map_nx; ix++){
			if(!P(pmap,ix,iy)){
			    int count=0;//count how many neighbors of lowest level
			    /*choose the minimum level of interpolation*/
			    int min_level=INT_MAX, min_jx=0, min_jy=0, min_sep=INT_MAX;
			    real min_rad=INFINITY;
			    for(int jy=MAX(0, iy-1); jy<MIN(map_ny, iy+2); jy++){
				for(int jx=MAX(0, ix-1); jx<MIN(map_nx, ix+2); jx++){
				    int sep=(abs(iy-jy)+abs(ix-jx));
				    int iphi;
				    if((sep==1 || sep==2) && (iphi=fabs(P(pmap,jx,jy)))){//edge
					iphi--;
					real rad=locx[iphi]*locx[iphi]+locy[iphi]*locy[iphi];//dist^2 to center
					if(P(level,jx,jy)<min_level){
					    min_level=P(level,jx,jy);
					    min_rad=rad;
					    min_sep=sep;
					    min_jy=jy;
					    min_jx=jx;
					    count=1;//reset to one
					}else if(P(level,jx,jy)==min_level){
					    count++;
					    if(sep<min_sep){
						min_sep=sep;
						min_jy=jy;
						min_jx=jx;
					    }else if(sep==min_sep && rad<min_rad){
						min_rad=rad;
						min_jy=jy;
						min_jx=jx;
					    }
					}
				    }
				}
			    }
			    if((min_level==cur_level && count==3) || min_level<cur_level){
				//if level satisfy or even if three neighbor with higher level.
				P(pmap,ix,iy)=-fabs(P(pmap,min_jx,min_jy));
				P(level,ix,iy)=min_level+1;
				num_found++;
			    }
			}
		    }
		}
	    }while(num_found);
	}
	lfree(level);
    }
    UNLOCK(maplock);
}
/*
  Embed in into dest according to map defined in loc->map. Arbitraty map cannot
  be used here for mapping.
*/
void loc_embed(map_t *dest, const loc_t *loc, const real *in){
    map_t *map=loc->map;
    if(!map){
	loc_create_map((loc_t*)loc);
	map=loc->map;
    }
    if(dest->nx!=map->nx || dest->ny!=map->ny){
	error("dest and map doesn't agree\n");
    }
    const real *pin=in-1;//iphi count from 1
    for(long i=0; i<map->nx*map->ny; i++){
	long iphi=fabs(map->p[i]);
	if(iphi){
	    dest->p[i]=pin[iphi];
	}else{
	    dest->p[i]=NAN;
	}
    }
}
/*
  A convenient wrapper for loc_embed to be called by matlab or python
*/
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

/*
  Embed in into dest according to map defined in loc->map. Arbitraty map cannot
  be used here for mapping.
*/
void loc_embed_add(map_t *dest, const loc_t *loc, const real *in){
    map_t *map=loc->map;
    if(!map){
	loc_create_map((loc_t*)loc);
	map=loc->map;
    }
    if(dest->nx!=map->nx || dest->ny!=map->ny){
	error("dest and map doesn't agree\n");
    }
    const real *pin=in-1;//iphi count from 1
    for(long i=0; i<map->nx*map->ny; i++){
	long iphi=fabs(map->p[i]);
	if(iphi){
	    dest->p[i]+=pin[iphi];
	}
    }
}
/*
  Embed in into dest according to map defined in loc->map for cells. Arbitraty map cannot
  be used here for mapping.
*/
void loc_embed_cell(dcell **pdest, const loc_t *loc, const dcell *in){
    map_t *map=loc->map;
    if(!map){
	loc_create_map((loc_t*)loc);
	map=loc->map;
    }
    if(!*pdest){
	*pdest=dcellnew(map->nx, map->ny);
    }
    dcell *dest=*pdest;
    if(dest->nx!=map->nx || dest->ny!=map->ny){
	error("dest and map doesn't agree\n");
    }
    for(long i=0; i<map->nx*map->ny; i++){
	long iphi=fabs(map->p[i]);
	if(iphi){
	    dest->p[i]=dref(P(in, iphi-1));
	}else{
	    dest->p[i]=0;
	}
    }
}
/*
  Embed in into dest according to map defined in loc->map. Arbitraty map cannot
  be used here for mapping.
*/
void loc_extract(dmat *dest, const loc_t *loc, map_t *in){
    map_t *map=loc->map;
    if(!map){
	error("map is null\n");
    }
    if(in->nx!=map->nx || in->ny!=map->ny){
	error("in and map doesn't agree: in is %ldx%ld, map is %ldx%ld\n",
	      in->nx, in->ny, map->nx, map->ny);
    }
    real *restrict pdest=dest->p-1;//iphi count from 1
    for(long i=0; i<map->nx*map->ny; i++){
	long iphi=map->p[i];
	if(iphi>0){
	    pdest[iphi]=in->p[i];
	}
    }
}
/**
   Convert a map to a loc that collects all positive entries. */
loc_t* map2loc(map_t *map, real thres){
    const real dx=map->dx;
    const real dy=map->dy;
    const real ox=map->ox;
    const real oy=map->oy;
    const long nx=map->nx;
    const long ny=map->ny;
    if(thres){
	thres*=dmax((dmat*)map);
    }
    long ix,iy;
    loc_t *loc=locnew(nx*ny, dx, dy);
    long count=0;
    for(iy=0; iy<ny; iy++){
	for(ix=0; ix<nx; ix++){
	    if(P(map,ix,iy)>thres){
		loc->locx[count]=ix*dx+ox;
		loc->locy[count]=iy*dy+oy;
		count++;
	    }
	}
    }
    if(!count){
	warning("map2loc: there are 0 points.\n");
	writebin(map, "map2loc");
	print_backtrace();
    }
    loc->locx=myrealloc(loc->locx,count,real);
    loc->locy=myrealloc(loc->locy,count,real);
    loc->nloc=count;
    loc->iac=map->iac;
    loc->ht=map->h;
    return loc;
}
/**
   Create 1 dimensional loc with given vector.
*/
loc_t *mk1dloc_vec(real *x, long nx){
    real dx=fabs((x[nx-1]-x[0])/(nx-1));
    loc_t *loc=locnew(nx, dx, dx);
    memcpy(loc->locx, x, sizeof(real)*nx);
    return loc;
}
/**
   Create 1 dimensional loc with origin at x0, sampling of dx, and nx numbers.
*/
loc_t *mk1dloc(real x0, real dx, long nx){
    loc_t *loc=locnew(nx, dx, dx);
    for(long ix=0; ix<nx; ix++){
	loc->locx[ix]=x0+dx*ix;
    }
    return loc;
}
/**
   Create a loc array that covers the map_t. Notice that it is different from
   map2loc which only covers valid (value>0) regions.
*/
loc_t *mksqloc_map(const map_t*map){
    return mksqloc(map->nx, map->ny, map->dx, map->dy, map->ox, map->oy);
}
/**
   Create a loc array of size nx*ny at sampling dx with origin at (nx/2, ny/2)
*/
loc_t *mksqloc_auto(long nx, long ny, real dx, real dy){
    return mksqloc(nx,ny,dx,dy,-nx/2*dx,-ny/2*dy);
}
/**
   Create a loc array contains coorindates in a square map of size nx*ny, with
   sampling dx, and at origin ox,oy */
loc_t *mksqloc(long nx, long ny, real dx, real dy, real ox, real oy){
    loc_t *loc=locnew(nx*ny, dx, dy);
    for(long iy=0; iy<ny; iy++){
	real y=iy*dy+oy;
	for(long ix=0; ix<nx; ix++){
	    loc->locx[ix+iy*nx]=ix*dx+ox;
	    loc->locy[ix+iy*nx]=y;
	}
    }
    return loc;
}

/**
   Create a loc array within diameter D, and inner diameter Din, with spacing dx.
 */
loc_t *mkannloc(real D, real Din, real dx, real thres){
    long nx=D/dx+1;
    map_t *xy=mapnew(nx, nx, dx, dx);
    mapcircle(xy, D/2, 1);
    if(Din>0 && Din<D){
	mapcircle(xy, Din/2, -1);
    }
    loc_t *loc=map2loc(xy, thres);
    mapfree(xy);
    return loc;
}
/**
   A convenient wrapp for dcircle.
 */
dmat *mkcirmap(long nx, long ny, real cx, real cy, real r){
    dmat *map=dnew(nx, ny);
    dcircle(map, cx, cy, 1,1, r, 1);
    return map;
}
/**
   Estimate the diameter of loc
*/
real loc_diam(const loc_t *loc){
    real R2max=0;
    for(long i=0; i<loc->nloc; i++){
	real R2=loc->locx[i]*loc->locx[i]+loc->locy[i]*loc->locy[i];
	if(R2max<R2){
	    R2max=R2;
	}
    }
    return sqrt(R2max)*2;
}
/**
   Find the point that closes to origin (0,0) in loc. Useful in single point
   piston constraint in creating reconstructor.
*/
int loccenter(const loc_t *loc){
    int ipix,jpix=0;
    real r2,r2min=INFINITY;
    
    for(ipix=0; ipix<loc->nloc; ipix++){
	r2=pow(loc->locx[ipix],2)+pow(loc->locy[ipix],2);
	if(r2<r2min){
	    jpix=ipix;
	    r2min=r2;
	}
    }
    return jpix;
}

/**
   assemble modes of piston/tip/tilt into Nx3 matrix M, and compute their
   inverse covariance matrix.  mcc=M'*(amp.*M) */
dmat *loc_mcc_ptt(const loc_t *loc, const real *amp){
    const int nmod=3;
    real *mod[nmod];
    dmat *mcc=dnew(nmod,nmod);
    mod[0]=NULL;
    mod[1]=loc->locx;
    mod[2]=loc->locy;
    for(int jmod=0; jmod<nmod; jmod++){
	for(int imod=jmod; imod<nmod; imod++){
	    real tmp=dvecdot(mod[imod],mod[jmod],amp,loc->nloc);
	    P(mcc,imod,jmod)=P(mcc,jmod,imod)=tmp;
	}
    }
    return mcc;
}
/**
   same as loc_mcc_ptt, except using pts instead of loc.
   this is used to build imcc for TT/F powfs ztilt imcc
*/
dcell *pts_mcc_ptt(const pts_t *pts, const real *amp){
    const int nmod=3;
    const int nsa=pts->nsa;
    dcell *mcc=dcellnew_same(nsa,1,nmod,nmod);
    for(int isa=0; isa<nsa; isa++){
	const real origy=pts->origy[isa];
	const real origx=pts->origx[isa];
	const real dx=pts->dx;
	const real dy=pts->dy;
	const real *ampi=amp+pts->nx*pts->nx*isa;
	dmat *ATA=mcc->p[isa];
	real a00=0,a01=0,a02=0,a11=0,a12=0,a22=0;
	for(int iy=0; iy<pts->nx; iy++){
	    real y=iy*dy+origy;
	    const real *ampx=ampi+iy*pts->nx;
	    for(int ix=0; ix<pts->nx; ix++){
		real x=ix*dx+origx;
		real a=ampx[ix];
		a00+=a;
		a01+=a*x;
		a02+=a*y;
		a11+=a*x*x;
		a12+=a*x*y;
		a22+=a*y*y;
	    }
	}
	P(ATA,0,0)=a00;
	P(ATA,1,1)=a11;
	P(ATA,2,2)=a22;
	P(ATA,0,1)=P(ATA,1,0)=a01;
	P(ATA,0,2)=P(ATA,2,0)=a02;
	P(ATA,1,2)=P(ATA,2,1)=a12;
    }
    return mcc;
}

/**
   evaluate piston/tip-tilt/ removed wavefront error.
   output coeffout in unit of radian like units.
*/
void loc_calc_ptt(real *rmsout, real *coeffout,
		  const loc_t *loc, const real ipcc, 
		  const dmat *imcc, const real *amp, const real *opd){
    assert(imcc->nx==imcc->ny && imcc->nx==3);
    const long nloc=loc->nloc;
    const real *restrict locx=loc->locx;
    const real *restrict locy=loc->locy;

    real tot=0;
    real coeff[3]={0,0,0};
    if(amp){
	for(long iloc=0; iloc<nloc; iloc++){
	    const real junk=opd[iloc]*amp[iloc];
	    coeff[0]+=junk;
	    coeff[1]+=junk*locx[iloc];
	    coeff[2]+=junk*locy[iloc];
	    tot+=junk*opd[iloc];
	}
    }else{
	for(long iloc=0; iloc<nloc; iloc++){
	    const real junk=opd[iloc];
	    coeff[0]+=junk;
	    coeff[1]+=junk*locx[iloc];
	    coeff[2]+=junk*locy[iloc];
	    tot+=junk*opd[iloc];
	}
    }
    if(coeffout){
	dmulvec3(coeffout, imcc, coeff);
    }
    if(rmsout){
	real pis=ipcc*coeff[0]*coeff[0];/*piston mode variance */
	real ptt=dwdot3(coeff, imcc, coeff);/*p/t/t mode variance. */
	rmsout[0]=tot-pis;/*PR */
	rmsout[1]=ptt-pis;/*TT */
	rmsout[2]=tot-ptt;/*PTTR */
    }
}
/**
   Calculate variance of OPD in modes defined in mod.
   
   Notice the difference in rmsout for loc_calc_ptt and this loc_calc_mode.
   in loc_calc_ptt, out[0] is PR, out[1] is TT, out[2] is PTTR
   in loc_calc_mod, out[0] is PR, out[1] is PTR, out[2] is PTTR, etc
       
   The mod is already orth-noramlized so that we have the
   following simple forms.

   coeffout is in unit of zernike!
*/
void loc_calc_mod(real *rmsout, real *coeffout,const dmat *mod,
		  const real *amp, real *opd){
  
    const int nmod=mod->ny;
    const int nloc=mod->nx;
    real tot=0;
    real val[nmod];
    memset(val, 0, sizeof(real)*nmod);
    for(long iloc=0; iloc<nloc; iloc++){
	real junk=opd[iloc]*amp[iloc];
	tot+=opd[iloc]*junk;
	for(long imod=0; imod<nmod; imod++){
	    val[imod]+=P(mod,iloc,imod)*junk;
	}
    }
    for(long imod=0; imod<nmod; imod++){
	tot-=val[imod]*val[imod];
	rmsout[imod]=tot;
    }
    if(coeffout){
	for(long imod=0; imod<nmod; imod++){
	    coeffout[imod]=val[imod];
	}
    }
}

/**
   Remove Piston/Tip/Tilt (in radian) from OPD
*/
void loc_remove_ptt(real *opd, const real *ptt, const loc_t *loc){
    real ptt1[3];
    ptt1[0]=-ptt[0];
    ptt1[1]=-ptt[1];
    ptt1[2]=-ptt[2];
    loc_add_ptt(opd, ptt1, loc);
}

/**
   Add Piston/Tip/Tilt from OPD
*/
void loc_add_ptt(real *opd, const real *ptt, const loc_t *loc){
    const long nloc=loc->nloc;
    const real *restrict locx=loc->locx;
    const real *restrict locy=loc->locy;
    for(long iloc=0; iloc<nloc; iloc++){
	opd[iloc]+=ptt[0]+ptt[1]*locx[iloc]+ptt[2]*locy[iloc];
    }
}
/**
   Compute zernike best fit for all subapertures. add result to out.  returns
   radians not zernike modes of tip/tilt. Used in wfsgrad */
void pts_ztilt(dmat **out, const pts_t *pts, const dcell *imcc,
	       const real *amp, const real *opd){
    const int nsa=pts->nsa;
    assert(imcc->nx==nsa && imcc->ny==1);
    if(!*out) *out=dnew(nsa*2,1);
    real *res=(*out)->p;
    for(int isa=0; isa<nsa; isa++){
	const real origy=pts->origy[isa];
	const real origx=pts->origx[isa];
	const real dx=pts->dx;
	const real dy=pts->dy;
	const real *ampi=amp+pts->nx*pts->nx*isa;
	const real *opdi=opd+pts->nx*pts->nx*isa;
	assert(imcc->p[isa]->nx==3 && imcc->p[isa]->ny==3);
        real coeff[3]={0,0,0};
	real a0=0,a1=0,a2=0;
	for(int iy=0; iy<pts->nx; iy++){
	    real y=iy*dy+origy;
	    const real *ampx=ampi+iy*pts->nx;
	    const real *opdx=opdi+iy*pts->nx;
	    for(int ix=0; ix<pts->nx; ix++){
		real x=ix*dx+origx;
		real tmp=ampx[ix]*opdx[ix];
		a0+=tmp;
		a1+=tmp*x;
		a2+=tmp*y;
	    }
	}
	coeff[0]=a0;
	coeff[1]=a1;
	coeff[2]=a2;
	real outp[3];
	dmulvec3(outp, imcc->p[isa], coeff);
	/*
	  2010-07-19: Was =, modified to += to conform to the convention.
	*/
	res[isa]+=outp[1];
	res[isa+nsa]+=outp[2];
    }
}
/**
   Gather information about the starting of each column in loc.
*/
void loc_create_stat_do(loc_t *loc){
    locstat_t *locstat=mycalloc(1,locstat_t);
    loc->stat=locstat;
    const real *locx=loc->locx;
    const real *locy=loc->locy;
    int nloc=loc->nloc;
    real dx=locstat->dx=loc->dx;
    real dy=locstat->dy=loc->dy;
    int ncolmax=(int)round((locy[nloc-1]-locy[0])/dy)+2;
    locstat->cols=mymalloc(ncolmax,locstatcol_t);
    int colcount=0;
    /*do first column separately. */
    int iloc=0;
    locstat->cols[colcount].pos=iloc;
    locstat->cols[colcount].xstart=locx[iloc];
    locstat->cols[colcount].ystart=locy[iloc];
    locstat->ymin=locstat->cols[colcount].ystart;
    locstat->xmin=locstat->cols[colcount].xstart;
    real xmax=locstat->cols[colcount].xstart;
    const real dythres=fabs(dy)*0.1;
    const real dxthres=fabs(dx)*0.1;
    colcount++;
    for(iloc=1; iloc<loc->nloc; iloc++){
	if(fabs(locy[iloc]-locy[iloc-1])>dythres /*a new column starts */
	   || fabs(locx[iloc]-locx[iloc-1]-dx)>dxthres){
	    locstat->cols[colcount].pos=iloc;
	    locstat->cols[colcount].xstart=locx[iloc];
	    locstat->cols[colcount].ystart=locy[iloc];
	    if(xmax < locx[iloc-1]){
		xmax = locx[iloc-1];
	    }
	    if(locstat->xmin>locstat->cols[colcount].xstart){
		locstat->xmin=locstat->cols[colcount].xstart;
	    }
	    colcount++;
	    if(colcount>=ncolmax){/*expand the memory. */
		ncolmax*=2;
		locstat->cols=myrealloc(locstat->cols, ncolmax,locstatcol_t);
	    }
	}
    }
    locstat->ncol=colcount;
    locstat->cols[colcount].pos=loc->nloc;
    locstat->cols[colcount].xstart=locx[loc->nloc-1];
    locstat->cols[colcount].ystart=locy[loc->nloc-1];
    if(xmax < locx[loc->nloc-1]){
	xmax = locx[loc->nloc-1];
    }
    locstat->cols=myrealloc(locstat->cols,(locstat->ncol+1),locstatcol_t);
    locstat->ny=(long)round((locy[loc->nloc-1]-locstat->ymin)/dy)+1;
    locstat->nx=(long)round((xmax-locstat->xmin)/dx)+1;
}
/**
   Create a gray pixel circular map in phi using coordinates defined in loc, center
   defined using cx, cy, radius of r, and value of val */
void loccircle(real *phi,loc_t *loc,real cx,real cy,real r,real val){
    /*cx,cy,r are in unit of true unit, as in loc */
    real dx=loc->dx;
    real dy=loc->dy;
    int nres=10;
    real cres=(nres-1.)/2.;
    real res=1./nres;
    real dxres=res*dx;
    real dyres=res*dy;
    real res2=res*res;
    real r2=r*r;
    real r2l=(r-dx*1.5)*(r-dy*1.5);
    real r2u=(r+dx*1.5)*(r+dy*1.5);
    real *locx=loc->locx;
    real *locy=loc->locy;
    real iix,iiy;
    for(int iloc=0; iloc<loc->nloc; iloc++){
	real x=locx[iloc];
	real y=locy[iloc];
	real r2r=pow(x-cx,2)+pow(y-cy,2);
	if(r2r<r2l) 
	    phi[iloc]+=val;
	else if(r2r<r2u){
	    long tot=0;
	    for(int jres=0; jres<nres; jres++){
		iiy=y+(jres-cres)*dyres;
		real rr2y=(iiy-cy)*(iiy-cy);
		for(int ires=0; ires<nres; ires++){
		    iix=x+(ires-cres)*dxres;
		    real rr2r=(iix-cx)*(iix-cx)+rr2y;
		    if(rr2r<=r2)
			tot++;
		}
	    }
	    phi[iloc]+=(real)tot*res2*val;
	}
    }
}
/**
   Create a gray pixel annular map in phi using coordinates defined in loc,
   center defined using cx, cy, radius of r, and value of val
*/
void locannular(real *phi,loc_t *loc,real cx,real cy,real r,real rin,real val){
    loccircle(phi,loc,cx,cy,r,val);
    if(rin>EPS){
	loccircle(phi,loc,cx,cy,rin,-val);
    }
}
/**
   Create a hard annular mask in phi.
*/
void locannularmask(real *phi,loc_t *loc,real cx,real cy,real r,real rin){
    /*apply the hard pupil mask of aper.d, using loc. not locm */
    /* 2011-07-13: changed from r^2 to (r+0.5*dx)^2*/
    real rr2max=pow(r+0.25*(loc->dx+loc->dy),2);
    real rr2min=MIN(rin*rin, pow(rin-0.25*(loc->dx+loc->dy),2));
    for(long iloc=0; iloc<loc->nloc; iloc++){
	real r2=pow(loc->locx[iloc]-cx,2)+pow(loc->locy[iloc]-cy,2);
	if(r2<rr2min || r2>rr2max){
	    phi[iloc]=0;
	}
    }
}
/**
   Create a gray pixel elliptical map in phi using coordinates defined in loc,
   center defined using cx, cy, radii of rx, ry, and value of val */
void locellipse(real *phi,loc_t *loc,real cx,real cy,
		real rx,real ry,real val){
    /*cx,cy,r are in unit of true unit, as in loc */
    real dx=loc->dx;
    real dy=loc->dy;
    int nres=10;
    real cres=(nres-1.)/2.;
    real res=1./nres;
    real dxres=res*dx;
    real dyres=res*dy;
    real res2=res*res;
    real rx1=1./rx;
    real ry1=1./ry;
    real rx12=rx1*rx1;
    real ry12=ry1*ry1;
    real r2=2.;
    real r2l=(rx-dx)*(rx-dx)*rx12+(ry-dy)*(ry-dy)*ry12;
    real r2u=(rx+dx)*(rx+dx)*rx12+(ry+dy)*(ry+dy)*ry12;
    real *locx=loc->locx;
    real *locy=loc->locy;
    real iix,iiy;
    for(int iloc=0; iloc<loc->nloc; iloc++){
	real x=locx[iloc];
	real y=locy[iloc];
	real r2r=pow(x-cx,2)*rx12+pow(y-cy,2)*ry12;
	if(r2r<r2l) 
	    phi[iloc]+=val;
	else if(r2r<r2u){
	    long tot=0;
	    for(int jres=0; jres<nres; jres++){
		iiy=y+(jres-cres)*dyres;
		real rr2y=(iiy-cy)*(iiy-cy)*ry12;
		for(int ires=0; ires<nres; ires++){
		    iix=x+(ires-cres)*dxres;
		    real rr2r=(iix-cx)*(iix-cx)*rx12+rr2y;
		    if(rr2r<r2)
			tot++;
		}
	    }
	    phi[iloc]+=tot*res2*val;
	}
    }
}
/**
   Remove points in loc that have zero value in amp and modify amp in the same
   time. Keep each row continuous if cont==1. Return in skipout the index of
   skipped points if skipout is not NULL.  */

void loc_reduce(loc_t *loc, dmat *amp, real thres, int cont, int **skipout){
    if(thres<=0) thres=EPS;
    int redo_stat=loc->stat?1:0;
    int nloc=loc->nloc; 
    int *skip=mycalloc(nloc,int);
    if(cont){/*make sure loc is continuous. */
	loc_create_stat(loc);
	locstat_t *locstat=loc->stat;
	int ncol=locstat->ncol;
	for(int icol=0; icol<ncol; icol++){
	    int pstart=locstat->cols[icol].pos;
	    int pend=locstat->cols[icol+1].pos;
	    for(int pos=pstart; pos<pend; pos++){
		if(amp->p[pos]<thres){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	    for(int pos=pend-1; pos>pstart-1; pos--){
		if(amp->p[pos]<thres){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	}
    }else{
	for(int iloc=0; iloc<nloc; iloc++){
	    if(amp->p[iloc]<thres){
		skip[iloc]=1;
	    }
	}
    }
    int count=0;
    for(int iloc=0; iloc<nloc; iloc++){
	loc->locx[count]=loc->locx[iloc];
	loc->locy[count]=loc->locy[iloc];
	amp->p[count]=amp->p[iloc];
	if(!skip[iloc]) count++;
    }
    locresize(loc, count);
    dresize(amp, count, 1);
    if(redo_stat){
	loc_create_stat(loc);
    }
    if(skipout){
	*skipout=skip;
    }else{
	free(skip);
    }
}
/**
   Remove uncoupled points in loc and modify spc in the same time. debugged on 2009-12-20. Not used often. using
   a dspcell, to compute the coupling which is modified accordingly.  */
void loc_reduce_spcell(loc_t *loc, dspcell *spc, int dim, int cont){
    int nloc=loc->nloc;
    dmat *sum=NULL;
    for(int isp=0; isp<spc->nx*spc->ny; isp++){
	if(spc->p[isp]){
	    dmat* sum0=dspsumabs(spc->p[isp],3-dim);
	    dadd(&sum,1,sum0,1);
	    dfree(sum0);
	}
    }
    if(!sum || sum->nx*sum->ny!=loc->nloc){
	error("Mismatch\n");
    }
    int *skip;
    loc_reduce(loc, sum, EPS, cont, &skip);
    dfree(sum);
    int count;
    if(dim==1){/*opd(loc)=sp*opd(locin); */
	count=0;
	int *map=mycalloc(nloc,int);
	for(int iloc=0; iloc<nloc; iloc++){
	    map[iloc]=count;
	    if(!skip[iloc]) count++;
	}
	for(int isp=0; isp<spc->nx*spc->ny;isp++){
	    dsp *sp=spc->p[isp];
	    if(!sp) continue;
	    sp->nx=count;
	    for(int iz=0; iz<sp->nzmax; iz++){
		sp->i[iz]=map[sp->i[iz]];
	    }
	}
	free(map);
    }else if(dim==2){
	for(int isp=0; isp<spc->nx*spc->ny;isp++){
	    dsp *sp=spc->p[isp];
	    if(!sp) continue;
	    count=0;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(!skip[iloc]){
		    sp->p[count]=sp->p[iloc];
		    count++;
		}
	    }
	    sp->p[count]=sp->p[nloc];
	    sp->ny=count;
	    sp->p=myrealloc(sp->p,(count+1),spint);
	}
    }
    free(skip);
}

/**
   Remove uncoupled points in loc and modify sp in the same time. debugged on 2009-12-20.  use single sparse
   matrix, which is modified accordingly.  */
void loc_reduce_sp(loc_t *loc, dsp *sp, int dim, int cont){
    int nloc=loc->nloc;
    if((dim==1 && nloc!=sp->nx) || (dim==2 && nloc!=sp->ny) || dim<0 || dim>2)
	error("Mismatch dimension\n");
    dmat* sum=dspsumabs(sp,3-dim);
    int *skip;
    loc_reduce(loc, sum, EPS, cont, &skip);
    dfree(sum);
    int count=0;
    if(dim==1){/*opd(loc)=sp*opd(locin); */
	int *map=mycalloc(nloc,int);
	for(int iloc=0; iloc<nloc; iloc++){
	    map[iloc]=count;
	    if(!skip[iloc]) count++;
	}
	sp->nx=count;
	for(int iz=0; iz<sp->nzmax; iz++){
	    sp->i[iz]=map[sp->i[iz]];
	}
	free(map);
    }else if(dim==2){
	for(int iloc=0; iloc<nloc; iloc++){
	    if(!skip[iloc]){
		sp->p[count]=sp->p[iloc];
		count++;
	    }
	}
	sp->p[count]=sp->p[nloc];
	sp->ny=count;
	sp->p=myrealloc(sp->p,(count+1),spint);
    }
    free(skip);
}

/**
   Add val amount of focus to opd. The unit is in radian like.
   Piston is removed.
*/
void loc_add_focus(real *opd, loc_t *loc, real val){
    if(fabs(val)<1.e-15) return;
    const real *restrict locx=loc->locx;
    const real *restrict locy=loc->locy;
    real piston=0;
    for(long iloc=0; iloc<loc->nloc; iloc++){
	real tmp=(locx[iloc]*locx[iloc]+locy[iloc]*locy[iloc])*val;
	opd[iloc]+=tmp;
	piston+=tmp;
    }
    piston=-piston/loc->nloc;
    for(long iloc=0; iloc<loc->nloc; iloc++){
	opd[iloc]+=piston;
    }
}
/**
   Create a dmat containing locx and locy in two columns
*/
dmat *loc2mat(loc_t *loc,int piston){
    dmat *out=NULL;
    real *ptt=NULL;
    if(piston){
	out=dnew(loc->nloc,3);
	for(long i=0; i<loc->nloc; i++){
	    out->p[i]=1;
	}
	ptt=out->p+loc->nloc;
    }else{
	out=dnew(loc->nloc,2);
	ptt=out->p;
    }
    memcpy(ptt,loc->locx,sizeof(real)*loc->nloc);
    memcpy(ptt+loc->nloc,loc->locy,sizeof(real)*loc->nloc);
    return out;
}

/**
   Convert pts_t to loc_t. Each entry in pts recors the lower
   left point in each subaperture. 
*/
loc_t *pts2loc(pts_t *pts){
    long nsa  = pts->nsa;
    long nx   = pts->nx;
    long nxsa = nx*nx;
    real dx = pts->dx;
    real dy = pts->dy;
    if(dy==0){
	dy=dx;
    }
    loc_t *loc = locnew(nsa*nxsa, dx, dy);
    for(int isa=0; isa<nsa; isa++){
	const real origx = pts->origx[isa];
	const real origy = pts->origy[isa];
	real *locx = loc->locx+isa*nxsa;
	real *locy = loc->locy+isa*nxsa;
	for(int jj=0; jj<nx; jj++){
	    const real yy=origy+(real)jj*dy;
	    for(int ii=0; ii<nx; ii++){
		locx[jj*nx+ii]=origx+(real)ii*dx;
		locy[jj*nx+ii]=yy;
	    }
	}
    }
    return loc;
}

/**
   Rotate the coordinates by theta (radian) CCW.
*/
void locrot(loc_t *loc, const real theta){
    const real ctheta=cos(theta);
    const real stheta=sin(theta);

    real *x=loc->locx;
    real *y=loc->locy;
    for(int i=0; i<loc->nloc; i++){
	real tmp=x[i]*ctheta-y[i]*stheta;
	y[i]=x[i]*stheta+y[i]*ctheta;
	x[i]=tmp;
    }
}
/**
   Determine the angle of rotation from loc1 to loc2. CCW is positive.
 */
real loc_angle(const loc_t *loc1, const loc_t *loc2){
    real mx1, mx2, my1, my2;
    locmean(&mx1, &my1, loc1);
    locmean(&mx2, &my2, loc2);
    real sina=0; 
    long count=0;
    const real thres=(loc1->dx+loc1->dy)*0.5;
    for(long i=0; i<loc1->nloc; i++){
	real x1=loc1->locx[i]-mx1;
	real y1=loc1->locy[i]-my1;
	real x2=loc2->locx[i]-mx2;
	real y2=loc2->locy[i]-my2;
	real lcross=x1*y2-x2*y1;
	real lamp=sqrt((x1*x1+y1*y1)*(x2*x2+y2*y2));
	if(lamp>thres){
	    sina+=lcross/lamp;
	    count++;
	}
    }
    sina/=count;
    return asin(sina);
}
/**
   Stretch the coordinate by frac along theta (radian) CCW.
*/
void locstretch(loc_t *loc, const real theta, const real frac){
    locrot(loc, -theta);
    for(int i=0; i<loc->nloc; i++){
	loc->locx[i]*=frac;
    }
    locrot(loc, theta);
}

/**
   duplicate a loc_t object; */
loc_t *locdup(loc_t *loc){
    loc_t *res=locnew(loc->nloc, loc->dx, loc->dy);
    memcpy(res->locx,loc->locx,sizeof(real)*loc->nloc);
    memcpy(res->locy,loc->locy,sizeof(real)*loc->nloc);
    res->iac=loc->iac;
    res->ht=loc->ht;
    return res;
}
/**
   Compute average of locx, locy*/
void locmean(real *xm, real *ym, const loc_t *loc){
    real x=0, y=0;
    for(long i=0; i<loc->nloc; i++){
	x+=loc->locx[i];
	y+=loc->locy[i];
    }
    *xm=x/loc->nloc;
    *ym=y/loc->nloc;
}

/**Parse string representation of polynominal to array representation.*/
static dmat *parse_poly(const char *_ps){
    char *ps=(char *)_ps;
    int ncx=5;
    dmat *cx=dnew(3, ncx);
    int icx=0;
    char *endptr;
    while(ps[0] && ps[0]!=';'){
	P(cx,0,icx)=strtod(ps, &endptr);//coefficient
	if(ps==endptr){
	    if(ps[0]=='-'){
		ps++;
		P(cx,0,icx)=-1;
	    }else if(ps[0]=='+'){
		ps++;
		P(cx,0,icx)=1;
	    }
	    if(ps[0]=='x' || ps[0]=='y'){
		if(P(cx,0,icx)==0) P(cx,0,icx)=1;
	    }else{
		error("Unable to parse (%s). ps=(%s)\n", _ps, ps);
	    }
	}else{
	    ps=endptr;
	}
	P(cx,1,icx)=P(cx,2,icx)=0;
	while(ps[0]==' ') ps++;
	if(ps[0]=='*') ps++;
	while(ps[0]==' ') ps++;
        //parse each term
	while(ps[0]!=0 && ps[0]!='+' && ps[0]!='-' && ps[0]!=';'){
	    if(ps[0]=='x' || ps[0]=='y'){
		int iy=1;
		if(ps[0]=='y'){
		    iy=2;
		}
		ps++;
		if(ps[0]=='^'){
		    ps++;
		    P(cx,iy,icx)+=strtol(ps, &endptr, 10);
		    if(ps==endptr){
			error("Unable to parse %s\n", _ps);
		    }
		    ps=endptr;
		}else{
		    P(cx,iy,icx)++;
		}
	    }else{
		P(cx,0,icx)*=strtod(ps, &endptr);
		if(ps==endptr){
		    error("Unable to parse %s\n", _ps);
		}
		ps=endptr;
	    }
	    while(ps[0]==' ') ps++;
	}
	icx++;
	if(ps[0]==';' || ps[0]==0){
	    break;
	}else if(ps[0]=='+' || ps[0]=='-'){
	    if(ps[0]=='+') ps++;
	    if(icx>=ncx){
		ncx*=2;
		dresize(cx, 3, ncx);
	    }
	}else{
	    error("Unable to parse (%s)\n", ps);
	}
    }
    ncx=icx;
    dresize(cx, 3, ncx);
    return cx;
    
}
static loc_t *loctransform_do(const loc_t *loc, const dmat *cx, const dmat *cy){
    loc_t *locm=locnew(loc->nloc, loc->dx, loc->dy);
    real *restrict xm=locm->locx;
    real *restrict ym=locm->locy;
    const real *restrict x=loc->locx;
    const real *restrict y=loc->locy;
    int np=0;
    int nonint=0;
    for(int ic=0; ic<cx->ny; ic++){
	int ocx=(int)round(P(cx,1,ic));
	int ocy=(int)round(P(cx,2,ic));
	if((P(cx,1,ic)-ocx)>EPS || (P(cx,2,ic)-ocy)>EPS){//not integer
	    nonint=1;
	}
	if(ocx>np) np=ocx;
	if(ocy>np) np=ocy;
    }
    for(int ic=0; ic<cy->ny; ic++){
	int ocx=(int)round(P(cy,1,ic));
	int ocy=(int)round(P(cy,2,ic));
	if((P(cy,1,ic)-ocx)>EPS || (P(cy,2,ic)-ocy)>EPS){//not integer
	    nonint=1;
	}
	if(ocx>np) np=ocx;
	if(ocy>np) np=ocy;
    }
    np++; //max order of x or y +1
    for(long iloc=0; iloc<loc->nloc; iloc++){
	if(nonint){
	    for(long ic=0; ic<cx->ny; ic++){
		xm[iloc]+=P(cx,0,ic)*pow(x[iloc],P(cx,1,ic))*pow(y[iloc],P(cx,2,ic));
	    }
		
	    for(long ic=0; ic<cy->ny; ic++){
		ym[iloc]+=P(cy,0,ic)*pow(x[iloc],P(cy,1,ic))*pow(y[iloc],P(cy,2,ic));
	    }
	}else{/*faster method for integer powers (>10x speed up). */
	    real xp[np], yp[np];
	    xp[0]=1; yp[0]=1;
	    for(long ip=1; ip<np; ip++){
		xp[ip]=xp[ip-1]*x[iloc];
		yp[ip]=yp[ip-1]*y[iloc];
	    }

	    for(long ic=0; ic<cx->ny; ic++){
		xm[iloc]+=P(cx,0,ic)*xp[(int)P(cx,1,ic)]*yp[(int)P(cx,2,ic)];
	    }

	    for(long ic=0; ic<cy->ny; ic++){
		ym[iloc]+=P(cy,0,ic)*xp[(int)P(cy,1,ic)]*yp[(int)P(cy,2,ic)];
	    }
	}
    }
    return locm;
}
loc_t *loctransform2(const loc_t *loc, const dmat *coeff){
    dmat *cx=dnew(3, coeff->nx);
    dmat *cy=dnew(3, coeff->nx);
    for(int ic=0; ic<cx->ny; ic++){
	P(cx, 0, ic)=P(coeff, ic, 2);//coeff for x
	P(cy, 0, ic)=P(coeff, ic, 3);//coeff for y
	P(cx, 1, ic)=P(cy, 1, ic)=P(coeff, ic, 0);//power for x
	P(cx, 2, ic)=P(cy, 2, ic)=P(coeff, ic, 1);//power for y
    }
    loc_t* locm=loctransform_do(loc, cx, cy);
    dfree(cx);
    dfree(cy);
    return locm;
}
/**
   Transform coordinate. loc (input) contains x,y; locm (return) contains xm, ym.
   polyn contains the formula
   xm{ip}=\{sum}_{ic}(coeff[0](0,ic)*pow(x,coeff[0](1,ic))*pow(y,coeff[0](2,ic)))
   ym{ip}=\{sum}_{ic}(coeff[1](0,ic)*pow(x,coeff[1](1,ic))*pow(y,coeff[1](2,ic)))
   was using string input. New scheme uses bin files to preserve precision.
*/
loc_t *loctransform(const loc_t *loc, const char *polycoeff){
    /*Test whether the transform is pure shift. */
    //Parse from string to 3xn array
    int input_type;
    if(check_suffix(polycoeff, ".bin")){
	input_type=1;
    }else if(*((const uint32_t*)polycoeff)==M_REAL){
	input_type=2;
    }else{
	input_type=0;
    }
    loc_t *locm=0;
    if(input_type>0){
	dmat *coeff=0;
	if(input_type==1){
	    coeff=dread("%s", polycoeff);
	}else if(input_type==2){
	    coeff=(dmat*)polycoeff;
	}
	if(coeff->ny==4){
	    locm=loctransform2(loc, coeff);

	}else{
	    error("coeff is in wrong format\n");
	}
	if(input_type==1){
	    dfree(coeff);
	}
    }else{
	char *polyn=strdup(polycoeff);
	char *px=polyn;
	char *py=strchr(polyn, ';');
	if(py==polyn || !py) error("Wrong format '%s'\n", polyn);
	py[0]='\0';
	py++;
	//info("polyn=(%s)\npx=(%s)\npy=(%s)\n", _polyn, px, py);
	//Now parse the strings.
	dmat *cx=parse_poly(px);
	dmat *cy=parse_poly(py);
	free(polyn); polyn=0;
    
	locm=loctransform_do(loc, cx, cy);
	dfree(cx);
	dfree(cy);
    }
    return locm;
}

/**
   Shift a loc coordinate
*/
loc_t *locshift(const loc_t *loc, real sx, real sy){
    loc_t *loc2=locnew(loc->nloc, loc->dx, loc->dy);
    for(long iloc=0; iloc<loc->nloc; iloc++){
	loc2->locx[iloc]=loc->locx[iloc]+sx;
	loc2->locy[iloc]=loc->locy[iloc]+sy;
    }
    return loc2;
}
/**
   Compute the size of a map that is used to build loc or can fully contain loc
   (embed)
*/
void loc_nxny(long *nx, long *ny, const loc_t *loc){
    real xmax, xmin, ymax, ymin;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    *nx=(long)round((xmax-xmin)/loc->dx)+1;
    *ny=(long)round((ymax-ymin)/loc->dy)+1;
}
/**
   Resize a loc_t by shrinking.
*/
void locresize(loc_t *loc, long nloc){
    if(!loc) return;
    loc_free_map(loc);
    loc_free_stat(loc);
    loc->locx=myrealloc(loc->locx,nloc,real);
    loc->locy=myrealloc(loc->locy,nloc,real);
    loc->nloc=nloc;
}

/**
   Embeding an OPD defined on loc to another array. 
   Do the embeding using locstat to have best speed.
   reverse = 0 : from oin to out: out=out*alpha+in*beta
   reverse = 1 : from out to oin: in=in*beta+out*alpha
*/
#define LOC_EMBED_DEF(X, T, R, OPX)					\
    void X(embed_locstat)(X(mat) **restrict out, real alpha,		\
			  loc_t *restrict loc,				\
			  R *restrict oin, real beta, int reverse){	\
	locstat_t *restrict locstat=loc->stat;				\
	if(!*out){							\
	    if(reverse == 0){						\
		*out=X(new)(locstat->nx, locstat->ny);			\
	    }else{							\
		error("For reverse embedding the array needs to be non-empty\n"); \
	    }								\
	}else{								\
	    if((*out)->nx < locstat->nx || (*out)->ny < locstat->ny){	\
		error("Preallocated array %ldx%ld is too small, we need %ldx%ld\n", \
		      (*out)->nx, (*out)->ny, locstat->nx, locstat->ny); \
	    }								\
	}								\
	X(mat*) p=*out;							\
	real dx1=1./locstat->dx;					\
	long xoff0=((*out)->nx - locstat->nx +1)/2;			\
	long yoff0=((*out)->ny - locstat->ny +1)/2;			\
									\
	for(long icol=0; icol<locstat->ncol; icol++){			\
	    long xoff=(long)round((locstat->cols[icol].xstart-locstat->xmin)*dx1); \
	    long yoff=(long)round((locstat->cols[icol].ystart-locstat->ymin)*dx1); \
	    long pos1=locstat->cols[icol].pos;				\
	    long pos2=locstat->cols[icol+1].pos;			\
	    T *restrict dest=PP(p,xoff+xoff0,yoff+yoff0);		\
	    if(!reverse){						\
		if(oin){						\
		    const R *restrict oin2=oin+pos1;			\
		    if(fabs(alpha)>EPS){				\
			for(long ix=0; ix<pos2-pos1; ix++){		\
			    dest[ix]=dest[ix]*alpha+oin2[ix]*beta;	\
			}						\
		    }else{						\
			for(long ix=0; ix<pos2-pos1; ix++){		\
			    dest[ix]=oin2[ix]*beta;			\
			}						\
		    }							\
		}else{							\
		    if(fabs(alpha)>EPS){				\
			for(long ix=0; ix<pos2-pos1; ix++){		\
			    dest[ix]=dest[ix]*alpha+beta;		\
			}						\
		    }else{						\
			for(long ix=0; ix<pos2-pos1; ix++){		\
			    dest[ix]=beta;				\
			}						\
		    }							\
		}							\
	    }else{							\
		R *restrict oin2=oin+pos1;				\
		if(fabs(beta)>EPS){					\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
			oin2[ix]=oin2[ix]*beta+alpha*OPX(dest[ix]);	\
		    }							\
		}else{							\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
		    oin2[ix]=alpha*OPX(dest[ix]);			\
		}							\
	    }								\
	}								\
    }									\
}

#define OPX(A) (A)
LOC_EMBED_DEF(AOS_DMAT, real, real, OPX)
#undef OPX
#define OPX(A) creal(A)
LOC_EMBED_DEF(AOS_CMAT, comp, real, OPX)
#undef OPX

/**
   Determine dx and dy from data.
*/
void loc_dxdy(loc_t *out){
    const real tol=1e-7;

    real dxd=INFINITY, dyd=INFINITY;
    for(long i=0; i<out->nloc-1; i++){
	real dxi=fabs(out->locx[i+1]-out->locx[i]);
	real dyi=fabs(out->locy[i+1]-out->locy[i]);
	if(dxi>tol && dxi+tol<dxd){
	    dxd=dxi;
	}
	if(dyi>tol && dyi+tol<dyd){
	    dyd=dyi;
	}
    }
    if(out->dx==0 || isnan(out->dx)){
	out->dx=dxd;//use value derived from data
    }else if(fabs(out->dx-dxd)>tol && isfinite(dxd)){
	warning("Specified dx=%.15g doesn't agree with data: %.15g, replace.\n", out->dx, dxd);
	out->dx=dxd;
    }

    if(out->dy==0 || isnan(out->dy)){
	out->dy=dyd;
    }else if(fabs(out->dy-dyd)>tol && isfinite(dyd)){
	warning("Specified dy=%.15g doesn't agree with data: %.15g, replace.\n", out->dy, dyd);
	out->dy=dyd;
    }
}
/**
   Verify the magic, dimension and read in the loc_t by calling locreaddata2().
 */
loc_t *locreaddata(file_t *fp, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    if((header->magic&M_REAL)!=M_REAL){
	error("magic=%u. Expect %x\n", header->magic, M_REAL);
    }
    real dx=fabs(search_header_num(header->str,"dx"));
    real dy=fabs(search_header_num(header->str,"dy"));
    
    free(header->str);header->str=0;
    long nx=header->nx;
    long ny=header->ny;
    loc_t *out;
    if(nx==0 || ny==0){
	out=NULL;
    }else{
	out=locnew(nx, dx, dy);
	readvec(out->locx, M_REAL, header->magic, sizeof(real), nx, fp);
	readvec(out->locy, M_REAL, header->magic, sizeof(real), nx, fp);
	
	loc_dxdy(out);
    }
    return out;
}

void locwritedata(file_t *fp, const loc_t *loc){
    char str[120];
    snprintf(str,120,"dx=%.15g;\ndy=%.15g;iac=%.15g\n",loc->dx,loc->dy,loc->iac);
    header_t header={M_LOC, 0, 0, str};
    if(loc){
	header.nx=loc->nloc;
	header.ny=2;
    }
    write_header(&header,fp);
    if(loc){
	zfwrite(loc->locx, sizeof(real),loc->nloc,fp);
	zfwrite(loc->locy, sizeof(real),loc->nloc,fp);
    }
}
