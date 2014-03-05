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

#include <sys/mman.h>
#include <sys/file.h>
#include <unistd.h>
#include "../sys/sys.h"
#include "../math/mathdef.h"
#include "loc.h"
#include "draw.h"

/**
   Free pts_t data
*/
void ptsfree_do(pts_t *pts){
    if(pts){
	if(pts->origx){
	    free(pts->origx);
	    free(pts->origy);
	}
	if(pts->map){
	    free(pts->map->p);
	    free(pts->map);
	}
	free(pts);
    }
}

/**
   Free the MAP in loc_t
*/
void loc_free_map(loc_t *loc){
    if(loc->map){
	free(loc->map->p);
	free(loc->map);
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
    if(!loc) return;
    loc_free_stat(loc);
    loc_free_map(loc);
    if(loc->locx){
	free(loc->locx);
	free(loc->locy);
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
   Free map_t data
*/
void mapfree_do(map_t *map){
    dfree_do((dmat*)map, 0);
}

/**
   Free map_t array
 */
void maparrfree_do(map_t **map, int nmap){
    if(!map) return;
    for(int imap=0; imap<nmap; imap++){
	mapfree(map[imap]);
    }
    free(map);
}

/**
   Free rectmap_t data
*/
void rectmapfree_do(rectmap_t *map){
    dfree_do((dmat*)map, 0);
}
/**
   Create a loc with nloc elements.
*/
loc_t *locnew(long nloc,double dx, double dy){
    loc_t *loc=calloc(1, sizeof(loc_t));
    loc->id=M_LOC64;
    loc->locx=calloc(nloc, sizeof(double));
    loc->locy=calloc(nloc, sizeof(double));
    loc->nloc=nloc;
    loc->dx=dx;
    loc->dy=dy;
    return loc;
}
/**
   Create a pts with nsa, dsa, nx, dx
*/
pts_t *ptsnew(long nsa, double dsax, double dsay, long nx, double dx, double dy){
    pts_t *pts=calloc(1, sizeof(pts_t));
    pts->id=M_LOC64;
    pts->origx=calloc(nsa, sizeof(double));
    pts->origy=calloc(nsa, sizeof(double));
    pts->dsa=dsax;
    pts->dsay=dsay;
    pts->nx=nx;
    pts->dx=dx;
    pts->dy=dy;
    return pts;
}
/**
   Create an vector to embed OPD into square array for FFT purpose. oversize is
2 for fft.  */
long *loc_create_embed(long *nembed, const loc_t *loc, int oversize){
    double xmin,xmax,ymin,ymax;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    const double dx_in1=1./loc->dx;
    const double dy_in1=1./loc->dy;
    long nx=(long)round((xmax-xmin)*dx_in1)+1;
    long ny=(long)round((ymax-ymin)*dy_in1)+1;
    long nxy=(nx>ny?nx:ny)*oversize;/*minimum size */
    long mapn;
    if(*nembed<=0){
	mapn=nextfftsize(nxy);
	*nembed=mapn;
    }else{
	if(*nembed<(long)(nxy*0.6)){
	    error("Supplied nembed %ld is too small, recommend %ld\n",*nembed, nxy);
	}else if(*nembed<nxy){
	    warning("Supplied nembed %ld maybe too small, recommend %ld\n",*nembed, nxy);
	}
	mapn=*nembed;
    }
    xmin-=(mapn-nx)/2*loc->dx;
    ymin-=(mapn-ny)/2*loc->dy;
    long *embed=calloc(loc->nloc, sizeof(long));
    
    for(int iloc=0; iloc<loc->nloc; iloc++){
	long ix=(long)round((loc->locx[iloc]-xmin)*dx_in1);
	long iy=(long)round((loc->locy[iloc]-ymin)*dy_in1);
	embed[iloc]=ix+iy*mapn;
    }
    return embed;
}
/**
   Create a map for loc with padding of 1. 
*/
void loc_create_map(loc_t *loc){
    loc_create_map_npad(loc,1);
}
PNEW(maplock);
/**
   Create a map for loc so that we can obtain the index in
   loc by x,y coordinate. Useful in ray tracing (accphi.c)
*/
void loc_create_map_npad(loc_t *loc, int npad){
    LOCK(maplock);
    if(loc->map){
	if(loc->map->npad<npad){
	    loc_free_map(loc);
	}else{
	    UNLOCK(maplock);
	    return;
	}
    }
    if(loc->nloc==0){
	UNLOCK(maplock);
	return;
    }
    loc->map = calloc(1,sizeof(locmap_t));
    loc->map->npad = npad;/*just record the information. */
    double xmin,xmax,ymin,ymax;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    int map_nx, map_ny;
    /*padding the map. normally don't need. */
    if(npad>0){
	xmin-=npad*loc->dx;
	ymin-=npad*loc->dy;
	xmax+=npad*loc->dx;
	ymax+=npad*loc->dy;
    }
    const double dx_in1=1./loc->dx;
    const double dy_in1=1./loc->dy;
    map_nx=(int) ceil((xmax-xmin)*dx_in1)+1;
    map_ny=(int) ceil((ymax-ymin)*dy_in1)+1;
    loc->map->p=calloc(map_nx*map_ny,sizeof(long));
    loc->map->nx=map_nx;
    loc->map->ny=map_ny;
    loc->map->ox=xmin;
    loc->map->oy=ymin;
    long (*pmap)[map_nx]=(long(*)[map_nx])loc->map->p;
    const double *locx,*locy;
    locx=loc->locx;
    locy=loc->locy;
    int ix,iy;
    for(long iloc=0; iloc<loc->nloc; iloc++){
	ix=(int)round((locx[iloc]-xmin)*dx_in1);
	iy=(int)round((locy[iloc]-ymin)*dy_in1);
	pmap[iy][ix]=iloc+1;/*start from 1. */
    }
    UNLOCK(maplock);
}

/**
   Create an index array can put embed opd defined in loc from map2loc into
   a square array of size nx*ny
*/
long *map2embed(map_t *amp){
    return d2embed((dmat*)amp);
}
/**
   Convert a map to a loc that collects all positive entries. */
loc_t* map2loc(map_t *map){
    const double dx=map->dx;
    const double dy=map->dy;
    const double ox=map->ox;
    const double oy=map->oy;
    const long nx=map->nx;
    const long ny=map->ny;
    PDMAT((dmat*)map,map0);
    long ix,iy;
    loc_t *loc=locnew(nx*ny, dx, dy);
    long count=0;
    for(iy=0; iy<ny; iy++){
	for(ix=0; ix<nx; ix++){
	    if(map0[iy][ix]>0){
		loc->locx[count]=ix*dx+ox;
		loc->locy[count]=iy*dy+oy;
		count++;
	    }
	}
    }
    loc->locx=realloc(loc->locx, sizeof(double)*count);
    loc->locy=realloc(loc->locy, sizeof(double)*count);
    loc->nloc=count;
    return loc;
}
/**
   convert loc to map.
*/
map_t *loc2map(loc_t *loc){
    double xmax, xmin, ymax, ymin;
    dmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
    dmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
    int nx=ceil((xmax-xmin)/loc->dx)+1;
    int ny=ceil((ymax-ymin)/loc->dy)+1;
    map_t *map=mapnew(nx, ny, loc->dx, loc->dy, NULL);
    map->ox=xmin;
    map->oy=ymin;
    double dx_in1=1./loc->dx;
    double dy_in1=1./loc->dy;
    for(int iloc=0; iloc<loc->nloc; iloc++){
	int ix=(int)round((loc->locx[iloc]-xmin)*dx_in1);
	int iy=(int)round((loc->locy[iloc]-ymin)*dy_in1);
	map->p[ix+iy*nx]=1.;
    }	
    return map;
}
/**
   map is a square array of values of 0, or nonzero. 0 means that point is
   missing. the returned vector of long integers can be used to embed an
   irregular OPD into the square array of map->nx*map->ny. The OPDs only goes to
   locations where the value is nonzero.  
*/
long *d2embed(dmat *map){
    long *embed=malloc(sizeof(long)*map->nx*map->ny);
    long count=0;
    for(long i=0; i<map->nx*map->ny; i++){
	if(map->p[i]>0){
	    embed[count]=i;
	    count++;
	}
    }
    embed=realloc(embed, sizeof(long)*count);
    return embed;
}

/**
   Create 1 dimensional loc with given vector.
*/
loc_t *mk1dloc_vec(double *x, long nx){
    double dx=(x[nx-1]-x[0])/(nx-1);
    loc_t *loc=locnew(nx, dx, dx);
    memcpy(loc->locx, x, sizeof(double)*nx);
    return loc;
}
/**
   Create 1 dimensional loc with origin at x0, sampling of dx, and nx numbers.
*/
loc_t *mk1dloc(double x0, double dx, long nx){
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
loc_t *mksqloc_map(map_t*map){
    return mksqloc(map->nx, map->ny, map->dx, map->dy, map->ox, map->oy);
}
/**
   Create a loc array of size nx*ny at sampling dx with origin at (nx/2, ny/2)
 */
loc_t *mksqloc_auto(long nx, long ny, double dx, double dy){
    return mksqloc(nx,ny,dx,dy,-nx/2*dx,-ny/2*dy);
}
/**
   Create a loc array contains coorindates in a square map of size nx*ny, with
   sampling dx, and at origin ox,oy */
loc_t *mksqloc(long nx, long ny, double dx, double dy, double ox, double oy){
    
    loc_t *loc=locnew(nx*ny, dx, dy);
    long ix,iy;
    double (*locx)[nx]=(double(*)[nx])loc->locx;
    double (*locy)[nx]=(double(*)[nx])loc->locy;
    double y;
    for(iy=0; iy<ny; iy++){
	y=iy*dy+oy;
	for(ix=0; ix<nx; ix++){
	    locx[iy][ix]=ix*dx+ox;
	    locy[iy][ix]=y;
	}
    }
    return loc;
}
/**
   like mksqloc, but roated theta CCW. Useful in creating the
   pixel coordinates in a polar coordinate CCD.
*/
loc_t *mksqlocrot(long nx, long ny, double dx, double dy, double ox, double oy, double theta){
    loc_t *loc=locnew(nx*ny, dx, dy);
    long ix,iy;
    double (*locx)[nx]=(double(*)[nx])loc->locx;
    double (*locy)[nx]=(double(*)[nx])loc->locy;
    double y,x;
    double ct=cos(theta);
    double st=sin(theta);
    for(iy=0; iy<ny; iy++){
	y=iy*dy+oy;
	for(ix=0; ix<nx; ix++){
	    x=ix*dx+ox;
	    locx[iy][ix]=ct*x-st*y;
	    locy[iy][ix]=st*x+ct*y;
	}
    }
    return loc;
}
/**
   Find the point that closes to origin (0,0) in loc. Useful in single point
   piston constraint in creating reconstructor.
*/
int loccenter(loc_t *loc){
    int ipix,jpix=0;
    double r2,r2min=INFINITY;
    
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
dmat *loc_mcc_ptt(const loc_t *loc, const double *amp){
    const int nmod=3;
    double *mod[nmod];
    dmat *mcc=dnew(nmod,nmod);
    mod[0]=NULL;
    mod[1]=loc->locx;
    mod[2]=loc->locy;
    double (*ATA)[nmod]=(double(*)[nmod])mcc->p;
    for(int jmod=0; jmod<nmod; jmod++){
	for(int imod=jmod; imod<nmod; imod++){
	    double tmp=dotdbl(mod[imod],mod[jmod],amp,loc->nloc);
	    ATA[jmod][imod]=ATA[imod][jmod]=tmp;
	}
    }
    return mcc;
}
/**
   same as loc_mcc_ptt, except using pts instead of loc.
   this is used to build imcc for TT/F powfs ztilt imcc
*/
dcell *pts_mcc_ptt(const pts_t *pts, const double *amp){
    const int nmod=3;
    const int nsa=pts->nsa;
    dcell *mcc=dcellnew(nsa,1);
    for(int isa=0; isa<nsa; isa++){
	const double origy=pts->origy[isa];
	const double origx=pts->origx[isa];
	const double dx=pts->dx;
	const double dy=pts->dy;
	const double *ampi=amp+pts->nx*pts->nx*isa;
	mcc->p[isa]=dnew(nmod,nmod);
	double (*ATA)[nmod]=(double(*)[nmod])(mcc->p[isa]->p);
	double a00=0,a01=0,a02=0,a11=0,a12=0,a22=0;
	for(int iy=0; iy<pts->nx; iy++){
	    double y=iy*dy+origy;
	    const double *ampx=ampi+iy*pts->nx;
	    for(int ix=0; ix<pts->nx; ix++){
		double x=ix*dx+origx;
		double a=ampx[ix];
		a00+=a;
		a01+=a*x;
		a02+=a*y;
		a11+=a*x*x;
		a12+=a*x*y;
		a22+=a*y*y;
	    }
	}
	ATA[0][0]=a00;
	ATA[1][1]=a11;
	ATA[2][2]=a22;
	ATA[1][0]=ATA[0][1]=a01;
	ATA[2][0]=ATA[0][2]=a02;
	ATA[2][1]=ATA[1][2]=a12;
    }
    return mcc;
}

/**
   evaluate piston/tip-tilt/ removed wavefront error.
   output coeffout in unit of radian like units.
*/
void loc_calc_ptt(double *rmsout, double *coeffout,
		  const loc_t *loc, const double ipcc, 
		  const dmat *imcc, const double *amp, const double *opd){
    assert(imcc->nx==imcc->ny && imcc->nx==3);
    const long nloc=loc->nloc;
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;

    double tot=0;
    double coeff[3]={0,0,0};
    if(amp){
	for(long iloc=0; iloc<nloc; iloc++){
	    const double junk=opd[iloc]*amp[iloc];
	    coeff[0]+=junk;
	    coeff[1]+=junk*locx[iloc];
	    coeff[2]+=junk*locy[iloc];
	    tot+=junk*opd[iloc];
	}
    }else{
	for(long iloc=0; iloc<nloc; iloc++){
	    const double junk=opd[iloc];
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
	double pis=ipcc*coeff[0]*coeff[0];/*piston mode variance */
	double ptt=dwdot3(coeff, imcc, coeff);/*p/t/t mode variance. */
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
void loc_calc_mod(double *rmsout, double *coeffout,const dmat *mod,
		 const double *amp, double *opd){
  
    const int nmod=mod->ny;
    const int nloc=mod->nx;
    PDMAT(mod,pmod);
 
    double tot=0;
    double val[nmod];
    memset(val, 0, sizeof(double)*nmod);
    for(long iloc=0; iloc<nloc; iloc++){
	double junk=opd[iloc]*amp[iloc];
	tot+=opd[iloc]*junk;
	for(long imod=0; imod<nmod; imod++){
	    val[imod]+=pmod[imod][iloc]*junk;
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
void loc_remove_ptt(double *opd, const double *ptt, const loc_t *loc){
    double ptt1[3];
    ptt1[0]=-ptt[0];
    ptt1[1]=-ptt[1];
    ptt1[2]=-ptt[2];
    loc_add_ptt(opd, ptt1, loc);
}

/**
   Add Piston/Tip/Tilt from OPD
*/
void loc_add_ptt(double *opd, const double *ptt, const loc_t *loc){
    const long nloc=loc->nloc;
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    for(long iloc=0; iloc<nloc; iloc++){
	opd[iloc]+=ptt[0]+ptt[1]*locx[iloc]+ptt[2]*locy[iloc];
    }
}
/**
   Compute zernike best fit for all subapertures. add result to out.  returns
   radians not zernike modes of tip/tilt. Used in wfsgrad */
void pts_ztilt(dmat **out, const pts_t *pts, const dcell *imcc,
	       const double *amp, const double *opd){
    const int nsa=pts->nsa;
    assert(imcc->nx==nsa && imcc->ny==1);
    if(!*out) *out=dnew(nsa*2,1);
    double *res=(*out)->p;
    for(int isa=0; isa<nsa; isa++){
	const double origy=pts->origy[isa];
	const double origx=pts->origx[isa];
	const double dx=pts->dx;
	const double dy=pts->dy;
	const double *ampi=amp+pts->nx*pts->nx*isa;
	const double *opdi=opd+pts->nx*pts->nx*isa;
	assert(imcc->p[isa]->nx==3 && imcc->p[isa]->ny==3);
        double coeff[3]={0,0,0};
	double a0=0,a1=0,a2=0;
	for(int iy=0; iy<pts->nx; iy++){
	    double y=iy*dy+origy;
	    const double *ampx=ampi+iy*pts->nx;
	    const double *opdx=opdi+iy*pts->nx;
	    for(int ix=0; ix<pts->nx; ix++){
		double x=ix*dx+origx;
		double tmp=ampx[ix]*opdx[ix];
		a0+=tmp;
		a1+=tmp*x;
		a2+=tmp*y;
	    }
	}
	coeff[0]=a0;
	coeff[1]=a1;
	coeff[2]=a2;
	double outp[3];
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
    locstat_t *locstat=calloc(1, sizeof(locstat_t));
    loc->stat=locstat;
    const double *locx=loc->locx;
    const double *locy=loc->locy;
    int nloc=loc->nloc;
    double dx=locstat->dx=loc->dx;
    double dy=locstat->dy=loc->dy;
    int ncolmax=(int)round((locy[nloc-1]-locy[0])/dy)+2;
    locstat->cols=malloc(ncolmax*sizeof(locstatcol_t));
    int colcount=0;
    /*do first column separately. */
    int iloc=0;
    locstat->cols[colcount].pos=iloc;
    locstat->cols[colcount].xstart=locx[iloc];
    locstat->cols[colcount].ystart=locy[iloc];
    locstat->ymin=locstat->cols[colcount].ystart;
    locstat->xmin=locstat->cols[colcount].xstart;
    double xmax=locstat->cols[colcount].xstart;

    colcount++;
    for(iloc=1; iloc<loc->nloc; iloc++){
	if(fabs(locy[iloc]-locy[iloc-1])>1.e-12 /*a new column starts */
	   || fabs(locx[iloc]-locx[iloc-1]-dx)>1.e-12){
	    if(colcount>=ncolmax){/*expand the memory. */
		ncolmax*=2;
		locstat->cols=realloc(locstat->cols, ncolmax*sizeof(locstatcol_t));
	    }
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
	}
    }
    locstat->ncol=colcount;
    locstat->cols[colcount].pos=loc->nloc;
    locstat->cols[colcount].xstart=locx[loc->nloc-1];
    locstat->cols[colcount].ystart=locy[loc->nloc-1];
    if(xmax < locx[loc->nloc-1]){
	xmax = locx[loc->nloc-1];
    }
    locstat->nrow=(long)round((xmax-locstat->xmin)/dx)+1;
    locstat->cols=realloc(locstat->cols,(locstat->ncol+1)*sizeof(locstatcol_t));
}
/**
   Create a gray pixel circular map in phi using coordinates defined in loc, center
defined using cx, cy, radius of r, and value of val */
void loccircle(double *phi,loc_t *loc,double cx,double cy,double r,double val){
    /*cx,cy,r are in unit of true unit, as in loc */
    double dx=loc->dx;
    double dy=loc->dy;
    int nres=10;
    double cres=(nres-1.)/2.;
    double res=1./nres;
    double dxres=res*dx;
    double dyres=res*dy;
    double res2=res*res;
    double r2=r*r;
    double r2l=(r-dx*1.5)*(r-dy*1.5);
    double r2u=(r+dx*1.5)*(r+dy*1.5);
    double *locx=loc->locx;
    double *locy=loc->locy;
    double iix,iiy;
    for(int iloc=0; iloc<loc->nloc; iloc++){
	double x=locx[iloc];
	double y=locy[iloc];
	double r2r=pow(x-cx,2)+pow(y-cy,2);
	if(r2r<r2l) 
	    phi[iloc]+=val;
	else if(r2r<r2u){
	    long tot=0;
	    for(int jres=0; jres<nres; jres++){
		iiy=y+(jres-cres)*dyres;
		double rr2y=(iiy-cy)*(iiy-cy);
		for(int ires=0; ires<nres; ires++){
		    iix=x+(ires-cres)*dxres;
		    double rr2r=(iix-cx)*(iix-cx)+rr2y;
		    if(rr2r<=r2)
			tot++;
		}
	    }
	    phi[iloc]+=(double)tot*res2*val;
	}
    }
}
/**
   Create a gray pixel annular map in phi using coordinates defined in loc,
center defined using cx, cy, radius of r, and value of val
 */
void locannular(double *phi,loc_t *loc,double cx,double cy,double r,double rin,double val){
    loccircle(phi,loc,cx,cy,r,val);
    if(rin>EPS){
	loccircle(phi,loc,cx,cy,rin,-val);
    }
}
/**
   Create a hard annular mask in phi.
*/
void locannularmask(double *phi,loc_t *loc,double cx,double cy,double r,double rin){
    /*apply the hard pupil mask of aper.d, using loc. not locm */
    /* 2011-07-13: changed from r^2 to (r+0.5*dx)^2*/
    double rr2max=pow(r+0.25*(loc->dx+loc->dy),2);
    double rr2min=MIN(rin*rin, pow(rin-0.25*(loc->dx+loc->dy),2));
    for(long iloc=0; iloc<loc->nloc; iloc++){
	double r2=pow(loc->locx[iloc],2)+pow(loc->locy[iloc],2);
	if(r2<rr2min || r2>rr2max){
	    phi[iloc]=0;
	}
    }
}
/**
 Create a gray pixel elliptical map in phi using coordinates defined in loc,
center defined using cx, cy, radii of rx, ry, and value of val */
void locellipse(double *phi,loc_t *loc,double cx,double cy,
		double rx,double ry,double val){
    /*cx,cy,r are in unit of true unit, as in loc */
    double dx=loc->dx;
    double dy=loc->dy;
    int nres=10;
    double cres=(nres-1.)/2.;
    double res=1./nres;
    double dxres=res*dx;
    double dyres=res*dy;
    double res2=res*res;
    double rx1=1./rx;
    double ry1=1./ry;
    double rx12=rx1*rx1;
    double ry12=ry1*ry1;
    double r2=2.;
    double r2l=(rx-dx)*(rx-dx)*rx12+(ry-dy)*(ry-dy)*ry12;
    double r2u=(rx+dx)*(rx+dx)*rx12+(ry+dy)*(ry+dy)*ry12;
    double *locx=loc->locx;
    double *locy=loc->locy;
    double iix,iiy;
    for(int iloc=0; iloc<loc->nloc; iloc++){
	double x=locx[iloc];
	double y=locy[iloc];
	double r2r=pow(x-cx,2)*rx12+pow(y-cy,2)*ry12;
	if(r2r<r2l) 
	    phi[iloc]+=val;
	else if(r2r<r2u){
	    long tot=0;
	    for(int jres=0; jres<nres; jres++){
		iiy=y+(jres-cres)*dyres;
		double rr2y=(iiy-cy)*(iiy-cy)*ry12;
		for(int ires=0; ires<nres; ires++){
		    iix=x+(ires-cres)*dxres;
		    double rr2r=(iix-cx)*(iix-cx)*rx12+rr2y;
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

void loc_reduce(loc_t *loc, dmat *amp, int cont, int **skipout){
    int redo_stat=loc->stat?1:0;
    loc_free_map(loc);/*remove the internal map before touchlong loc. */
    loc_free_stat(loc);
    int nloc=loc->nloc; 
    int *skip=calloc(nloc,sizeof(int));
    if(cont){/*make sure loc is continuous. */
	loc_create_stat(loc);
	locstat_t *locstat=loc->stat;
	int ncol=locstat->ncol;
	for(int icol=0; icol<ncol; icol++){
	    int pstart=locstat->cols[icol].pos;
	    int pend=locstat->cols[icol+1].pos;
	    for(int pos=pstart; pos<pend; pos++){
		if(amp->p[pos]<EPS){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	    for(int pos=pend-1; pos>pstart-1; pos--){
		if(amp->p[pos]<EPS){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	}
	loc_free_stat(loc);
    }else{
	for(int iloc=0; iloc<nloc; iloc++){
	    if(amp->p[iloc]<EPS){
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
    loc->locx=realloc(loc->locx,sizeof(double)*count);
    loc->locy=realloc(loc->locy,sizeof(double)*count);
    loc->nloc=count;
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
   a spcell, to compute the coupling which is modified accordingly.  */
void loc_reduce_spcell(loc_t *loc, spcell *spc, int dim, int cont){
    loc_free_map(loc);
    loc_free_stat(loc);
    int nloc=loc->nloc;
    dmat *sum=NULL;
    for(int isp=0; isp<spc->nx*spc->ny; isp++){
	if(spc->p[isp]){
	    dmat* sum0=spsumabs(spc->p[isp],3-dim);
	    dadd(&sum,1,sum0,1);
	    dfree(sum0);
	}
    }
    if(!sum || sum->nx*sum->ny!=loc->nloc){
	error("Mismatch\n");
    }
    int *skip;
    loc_reduce(loc, sum, cont, &skip);
    dfree(sum);
    int count;
    if(dim==1){/*opd(loc)=sp*opd(locin); */
	count=0;
	int *map=calloc(nloc,sizeof(int));
	for(int iloc=0; iloc<nloc; iloc++){
	    map[iloc]=count;
	    if(!skip[iloc]) count++;
	}
	for(int isp=0; isp<spc->nx*spc->ny;isp++){
	    dsp *sp=spc->p[isp];
	    if(!sp) continue;
	    sp->m=count;
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
	    sp->n=count;
	    sp->p=realloc(sp->p,sizeof(long)*(count+1));
	}
    }
    free(skip);
}

/**
   Remove uncoupled points in loc and modify sp in the same time. debugged on 2009-12-20.  use single sparse
   matrix, which is modified accordingly.  */
void loc_reduce_sp(loc_t *loc, dsp *sp, int dim, int cont){
    loc_free_map(loc);/*remove the internal map before touchlong loc. */
    loc_free_stat(loc);
    int nloc=loc->nloc;
    if((dim==1 && nloc!=sp->m) || (dim==2 && nloc!=sp->n) || dim<0 || dim>2)
	error("Mismatch dimension\n");
    dmat* sum=spsumabs(sp,3-dim);
    int *skip;
    loc_reduce(loc, sum, cont, &skip);
    dfree(sum);
    int count=0;
    if(dim==1){/*opd(loc)=sp*opd(locin); */
	int *map=calloc(nloc,sizeof(int));
	for(int iloc=0; iloc<nloc; iloc++){
	    map[iloc]=count;
	    if(!skip[iloc]) count++;
	}
	sp->m=count;
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
	sp->n=count;
	sp->p=realloc(sp->p,sizeof(long)*(count+1));
    }
    free(skip);
}

/**
   Generating Zernike Rnm for radial order ir, azimuthal or im.  used by loc_zernike.
   */
static dmat *genRnm(dmat *locr, int ir, int im){
    if(ir<0 || im < 0|| im>ir || (ir-im)%2!=0)
	error("Invalid ir, im (%d, %d)\n", ir, im);
    
    const long nloc=locr->nx;
    dmat *Rnm=dnew(nloc,1);
    for(int s=0; s<=(ir-im)/2; s++){
	double coeff=pow(-1,s)*factorial(ir-s)/factorial(s)/factorial((ir+im)/2-s)/factorial((ir-im)/2-s);
	int power=ir-2*s;
	if(power==0){
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff;
	    }
	}else if(power==1){
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff*locr->p[iloc];
	    }
	}else{
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff*pow(locr->p[iloc],power);
	    }
	}
    }
    return Rnm;
}

/**
   Create Zernike modes on loc with radius R and radial order upto nr
   nr=0 is piston
   nr=1 is tip/tilt
   nr=2 is quadratic modes
   
*/
dmat* loc_zernike(loc_t *loc, double R, int nr){
    if(nr<0) error("Invalid nr\n");
    int nmod=(nr+1)*(nr+2)/2;
    const long nloc=loc->nloc;
    dmat *restrict MOD=dnew(nloc,nmod);
    
    dmat *restrict locr=dnew(nloc,1);
    dmat *restrict locs=dnew(nloc,1);
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    const double R1=1./R;
    for(long iloc=0; iloc<nloc; iloc++){
	locr->p[iloc]=sqrt(pow(locx[iloc],2)+pow(locy[iloc],2))*R1;
	locs->p[iloc]=atan2(locy[iloc], locx[iloc]);
    }
    int cmod=0;
    for(int ir=0; ir<=nr; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    dmat *Rnm=genRnm(locr, ir, im);
	    /*dwrite(Rnm,"Rnm_%d_%d",ir,im); */
	    if(im==0){
		double coeff=sqrt(ir+1.);
		double *restrict pmod=MOD->p+nloc*cmod;
		for(long iloc=0; iloc<nloc; iloc++){
		    pmod[iloc]=Rnm->p[iloc]*coeff;
		}
		cmod++;
	    }else{
		double coeff=sqrt(2*(ir+1.));
		double *restrict pmodc;
		double *restrict pmods;
		if((cmod+1) % 2 == 1){
		    pmods=MOD->p+nloc*cmod;
		    pmodc=MOD->p+nloc*(cmod+1);
		}else{
		    pmodc=MOD->p+nloc*cmod;
		    pmods=MOD->p+nloc*(cmod+1);
		}
		for(long iloc=0; iloc<nloc; iloc++){
		    pmods[iloc]=Rnm->p[iloc]*coeff*sin(im*locs->p[iloc]);
		    pmodc[iloc]=Rnm->p[iloc]*coeff*cos(im*locs->p[iloc]);
		}
		cmod+=2;
	    }
	    dfree(Rnm);
	}
    }
    dfree(locr);
    dfree(locs);
    return MOD;
}

/**
   Add val amount of focus to opd. The unit is in radian like.
   Piston is not removed.
*/
void loc_add_focus(double *opd, loc_t *loc, double val){
    if(fabs(val)<1.e-15) return;
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    for(long iloc=0; iloc<loc->nloc; iloc++){
	opd[iloc]+=(locx[iloc]*locx[iloc]+locy[iloc]*locy[iloc])*val;
    }
}
/**
   Create a dmat containing locx and locy in two columns
*/
dmat *loc2mat(loc_t *loc,int piston){
    dmat *out=NULL;
    double *ptt=NULL;
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
    memcpy(ptt,loc->locx,sizeof(double)*loc->nloc);
    memcpy(ptt+loc->nloc,loc->locy,sizeof(double)*loc->nloc);
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
    double dx = pts->dx;
    double dy = pts->dy;
    loc_t *loc = locnew(nsa*nxsa, dx, dy);
    for(int isa=0; isa<nsa; isa++){
	const double origx = pts->origx[isa];
	const double origy = pts->origy[isa];
	double *locx = loc->locx+isa*nxsa;
	double *locy = loc->locy+isa*nxsa;
	for(int jj=0; jj<nx; jj++){
	    const double yy=origy+(double)jj*dy;
	    for(int ii=0; ii<nx; ii++){
		locx[jj*nx+ii]=origx+(double)ii*dx;
		locy[jj*nx+ii]=yy;
	    }
	}
    }
    return loc;
}

/**
   Rotate the coordinates by theta CCW.
*/
void locrot(loc_t *loc, const double theta){
    const double ctheta=cos(theta);
    const double stheta=sin(theta);

    double *x=loc->locx;
    double *y=loc->locy;
    for(int i=0; i<loc->nloc; i++){
	double tmp=x[i]*ctheta-y[i]*stheta;
	y[i]=x[i]*stheta+y[i]*ctheta;
	x[i]=tmp;
    }
}

/**
   duplicate a loc_t object; */
loc_t *locdup(loc_t *loc){
    loc_t *res=locnew(loc->nloc, loc->dx, loc->dy);
    memcpy(res->locx,loc->locx,sizeof(double)*loc->nloc);
    memcpy(res->locy,loc->locy,sizeof(double)*loc->nloc);
    return res;
}
/**Parse string representation of polynominal to array representation.*/
static int parse_poly(double (**pcx)[3], const char *_ps){
    char *ps=(char *)_ps;
    int ncx=5;
    *pcx=malloc(3*ncx*sizeof(double));
    double (*cx)[3]=*pcx;
    int icx=0;
    char *endptr;
    while(ps[0] && ps[0]!=';'){
	cx[icx][0]=strtod(ps, &endptr);//coefficient
	if(ps==endptr){
	    if(ps[0]=='-'){
		ps++;
		cx[icx][0]=-1;
	    }else if(ps[0]=='+'){
		ps++;
		cx[icx][0]=1;
	    }
	    if(ps[0]=='x' || ps[0]=='y'){
		if(cx[icx][0]==0) cx[icx][0]=1;
	    }else{
		error("Unable to parse (%s). ps=(%s)\n", _ps, ps);
	    }
	}else{
	    ps=endptr;
	}
	cx[icx][1]=cx[icx][2]=0;
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
		    cx[icx][iy]+=strtol(ps, &endptr, 10);
		    if(ps==endptr){
			error("Unable to parse %s\n", _ps);
		    }
		    ps=endptr;
		}else{
		    cx[icx][iy]++;
		}
	    }else{
		cx[icx][0]*=strtod(ps, &endptr);
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
		cx=*pcx=realloc(*pcx, 3*ncx*sizeof(double));
	    }
	}else{
	    error("Unable to parse (%s)\n", ps);
	}
    }
    ncx=icx;
    cx=*pcx=realloc(*pcx, 3*ncx*sizeof(double));
    return ncx;
    
}
/**
   Transform coordinate. loc contains x,y; locm contains xm, ym.
   polyn contains the formula
   xm{ip}=\{sum}_{ic}(coeff[0](1,ic)*pow(x,coeff[0](2,ic))*pow(y,coeff[0](3,ic)))
   ym{ip}=\{sum}_{ic}(coeff[1](1,ic)*pow(x,coeff[1](2,ic))*pow(y,coeff[1](3,ic)))
 */
loc_t *loctransform(loc_t *loc, const char *_polyn){
    const double *restrict x=loc->locx;
    const double *restrict y=loc->locy;
    /*Test whether the transform is pure shift. */
    loc_t *locm=locnew(loc->nloc, loc->dx, loc->dy);
    double *restrict xm=locm->locx;
    double *restrict ym=locm->locy;
    //Parse from string to 3xn array
    double (*cx)[3]=0, (*cy)[3]=0;
    int ncx=0, ncy=0;
    {
	char *polyn=strdup(_polyn);
	char *px=polyn;
	char *py=strchr(polyn, ';');
	if(py==polyn || !py) error("Wrong format '%s'\n", polyn);
	py[0]='\0';
	py++;
	//info2("polyn=(%s)\npx=(%s)\npy=(%s)\n", _polyn, px, py);
	//Now parse the strings.
	ncx=parse_poly(&cx, px);
	ncy=parse_poly(&cy, py);
	free(polyn); polyn=0;px=0; py=0;
    }
    int np=0;
    int nonint=0;
    for(int ic=0; ic<ncx; ic++){
	int ocx=(int)round(cx[ic][1]);
	int ocy=(int)round(cx[ic][2]);
	if((cx[ic][1]-ocx)>EPS || (cx[ic][2]-ocy)>EPS){//not integer
	    nonint=1;
	}
	if(ocx>np) np=ocx;
	if(ocy>np) np=ocy;
    }
    for(int ic=0; ic<ncy; ic++){
	int ocx=(int)round(cy[ic][1]);
	int ocy=(int)round(cy[ic][2]);
	if((cy[ic][1]-ocx)>EPS || (cy[ic][2]-ocy)>EPS){//not integer
	    nonint=1;
	}
	if(ocx>np) np=ocx;
	if(ocy>np) np=ocy;
    }
    np++; //max order of x or y +1
    for(long iloc=0; iloc<loc->nloc; iloc++){
	if(nonint){
	    for(long ic=0; ic<ncx; ic++){
		xm[iloc]+=cx[ic][0]*pow(x[iloc],cx[ic][1])*pow(y[iloc],cx[ic][2]);
	    }
		
	    for(long ic=0; ic<ncy; ic++){
		ym[iloc]+=cy[ic][0]*pow(x[iloc],cy[ic][1])*pow(y[iloc],cy[ic][2]);
	    }
	}else{/*faster method for integer powers (>10x speed up). */
	    double xp[np], yp[np];
	    xp[0]=1; yp[0]=1;
	    for(long ip=1; ip<np; ip++){
		xp[ip]=xp[ip-1]*x[iloc];
		yp[ip]=yp[ip-1]*y[iloc];
	    }

	    for(long ic=0; ic<ncx; ic++){
		xm[iloc]+=cx[ic][0]*xp[(int)cx[ic][1]]*yp[(int)cx[ic][2]];
	    }

	    for(long ic=0; ic<ncy; ic++){
		ym[iloc]+=cy[ic][0]*xp[(int)cy[ic][1]]*yp[(int)cy[ic][2]];
	    }
	}
    }
    return locm;
}
/**
   Shift a loc coordinate
*/
loc_t *locshift(const loc_t *loc, double sx, double sy){
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
    double xmax, xmin, ymax, ymin;
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
    loc->locx=realloc(loc->locx, sizeof(double)*nloc);
    loc->locy=realloc(loc->locy, sizeof(double)*nloc);
    loc->nloc=nloc;
}
/**
   create a new map_t object.
*/
map_t *mapnew(long nx, long ny, double dx, double dy, double *p){
    map_t *map=realloc(dnew_data(nx, ny, p),sizeof(map_t));
    map->h=0;
    map->dx=dx;
    map->dy=dy;
    map->ox=-map->nx/2*map->dx;
    map->oy=-map->ny/2*map->dy;
    map->vx=0;
    map->vy=0;
    map->cubic=0;
    map->iac=0;
    return map;
}
/**
   ceate a new map_t object from existing one. P is left empty.
 */
map_t *mapnew2(map_t* A){
    map_t *map=realloc(dnew_data(A->nx, A->ny, NULL),sizeof(map_t));
    map->h=A->h;
    map->dx=A->dx;
    map->dy=A->dy;
    map->ox=A->ox;
    map->oy=A->oy;
    map->vx=A->vx;
    map->vy=A->vy;
    map->cubic=A->cubic;
    map->iac=A->iac;
    return map;
}
/**
   Create a circular aperture on map_t.
*/
void mapcircle(map_t *map, double r, double val){
    dcircle((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r, val);
}
/**
   Create a circular aperture on map_t.
*/
void mapcircle_symbolic(map_t *map, double r){
    dcircle_symbolic((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r);
}
/**
   Find the inner and outer diameter of an amplitude map contained in map_t.
*/
void map_d_din(map_t *map, double *d, double *din){
    PDMAT(map, p);
    double r2min=INFINITY, r2max=0;
    for(long iy=0; iy<map->ny; iy++){
	double y=iy*map->dy+map->oy;
	for(long ix=0; ix<map->nx; ix++){
	    if(p[iy][ix]>EPS){
		double x=ix*map->dx+map->ox;
		double r2=x*x+y*y;
		if(r2>r2max) r2max=r2;
		if(r2<r2min) r2min=r2;
	    }
	}
    }
    *d=sqrt(r2max)*2;
    *din=sqrt(r2min)*2;
}

/**
   Embeding an OPD defined on loc to another array. *out is dmat or cmat depends
   on iscomplex. Do the embeding using locstat to have best speed.
   reverse = 0 : from oin to out: out=out*alpha+in*beta
   reverse = 1 : from out to oin: in=in*beta+out*alpha
*/
#define LOC_EMBED_DEF(X, T, R, REAL)					\
void X(embed_locstat)(X(mat) **restrict out, double alpha,		\
		      loc_t *restrict loc,				\
		      R *restrict oin, double beta, int reverse){	\
    locstat_t *restrict locstat=loc->stat;				\
    if(!*out){								\
	if(reverse == 0){						\
	    *out=X(new)(locstat->nrow, locstat->ncol);			\
	}else{								\
	    error("For reverse embedding the array needs to be non-empty\n"); \
	}								\
    }else{								\
	if((*out)->nx < locstat->nrow || (*out)->ny < locstat->ncol){	\
	    error("Preallocated array %ldx%ld is too small, we need %ldx%ld\n",	\
		  (*out)->nx, (*out)->ny, locstat->nrow, locstat->ncol); \
	}								\
    }									\
    T (*restrict p)[(*out)->nx]=(void *)(*out)->p;			\
    double dx1=1./locstat->dx;						\
    long xoff0=((*out)->nx - locstat->nrow +1)/2;			\
    long yoff0=((*out)->ny - locstat->ncol +1)/2;			\
									\
    for(long icol=0; icol<locstat->ncol; icol++){			\
	long xoff=(long)round((locstat->cols[icol].xstart-locstat->xmin)*dx1); \
	long yoff=(long)round((locstat->cols[icol].ystart-locstat->ymin)*dx1); \
	long pos1=locstat->cols[icol].pos;				\
	long pos2=locstat->cols[icol+1].pos;				\
	T *restrict dest=&p[yoff+yoff0][xoff+xoff0];			\
	if(!reverse){							\
	    if(oin){							\
		const R *restrict oin2=oin+pos1;			\
		if(fabs(alpha)>EPS){					\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
			dest[ix]=dest[ix]*alpha+oin2[ix]*beta;		\
		    }							\
		}else{							\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
			dest[ix]=oin2[ix]*beta;				\
		    }							\
		}							\
	    }else{							\
		if(fabs(alpha)>EPS){					\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
			dest[ix]=dest[ix]*alpha+beta;			\
		    }							\
		}else{							\
		    for(long ix=0; ix<pos2-pos1; ix++){			\
			dest[ix]=beta;					\
		    }							\
		}							\
	    }								\
	}else{								\
	    R *restrict oin2=oin+pos1;					\
	    if(fabs(beta)>EPS){						\
		for(long ix=0; ix<pos2-pos1; ix++){			\
		    oin2[ix]=oin2[ix]*beta+alpha*REAL(dest[ix]);	\
		}							\
	    }else{							\
		for(long ix=0; ix<pos2-pos1; ix++){			\
		    oin2[ix]=alpha*REAL(dest[ix]);			\
		}							\
	    }								\
	}								\
    }									\
}

#define REAL(A) (A)
LOC_EMBED_DEF(AOS_DMAT, double, double, REAL)
#undef REAL
#define REAL(A) creal(A)
LOC_EMBED_DEF(AOS_CMAT, dcomplex, double, REAL)
#undef REAL
