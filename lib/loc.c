/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
#include "common.h"
#include "misc.h"
#include "loc.h"
#include "cell.h"
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"
#include "bin.h"
#include "draw.h"
#include "matbin.h"
#include "mathmisc.h"
/**
   \file loc.c
   This file defines functions relates to PTS_T, LOC_T, MAP_T, etc.
 */
/**
   Free PTS_T data
*/
void ptsfree_do(PTS_T *pts){
    if(pts){
	if(pts->origx){
	    free(pts->origx);
	    free(pts->origy);
	}
	if(pts->map && pts->map->loc==(LOC_T*)pts){
	    free(pts->map->p);
	    free(pts->map);
	}
	if(pts->area){
	    free(pts->area);
	}
	free(pts);
    }
}
/**
   Free LOCSTAT_T data
*/
void locstatfree_do(LOCSTAT_T *locstat){
    free(locstat->cols);
    free(locstat);
}
/**
   Free the MAP in LOC_T
*/
void loc_free_map(LOC_T *loc){
    if(loc->map && loc->map->loc==loc){
	free(loc->map->p);
	loc->map->p=NULL;
	loc->map->loc=NULL;
	free(loc->map);
	loc->map=NULL;
    }
}
/**
   Free LOC_T data
 */
void locfree_do(LOC_T *loc){
    if(!loc) return;
    if(loc->locx){
	free(loc->locx);
	free(loc->locy);
	loc_free_map(loc);
    }
    free(loc);
}

/**
   Free LOC_T array
*/
void locarrfree_do(LOC_T **loc, int nloc){
    if(!loc) return;
    for(int iloc=0; iloc<nloc; iloc++){
	locfree(loc[iloc]);
    }
    free(loc);
}
/**
   Free MAP_T data
*/
void sqmapfree_do(MAP_T *map){
    if (!map) return;
    if(map->shm){
#if USE_POSIX_SHM
	shm_unmap(map->p, map->shm);
#else
	error("Fatal Error. Should never happen\n");
#endif
	free(map);
    }else{
	free(map->p);
	map->p=NULL;
	free(map);
    }
}

/**
   Free MAP_T array
 */
void sqmaparrfree_do(MAP_T **map, int nmap){
    if(!map) return;
    for(int imap=0; imap<nmap; imap++){
	sqmapfree(map[imap]);
    }
    free(map);
}

/**
   Free RECTMAP_T data
*/
void rectmapfree_do(RECTMAP_T *map){
    if (!map) return;
    free(map->p);
    free(map);
}

/**
   Create an vector to embed OPD into square array for FFT purpose.
*/
int *loc_create_embed(int *nembed, const LOC_T *loc){
    double xmin,xmax,ymin,ymax,tmp;
    maxmindbl(loc->locx, loc->nloc, &xmax, &xmin,&tmp);
    maxmindbl(loc->locy, loc->nloc, &ymax, &ymin,&tmp);
    const double dx_in1=1./loc->dx;
    int nx=(int)round((xmax-xmin)*dx_in1)+1;
    int ny=(int)round((ymax-ymin)*dx_in1)+1;
    int nxy=(nx>ny?nx:ny)*2;//minimum size
    int mapn;
    if(*nembed<=0){
	mapn=(int)exp2(ceil(log2(nxy)));
	*nembed=mapn;
    }else{
	if(*nembed<nxy){
	    error("Supplied nembed %d is not valid, need at least %d\n",*nembed,nxy);
	}
	mapn=*nembed;
    }
    info2("Science PSF is using grid size of %d\n",mapn);
    xmin-=(mapn-nx)/2*loc->dx;
    ymin-=(mapn-ny)/2*loc->dx;
    int *embed=calloc(loc->nloc, sizeof(int));
    
    for(int iloc=0; iloc<loc->nloc; iloc++){
	int ix=(int)round((loc->locx[iloc]-xmin)*dx_in1);
	int iy=(int)round((loc->locy[iloc]-ymin)*dx_in1);
	embed[iloc]=ix+iy*mapn;
    }
    return embed;
}
/**
   Create a map for loc with padding of 1. 
*/
void loc_create_map(LOC_T *loc){
    loc_create_map_npad(loc,1);
}
PNEW(maplock);
/**
   Create a map for loc so that we can obtain the index in
   loc by x,y coordinate. Useful in ray tracing (accphi.c)
*/
void loc_create_map_npad(LOC_T *loc, int npad){
    LOCK(maplock);
    if(loc->map){
	if(loc->map->loc==loc)//already have a map that is valid.
	    if(loc->map->nloc!=loc->nloc ||loc->map->npad<npad){
		loc_free_map(loc);
	    }else{
		UNLOCK(maplock);
		return;
	    }
	else
	    loc->map=NULL;
    }
    loc->map = calloc(1,sizeof(LOCMAP_T));
    loc->map->loc = loc;
    loc->map->nloc= loc->nloc;
    loc->map->npad = npad;//just record the information.
    double xmin,xmax,ymin,ymax,tmp;
    maxmindbl(loc->locx, loc->nloc, &xmax, &xmin,&tmp);
    maxmindbl(loc->locy, loc->nloc, &ymax, &ymin,&tmp);
    int map_nx, map_ny;
    //padding the map. normally don't need.
    if(npad>0){
	xmin-=npad*loc->dx;
	ymin-=npad*loc->dx;
	xmax+=npad*loc->dx;
	ymax+=npad*loc->dx;
    }
    const double dx_in1=1./loc->dx;
    map_nx=(int) ceil((xmax-xmin)*dx_in1)+1;
    map_ny=(int) ceil((ymax-ymin)*dx_in1)+1;
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
	iy=(int)round((locy[iloc]-ymin)*dx_in1);
	pmap[iy][ix]=iloc+1;//start from 1.
    }
    UNLOCK(maplock);
}

/**
   convert a sqmap to a loc object. only positive numbers are
   converted.
*/
LOC_T* sqmap2loc(MAP_T *amp){
    return map2loc(amp->dx,amp->nx,amp->ny,amp->ox,amp->oy,amp->p);
}
/**
   Create an index array can put embed opd defined in loc from map2loc into
   a square array of size nx*ny
*/
long *sqmap2embed(MAP_T *amp){
    return map2embed(amp->nx,amp->ny,amp->p);
}
/**
   Convert a map to a loc that collects all positive entries. */
LOC_T* map2loc(double dx, long nx, long ny, 
	       double ox, double oy, double *map){
    double (*map0)[nx]=(double(*)[nx]) map;
    long ix,iy;
    LOC_T *loc=calloc(1,sizeof(LOC_T));;
    loc->locx=malloc(sizeof(double)*nx*ny);
    loc->locy=malloc(sizeof(double)*nx*ny);

    long count=0;
    for(iy=0; iy<ny; iy++){
	for(ix=0; ix<nx; ix++){
	    if(map0[iy][ix]>0){
		loc->locx[count]=ix*dx+ox;
		loc->locy[count]=iy*dx+oy;
		count++;
	    }
	}
    }
    loc->locx=realloc(loc->locx, sizeof(double)*count);
    loc->locy=realloc(loc->locy, sizeof(double)*count);
    loc->nloc=count;
    loc->dx=dx;
    loc->map=NULL;
    return loc;
}
/**
   Create an index array that can be used to embed opd defined in loc from
   map2loc into a square array of size nx*ny */
long *map2embed(long nx, long ny, double *map){
    long *embed=malloc(sizeof(long)*nx*ny);
    long count=0;
    for(long i=0; i<nx*ny; i++){
	if(map[i]>0){
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
LOC_T *mk1dloc_vec(double *x, long nx){
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    loc->nloc=nx;
    loc->dx=(x[nx-1]-x[0])/(nx-1);
    loc->locx=malloc(sizeof(double)*nx);
    loc->locy=calloc(sizeof(double),nx);
    memcpy(loc->locx, x, sizeof(double)*nx);
    return loc;
}
/**
   Create 1 dimensional loc with origin at x0, sampling of dx, and nx numbers.
*/
LOC_T *mk1dloc(double x0, double dx, long nx){
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    loc->nloc=nx;
    loc->dx=dx;
    loc->locx=malloc(sizeof(double)*nx);
    loc->locy=calloc(sizeof(double),nx);//initialize to zero.
    for(long ix=0; ix<nx; ix++){
	loc->locx[ix]=x0+dx*ix;
    }
    return loc;
}
/**
   Create a loc array that covers the MAP_T. Notice that it is different from
   sqmap2loc which only covers valid (value>0) regions.
 */
LOC_T *mksqloc_map(MAP_T*map){
    return mksqloc(map->nx, map->ny, map->dx, map->ox, map->oy);
}
/**
   Create a loc array of size nx*ny at sampling dx with origin at (nx/2, ny/2)
 */
LOC_T *mksqloc_auto(long nx, long ny, double dx){
    return mksqloc(nx,ny,dx,-nx/2*dx,-ny/2*dx);
}
/**
   Create a loc array contains coorindates in a square map of size nx*ny, with
   sampling dx, and at origin ox,oy */
LOC_T *mksqloc(long nx, long ny, double dx, double ox, double oy){
    
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    loc->nloc=nx*ny;
    loc->dx=dx;
    loc->locx=malloc(sizeof(double)*loc->nloc);
    loc->locy=malloc(sizeof(double)*loc->nloc);
    long ix,iy;
    double (*locx)[nx]=(double(*)[nx])loc->locx;
    double (*locy)[nx]=(double(*)[nx])loc->locy;
    double y;
    for(iy=0; iy<ny; iy++){
	y=iy*dx+oy;
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
LOC_T *mksqlocrot(long nx, long ny, double dx, double ox, double oy, double theta){
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    loc->nloc=nx*ny;
    loc->dx=dx;
    loc->locx=malloc(sizeof(double)*loc->nloc);
    loc->locy=malloc(sizeof(double)*loc->nloc);
    long ix,iy;
    double (*locx)[nx]=(double(*)[nx])loc->locx;
    double (*locy)[nx]=(double(*)[nx])loc->locy;
    double y,x;
    double ct=cos(theta);
    double st=sin(theta);
    for(iy=0; iy<ny; iy++){
	y=iy*dx+oy;
	for(ix=0; ix<nx; ix++){
	    x=ix*dx+ox;
	    locx[iy][ix]=ct*x-st*y;
	    locy[iy][ix]=st*x+ct*y;
	}
    }
    return loc;
}
/**
   Create rough circular grid, within diameter of dcir. Not used often.
 */
LOC_T *mkcirloc(double dcir, double dx){
    double rmax2,xmax2;
    int N,i,j,count;
    N=iceil(dcir/dx/2+5)*2;
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    loc->locx=malloc(sizeof(double)*N*N);
    loc->locy=malloc(sizeof(double)*N*N);
    rmax2=(dcir/dx/2.+1*1.414);
    rmax2*=rmax2;
    count=0;
   
    for(j=-N/2;j<N/2;j++){
	xmax2=rmax2-j*j;
	/*this method is good. there are no holes*/
	for(i=-N/2;i<N/2;i++){
	    if(i*i<xmax2){
		loc->locx[count]=i*dx;
		loc->locy[count]=j*dx;
		count++;
	    }
	}
    }
    loc->locx = realloc(loc->locx, sizeof(double)*count);
    loc->locy = realloc(loc->locy, sizeof(double)*count);
    loc->nloc = count;
    loc->dx   = dx;
    return loc;
}

/**
   Create LOC map when amplitude map is specified, within a maximum
   diameter. The amplitude map must be defined in the same sampling of dx.

   - cropmap=0: LOC is purely defined by the amplitude map. dcir is irrelevant in this case.
   - cropmap=1: LOC is forced to be within dcir. ampin is cropped.
   
   cstat records the starting point for each row to speed up accphi.  */
LOC_T* mkcirloc_amp(double** ampout,  /**<[out] amplitude map defined on loc*/
		    LOCSTAT_T *cstat, /**<[out] statistics on loc*/
		    MAP_T* ampin,     /**<[in/out] the amplitude that defines
					 the circular aperture, modified to
					 widhtin dtel if cropmap=1*/
		    double dcir,      /**<[in] maximum diameter of the circular loc*/
		    double dx,        /**<[in] sampling of output loc*/
		    int cropamp       /**<[in] do we zero ampin points if it outside of dcir*/
		    ){
    LOC_T *loc=calloc(1, sizeof(LOC_T));
    double *restrict ampp = ampin->p;
    int count,colcount,count0;
    double rmax2=(dcir/dx/2.+1*1.414);
    rmax2*=rmax2;

    if(ampin->nx*dx<dcir || ampin->ny*dx<dcir){
	error("maximum diameter is %g m, amplitude is %g mx%g m, "
	      "with sampling of 1/%g m\n",dcir,
	      ampin->nx*dx,ampin->ny*dx,1/dx);
    }else if(ampin->nx*dx>2*dcir){
	warning("maximum diameter is %g m, amplitude map is %gx%g m, "
		"with sampling of 1/%g m\n"
		"\t\t\tWill crop the amplitude map to maximum diameter of %g\n",
		dcir,ampin->nx*dx,ampin->ny*dx,1/dx,dcir);
    }
    
    loc->locx=malloc(sizeof(double)*ampin->nx*ampin->ny);
    loc->locy=malloc(sizeof(double)*ampin->nx*ampin->ny);
    *ampout  =malloc(sizeof(double)*ampin->nx*ampin->ny);
    count=0;   //count number of points.
    colcount=0;//count number of columns (y)

    int ncolstat=ampin->ny+1;
    if(cstat){
	cstat->cols=calloc(ncolstat, sizeof(LOCSTATCOL_T));
    }
    int offsety=ampin->ny/2;
    int offsetx=ampin->nx/2;
    /*scan through the input maps*/
    for(int iy=0; iy<ampin->ny; iy++){
	int jy=iy-offsety;
	int imin, imax, iminnotfound;
	double xmax2=rmax2-jy*jy;//max limit of x coord
	count0=count;//starting number of this column

	imin=INT_MAX;
	imax=-INT_MAX;
	iminnotfound=1;
	int offset2=iy*ampin->nx;
	/*we find nonzero min and max index for each column.*/
	for(int ix=0; ix<ampin->nx; ix++){
	    if(ampp[ix+offset2]>0){
		if(iminnotfound){
		    imin=ix;
		    iminnotfound=0;
		}
		imax=ix;
	    }
	}
	for(int ix=imin;ix<imax+1;ix++){
	    /*use amp only will create locs with holes it. the method
	      using stat to do accphi will be in trouble*/
	    int jx=ix-offsetx;
	    if(!cropamp || jx*jx<xmax2){
		loc->locx[count]=jx*dx;
		loc->locy[count]=jy*dx;
		(*ampout)[count]=ampp[ix+offset2];
		count++;
	    }else{
		ampp[ix+offset2]=0;/*zero out the amplitude map*/
	    }
	}
	if(cstat){
	    if(count>count0){
		/*row valid*/
		if(colcount>=ncolstat){
		    ncolstat*=2;
		    cstat->cols=realloc(cstat->cols, sizeof(LOCSTATCOL_T)*ncolstat);
		}
		(cstat->cols)[colcount].xstart=loc->locx[count0];
		(cstat->cols)[colcount].ystart=loc->locy[count0];
		(cstat->cols)[colcount].pos=count0;
		colcount++;
	    }
	}
    }
    /*append an addition row which marks the boundary*/
    if(cstat){
	cstat->cols[colcount].pos=count;/*record last*/
	cstat->cols = realloc(cstat->cols, sizeof(LOCSTATCOL_T)*(colcount+1));
	cstat->ncol = colcount;
	cstat->dx   = dx;
    }
    loc->locx = realloc(loc->locx, sizeof(double)*count);
    loc->locy = realloc(loc->locy, sizeof(double)*count);
    loc->nloc = count;
    loc->dx   = dx;
    *ampout=realloc(*ampout, sizeof(double)*count);
    return loc;
}


/**
   Find the point that closes to origin (0,0) in loc. Useful in single point
   piston constraint in creating reconstructor.
*/
int loccenter(LOC_T *loc){
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
dmat *loc_mcc_ptt(const LOC_T *loc, const double *amp){
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
dcell *pts_mcc_ptt(const PTS_T *pts, const double *amp){
    const int nmod=3;
    const int nsa=pts->nsa;
    dcell *mcc=dcellnew(nsa,1);
    for(int isa=0; isa<nsa; isa++){
	const double origy=pts->origy[isa];
	const double origx=pts->origx[isa];
	const double dx=pts->dx;
	const double *ampi=amp+pts->nx*pts->nx*isa;
	mcc->p[isa]=dnew(nmod,nmod);
	double (*ATA)[nmod]=(double(*)[nmod])(mcc->p[isa]->p);
	double a00=0,a01=0,a02=0,a11=0,a12=0,a22=0;
	for(int iy=0; iy<pts->nx; iy++){
	    double y=iy*dx+origy;
	    const double *ampx=ampi+iy*pts->nx;
	    for(int ix=0; ix<pts->nx; ix++){
		double x=ix*dx+origx;
		a00+=ampx[ix];
		a01+=ampx[ix]*x;
		a02+=ampx[ix]*y;
		a11+=ampx[ix]*x*x;
		a12+=ampx[ix]*x*y;
		a22+=ampx[ix]*y*y;
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
		  const LOC_T *loc, const double ipcc, 
		  const dmat *imcc, const double *amp, const double *opd){
    assert(imcc->nx==imcc->ny && imcc->nx==3);
    const long nloc=loc->nloc;
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;

    double tot=0;
    double coeff[3]={0,0,0};
    if(amp){
	if(rmsout){
	    for(long iloc=0; iloc<nloc; iloc++){
		const double junk=opd[iloc]*amp[iloc];
		coeff[0]+=junk;
		coeff[1]+=junk*locx[iloc];
		coeff[2]+=junk*locy[iloc];
		tot+=junk*opd[iloc];
	    }
	}else{
	    for(long iloc=0; iloc<nloc; iloc++){
		const double junk=opd[iloc]*amp[iloc];
		coeff[0]+=junk;
		coeff[1]+=junk*locx[iloc];
		coeff[2]+=junk*locy[iloc];
		//tot+=junk*opd[iloc];
	    }
	}
    }else{
	if(rmsout){
	    for(long iloc=0; iloc<nloc; iloc++){
		const double junk=opd[iloc];
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
		//tot+=junk*opd[iloc];
	    }
	}
    }
    if(coeffout){
	dmulvec3(coeffout, imcc, coeff);
    }
    if(rmsout){
	double pis=ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, imcc, coeff);
	rmsout[0]=tot-pis;//PR
	rmsout[1]=ptt-pis;//TT
	rmsout[2]=tot-ptt;//PTTR
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
void loc_remove_ptt(double *opd, const double *ptt, const LOC_T *loc){
    double ptt1[3];
    ptt1[0]=-ptt[0];
    ptt1[1]=-ptt[1];
    ptt1[2]=-ptt[2];
    loc_add_ptt(opd, ptt1, loc);
}

/**
   Add Piston/Tip/Tilt from OPD
*/
void loc_add_ptt(double *opd, const double *ptt, const LOC_T *loc){
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
void pts_ztilt(double *out, const PTS_T *pts, const dcell *imcc,
	       const double *amp, const double *opd){
    const int nsa=pts->nsa;
    assert(imcc->nx==nsa && imcc->ny==1);
    double (*res)[nsa]=(double(*)[nsa])out;
    for(int isa=0; isa<nsa; isa++){
	const double origy=pts->origy[isa];
	const double origx=pts->origx[isa];
	const double dx=pts->dx;
	const double *ampi=amp+pts->nx*pts->nx*isa;
	const double *opdi=opd+pts->nx*pts->nx*isa;
	assert(imcc->p[isa]->nx==3 && imcc->p[isa]->nx==3);
        double coeff[3]={0,0,0};
	double a0=0,a1=0,a2=0;
	for(int iy=0; iy<pts->nx; iy++){
	    double y=iy*dx+origy;
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
	res[0][isa]+=outp[1];
	res[1][isa]+=outp[2];
    }
}
/**
   Gather information about the starting of each column in loc.
 */
LOCSTAT_T *mklocstat(const LOC_T *loc){
    LOCSTAT_T *locstat=calloc(1, sizeof(LOCSTAT_T));
    const double *locx=loc->locx;
    const double *locy=loc->locy;
    int nloc=loc->nloc;
    double dx=locstat->dx=loc->dx;
    int ncolmax=(int)round((locy[nloc-1]-locy[0])/dx)+2;
    locstat->cols=malloc(ncolmax*sizeof(LOCSTATCOL_T));
    int colcount=0;
    int iloc=0;
    locstat->cols[colcount].pos=iloc;
    locstat->cols[colcount].xstart=locx[iloc];
    locstat->cols[colcount].ystart=locy[iloc];
    colcount++;
    for(iloc=1; iloc<loc->nloc; iloc++){
	if(fabs(locy[iloc]-locy[iloc-1])>1.e-12 
	   || fabs(locx[iloc]-locx[iloc-1]-dx)>1.e-12){
	    if(colcount>=ncolmax){
		ncolmax*=2;
		locstat->cols=realloc(locstat->cols, ncolmax*sizeof(LOCSTATCOL_T));
	    }
	    locstat->cols[colcount].pos=iloc;
	    locstat->cols[colcount].xstart=locx[iloc];
	    locstat->cols[colcount].ystart=locy[iloc];	    
	    colcount++;
	}
    }
    locstat->ncol=colcount;
    locstat->cols[colcount].pos=loc->nloc;
    locstat->cols[colcount].xstart=locx[loc->nloc-1];
    locstat->cols[colcount].ystart=locy[loc->nloc-1];
    locstat->cols=realloc(locstat->cols,(locstat->ncol+1)*sizeof(LOCSTATCOL_T));
    return locstat;
}
/**
   Create a gray pixel circular map in phi using coordinates defined in loc, center
defined using cx, cy, radius of r, and value of val */
void loccircle(double *phi,LOC_T *loc,double cx,double cy,double r,double val){
    //cx,cy,r are in unit of true unit, as in loc
    double dx=loc->dx;
    int nres=10;
    double cres=(nres-1.)/2.;
    double res=1./nres;
    double dxres=res*dx;
    double res2=res*res;
    double r2=r*r;
    double r2l=(r-dx*1.5)*(r-dx*1.5);
    double r2u=(r+dx*1.5)*(r+dx*1.5);
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
	    for(int jres=0; jres<10; jres++){
		iiy=y+(jres-cres)*dxres;
		double rr2y=(iiy-cy)*(iiy-cy);
		for(int ires=0; ires<nres; ires++){
		    iix=x+(ires-cres)*dxres;
		    double rr2r=(iix-cx)*(iix-cx)+rr2y;
		    if(rr2r<r2)
			tot++;
		}
	    }
	    phi[iloc]+=tot*res2*val;
	}
    }
}
/**
 Create a gray pixel elliptical map in phi using coordinates defined in loc,
center defined using cx, cy, radii of rx, ry, and value of val */
void locellipse(double *phi,LOC_T *loc,double cx,double cy,
		double rx,double ry,double val){
    //cx,cy,r are in unit of true unit, as in loc
    double dx=loc->dx;
    int nres=10;
    double cres=(nres-1.)/2.;
    double res=1./nres;
    double dxres=res*dx;
    double res2=res*res;
    double rx1=1./rx;
    double ry1=1./ry;
    double rx12=rx1*rx1;
    double ry12=ry1*ry1;
    double r2=2.;
    double r2l=(rx-dx)*(rx-dx)*rx12+(ry-dx)*(ry-dx)*ry12;
    double r2u=(rx+dx)*(rx+dx)*rx12+(ry+dx)*(ry+dx)*ry12;
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
		iiy=y+(jres-cres)*dxres;
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
   Remove uncoupled points in loc. debugged on 2009-12-20. Not used often. using
   a spcell, to compute the coupling which is modified accordingly.  */
void loc_reduce_spcell(LOC_T *loc, spcell *spc, int dim, int cont){
    int nloc=loc->nloc;
    int32_t *skip=calloc(nloc,sizeof(int32_t));
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
    if(cont){//make sure loc is continuous.
	LOCSTAT_T *locstat=mklocstat(loc);
	int ncol=locstat->ncol;
	for(int icol=0; icol<ncol; icol++){
	    int pstart=locstat->cols[icol].pos;
	    int pend=locstat->cols[icol+1].pos;
	    for(int pos=pstart; pos<pend; pos++){
		if(sum->p[pos]<1.e-200){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	    for(int pos=pend-1; pos>pstart-1; pos--){
		if(sum->p[pos]<1.e-200){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	}
	locstatfree(locstat);
    }else{
	for(int iloc=0; iloc<nloc; iloc++){
	    if(sum->p[iloc]<1.e-200){
		skip[iloc]=1;
	    }
	}
    }
    //writeint32(skip,nloc,1,"skip");
    int count=0;
    for(int iloc=0; iloc<nloc; iloc++){
	loc->locx[count]=loc->locx[iloc];
	loc->locy[count]=loc->locy[iloc];
	if(!skip[iloc]) count++;
    }
    loc->locx=realloc(loc->locx,sizeof(double)*count);
    loc->locy=realloc(loc->locy,sizeof(double)*count);
    loc->nloc=count;
    if(dim==1){//opd(loc)=sp*opd(locin);
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
    dfree(sum);
}

/**
   Remove uncoupled points in loc. debugged on 2009-12-20.  use single sparse
   matrix, which is modified accordingly.  */
void loc_reduce_sp(LOC_T *loc, dsp *sp, int dim, int cont){
    loc_free_map(loc);//remove the internal map before touchlong loc.
    int nloc=loc->nloc;
    if((dim==1 && nloc!=sp->m) || (dim==2 && nloc!=sp->n) || dim<0 || dim>2)
	error("Mismatch dimension\n");
    int32_t *skip=calloc(nloc,sizeof(int32_t));
    dmat* sum=spsumabs(sp,3-dim);
    if(cont){//make sure loc is continuous.
	LOCSTAT_T *locstat=mklocstat(loc);
	int ncol=locstat->ncol;
	for(int icol=0; icol<ncol; icol++){
	    int pstart=locstat->cols[icol].pos;
	    int pend=locstat->cols[icol+1].pos;
	    for(int pos=pstart; pos<pend; pos++){
		if(sum->p[pos]<1.e-200){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	    for(int pos=pend-1; pos>pstart-1; pos--){
		if(sum->p[pos]<1.e-200){
		    skip[pos]=1;
		}else{
		    break;
		}
	    }
	}
	locstatfree(locstat);
    }else{
	for(int iloc=0; iloc<nloc; iloc++){
	    if(sum->p[iloc]<1.e-200){
		skip[iloc]=1;
	    }
	}
    }
    writeint32(skip,nloc,1,"skip");
    int count=0;
    for(int iloc=0; iloc<nloc; iloc++){
	loc->locx[count]=loc->locx[iloc];
	loc->locy[count]=loc->locy[iloc];
	if(!skip[iloc]) count++;
    }
    loc->locx=realloc(loc->locx,sizeof(double)*count);
    loc->locy=realloc(loc->locy,sizeof(double)*count);
    loc->nloc=count;
    if(dim==1){//opd(loc)=sp*opd(locin);
	count=0;
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
    free(skip);
    dfree(sum);
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
	double coeff=pow(-1,s)*(double)factorial(ir-s)
	    /(double)(factorial(s)*factorial((ir+im)/2-s)*factorial((ir-im)/2-s));
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
dmat* loc_zernike(LOC_T *loc, double R, int nr){
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
	    //dwrite(Rnm,"Rnm_%d_%d",ir,im);
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
void loc_add_focus(double *opd, LOC_T *loc, double val){
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
dmat *loc2mat(LOC_T *loc,int piston){
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
   Convert PTS_T to LOC_T. Each entry in pts recors the lower
   left point in each subaperture. 
*/
LOC_T *pts2loc(PTS_T *pts){
    LOC_T*loc = calloc(1, sizeof(LOC_T));
    long nsa  = pts->nsa;
    long nx   = pts->nx;
    long nxsa = nx*nx;
    loc->locx = malloc(sizeof(double)*nsa*nxsa);
    loc->locy = malloc(sizeof(double)*nsa*nxsa);
    loc->nloc = nsa*nxsa;
    double dx = loc->dx = pts->dx;
    for(int isa=0; isa<nsa; isa++){
	const double origx = pts->origx[isa];
	const double origy = pts->origy[isa];
	double *locx = loc->locx+isa*nxsa;
	double *locy = loc->locy+isa*nxsa;
	for(int jj=0; jj<nx; jj++){
	    const double yy=origy+(double)jj*dx;
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
void locrot(LOC_T *loc, const double theta){
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
   duplicate a LOC_T object; */
LOC_T *locdup(LOC_T *loc){
    LOC_T *new=calloc(1,sizeof(LOC_T));
    new->locx=malloc(loc->nloc*sizeof(double));
    new->locy=malloc(loc->nloc*sizeof(double));
    new->nloc=loc->nloc;
    new->dx=loc->dx;
    memcpy(new->locx,loc->locx,sizeof(double)*loc->nloc);
    memcpy(new->locy,loc->locy,sizeof(double)*loc->nloc);
    return new;
}

/**
   Transform coordinate. loc contains x,y; locm contains xm, ym.
   xm{ip}=\{sum}_{ic}(coeff[0](1,ic)*pow(x,coeff[0](2,ic))*pow(y,coeff[0](3,ic)))
   ym{ip}=\{sum}_{ic}(coeff[1](1,ic)*pow(x,coeff[1](2,ic))*pow(y,coeff[1](3,ic)))
*/
LOC_T *loctransform(LOC_T *loc, dmat **coeff){

    if(coeff[0]->nx!=3 || coeff[1]->nx!=3){
	error("Coeff is in wrong format\n");
    }
    PDMAT(coeff[0],cx);
    PDMAT(coeff[1],cy);

    LOC_T *locm=calloc(1, sizeof(LOC_T));
    locm->nloc=loc->nloc;
    locm->dx=loc->dx;
    locm->locx=calloc(locm->nloc, sizeof(double));
    locm->locy=calloc(locm->nloc, sizeof(double));
    const double *restrict x=loc->locx;
    const double *restrict y=loc->locy;
    double *restrict xm=locm->locx;
    double *restrict ym=locm->locy;

    for(long iloc=0; iloc<loc->nloc; iloc++){
	for(long ic=0; ic<coeff[0]->ny; ic++){
	    xm[iloc]+=cx[ic][0]*pow(x[iloc],cx[ic][1])*pow(y[iloc],cx[ic][2]);
	}
	for(long ic=0; ic<coeff[1]->ny; ic++){
	    ym[iloc]+=cy[ic][0]*pow(x[iloc],cy[ic][1])*pow(y[iloc],cy[ic][2]);
	}
    }
    return locm;
}

/**
   Compute the size of a map that is used to build loc or can fully contain loc
   (embed)
*/
void loc_nxny(long *nx, long *ny, const LOC_T *loc){
    double xmax, xmin, ymax, ymin;
    maxmindbl(loc->locx, loc->nloc, &xmax, &xmin, NULL);
    maxmindbl(loc->locy, loc->nloc, &ymax, &ymin, NULL);
    *nx=(long)round((xmax-xmin)/loc->dx)+1;
    *ny=(long)round((ymax-ymin)/loc->dx)+1;
}
