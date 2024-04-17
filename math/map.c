/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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



#include "map.h"
#define check_nonempty(A) (A?(A->nx?(A->ny?1:0):0):0)
static void map_keywords(map_t *map){
	if(!map || map->keywords) return;
	char keywords[1024];
	snprintf(keywords, sizeof(keywords), "dx=%g; dy=%g; ox=%g; oy=%g;", map->dx, map->dy, map->ox, map->oy);
	map->keywords=strdup(keywords);
}
/**
   create a new map_t object.
*/
map_t* mapnew(long nx, long ny, real dx, real dy){
	map_t* map=mycalloc(1, map_t);
	dinit((dmat**)&map, nx, ny);
	map->h=0;
	map->dx=dx;
	map->dy=dy;
	map->ox=-map->nx/2*map->dx;
	map->oy=-map->ny/2*map->dy;
	map->vx=0;
	map->vy=0;
	map->iac=0;
	map_keywords(map);
	return map;
}
/**
   ceate a new map_t object from existing one. P is left empty.
*/
map_t* mapnew2(map_t* A){
	if(!check_nonempty(A)) return NULL;
	map_t* map=mycalloc(1, map_t);
	dinit((dmat**)&map, A->nx, A->ny);
	map->h=A->h;
	map->dx=A->dx;
	map->dy=A->dy;
	map->ox=A->ox;
	map->oy=A->oy;
	map->vx=A->vx;
	map->vy=A->vy;
	map->iac=A->iac;
	map_keywords(map);
	return map;
}
map_t* mapref(map_t* in){
	if(!check_nonempty(in)) return NULL;
	map_t* out=mycalloc(1, map_t);
	memcpy(out, in, sizeof(map_t));
	if(in->keywords) out->keywords=strdup(in->keywords);
	out->mem=mem_ref(in->mem);
	out->fp=NULL;
	out->fft=NULL;
	return out;
}

/**
   Create a circular aperture on map_t.
*/
void mapcircle(map_t* map, real r, real val){
	if(!check_nonempty(map)) return;
	dcircle((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r, val);
}

/**
   Create a circular aperture on map_t.
*/
void mapcircle_symbolic(map_t* map, real r){
	if(!check_nonempty(map)) return;
	dcircle_symbolic((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r);
}

/**
   Find the inner and outer diameter of an amplitude map contained in map_t.
*/
void map_d_din(map_t* map, real* d, real* din){
	if(!check_nonempty(map)||!d||!din) return;
	real r2min=INFINITY, r2max=0;
	for(long iy=0; iy<map->ny; iy++){
		real y=iy*map->dy+map->oy;
		for(long ix=0; ix<map->nx; ix++){
			if(P(map, ix, iy)>EPS){
				real x=ix*map->dx+map->ox;
				real r2=x*x+y*y;
				if(r2>r2max) r2max=r2;
				if(r2<r2min) r2min=r2;
			}
		}
	}
	*d=sqrt(r2max)*2;
	*din=sqrt(r2min)*2;
}

/**
   create a metapupil map, with size nx*ny, origin at (ox,oy), sampling of dx, dy,
   height of ht, that can cover all the directions specified in dirs

   offset: distance in pixel from the point closest to the origin to origin (right side).
   0: there is a point on the origin.
   1/2: the closest point to origin is 1/2 pixel.

   pad!=0: round nx, ny to power of 2.  */

void create_metapupil(map_t** mapout,/**<[out] map*/
	long* nxout,  /**<[out] nx*/
	long* nyout,  /**<[out] ny*/
	dmat* dirs,   /**<[in] All Directions to cover (thetax, thetay, hs, hc)*/
	real D,     /**<[in] Diameter (meter)*/
	real ht0,   /**<[in] Conjugation Height (meter)*/
	real dx,    /**<[in] Sampling along x (meter)*/
	real dy,    /**<[in] Sampling along y (meter)*/
	real offset,/**<[in] Fractional offset of point closet from origin. between [0, 1)*/
	real guard, /**<[in] Width of guard area, in meter*/
	long ninx,    /**<[in] Suggested size along x*/
	long niny,    /**<[in] Suggested size along y*/
	int pad,      /**<[in] Increase nx, ny to power of 2*/
	int square    /**<[in] Full square/rectangular grid*/
){
	if(!dirs) return;
	const real R=D/2;
	real minx=INFINITY, miny=INFINITY, maxx=-INFINITY, maxy=-INFINITY;
	if(dirs->nx<4||dirs->ny<=0){
		error("dirs should have no less than 3 rows and positive number of cols.\n");
	}
	for(int idir=0; idir<dirs->ny; idir++){
		real ht=ht0-P(dirs,3, idir);//hc
		real RR=(1.-ht/P(dirs, 2, idir))*R+guard;
		real sx1=(P(dirs, 0, idir)*ht)-RR;
		real sx2=(P(dirs, 0, idir)*ht)+RR;
		real sy1=(P(dirs, 1, idir)*ht)-RR;
		real sy2=(P(dirs, 1, idir)*ht)+RR;
		//Need to work when ht<0;
		if(sx1<minx) minx=sx1;
		if(sx1>maxx) maxx=sx1;
		if(sx2<minx) minx=sx2;
		if(sx2>maxx) maxx=sx2;
		if(sy1<miny) miny=sy1;
		if(sy1>maxy) maxy=sy1;
		if(sy2<miny) miny=sy2;
		if(sy2>maxy) maxy=sy2;
	}
	if(square){//if square, also make symmetric.
		maxx=MAX(fabs(minx), fabs(maxx));
		minx=-maxx;
		maxy=MAX(fabs(miny), fabs(maxy));
		miny=-maxy;
	}
	/*ajust central point offset*/
	{
		offset=offset-floor(offset);//between 0 and 1
		real mind=minx/dx;
		real adjust=mind-floor(mind)-offset;
		if(adjust<0){
			adjust++;
		}
		minx-=adjust*dx;
		mind=miny/dy;
		adjust=mind-floor(mind)-offset;
		if(adjust<0){
			adjust++;
		}
		miny-=adjust*dy;
	}
	real ox=minx;
	real oy=miny;
	long nx=ceil((maxx-ox)/dx)+1;
	long ny=ceil((maxy-oy)/dy)+1;

	/*Make it square */
	if(square){
		ny=nx=(nx<ny)?ny:nx;
	}
	if(pad){/*pad to power of 2 */
		ninx=1<<iceil(log2((real)nx));
		niny=1<<iceil(log2((real)ny));
	}
	if(ninx>1){
		if(ninx<nx) warning("ninx=%ld is too small. need %ld\n", ninx, nx);
		ox=ox-(ninx-nx)/2*dx;
		nx=ninx;
	}
	if(niny>1){
		if(niny<ny)  warning("niny=%ld is too small. need %ld\n", niny, ny);
		oy=oy-(niny-ny)/2*dy;
		ny=niny;
	}
	if(nxout)
		*nxout=nx;
	if(nyout)
		*nyout=ny;
	if(mapout){
		*mapout=mapnew(nx, ny, dx, dy);
		(*mapout)->ox=ox;
		(*mapout)->oy=oy;
		(*mapout)->h=ht0;
		dmat* dmap=(dmat*)(*mapout);
		if(square){/**Only want square grid*/
			dset(dmap, 1);
		} else{/*Want non square grid*/
			for(int idir=0; idir<dirs->ny; idir++){
				real ht=ht0-P(dirs, 3, idir);
				real sx=-ox+(P(dirs, 0, idir)*ht);
				real sy=-oy+(P(dirs, 1, idir)*ht);
				real RR=R*(1.-ht/P(dirs, 2, idir))+guard;
				dcircle_symbolic(dmap, sx, sy, dx, dy, RR);
			}
			for(int i=0; i<nx*ny; i++){
				P(dmap,i)=(P(dmap,i))>1.e-15?1:0;
			}
		}
	}
}

/**
   convert a dmat to map_t.

   keywords are used to specify parameters:
   - the sampling (dx or D; strong recommended; default to 1/64; also specify dy if spacing is different.),
   - layer heiht (default to 0),
   - frozen flow speed (vx, vy)
   - origin (ox, oy; will use nx/2*dx or ny/2*dy if not specified)
   - for DM grid only: inter-actuator-coupling (iac)
*/
map_t* d2map(const dmat* in){
	if(!check_nonempty(in)) return NULL;
	dmat *tmp=dref(in);
	map_t* map=myrealloc(tmp, 1, map_t);
	memset((char*)map+sizeof(dmat), 0, sizeof(map_t)-sizeof(dmat));
	char* keywords=in->keywords;
	map->iac=0;
	map->ox=search_keyword_num(keywords, "ox");
	map->oy=search_keyword_num(keywords, "oy");
	map->dx=search_keyword_num(keywords, "dx");
	map->dy=search_keyword_num(keywords, "dy");
	map->h=search_keyword_num_default(keywords, "h", 0);
	map->vx=search_keyword_num_default(keywords, "vx", 0);
	map->vy=search_keyword_num_default(keywords, "vy", 0);
	real D=search_keyword_num(keywords, "D");
	real offset=search_keyword_num_default(keywords, "offset",0);

	if(isnan(map->dx)){
		if(isnan(D)){
			map->dx=1./64.;
			warning_once("dx and D are not specified. Set dx, dy to 1/%g.\n", 1./map->dx);
		}else{
			map->dx=D/map->nx;
			dbg2("dx is not specified, but D is. Set dx, dy to 1/%g.\n", 1/map->dx);
		}
	}
	if(isnan(map->dy)){
		map->dy=map->dx;
	}
	if(isnan(map->ox)){
		map->ox=(-map->nx/2+offset)*map->dx;
		dbg2("ox is not specified. set to %g\n", map->ox);
	}
	if(isnan(map->oy)){
		map->oy=(-map->ny/2+offset)*map->dy;
		dbg2("oy is not specified. set to %g\n", map->oy);
	}
	return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
mapcell* dcell2map(const dcell* in){
	if(!check_nonempty(in)) return NULL;
	mapcell* map=(mapcell*)cellnew(in->nx, in->ny);
	for(long i=0; i<in->nx*in->ny; i++){
		if(!P(in,i)->keywords&&in->keywords){
			P(in,i)->keywords=strdup(in->keywords);
		}
		P(map,i)=d2map(P(in,i));
	}
	return map;
}


/**
    convert a dmat to rmap_t.

  	keywords are used to specify parameters:
    - the sampling (dx or D; strong recommended; default to 1/64; also specify dy if spacing is different.),
    - origin (ox, oy; will use nx/2*dx or ny/2*dy if not specified)
    - txdeg, dydeg: tilt of surface in degree wrt beam. only 1 can be less than 90 degree.
	- ftel:  focal length of the telescope
	- fexit: distance from exit pupil to focus
	- fsurf  distance from surface to focu
*/
rmap_t* d2rmap(const dmat* in){
	if(!check_nonempty(in)) return NULL;
	dmat *tmp=dref(in);
	rmap_t* map=myrealloc(tmp, 1, rmap_t);
	memset((char*)map+sizeof(dmat), 0, sizeof(rmap_t)-sizeof(dmat));
	if(!in->keywords){
		error("this dmat has no header\n");
	}
	const char *keywords=in->keywords;
	map->ox=search_keyword_num(keywords, "ox");
	map->oy=search_keyword_num(keywords, "oy");
	map->dx=search_keyword_num(keywords, "dx");
	map->dy=search_keyword_num(keywords, "dy");
	map->txdeg=search_keyword_num(keywords, "txdeg");
	map->tydeg=search_keyword_num(keywords, "tydeg");
	map->ftel=search_keyword_num(keywords, "ftel");
	map->fexit=search_keyword_num(keywords, "fexit");
	map->fsurf=search_keyword_num(keywords, "fsurf");
	if(isnan(map->ox)||isnan(map->oy)||isnan(map->dx)||isnan(map->dy)||isnan(map->fsurf)){
		error("header is needed to convert dmat to map_t\n");
	}
	return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
rmap_t** dcell2rmap(int* nlayer, const dcell* in){
	if(!check_nonempty(in)) return NULL;
	*nlayer=in->nx*in->ny;
	rmap_t** map=mycalloc(in->nx*in->ny, rmap_t*);
	for(long i=0; i<in->nx*in->ny; i++){
		map[i]=d2rmap(P(in,i));
	}
	return map;
}


/**
	convert fields to keyword string
*/
void map_header(map_t* map){
	if(map&&!map->keywords){
		char keywords[1024];
		snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
			map->ox, map->oy, map->dx, map->dy, map->h, map->vx, map->vy);
		map->keywords=strdup(keywords);
	}
}

/**
 * convert fields to keyword string
*/
void rmap_header(rmap_t* map){
	if(map&&!map->keywords){
		char keywords[1024];
		snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\ntxdeg=%.15g\ntydeg=%.15g\nftel=%.15g\nfexit=%.15g\nfsurf=%.15g\n",
			map->ox, map->oy, map->dx, map->dy, map->txdeg, map->tydeg, map->ftel, map->fexit, map->fsurf);
		map->keywords=strdup(keywords);
	}
}
