/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define check_map(A) (A?(A->nx?(A->ny?1:0):0):0)
/**
   create a new map_t object.
*/
map_t* mapnew(long nx, long ny, real dx, real dy){
	dmat *tmp=dnew_do(nx, ny, NULL, 0);
	map_t* map=myrealloc(tmp, 1, map_t);
	map->h=0;
	map->dx=dx;
	map->dy=dy;
	map->ox=-map->nx/2*map->dx;
	map->oy=-map->ny/2*map->dy;
	map->vx=0;
	map->vy=0;
	map->iac=0;
	return map;
}
/**
   ceate a new map_t object from existing one. P is left empty.
*/
map_t* mapnew2(map_t* A){
	if(!check_map(A)) return NULL;
	dmat *tmp=dnew_do(A->nx, A->ny, NULL, 0);
	map_t* map=myrealloc(tmp, 1, map_t);
	map->h=A->h;
	map->dx=A->dx;
	map->dy=A->dy;
	map->ox=A->ox;
	map->oy=A->oy;
	map->vx=A->vx;
	map->vy=A->vy;
	map->iac=A->iac;
	return map;
}
map_t* mapref(map_t* in){
	if(!check_map(in)) return NULL;
	map_t* out=mycalloc(1, map_t);
	memcpy(out, in, sizeof(map_t));
	out->mem=mem_ref(in->mem);
	return out;
}

/**
   Create a circular aperture on map_t.
*/
void mapcircle(map_t* map, real r, real val){
	if(!check_map(map)) return;
	dcircle((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r, val);
}
/**
   Create a circular aperture on map_t.
*/
void mapcircle_symbolic(map_t* map, real r){
	if(!check_map(map)) return;
	dcircle_symbolic((dmat*)map, (-map->ox), (-map->oy), map->dx, map->dy, r);
}


/**
   Find the inner and outer diameter of an amplitude map contained in map_t.
*/
void map_d_din(map_t* map, real* d, real* din){
	if(!check_map(map)||!d||!din) return;
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
*/
map_t* d2map(const dmat* in){
	if(!check_map(in)) return NULL;
	dmat *tmp=dref(in);
	map_t* map=myrealloc(tmp, 1, map_t);
	memset((char*)map+sizeof(dmat), 0, sizeof(map_t)-sizeof(dmat));
	char* header=in->header;
	map->iac=0;
	map->ox=search_header_num(header, "ox");
	map->oy=search_header_num(header, "oy");
	map->dx=search_header_num(header, "dx");
	map->dy=search_header_num(header, "dy");
	map->h=search_header_num(header, "h");
	map->vx=search_header_num(header, "vx");
	map->vy=search_header_num(header, "vy");
	if(isnan(map->dx)){
		warning_once("dx is not specified in header, set to 1/64.\n");
		map->dx=1./64.;
	}
	if(isnan(map->dy)){
		map->dy=map->dx;
	}
	if(isnan(map->ox)||isnan(map->oy)){
		map->ox=-map->nx/2*map->dx;
		map->oy=-map->ny/2*map->dy;
	}
	if(isnan(map->h)) map->h=0;
	if(isnan(map->vx)) map->vx=0;
	if(isnan(map->vy)) map->vy=0;
	return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
mapcell* dcell2map(const dcell* in){
	if(!check_map(in)) return NULL;
	mapcell* map=(mapcell*)cellnew(in->nx, in->ny);
	for(long i=0; i<in->nx*in->ny; i++){
		if(!P(in,i)->header&&in->header){
			P(in,i)->header=strdup(in->header);
		}
		P(map,i)=d2map(P(in,i));
	}
	return map;
}


/**
   convert a dmat to rmap_t.
*/
rmap_t* d2rmap(const dmat* in){
	if(!check_map(in)) return NULL;
	dmat *tmp=dref(in);
	rmap_t* map=myrealloc(tmp, 1, rmap_t);
	memset((char*)map+sizeof(dmat), 0, sizeof(rmap_t)-sizeof(dmat));
	char* header=in->header;
	if(!in->header){
		error("this dmat has no header\n");
	}
	map->ox=search_header_num(header, "ox");
	map->oy=search_header_num(header, "oy");
	map->dx=search_header_num(header, "dx");
	map->dy=search_header_num(header, "dy");
	map->txdeg=search_header_num(header, "txdeg");
	map->tydeg=search_header_num(header, "tydeg");
	map->ftel=search_header_num(header, "ftel");
	map->fexit=search_header_num(header, "fexit");
	map->fsurf=search_header_num(header, "fsurf");
	if(isnan(map->ox)||isnan(map->oy)||isnan(map->dx)||isnan(map->dy)||isnan(map->fsurf)){
		error("header is needed to convert dmat to map_t\n");
	}
	return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
rmap_t** dcell2rmap(int* nlayer, const dcell* in){
	if(!check_map(in)) return NULL;
	*nlayer=in->nx*in->ny;
	rmap_t** map=mycalloc(in->nx*in->ny, rmap_t*);
	for(long i=0; i<in->nx*in->ny; i++){
		map[i]=d2rmap(P(in,i));
	}
	return map;
}


//void mapwritedata(file_t *fp, map_t *map){
void map_header(map_t* map){
	if(map&&!map->header){
		char header[1024];
		snprintf(header, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
			map->ox, map->oy, map->dx, map->dy, map->h, map->vx, map->vy);
		map->header=strdup(header);
	}
}

//void rmapwritedata(file_t *fp, rmap_t *map){
void rmap_header(rmap_t* map){
	if(map&&!map->header){
		char header[1024];
		snprintf(header, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\ntxdeg=%.15g\ntydeg=%.15g\nftel=%.15g\nfexit=%.15g\nfsurf=%.15g\n",
			map->ox, map->oy, map->dx, map->dy, map->txdeg, map->tydeg, map->ftel, map->fexit, map->fsurf);
		map->header=strdup(header);
	}
}
