/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
 * @brief Convert a dmat to map object. Ensure that the type is properly initialized
 * 
 * @param A 	input dmat object. will be changed in place.
 * @return map_t* 
 */
map_t* map_convert(dmat *A){
	if(!A || A->id!=M_REAL) return NULL;
	map_t* map=myrecalloc(A, dmat, map_t);
	map->id=M_MAP;
	map->make_keywords=map_make_keywords;
	if(map->keywords) map_parse_keywords(map); 
	return map;
}

/**
 * @brief Convert a dmat to rmap object. Ensure that the type is properly initialized
 * 
 * @param A 	input dmat object. will be changed in place.
 * @return map_t* 
 */
rmap_t* rmap_convert(dmat* A){
	if(!A || A->id!=M_REAL) return NULL;
	rmap_t* map=myrecalloc(A, dmat, rmap_t);
	map->id=M_RMAP;
	map->make_keywords=rmap_make_keywords;
	if(map->keywords) rmap_parse_keywords(map);
	return map;
}
/**
   create a new map_t object.
*/
map_t* mapnew(long nx, long ny, real dx, real dy){
	map_t* map=map_convert(dnew(nx, ny));
	map->ht=0;
	map->dx=dx;
	map->dy=dy;
	map->ox=-(map->nx/2)*map->dx;
	map->oy=-(map->ny/2)*map->dy;
	map->vx=0;
	map->vy=0;
	map->iac=0;
	map->dratio=0;
	return map;
}
/**
   ceate a new map_t object from existing one. P is left empty.
*/
map_t* mapnew2(map_t* A){
	if(isempty(A)) return NULL;
	map_t* map=mapnew(A->nx, A->ny, A->dx, A->dy);
	map->ht=A->ht;
	map->ox=A->ox;
	map->oy=A->oy;
	map->vx=A->vx;
	map->vy=A->vy;
	map->iac=A->iac;
	map->dratio=A->dratio;
	return map;
}
map_t* mapref(const map_t* in){
	if(isempty(in)) return NULL;
	map_t* out=map_convert(dref(in->dmat));
	memcpy((char*)out+sizeof(dmat), (char*)in+sizeof(dmat), sizeof(map_t)-sizeof(dmat));
	return out;
}
map_t* mapdup(const map_t* in){
	if(isempty(in)) return NULL;
	map_t* out=map_convert(ddup(in->dmat));
	memcpy((char*)out+sizeof(dmat), (char*)in+sizeof(dmat), sizeof(map_t)-sizeof(dmat));
	return out;
}

/**
   Create a circular aperture on map_t.
*/
void mapcircle(map_t* map, real r, real val){
	if(isempty(map)) return;
	dcircle(DMAT(map), (-map->ox), (-map->oy), map->dx, map->dy, r, val);
}

/**
   Find the inner and outer diameter of an amplitude map contained in map_t.
*/
void map_d_din(map_t* map, real* d, real* din){
	if(isempty(map)||!d||!din) return;
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
	const real R=D/2;
	real ox=0, oy=0;
	long nx, ny;
	{
		real minx=-R-guard, miny=-R-guard, maxx=R+guard, maxy=R+guard;
		if(dirs && ht0!=0){
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
			//info("range is %g, %g; %g, %g\n", minx, maxx, miny, maxy);
		}else if(ht0!=0){
			error("if ht0 is 0, dirs must be set\n");
		}
		ox=floor(minx/dx)*dx;//round
		nx=ceil((maxx-ox)/dx)+1;
		oy=floor(miny/dx)*dx;
		ny=ceil((maxy-oy)/dy)+1;
		//info("ox=%g, oy=%g, nx=%ld, ny=%ld\n", ox, oy, nx, ny);
	}
	if(square){//make it a square grid.
		ox=oy=MIN(ox,oy);
		nx=ny=MAX(nx,ny);
	}
	/*ajust central point offset*/
	offset=offset-floor(offset);//between 0 and 1
	if(offset){
		ox-=offset*dx; nx++;
		oy-=offset*dy; ny++;
	}
	if(pad){/*pad to power of 2 */
		long nx2=1<<iceil(log2((real)nx));
		long ny2=1<<iceil(log2((real)ny));
		ox-=(nx2-nx)/2*dx;
		oy-=(ny2-ny)/2*dx;
		nx=nx2;
		ny=ny2;
	}
	if(ninx){
		if(ninx<nx) warning("ninx=%ld is too small. need %ld\n", ninx, nx);
		ox-=(ninx-nx)/2*dx;
		nx=ninx;
	}
	if(niny){
		if(niny<ny)  warning("niny=%ld is too small. need %ld\n", niny, ny);
		oy-=(niny-ny)/2*dy;
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
		(*mapout)->ht=ht0;
		dmat* dmap=DMAT((*mapout));
		if(square){/**Only want square grid*/
			dset(dmap, 1);
		} else{/*Want non square grid*/
			if(dirs){
				for(int idir=0; idir<dirs->ny; idir++){
					real ht=ht0-P(dirs, 3, idir);
					real sx=-ox+(P(dirs, 0, idir)*ht);
					real sy=-oy+(P(dirs, 1, idir)*ht);
					real RR=R*(1.-ht/P(dirs, 2, idir))+guard;
					dcircle(dmap, sx, sy, dx, dy, RR, 1);
				}
			}else{
				dcircle(dmap, -ox, -oy, dx, dy, R, 1);
			}
			for(int i=0; i<nx*ny; i++){
				P(dmap,i)=(P(dmap,i))>1.e-15?1:0;
			}
		}
	}
}

/**
   parse keywords for map_t. 
   keywords are used to specify parameters:
   - the sampling (dx or D; strong recommended; default to 1/64; also specify dy if spacing is different.),
   - layer heiht (default to 0),
   - frozen flow speed (vx, vy)
   - origin (ox, oy; will use nx/2*dx or ny/2*dy if not specified)
   - for DM grid only: inter-actuator-coupling (iac)
*/
void map_parse_keywords(map_t* map){
	if(!map || map->id!=M_MAP || !map->keywords) return;
	const char* keywords=map->keywords;
	map->ox=search_keyword_num(keywords, "ox");
	map->oy=search_keyword_num(keywords, "oy");
	map->dx=search_keyword_num(keywords, "dx");
	map->dy=search_keyword_num(keywords, "dy");
	map->ht=search_keyword_num_default(keywords, "h", 0);
	map->vx=search_keyword_num_default(keywords, "vx", 0);
	map->vy=search_keyword_num_default(keywords, "vy", 0);
	map->iac=search_keyword_num_default(keywords, "iac", 0);
	map->dratio=search_keyword_num_default(keywords, "dratio", 0);
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
}

/**
    Parse keywords for rmap_t.

  	keywords are used to specify parameters:
    - the sampling (dx or D; strong recommended; default to 1/64; also specify dy if spacing is different.),
    - origin (ox, oy; will use nx/2*dx or ny/2*dy if not specified)
    - txdeg, dydeg: tilt of surface in degree wrt beam. only 1 can be less than 90 degree.
	- ftel:  focal length of the telescope
	- fexit: distance from exit pupil to focus
	- fsurf  distance from surface to focu
*/
void rmap_parse_keywords(rmap_t* map){
	if(!map || map->id!=M_RMAP || !map->keywords) return;
	const char *keywords=map->keywords;
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
		error("valid header is needed to construct rmap_t\n");
	}
}

/**
	convert fields to keyword string
*/
void map_make_keywords(cell* p){
	if(!p || p->id!=M_MAP) return;
	map_t *map=(map_t *)p;
	if(map->keywords) free(map->keywords);
	char keywords[1024];
	snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		map->ox, map->oy, map->dx, map->dy, map->ht, map->vx, map->vy);
	map->keywords=strdup(keywords);
}

/**
 * convert fields to keyword string
*/
void rmap_make_keywords(cell* p){
	if(!p || p->id!=M_RMAP) return;
	rmap_t* map=(rmap_t*)p;
	if(map->keywords) free(map->keywords);
	char keywords[1024];
	snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\ntxdeg=%.15g\ntydeg=%.15g\nftel=%.15g\nfexit=%.15g\nfsurf=%.15g\n",
		map->ox, map->oy, map->dx, map->dy, map->txdeg, map->tydeg, map->ftel, map->fexit, map->fsurf);
	map->keywords=strdup(keywords);
}

/**
   blend atm1 with atm2 according to wind direction angle and required
   overlapping region of at least overx*overy. Not used.
*/

void map_blend(map_t* atm1, map_t* atm2, long overx, long overy){
	const long nx=NX(atm1);
	const long ny=NY(atm1);
	int ca=0;
	if(atm1->vx>EPS){
		ca=-1;/*reverse sign of vx */
	} else if(atm1->vx<-EPS){
		ca=1;
	}
	int sa=0;
	if(atm1->vy>EPS){
		sa=-1;/*reverse sign of vy */
	} else if(atm1->vy<-EPS){
		sa=1;
	}
	long rr;
	long offx=nx-overx;
	long offy=(ny-overy)*nx;
	if(ca==0){/*along y. */
		rr=ny-overy;/*distance between the origins. */
		atm2->oy=atm1->oy+rr*sa*atm1->dx;
		atm2->ox=atm1->ox;
		real wty=sa<0?1:0;
		real* p1=P(atm1)+(1-(long)wty)*offy;
		real* p2=P(atm2)+(long)wty*offy;
		real overyd=(real)overy;
		for(long iy=0; iy<overy; iy++){
			real wt1=fabs(wty-(real)(iy+1)/overyd);
			for(long ix=0; ix<nx; ix++){
				p1[ix+iy*nx]=(1-wt1)*p1[ix+iy*nx]+wt1*p2[ix+iy*nx];
				p2[ix+iy*nx]=p1[ix+iy*nx];
			}
		}
	} else if(sa==0){
		rr=nx-overx;/*distance between the origins. */
		atm2->ox=atm1->ox+rr*ca*atm1->dx;
		atm2->oy=atm1->oy;
		real wtx=ca<0?1:0;
		real* p1=P(atm1)+(1-(long)wtx)*offx;
		real* p2=P(atm2)+(long)wtx*offx;
		real wts[overx];
		real overxd=(real)overx;
		for(long ix=0; ix<overx; ix++){
			wts[ix]=fabs(wtx-(real)(ix+1)/overxd);
		}
		for(long iy=0; iy<ny; iy++){
			for(long ix=0; ix<overx; ix++){
				p1[ix+iy*nx]=(1-wts[ix])*p1[ix+iy*nx]+wts[ix]*p2[ix+iy*nx];
				p2[ix+iy*nx]=p1[ix+iy*nx];
			}
		}
	} else{
		error("We do not support this wind direction: ca=%d, sa=%d\n", ca, sa);
	}
}
