/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "../math/mathdef.h"
#include "loc.h"
#include "map.h"


/**
   Free the MAP in loc_t

   @param[in,out] loc The loc whose map is freed and set to NULL.
*/
void loc_free_map(loc_t* loc){
	if(loc&&loc->map){
		mapfree(loc->map);
		loc->map=NULL;
	}
}
/**
   Free the stat in loc_t

   @param[in,out] loc The loc whose stat is freed and set to NULL.
*/
void loc_free_stat(loc_t* loc){
	if(loc&&loc->stat){
		free(loc->stat->cols);
		free(loc->stat);
		loc->stat=NULL;
	}
}
/**
   Free the internal data (stat and map) of loc_t. Used as the deinit callback
   of loc_t; the coordinate memory is freed by the dmat machinery.

   @param[in,out] p The loc_t (as cell). Ignored if not a loc_t.
*/
void locfree_do(cell* p){
	if(p && p->id==M_LOC){
		loc_t *loc=(loc_t*)p;
		loc_free_stat(loc);
		loc_free_map(loc);
	}
}
/**
   Convert dmat to loc in place without parsing keywords. dx and dy are set to NAN.

   @param[in,out] A nloc*2 matrix of coordinates. It is reallocated in place and must no longer be used.
   @return The converted loc_t, which takes over A; NULL if A is NULL.
*/
static loc_t *dmat2loc_inplace(dmat*A){
	if(!A) return NULL;
	loc_t* loc=myrecalloc(A, dmat, loc_t);
	loc->id=M_LOC;
	loc->locy=loc->locx+loc->nloc;
	loc->dmat->deinit=locfree_do;
	loc->dmat->make_keywords=loc_keywords;
	loc->dx=NAN;
	loc->dy=NAN;
	loc->ht=0;
	loc->iac=0;
	loc->dratio=0;
	loc->aoi=0;
	return loc;
}
/**
   Parse dx, dy, ht, iac, dratio, aoi from the keywords of loc and determine
   dx, dy from the data if needed.

   @param[in,out] loc The loc whose fields are set.
*/
static void loc_parse_header(loc_t* loc){
	if(!loc) return;
	loc->dx=search_keyword_num(loc->keywords, "dx");
	loc->dy=search_keyword_num_default(loc->keywords, "dy", loc->dx);
	loc->ht=search_keyword_num_default(loc->keywords, "ht", 0);
	loc->iac=search_keyword_num_default(loc->keywords, "iac", 0);
	loc->dratio=search_keyword_num_default(loc->keywords, "dratio", 0);
	loc->aoi=search_keyword_num_default(loc->keywords, "aoi", 0);
	loc_dxdy(loc);
}
/**
   Convert a dmat in place to loc_t object. Keywords are parsed into the fields.

   @param[in,out] A nloc*2 input array of coordinates. It is reallocated in place and must no longer be used.
   @return The converted loc_t, which takes over A; NULL if A is NULL or not a real matrix.
*/
loc_t* loc_convert(dmat* A){
	if(!A || A->id!=M_REAL) return NULL;
	loc_t *loc=dmat2loc_inplace(A);
	loc_parse_header(loc);
	return loc;
}
/**
   loc_convert() with the signature required by cellconvert().

   @param[in,out] A Input dmat object as cell*. See loc_convert().
   @return The converted loc_t as cell*.
*/
cell* loc_convert_cell(cell* A){
	return (cell*)loc_convert(dmat_cast(A));
}
/**
   Create a loc with nloc elements.

   @param[in] nloc Number of points.
   @param[in] dx Sampling along x (m).
   @param[in] dy Sampling along y (m). Uses dx if 0.
   @return Newly allocated loc with zero coordinates; the caller frees it with locfree().
*/
loc_t* locnew(long nloc, real dx, real dy){
	loc_t* loc=dmat2loc_inplace(dnew(nloc, 2));
	loc->dx=dx;
	loc->dy=dy?dy:dx;
	return loc;
}
/**
   Reference an existing loc. map and stat are not referenced.

   @param[in] in The loc to reference.
   @return Newly allocated loc sharing coordinates with in; the caller frees it with locfree(). NULL if in is NULL.
*/
loc_t* locref(loc_t* in){
	if(!in) return NULL;
	loc_t* loc=dmat2loc_inplace(dref(in->dmat));
	loc->dx=in->dx;
	loc->dy=in->dy;
	loc->ht=in->ht;
	loc->iac=in->iac;
	loc->dratio=in->dratio;
	loc->aoi=in->aoi;
	return loc;
}
/**
   Create a pts with nsa, dsa, nx, dx

   @param[in] nsa Number of subapertures.
   @param[in] dsax Subaperture size along x (m).
   @param[in] dsay Subaperture size along y (m).
   @param[in] nxsa Number of points per subaperture along x.
   @param[in] nysa Number of points per subaperture along y.
   @param[in] dx Sampling of points within a subaperture along x (m).
   @param[in] dy Sampling of points within a subaperture along y (m).
   @return Newly allocated pts; the caller frees it with ptsfree().
*/
pts_t* ptsnew(long nsa, real dsax, real dsay, long nxsa, long nysa, real dx, real dy){
	pts_t* pts=myrecalloc(locnew(nsa, dsax, dsay), loc_t, pts_t);
	pts->nxsa=nxsa;
	pts->nysa=nysa;
	pts->dx=dx;
	pts->dy=dy;
	if(dx!=dy || dsax != dsay || nxsa!=nysa){
		warning("please double check non-square grid\n");
	}
	return pts;
}
/**
   apply hashlittle() to loc_t.

   @param[in] loc The loc to hash.
   @param[in] key Initial hash key.
   @return Updated hash key; key unchanged if loc is NULL.
*/
uint32_t lochash(const loc_t* loc, uint32_t key){
	if(loc){
		key=dhash(loc->dmat, key);
	}
	return key;
}
/**
   Create an vector to embed OPD into a square array for FFT purpose. 
   This is reverse operation from loc_create_map which maps 2d array indexing info loc indexing.

   @param[in,out] nembed Size of square array. Computed and set if <=0 on input, otherwise used as is.
   @param[in] loc Irregular grid.
   @param[in] oversize Ratio of oversizing. 2 for FFT to achieve Nyquist sampling.
   @param[in] fftpad Pad to size that is fast for FFT.
   @return Newly allocated vector of the same length of loc containing 1d index into the square array; the caller frees it.
*/
lmat* loc_create_embed(long* nembed, const loc_t* loc, real oversize, int fftpad){
	if(!loc) return NULL;
	real xmin, xmax, ymin, ymax;
	dvecmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
	dvecmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
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
	} else{
		if(*nembed<(long)(nxy*0.6)){
			error("Supplied nembed %ld is too small, recommend %ld\n", *nembed, nxy);
			return NULL;
		} else if(*nembed<nxy){
			warning("Supplied nembed %ld is too small, recommend %ld\n", *nembed, nxy);
		}
		nxy=*nembed;
	}
	xmin-=(nxy-nx+1)/2*loc->dx;//odd -> pad on left
	ymin-=(nxy-ny+1)/2*loc->dy;
	lmat* embed=lnew(loc->nloc, 1);
	for(int iloc=0; iloc<loc->nloc; iloc++){
		long ix=(long)round((loc->locx[iloc]-xmin)*dx_in1);
		long iy=(long)round((loc->locy[iloc]-ymin)*dy_in1);
		P(embed,iloc)=ix+iy*nxy;
	}
	return embed;
}
/**
   Create a map for loc without padding. See loc_create_map_npad().

   @param[in] loc The loc. Its map is created (the struct is modified despite const).
*/
void loc_create_map(const loc_t* loc){
	loc_create_map_npad((loc_t*)loc, 0, 0, 0);
}
PNEW(maplock);
DEF_ENV_FLAG(LOC_MAP_EXTEND, 1);

/**
   Create a map for loc so that we can obtain the index in loc by x,y
   coordinate. Useful in ray tracing (accphi.c). If extend is true, the invalid
   region in the map is replaced with closest valid point. This extends the
   coverage of the grid automatically.

   The map stores loc index plus 1 at valid points, and negative index plus 1
   of the neighbor for extended points. Does nothing if the map already exists.

   @param[in] loc The loc. Its map and npad are set (the struct is modified despite const).
   @param[in] npad Number of points to pad on each side.
   @param[in] nx Size of the map along x. 0 to determine automatically.
   @param[in] ny Size of the map along y. 0 to determine automatically.
*/
void loc_create_map_npad(const loc_t* loc, int npad, int nx, int ny){
	if(!loc) return;
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
	real xmin, xmax, ymin, ymax;
	dvecmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
	dvecmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
	/*padding the map. normally don't need. */
	if(npad>0){
		xmin-=npad*fabs(loc->dx);
		ymin-=npad*fabs(loc->dy);
		xmax+=npad*fabs(loc->dx);
		ymax+=npad*fabs(loc->dy);
	}
	const real dx_in1=1./fabs(loc->dx);
	const real dy_in1=1./fabs(loc->dy);
	long map_nx=(long)round((xmax-xmin)*dx_in1)+1;
	long map_ny=(long)round((ymax-ymin)*dy_in1)+1;
	if(nx&&ny){
		if(map_nx>nx||map_ny>ny){
			error("Specified size %dx%d is too small, need at least %ldx%ld\n",
				nx, ny, map_nx, map_ny);
			UNLOCK(maplock);
			return;
		}
		//the guaranteed padding on all sides is the requested npad plus the minimum extra padding along x and y.
		npad+=(int)MIN((nx-map_nx)/2, (ny-map_ny)/2);
		xmin-=((nx-map_nx)/2)*loc->dx;
		ymin-=((ny-map_ny)/2)*loc->dy;
		map_nx=nx;
		map_ny=ny;
	}
	((loc_t*)loc)->npad=npad;/*record the information. */
	((loc_t*)loc)->map=mapnew(map_nx, map_ny, loc->dx, loc->dy);
	loc->map->iac=loc->iac;
	loc->map->dratio=loc->dratio;
	loc->map->aoi=loc->aoi;
	loc->map->ox=xmin;
	loc->map->oy=ymin;
	dmat* pmap=DMAT(loc->map);
	const real* locx=loc->locx;
	const real* locy=loc->locy;
	for(long iloc=0; iloc<loc->nloc; iloc++){
		int ix=(int)round((locx[iloc]-xmin)*dx_in1);
		int iy=(int)round((locy[iloc]-ymin)*dy_in1);
		P(pmap, ix, iy)=iloc+1;/*start from 1. */
	}
	if(LOC_MAP_EXTEND&&loc->nloc<map_nx*map_ny){
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
		lmat* level=lnew(map_nx, map_ny);
		for(int iy=0; iy<map_ny; iy++){
			for(int ix=0; ix<map_nx; ix++){
				if(P(pmap, ix, iy)){
					P(level, ix, iy)=0;//level of existence.
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
						if(!P(pmap, ix, iy)){
							int count=0;//count how many neighbors of lowest level
							/*choose the minimum level of interpolation*/
							int min_level=INT_MAX, min_jx=0, min_jy=0, min_sep=INT_MAX;
							real min_rad=INFINITY;
							for(int jy=MAX(0, iy-1); jy<MIN(map_ny, iy+2); jy++){
								for(int jx=MAX(0, ix-1); jx<MIN(map_nx, ix+2); jx++){
									int sep=(abs(iy-jy)+abs(ix-jx));
									int iphi;
									if((sep==1||sep==2)&&(iphi=fabs(P(pmap, jx, jy)))){//edge
										iphi--;
										real rad=locx[iphi]*locx[iphi]+locy[iphi]*locy[iphi];//dist^2 to center
										if(P(level, jx, jy)<min_level){
											min_level=(int)P(level, jx, jy);
											min_rad=rad;
											min_sep=sep;
											min_jy=jy;
											min_jx=jx;
											count=1;//reset to one
										} else if(P(level, jx, jy)==min_level){
											count++;
											if(sep<min_sep){
												min_sep=sep;
												min_jy=jy;
												min_jx=jx;
											} else if(sep==min_sep&&rad<min_rad){
												min_rad=rad;
												min_jy=jy;
												min_jx=jx;
											}
										}
									}
								}
							}
							if((min_level==cur_level&&count==3)||min_level<cur_level){
							//if level satisfy or even if three neighbor with higher level.
								P(pmap, ix, iy)=-fabs(P(pmap, min_jx, min_jy));
								P(level, ix, iy)=min_level+1;
								num_found++;
							}
						}
					}
				}
			} while(num_found);
		}
		lfree(level);
	}
	UNLOCK(maplock);
}
#define CALC_EMBED_OFFSET(doffx,moffx,mx,dnx,mnx)\
if(dnx>=mnx){/*output array is larger*/\
	doffx=(dnx-mnx)>>1;\
	moffx=0;\
	mx=mnx;\
} else{/*map is larger*/\
	doffx=0;\
	moffx=(mnx-dnx)>>1;\
	mx=dnx;\
}
/**
  Embed in into dest according to map defined in loc->map. The two arrays are assumed to be concentric.

  @param[in,out] _dest Destination 2-d array. Resized to the map size without padding if empty. Points not on loc are set to NAN unless add is set.
  @param[in] loc The loc. Its map is created if not yet.
  @param[in] _in Values defined on loc.
  @param[in] expanded Also fill the extended (ghost) points of the map.
  @param[in] add Add to dest instead of overwriting.
*/
void loc_embed(anydmat _dest, const loc_t* loc, const_anydmat _in, int expanded, int add){
	if(!loc||!_dest.dm||!_in.dm) return;
	if(!loc->map){
		loc_create_map((loc_t*)loc);
	}
	map_t *map=loc->map;
	dmat *dest=dmat_cast(_dest.dm);
	dmat *in=dmat_cast(_in.dm);
	if(!dest || !in){
		error("Needs implementation for these types: %u %u", _dest.dm?_dest.dm->id:0, _in.dm?_in.dm->id:0);
		return;
	}
	if(!dest->nx||!dest->ny){
		dresize(dest, map->nx-loc->npad*2, map->ny-loc->npad*2);
	}
	long doffx,doffy,moffx,moffy,mx,my;
	CALC_EMBED_OFFSET(doffx, moffx, mx, dest->nx, map->nx);
	CALC_EMBED_OFFSET(doffy, moffy, my, dest->ny, map->ny);
	if(!add) dset(dest, NAN);
	const real *pin=P(in)-1;//iphi count from 1
	for(int iy=0; iy<my; iy++){
		for(int ix=0; ix<mx; ix++){
			long iphi=P(map, ix+moffx, iy+moffy);
			if(expanded && iphi<0) iphi=-iphi;
			if(iphi>0){
				P(dest, ix+doffx, iy+doffy)=add?(P(dest, ix+doffx, iy+doffy)+pin[iphi]):pin[iphi];
			}
		}
	}
}
/**
  A convenient wrapper for loc_embed to be called by matlab or python
  @param[in] loc The location coordinate. Its map is created if not yet.
  @param[in] arr The points defined on loc to be embeded to 2-d array. May have multiple vectors.
  @return Newly allocated cell of embeded 2-d arrays, one per vector; the caller frees it.
*/
dcell* loc_embed2(const loc_t* loc, const dmat* arr){
	if(!loc||!arr) return NULL;
	if(!loc->map){
		loc_create_map((loc_t*)loc);
	}
	long arrnx=arr->nx;
	long arrny=arr->ny;
	if(arr->nx==1&&arr->ny!=1){
		arrnx=arr->ny;
		arrny=1;
	}
	long nx=arrnx/loc->nloc;
	long ny=arrny;
	if(nx*loc->nloc!=arrnx){
		error("arr has wrong dimension: %ldx%ld. loc length is %ld \n", arrnx, arrny, loc->nloc);
		return NULL;
	}
	dcell* dest=dcellnew(nx, ny);
	dmat *tmp=dnew_do(loc->nloc,1,P(arr),0);//weak reference
	for(long ix=0; ix<nx*ny; ix++){
		P(dest, ix)=dnew(loc->map->nx, loc->map->ny);
		tmp->p=&P(arr, ix*loc->nloc);
		loc_embed(P(dest, ix), loc, tmp, 0, 0);
	}
	dfree(tmp);
	return dest;
}

/**
  Embed in into dest according to map defined in loc->map for cells. Arbitraty map cannot
  be used here for mapping.

  @param[in,out] pdest Destination cell array of the map size. Allocated if NULL. Each element references the corresponding element of in, or is NULL.
  @param[in] loc The loc. Its map is created if not yet.
  @param[in] in Cell array defined on loc.
*/
void loc_embed_cell(dcell** pdest, const loc_t* loc, const dcell* in){
	if(!pdest||!loc||!in) return;
	if(!loc->map){
		loc_create_map((loc_t*)loc);
	}
	map_t *map=loc->map;
	if(!*pdest){
		*pdest=dcellnew(map->nx, map->ny);
	}
	dcell* dest=*pdest;
	if(dest->nx!=map->nx||dest->ny!=map->ny){
		error("dest and map doesn't agree\n");
		return;
	}
	for(long i=0; i<map->nx*map->ny; i++){
		long iphi=P(map,i);
		if(iphi>0){//do not embed for extra padding.
			P(dest,i)=dref(P(in, iphi-1));
		} else{
			P(dest,i)=0;
		}
	}
}
/**
  Extract values from in into dest according to map defined in loc->map.
  Arbitraty map cannot be used here for mapping.

  @param[out] dest Values on loc. Must be preallocated with nloc elements.
  @param[in] loc The loc. Its map must exist.
  @param[in] in 2-d array with the same dimension as loc->map.
*/
void loc_extract(dmat* dest, const loc_t* loc, map_t* in){
	if(!dest||!loc||!in) return;
	map_t* map=loc->map;
	if(!map){
		error("map is null\n");
		return;
	}
	if(in->nx!=map->nx||in->ny!=map->ny){
		error("in and map doesn't agree: in is %ldx%ld, map is %ldx%ld\n",
			in->nx, in->ny, map->nx, map->ny);
		return;
	}
	real* restrict pdest1=P(dest)-1;//iphi count from 1
	for(long i=0; i<map->nx*map->ny; i++){
		long iphi=P(map,i);
		if(iphi>0){
			pdest1[iphi]=P(in,i);
		}
	}
}
/**
   Convert a map to a loc that collects all positive entries while keep each row continuous

   @param[in] map The map.
   @param[in] thres Threshold relative to the maximum of map. Points above it are valid.
   @return Newly allocated loc; the caller frees it with locfree(). NULL if map is NULL.
*/
loc_t* loc_from_map(map_t* map, real thres){
	if(!map) return NULL;
	const real dx=map->dx;
	const real dy=map->dy;
	const real ox=map->ox;
	const real oy=map->oy;
	const long nx=map->nx;
	const long ny=map->ny;
	if(thres){
		thres*=dmax(map->dmat);
	}
	long ix, iy;
	loc_t* loc=locnew(nx*ny, dx, dy);
	long count=0;
	for(iy=0; iy<ny; iy++){
		long istart=-1,iend=-1;
		for(ix=0; ix<nx; ix++){
			if(P(map, ix, iy)>thres){
				istart=ix;
				break;
			}
		}
		for(ix=nx-1; ix>-1; ix--){
			if(P(map, ix, iy)>thres){
				iend=ix;
				break;
			}
		}
		if(istart!=-1){
			for(ix=istart; ix<=iend; ix++){
				loc->locx[count]=ix*dx+ox;
				loc->locy[count]=iy*dy+oy;
				count++;
			}
		}
	}
	if(!count){
		warning("loc_from_map: there are 0 points.\n");
		print_backtrace();
	}
	locresize(loc, count);
	loc->iac=map->iac;
	loc->dratio=map->dratio;
	loc->aoi=map->aoi;
	loc->ht=map->ht;
	return loc;
}
/**
   Create 1 dimensional loc with given vector.

   @param[in] x Coordinates along x.
   @param[in] nx Number of points.
   @return Newly allocated loc with y set to 0; the caller frees it with locfree(). NULL if x is NULL or nx is 0.
*/
loc_t* mk1dloc_vec(real* x, long nx){
	if(!x||!nx) return NULL;
	real dx=fabs((x[nx-1]-x[0])/(nx-1));
	loc_t* loc=locnew(nx, dx, dx);
	memcpy(loc->locx, x, sizeof(real)*nx);
	return loc;
}
/**
   Create 1 dimensional loc with origin at x0, sampling of dx, and nx numbers.

   @param[in] x0 Origin (m).
   @param[in] dx Sampling (m).
   @param[in] nx Number of points.
   @return Newly allocated loc with y set to 0; the caller frees it with locfree().
*/
loc_t* mk1dloc(real x0, real dx, long nx){
	loc_t* loc=locnew(nx, dx, dx);
	for(long ix=0; ix<nx; ix++){
		loc->locx[ix]=x0+dx*ix;
	}
	return loc;
}
/**
   Create a loc array that covers the map_t. Notice that it is different from
   loc_from_map which only covers valid (value>0) regions.

   @param[in] map The map.
   @return Newly allocated loc; the caller frees it with locfree(). NULL if map is NULL.
*/
loc_t* mksqloc_map(const map_t* map){
	if(!map) return NULL;
	return mksqloc(map->nx, map->ny, map->dx, map->dy, map->ox, map->oy);
}
/**
   Create a loc array of size nx*ny at sampling dx with origin at (nx/2, ny/2)

   @param[in] nx Number of points along x.
   @param[in] ny Number of points along y.
   @param[in] dx Sampling along x (m).
   @param[in] dy Sampling along y (m).
   @return Newly allocated loc; the caller frees it with locfree().
*/
loc_t* mksqloc_auto(long nx, long ny, real dx, real dy){
	return mksqloc(nx, ny, dx, dy, -nx/2*dx, -ny/2*dy);
}
/**
   Create a loc array contains coorindates in a square map of size nx*ny, with
   sampling dx, and at origin ox,oy

   @param[in] nx Number of points along x.
   @param[in] ny Number of points along y.
   @param[in] dx Sampling along x (m).
   @param[in] dy Sampling along y (m).
   @param[in] ox Coordinate of the first point along x (m).
   @param[in] oy Coordinate of the first point along y (m).
   @return Newly allocated loc; the caller frees it with locfree().
*/
loc_t* mksqloc(long nx, long ny, real dx, real dy, real ox, real oy){
	loc_t* loc=locnew(nx*ny, dx, dy);
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

   @param[in] D Outer diameter (m).
   @param[in] Din Inner diameter (m).
   @param[in] dx Sampling (m).
   @param[in] thres Threshold of the gray pixel amplitude, relative to maximum, for a point to be included.
   @return Newly allocated loc; the caller frees it with locfree().
*/
loc_t* mkannloc(real D, real Din, real dx, real thres){
	long nx=MAX(D/dx+1, 1);
	map_t* xy=mapnew(nx, nx, dx, dx);
	mapcircle(xy, D/2, 1);
	if(Din>0&&Din<D){
		mapcircle(xy, Din/2, -1);
	}
	loc_t* loc=loc_from_map(xy, thres);
	mapfree(xy);
	return loc;
}

/**
   Estimate the diameter of loc.
   2021-09-03: corrected calculation for square grid. Use the inscribed diameter but not circumcircle.

   @param[in] loc The loc.
   @return Maximum extent along x or y (m); 0 if loc is NULL.
*/
real loc_diam(const loc_t* loc){
	if(!loc) return 0;
	real xmin, xmax, ymin, ymax;
	dvecmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
	dvecmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
	return MAX(xmax-xmin, ymax-ymin);
}
/**
   Find the point that closes to origin (0,0) in loc. Useful in single point
   piston constraint in creating reconstructor.

   @param[in] loc The loc.
   @return Index of the point closest to origin; 0 if loc is NULL.
*/
int loccenter(const loc_t* loc){
	if(!loc) return 0;
	int ipix, jpix=0;
	real r2, r2min=INFINITY;

	for(ipix=0; ipix<loc->nloc; ipix++){
		r2=pow(loc->locx[ipix], 2)+pow(loc->locy[ipix], 2);
		if(r2<r2min){
			jpix=ipix;
			r2min=r2;
		}
	}
	return jpix;
}

/**
   assemble modes of piston/tip/tilt into Nx3 matrix M, and compute their
   inverse covariance matrix.  mcc=M'*(amp.*M)

   @param[in] loc The loc.
   @param[in] amp Amplitude weighting. NULL for uniform.
   @return Newly allocated 3x3 matrix mcc (not inverted); the caller frees it. NULL if loc is NULL.
*/
dmat* loc_mcc_ptt(const loc_t* loc, const real* amp){
	if(!loc) return NULL;
	const int nmod=3;
	real* mod[nmod];
	dmat* mcc=dnew(nmod, nmod);
	mod[0]=NULL;
	mod[1]=loc->locx;
	mod[2]=loc->locy;
	for(int jmod=0; jmod<nmod; jmod++){
		for(int imod=jmod; imod<nmod; imod++){
			real tmp=dvecdot(mod[imod], mod[jmod], amp, loc->nloc);
			P(mcc, imod, jmod)=P(mcc, jmod, imod)=tmp;
		}
	}
	return mcc;
}
/**
   same as loc_mcc_ptt, except using pts instead of loc.
   this is used to build imcc for TT/F powfs ztilt imcc

   @param[in] pts The subaperture points.
   @param[in] amp Amplitude weighting, nxsa*nysa per subaperture.
   @return Newly allocated cell of 3x3 matrices, one per subaperture (not inverted); the caller frees it.
*/
dcell* pts_mcc_ptt(const pts_t* pts, const real* amp){
	if(!pts) return NULL;
	const int nmod=3;
	const long nsa=pts->nsa;
	dcell* mcc=dcellnew_same(nsa, 1, nmod, nmod);
	for(long isa=0; isa<nsa; isa++){
		const real origy=pts->origy[isa];
		const real origx=pts->origx[isa];
		const real dx=pts->dx;
		const real dy=pts->dy;
		const real* ampi=amp+pts->nxsa*pts->nysa*isa;
		dmat* ATA=P(mcc,isa);
		double a00=0, a01=0, a02=0, a11=0, a12=0, a22=0;//use double for accumulation
		for(int iy=0; iy<pts->nysa; iy++){
			real y=iy*dy+origy;
			const real* ampx=ampi+iy*pts->nxsa;
			for(int ix=0; ix<pts->nxsa; ix++){
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
		P(ATA, 0, 0)=a00;
		P(ATA, 1, 1)=a11;
		P(ATA, 2, 2)=a22;
		P(ATA, 0, 1)=P(ATA, 1, 0)=a01;
		P(ATA, 0, 2)=P(ATA, 2, 0)=a02;
		P(ATA, 1, 2)=P(ATA, 2, 1)=a12;
	}
	return mcc;
}

/**
   evaluate piston/tip-tilt/ removed wavefront error.
   output coeffout in unit of radian like units.

   @param[out] rmsout Variance of 3 elements: PR (piston removed), TT (tip/tilt), and PTTR (piston/tip/tilt removed). Skipped if NULL.
   @param[out] coeffout Piston, tip, and tilt coefficients. Skipped if NULL.
   @param[in] loc The loc.
   @param[in] ipcc Inverse of piston mode weight. Computed internally if imcc is NULL.
   @param[in] imcc Inverse of 3x3 mode covariance from loc_mcc_ptt(). Computed internally if NULL.
   @param[in] amp Amplitude weighting. NULL for uniform.
   @param[in] opd OPD on loc.
   @param[in] stride Stride between consecutive OPD elements.
*/
void loc_calc_ptt_stride(real* rmsout, real* coeffout, const loc_t* loc, real ipcc, const dmat* imcc, const real* amp,
	const real* opd, long stride){
	if(!loc||!opd) return;
	dmat *imcc_internal=NULL;
	if(imcc){
		assert(imcc->nx==imcc->ny&&imcc->nx==3);
	}else{
		imcc=imcc_internal=loc_mcc_ptt(loc, amp);
		ipcc=1./P(imcc_internal,0,0);
		dsvd_pow(imcc_internal, -1);
	}
	const long nloc=loc->nloc;
	const real* restrict locx=loc->locx;
	const real* restrict locy=loc->locy;

	real tot=0;
	double coeff[3]={0,0,0};//use double for accumulation
	if(amp){
		long iopd=0;
		for(long iloc=0; iloc<nloc; iloc++){
			const real junk=opd[iopd]*amp[iloc];
			coeff[0]+=junk;
			coeff[1]+=junk*locx[iloc];
			coeff[2]+=junk*locy[iloc];
			tot+=junk*opd[iopd];
			iopd+=stride;
		}
	} else{
		long iopd=0;
		for(long iloc=0; iloc<nloc; iloc++){
			const real junk=opd[iopd];
			coeff[0]+=junk;
			coeff[1]+=junk*locx[iloc];
			coeff[2]+=junk*locy[iloc];
			tot+=junk*opd[iopd];
			iopd+=stride;
		}
	}
	real coeff2[3]={coeff[0],coeff[1],coeff[2]};
	if(coeffout){
		dmulvec3(coeffout, imcc, coeff2);
	}
	if(rmsout){
		real pis=ipcc*coeff[0]*coeff[0];/*piston mode variance */
		real ptt=dwdot(coeff2, imcc, coeff2);/*p/t/t mode variance. */
		rmsout[0]=tot-pis;/*PR */
		rmsout[1]=ptt-pis;/*TT */
		rmsout[2]=tot-ptt;/*PTTR */
	}
	dfree(imcc_internal);
}
/**
   Same as loc_calc_ptt_stride() with stride 1. rmsout[0..2] receives PR, TT,
   and PTTR variance; coeffout (optional) receives the piston/tip/tilt
   coefficients. ipcc and imcc are computed internally if imcc is NULL.

   @param[out] rmsout Variance of PR, TT, and PTTR. Skipped if NULL.
   @param[out] coeffout Piston, tip, and tilt coefficients. Skipped if NULL.
   @param[in] loc The loc.
   @param[in] ipcc Inverse of piston mode weight.
   @param[in] imcc Inverse of 3x3 mode covariance. Computed internally if NULL.
   @param[in] amp Amplitude weighting. NULL for uniform.
   @param[in] opd OPD on loc.
*/
void loc_calc_ptt(real *rmsout, real *coeffout,
	const loc_t *loc, const real ipcc,
	const dmat *imcc, const real *amp, const real *opd){
	loc_calc_ptt_stride(rmsout, coeffout, loc, ipcc, imcc, amp, opd, 1);
}
/**
   Calculate variance of OPD in modes defined in mod.

   Notice the difference in rmsout for loc_calc_ptt and this loc_calc_mode.
   in loc_calc_ptt, out[0] is PR, out[1] is TT, out[2] is PTTR
   in loc_calc_mod, out[0] is PR, out[1] is PTR, out[2] is PTTR, etc

   The mod is already orth-noramlized so that we have the following simple forms.

   coeffout is in unit of zernike!

   @param[out] rmsout Residual variance after removing the first imod+1 modes, one per mode.
   @param[out] coeffout Mode coefficients. Skipped if NULL.
   @param[in] mod Orthonormal modes, nloc*nmod.
   @param[in] amp Amplitude weighting. Must not be NULL.
   @param[in] opd OPD on loc.
*/
void loc_calc_mod(real* rmsout, real* coeffout, const dmat* mod, const real* amp, real* opd){
	if(!mod||!opd||!amp) return;
	const long nmod=mod->ny;
	const long nloc=mod->nx;
	double tot=0; //use double for accumulation
	double val[nmod];
	memset(val, 0, sizeof(double)*nmod);
	//val=(Mod*Amp*Opd)
	for(long iloc=0; iloc<nloc; iloc++){
		real junk=opd[iloc]*amp[iloc];
		tot+=opd[iloc]*junk;
		for(long imod=0; imod<nmod; imod++){
			val[imod]+=P(mod, iloc, imod)*junk;
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
   Subtract Piston/Tip/Tilt (in radian) from OPD

   @param[in,out] opd OPD on loc.
   @param[in] ptt Piston, tip, and tilt coefficients.
   @param[in] loc The loc.
*/
void loc_sub_ptt(dmat* opd, const real* ptt, const loc_t* loc){
	if(!opd||!ptt||!loc) return;
	real ptt1[3];
	ptt1[0]=-ptt[0];
	ptt1[1]=-ptt[1];
	ptt1[2]=-ptt[2];
	loc_add_ptt(opd, ptt1, loc);
}
/**
   Add Piston/Tip/Tilt to OPD

   @param[in,out] opd OPD on loc.
   @param[in] stride Stride between consecutive OPD elements.
   @param[in] ptt Piston, tip, and tilt coefficients.
   @param[in] loc The loc.
*/
void loc_add_ptt_stride(real *opd, long stride, const real *ptt, const loc_t *loc){
	if(!opd||!ptt||!loc) return;
	const long nloc=loc->nloc;
	const real *restrict locx=loc->locx;
	const real *restrict locy=loc->locy;
	OMP_TASK_FOR(4)
	for(long iloc=0; iloc<nloc; iloc++){
		opd[iloc*stride]+=ptt[0]+ptt[1]*locx[iloc]+ptt[2]*locy[iloc];
	}
}
/**
   Add Piston/Tip/Tilt to OPD

   @param[in,out] opd OPD on loc.
   @param[in] ptt Piston, tip, and tilt coefficients.
   @param[in] loc The loc.
*/
void loc_add_ptt(dmat* opd, const real* ptt, const loc_t* loc){
	if(!opd||!ptt||!loc) return;
	if(loc->nloc!=opd->nx*opd->ny){
		error("Invalid dimensions. loc has %ld, opd has %ldx%ld\n", loc->nloc, opd->nx, opd->ny);
		return;
	}
	loc_add_ptt_stride(P(opd), 1, ptt, loc);
	/*const long nloc=loc->nloc;
	const real* restrict locx=loc->locx;
	const real* restrict locy=loc->locy;
OMP_TASK_FOR(4)
	for(long iloc=0; iloc<nloc; iloc++){
		P(opd,iloc)+=ptt[0]+ptt[1]*locx[iloc]+ptt[2]*locy[iloc];
	}*/
}
/**
   Project piston/tip/tilt from opd which may have multiple independent columns.
   If cov is set, opd should be symmetric (covariance) and ptt is removed from left and right side.

   @param[in,out] opd OPD with nloc rows. Modified in place.
   @param[in] loc The loc.
   @param[in] amp Amplitude weighting. NULL for uniform.
   @param[in] imcc Inverse of 3x3 mode covariance. Computed internally if NULL.
   @param[in] cov Also remove from the rows, treating opd as a covariance matrix.
*/
void loc_remove_ptt(dmat *opd, const loc_t *loc, const real *amp, const dmat *imcc, int cov){
	if(NX(opd)!=loc->nloc){
		error("loc_remove_ptt: opd should have rows equal to loc->nloc and amplitude map size.\n");
		return;
	}
	real ptt[3];
	dmat *imcc_internal=NULL;
	if(!imcc){
		imcc=imcc_internal=loc_mcc_ptt(loc, amp);
		dsvd_pow(imcc_internal, -1);//more resilient to numerical errors than invspd
	}
	for(int ic=0; ic<NY(opd); ic++){
		loc_calc_ptt_stride(NULL, ptt, loc, 0, imcc, amp, PCOL(opd, ic),1);
		ptt[0]=-ptt[0]; ptt[1]=-ptt[1]; ptt[2]=-ptt[2];
		loc_add_ptt_stride(PCOL(opd,ic), 1, ptt, loc);
	}
	if(cov){
		if(NX(opd)!=NY(opd)){
			error("when cov is set, opd shall be symmetric\n");
		}else{
			for(int ic=0; ic<NX(opd); ic++){
				loc_calc_ptt_stride(NULL, ptt, loc, 0, imcc, amp, &P(opd, ic, 0), NX(opd));
				ptt[0]=-ptt[0]; ptt[1]=-ptt[1]; ptt[2]=-ptt[2];
				loc_add_ptt_stride(&P(opd, ic, 0), NX(opd), ptt, loc);
			}
		}
	}
	if(imcc_internal) dfree(imcc_internal);
}

/**
   Remove focus mode from gradients. Do not consider noise weighting

   @param[in,out] grad Gradients (all x then all y) on saloc. Focus is removed if factor is not 0.
   @param[in] saloc Subaperture locations.
   @param[in] factor Fraction of focus to remove. 0 to only estimate.
   @return Focus coefficient a of a*(x^2+y^2) estimated from the input gradients.
*/
real loc_remove_focus_grad(dmat *grad, const loc_t *saloc, real factor){
	const dmat *mode=saloc->dmat;
	real d1=dvecdot(P(mode), P(grad), NULL, PN(grad));
	real d2=dvecdot(P(mode), P(mode), NULL, PN(grad));
	real focus=d1/d2;
	if(factor) dadd(&grad, 1, mode, -focus*factor);
	//multiply by 0.5 because focus has gradient of 2*locx, 2*locy
	return focus*0.5;
}
/**
   Compute zernike best fit for all subapertures. add result to out.  returns
   radians not zernike modes of tip/tilt. Used in wfsgrad

   @param[in,out] out Gradients, x for all subapertures followed by y. Allocated if NULL; otherwise accumulated into.
   @param[in] pts The subaperture points.
   @param[in] imcc Inverse of 3x3 mode covariance for each subaperture.
   @param[in] amp Amplitude, nxsa*nysa per subaperture.
   @param[in] opd OPD, nxsa*nysa per subaperture.
*/
void pts_ztilt(dmat** out, const pts_t* pts, const dcell* imcc, const real* amp, const real* opd){
	if(!pts||!imcc||!opd) return;
	const long nsa=pts->nsa;
	assert(imcc->nx==nsa&&imcc->ny==1);
	if(!*out) *out=dnew(nsa*2, 1);
	real* res=P(*out);
	for(long isa=0; isa<nsa; isa++){
		const real origy=pts->origy[isa];
		const real origx=pts->origx[isa];
		const real dx=pts->dx;
		const real dy=pts->dy;
		const real* ampi=amp+pts->nxsa*pts->nysa*isa;
		const real* opdi=opd+pts->nxsa*pts->nysa*isa;
		assert(P(imcc,isa)->nx==3&&P(imcc,isa)->ny==3);
		real coeff[3]={0,0,0};
		double a0=0, a1=0, a2=0;
		for(int iy=0; iy<pts->nysa; iy++){
			real y=iy*dy+origy;
			const real* ampx=ampi+iy*pts->nxsa;
			const real* opdx=opdi+iy*pts->nxsa;
			for(int ix=0; ix<pts->nxsa; ix++){
				if(ampx[ix]){
					real x=ix*dx+origx;
					real tmp=ampx[ix]*opdx[ix];
					a0+=tmp;
					a1+=tmp*x;
					a2+=tmp*y;
				}
			}
		}
		//info("a0=%g, a1=%g, a2=%g\n", a0, a1, a2);
		coeff[0]=a0;
		coeff[1]=a1;
		coeff[2]=a2;
		real outp[3];
		dmulvec3(outp, P(imcc,isa), coeff);
		/*
		  2010-07-19: Was =, modified to += to conform to the convention.
		*/
		res[isa]+=outp[1];
		res[isa+nsa]+=outp[2];
	}
}
/**
   Gather information about the starting of each column in loc.

   @param[in,out] loc The loc. Its stat is allocated and filled.
*/
void loc_create_stat_do(loc_t* loc){
	if(!loc) return;
	locstat_t* locstat=mycalloc(1, locstat_t);
	loc->stat=locstat;
	const real* locx=loc->locx;
	const real* locy=loc->locy;
	long nloc=loc->nloc;
	real dx=locstat->dx=loc->dx;
	real dy=locstat->dy=loc->dy;
	long ncolmax=(long)round((locy[nloc-1]-locy[0])/dy)+2;
	locstat->cols=mymalloc(ncolmax, locstatcol_t);
	long colcount=0;
	/*do first column separately. */
	long iloc=0;
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
			||fabs(locx[iloc]-locx[iloc-1]-dx)>dxthres){
			locstat->cols[colcount].pos=iloc;
			locstat->cols[colcount].xstart=locx[iloc];
			locstat->cols[colcount].ystart=locy[iloc];
			if(xmax<locx[iloc-1]){
				xmax=locx[iloc-1];
			}
			if(locstat->xmin>locstat->cols[colcount].xstart){
				locstat->xmin=locstat->cols[colcount].xstart;
			}
			colcount++;
			if(colcount>=ncolmax){/*expand the memory. */
				ncolmax*=2;
				locstat->cols=myrealloc(locstat->cols, ncolmax, locstatcol_t);
			}
		}
	}
	locstat->ncol=colcount;
	locstat->cols[colcount].pos=loc->nloc;
	locstat->cols[colcount].xstart=locx[loc->nloc-1];
	locstat->cols[colcount].ystart=locy[loc->nloc-1];
	if(xmax<locx[loc->nloc-1]){
		xmax=locx[loc->nloc-1];
	}
	locstat->cols=myrealloc(locstat->cols, (locstat->ncol+1), locstatcol_t);
	locstat->ny=(long)round((locy[loc->nloc-1]-locstat->ymin)/dy)+1;
	locstat->nx=(long)round((xmax-locstat->xmin)/dx)+1;
}
/**
   Add a gray pixel circular map in phi using coordinates defined in loc, center
   defined using cx, cy, radius of r, and value of val

   @param[in,out] phi Values on loc to add to.
   @param[in] loc The loc.
   @param[in] cx Center along x (m).
   @param[in] cy Center along y (m).
   @param[in] r Outer radius (m).
   @param[in] rin Inner radius (m). The inner circle is subtracted if >0.
   @param[in] val Value to add inside the annulus.
*/
void loc_circle_add(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val){
	if(!phi||!loc || PN(phi)!=loc->nloc){
		warning("Invalid usage: requires phi, loc and compatible dimension.\n");
		return;
	}
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
	real r2l=(r>dx*1.5&&r>dy*1.5)?(r-dx*1.5)*(r-dy*1.5):0;
	real r2u=(r+dx*1.5)*(r+dy*1.5);
	real* locx=loc->locx;
	real* locy=loc->locy;
	real iix, iiy;
	for(int iloc=0; iloc<loc->nloc; iloc++){
		real x=locx[iloc];
		real y=locy[iloc];
		real r2r=pow(x-cx, 2)+pow(y-cy, 2);
		if(r2r<r2l)
			P(phi,iloc)+=val;
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
			P(phi,iloc)+=(real)tot*res2*val;
		}
	}
	if(rin>0){
		loc_circle_add(phi, loc, cx, cy, rin, 0, -val);
	}
}
/**
   Apply an annular mask in phi.

   @param[in,out] phi Values on loc multiplied by the mask.
   @param[in] loc The loc.
   @param[in] cx Center along x (m).
   @param[in] cy Center along y (m).
   @param[in] r Outer radius (m).
   @param[in] rin Inner radius (m).
   @param[in] val Value of the mask inside the annulus.
*/
void loc_circle_mul(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val){
	dmat *mask=dnew(loc->nloc,1);
	loc_circle_add(mask, loc, cx, cy, r, rin, val);
	dcwm(phi, mask);
	dfree(mask);
}
/**
   Create a gray pixel elliptical map in phi using coordinates defined in loc,
   center defined using cx, cy, radii of rx, ry, and value of val

   @param[in,out] phi Values on loc to add to.
   @param[in] loc The loc.
   @param[in] cx Center along x (m).
   @param[in] cy Center along y (m).
   @param[in] rx Radius along x (m).
   @param[in] ry Radius along y (m).
   @param[in] val Value to add inside the ellipse.
*/
void loc_ellipse_add(dmat* phi, const loc_t* loc, real cx, real cy, real rx, real ry, real val){
	if(!phi||!loc) return;
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
	real r2=1.;//normalized radius of the ellipse
	real dr=MAX(fabs(dx)*rx1, fabs(dy)*ry1);//normalized pixel size
	real r2l=dr<1?(1-dr)*(1-dr):0;
	real r2u=(1+dr)*(1+dr);
	real* locx=loc->locx;
	real* locy=loc->locy;
	real iix, iiy;
	for(int iloc=0; iloc<loc->nloc; iloc++){
		real x=locx[iloc];
		real y=locy[iloc];
		real r2r=pow(x-cx, 2)*rx12+pow(y-cy, 2)*ry12;
		if(r2r<r2l)
			P(phi,iloc)+=val;
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
			P(phi,iloc)+=tot*res2*val;
		}
	}
}
/**
   Remove points in loc that have zero value in amp and modify amp in the same
   time. Keep each row continuous if cont==1. Return in skipout the index of
   skipped points if skipout is not NULL.

   @param[in,out] loc The loc. Reduced in place.
   @param[in,out] amp Amplitude on loc. Reduced in place.
   @param[in] thres Points with amp below thres are removed. EPS is used if <=0.
   @param[in] cont 0: remove any point below thres. 1: only remove points at the ends of each column. Other: within each column, switch between keep and skip only after 4 consecutive points of the same kind.
   @param[out] skipout Newly allocated array of nloc flags marking removed points; the caller frees it. Skipped if NULL.
*/
void loc_reduce(loc_t* loc, dmat* amp, real thres, int cont, int** skipout){
	if(!loc||!amp) return;
	if(thres<=0) thres=EPS;
	int redo_stat=loc->stat?1:0;
	long nloc=loc->nloc;
	int* skip=mycalloc(nloc, int);
	if(cont){/*make sure loc is continuous. */
		loc_create_stat(loc);
		locstat_t* locstat=loc->stat;
		long ncol=locstat->ncol;
		for(long icol=0; icol<ncol; icol++){
			long pstart=locstat->cols[icol].pos;
			long pend=locstat->cols[icol+1].pos;
			if(cont==1){//original method
				for(long pos=pstart; pos<pend; pos++){
					if(P(amp, pos)<thres){
						skip[pos]=1;
					} else{
						break;
					}
				}
				for(long pos=pend-1; pos>pstart-1; pos--){
					if(P(amp, pos)<thres){
						skip[pos]=1;
					} else{
						break;
					}
				}
			}else{//new method
			int nkeep=0, nskip=0;//count continuous keep/skip
			int lastskip=1;
			for(long pos=pstart; pos<pend; pos++){
				if(P(amp,pos)<thres){
					nskip++;
					nkeep=0;
				}else{
					nkeep++;
					nskip=0;
				}
				if(nkeep==4){//continuous keep
					for(int it=1; it<nkeep; it++){
						skip[pos-it]=0;
					}
					lastskip=0;
				}
				if(nskip==4){//continuous skip
					for(int it=1; it<nskip; it++){
						skip[pos-it]=1;
					}
					lastskip=1;
				}
				skip[pos]=lastskip;
			}
			}
		}
	} else{
		for(long iloc=0; iloc<nloc; iloc++){
			if(P(amp,iloc)<thres){
				skip[iloc]=1;
			}
		}
	}
	long count=0;
	for(long iloc=0; iloc<nloc; iloc++){
		loc->locx[count]=loc->locx[iloc];
		loc->locy[count]=loc->locy[iloc];
		P(amp,count)=P(amp,iloc);
		if(!skip[iloc]) count++;
	}
	locresize(loc, count);
	dresize(amp, count, 1);
	if(redo_stat){
		loc_create_stat(loc);
	}
	if(skipout){
		*skipout=skip;
	} else{
		free(skip);
	}
}
/**
   Remove uncoupled points in loc and modify spc in the same time. debugged on 2009-12-20. Not used often. using
   a dspcell, to compute the coupling which is modified accordingly.

   @param[in,out] loc The loc. Reduced in place.
   @param[in,out] spc Sparse matrices. The dimension along loc is reduced in place.
   @param[in] dim 1: loc is along rows of spc. 2: loc is along columns of spc.
   @param[in] cont See loc_reduce().
*/
void loc_reduce_spcell(loc_t* loc, dspcell* spc, int dim, int cont){
	if(!loc||!spc) return;
	long nloc=loc->nloc;
	dmat* sum=NULL;
	for(int isp=0; isp<spc->nx*spc->ny; isp++){
		if(P(spc,isp)){
			dmat* sum0=dspsumabs(P(spc,isp), 3-dim);
			dadd(&sum, 1, sum0, 1);
			dfree(sum0);
		}
	}
	if(!sum||sum->nx*sum->ny!=loc->nloc){
		error("Mismatch\n");
		dfree(sum);
		return;
	}
	int* skip;
	loc_reduce(loc, sum, EPS, cont, &skip);
	dfree(sum);
	long count;
	if(dim==1){/*opd(loc)=sp*opd(locin); */
		count=0;
		spint* map=mycalloc(nloc, spint);
		for(long iloc=0; iloc<nloc; iloc++){
			map[iloc]=(spint)count;
			if(!skip[iloc]) count++;
		}
		for(int isp=0; isp<spc->nx*spc->ny;isp++){
			dsp* sp=P(spc,isp);
			if(!sp) continue;
			sp->nx=count;
			for(long iz=0; iz<sp->pp[sp->ny]; iz++){
				sp->pi[iz]=map[sp->pi[iz]];
			}
		}
		free(map);
	} else if(dim==2){
		for(int isp=0; isp<spc->nx*spc->ny;isp++){
			dsp* sp=P(spc,isp);
			if(!sp) continue;
			count=0;
			for(long iloc=0; iloc<nloc; iloc++){
				if(!skip[iloc]){
					sp->pp[count]=sp->pp[iloc];
					count++;
				}
			}
			sp->pp[count]=sp->pp[nloc];
			sp->ny=count;
			sp->pp=myrealloc(sp->pp, (count+1), spint);
		}
	}
	free(skip);
}

/**
   Remove uncoupled points in loc and modify sp in the same time. debugged on 2009-12-20.  use single sparse
   matrix, which is modified accordingly.

   @param[in,out] loc The loc. Reduced in place.
   @param[in,out] sp Sparse matrix. The dimension along loc is reduced in place.
   @param[in] dim 1: loc is along rows of sp. 2: loc is along columns of sp.
   @param[in] cont See loc_reduce().
*/
void loc_reduce_sp(loc_t* loc, dsp* sp, int dim, int cont){
	if(!loc||!sp) return;
	long nloc=loc->nloc;
	if((dim==1&&nloc!=sp->nx)||(dim==2&&nloc!=sp->ny)||dim<1||dim>2){
		error("Mismatch dimension\n");
		return;
	}
	dmat* sum=dspsumabs(sp, 3-dim);
	int* skip=NULL;
	loc_reduce(loc, sum, EPS, cont, &skip);
	dfree(sum);
	long count=0;
	if(dim==1){/*opd(loc)=sp*opd(locin); */
		spint* map=mycalloc(nloc, spint);
		for(long iloc=0; iloc<nloc; iloc++){
			map[iloc]=(spint)count;
			if(!skip[iloc]) count++;
		}
		sp->nx=count;
		for(long iz=0; iz<sp->pp[sp->ny]; iz++){
			sp->pi[iz]=map[sp->pi[iz]];
		}
		free(map);
	} else if(dim==2){
		for(long iloc=0; iloc<nloc; iloc++){
			if(!skip[iloc]){
				sp->pp[count]=sp->pp[iloc];
				count++;
			}
		}
		sp->pp[count]=sp->pp[nloc];
		sp->ny=count;
		sp->pp=myrealloc(sp->pp, (count+1), spint);
	}
	free(skip);
}

/**
   Add val amount of focus to opd. The unit is in radian like.
   Piston is removed. The focus is centered at ox, oy

   @param[in,out] opd OPD on loc. Modified in place despite const.
   @param[in] loc The loc.
   @param[in] val Coefficient a of focus a*(x^2+y^2).
   @param[in] ox Center along x (m).
   @param[in] oy Center along y (m).
*/
void loc_add_focus_offset(const dmat* opd, const loc_t* loc, const real val, const real ox, const real oy){
	if(!opd||!loc||!val) return;
	if(loc->nloc!=opd->nx*opd->ny){
		error("Invalid dimensions. loc has %ld, opd has %ldx%ld\n", loc->nloc, opd->nx, opd->ny);
		return;
	}
	const real* restrict locx=loc->locx;
	const real* restrict locy=loc->locy;
	real piston=0;
	for(long iloc=0; iloc<loc->nloc; iloc++){
		const real x=locx[iloc]-ox;
		const real y=locy[iloc]-oy;
		real tmp=(x*x+y*y)*val;
		P(opd,iloc)+=tmp;
		piston+=tmp;
	}
	piston=-piston/loc->nloc;
	for(long iloc=0; iloc<loc->nloc; iloc++){
		P(opd,iloc)+=piston;
	}
}
/**
   Add val amount of focus (radian like unit) centered at the origin to opd.
   Piston is removed. See loc_add_focus_offset().

   @param[in,out] opd OPD on loc. Modified in place despite const.
   @param[in] loc The loc.
   @param[in] val Coefficient a of focus a*(x^2+y^2).
*/
void loc_add_focus(const dmat* opd, const loc_t* loc, const real val){
	loc_add_focus_offset(opd, loc, val, 0, 0);
}
/**
   Create a dmat containing locx and locy in two columns

   @param[in] loc The loc.
   @param[in] piston Prepend a column of 1 for piston.
   @return Newly allocated nloc*2 (or nloc*3 with piston) matrix; the caller frees it. NULL if loc is NULL.
*/
dmat* loc2mat(loc_t* loc, int piston){
	if(!loc) return NULL;
	dmat* out=NULL;
	real* ptt=NULL;
	if(piston){
		out=dnew(loc->nloc, 3);
		for(long i=0; i<loc->nloc; i++){
			P(out,i)=1;
		}
		ptt=PCOL(out, 1);
	} else{
		out=dnew(loc->nloc, 2);
		ptt=P(out);
	}
	memcpy(ptt, loc->locx, sizeof(real)*loc->nloc);
	memcpy(ptt+loc->nloc, loc->locy, sizeof(real)*loc->nloc);
	return out;
}

/**
   Convert pts_t to loc_t. Each entry in pts recors the lower
   left point in each subaperture.

   @param[in] pts The subaperture points.
   @return Newly allocated loc of all points in all subapertures; the caller frees it with locfree(). NULL if pts is NULL.
*/
loc_t* pts2loc(pts_t* pts){
	if(!pts) return NULL;
	long nsa=pts->nsa;
	long npsa=pts->nxsa*pts->nysa;
	real dx=pts->dx;
	real dy=pts->dy;
	if(dy==0){
		dy=dx;
	}
	loc_t* loc=locnew(nsa*npsa, dx, dy);
	for(int isa=0; isa<nsa; isa++){
		const real origx=pts->origx[isa];
		const real origy=pts->origy[isa];
		real* locx=loc->locx+isa*npsa;
		real* locy=loc->locy+isa*npsa;
		for(int jj=0; jj<pts->nysa; jj++){
			const real yy=origy+(real)jj*dy;
			for(int ii=0; ii<pts->nxsa; ii++){
				locx[jj*pts->nxsa+ii]=origx+(real)ii*dx;
				locy[jj*pts->nxsa+ii]=yy;
			}
		}
	}
	return loc;
}

/**
   Rotate the coordinate by theta (radian) CCW which is equivalent to rotate the points by theta -CCW.

   @param[in,out] loc The loc rotated in place.
   @param[in] theta Rotation angle (radian).
*/
void locrot(loc_t* loc, const real theta){
	if(!loc||!theta) return;
	const real ctheta=cos(theta);
	const real stheta=sin(theta);

	real* x=loc->locx;
	real* y=loc->locy;
	for(int i=0; i<loc->nloc; i++){
		real tmp=x[i]*ctheta+y[i]*stheta;
		y[i]=-x[i]*stheta+y[i]*ctheta;
		x[i]=tmp;
	}
}
/**
   Determine the angle of rotation from loc1 to loc2. CCW is positive.

   @param[in] loc1 The first loc.
   @param[in] loc2 The second loc with the same number of points as loc1.
   @return Rotation angle (radian); 0 if either is NULL.
*/
real loc_angle(const loc_t* loc1, const loc_t* loc2){
	if(!loc1||!loc2) return 0;
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
	if(count) sina/=count;
	return asin(sina);
}
/**
   Stretch the coordinate by frac along theta (radian) CCW.

   @param[in,out] loc The loc stretched in place.
   @param[in] theta Direction of the stretch (radian).
   @param[in] frac Stretch factor.
*/
void locstretch(loc_t* loc, const real theta, const real frac){
	if(!loc) return;
	locrot(loc, theta);
	for(int i=0; i<loc->nloc; i++){
		loc->locx[i]*=frac;
	}
	locrot(loc, -theta);
}

/**
   duplicate a loc_t object; map and stat are not copied.

   @param[in] loc The loc to copy.
   @return Newly allocated loc; the caller frees it with locfree(). NULL if loc is NULL.
*/
loc_t* locdup(loc_t* loc){
	if(!loc) return NULL;
	loc_t* res=locnew(loc->nloc, loc->dx, loc->dy);
	memcpy(res->locx, loc->locx, sizeof(real)*loc->nloc);
	memcpy(res->locy, loc->locy, sizeof(real)*loc->nloc);
	res->iac=loc->iac;
	res->dratio=loc->dratio;
	res->aoi=loc->aoi;
	res->ht=loc->ht;
	return res;
}
/**
   Compute average of locx, locy

   @param[out] xm Average of locx.
   @param[out] ym Average of locy.
   @param[in] loc The loc.
*/
void locmean(real* xm, real* ym, const loc_t* loc){
	if(!xm||!ym||!loc) return;
	real x=0, y=0;
	for(long i=0; i<loc->nloc; i++){
		x+=loc->locx[i];
		y+=loc->locy[i];
	}
	*xm=x/loc->nloc;
	*ym=y/loc->nloc;
}

/**
   Parse string representation of polynominal to array representation.

   @param[in] _ps Polynomial string in x and y, such as "0.1+x+0.2*x^2*y".
   @return Newly allocated 3*n matrix with coefficient, power of x, and power of y for each term; the caller frees it. NULL if _ps is NULL or cannot be parsed.
*/
dmat* parse_poly(const char* _ps){
	if(!_ps) return NULL;
	if(strchr(_ps, ';')){
		warning("There should not be any ';' :%s", _ps);
	}
	char* ps=(char*)_ps;
	int ncx=5;
	dmat* cx=dnew(3, ncx);
	int icx=0;
	char* endptr;
	while(ps[0]&&ps[0]!=';'){
		P(cx, 0, icx)=strtod(ps, &endptr);//coefficient
		if(ps==endptr){
			if(ps[0]=='-'){
				ps++;
				P(cx, 0, icx)=-1;
			} else if(ps[0]=='+'){
				ps++;
				P(cx, 0, icx)=1;
			}
			if(ps[0]=='x'||ps[0]=='y'){
				if(P(cx, 0, icx)==0) P(cx, 0, icx)=1;
			} else{
				error("Unable to parse (%s). ps=(%s)\n", _ps, ps);
				goto fail;
			}
		} else{
			ps=endptr;
		}
		P(cx, 1, icx)=P(cx, 2, icx)=0;
		while(ps[0]==' ') ps++;
		if(ps[0]=='*') ps++;
		while(ps[0]==' ') ps++;
		//parse each term
		while(ps[0]!=0&&ps[0]!='+'&&ps[0]!='-'&&ps[0]!=';'){
			if(ps[0]=='x'||ps[0]=='y'){
				int iy=1;
				if(ps[0]=='y'){
					iy=2;
				}
				ps++;
				if(ps[0]=='^'){
					ps++;
					P(cx, iy, icx)+=strtol(ps, &endptr, 10);
					if(ps==endptr){
						error("Unable to parse %s\n", _ps);
						goto fail;
					}
					ps=endptr;
				} else{
					P(cx, iy, icx)++;
				}
			} else{
				P(cx, 0, icx)*=strtod(ps, &endptr);
				if(ps==endptr){
					error("Unable to parse %s\n", _ps);
					goto fail;
				}
				ps=endptr;
			}
			while(ps[0]==' '||ps[0]=='*') ps++;//factors may be separated by *
		}
		icx++;
		if(ps[0]==';'||ps[0]==0){
			break;
		} else if(ps[0]=='+'||ps[0]=='-'){
			if(ps[0]=='+') ps++;
			if(icx>=ncx){
				ncx*=2;
				dresize(cx, 3, ncx);
			}
		} else{
			error("Unable to parse (%s)\n", ps);
			goto fail;
		}
	}
	ncx=icx;
	dresize(cx, 3, ncx);
	return cx;
fail:
	dfree(cx);
	return NULL;

}
/**
   Transform coordinate with polynomial coefficients for x and y.

   @param[in] loc Input loc.
   @param[in] cx 3*n matrix of coefficient, power of x, and power of y for the new x.
   @param[in] cy 3*n matrix of coefficient, power of x, and power of y for the new y.
   @return Newly allocated loc; the caller frees it with locfree(). NULL if any input is NULL.
*/
static loc_t* loctransform_do(const loc_t* loc, const dmat* cx, const dmat* cy){
	if(!loc||!cx||!cy) return NULL;
	loc_t* locm=locnew(loc->nloc, loc->dx, loc->dy);
	real* restrict xm=locm->locx;
	real* restrict ym=locm->locy;
	const real* restrict x=loc->locx;
	const real* restrict y=loc->locy;
	int np=0;
	int nonint=0;
	for(int ic=0; ic<cx->ny; ic++){
		int ocx=(int)round(P(cx, 1, ic));
		int ocy=(int)round(P(cx, 2, ic));
		if(fabs(P(cx, 1, ic)-ocx)>EPS||fabs(P(cx, 2, ic)-ocy)>EPS||ocx<0||ocy<0){//not integer or negative
			nonint=1;
		}
		if(ocx>np) np=ocx;
		if(ocy>np) np=ocy;
	}
	for(int ic=0; ic<cy->ny; ic++){
		int ocx=(int)round(P(cy, 1, ic));
		int ocy=(int)round(P(cy, 2, ic));
		if(fabs(P(cy, 1, ic)-ocx)>EPS||fabs(P(cy, 2, ic)-ocy)>EPS||ocx<0||ocy<0){//not integer or negative
			nonint=1;
		}
		if(ocx>np) np=ocx;
		if(ocy>np) np=ocy;
	}
	np++; //max order of x or y +1
	for(long iloc=0; iloc<loc->nloc; iloc++){
		if(nonint){
			for(long ic=0; ic<cx->ny; ic++){
				xm[iloc]+=P(cx, 0, ic)*pow(x[iloc], P(cx, 1, ic))*pow(y[iloc], P(cx, 2, ic));
			}

			for(long ic=0; ic<cy->ny; ic++){
				ym[iloc]+=P(cy, 0, ic)*pow(x[iloc], P(cy, 1, ic))*pow(y[iloc], P(cy, 2, ic));
			}
		} else{/*faster method for integer powers (>10x speed up). */
			real xp[np], yp[np];
			xp[0]=1; yp[0]=1;
			for(long ip=1; ip<np; ip++){
				xp[ip]=xp[ip-1]*x[iloc];
				yp[ip]=yp[ip-1]*y[iloc];
			}

			for(long ic=0; ic<cx->ny; ic++){
				xm[iloc]+=P(cx, 0, ic)*xp[(int)P(cx, 1, ic)]*yp[(int)P(cx, 2, ic)];
			}

			for(long ic=0; ic<cy->ny; ic++){
				ym[iloc]+=P(cy, 0, ic)*xp[(int)P(cy, 1, ic)]*yp[(int)P(cy, 2, ic)];
			}
		}
	}
	return locm;
}
/**
   Transform coordinate with coefficients. See loctransform()

   @param[in] loc Input loc_t.
   @param[in] coeff n*4 matrix with columns power of x, power of y, coefficient for new x, and coefficient for new y.
   @return Newly allocated loc; the caller frees it with locfree(). NULL if any input is NULL.
*/
loc_t* loctransform2(const loc_t* loc, const dmat* coeff){
	if(!loc||!coeff) return NULL;
	dmat* cx=dnew(3, coeff->nx);
	dmat* cy=dnew(3, coeff->nx);
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
   Whether p points to a dmat (cast to char*) rather than a string. The bytes
   are compared one at a time to avoid an unaligned load and reading past the
   end of a short string.

   @param[in] p dmat or string.
   @return 1 if p starts with the dmat magic number.
*/
static int is_dmat_id(const char *p){
	const unsigned char *b=(const unsigned char *)p;
	const uint32_t id=M_REAL;
	for(int i=0; i<4; i++){
		if(b[i]!=((id>>(8*i))&0xFF)) return 0;
	}
	return 1;
}
/**
   Transform coordinate. loc (input) contains x,y; locm (return) contains xm, ym.
   polyn contains the formula
   \f[
   xm{ip}=\sum_{ic}(coeff[0](0,ic)*pow(x,coeff[0](1,ic))*pow(y,coeff[0](2,ic)))
   ym{ip}=\sum_{ic}(coeff[1](0,ic)*pow(x,coeff[1](1,ic))*pow(y,coeff[1](2,ic)))
   \f]
   was using string input. New scheme uses bin files to preserve precision.

   @param[in] loc Input loc.
   @param[in] polycoeff Either a .bin file name of n*4 coefficients, a dmat cast to char*, or a string "polyx;polyy" (see parse_poly()).
   @return Newly allocated loc; the caller frees it with locfree(). NULL if any input is NULL.
*/
loc_t* loctransform(const loc_t* loc, const char* polycoeff){
	if(!loc||!polycoeff) return NULL;
	/*Test whether the transform is pure shift. */
	//Parse from string to 3xn array
	int input_type;
	if(check_suffix(polycoeff, ".bin")){//file input
		input_type=1;
	} else if(is_dmat_id(polycoeff)){//dmat input
		input_type=2;
	} else{
		input_type=0;//string input as formula
	}
	loc_t* locm=0;
	if(input_type>0){
		dmat* coeff=0;
		if(input_type==1){
			coeff=dread("%s", polycoeff);
		} else if(input_type==2){
			coeff=dmat_cast((cell*)polycoeff);
		}
		if(coeff->ny==4){
			locm=loctransform2(loc, coeff);

		} else{
			error("coeff is in wrong format\n");
		}
		if(input_type==1){
			dfree(coeff);
		}
	} else{
		char* polyn=strdup(polycoeff);
		char* px=polyn;
		char* py=strchr(polyn, ';');
		if(py==px||!py){
			error("Wrong format. There should be a ';': {%s}\n", polycoeff);
			free(polyn);
			return NULL;
		}
		py[0]='\0';
		py++;
		//info("polyn=(%s)\npx=(%s)\npy=(%s)\n", _polyn, px, py);
		//Now parse the strings.
		dmat* cx=parse_poly(px);
		dmat* cy=parse_poly(py);
		free(polyn); polyn=0;

		locm=loctransform_do(loc, cx, cy);
		dfree(cx);
		dfree(cy);
	}
	return locm;
}

/**
   Shift a loc coordinate by sx and sy.

   @param[in,out] loc The loc shifted in place.
   @param[in] sx Shift along x (m).
   @param[in] sy Shift along y (m).
*/
void locshift(loc_t* loc, real sx, real sy){
	if(!loc || (!sx && !sy)) return;
	for(long iloc=0; iloc<loc->nloc; iloc++){
		loc->locx[iloc]+=sx;
		loc->locy[iloc]+=sy;
	}
}
/**
   Compute the size of a map that is used to build loc or can fully contain loc
   (embed)

   @param[out] nx Size along x.
   @param[out] ny Size along y.
   @param[in] loc The loc.
*/
void loc_nxny(long* nx, long* ny, const loc_t* loc){
	if(!loc) return;
	real xmax, xmin, ymax, ymin;
	dvecmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
	dvecmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
	*nx=(long)round((xmax-xmin)/loc->dx)+1;
	*ny=(long)round((ymax-ymin)/loc->dy)+1;
}
/**
   Resize a loc_t by shrinking. map and stat are freed.

   @param[in,out] loc The loc resized in place.
   @param[in] nloc New number of points.
*/
void locresize(loc_t* loc, long nloc){
	if(!loc || loc->nloc==nloc) return;
	loc_free_map(loc);
	loc_free_stat(loc);
	dresize(loc->dmat, nloc, 2);
	loc->locy=loc->locx+nloc;
}


#define LOC_EMBED_DEF(X, T, R, OPX)					\
    void X(embed_locstat)(X(mat) **restrict out, real alpha,		\
			  loc_t *restrict loc, R *restrict oin, real beta, int reverse){	\
	if(!loc ||!oin) return ;\
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
	real dy1=1./locstat->dy;					\
	long xoff0=((*out)->nx - locstat->nx +1)/2;			\
	long yoff0=((*out)->ny - locstat->ny +1)/2;			\
									\
	for(long icol=0; icol<locstat->ncol; icol++){			\
	    long xoff=(long)round((locstat->cols[icol].xstart-locstat->xmin)*dx1); \
	    long yoff=(long)round((locstat->cols[icol].ystart-locstat->ymin)*dy1); \
	    long pos1=locstat->cols[icol].pos;				\
	    long pos2=locstat->cols[icol+1].pos;			\
	    T *restrict dest=&P(p,xoff+xoff0,yoff+yoff0);		\
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
/**
   Embeding an OPD defined on loc to another array.
   Do the embeding using locstat to have best speed.
   reverse = 0 : from oin to out: out=out*alpha+in*beta
   reverse = 1 : from out to oin: in=in*beta+out*alpha
*/
LOC_EMBED_DEF(AOS_DMAT, real, real, OPX)
#undef OPX
#define OPX(A) creal(A)
/**
   Embeding an OPD defined on loc to another array.
   Do the embeding using locstat to have best speed.
   reverse = 0 : from oin to out: out=out*alpha+in*beta
   reverse = 1 : from out to oin: in=in*beta+out*alpha
*/
LOC_EMBED_DEF(AOS_CMAT, comp, real, OPX)
#undef OPX

/**
   Determine dx and dy from data.

   @param[in,out] out The loc. dx and dy are set from data if not set or disagree with data.
*/
void loc_dxdy(loc_t* out){
	if(!out) return;
	const real tol=1e-7;

	real dxd=INFINITY, dyd=INFINITY;
	for(long i=0; i<out->nloc-1; i++){
		real dxi=fabs(out->locx[i+1]-out->locx[i]);
		real dyi=fabs(out->locy[i+1]-out->locy[i]);
		if(dxi>tol&&dxi+tol<dxd){
			dxd=dxi;
		}
		if(dyi>tol&&dyi+tol<dyd){
			dyd=dyi;
		}
	}
	if(out->dx==0||isnan(out->dx)){
		out->dx=dxd;//use value derived from data
	} else if(fabs(out->dx-dxd)>tol&&!isinf(dxd)){
		warning("Specified dx=%.15g doesn't agree with data: %.15g, replace.\n", out->dx, dxd);
		out->dx=dxd;
	}

	if(out->dy==0||isnan(out->dy)){
		out->dy=dyd;
	} else if(fabs(out->dy-dyd)>tol&&!isinf(dyd)){
		warning("Specified dy=%.15g doesn't agree with data: %.15g, replace.\n", out->dy, dyd);
		out->dy=dyd;
	}
}

/**
   Generate the keywords string (dx, dy, ht, iac, dratio, aoi) of loc_t for
   saving to file. Used as the make_keywords callback of loc_t.

   @param[in,out] p The loc_t (as cell) whose keywords are regenerated. Ignored if not a loc_t.
*/
void loc_keywords(cell *p){
	if(!p) return;
	if(p->id==M_LOC){
		loc_t *loc=(loc_t*)p;
		if(loc->keywords) free(loc->keywords);
		char str[256];
		snprintf(str, sizeof(str), "dx=%.15g;\ndy=%.15g;ht=%.15g\n;iac=%.15g;\ndratio=%.15g;\naoi=%.15g\n", loc->dx, loc->dy, loc->ht, loc->iac, loc->dratio, loc->aoi);
		loc->keywords=strdup(str);
	}
}

/**
   Parse the input dead actuator location to actuator indices based on aloc.
   2015-03-30: build a mask for dead actuators based on coordinate.

   @param[in] aloc Actuator locations. Its map is created if not yet.
   @param[in,out] dead n*2 or n*3 matrix containing dead actuator coordinates (m) and optionally the stuck value (m). A row vector is reshaped in place.
   @return Newly allocated nloc vector; nonzero marks a dead actuator with value 1 or the stuck value in nm. The caller frees it.
*/
lmat* loc_coord2ind(loc_t* aloc, dmat* dead){
	if(NY(dead)==1 && NX(dead)>1){
		dead->ny=dead->nx;
		dead->nx=1;
	}
	if(NY(dead)!=2&&NY(dead)!=3){
		error("dead must contain 2 or 3 columns of data. %ldx%ld\n", NX(dead), NY(dead));
		return NULL;
	}
	loc_create_map(aloc);
	map_t* map=aloc->map;
	real ox=aloc->map->ox;
	real oy=aloc->map->oy;
	real dx1=1./aloc->dx;
	real dy1=1./aloc->dy;
	dmat* ps=dead/*PDMAT*/;
	lmat* out=lnew(aloc->nloc, 1);
	for(long jact=0; jact<NX(dead); jact++){
		long mapx=(long)round((P(ps, jact, 0)-ox)*dx1);
		long mapy=(long)round((P(ps, jact, 1)-oy)*dy1);
		long iact=loc_map_get(map, mapx, mapy)-1;
		if(iact>=0){
			P(out, iact)=(NY(dead)==3?P(ps, jact, 2)*1e9:1);//integer in nm.
		}
	}
	return out;
}
