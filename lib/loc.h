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
#ifndef AOS_LIB_LOC_H
#define AOS_LIB_LOC_H
#include "../math/mathdef.h"
#include "map.h"
/**
   \file loc.h
   This file defines functions relates to pts_t, loc_t, map_t, etc.
 */
/**
   Store starting x,y for each col
*/
typedef struct locstatcol_t{
    real xstart; /**<starting x of this column*/
    real ystart; /**<starting y of this column*/
    long pos;    /**<starting index of this column*/
}locstatcol_t;

/**
   Stores array of locstatcol_t

*/
typedef struct locstat_t{
    locstatcol_t* cols; /**<Information about each column*/
    real dx;          /**<Sampling of the grid along x*/
    real dy;          /**<Sampling of the grid along y*/
    real xmin;        /**<Minimum x*/
    real ymin;        /**<Minimum y*/
    long ncol;        /**<Number of consecutive columns found*/
    long nx, ny;       /**<Size for embedding*/
}locstat_t;

/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
struct loc_t{
    union{
        dmat dmat[1];
        struct{//must be the same as dmat
            M_ID id;
            real* locx;  /**< x coordinates of each point [allocates memory for both locx and locy]*/
            long nloc;   /**< number of points*/
            long two;    /**<Constant 2. to be binary compatble with rmat*/
            char *keywords; /**<The keywords as a string*/
		};
	};
	
	real* locy;  /**< y coordinates of each point*/
	locstat_t* stat;/**<points to column statistics*/
	map_t* map;    /**< point to the map used for identifying neihboring points.*/
	real dx;     /**< Sampling along x*/
	real dy;     /**< Sampling along y*/
	real ht;     /**< Conjugation height of the loc grid.*/
	real iac;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
	real dratio;  /**<M1 to DM demagnification ratio (>1 for shrinking beam)*/
	real aoi;	  /**<Angle of incidence in radian for ray tracing input */
	int npad;      /*padding when create map*/
};
/**
   low left point of each subaperture.

   don't change the leading elements. so that pts_t can be used as loc_t.
*/
struct pts_t{
	union{
		loc_t loc[1];
		struct{
			union{
				dmat dmat[1];
				struct{
					M_ID id;
					real* origx;    /**<The x origin of each subaperture. Contains memory for origy*/
					long nsa;       /**<number of subapertures.*/
					long two;       /**<Constant 2. to be binary compatble with rmat*/
					char *keywords; /**<The keywords as a string*/
				};
			};
	        real* origy; 		/**<The y origin of each subaperture*/
            locstat_t* stat;	/**<padding so that we can be casted to loc_t*/
            map_t* map;    		/**<treat pts_t as loc_t and compute the MAP*/
            union{
                real dsa;    /**<side length of subaperture*/
                real dsax;   /**<side length of subaperture*/
            };
            real dsay;   	 /**<side length of subaperture. */
            real dummy_ht; 	 /**<Place holder*/
            real dummy_iac;  /**<Place holder*/
			real dratio;  	 /**<M1 to DM demagnification ratio (>1 for shrinking beam)*/
			real aoi;	     /**<Angle of incidence in radian for ray tracing input */
            int npad;      	 /**<padding when create map*/
		};
	};
	//before this line same memory layout as loc_t
	//the following are specific to pts_t
	long nxsa;       /**<number of cols per subaperture*/
	long nysa;       /**<number of rows per subaperture*/
	real dx;     /**<sampling of points in each subaperture*/
	real dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
};
//loc_t can be used as dmat with nx=nloc and ny=2
TYPE_ASSERT_SAME(loc_t, id, dmat, id);
TYPE_ASSERT_SAME(loc_t, locx, dmat, p);
TYPE_ASSERT_SAME(loc_t, nloc, dmat, nx);
TYPE_ASSERT_SAME(loc_t, two, dmat, ny);
TYPE_ASSERT_SAME(loc_t, keywords, dmat, keywords);
TYPE_ASSERT(offsetof(loc_t, locy)>=sizeof(dmat));
//pts_t can be used as loc_t
TYPE_ASSERT_SAME(pts_t, origx, loc_t, locx);
TYPE_ASSERT_SAME(pts_t, nsa, loc_t, nloc);
TYPE_ASSERT_SAME(pts_t, origy, loc_t, locy);
TYPE_ASSERT_SAME(pts_t, stat, loc_t, stat);
TYPE_ASSERT_SAME(pts_t, map, loc_t, map);
TYPE_ASSERT_SAME(pts_t, dsax, loc_t, dx);
TYPE_ASSERT_SAME(pts_t, dsay, loc_t, dy);
TYPE_ASSERT_SAME(pts_t, dummy_ht, loc_t, ht);
TYPE_ASSERT_SAME(pts_t, dummy_iac, loc_t, iac);
TYPE_ASSERT_SAME(pts_t, dratio, loc_t, dratio);
TYPE_ASSERT_SAME(pts_t, aoi, loc_t, aoi);
TYPE_ASSERT_SAME(pts_t, npad, loc_t, npad);
TYPE_ASSERT(offsetof(pts_t, nxsa)>=sizeof(loc_t));

#define LOC(A)  ((A)?(A)->loc:(loc_t*)0)
loc_t* loc_convert(dmat* A);
cell* loc_convert_cell(cell* A);
#define locwrite(out, A...) writebin(out, A)
#define locread(A...)    loc_convert(dread(A))
#define loccellread(A...) (loccell*)cellconvert(readbin_id(M_LOC, 1, A), loc_convert_cell)
#define loccellnew (loccell*)cellnew
#define locccellnew (locccell*)cellnew
#define loccellfree(A) cellfree(A) //if(A){cellfree_do((cell*)A); A=NULL;};
#define locfree(A) cellfree(A) //if(A){cellfree_do((cell*)A);A=NULL;}
#define ptsfree(A) cellfree(A) //if(A){cellfree_do((cell*)A);A=NULL;}
//#define loccellwrite(A,B...) cellwrite((cell*)A,B)
lmat* loc_create_embed(long* nembed, const loc_t* loc, real oversize, int fftpad);
void loc_create_map_npad(const loc_t* loc, int npad, int nx, int ny);
void loc_create_map(const loc_t* loc);
/**
   Obtain an entry in the map, with boundary checking enforced. For the map
   created by loc_create_map(), the value is the loc index plus 1, and 0 means
   no point or out of range.

   @param[in] map Map to read.
   @param[in] ix Index along x.
   @param[in] iy Index along y.
   @return Value of the map at (ix, iy), converted to long; 0 if out of range.
*/
static inline long loc_map_get(map_t* map, long ix, long iy){
	if(ix>=0&&ix<map->nx&&iy>=0&&iy<map->ny){
		return (long)map->p[ix+iy*map->nx];
	} else{
		return 0;
	}
}
uint32_t lochash(const loc_t* loc, uint32_t key);
void loc_embed(anydmat _dest, const loc_t *loc, const_anydmat in, int expanded, int add);
dcell* loc_embed2(const loc_t* loc, const dmat* arr);
void loc_embed_cell(dcell** dest, const loc_t* loc, const dcell* in);
void loc_extract(dmat* dest, const loc_t* loc, map_t* in);
loc_t* loc_from_map(map_t* amp, real thres);

void loc_free_map(loc_t* loc);
void locfree_do(cell* loc);

real loc_diam(const loc_t* loc);
int loccenter(const loc_t* loc);
loc_t* locnew(long nloc, real dx, real dy);
loc_t* locref(loc_t* in);
pts_t* ptsnew(long nsa, real dsax, real dsay, long nxsa, long nysa, real dx, real dy);
void loc_calc_ptt(real* out, real* coeffout,
	const loc_t* loc, real ipcc,
	const dmat* imcc, const real* amp, const real* opd);
void loc_calc_mod(real* out, real* coeffout,
	const dmat* mod, const real* amp, real* opd);
dmat* loc_mcc_ptt(const loc_t* loc, const real* amp);
dcell* pts_mcc_ptt(const pts_t* pts, const real* amp);
void loc_sub_ptt(dmat* opd, const real* ptt, const loc_t* loc);
void loc_add_ptt(dmat* opd, const real* ptt, const loc_t* loc);
void loc_remove_ptt(dmat *opd, const loc_t *loc, const real *amp, const dmat *imcc, int both);
real loc_remove_focus_grad(dmat *grad, const loc_t *saloc, real factor);
void pts_ztilt(dmat** out, const pts_t* pts, const dcell* imcc,
	const real* amp, const real* opd);
loc_t* mk1dloc_vec(real* x, long nx);
loc_t* mk1dloc(real x0, real dx, long nx);
loc_t* mksqloc_auto(long nx, long ny, real dx, real dy);
loc_t* mksqloc_map(const map_t* map);
loc_t* mksqloc(long nx, long ny, real dx, real dy, real ox, real oy);
loc_t* mkannloc(real D, real Din, real dx, real thres);

void loc_create_stat_do(loc_t* loc);
#define loc_create_stat(loc) if(!loc->stat) loc_create_stat_do(loc);
void loc_free_stat(loc_t* loc);
void loc_circle_add(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val);
void loc_circle_mul(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val);
void loc_ellipse_add(dmat* phi, const loc_t* loc, real cx, real cy,real rx, real ry, real val);
void loc_reduce(loc_t* loc, dmat* amp, real thres, int cont, int** skipout);
void loc_reduce_spcell(loc_t* loc, dspcell* sp, int dim, int cont);
void loc_reduce_sp(loc_t* loc, dsp* sp, int dim, int cont);

void loc_add_focus(const dmat* opd, const loc_t* loc, real val);
void loc_add_focus_offset(const dmat* opd, const loc_t* loc, real val, real ox, real oy);
dmat* loc2mat(loc_t* loc, int piston);
loc_t* pts2loc(pts_t* pts);
void locrot(loc_t* loc, const real theta);
real loc_angle(const loc_t* loc1, const loc_t* loc2);
void locstretch(loc_t* loc, const real theta, const real frac);
loc_t* locdup(loc_t* loc);
void locmean(real* xm, real* ym, const loc_t* loc);
dmat *parse_poly(const char *_ps);
loc_t* loctransform(const loc_t* loc, const char* ps);
loc_t* loctransform2(const loc_t* loc, const dmat* coeff);
void locshift(loc_t* loc, real sx, real sy);
void loc_nxny(long* nx, long* ny, const loc_t* loc);

void locresize(loc_t* loc, long nloc);
#define ptsresize(pts, nsa) if(pts) locresize(pts->loc, nsa)
void dembed_locstat(dmat** out, real alpha, loc_t* loc, real* oin, real beta, int reverse);
void cembed_locstat(cmat** out, real alpha, loc_t* loc, real* oin, real beta, int reverse);
void loc_dxdy(loc_t* loc);
loc_t* locreaddata(file_t* fp, header_t* header);
void loc_keywords(cell *loc);
lmat* loc_coord2ind(loc_t *aloc, dmat *dead);
#endif
