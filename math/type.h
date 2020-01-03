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

#ifndef AOS_MATARRH_TYPE_H
#define AOS_MATARRH_TYPE_H
#include "numtype.h"
/**
   \file type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.

   Don't use ulong for dimensions because subtracting a bigger ulong from a
   smaller ulong overflows.  */

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;

/*
  We use pointers for reference counter because different array
  may use the same pointer, but with different nx or ny
  partition. */

#define ARR(T)								\
    uint32_t id;   /**< to identify the array type. Must be the first element*/	\
    T *restrict p; /**<The data pointer*/				\
    long nx;       /**< number of rows */				\
    long ny;       /**< number of columns */				\
    char *header;  /**<The header*/					\
    struct fft_t *fft					

#define MATARR(T)				\
    ARR(T);					\
    struct mem_t *mem /**< Memory management*/	\
    

#define CELLARR(T)struct{			\
	ARR(T);					\
	T m; /*store continuous data*/		\
    }

#define SPMATARR(T) struct{						\
	uint32_t id;         /**<to identify the array type*/		\
	T *restrict x;       /**< numerical values, size nzmax */	\
	long nx;             /**< number of rows */			\
	long ny;             /**< number of columns */			\
	char *header;        /**<header*/				\
	long nzmax ;         /**< maximum number of entries */		\
        spint *restrict p ;  /**< column pointers (size n+1) or col indices (size nzmax) when nz!=-1 */ \
	spint *restrict i ;  /**< row indices, size nzmax */		\
	int *nref;           /**< reference counting like dmat */	\
    }

typedef struct{MATARR(float);} smat;
typedef struct{MATARR(fcomplex);} zmat;
#ifdef USE_DOUBLE
typedef struct{MATARR(real);} dmat;/*a real matrix object contains 2-d array of real numbers*/
typedef struct{MATARR(comp);} cmat;
#else
typedef smat dmat;
typedef zmat cmat;
#endif
typedef struct{MATARR(long);} lmat;
typedef struct{MATARR(int);} imat;

typedef SPMATARR(float) ssp;
typedef SPMATARR(fcomplex) zsp;
#ifdef USE_DOUBLE
typedef SPMATARR(real) dsp;
typedef SPMATARR(comp) csp;
#else
typedef ssp dsp;
typedef zsp csp;
#endif

/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat
*/
typedef struct map_t{
    /*The OPD, takes the same form of dmat so can be casted. */
    MATARR(real);
    real ox;      /**<Origin in x*/
    real oy;      /**<Origin in y*/
    real dx;      /**<Sampling along x*/
    real dy;      /**<Sampling along y*/
    real h;       /**<Heigh conjugation of this surface*/
    real vx;      /**Wind velocity. Useful for atmospheric grid*/
    real vy;      /**Wind velocity. Useful for atmospheric grid*/
    real iac;     /**<Inter-actuator coupling. >0: use cubic influence function*/
} map_t;

/**
   Map with different x/y sampling. Can be cased to dmat
*/
typedef struct rmap_t{
    MATARR(real);
    real ox;      /**<Origin in x*/
    real oy;      /**<Origin in y*/
    real dx;      /**<Sampling along x (first dimension)*/
    real dy;      /**<Sampling along y (second dimension)*/
    real txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    real tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    real ftel;    /**<Effective focal length of the telescope*/
    real fexit;   /**<The distance between the exit pupil and the focus*/
    real fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
}rmap_t;

/**
   Store starting x,y for each col
*/
typedef struct locstatcol_t{
    real xstart; /**<starting x of this column*/
    real ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
}locstatcol_t;

/**
   Stores array of locstatcol_t

*/
typedef struct locstat_t{
    locstatcol_t *cols; /**<Information about each column*/
    real dx;          /**<Sampling of the grid along x*/
    real dy;          /**<Sampling of the grid along y*/
    real xmin;        /**<Minimum x*/
    real ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nx,ny;       /**<Size for embedding*/
}locstat_t;

/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
typedef struct loc_t{
    uint32_t id;
    real *locx;  /**< x coordinates of each point*/
    real *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    real dx;     /**< Sampling along x*/
    real dy;     /**< Sampling along y*/
    real ht;     /**< Conjugation height of the loc grid.*/
    real iac;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
    locstat_t *stat;/**<points to column statistics*/
    map_t *map;    /**< point to the map used for identifying neihboring points.*/
    int npad;      /*padding when create map*/
    int *nref;       /**<Reference counting*/
}loc_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be used as loc_t.
*/
typedef struct pts_t{
    uint32_t id;
    real *origx; /**<The x origin of each subaperture*/
    real *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    union{
	real dsa;    /**<side length of subaperture*/
	real dsax;   /**<side length of subaperture*/
    };
    real dsay;   /**<side length of subaperture*/
    real dummy1; /**<Place holder*/
    real dummy2;  /**<Place holder*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    map_t *map;    /**<treat pts_t as loc_t and compute the MAP*/
    int npad;      /*padding when create map*/
    int *nref;     /**<Reference counting*/
    int nx;        /**<number of cols per subaperture*/
    int ny;        /**<number of rows per subaperture*/
    real dx;     /**<sampling of points in each subaperture*/
    real dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
}pts_t;

typedef CELLARR(cmat*) ccell;
typedef CELLARR(zmat*) zcell;
typedef CELLARR(dmat*) dcell;
typedef CELLARR(smat*) scell;
typedef CELLARR(lmat*) lcell;

typedef CELLARR(dsp*) dspcell;
typedef CELLARR(ssp*) sspcell;
typedef CELLARR(csp*) cspcell;
typedef CELLARR(zsp*) zspcell;

typedef CELLARR(ccell*) cccell;
typedef CELLARR(zcell*) zccell;
typedef CELLARR(dcell*) dccell;
typedef CELLARR(scell*) sccell;
typedef CELLARR(lcell*) iccell;

typedef CELLARR(cccell*) ccccell;
typedef CELLARR(zccell*) zcccell;
typedef CELLARR(dccell*) dcccell;
typedef CELLARR(sccell*) scccell;
typedef CELLARR(iccell*) icccell;

typedef CELLARR(map_t*) mapcell;
typedef CELLARR(rmap_t*) rmapcell;
typedef CELLARR(loc_t*) loccell;

typedef CELLARR(mapcell*) mapccell;
typedef CELLARR(rmapcell*) rmapccell;
typedef CELLARR(loccell*) locccell;

typedef struct cell{
    ARR(struct cell*);
    struct cell *m;
}cell;
#undef ARR
#undef CELLARR
#undef MATARR
static inline int iscell(const void *id){
    const uint32_t magic=*((const uint32_t*)id);
    return magic==MCC_ANY;
    //return (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420);
}
/*A method to simulate operator overloading for indexing arrys*/
#if DEBUG
static inline long index_1d(long i, long nx, long ny){
    if(i<0 || i>=nx*ny){
	error("%ld is out of range for (%ld,%ld) array\n", i, nx, ny);
    }
    return i;
}
static inline long index_2d(long ix, long iy, long nx, long ny){
    if(ix<0 || ix>=nx || iy<0 || iy>=ny){
	error("(%ld,%ld) is out of range for (%ld,%ld) array\n", ix, iy, nx, ny);
    }
    return ix+iy*nx;
}
#else
#define index_1d(i,    nx,ny) (i)
#define index_2d(ix,iy,nx,ny) ((ix)+(iy)*(nx))
#endif
    
#define P0(A)  ((A)->p)

#define P1(A,i)     (A)->p[index_1d((i),        (A)->nx, (A)->ny)]
#define P2(A,ix,iy) (A)->p[index_2d((ix), (iy), (A)->nx, (A)->ny)]

#define PP1(A,i)     ((A)->p+index_1d((i),        (A)->nx, (A)->ny))
#define PP2(A,ix,iy) ((A)->p+index_2d((ix), (iy), (A)->nx, (A)->ny))

#define P_GET(_0,_1,_2,_3,NAME,...) NAME
#define P(...) P_GET(_0,__VA_ARGS__,P2,P1,P0)(__VA_ARGS__)
#define PP(...) P_GET(_0,__VA_ARGS__,PP2,PP1,P0)(__VA_ARGS__)
#define PCOL(A,iy) ((A)->p+(iy)*(A)->nx)

//Define indexing using wrapping. See wrap()
#define PR(A,ix,iy)   P2((A), wrap((ix), (A)->nx), wrap((iy), (A)->ny))
#define PPR(A,ix,iy) PP2((A), wrap((ix), (A)->nx), wrap((iy), (A)->ny))
#define PCOLR(A,iy) ((A)->p+wrap(iy, A->ny)*(A)->nx)

//Return Number of elements
#define PN(A) ((A)->nx*(A)->ny)

#endif
