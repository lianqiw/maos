/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

typedef struct{MATARR(double);} dmat;/*a double matrix object contains 2-d array of double numbers*/
typedef struct{MATARR(float);} smat;
typedef struct{MATARR(dcomplex);} cmat;
typedef struct{MATARR(fcomplex);} zmat;
typedef struct{MATARR(long);} lmat;

typedef SPMATARR(double) dsp;
typedef SPMATARR(float) ssp;
typedef SPMATARR(dcomplex) csp;
typedef SPMATARR(fcomplex) zsp;



/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat
*/
typedef struct map_t{
    /*The OPD, takes the same form of dmat so can be casted. */
    MATARR(double);
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x*/
    double dy;      /**<Sampling along y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx;      /**Wind velocity. Useful for atmospheric grid*/
    double vy;      /**Wind velocity. Useful for atmospheric grid*/
    double iac;     /**<Inter-actuator coupling. >0: use cubic influence function*/
} map_t;

/**
   Map with different x/y sampling. Can be cased to dmat
*/
typedef struct rmap_t{
    MATARR(double);
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x (first dimension)*/
    double dy;      /**<Sampling along y (second dimension)*/
    double txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel;    /**<Effective focal length of the telescope*/
    double fexit;   /**<The distance between the exit pupil and the focus*/
    double fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
}rmap_t;

/**
   Store starting x,y for each col
*/
typedef struct locstatcol_t{
    double xstart; /**<starting x of this column*/
    double ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
}locstatcol_t;

/**
   Stores array of locstatcol_t

*/
typedef struct locstat_t{
    locstatcol_t *cols; /**<Information about each column*/
    double dx;          /**<Sampling of the grid along x*/
    double dy;          /**<Sampling of the grid along y*/
    double xmin;        /**<Minimum x*/
    double ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nx,ny;       /**<Size for embedding*/
}locstat_t;

/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
typedef struct loc_t{
    uint32_t id;
    double *locx;  /**< x coordinates of each point*/
    double *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling along x*/
    double dy;     /**< Sampling along y*/
    double ht;     /**< Conjugation height of the loc grid.*/
    double iac;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
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
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    union{
	double dsa;    /**<side length of subaperture*/
	double dsax;   /**<side length of subaperture*/
    };
    double dsay;   /**<side length of subaperture*/
    double dummy1; /**<Place holder*/
    double dummy2;  /**<Place holder*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    map_t *map;    /**<treat pts_t as loc_t and compute the MAP*/
    int npad;      /*padding when create map*/
    int *nref;     /**<Reference counting*/
    int nx;        /**<number of cols per subaperture*/
    int ny;        /**<number of rows per subaperture*/
    double dx;     /**<sampling of points in each subaperture*/
    double dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
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

/*A method to simulate operator overloading for indexing arrys*/
#if DEBUG
static inline void assert_1d(long i, long nx, long ny){
    if(i<0 || i>=nx*ny){
	error("%ld is out of range for (%ld,%ld) array\n", i, nx, ny);
    }
}
static inline void assert_2d(long ix, long iy, long nx, long ny){
    if(ix<0 || ix>=nx || iy<0 || iy>=ny){
	error("(%ld,%ld) is out of range for (%ld,%ld) array\n", ix, iy, nx, ny);
    }
}
#define P1(A,i) ((A)->p[assert_1d((i), (A)->nx, (A)->ny),(i)])
#define P2(A,ix,iy) ((A)->p[assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)])
//#define PP1(A,i) ((A)->p+(assert_1d(i, (A)->nx, (A)->ny),(i)))
//#define PP2(A,ix,iy) ((A)->p+(assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)))
#else
#define P1(A,i) ((A)->p[(i)])
#define P2(A,ix,iy) ((A)->p[(ix)+(A)->nx*(iy)])
#endif
#define P0(A) ((A)->p)
#define PP0(A) ((A)->p)
#define PP1(A,i) ((A)->p+(i))
#define PP2(A,ix,iy) ((A)->p+(ix)+(A)->nx*(iy))
//#endif
//#define P0(A) _Pragma("#error Invalid use. Use P(A,i) or P(A,ix,iy)\n")
//#define PP0(A) _Pragma("#error Invalid use. Use PP(A,i) or PP(A,ix,iy)\n")
#define P_GET(_0,_1,_2,_3,NAME,...) NAME
#define P(...) P_GET(_0,__VA_ARGS__,P2,P1,P0)(__VA_ARGS__)
#define PP(...) P_GET(_0,__VA_ARGS__,PP2,PP1,PP0)(__VA_ARGS__)
#define PCOL(A,iy) ((A)->p+(iy)*(A)->nx)

//Define indexing using wrapping. See wrap()
//#define P1R(A,i) _Pragma("#error Invalid use. Use PR(A,i,j)")
#define PR(A,ix,iy) P2(A, wrap(ix, A->nx), wrap(iy, A->ny))
//#define PR(...) P_GET(_0,__VA_ARGS__,P2R,P1R,P1R,P1R)(__VA_ARGS__)
//#define PP1R(A,i) _Pragma("#error Invalid use. Use PPR(A,i,j)")
#define PPR(A,ix,iy) PP2(A, wrap(ix, A->nx), wrap(iy, A->ny))
//#define PPR(...) P_GET(_0,__VA_ARGS__,PP2R,PP1R,PP1R,PP1R)(__VA_ARGS__)
#define PCOLR(A,iy) ((A)->p+wrap(iy, A->ny)*(A)->nx)
#endif
