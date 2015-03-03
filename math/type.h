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

#ifndef AOS_MATH_TYPE_H
#define AOS_MATH_TYPE_H
#include "../sys/sys.h"
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

#define ARR(T)						\
    long id;       /**< to identify the array type*/	\
    T *restrict p; /**<The data pointer*/		\
    long nx;       /**< number of rows */		\
    long ny;       /**< number of columns */		\
    char *header;  /**<The header*/			\
    struct mmap_t *mmap;/**< not NULL if mmaped.*/	\
    int *nref; /**< reference count */			\
    struct fft_t *fft					

#define MAT(T) struct{ \
	ARR(T);	       \
    }

#define CELL(T) struct{		      \
	ARR(T);			      \
	T m;/*store continuous data*/ \
    }

#define SPMAT(T) struct{						\
	long id;/**<to identify the array type*/			\
	T *restrict x;       /**< numerical values, size nzmax */	\
	union{long m;long nx;};	          /**< number of rows */	\
	union{long n;long ny;};	          /**< number of columns */	\
	char *header;         /**<header*/				\
	long nzmax ;          /**< maximum number of entries */		\
        spint *restrict p ;   /**< column pointers (size n+1) or col indices (size nzmax) when nz!=-1 */ \
	spint *restrict i ;   /**< row indices, size nzmax */		\
	long nz ;             /**< number of entries in triplet matrix, -1 for compressed-col */ \
	int *nref;           /**< reference counting like dmat */	\
    }

typedef MAT(double) dmat;/*a double matrix object contains 2-d array of double numbers*/
typedef MAT(float) smat;
typedef MAT(dcomplex) cmat;
typedef MAT(fcomplex) zmat;
typedef MAT(long) lmat;

typedef SPMAT(double) dsp;
typedef SPMAT(float) ssp;
typedef SPMAT(dcomplex) csp;
typedef SPMAT(fcomplex) zsp;



/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat
*/
typedef struct map_t{
    /*The OPD, takes the same form of dmat so can be casted. */
    ARR(double);
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x*/
    double dy;      /**<Sampling along y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx;      /**Wind velocity. Useful for atmospheric grid*/
    double vy;      /**Wind velocity. Useful for atmospheric grid*/
    int cubic;
    double iac;
} map_t;

/**
   Map with different x/y sampling. Can be cased to dmat
*/
typedef struct rmap_t{
    ARR(double);
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
    long id;
    double *locx;  /**< x coordinates of each point*/
    double *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling along x*/
    double dy;     /**< Sampling along y*/ 
    locstat_t *stat;/**<points to column statistics*/
    map_t *map;    /**< point to the map used for identifying neihboring points.*/
    int npad;      /*padding when create map*/
    int ref;       /**<Data is referenced*/
}loc_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be
   cast to loc_t
*/
typedef struct pts_t{
    long id;
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    union{
	double dsa;    /**<side length of subaperture*/
	double dsax;   /**<side length of subaperture*/
    };
    double dsay;   /**<side length of subaperture*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    map_t *map;    /**<treat pts_t as loc_t and compute the MAP*/
    int npad;      /*padding when create map*/
    int nx;        /**<number of cols per subaperture*/
    int ny;        /**<number of rows per subaperture*/
    double dx;     /**<sampling of points in each subaperture*/
    double dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
}pts_t;

typedef CELL(cmat*) ccell;
typedef CELL(zmat*) zcell;
typedef CELL(dmat*) dcell;
typedef CELL(smat*) scell;
typedef CELL(lmat*) lcell;

typedef CELL(dsp*) dspcell;
typedef CELL(ssp*) sspcell;
typedef CELL(csp*) cspcell;
typedef CELL(zsp*) zspcell;

typedef CELL(ccell*) cccell;
typedef CELL(zcell*) zccell;
typedef CELL(dcell*) dccell;
typedef CELL(scell*) sccell;
typedef CELL(lcell*) iccell;

typedef CELL(cccell*) ccccell;
typedef CELL(zccell*) zcccell;
typedef CELL(dccell*) dcccell;
typedef CELL(sccell*) scccell;
typedef CELL(iccell*) icccell;

typedef CELL(map_t*) mapcell;
typedef CELL(rmap_t*) rmapcell;
typedef CELL(loc_t*) loccell;

typedef CELL(mapcell*) mapccell;
typedef CELL(rmapcell*) rmapccell;
typedef CELL(loccell*) locccell;

typedef struct cell{
    ARR(struct cell*);
    struct cell *m;
}cell;

/*A method to simulate operator overloading for indexing arrys*/
#define IND0(A) error("Invalid use. Use IND(A,i) or IND(A,ix,iy)\n");
#define IND1(A,i) A->p[i]
#define IND2(A,ix,iy) A->p[ix+A->nx*iy]
#define IND_GET(_0,_1,_2,_3,NAME,...) NAME
#define IND(...) IND_GET(_0,__VA_ARGS__,IND2,IND1,IND0,IND0)(__VA_ARGS__)

#endif
