/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_TYPE_H
#define AOS_TYPE_H
#include "common.h"
#include "misc.h"
#include "bin.h"

/**
   \file lib/type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.

   Don't use ulong for dimensions because subtracting a bigger ulong from a
   smaller ulong overflows.  */


/*
  Separate definition of struct with typedef. Put definition of struct in a
  private file and typedef in a public interface to hide the struct from the
  user */
typedef enum{
    M_DMAT=1,
    M_CMAT,
    M_DSP,
    M_CSP,
    M_CELL,
}M_TYPE;
typedef void (*vtbl_write)(void *A, const char *format, ...);
typedef void (*vtbl_writedata)(file_t *fp, void *A);
typedef void *(*vtbl_read)(const char *format, ...);
typedef void *(*vtbl_readdata)(file_t *fp);
/**
   A virtual table that contails pointers to functions.
*/
typedef struct vtbl{
    M_TYPE type;/**<Denoting the type of this data*/
    vtbl_write write;
    vtbl_writedata writedata;
    vtbl_read read;
    vtbl_readdata readdata;
}vtbl;
extern vtbl dmat_vtbl;
extern vtbl cmat_vtbl;
extern vtbl cell_vtbl;
extern vtbl dsp_vtbl;
extern vtbl csp_vtbl;
extern vtbl map_vtbl;
extern vtbl rectmap_vtbl;
extern vtbl loc_vtbl;
extern vtbl pts_vtbl;
#define MAT(T) \
    vtbl *vtbl;    /**<The virtual function table*/\
    T *restrict p; /**<The data pointer*/\
    long nx;       /**< number of rows */\
    long ny;       /**< number of columns */\
    char *header;  /**<The header*/\
    struct mmap_t *mmap;/**< not NULL if mmaped. unmap the mmaped memory*/
/**
   a double matrix object contains 2-d array of double numbers
*/
typedef struct dmat{
    MAT(double)
    long *nref; /**< reference count */
}dmat;

/**
   a double complex matrix object contains 2-d arrays of
   double complex numbers. */
typedef struct cmat{
    MAT(dcomplex)
    long *nref;  /**< reference count */
    struct fft_t *fft;
}cmat;

/**
   an 2-d block matrix of cmat.
 */
typedef struct ccell{
    MAT(cmat*)
}ccell;
/**
   an 2-d block matrix of dmat.
 */
typedef struct dcell{
    MAT(dmat*)
}dcell;

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;
/**
   a sparse array of double numbers stored in
   compressed column format, i.e. MATLAB format */
typedef struct dsp{
    vtbl *vtbl;           /**<The virtual function table*/
    double *restrict x ;  /**< numerical values, size nzmax */
    long m ;	          /**< number of rows */
    long n ;	          /**< number of columns */
    char *header;         /**<header*/
    long nzmax ;          /**< maximum number of entries */
    spint *restrict p ;   /**< column pointers (size n+1) or col indlces (size nzmax) when nz!=-1 */
    spint *restrict i ;   /**< row indices, size nzmax */
    long nz ;             /**< number of entries in triplet matrix, -1 for compressed-col */
    long *nref;           /**< reference counting like dmat */
}dsp;

/**
   a sparse array of double complex numbers stored in
   compressed column format */
typedef struct csp{
    vtbl *vtbl; /**<The virtual function table*/
    dcomplex *x;/**< numerical values, size nzmax */
    long m ;	/**< number of rows */
    long n ;	/**< number of columns */
    char *header;/**<header*/
    long nzmax ;/**< maximum number of entries */
    spint *p ;   /**< column pointers (size n+1) or col indlces (size nzmax) */
    spint *i ;   /**< row indices, size nzmax */
    long nz ;   /**< # of entries in triplet matrix, -1 for compressed-col */
    long *nref; /**< reference counting like cmat*/
}csp;

/**
   an 2-d array of sparse.
 */
typedef struct spcell{
    MAT(dsp*)
}spcell;

/**
   an 2-d array of csp.
 */
typedef struct cspcell{
    MAT(csp*)
}cspcell;


/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat
*/
typedef struct map_t{
    //The OPD, takes the same form of dmat so can be casted.
    MAT(double)
    long *nref; /**< reference count */
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x and y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx;      /**Wind velocity. Useful for atmospheric grid*/
    double vy;      /**Wind velocity. Useful for atmospheric grid*/
#if USE_POSIX_SHM == 1
    //long shm;       /**Records the length of the memory mmaped of positive. -1 means it is part of shared shm.*/
#endif
} map_t;

/**
   Map with different x/y sampling. Can be cased to dmat
*/
typedef struct rectmap_t{
    MAT(double)
    long *nref;     /**< reference count */
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x (first dimension)*/
    double dy;      /**<Sampling along y (second dimension)*/
    double txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel;    /**<Effective focal length of the telescope*/
    double fexit;   /**<The distance between the exit pupil and the focus*/
    double fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
}rectmap_t;

/**
   map of locs
*/
typedef struct locmap_t{
    long *restrict p;       /**<The map, of size nx*ny*/
    int nx;        /**<Number of points along x*/
    int ny;        /**<Number of points along y*/
    double ox;     /**<Origin of the map along x*/
    double oy;     /**<Origin of the map along y*/
    int npad;      /**<Padding along the boundary. just for checking. no need in computation.*/
}locmap_t;

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
    double dx;          /**<Sampling of the grid*/
    double xmin;        /**<Minimum x*/
    double ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nrow;        /**<Maximum number of rows*/
}locstat_t;
/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
typedef struct loc_t{
    double *restrict locx;  /**< x coordinates of each point*/
    double *restrict locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling*/
    locmap_t *restrict map; /**< point to the map used for identifying neihboring points.*/
    locstat_t *stat;/**<points to column statistics*/
}loc_t;

/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be
   cast to loc_t
*/
typedef struct pts_t{
    double *restrict origx; /**<The x origin of each subaperture*/
    double *restrict origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    double dsa;    /**<side length of subaperture*/
    locmap_t *restrict map; /**<treat pts_t as loc_t and compute the MAP*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    double dx;     /**<sampling of points in each subaperture*/
    int nx;        /**<number points in each col or row per subaperture*/
}pts_t;

/**
   Cell of anything
*/
typedef struct cell{
    MAT(void*)
}cell;

#define AOS_CMAT(A) c##A
#define AOS_CSP(A) c##A
#define AOS_DMAT(A) d##A
#define AOS_DSP(A) A

#endif
