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


#ifndef AOS_MATARRH_TYPE_H
#define AOS_MATARRH_TYPE_H
#include "numtype.h"
/**
   \file type.h

   Defines the math data types like dmat, cmat, dcell, ccell, dsp, csp data
   types.

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
    M_ID id;   /**< to identify the array type. Must be the first element*/	\
    T *restrict p; /**<The data pointer*/				\
    long nx;       /**< number of rows */				\
    long ny;       /**< number of columns */				\
    char *keywords;/**<The keywords as a string*/					\
    file_t *fp;    /**<Opened file, to be saved to when called or freed*/\
    struct fft_t* fft					\

typedef struct cell{
    ARR(struct cell*);
    struct cell* m;
}cell;

//Make sure the memory layout of CELLDEF is identical to cell
//base[1] is conveniently used to return a pointer as cell without casting which is error prone
#define CELLDEF(T,S) typedef struct S{		\
    union{  \
            cell cell[1];\
            struct{\
	            ARR(struct T*);				\
	            struct T* m; /*continuous data*/	\
            };\
        };\
    }S


#define MATARR(T)				\
    union{\
        cell cell[1];\
        struct{\
            ARR(T);					\
            mem_t *mem; /**< Memory management*/	\
            async_t *async; /**<async io*/\
        };\
    }\

#define MATDEF(T,S) typedef struct S{MATARR(T);} S

#define SPMATDEF(T,S) typedef struct S{					        \
    union{\
        cell cell[1];     /**<Use array to get pointer easily*/  \
        struct{\
            M_ID id;             /**<to identify the array type*/		\
            T *restrict px;      /**<numerical values, size nzmax */	\
            long nx;             /**<number of rows */			        \
            long ny;             /**<number of columns */			    \
            char *keywords;      /**<the keywords as a string*/         \
            char *fn;            /**<The file, to be saved upon free*/  \
            long nzmax;          /**<maximum number of entries */		\
            spint *restrict pp;  /**<col indices (size nzmax)  */       \
            spint *restrict pi;  /**<row indices, size nzmax */		    \
            unsigned int *nref;           /**<reference counting for px, pp, pi*/\
        };\
    };/*todo: padding base to the maximum size. Make sure to update python ctypes definition as well*/\
}S


MATDEF(float, smat);
MATDEF(fcomplex, zmat);

MATDEF(real, dmat);
MATDEF(comp, cmat);

MATDEF(long, lmat);
MATDEF(int, imat);

SPMATDEF(float, ssp);
SPMATDEF(fcomplex, zsp);

SPMATDEF(real, dsp);
SPMATDEF(comp, csp);



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
   TODO: use rmat interface of mem_t?
*/
typedef struct loc_t{
    union{
        cell cell[1];
        dmat dmat[1];
        struct{
            M_ID id;
            real* locx;  /**< x coordinates of each point [allocates memory for both locx and locy]*/
            long nloc;   /**< number of points*/
            long two;    /**<Constant 2. to be binary compatble with rmat*/
            real* locy;  /**< y coordinates of each point*/
            locstat_t* stat;/**<points to column statistics*/
            map_t* map;    /**< point to the map used for identifying neihboring points.*/
            unsigned int *nref;    /**<Reference counting*/
            real dx;     /**< Sampling along x*/
            real dy;     /**< Sampling along y*/
            real ht;     /**< Conjugation height of the loc grid.*/
            real iac;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
            int npad;      /*padding when create map*/
        };
    };
}loc_t;
/**
   low left point of each subaperture.

   don't change the leading elements. so that pts_t can be used as loc_t.
*/
typedef struct pts_t{
    union{
        cell cell[1];
        dmat dmat[1];
        loc_t loc[1];
        struct{
            M_ID id;
            real* origx; /**<The x origin of each subaperture. Contains memory for origy*/
            long nsa;      /**<number of subapertures.*/
            long two;
            real* origy; /**<The y origin of each subaperture*/
            locstat_t* stat;/**<padding so that we can be casted to loc_t*/
            map_t* map;    /**<treat pts_t as loc_t and compute the MAP*/
            unsigned int *nref;     /**<Reference counting*/
            union{
                real dsa;    /**<side length of subaperture*/
                real dsax;   /**<side length of subaperture*/
            };
            real dsay;   /**<side length of subaperture*/
            real dummy1; /**<Place holder*/
            real dummy2;  /**<Place holder*/
            int npad;      /*padding when create map*/
            //before this line same memory layout as loc_t
            //the following are specific to pts_t
            int nxsa;        /**<number of cols per subaperture*/
            int nysa;        /**<number of rows per subaperture*/
            real dx;     /**<sampling of points in each subaperture*/
            real dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
        };
    };
}pts_t;

CELLDEF(cmat, ccell);
CELLDEF(zmat, zcell);
CELLDEF(dmat, dcell);
CELLDEF(smat, scell);
CELLDEF(lmat, lcell);

CELLDEF(dsp, dspcell);
CELLDEF(ssp, sspcell);
CELLDEF(csp, cspcell);
CELLDEF(zsp, zspcell);

CELLDEF(ccell, cccell);
CELLDEF(zcell, zccell);
CELLDEF(dcell, dccell);
CELLDEF(scell, sccell);
CELLDEF(lcell, iccell);

CELLDEF(cccell, ccccell);
CELLDEF(zccell, zcccell);
CELLDEF(dccell, dcccell);
CELLDEF(sccell, scccell);
CELLDEF(iccell, icccell);

CELLDEF(map_t, mapcell);
CELLDEF(rmap_t, rmapcell);
CELLDEF(loc_t, loccell);

CELLDEF(mapcell, mapccell);
CELLDEF(rmapcell, rmapccell);
CELLDEF(loccell, locccell);

//CELLDEF(cell, cell);
/*
    Reshape array without altering the number of elements.
*/
#define reshape(in, nx_, ny_) \
({\
  int ans=0;\
  if(in){\
    long nx__=nx_;/*preserve input value*/\
    long ny__=ny_;\
    if(PN(in)==nx__*ny__){\
        in->nx=nx__;\
        in->ny=ny__;\
    }else if(ny__==0 && nx__>0 && PN(in)%nx__==0){ /* change number of rows*/\
        in->ny=PN(in)/nx__;\
        in->nx=nx__;\
    }else if(nx__==0 && ny__==1){ /* change to vector */\
        in->nx=PN(in);\
        in->ny=1;\
    } else{\
        error("Must not change number of elements, aborted.\n");\
        ans=-1;\
    }\
  };ans;\
})
#undef ARR
#undef CELLARR
#undef MATARR
static inline int iscell(const void* id){
    return id?(*(const M_ID*)id==MCC_ANY):0;
    //return (((magic)&0x6410)==0x6410 || ((magic)&0x6420) == 0x6420);
}
/*A method to simulate operator overloading for indexing arrys*/
#if DEBUG
static inline long index_1d(long i, long nx, long ny){
    if(i<0 || i>=nx*ny){
        error("Index %ld is out of range for (%ld,%ld) array\n", i, nx, ny);
    }
    return i;
}
static inline long index_2d(long ix, long iy, long nx, long ny){
    if(ix<0 || ix>=nx || iy<0 || iy>=ny){
        error("Index (%ld,%ld) is out of range for (%ld,%ld) array\n", ix, iy, nx, ny);
    }
    return ix+iy*nx;
}
static inline long index_col(long iy, long nx, long ny){
    if(iy<0 || iy>=ny){
        error("Column %ld is out of range for (%ld,%ld) array\n", iy, nx, ny);
    }
    return iy*nx;
}
#else
#define index_1d(i,nx,ny) (i)
#define index_2d(ix,iy,nx,ny) ((ix)+(iy)*(nx))
#define index_col(iy,nx,ny) ((iy)*(nx))
#endif
/**
    \def P()
    Used to obtain elements of an array \c P(A,i) or \c P(A,i,j). Returns pointer if no index is given \c P(A).

    \def PCOL(A,iy)
    Used to obtain pointer to columns iy of array A

    \def PR(A,ix,iy)
    Calls P() with wrapping of index

    \def PCOLR(A, iy)
    Calls PCOL(A, iy) with wrapping of index

    \def PN(A)
    Return number of elements

    \def NX(A)
    Return number of rows (inner dimension)

    \def NY(A)
    Return number of cols (outer dimension)
 */

#define P0(A)       (A)->p //cannot do test here, bothers memcpy
#define P1(A,i)     (A)->p[index_1d((i),        (A)->nx, (A)->ny)]
#define P2(A,ix,iy) (A)->p[index_2d((ix), (iy), (A)->nx, (A)->ny)]
#define P3(Ac,icx,ix,iy) P2(P1(Ac,icx),ix,iy)
#define P4(Ac,icx,icy,ix,iy) P2(P2(Ac,icx,icy),ix,iy)

#define P_GET5(_0,_1,_2,_3,_4,_5,NAME,...) NAME
#define P(...) P_GET5(_0,__VA_ARGS__,P4,P3,P2,P1,P0)(__VA_ARGS__)
#define PCOL(A,iy) ((A)->p+index_col((iy), (A)->nx, (A)->ny))

//Define indexing using wrapping. 
#define P_GET3(_0,_1,_2,_3,NAME,...) NAME
#define PR1(A,ix)    P1((A), (((A)->nx==1 && (A)->ny==1)?0:ix))
#define PR2(A,ix,iy) P2((A), ((A)->nx==1?0:ix), ((A)->ny==1?0:iy))
#define PR(...) P_GET3(_0,__VA_ARGS__,PR2,PR1,P0)(__VA_ARGS__)
#define PCOLR(A,iy) PCOL((A),((A)->ny==1?0:iy))

//Return Number of elements
#define PN0(A)       ((A)?((A)->nx*(A)->ny):0)
#define NX0(A)       ((A)?(A)->nx:0)
#define NY0(A)       ((A)?(A)->ny:0)

//Do not do checking which slows done loops a lot
#define PN1(A,i)     PN0(P1(A,i))
#define PN2(A,ix,iy) PN0(P2(A,ix,iy))
#define PN(...)      P_GET3(_0,__VA_ARGS__,PN2,PN1,PN0)(__VA_ARGS__)

#define NX1(A,i)     NX0(P1(A,i))
#define NX2(A,ix,iy) NX0(P2(A,ix,iy))
#define NX(...)      P_GET3(_0,__VA_ARGS__,NX2,NX1,NX0)(__VA_ARGS__)

#define NY1(A,i)     NY0(P1(A,i))
#define NY2(A,ix,iy) NY0(P2(A,ix,iy))
#define NY(...)      P_GET3(_0,__VA_ARGS__,NY2,NY1,NY0)(__VA_ARGS__)

//Check non-empty array
#define NE0(A)       ((A) && (A)->nx && (A)->ny)
#define NE1(A,i)     ((A) && NE0(P(A,i)))
#define NE2(A,ix,iy) ((A) && NE0(P(A,ix,iy)))
#define NE(...)      P_GET3(_0,__VA_ARGS__,NE2,NE1,NE0)(__VA_ARGS__)
#endif
