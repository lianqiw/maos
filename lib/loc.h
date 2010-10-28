/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LOC_H
#define AOS_LOC_H
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"


/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y */
typedef struct MAP_T{
    double *p;      /**<The OPD*/
    long nx;        /**<Number of points along x*/
    long ny;        /**<Number of points along y*/
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x and y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx, vy;  /**Wind velocity. Useful for atmospheric grid*/
    long shm;       /**Records the length of the memory mmaped of positive. -1 means it is part of shared shm.*/
} MAP_T;

/**
   Map with different x/y sampling.
*/
typedef struct RECTMAP_T{
    double *p;      /**<The OPD*/
    long nx;        /**<Number of points along x*/
    long ny;        /**<Number of points along y*/
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x (first dimension)*/
    double dy;      /**<Sampling along y (second dimension)*/
    double txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel;    /**<Effective focal length of the telescope*/
    double fexit;   /**<The distance between the exit pupil and the focus*/
    double fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
}RECTMAP_T;

/**
   map of locs
*/
typedef struct LOCMAP_T{
    struct LOC_T *loc;
    int nloc;      /**<record nloc. to check against loc to see if loc has changed.*/
    long *p;       /**<The map, of size nx*ny*/
    double ox;     /**<Origin of the map along x*/
    double oy;     /**<Origin of the map along y*/
    int nx;        /**<Number of points along x*/
    int ny;        /**<Number of points along y*/
    int npad;      /**<Padding along the boundary. just for checking. no need in computation.*/
}LOCMAP_T;

/**
   Struct for loc like plocs, xloc, aloc etc.
*/
typedef struct LOC_T{
    double *locx;  /**< x coordinates of each point*/
    double *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling*/
    LOCMAP_T *map; /**< point to the map used for identifying neihboring points.*/
}LOC_T;

/**
   Store starting x,y for each col
*/
typedef struct LOCSTATCOL_T{
    double xstart; /**<starting x of this column*/
    double ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
}LOCSTATCOL_T;

/**
   Stores array of LOCSTATCOL_T
*/
typedef struct LOCSTAT_T{
    LOCSTATCOL_T *cols; /**<Information about each column*/
    long   ncol;        /**<Number of consecutive columns found*/
    double dx;          /**<Sampling of the grid*/
}LOCSTAT_T;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that PTS_T can be
   cast to LOC_T

   2009-12-08: Bug found: I set the origx, origy to the lower
   left point in the subaperture instead of the true
   edge. This is good for raytracing. But I used it as the
   saloc orig for mkg, which is really wrong.  use
   powfs->saloc as the orig of subapertures. PTS just for ray
   tracing.
*/
typedef struct PTS_T{
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    double dsa;    /**<side length of subaperture*/
    LOCMAP_T *map; /**<treat PTS_T as LOC_T and compute the MAP*/
    double *area;  /**<subaperture area, sum(amp^2)*/
    double dx;     /**<sampling of points in each subaperture*/
    int nx;        /**<number of col and row of points per subaperture*/
}PTS_T;
int *loc_create_embed(int *nembed, const LOC_T *loc);
void loc_create_map_npad(LOC_T *loc, int npad);
void loc_create_map(LOC_T *loc);

LOC_T * map2loc(double dx, long nx, long ny, 
		double ox, double oy, double *map);
long *map2embed(long nx, long ny, double *map);
LOC_T * sqmap2loc(MAP_T *amp);
long* sqmap2embed(MAP_T *amp);
void rectmapfree_do(RECTMAP_T *map);
#define rectmapfree(A) ({rectmapfree_do(A);A=NULL;})
void sqmapfree_do(MAP_T *map);
#define sqmapfree(A) ({sqmapfree_do(A);A=NULL;})
void sqmaparrfree_do(MAP_T **map, int nmap);
#define sqmaparrfree(A,B) ({sqmaparrfree_do(A,B);A=NULL;})
void loc_free_map(LOC_T *loc);

void locfree_do(LOC_T *loc);
#define locfree(A) ({locfree_do(A);A=NULL;})
void ptsfree_do(PTS_T *pts);
#define ptsfree(A) ({ptsfree_do(A);A=NULL;})
void locarrfree_do(LOC_T **loc, int nloc);
#define locarrfree(A,B) ({locarrfree_do(A,B);A=NULL;})

void locstatfree_do(LOCSTAT_T *locstat);
#define locstatfree(A) ({locstatfree_do(A);A=NULL;})

int loccenter(LOC_T *loc);

void loc_calc_ptt(double *out, double *coeffout, 
		  const LOC_T *loc, const double ipcc, 
	       const dmat *imcc, const double *amp, const double *opd);
void loc_calc_mod(double *out, double *coeffout, 
		 const dmat *mod, const double *amp, double *opd);
dmat *loc_mcc_ptt(const LOC_T *loc, const double *amp);
dcell *pts_mcc_ptt(const PTS_T *pts, const double *amp);
void loc_remove_ptt(double *opd, const double *ptt, const LOC_T *loc);
void loc_add_ptt(double *opd, const double *ptt, const LOC_T *loc);
void pts_ztilt(double *out, const PTS_T *pts, const dcell *imcc,
	       const double *amp, const double *opd);
LOC_T *mk1dloc_vec(double *x, long nx);
LOC_T *mk1dloc(double x0, double dx, long nx);
LOC_T *mksqloc_auto(long nx, long ny, double dx);
LOC_T *mksqloc_map(MAP_T*map);
LOC_T *mksqloc(long nx, long ny, double dx, double ox, double oy);
LOC_T *mksqlocrot(long nx, long ny, double dx, 
		  double ox, double oy, double theta);
LOC_T *mkcirloc(double d, double dx);
LOC_T *mkcirloc_amp(double** restrict ampout, LOCSTAT_T *cstat, 
		    MAP_T* ampin, double dtel, double dx, int cropamp);
LOCSTAT_T *mklocstat(const LOC_T *loc);
void loccircle(double *phi,LOC_T *loc,double cx,double cy,double r,double val);
void locellipse(double *phi,LOC_T *loc,double cx,double cy,
		double rx,double ry,double val);
void loc_reduce_spcell(LOC_T *loc, spcell *sp, int dim, int cont);
void loc_reduce_sp(LOC_T *loc, dsp *sp, int dim, int cont);

dmat* loc_zernike(LOC_T *loc, double R, int nr);
void loc_add_focus(double *opd, LOC_T *loc, double val);
dmat *loc2mat(LOC_T *loc,int piston);
LOC_T *pts2loc(PTS_T *pts);
void locrot(LOC_T *loc, const double theta);
LOC_T *locdup(LOC_T *loc);
LOC_T *loctransform(LOC_T *loc, dmat **coeff);
void loc_nxny(long *nx, long *ny, const LOC_T *loc);
#endif
