/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
typedef struct map_t{
    double *p;      /**<The OPD*/
    long nx;        /**<Number of points along x*/
    long ny;        /**<Number of points along y*/
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x and y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx, vy;  /**Wind velocity. Useful for atmospheric grid*/
    long shm;       /**Records the length of the memory mmaped of positive. -1 means it is part of shared shm.*/
} map_t;

/**
   Map with different x/y sampling.
*/
typedef struct rectmap_t{
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
}rectmap_t;

/**
   map of locs
*/
typedef struct locmap_t{
    struct loc_t *loc;
    int nloc;      /**<record nloc. to check against loc to see if loc has changed.*/
    long *p;       /**<The map, of size nx*ny*/
    double ox;     /**<Origin of the map along x*/
    double oy;     /**<Origin of the map along y*/
    int nx;        /**<Number of points along x*/
    int ny;        /**<Number of points along y*/
    int npad;      /**<Padding along the boundary. just for checking. no need in computation.*/
}locmap_t;

/**
   Struct for loc like plocs, xloc, aloc etc.
*/
typedef struct loc_t{
    double *locx;  /**< x coordinates of each point*/
    double *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling*/
    locmap_t *map; /**< point to the map used for identifying neihboring points.*/
}loc_t;

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
    long   ncol;        /**<Number of consecutive columns found*/
    double dx;          /**<Sampling of the grid*/
}locstat_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be
   cast to loc_t

   2009-12-08: Bug found: I set the origx, origy to the lower
   left point in the subaperture instead of the true
   edge. This is good for raytracing. But I used it as the
   saloc orig for mkg, which is really wrong.  use
   powfs->saloc as the orig of subapertures. PTS just for ray
   tracing.
*/
typedef struct pts_t{
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    double dsa;    /**<side length of subaperture*/
    locmap_t *map; /**<treat pts_t as loc_t and compute the MAP*/
    double *area;  /**<subaperture area, sum(amp^2)*/
    double dx;     /**<sampling of points in each subaperture*/
    int nx;        /**<number of col and row of points per subaperture*/
}pts_t;
int *loc_create_embed(int *nembed, const loc_t *loc);
void loc_create_map_npad(loc_t *loc, int npad);
void loc_create_map(loc_t *loc);

loc_t * map2loc(double dx, long nx, long ny, 
		double ox, double oy, double *map);
long *map2embed(long nx, long ny, double *map);
loc_t * sqmap2loc(map_t *amp);
long* sqmap2embed(map_t *amp);
void rectmapfree_do(rectmap_t *map);
#define rectmapfree(A) ({rectmapfree_do(A);A=NULL;})
void sqmapfree_do(map_t *map);
#define sqmapfree(A) ({sqmapfree_do(A);A=NULL;})
void sqmaparrfree_do(map_t **map, int nmap);
#define sqmaparrfree(A,B) ({sqmaparrfree_do(A,B);A=NULL;})
void loc_free_map(loc_t *loc);

void locfree_do(loc_t *loc);
#define locfree(A) ({locfree_do(A);A=NULL;})
void ptsfree_do(pts_t *pts);
#define ptsfree(A) ({ptsfree_do(A);A=NULL;})
void locarrfree_do(loc_t **loc, int nloc);
#define locarrfree(A,B) ({locarrfree_do(A,B);A=NULL;})

void locstatfree_do(locstat_t *locstat);
#define locstatfree(A) ({locstatfree_do(A);A=NULL;})

int loccenter(loc_t *loc);
loc_t *locnew(long nloc);
void loc_calc_ptt(double *out, double *coeffout, 
		  const loc_t *loc, const double ipcc, 
	       const dmat *imcc, const double *amp, const double *opd);
void loc_calc_mod(double *out, double *coeffout, 
		 const dmat *mod, const double *amp, double *opd);
dmat *loc_mcc_ptt(const loc_t *loc, const double *amp);
dcell *pts_mcc_ptt(const pts_t *pts, const double *amp);
void loc_remove_ptt(double *opd, const double *ptt, const loc_t *loc);
void loc_add_ptt(double *opd, const double *ptt, const loc_t *loc);
void pts_ztilt(double *out, const pts_t *pts, const dcell *imcc,
	       const double *amp, const double *opd);
loc_t *mk1dloc_vec(double *x, long nx);
loc_t *mk1dloc(double x0, double dx, long nx);
loc_t *mksqloc_auto(long nx, long ny, double dx);
loc_t *mksqloc_map(map_t*map);
loc_t *mksqloc(long nx, long ny, double dx, double ox, double oy);
loc_t *mksqlocrot(long nx, long ny, double dx, 
		  double ox, double oy, double theta);
loc_t *mkcirloc(double d, double dx);
loc_t *mkcirloc_amp(double** restrict ampout, locstat_t *cstat, 
		    map_t* ampin, double dtel, double dx, int cropamp);
locstat_t *mklocstat(const loc_t *loc);
void loccircle(double *phi,loc_t *loc,double cx,double cy,double r,double val);
void locellipse(double *phi,loc_t *loc,double cx,double cy,
		double rx,double ry,double val);
void loc_reduce_spcell(loc_t *loc, spcell *sp, int dim, int cont);
void loc_reduce_sp(loc_t *loc, dsp *sp, int dim, int cont);

dmat* loc_zernike(loc_t *loc, double R, int nr);
void loc_add_focus(double *opd, loc_t *loc, double val);
dmat *loc2mat(loc_t *loc,int piston);
loc_t *pts2loc(pts_t *pts);
void locrot(loc_t *loc, const double theta);
loc_t *locdup(loc_t *loc);
loc_t *loctransform(loc_t *loc, dmat **coeff);
void loc_nxny(long *nx, long *ny, const loc_t *loc);
map_t *mapnew(long nx, long ny, double dx, double *p);
void mapcircle(map_t *map, double r, double val);
#endif
