#ifndef AOS_LIB_TYPE_H
#define AOS_LIB_TYPE_H
#include "../math/type.h"

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
typedef struct rectmap_t{
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
}rectmap_t;

/**
   map of locs. Convert any coordinate (x,y) to corresponding iloc that has the
   coordinate.  */
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
    double *restrict locx;  /**< x coordinates of each point*/
    double *restrict locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling along x*/
    double dy;     /**< Sampling along y*/ 
    locmap_t *restrict map; /**< point to the map used for identifying neihboring points.*/
    locstat_t *stat;/**<points to column statistics*/
}loc_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be
   cast to loc_t
*/
typedef struct pts_t{
    long id;
    double *restrict origx; /**<The x origin of each subaperture*/
    double *restrict origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    union{
	double dsa;    /**<side length of subaperture*/
	double dsax;   /**<side length of subaperture*/
    };
    double dsay;   /**<side length of subaperture*/
    locmap_t *restrict map; /**<treat pts_t as loc_t and compute the MAP*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    double dx;     /**<sampling of points in each subaperture*/
    double dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
    int nx;        /**<number points in each col or row per subaperture*/
}pts_t;

#endif
