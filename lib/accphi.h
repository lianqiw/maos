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

#ifndef AOS_ACCPHI_H
#define AOS_ACCPHI_H


#ifndef AOS_LIB_THREAD_H
typedef struct thread_t thread_t;
#endif
#include "loc.h"

/**
   Unified data structure for automatic selection of propagation.
*/
typedef struct PROPDATA_T{
    /*Input.  */
    map_t *mapin;
    /*Or */
    loc_t *locin; const double *phiin;
    
    /*Output */
    map_t *mapout;
    /*Or */
    double *phiout;
    /*Combined with */
    const pts_t *ptsout;
    /*Or */
    const loc_t *locout; const double *ampout;
    /*Or  */
    const locstat_t *ostat;

    /*Constant displacement */
    double displacex0, displacey0;
    /*Time step dependent displacement */
    double displacex1, displacey1;
    /*scale of coordinate */
    double scale;
    /*scale of value: */
    double alpha;
    /*Options: */
    int cubic; double cubic_iac; /*for cubic interpolation. */
    int wrap;
    int nooptim;/*disable optim. */
    int index;
}PROPDATA_T;

void prop(thread_t *data);/*A unified wrapper */
void prop_index(PROPDATA_T *propdata);/*A unified wrapper */

#define ARGIN_GRID						\
    const map_t *mapin /**<[in] OPD defind on a square grid*/

#define ARGIN_NONGRID							\
    loc_t *locin,     /**<[in] Coordinate of iregular source grid*/	\
	const double *phiin/**<[in] Input OPD defined in locin*/

#define ARGOUT_LOC							\
    const loc_t *locout, /**<[in] Coordinate of irregular output grid*/	\
	const double *ampout,/**<[in] Amplitude defined on locout. skip point of amp is 0*/ \
	double* phiout      /**<[in,out] Output OPD defined in locout*/
#define ARGOUT_PTS							\
    const pts_t *pts,   /**<[in] coordinate of destination grid*/	\
	const double *ampout,/**<[in] Amplitude defined on locout. skip point of amp is 0*/ \
	double *phiout     /**<[in,out] OPD defined on locout*/

#define ARGOUT_MAP\
    map_t *mapout    /**<[in,out] Output OPD defined in a square grid*/

#define ARGOUT_STAT							\
    const locstat_t *ostat,/*<[in] statics of columns in a loc_t*/	\
	double *phiout /**<[in, out] Output OPD defined on ostat*/

#define ARG_PROP							\
    const double alpha,     /**<[in] scaling of OPD*/			\
	double displacex, /**<[in] displacement of the ray */	\
	double displacey, /**<[in] displacement of the ray */	\
	const double scale     /**<[in] scaling of the beam diameter (cone)*/

void prop_grid_map(ARGIN_GRID, ARGOUT_MAP, ARG_PROP, int wrap, long start, long end);
void prop_grid_pts(ARGIN_GRID, ARGOUT_PTS, ARG_PROP, int wrap, long sastart, long saend);
void prop_grid    (ARGIN_GRID, ARGOUT_LOC, ARG_PROP, int wrap, long start, long end);

void prop_nongrid(ARGIN_NONGRID, ARGOUT_LOC, ARG_PROP, long start, long end);
void prop_nongrid_map(ARGIN_NONGRID, ARGOUT_MAP, ARG_PROP, long start, long end);
void prop_nongrid_pts(ARGIN_NONGRID, ARGOUT_PTS, ARG_PROP, long start, long end);
/*
  A few cubic propagations.
 */
void prop_grid_cubic    (ARGIN_GRID, ARGOUT_LOC, ARG_PROP, double cubic_iac, long start, long end);
void prop_grid_pts_cubic(ARGIN_GRID, ARGOUT_PTS, ARG_PROP, double cubic_iac, long start, long end);
void prop_grid_map_cubic(ARGIN_GRID, ARGOUT_MAP, ARG_PROP, double cubic_iac, long start, long end);
void prop_nongrid_cubic    (ARGIN_NONGRID, ARGOUT_LOC, ARG_PROP, double cubic_iac, long start, long end);
void prop_nongrid_pts_cubic(ARGIN_NONGRID, ARGOUT_PTS, ARG_PROP, double cubic_iac, long start, long end);
void prop_nongrid_map_cubic(ARGIN_NONGRID, ARGOUT_MAP, ARG_PROP, double cubic_iac, long start, long end);

/**
   Propagate OPD defines on grid mapin to locstat_t that stores the starting
   point of each column.  alpha is the scaling of data. displacex, displacy is
   the displacement of the center of the beam on the input grid.  scale is the
   cone effect.*/
void prop_grid_stat(const map_t *mapin, const locstat_t *ostat, 
		    double *phiout, double alpha,
		    double displacex, double displacey,
		    double scale, int wrap,
		    long colstart, long colend);
/**
   Do the reverse propagation of prop_grid_stat. If prop_grid_stat does y+=H*x;
This just does x+=H'*y; */
void prop_grid_stat_transpose(map_t *mapin, const locstat_t *ostat, 
			      const double *phiout, double alpha,
			      double displacex, double displacey,
			      double scale, int wrap,
			      long colstart, long colend);
/*
  the following routine is used to du down sampling by doing binning ray tracing.
  locout is coarse sampling, locint is fine sampling.
  phiin->phiout.
*/
void prop_nongrid_bin(const loc_t *locin,
		      const double* phiin,
		      loc_t *locout, 
		      const double *ampout,
		      double* phiout,
		      ARG_PROP);
#endif
