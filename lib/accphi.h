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

#ifndef AOS_LIB_ACCPHI_H
#define AOS_LIB_ACCPHI_H
#include "../math/mathdef.h"

/**
   \file accphi.h

   Contains ray tracing routines optimized for different input/output
   formats. Notice that the OPDs are accumulated.
*/
/**
   Unified data structure for automatic selection of propagation.
*/
typedef struct PROPDATA_T{
    /*Input.  */
    map_t *mapin;
    /*Or */
    loc_t *locin; const real *phiin;
    
    /*Output */
    map_t *mapout;
    /*Or */
    real *phiout;
    /*Combined with */
    const pts_t *ptsout;
    /*Or */
    const loc_t *locout; 
    /*Or  */
    const locstat_t *ostat;

    /*Constant displacement */
    real displacex0, displacey0;
    /*Time step dependent displacement */
    real displacex1, displacey1;
    /*scale of coordinate */
    real scale;
    /*scale of value: */
    real alpha;
    int wrap;
    int nooptim;/*disable optim. */
    int index;
}PROPDATA_T;
void prop(thread_t *data);/*A unified wrapper */

#define ARGIN_GRID						\
    const map_t *mapin /**<[in] OPD defind on a square grid*/
#define ARGIN_GRID2 mapin
#define ARGIN_NONGRID							\
    const loc_t *locin,     /**<[in] Coordinate of iregular source grid*/	\
    const real *phiin/**<[in] Input OPD defined in locin*/
#define ARGIN_NONGRID2 locin, phiin
#define ARGOUT_LOC							\
    const loc_t *locout, /**<[in] Coordinate of irregular output grid*/	\
    real* phiout      /**<[in,out] Output OPD defined in locout*/
#define ARGOUT_LOC2 locout, phiout
#define ARGOUT_PTS							\
    const pts_t *pts,   /**<[in] coordinate of destination grid*/	\
    real *phiout     /**<[in,out] OPD defined on locout*/
#define ARGOUT_PTS2 pts, phiout
#define ARGOUT_MAP							\
    map_t *mapout    /**<[in,out] Output OPD defined in a square grid*/
#define ARGOUT_MAP2 mapout
#define ARGOUT_STAT							\
    const locstat_t *ostat,/**<[in] statics of columns in a loc_t*/	\
    real *phiout /**<[in, out] Output OPD defined on ostat*/
#define ARGOUT_STAT2 ostat, phiout
#define ARG_PROP							\
    const real alpha,     /**<[in] scaling of OPD*/			\
    real displacex, /**<[in] displacement of the ray */		\
    real displacey, /**<[in] displacement of the ray */		\
    const real scale     /**<[in] scaling of the beam diameter (cone)*/
#define ARG_PROP2 alpha, displacex, displacey, scale

#define ARG_INDEX                               \
    long start, /**<[in] start index*/           \
    long end /**<[ind] end index*/

#define ARG_PROP_WRAP                           \
    ARG_PROP,                                   \
    int wrap, /**<[in] warp around */           \
    ARG_INDEX

#define ARG_PROP_NOWRAP                         \
    ARG_PROP,                                   \
    ARG_INDEX

#define ARG_PROP_CUBIC                                          \
    ARG_PROP,                                                           \
    real cubic_iac, /**<[in] inter-actuator coupling coeffcient*/       \
    ARG_INDEX

void prop_grid     (ARGIN_GRID, ARGOUT_LOC, ARG_PROP_WRAP);
void prop_grid_map (ARGIN_GRID, ARGOUT_MAP, ARG_PROP_WRAP);
void prop_grid_pts (ARGIN_GRID, ARGOUT_PTS, ARG_PROP_WRAP);
void prop_grid_stat(ARGIN_GRID, ARGOUT_STAT,ARG_PROP_WRAP);

void prop_nongrid    (ARGIN_NONGRID, ARGOUT_LOC, ARG_PROP_NOWRAP);
void prop_nongrid_map(ARGIN_NONGRID, ARGOUT_MAP, ARG_PROP_NOWRAP);
void prop_nongrid_pts(ARGIN_NONGRID, ARGOUT_PTS, ARG_PROP_NOWRAP);

/*
  A few cubic spline propagations.
*/
void prop_grid_cubic     (ARGIN_GRID, ARGOUT_LOC, ARG_PROP_CUBIC);
void prop_grid_map_cubic (ARGIN_GRID, ARGOUT_MAP, ARG_PROP_CUBIC);
void prop_grid_pts_cubic (ARGIN_GRID, ARGOUT_PTS, ARG_PROP_CUBIC);
void prop_grid_stat_cubic(ARGIN_GRID, ARGOUT_STAT,ARG_PROP_CUBIC);

void prop_nongrid_cubic    (ARGIN_NONGRID, ARGOUT_LOC, ARG_PROP_CUBIC);
void prop_nongrid_pts_cubic(ARGIN_NONGRID, ARGOUT_PTS, ARG_PROP_CUBIC);
void prop_nongrid_map_cubic(ARGIN_NONGRID, ARGOUT_MAP, ARG_PROP_CUBIC);

void prop_grid_map_transpose(map_t *mapin, const map_t *mapout, 
			     real alpha, real displacex, real displacey, real scale,
			     int wrap, long start, long end);
/**
   Do the reverse propagation of prop_grid_stat. If prop_grid_stat does y+=H*x;
   This just does x+=H'*y; */
void prop_grid_stat_transpose(map_t *mapin, const locstat_t *ostat, 
			      const real *phiout, real alpha,
			      real displacex, real displacey,
			      real scale, int wrap,
			      long colstart, long colend);
/*
  the following routine is used to do down sampling by doing binning ray tracing.
  locout is coarse sampling, locint is fine sampling.
  phiin->phiout.
*/
void prop_nongrid_bin(const loc_t *locin,
		      const real* phiin,
		      loc_t *locout, 
		      real* phiout,
		      ARG_PROP);

#endif
