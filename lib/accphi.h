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

#ifndef AOS_ACCPHI_H
#define AOS_ACCPHI_H
#include "loc.h"
#include "common.h"
/*Unified data structure for automatic selection of propagation.*/
typedef struct PROPDATA_T{
    //Input. 
    MAP_T *mapin;
    //Or
    LOC_T *locin; const double *phiin;
    
    //Output
    MAP_T *mapout;
    //Or
    double *phiout;
    //Combined with
    const PTS_T *ptsout;
    //Or
    const LOC_T *locout; const double *ampout;
    //Or 
    const LOCSTAT_T *ostat;

    //Displacement, scale of coordinate, and alpha of value:
    double displacex, displacey, scale, alpha;
    //Options:
    int cubic; double cubic_iac; //for cubic interpolation.
    int wrap;
    int nooptim;//disable optim.
}PROPDATA_T;
/*Structure that wraps PROPDATA_T for threading */
/*typedef struct PROP_T{ //moved to thread_t
    PROPDATA_T *propdata;
    long start;//starting point
    long end;//end point
    long step;
}PROP_T;
*/
void prop(thread_t *data);//A unified wrapper
void prop_grid_grid(const MAP_T *mapin, MAP_T *mapout,
		    double alpha,
		    double displacex, double displacey, 
		    double scale, int wrap );
void prop_grid_pts(const MAP_T *mapin, const PTS_T *pts, 
		   double *phiout0, double alpha,
		   double displacex, double displacey, 
		   double scale, int wrap, 
		   long sastart, long saend, int optim);
void prop_grid(const MAP_T *mapin, const LOC_T *locout, 
	       double *phiout, double alpha,
	       double displacex, double displacey,
	       double scale, int wrap,
	       long start, long end);
void prop_grid_stat(const MAP_T *mapin, const LOCSTAT_T *ostat, 
		    double *phiout0, double alpha,
		    double displacex, double displacey,
		    double scale, int wrap,
		    long colstart, long colend, int optim);
void prop_nongrid(LOC_T *locin, const double* phiin, 
		  const LOC_T *locout,const double *amp,
		  double* phiout, double alpha,
		  double displacex, double displacey,
		  double scale, long start, long end);
void prop_nongrid_map(LOC_T *locin, const double *phiin,
		      MAP_T *mapout, double alpha,
		      double displacex, double displacey,
		      double scale, long start, long end, long step);
void prop_nongrid_pts(LOC_T *locin, const double *phiin,
		      const PTS_T *pts,const double *ampout,
		      double *phiout, double alpha,
		      double displacex, double displacey,
		      double scale, long start, long end, long step);
void prop_nongrid_cubic(LOC_T *locin, const  double* phiin, 
			const LOC_T *locout, const double *amp,
			double* phiout,  double alpha,
			double displacex, double displacey,
			double scale, double cubic_iac, long start, long end);
void prop_nongrid_pts_cubic(LOC_T *locin, const double* phiin, 
			    const PTS_T *pts, const double *ampout, 
			    double* phiout, double alpha,
			    double displacex, double displacey,
			    double scale, double cubic_iac, 
			    long start, long end, long step);
void prop_nongrid_map_cubic(LOC_T *locin, const double* phiin, 
			    MAP_T* mapout, double alpha,
			    double displacex, double displacey,
			    double scale, double cubic_iac, 
			    long start, long end, long step);
/*
  the following routine is used to du down sampling by doing *reverse* ray tracing.
  locin is coarse sampling, locout is fine sampling. phiin is the unknown
*/
void prop_nongrid_reverse(LOC_T *locin, double* phiin, 
			  const LOC_T *locout,
			  const double *ampout, 
			  const double* phiout, double alpha,
			  double displacex, double displacey,
			  double scale);
#endif
