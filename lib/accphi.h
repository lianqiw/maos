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
    //Input. 
    map_t *mapin;
    //Or
    loc_t *locin; const double *phiin;
    
    //Output
    map_t *mapout;
    //Or
    double *phiout;
    //Combined with
    const pts_t *ptsout;
    //Or
    const loc_t *locout; const double *ampout;
    //Or 
    const locstat_t *ostat;

    //Constant displacement
    double displacex0, displacey0;
    //Time step dependent displacement
    double displacex1, displacey1;
    //scale of coordinate
    double scale;
    //scale of value:
    double alpha;
    //Options:
    int cubic; double cubic_iac; //for cubic interpolation.
    int wrap;
    int nooptim;//disable optim.
    int index;
}PROPDATA_T;

void prop(thread_t *data);//A unified wrapper
void prop_index(PROPDATA_T *propdata);//A unified wrapper
void prop_grid_grid(const map_t *mapin, map_t *mapout,
		    double alpha,
		    double displacex, double displacey, 
		    double scale, int wrap);
void prop_grid_pts(const map_t *mapin, const pts_t *pts, 
		   double *phiout0, double alpha,
		   double displacex, double displacey, 
		   double scale, int wrap, 
		   long sastart, long saend);
void prop_grid(const map_t *mapin, const loc_t *locout, 
	       double *phiout, double alpha,
	       double displacex, double displacey,
	       double scale, int wrap,
	       long start, long end);
void prop_grid_stat(const map_t *mapin, const locstat_t *ostat, 
		    double *phiout0, double alpha,
		    double displacex, double displacey,
		    double scale, int wrap,
		    long colstart, long colend);
void prop_nongrid(loc_t *locin, const double* phiin, 
		  const loc_t *locout,const double *amp,
		  double* phiout, double alpha,
		  double displacex, double displacey,
		  double scale, long start, long end);
void prop_nongrid_map(loc_t *locin, const double *phiin,
		      map_t *mapout, double alpha,
		      double displacex, double displacey,
		      double scale, long start, long end);
void prop_nongrid_pts(loc_t *locin, const double *phiin,
		      const pts_t *pts,const double *ampout,
		      double *phiout, double alpha,
		      double displacex, double displacey,
		      double scale, long start, long end);
void prop_nongrid_cubic(loc_t *locin, const  double* phiin, 
			const loc_t *locout, const double *amp,
			double* phiout,  double alpha,
			double displacex, double displacey,
			double scale, double cubic_iac, long start, long end);
void prop_nongrid_pts_cubic(loc_t *locin, const double* phiin, 
			    const pts_t *pts, const double *ampout, 
			    double* phiout, double alpha,
			    double displacex, double displacey,
			    double scale, double cubic_iac, 
			    long start, long end);
void prop_nongrid_map_cubic(loc_t *locin, const double* phiin, 
			    map_t* mapout, double alpha,
			    double displacex, double displacey,
			    double scale, double cubic_iac, 
			    long start, long end);
/*
  the following routine is used to du down sampling by doing *reverse* ray tracing.
  locin is coarse sampling, locout is fine sampling. phiin is the unknown
*/
void prop_nongrid_reverse(loc_t *locin, double* phiin, 
			  const loc_t *locout,
			  const double *ampout, 
			  const double* phiout, double alpha,
			  double displacex, double displacey,
			  double scale);
#endif
