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

#ifndef AOS_LIB_LOC_H
#define AOS_LIB_LOC_H
#include "../math/mathdef.h"
#include "type.h"
/**
   \file loc.h
   This file defines functions relates to pts_t, loc_t, map_t, etc.
 */

long *loc_create_embed(long *nembed, const loc_t *loc, int oversize);
void loc_create_map_npad(loc_t *loc, int npad);
void loc_create_map(loc_t *loc);
void loc_embed(map_t *dest, const loc_t *loc, const double *in);
void loc_embed_add(map_t *dest, const loc_t *loc, const double *in);
void loc_extract(dmat *dest, const loc_t *loc, map_t *in);
loc_t * map2loc(map_t *amp);
void rectmapfree_do(rectmap_t *map);
#define rectmapfree(A) ({rectmapfree_do(A);A=NULL;})
void mapfree_do(map_t *map);
#define mapfree(A) ({mapfree_do(A);A=NULL;})
void maparrfree_do(map_t **map, int nmap);
#define maparrfree(A,B) ({maparrfree_do(A,B);A=NULL;})
void loc_free_map(loc_t *loc);

void locfree_do(loc_t *loc);
#define locfree(A) ({locfree_do(A);A=NULL;})
void ptsfree_do(pts_t *pts);
#define ptsfree(A) ({ptsfree_do(A);A=NULL;})
void locarrfree_do(loc_t **loc, int nloc);
#define locarrfree(A,B) ({locarrfree_do(A,B);A=NULL;})

int loccenter(loc_t *loc);
loc_t *locnew(long nloc, double dx, double dy);
pts_t *ptsnew(long nsa, double dsax, double dsay, long nx, double dx, double dy);
void loc_calc_ptt(double *out, double *coeffout, 
		  const loc_t *loc, const double ipcc, 
	       const dmat *imcc, const double *amp, const double *opd);
void loc_calc_mod(double *out, double *coeffout, 
		 const dmat *mod, const double *amp, double *opd);
dmat *loc_mcc_ptt(const loc_t *loc, const double *amp);
dcell *pts_mcc_ptt(const pts_t *pts, const double *amp);
void loc_remove_ptt(double *opd, const double *ptt, const loc_t *loc);
void loc_add_ptt(double *opd, const double *ptt, const loc_t *loc);
void pts_ztilt(dmat **out, const pts_t *pts, const dcell *imcc,
	       const double *amp, const double *opd);
loc_t *mk1dloc_vec(double *x, long nx);
loc_t *mk1dloc(double x0, double dx, long nx);
loc_t *mksqloc_auto(long nx, long ny, double dx, double dy);
loc_t *mksqloc_map(map_t*map);
loc_t *mksqloc(long nx, long ny, double dx, double dy, double ox, double oy);
loc_t *mksqlocrot(long nx, long ny, double dx, double dy,
		  double ox, double oy, double theta);
void loc_create_stat_do(loc_t *loc);
#define loc_create_stat(loc) if(!loc->stat) loc_create_stat_do(loc);
void loc_free_stat(loc_t *loc);
void loccircle(double *phi,loc_t *loc,double cx,double cy,double r,double val);
void locannular(double *phi,loc_t *loc,double cx,double cy,double r,double rin,double val);
void locannularmask(double *phi,loc_t *loc,double cx,double cy,double r,double rin);
void locellipse(double *phi,loc_t *loc,double cx,double cy,
		double rx,double ry,double val);
void loc_reduce(loc_t *loc, dmat *amp, int cont, int **skipout);
void loc_reduce_spcell(loc_t *loc, spcell *sp, int dim, int cont);
void loc_reduce_sp(loc_t *loc, dsp *sp, int dim, int cont);

dmat* loc_zernike(loc_t *loc, double R, int nr);
void loc_add_focus(double *opd, loc_t *loc, double val);
dmat *loc2mat(loc_t *loc,int piston);
loc_t *pts2loc(pts_t *pts);
void locrot(loc_t *loc, const double theta);
loc_t *locdup(loc_t *loc);
loc_t *loctransform(loc_t *loc, const char *ps);
loc_t *locshift(const loc_t *loc, double sx, double sy);
void loc_nxny(long *nx, long *ny, const loc_t *loc);
map_t *mapnew(long nx, long ny, double dx, double dy, double *p);
map_t *mapnew2(map_t *A);
void mapcircle(map_t *map, double r, double val);
void mapcircle_symbolic(map_t *map, double r);
void map_d_din(map_t *map, double *d, double *din);
void locresize(loc_t *loc, long nloc);
#define ptsresize(pts, nsa) locresize((loc_t*)pts, nsa)
void dembed_locstat(dmat **out, double alpha, loc_t *loc, double *oin, double beta, int reverse);
void cembed_locstat(cmat **out, double alpha, loc_t *loc, double *oin, double beta, int reverse);
#endif
