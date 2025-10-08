/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "mathdef.h"
#include "type.h"
/**
   \file loc.h
   This file defines functions relates to pts_t, loc_t, map_t, etc.
 */
loc_t* loc_convert(dmat* A);
#define locwrite(out, A...) writebin(out, A)
#define locread(A...)    loc_convert(dread(A))
#define loccellread(A...) (loccell*)cellconvert(readbin_id(M_LOC, 1, A), (cell*(*)(cell*))loc_convert)
#define loccellnew (loccell*)cellnew
#define locccellnew (locccell*)cellnew
#define loccellfree(A) cellfree(A) //if(A){cellfree_do((cell*)A); A=NULL;};
#define locfree(A) cellfree(A) //if(A){cellfree_do((cell*)A);A=NULL;}
#define ptsfree(A) cellfree(A) //if(A){cellfree_do((cell*)A);A=NULL;}
//#define loccellwrite(A,B...) cellwrite((cell*)A,B)
lmat* loc_create_embed(long* nembed, const loc_t* loc, real oversize, int fftpad);
void loc_create_map_npad(const loc_t* loc, int npad, int nx, int ny);
void loc_create_map(const loc_t* loc);
/*Obtain an entry in the map, with boundary checking enforced*/
static inline long loc_map_get(map_t* map, long ix, long iy){
	if(ix>=0&&ix<map->nx&&iy>=0&&iy<map->ny){
		return (long)map->p[ix+iy*map->nx];
	} else{
		return 0;
	}
}
uint32_t lochash(const loc_t* loc, uint32_t key);
void loc_embed_do(anydmat _dest, const loc_t *loc, const_anydmat in, int add);
static inline void loc_embed(anydmat dest, const loc_t* loc, const_anydmat in){
  loc_embed_do(dest, loc, in, 0);
}
static inline void loc_embed_add(anydmat dest, const loc_t *loc, const_anydmat in){
  loc_embed_do(dest, loc, in, 1);
}
dcell* loc_embed2(const loc_t* loc, const dmat* arr);
void loc_embed_cell(dcell** dest, const loc_t* loc, const dcell* in);
void loc_extract(dmat* dest, const loc_t* loc, map_t* in);
loc_t* loc_from_map(map_t* amp, real thres);

void loc_free_map(loc_t* loc);
void locfree_do(cell* loc);

real loc_diam(const loc_t* loc);
int loccenter(const loc_t* loc);
loc_t* locnew(long nloc, real dx, real dy);
loc_t* locref(loc_t* in);
pts_t* ptsnew(long nsa, real dsax, real dsay, long nxsa, long nysa, real dx, real dy);
void loc_calc_ptt(real* out, real* coeffout,
	const loc_t* loc, real ipcc,
	const dmat* imcc, const real* amp, const real* opd);
void loc_calc_mod(real* out, real* coeffout,
	const dmat* mod, const real* amp, real* opd);
dmat* loc_mcc_ptt(const loc_t* loc, const real* amp);
dcell* pts_mcc_ptt(const pts_t* pts, const real* amp);
void loc_sub_ptt(dmat* opd, const real* ptt, const loc_t* loc);
void loc_add_ptt(dmat* opd, const real* ptt, const loc_t* loc);
void loc_remove_ptt(dmat *opd, const loc_t *loc, const real *amp, const dmat *imcc, int both);
real loc_remove_focus_grad(dmat *grad, const loc_t *saloc, real factor);
void pts_ztilt(dmat** out, const pts_t* pts, const dcell* imcc,
	const real* amp, const real* opd);
loc_t* mk1dloc_vec(real* x, long nx);
loc_t* mk1dloc(real x0, real dx, long nx);
loc_t* mksqloc_auto(long nx, long ny, real dx, real dy);
loc_t* mksqloc_map(const map_t* map);
loc_t* mksqloc(long nx, long ny, real dx, real dy, real ox, real oy);
loc_t* mkannloc(real D, real Din, real dx, real thres);

void loc_create_stat_do(loc_t* loc);
#define loc_create_stat(loc) if(!loc->stat) loc_create_stat_do(loc);
void loc_free_stat(loc_t* loc);
void loc_circle_add(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val);
void loc_circle_mul(dmat* phi, const loc_t* loc, real cx, real cy, real r, real rin, real val);
void loc_ellipse_add(dmat* phi, const loc_t* loc, real cx, real cy,real rx, real ry, real val);
void loc_reduce(loc_t* loc, dmat* amp, real thres, int cont, int** skipout);
void loc_reduce_spcell(loc_t* loc, dspcell* sp, int dim, int cont);
void loc_reduce_sp(loc_t* loc, dsp* sp, int dim, int cont);

void loc_add_focus(const dmat* opd, const loc_t* loc, real val);
dmat* loc2mat(loc_t* loc, int piston);
loc_t* pts2loc(pts_t* pts);
void locrot(loc_t* loc, const real theta);
real loc_angle(const loc_t* loc1, const loc_t* loc2);
void locstretch(loc_t* loc, const real theta, const real frac);
loc_t* locdup(loc_t* loc);
void locmean(real* xm, real* ym, const loc_t* loc);
dmat *parse_poly(const char *_ps);
loc_t* loctransform(const loc_t* loc, const char* ps);
loc_t* loctransform2(const loc_t* loc, const dmat* coeff);
void locshift(loc_t* loc, real sx, real sy);
void loc_nxny(long* nx, long* ny, const loc_t* loc);

void locresize(loc_t* loc, long nloc);
#define ptsresize(pts, nsa) if(pts) locresize(pts->loc, nsa)
void dembed_locstat(dmat** out, real alpha, loc_t* loc, real* oin, real beta, int reverse);
void cembed_locstat(cmat** out, real alpha, loc_t* loc, real* oin, real beta, int reverse);
void loc_dxdy(loc_t* loc);
loc_t* locreaddata(file_t* fp, header_t* header);
void loc_keywords(cell *loc);
lmat* loc_coord2ind(loc_t *aloc, dmat *dead);
#endif
