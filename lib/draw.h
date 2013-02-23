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

#ifndef __DRAW_H
#define __DRAW_H

#include "dmat.h"
#include "cmat.h"
#include "loc.h"
extern int DRAW_ID;/*number to use for drawdaemon, usually PID. */
extern int DRAW_DIRECT;
void draw_helper(void);
int draw_add(int fd);
void imagesc(char *fig, long nx, long ny, const double *limit, const double *zlimit,
	     const void *p, const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...) CHECK_ARG(10);
void imagesc_cmp_ri(char *fig, long nx, long ny, const double *limit, const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel, 
		    const char *format,...) CHECK_ARG(10);
void imagesc_cmp_ap(char *fig, long nx, long ny, const double *limit, const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...) CHECK_ARG(10);
void imagesc_cmp_abs(char *fig, long nx, long ny, const double *limit,const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		     const char *format,...) CHECK_ARG(10);
void plot_points(char *fig, long ngroup, loc_t **loc, dcell *dc,
		 const int32_t *style, const double *limit, const char *xylog, int ncir, double(*pcir)[4], 
		 char **legend, const char *title, const char *xlabel, const char *ylabel,
		 const char *format,...) CHECK_ARG(14);
void ddraw(char *fig, const dmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...) CHECK_ARG(8);
void cdraw(char *fig, const cmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char*format,...) CHECK_ARG(8);
void cdrawabs(char *fig, const cmat *A, double *xylim, double *zlim,
	      const char *title, const char *xlabel, const char *ylabel,
	      const char *format,...) CHECK_ARG(8);
void cdrawri(char *fig, const cmat *A, double *xylim, double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char*format,...) CHECK_ARG(8);

void drawmap(char *fig, const map_t *map,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...) CHECK_ARG(7);
void drawloc(char *fig, loc_t *loc,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...) CHECK_ARG(7);
void drawopd(char *fig, loc_t *loc, const double *opd,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel, 
	     const char* format,...) CHECK_ARG(8);
void drawopdamp(char *fig, loc_t *loc, const double *opd, const double *amp, double *zlim,
		const char *title, const char *xlabel, const char *ylabel,
		const char* format,...) CHECK_ARG(9);
enum{
    DRAW_START=0, /*Mark the starting of data stream. */
    DRAW_DATA,
    DRAW_SHM,
    DRAW_POINTS,
    DRAW_STYLE,
    DRAW_CIRCLE,
    DRAW_LIMIT,
    DRAW_FIG,
    DRAW_NAME,
    DRAW_TITLE,
    DRAW_XLABEL,
    DRAW_YLABEL,
    DRAW_ZLIM,
    DRAW_LEGEND,/*legend */
    DRAW_XYLOG,/*logscale long x, y*/
    DRAW_FIGFN,
    DRAW_PAUSE,
    DRAW_RESUME,
    DRAW_END=100
};

#endif
