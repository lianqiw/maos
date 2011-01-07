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

#ifndef __DRAW_H
#define __DRAW_H
#include "misc.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
extern int DRAW_ID;//number to use for drawdaemon, usually PID.
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
void plot_points(char *fig, long ngroup, loc_t **loc, dcell *cell,
		 const int32_t *style, const double *limit, int ncir, double(*pcir)[4], 
		 char **legend, const char *title, const char *xlabel, const char *ylabel,
		 const char *format,...) CHECK_ARG(13);
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

//set to one to disable drawing.
extern int disable_draw;
#define DRAW_GRAY   0x0
#define DRAW_COLOR  0x1
enum{
    FIFO_START=0, //Mark the starting of data stream.
    FIFO_DATA,
    FIFO_SHM,
    FIFO_POINTS,
    FIFO_STYLE,
    FIFO_CIRCLE,
    FIFO_LIMIT,
    FIFO_FIG,
    FIFO_NAME,
    FIFO_TITLE,
    FIFO_XLABEL,
    FIFO_YLABEL,
    FIFO_ZLIM,
    FIFO_LEGEND,//legend
    FIFO_END=100
};

#endif
