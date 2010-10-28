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

#ifndef __DRAW_H
#define __DRAW_H
#include "misc.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
#if USE_DAEMON==0
#define plot_psf(A...)
#define imagesc(A...)
#define imagesc_long(A...)
#define imagesc_cmp(A...)
#define imagesc_cmp_ri(A...)
#define imagesc_cmp_ap(A...)
#define plot_pts(A...)
#define ddraw(A...)
#define cdraw(A...)
#define cdrawabs(A...)
#define cdrawri(A...)
#define cdrawap(A...)
#define drawmap(A...)
#define drawloc(A...)
#define drawopd(A...)
#define drawopdamp(A...)
#else
void imagesc(char *fig, long nx, long ny, const double *limit,
	     const double *p, int color,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...) CHECK_ARG(10);
void imagesc_long(char *fig, long nx, long ny,const double *limit,
		  const long *p, int color, 
		  const char *title, const char *xlabel, const char *ylabel,
		  const char *format,...) CHECK_ARG(10);
void imagesc_cmp(char *fig, long nx, long ny, const double *limit,
		 const dcomplex *p, int color, 
		 const char *title, const char *xlabel, const char *ylabel,
		 const char *format,...) CHECK_ARG(10);
void imagesc_cmp_ri(char *fig, long nx, long ny, const double *limit,
		    const dcomplex *p, int color,
		    const char *title, const char *xlabel, const char *ylabel, 
		    const char *format,...) CHECK_ARG(10);
void imagesc_cmp_ap(char *fig, long nx, long ny, const double *limit,
		    const dcomplex *p, int color, 
		    const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...) CHECK_ARG(10);
void plot_coord(char *fig, long npts, const double *ptsx, const double *ptsy, 
		const long *style, const double *limit,
		int ncir, const void *pcir, 
		const char *title, const char *xlabel, const char *ylabel,
		const char *format,...) CHECK_ARG(12);
void ddraw(char *fig, const dmat *A, 
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...) CHECK_ARG(6);
void cdraw(char *fig, const cmat *A, 
	   const char *title, const char *xlabel, const char *ylabel,
	   const char*format,...) CHECK_ARG(6);
void cdrawabs(char *fig, const cmat *A, 
	      const char *title, const char *xlabel, const char *ylabel,
	      const char *format,...) CHECK_ARG(6);
void cdrawri(char *fig, const cmat *A, 
	     const char *title, const char *xlabel, const char *ylabel,
	     const char*format,...) CHECK_ARG(6);

void drawmap(char *fig, const MAP_T *map, 
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...) CHECK_ARG(6);
void drawloc(char *fig, LOC_T *loc, 
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...) CHECK_ARG(6);
void drawopd(char *fig, LOC_T *loc, const double *opd, 
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...) CHECK_ARG(7);
void drawopdamp(char *fig, LOC_T *loc, const double *opd, const double *amp, 
		const char *title, const char *xlabel, const char *ylabel,
		const char* format,...) CHECK_ARG(8);

//set to one to disable drawing.
extern int disable_draw;
#define DRAW_GRAY   0x0
#define DRAW_COLOR  0x1
#endif
double* bin(long nx, long ny, double *p, int factor);
void dbl2gray(long nx, long ny, const double *restrict p, 
	      unsigned char *pc, double *info);
void cmp2gray(long nx, long ny, const dcomplex *restrict p, 
	      unsigned char *pc, double *info);
void long2gray(long nx, long ny, const long *p, 
	       unsigned char *pc, double *info);
void dbl2rgb(long nx, long ny, const double *p, 
	     int *pi, double *info);
void cmp2rgb(long nx, long ny, const dcomplex *restrict p, 
	     int *pi, double *info);
void long2rgb(long nx, long ny, const long *p, 
	      int *pi, double *info);

#endif
