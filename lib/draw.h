/*
  Copyright 2009-2026 Lianqi Wang
  
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

#include "../math/mathdef.h"
/**
   \file draw.h
   Contains functions for data visualization.
*/
extern int draw_id;/*number to use for drawdaemon, usually PID. */
extern int draw_direct;
extern int draw_disabled;
extern int draw_single;
void draw_helper(void);
int draw_add(int fd);
void draw_final(int reuse);
int draw_current(const char* fig, const char* fn);
int draw_current_format(const char *fig, const char *format, ...) CHECK_ARG(2);

typedef struct {
  //For 2-d images
  union{
    const dmat* image;//2-d image. real valued. 
    const cmat* cimage;//2-d image. complex value. ctype determines how to convert to real value
  };
  const real zlim[2];//for 2-d image. value limit min, max 
  const int zlog;//for 2-d image. Use log10 of value.
  const int ctype; //how to convert cimage to image. 0: abs, 1: phase. 2: real, 3: imaginary.
  //For lines
  dmat** dd;//array of dmat
  loc_t** loc;//2-d grid points to plot.
  const dcell* dc;//dcell. takes priority over dd
  int ngroup;///number of lines. 0: automatic
  int maxlen;///maximum nx to plot
  const real dx;//spacing of x if only y is provided.
  const dmat *dxs;//spacing of x if only y is provided. one value per line. (if dx is not specified)
  const char* xylog;//for lines, take log along x (xylog[0]=='y') and/or y (xylog[1]=='n')
  const int32_t* style;//line style
  const char* const* const legend;//legends for lines
  //common for 2-d image or lines
  const dmat* cir;//circles to draw
  const real* limit;//manual xmin,xmax,ymin,ymax
  int always;//1: always draw (do not check for active panel)
} plot_opts;

int draw(const char* fig, const plot_opts opts, const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(6);

int drawmap(const char *fig, const map_t *map, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(7);
int drawloc(const char *fig, loc_t *loc, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(7);
int drawopd(const char *fig, loc_t *loc, const dmat *opd, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(8);
int drawgrad(const char* fig, loc_t* saloc, const dmat *saa, const dmat* grad, int grad2opd, int trs, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(11);
int drawopdamp(const char *fig, loc_t *loc, const dmat *opd, const dmat *amp, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(9);
int drawints(const char *fig, const loc_t *saloc, const dcell *ints, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(8);
#endif
