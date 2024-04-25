/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
int plot_empty(int sock_draw,const char *fig,const char *fn);
typedef struct {
  union{
    const dmat* image;//2d image. can be real or complex. set ctype for complex.
    const cmat* cimage;
  };
  int ngroup;///number of lines
  int maxlen;///limit maximum nx
  loc_t** loc;
  const dcell* dc;
  const int32_t* style;
  const real zlim[2];//z limit
  const real* limit;//x and y limit
  const int zlog;//z log.
  const int ctype; //how to convert cimage to image. 0: abs, 1: phase. 2: real, 3: imaginary.
  const char* xylog;
  const dmat* cir;
  const char* const* const legend;
} plot_opts;

int draw(const char* fig, plot_opts opts, const char* title, const char* xlabel, const char* ylabel,
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
int drawgrad(const char* fig, loc_t* saloc, const dmat* grad, int grad2opd, int trs, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(10);
int drawopdamp(const char *fig, loc_t *loc, const dmat *opd, const dmat *amp, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(9);
int drawints(const char *fig, const loc_t *saloc, const dcell *ints, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(8);
#endif
