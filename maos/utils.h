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

#ifndef __AOS_UTILS_H
#define __AOS_UTILS_H
#include "maos.h"
#include <signal.h>
#if USE_STATIC
extern char _binary____config_tar_gz_start;
extern char _binary____config_tar_gz_end;
#endif
extern SIM_T *cursimu;//Used for the signal handler only.
void addnoise(dmat *A, rand_t* rstat, 
	      const double bkgrnd, const double pcalib, 
	      const double *bkgrnd2, const double pcalib2,
	      const double rne);

void plotloc(char *fig, const PARMS_T *parms, 
	     loc_t *loc, double ht, char *format,...);
void rename_file(int sig);
void maos_signal_handler(int sig);
void wait_cpu(int nthread);
ARG_T* parse_args(int argc, char **argv);
cmat *strehlcomp(const dmat *iopdevl, const double *amp, const double wvl);
ccell *psfcomp(const dmat *iopdevl, const double *restrict amp,
	       int **embeds, const int *nembeds, const int *psfsize,
	       const int nwvl, const double *wvl);
void psfcomp_iwvl(thread_t *tdata);
void embed_in(double *out, const double *in, long nin, long *embed);
void embed_out(const double *out, double *in, long nin, long *embed);
void embedc_in(dcomplex *out, const double *in, long nin, long *embed);
void embedc_out(const dcomplex *out, double *in, long nin, long *embed);
int lock_seeds(PARMS_T *parms);
double calc_aniso(double r0, int nht, double *ht, double *wt);
double calc_aniso2(double r0, int nht, double *ht, double *wt, double hc1, double hc2);

void shift_inte(int nap, double *ap, dcell **inte);
#endif
