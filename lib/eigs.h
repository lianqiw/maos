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

#ifndef AOS_EIGS_H
#define AOS_EIGS_H
#include "dsp.h"
#include "dmat.h"
#include "cmat.h"
extern void dsaupd_(int *ido, char *bmat, int *n, char *which,
		    int *nev, double *tol, double *resid, int *ncv,
		    double *v, int *ldv, int *iparam, int *ipntr,
		    double *workd, double *workl, int *lworkl,
		    int *info);
extern void dseupd_(int *rvec, char *All, int *select0, double *d,
		    double *z, int *ldz, double *sigma,
		    char *bmat, int *n, char *which, int *nev,
		    double *tol, double *resid, int *ncv, double *v,
		    int *ldv, int *iparam, int *ipntr, double *workd,
		    double *workl, int *lworkl, int *ierr);
typedef enum ATYPE{
    A_SPARSE=0,
    A_SPCELL,
}ATYPE;
void eigs(dmat **eig, dmat **eigv, const void *A, ATYPE Atype, 
	  const int neig, char *eigtype);
double spmaxeig(const dsp *A);
double spmineig(const dsp *A);
double spcellmaxeig(const spcell *A);
double spcellmineig(const spcell *A);
void spcelltikcr(spcell *A, double thres);
void sptikcr(dsp *A, double thres);
#endif
