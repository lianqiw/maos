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

#ifndef AOS_LIB_CELLARR_H
#define AOS_LIB_CELLARR_H
#include "type.h"
/**
   \file cellarr.h
   cellarr is an object used to write arrays of dcell or ccell into file.
   Mainly used to output PSF into files.
*/
/*
   used to save array of dmat, cmat, ccell or dcell. mainly used to save
psfout. No user modifiable entries.  */
typedef struct cellarr{
    file_t *fp;     /**<pointer to file*/
    long cur;       /**<Current element*/
    long tot;       /**<Total number of elements*/
}cellarr;
cellarr* cellarr_init(long nx, long ny, const char*format,...) CHECK_ARG(3);
void cellarr_push(cellarr *ca, int i, const void *A);
void cellarr_dcell(cellarr *ca, int i, const dcell *A);
void cellarr_scell(cellarr *ca, int i, const scell *A);
void cellarr_ccell(cellarr *ca, int i, const ccell *A);
void cellarr_zcell(cellarr *ca, int i, const zcell *A);
void cellarr_dmat(cellarr *ca, int i, const dmat *A);
void cellarr_smat(cellarr *ca, int i, const smat *A);
void cellarr_cmat(cellarr *ca, int i, const cmat *A);
void cellarr_zmat(cellarr *ca, int i, const zmat *A);
void cellarr_close(cellarr *ca);
void cellarr_close_n(cellarr **ca, int nc);
#endif
