/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_CELLARR_H
#define AOS_CELLARR_H
#include "bin.h"
#include "type.h"
/**
   used to save array of dmat, cmat, ccell or dcell. mainly used to save
psfout. No user modifiable entries.  */
typedef struct cellarr{
    const char *fn; /**<file name*/
    file_t *fp;     /**<pointer to file*/
    long cur;       /**<Current element*/
    long tot;       /**<Total number of elements*/
}cellarr;
cellarr* cellarr_init(int tot,const char*format,...) CHECK_ARG(2);
void cellarr_dcell(cellarr *ca, const dcell *A);
void cellarr_ccell(cellarr *ca, const ccell *A);
void cellarr_dmat(cellarr *ca, const dmat *A);
void cellarr_cmat(cellarr *ca, const cmat *A);
void cellarr_close(cellarr *ca);
void cellarr_close_n(cellarr **ca, int nc);
#endif
