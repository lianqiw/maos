/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_POWFS_H
#define AOS_POWFS_H
#include "common.h"
POWFS_T * setup_powfs_init(const PARMS_T *parms, APER_T *aper);
void setup_powfs_phy(const PARMS_T *parms,  POWFS_T *powfs);
void setup_powfs_calib(const PARMS_T *parms, POWFS_T *powfs, loccell *aloc, dcell *dm_ncpa);
void free_powfs_unused(const PARMS_T *parms, POWFS_T *powfs);
void free_powfs(const PARMS_T *parms, POWFS_T *powfs);
void test_powfs(const PARMS_T *parms, POWFS_T *powfs);
void setup_powfs_etf(POWFS_T *powfs, const PARMS_T *parms, 
		     int ipowfs, int mode, int istep);
void wfspupmask(const PARMS_T *parms, loc_t *loc, dmat *amp, int iwfs);
void pywfs_setup(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs);
void pywfs_free(PYWFS_T *pywfs);
void pywfs_grad(dmat **pgrad, const PYWFS_T *pywfs, const dmat *ints);
void pywfs_fft(dmat **ints, const PYWFS_T *pywfs, const dmat *opd);
dsp* pywfs_mkg(const PYWFS_T *pywfs, const loc_t* ploc, int cubic, double iac);
#endif
