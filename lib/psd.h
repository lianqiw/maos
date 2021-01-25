/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_LIB_PSD_H
#define AOS_LIB_PSD_H
/**
   \file psd.h
   Compute the PSD from a sequence.
*/
dmat *psd1d(const dmat *v, long nseg );
dmat *psd1dt(const dmat *v, long nseg, real dt);
dmat *psdinterp1(const dmat *psdin, const dmat *fnew, int uselog);
dmat *psd_vibid(const dmat *psdin);
dmat *psdt2s(const dmat *psdt, real vmean);
dmat *psds2t(const dmat *psdt, real vmean);
real psd_inte(const real *nu, const real *psd, long n);
real psd_inte2(const dmat *psdin);
dmat* psd2time(const dmat *psdin, rand_t *rstat, real dt, int nstep);
dmat* add_psd(const dmat *psd1, const dmat *psd2, real scale2);
void add_psd2(dmat **out, const dmat *in, real scale);
void psd_sum(dmat *psd, real scale);
#endif
