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

double dtrapz(const dmat *x, const dmat *y);
dmat *psdinterp1(const dmat *psdin, const dmat *fnew, int uselog);
dmat* zernike(loc_t *loc, double D, int rmin, int rmax, int onlyr);
dmat *zernike_cov_kolmogorov(int nr);
dmat *diag_mod_cov(dmat *mz, dmat *cov);
dmat *KL_kolmogorov(loc_t *loc, double D, int nr, int nr2);
dmat *psd_vibid(const dmat *psdin);
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit);
cell* readbin(const char *fn);
void  writebin(const cell* dc, const char* fn);
dmat* psd2time(const dmat *psdin, rand_t *seed, double dt, int nstepin);
dmat *psdt2s(const dmat *psdt, double vmean);
dmat *psds2t(const dmat *psdt, double vmean);
