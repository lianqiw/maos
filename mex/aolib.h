/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

dmat *zernike_Rnm(const dmat *locr, int ir, int im);
dmat* zernike(const loc_t *loc, double D, int rmin, int rmax, int onlyr);
dmat *zernike_cov_kolmogorov(int nr);
dmat *cov_vonkarman(const loc_t *loc, const dmat *modz, double L0);
dmat *cov_diagnolize(const dmat *mod, const dmat *cov);
dmat *KL_vonkarman(const loc_t *loc, int nmod, double L0);
dmat* turbcov(dmat *r, double rmax, double r0, double L0);

cell* readbin(const char *fn);
void  writebin(const cell* dc, const char* fn);

dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit);
double dtrapz(const dmat *x, const dmat *y);
dmat *psdinterp1(const dmat *psdin, const dmat *fnew, int uselog);
dmat *psd_vibid(const dmat *psdin);
dmat* psd2time(const dmat *psdin, rand_t *seed, double dt, int nstepin);
dmat *psdt2s(const dmat *psdt, double vmean);
dmat *psds2t(const dmat *psdt, double vmean);
dmat *psd1d(const dmat *v, long nseg);
dmat *psd1dt(const dmat *v, long nseg, double dt);
dsp * mkh(loc_t *locin, loc_t *locout, double displacex, double displacey,double scale);
dsp * mkh_cubic(loc_t *locin, loc_t *locout, double displacex, double displacey,double scale, double cubic_iac);
void dsvd_pow(dmat *A, double power, double thres);
cell* dcellmm2(const cell *A, const cell *B, const char*trans);
dsp *mkg(loc_t* xloc, loc_t *ploc, dmat *amp, loc_t *saloc, double scale, double dispx, double dispy, int do_partial);
dmat *sho_filter(const dmat *xi, double dt, double f0, double zeta);
