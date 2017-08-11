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
/*
  This file is parsed by aolib.py. The content here does not necessarily adhere to C standard.
 */
dmat *zernike_Rnm(const dmat *locr, int ir, int im);
dmat* zernike(const loc_t *loc, double D, int rmin, int rmax, int onlyr);
dmat *zernike_cov_kolmogorov(int nr);
dmat *cov_vonkarman(const loc_t *loc, const dmat *modz, double L0);
dmat *cov_diagnolize(const dmat *mod, const dmat *cov);
dmat *karhunen_loeve=KL_vonkarman(const loc_t *loc, int nmod, double L0);
dmat* turbcov(dmat *r, double rmax, double r0, double L0);

cell* readbin(const char *fn);
void  writebin(const cell* dc, const char* fn);

dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, int vibid);
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
void svdpow=dsvd_pow(dmat *A, double power, double thres);
cell* cellmm=dcellmm2(const cell *A, const cell *B, const char*trans);
dsp *mkg(loc_t* xloc, loc_t *ploc, dmat *amp, loc_t *saloc, double scale, double dispx, double dispy, int do_partial);
dmat *sho_filter(const dmat *xi, double dt, double f0, double zeta);
cell *genotf=genotfmex(loc_t *loc, const dmat *amp, const dmat *opdbias, const dmat *area, double thres, double wvl, double dtheta, const dmat *cov, double r0, double l0, long ncompx, long ncompy, long nsa, long pttr);
dmat *calcenc=denc(dmat *psf, dmat *dvec, int type, int nthread);
dmat *m3proj=m3projmex(dmat *mapin_0, char *header, loc_t *locout, double thetax, double thetay, double hs);
dspcell *act_extrap(loccell *aloc, const dcell *actcplc, double thres);
dmat *mkcirmap(long nx, long ny, double cx, double cy, double r);
dmat *psd1d(dmat* data, long nseg);
dmat *psd1dt(dmat* data, long nseg, double dt);
dmat* mk2dcov(loc_t *loc, const dmat *amp, double ampthres, const dmat *cov, int norm);
void mkw=mkw_circular(loc_t *loc,double cx, double cy, double r,dsp **W0,dmat **W1);
dcell* servo_optim(const dmat *psdin,  double dt, long dtrat, double pmargin, const dmat* sigma2n, int servo_type);
dmat* servo_test(dmat *input, double dt, int dtrat, dmat *sigma2n, dmat *gain);
double servo_residual(double *noise_amp, const dmat *psdin, double dt, long dtrat, const dmat *gain, int servo_type);
kalman_t* sde_kalman(const dmat *coeff, double dthi, const dmat* dtrat, const dcell *Gwfs, const dcell *Rwfs, const dmat *Proj);
dmat *kalman_test(kalman_t *kalman, dmat *input);
dspcell *slaving(loccell *aloc, const dcell *actcpl, dcell *NW, const lcell *actstuck, const lcell *actfloat, double thres, double scl);
dmat *mtch(dmat **nea, const dmat *i0, const dmat *gx, const dmat *gy, const dmat *bkbrnd2, const dmat *bkgrnd2c, double bkgrnd, double bkgrndc, double rne, double pixthetax, double pixthetay, double pixrot, int radgx, int cr);
dmat *mtch2(dmat **nea, const dmat *i0, const dmat *gx, const dmat *gy, int cr);
dmat *sdepsd(const dmat *ff, const dmat *coeff);
