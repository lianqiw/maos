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
/*
  This file is parsed by lib2mex.py. The content here does not necessarily adhere to C standard.
 */
dmat* zernike_Rnm(const dmat* locr, int ir, int im);
dmat* zernike(const loc_t* loc, double D, int rmin, int rmax, int onlyr);
dmat* zernike_cov_kolmogorov(int nr);
dmat* cov_vonkarman(const loc_t* loc, const dmat* modz, double L0);
dmat* cov_diagnolize(const dmat* mod, const dmat* cov);
dmat* karhunen_loeve=KL_vonkarman(const loc_t*loc, int nmod, double L0);
dmat* turbcov(dmat* r, double rmax, double r0, double L0);
cell* read=read_by_id(0, -1, const char* fn);
cell* readsock(int sock);
void  write=write_by_id(const cell*dc, 0, const char* fn);
void  write_header=writebin_header(const cell*dc, const char* header, const char* fn);
void  writesock(const cell* dc, int sock);
dmat* sde_fit(const dmat* psdin, const dmat* coeff0, double tmax_fit, int vibid);
double dtrapz(const dmat* x, const dmat* y);
dmat* psdinterp1(const dmat* psdin, const dmat* fnew, int uselog);
dmat* psd_vibid(const dmat* psdin);
dmat* psd2time(const dmat* psdin, rand_t* seed, double dt, int nstepin);
dmat* psdt2s(const dmat* psdt, double vmean);
dmat* psds2t(const dmat* psdt, double vmean);
dmat* psd1d(const dmat* v, long nseg);
dmat* psd1dt(const dmat* v, long nseg, double dt);
dsp* mkh(loc_t* locin, loc_t* locout, double displacex, double displacey, double scale);
dsp* mkh_cubic(loc_t* locin, loc_t* locout, double displacex, double displacey, double scale, double cubic_iac);
void svdpow=dsvd_pow(dmat*A, double power, double thres, double tikcr);
void cellmm=dcellmm_any(cell **C, const cell*A, const cell*B, const char* trans, real alpha);
dsp* mkg(loc_t* xloc, loc_t* ploc, dmat* amp, loc_t* saloc, double scale, double dispx, double dispy, int do_partial);
dmat* sho_filter(const dmat* xi, double dt, double f0, double zeta);
void genotf(ccell**potf,loc_t*loc, const dmat*amp, const dmat*opdbias, const dmat*area, double thres, double wvl, const dmat*cov, double r0, double l0, long npsfx, long npsfy, long nsa, long pttr);
dmat* calcenc=denc(dmat*psf, dmat*dvec, int type, int nthread);
dmat* m3proj=m3proj2(dmat*mapin_0, char* header, loc_t*locout, double thetax, double thetay, double hs);
dspcell* act_extrap(loccell* aloc, const dcell* actcplc, double thres);
dmat* mkcirmap(long nx, long ny, double cx, double cy, double r);
dmat* psd1d(dmat* data, long nseg);
dmat* psd1dt(dmat* data, long nseg, double dt);
dmat* mk2dcov(loc_t* loc, const dmat* amp, double ampthres, const dmat* cov, int norm);
void mkw=mkw_circular(loc_t*loc, double cx, double cy, double r, dsp**W0, dmat**W1);
dmat* servo_rej2ol(const dmat* psdcl, double dt, long dtrat, double al, double gain, double sigma2n);
dcell* servo_optim(const dmat* psdin, double dt, long dtrat, long al, double pmargin, const dmat* sigma2n, int servo_type);
dmat* servo_test(dmat* input, double dt, int dtrat, dmat* sigma2n, dmat* gain);
double servo_residual(double* noise_amp, const dmat* psdin, double dt, long dtrat, long al, const dmat* gain, int servo_type);
kalman_t* sde_kalman(const dmat* coeff, double dthi, const lmat* dtrat, const dcell* Gwfs, const dcell* Rwfs, const dmat* Proj);
dmat* kalman_test(kalman_t* kalman, dmat* input);
dspcell* slaving(loccell* aloc, const dcell* actcpl, const lcell* actstuck, const lcell* actfloat, double thres, double scale, int mode);
void mtch(dmat** mtche, dmat** nea, const dmat* i0, const dmat* gx, const dmat* gy, const dmat* qe, const dmat* bkbrnd2, const dmat* bkgrnd2c, double bkgrnd, double bkgrndc, double rne, double pixthetax, double pixthetay, double pixrot, int radgx, int cr);
void mtch2(dmat** mtche, dmat** nea, const dmat* i0, const dmat* gx, const dmat* gy, int cr);
dmat* sde_psd2(const dmat* ff, const dmat* coeff);
cn2est_t* cn2est=cn2est_all(const dmat*wfspair, dmat*wfstheta, const loc_t*saloc,
	const dmat*saa, const double saat, const dmat*hs,
	const dmat*htrecon, int keepht, double l0, dcell*grad);

dtf_t* mkdtf(const dmat* wvls, double dxsa, double embfac, long notfx, long notfy, long pixpsax, long pixpsay, double pixthetax, double pixthetay, const dmat* pixoffx, const dmat* pixoffy, double pixblur, const dcell* pixrot);
dmat* add_psd(const dmat* psd1, const dmat* psd2, double scale2);
dmat* derive_by_fft(const dmat* i0, double theta);
void cure(dmat** phi, const dmat* gx, const dmat* gy, double dx);
void cure1d(dmat** phi, const dmat* gx, const dmat* gy, double dx);
void cure_loc(dmat** phix, const dmat* grad, const loc_t* gloc);
dcell* loc_embed=loc_embed2(loc_t*loc, dmat*arr);
double loc_angle(const loc_t* loc1, const loc_t* loc2);
void locrot(loc_t* loc, const double theta);
loc_t* mksqloc(long nx, long ny, double dx, double dy, double ox, double oy);
dmat* loc_calib(const dsp* GA, const loc_t* aloc, const loc_t* saloc,
	double dispx, double dispy, double scale, int maxorder);
dmat* poly2fit(const dmat* in, const dmat* out, int maxorder);
loc_t* loctransform2(const loc_t* loc, const dmat* coeff);
dsp* chol_factorize2(long** Cp, const dsp* A_in);
dmat* hyst_test(double hysteresis, double hyst_alpha, double hyst_stroke, const dmat* dmcmd);
void addpath(const char* path);
void rmpath(const char* path);
void fresnel_prop(cmat** out, double* dxout, const cmat* in, double dxin, double wvl, double z, double scale, int method);
double dgauss_width(dmat* A, double thres);
void dgauss_fit(double* mr, double* ma, double* mb, double* theta, dmat* A, double thres);
dmat* dpinv2(dmat*A, cell *wt, real svdthres, real tikcr);
dcell* dcellpinv2(dcell *A, cell *wt, real svdthres, real tikcr);
real loc_remove_focus_grad(dmat *grad, const loc_t *saloc, real factor);