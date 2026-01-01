/*
  Copyright 2009-2026 Lianqi Wang
  
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
/*dmat* karhunen_loeve=KL_vonkarman(dmat **cov, const loc_t*loc, double L0);*/
cell* read=readbin(const char* fn);
cell* readsock(int sock);
void  write=writebin(const cell*dc, const char* fn);
void  write_header=writebin_header(const cell*dc, const char* keywords, const char* fn);
void  writesock(const cell* dc, int sock);
double dtrapz(const dmat* x, const dmat* y);
void dsvd_pow(dmat*A, double power);
void dsvd_pow2(dmat*A, double power, real thres1, real thres2);
void dcellmm(cell **C, const cell*A, const cell*B, const char* trans, double alpha);
dmat* denc(dmat* psf, dmat* dvec, int type, int nthread);
void dcircle(dmat* A,double cx, double cy, double dx, double dy, double r, double val);
void drectangle(dmat *A, double cx, double cy, double rx, double ry, double theta, double val);
dsp* chol_factorize2(lmat** Cp, const dsp* A_in);
void addpath(const char* path);
void rmpath(const char* path);
double dgauss_width(dmat* A, double thres);
void dgauss_fit(double* mr, double* ma, double* mb, double* theta, dmat* A, double thres);
dmat* dpinv(dmat*A, cell *wt);
dmat* dpinv2(dmat*A, cell *wt, double thres1, double thres2);
dcell* dcellpinv(dcell *A, cell *wt);
dcell* dcellpinv2(dcell *A, cell *wt, double thres1, double thres2);
void dembed(dmat*out, dmat*in, double theta);
void dshift2center(dmat *A, real offx, real offy);
dmat *dcell2m(dcell*A);
dsp *dspaddI(dsp*A, double alpha);
