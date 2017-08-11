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
cell *genotfmex(loc_t *loc, const dmat *amp, const dmat *opdbias, const dmat *area, double thres, double wvl, double dtheta, const dmat *cov, double r0, double l0, long ncompx, long ncompy, long nsa, long pttr);
dmat *m3projmex(dmat *mapin_0, char *header, loc_t *locout, double thetax, double thetay, double hs);
dmat *mkcirmap(long nx, long ny, double cx, double cy, double r);
dmat *mtch2(dmat **nea, const dmat *i0, const dmat *gx, const dmat *gy, int cr);
dmat *sdepsd(const dmat *ff, const dmat *coeff);
