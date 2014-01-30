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
#ifndef SKYC_SETUP_STARS_H
#define SKYC_SETUP_STARS_H
STAR_S *setup_star(int *nstarout,SIM_S* simu, dmat *stars, int seed);
void free_istar(STAR_S *star, const PARMS_S *parms);
void free_star(STAR_S *star, int nstar, const PARMS_S *parms);
long setup_star_read_wvf(STAR_S *star, int nstar, const PARMS_S *parms, int seed);
long setup_star_read_ztilt(STAR_S *star, int nstar, const PARMS_S *parms, int seed);
void free_pistat(PISTAT_S *pistat, int npistat, const PARMS_S *parms);
#endif
