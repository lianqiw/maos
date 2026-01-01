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
#ifndef SKYC_SETUP_STARS_H
#define SKYC_SETUP_STARS_H
star_s *setup_star(int *nstarout,sim_s* simu, dmat *stars, int seed);
void free_istar(star_s *star, const parms_s *parms);
void free_star(star_s *star, int nstar, const parms_s *parms);
long setup_star_read_wvf(star_s *star, int nstar, const parms_s *parms, int seed);
long setup_star_read_ztilt(star_s *star, int nstar, const parms_s *parms, int seed);
void free_pistat(pistat_s *pistat, int npistat, const parms_s *parms);
#endif
