/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef SKYC_SIM_H
#define SKYC_SIM_H
real calc_rms(const dmat *mideal, const dmat *mcc, int istep0);
void ngsmod2wvf(cmat *wvf, real wvl, const dmat *modm,
		const powfs_s *powfs, int isa, real thetax, real thetay, 
		const parms_s *parms);
dmat* physim(dmat **mres, const dmat *mideal, const dmat *mideal_oa, real ngsol,
		 aster_s *aster, const powfs_s *powfs, const parms_s *parms, 
		 int idtrat, int noisy, int phystart);
void skysim_save(const sim_s *simu, const aster_s *aster, const real *ipres, int selaster, int seldtrat, int isky);
#endif
