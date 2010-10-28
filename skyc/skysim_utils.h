/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
void ngsmod2wvf(cmat *wvf, double wvl, const dmat *modm,
		const LOC_T *cloc, double fpc, double thetax, double thetay, 
		const PARMS_S *parms);
dcell* skysim_ztilt(dmat *mideal, ASTER_S *aster, const PARMS_S *parms);
dmat* skysim_phy(dmat **mres,SIM_S *simu, ASTER_S *aster, POWFS_S *powfs, const PARMS_S *parms, 
		 int idtrat, int noisy, int demotettf);
void skysim_save(SIM_S *simu, ASTER_S *aster, double *ipres, int selaster, int seldtrat, int isky);
#endif
