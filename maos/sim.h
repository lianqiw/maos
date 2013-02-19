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

#ifndef AOS_SIM_H
#define AOS_SIM_H
void sim(const PARMS_T *parms,  POWFS_T *powfs, APER_T *aper,  RECON_T *recon);
void wfsgrad(SIM_T *simu);
void maxapriori(double *g, dmat *ints, const PARMS_T *parms, 
		const POWFS_T *powfs, int iwfs, int isa, int noisy,
		double bkgrnd, double rne);
void wfslinearity(const PARMS_T *parms, POWFS_T *powfs, const int iwfs);
void perfevl_ievl(thread_t *info);
void perfevl(SIM_T *simu);
void prep_cachedm(SIM_T *simu);
void calc_cachedm(SIM_T *simu);
void filter_cl(SIM_T *simu);
void filter_ol(SIM_T *simu);
void filter(SIM_T *simu);
void update_dm(SIM_T *simu);
void hysterisis(HYST_T **hyst, dcell *dmreal, const dcell *dmcmd);
void wfsints(thread_t *thread_data);
void wfsgrad_iwfs(thread_t *info);
void addlow2dm(dcell **dmval, const SIM_T *simu, 
	       const dcell *low_val, double gain);
#define DEBUG_POWFSINTS 0
#endif
