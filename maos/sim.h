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
void perfevl_ievl(thread_t *info);
void perfevl(SIM_T *simu);
void prep_cachedm(SIM_T *simu);
void calc_cachedm(SIM_T *simu);
void filter_dm(SIM_T *simu);
void update_dm(SIM_T *simu);
void wfsgrad(SIM_T *simu);
void wfsints(thread_t *thread_data);
void wfsgrad_iwfs(thread_t *info);
void wfsgrad_post(thread_t *info);
void addlow2dm(dcell **dmval, const SIM_T *simu, 
	       const dcell *low_val, double gain);

#endif
