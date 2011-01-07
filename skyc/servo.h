/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef SKYC_SERVO_H
#define SKYC_SERVO_H
#include "skyc.h"
#include "types.h"
dmat* servo_typeII_optim(double *ress, double *resn, const dmat *psdin,
			 double fs, double lgsdt, double sigman);
double servo_typeII_residual(const dmat *gain, const dmat *psdin, double fs, double lgsdt);
void servo_typeII_filter(SERVO_S *st, dmat *merr, double dtngs, const dmat *gain);
void servo_typeI_filter(SERVO_S *st, dmat *merr, double gain);
dmat *psd2temp(dmat *psdin, double dt, double N, rand_t* rstat);
dmat* servo_typeII_test(dmat *mideal, dmat *gain, double dtngs, int dtrat);
void servo_free(SERVO_S *st);
cmat *servo_typeII_Hol(const dmat *gain, double fs, double lgsdt);
double psd_intelog(double *nu, double *psd, long n);
double psd_intelog2(dmat *psdin);
#endif
