/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef SKYC_NAFOCUS_H
#define SKYC_NAFOCUS_H
/**
   Compute sodium power spectrum density. alpha, beta are the parameters of the
   sodium power spectrum obtained by fitting measurement data at UBC. */
INLINE double nafocus_NaPSD(double nu, double alpha, double beta2){
    return beta2*pow(nu,alpha);/*we don't divide 2pi */
}
double nafocus_residual(double fs, double tau, double f_zc, double zeta, 
			double D, double hs, double alpha, double beta);
dmat* nafocus_time(double alpha, double beta, double dt, long nstep, rand_t *rstat);
#endif
