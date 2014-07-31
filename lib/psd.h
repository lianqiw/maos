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
#ifndef AOS_LIB_PSD_H
#define AOS_LIB_PSD_H
/**
   \file psd.h
   Compute the PSD from a sequence.
*/
dmat *psd1d(dmat *v, long lseg );
dmat *psd1dt(dmat *v, long lseg, double dt);
dmat *psdinterp1(const dmat *psdin, const dmat *fnew);
dmat *psd_vibid(const dmat *psdin);
#endif
