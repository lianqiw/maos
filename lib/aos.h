/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_AOS_H
#define AOS_LIB_AOS_H

#include "../math/mathdef.h"
#include "readcfg.h"
#include "draw.h"
#include "accphi.h"
#include "mkz.h"
#include "mkg.h"
#include "mkh.h"
#include "mkw.h"
#include "turbulence.h"
#include "proj.h"
#include "laplacian.h"
#include "pcg.h"
#include "muv.h"
#include "genotf.h"
#include "fractal.h"
#include "stfun.h"
#include "servo.h"
#include "hyst.h"
#include "slaving.h"
#include "mkdtf.h"
#include "psd.h"
#include "kalman.h"
#include "cn2est.h"
#include "locfft.h"
#include "zernike.h"
#include "misc.h"
#define SEC2RAD 4.848136811095360e-06 //arcsec in unit of radian
#define RAD2SEC 206264.8062470964 //radian in unit of arcsec
#define RAD2MAS 2.062648062470964e+08 //radian to milli-arcsecond
#endif
