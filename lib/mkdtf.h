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
#ifndef AOS_LIB_DTF_H
#define AOS_LIB_DTF_H
#include "type.h"
#include "loc.h"
/**
   \file mkdtf.h
   Routine to generate detector transfer function.
*/
void mkdtf(ccell **pnominal, /**<[out] to be multiplied to the OTF*/
	   spcell **psi,     /**<[out] to be applied after IFFT of the final OTF*/
	   int ncompx,       /**<[in] size of OTF FFT*/
	   int ncompy,       /**<[in] size of OTF FFT*/
	   double dtheta,    /**<[in] sampling of PSF*/
	   int pixpsax,      /**<[in] number of pixels along x dimension*/
	   int pixpsay,      /**<[in] number of pixels along y dimension*/
	   double pixthetax, /**<[in] size of pixel along x dimension*/
	   double pixthetay, /**<[in] size of pixel along y dimension*/
	   double pixoffx,   /**<[in] offset of the image from the center of detector.*/
	   double pixoffy,   /**<[in] offset of the image from the center of detector*/
	   double blurx,     /**<[in] blurring as a percentage of pixel*/
	   double blury,     /**<[in] blurring as a percentage of pixel*/
	   double wvl,       /**<[in] the wavelength in meter*/
	   dmat* theta       /**<[in] angle of rotation of each subaps for polar ccd. NULL for  geometry.*/
	   );
#endif
