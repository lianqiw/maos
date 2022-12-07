/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

/**
   \file locfft.h

   Routines used to embed OPDs defined on loc_t to square array with embeding
   factor and compute PSF using FFT.
*/

typedef struct{
    loc_t *loc; /*reference to the grid*/
    const dmat *amp;  /*reference to the amplitude map*/
    dmat *wvl;  /*reference to the wavelength*/
    lmat *nembed;/**<size of embedding array (square)*/
    lcell *embed; /**<embedding index*/
    real ampsum;/**<sum(amp)*/
    real ampnorm;/**<sum(amp.*amp)*/
    dcell *fieldmask;/**<Masking the PSF in fourier domain*/
    real fieldstop;
}locfft_t;

locfft_t *locfft_init(loc_t *loc, const dmat *amp, const dmat *wvl, const lmat *fftsize, real oversize, real fieldstop); 
void locfft_free(locfft_t *locfft);
void locfft_psf(ccell **psfs, const locfft_t *locfft, const dmat *opd, const lmat *psfsize, int sum2one);
void locfft_fieldstop(const locfft_t *locfft, dmat *opd, const dmat *wvlwts);
void fresnel_prop(cmat **out, real *dxout, const cmat *in, real dxin, real wvl, real z, real scale, int method);