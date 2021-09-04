/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include "../math/mathdef.h"
/**
   \file mkdtf.h
   Routine to generate detector transfer function and elongation transfer function due to sodium layer.
*/

/**
   contains the data associated with a detector transfer function for a
   subaperture. The PSF is computed as
   \f$\textrm{PSF}=\frac{1}{N^2\sum(\textrm{amp}^2)}|\textrm{fftshift}\mathcal{F}[A
   \exp(-\frac{2\pi}{\lambda}\textrm{opd})]|^2\f$.  The subaperture image is
   computed as
   \f$I=\textrm{si}*\mathcal{F}^{-1}[\mathcal{F}[\textrm{PSF}\times\textrm{nominal}]]\f$
*/
typedef struct dtf_t{
    ccell* nominal;    /**<The FFT of the pixel functions, don't apply if etf exists .*/
    dspcell* si;       /**<The pixel selection*/
    real wvl;          /**<Wavelength*/
    real dtheta;       /**<Sampling of PSF*/
    cmat* Ux;          /**<Special frequency vector along x*/
    cmat* Uy;          /**<Special frequency vector along y*/
    real dxsa;         /**<Subaperture size*/
    real pixthetax;    /**<Pixel size along x (radial)*/
    real pixthetay;    /**<Pixel size along y (radial*/
    int notfx;         /**<FFT size along x*/
    int notfy;         /**<FFT size along y*/
    int pixpsax;       /**<Number of pixels along x*/
    int pixpsay;       /**<Number of pixels along y*/
    int radpix;        /**<1: Pixels are along radial/azimuthal direction*/
    int nwvl;          /**<Number of dtf_t*/
}dtf_t;

typedef struct etf_t{
    ccell* etf;          /**<Store the 2D ETF. nominal is fused in always.*/
    double hs;           /**<Guide star height*/
    int icol;            /**<Store the column index*/
    int nwvl;            /**<Number of dtf_t*/
}etf_t;

dtf_t* mkdtf(const dmat* wvls, /**<List of wavelength*/
    real dxsa,        /**<Subaperture size*/
    real embfac,      /**<Embedding factor (2)*/
    long notfx,       /**<FFT size along x*/
    long notfy,       /**<FFT size along y*/
    long pixpsax,     /**<Number of pixels along x(r)*/
    long pixpsay,     /**<Number of pixels along y(a)*/
    real pixthetax,   /**<Pixel size along x (r)*/
    real pixthetay,   /**<Pixel size along y (a)*/
    const dmat* pixoffx,  /**<offset of image center from center of pixel array, along x or radial*/
    const dmat* pixoffy,  /**<offset of image center from center of pixel array, along y or azimuthal*/
    real pixblur,     /**<Pixel blur sigma(fraction of pixel)*/
    const dcell* pixrot /**<Rotation angle of pixels islands in each subaperture. for polar coordinate only*/
);
etf_t* mketf(const dtf_t* dtfs,  /**<The dtfs*/
    const dcell* sodium, /**<The sodium profile. First column is coordinate.*/
    int icol,     /**<Which sodium profile to use*/
    const dcell* srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
    const dcell* srsa,  /**<Subaperture to LLT distance*/
    real hs,      /**<Guide star focus range*/
    real htel,    /**<Telescope altitude*/
    real za_rad, /**<Telescope zenith angle in radian*/
    int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
);
void dtf_free_do(dtf_t* dtfs);
void etf_free_do(etf_t* etfs);
/**frees dtf_t */
#define dtf_free(dtfs) if(dtfs){dtf_free_do(dtfs); dtfs=NULL;}
/**frees etf_t */
#define etf_free(etfs) if(etfs){etf_free_do(etfs); etfs=NULL;}
dmat* smooth(const dmat* profile, real dxnew);
dcell* smooth_cell(const dcell* profile, real dxnew);
#endif
