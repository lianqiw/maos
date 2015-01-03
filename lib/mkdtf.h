/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
typedef struct DTF_T{
    ccell *nominal;      /**<The FFT of the pixel functions*/
    dspcell *si;         /**<The pixel selection*/
    double wvl;          /**<Wavelength*/
    double dtheta;       /**<Sampling of PSF*/
    cmat *Ux;            /**<Special frequency vector along x*/
    cmat *Uy;            /**<Special frequency vector along y*/
    double dxsa;         /**<Subaperture size*/
    long ncompx;         /**<FFT size along x*/
    long ncompy;         /**<FFT size along y*/
    int radpix;          /**<1: Pixels are along radial/azimuthal direction*/
    int radrot;          /**<For radial format CCD, rotate PSF/OTF into r/a coord. uses less memory*/
    int fused;           /**<Whether the DTF has been fused to ETF*/
}DTF_T;

typedef struct ETF_T{
    ccell *p1;           /**<Store the ETF along radial direction when radrot==1*/
    ccell *p2;           /**<Store the 2D ETF when radrot==0*/
}ETF_T;

DTF_T *mkdtf(dmat *wvls, /**<List of wavelength*/
	     double dxsa,/**<Subaperture size*/
	     double embfac,/**<Embedding factor (2)*/
	     long ncompx,/**<FFT size along x*/
	     long ncompy,/**<FFT size along y*/
	     long pixpsax,/**<Number of pixels along x(r)*/
	     long pixpsay,/**<Number of pixels along y(a)*/
	     double pixthetax,/**<Pixel size along x (r)*/
	     double pixthetay,/**<Pixel size along y (a)*/
	     double pixoffx,  /**<offset of image center from center of detector*/
	     double pixoffy,  /**<offset of image center from center of detector*/
	     double pixblur,  /**<Pixel blur (fraction of pixel)*/
	     dcell *srot, /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     int radpix,  /**<1: Pixels are along radial/azimuthal direction*/
	     int radrot  /**<For radial format CCD, rotate PSF/OTF into r/a coord. uses less memory*/
    );
ETF_T *mketf(DTF_T *dtfs,  /**<The dtfs*/
	     double hs,    /**<Guide star focus range*/
	     dcell *sodium, /**<The sodium profile. First column is coordinate.*/
	     int icol,     /**<Which sodium profile to use*/
	     int nwvl,     /**<Number of wavelength*/
	     dcell *srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     dcell *srsa,  /**<Subaperture to LLT distance*/
	     double za,    /**<Zenith angle*/
	     int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
    );
void dtf_free_do(DTF_T *dtfs, int nwvl);
void etf_free_do(ETF_T *etfs, int nwvl);
#define dtf_free(dtfs, nwvl) ({dtf_free_do(dtfs, nwvl); dtfs=0;})
#define etf_free(etfs, nwvl) ({etf_free_do(etfs, nwvl); etfs=0;})
dmat* smooth(dmat *profile, double dxnew);
#endif
