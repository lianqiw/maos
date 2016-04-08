/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_PYWFS_H
#define AOS_PYWFS_H
    
/**
   Parameters used by Pyramid WFS
 */
typedef struct PYWFS_T{
    double hs;         /**<Height of guide star*/
    double poke;       /**<How much to poke for mkg*/
    double modulate;   /**<Amount of modulation in radian*/
    int modulpos;      /**<Number of positions per modulation cycle*/
    int iwfs0;         /**<First iwfs for this powfs*/
    dmat *wvlwts;      /**<parms->powfs.wvlwts*/
    loc_t *loc;        /**<Pupil plane grid*/
    dmat  *amp;        /**<Pupil plane amplitude map*/
    locfft_t *locfft;  /**<First fft to form PSF*/
    ccell *pyramid;    /**<OPD of pyramid. Angular size of clear aperture is different*/
    cmat *nominal;     /**<For sampling results onto detector*/
    dspcell *si;       /**<For sampling results onto detector*/
    double gain;       /**<Optical gain of PYWFS*/
    dmat *saa;         /**<Subaperture area. Average is one*/
    dmat *gradoff;     /**<Gradient of a flat wavefront*/
    dmat *GTT;          /**<Response of TT mode.*/
}PYWFS_T;

void pywfs_setup(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs);
void pywfs_free(PYWFS_T *pywfs);
void pywfs_grad(dmat **pgrad, const PYWFS_T *pywfs, const dmat *ints);
void pywfs_fft(dmat **ints, const PYWFS_T *pywfs, const dmat *opd);
dmat* pywfs_mkg(const PYWFS_T *pywfs, const loc_t* ploc, const dmat *mod, double displacex,  double displacey, double scale);
dmat *pywfs_tt(const PYWFS_T *pywfs);

#endif
