/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "../math/mathdef.h"
#include "locfft.h"
#include "mkh.h"
/**
 * Pyramid WFS configurations. //Todo: move to pywfs.h after pywfs.h is moved to lib folder
 * */
typedef struct pywfs_cfg_t{
	int order;		 /**<Number of pixels per pupil in each dimension.*/
	real siglev;     /**<Nominal siglev per subaperture*/
	dmat *wvl; 		 /**<wavelength*/

	//optional parameters if loc is provided
	real dx;		 /**<OPD sampling. Use loc->dx if not present*/
	real D;		     /**<Aperture diameter. Use loc_diam if not present*/

	//optional parameters. default values are used if omitted (0).
	int nside;       /**<Number of sides. can be 4 (default) or 3. 4-side is traditional pyramid while 3-side is much easier to make.*/
	int raw;         /**<1: use normalized ints of each sub-pupil as gradient. 0 (default): use difference between sub-pupils*/
	real dsa;	     /**<Subaperture size (default is D/order)*/
	real saat;		 /**<drop subapertures below this fraction (default is 0)*/
	real hs;         /**<Height of guide star (default is inf)*/
	real hc;         /**<Conjugation height of WFS pupil (default is 0)*/
	real modulate;   /**<Amount of modulation in radian (default is 5 lambda/D)*/
	int modulpos;    /**<Number of positions per modulation cycle per ring(default is 32)*/
	int modulring;   /**<Number of rings within the maximum radius to modulate (default is 1, use more to simulate extended object)*/
	int sigmatch;    /**<Scale gradients by matching intensity (1: locally, 2 (default): globally).*/
	real poke;       /**<How much to poke for mkg. (default is 1e-7 m)*/
	dmat *wvlwts;    /**<Wavelength weights*/
	real pixblur;	 /**<pixel blur (default is 0.3) */
	real fieldstop;  /**<size of field stop in arcsec (default is none).*/
	//for testing
	int modulpos_i;  /**<For testing: Index of modulate position to use between 0 and modulpos*/
	//The following are used to simulate the implementation error.
	dmat *psx;       /**<pyramid WFS pupil shift along x (in pixel). pupil ordering: -x+y, +x+y, -x-y, +x-y.*/
	dmat *psy;       /**<pyramid WFS pupil shift along y (in pixel).*/
	real flate;      /**<pyramid flat edge angular width (in arcsec)*/
	real flatv;      /**<pyramid flat vertex angular width (in arcsec)*/
	real pupelong;   /**<pyramid pupil (detector) elongation ratio (long axis / short axis - 1).*/
}pywfs_cfg_t;

/**
  Parameters used by Pyramid WFS. 
 */
typedef struct pywfs_t{
  const pywfs_cfg_t *cfg; /**<Configuration. do not try to free.*/
  int iwfs0;         /**<First iwfs for this powfs*/
  int gpu;           /**<Whether GPU can be used*/
  loc_t *loc;        /**<Pupil plane grid*/
  dmat *amp;        /**<Pupil plane amplitude map*/
  loc_t *saloc;	     /**<powfs.saloc*/
  locfft_t *locfft;  /**<First fft to form PSF*/
  ccell *pyramid;    /**<OPD of pyramid. Angular size of clear aperture is different*/
  cmat *nominal;     /**<For sampling results onto detector*/
  dspcell *si;       /**<For sampling results onto detector*/
  dmat *sioff;       /**<Offset for si*/
  dmat *saa;         /**<Subaperture area. Average is one*/
  dmat *gradoff;     /**<Gradient of a flat wavefront*/
  dmat *GTT;          /**<Response of TT mode.*/
  dmat *pupilshift;  /**<Pupil shift. 4x2.*/
  loccell *msaloc;   /**<Mishaped saloc of each sub-pupil due to optical effects*/
  dmat *opdadd;      /**<Aberrations along the WFS path (on locfft grid)*/
  real gain;         /**<Optical gain of PYWFS*/
}pywfs_t;
pywfs_t* pywfs_new(pywfs_cfg_t *pycfg, loc_t *loc, const dmat *amp);
#define pywfs_ng(pycfg) ((pycfg && (pycfg->raw||pycfg->nside<3))?pycfg->nside:2)
void pywfs_free(pywfs_t *pywfs);
void pycfg_free(pywfs_cfg_t *pycfg);
void pywfs_grad(dmat **pgrad, const pywfs_t *pywfs, const dmat *ints);
void pywfs_ints(dmat **ints, const pywfs_t *pywfs, const dmat *opd, real siglev);
dmat* pywfs_mkg(pywfs_t *pywfs, const loc_t* ploc, const char *distortion, 
		const dmat *mod, const dmat *opdadd, real displacex,  real displacey);
dmat *pywfs_tt(const pywfs_t *pywfs);
void pywfs_simu(dmat **ints, dmat **grad, pywfs_cfg_t *pycfg, int order, dmat *wvl, real siglev, loc_t *loc, const dmat *amp, const dmat *opd);
void pywfs_gain_calibrate(pywfs_t *pywfs, const dmat *grad, real r0);
void pywfs_test(pywfs_t *pywfs);
#endif
