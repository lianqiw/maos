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
#ifndef AOS_LIB_RAYLEIGH_H
#define AOS_LIB_RAYLEIGH_H
/**
	\file rayleigh.h
	Compute rayleigh backscatter for WFS subaperture
*/
/**
 * @struct rayleigh_t
 * @brief Parameters for the rayleigh() scattering simulation.
 */
typedef struct {
    loc_t *saloc;  /**< Subaperture lower left corner location. */
    dmat *thetax;  /**< LGS on sky angle wrt optical axis (x-component). */
    dmat *thetay;  /**< LGS on sky angle wrt optical axis (y-component). */
    dmat *llt_ox;  /**< LLT location for each LGS (x-coordinate). */
    dmat *llt_oy;  /**< LLT location for each LGS (y-coordinate). */
    dmat *tau;     /**< n*3 matrix: [height, total throughput, scattering fraction]. */
	//Always use double for stable Python interface
    double dsa;            /**< Subaperture size. */
    double sbeam;          /**< Sigma (knee width) of the laser beam Gaussian profile. */
    double dbeam;          /**< Full width of laser beam as clipped by the LLT. */
    double hs;             /**< Altitude/Height of the Laser Guide Star (LGS). */
    double dtheta;         /**< Angular resolution for PSF sampling. */
    int npsf;            /**< Grid size for the PSF (e.g., npsf x npsf). */
} rayleigh_t;
dcell *rayleigh(rayleigh_t *cfg);
void rayleigh_free(rayleigh_t *cfg);
rayleigh_t *rayleigh_setup(real D);
#endif
