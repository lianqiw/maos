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
#ifndef AOS_MAOS_CN2EST_H
#define AOS_MAOS_CN2EST_H

/**
   contains the data related to Cn2 Estimation for each WFS pair.
*/
typedef struct CN2PAIR_T{
    int (**sapair)[2];/**<record pair of subapertures to use for each separation for each wfs pair*/
    int *nsapair;   /**<number of subaperture pairs for each separation for each wfs pair*/
    double dtheta;  /**<separation between the stars*/
    double beta;    /**<angle of separation between the stars*/
    int xstep;      /**<separation step of subapertures along x*/
    int ystep;      /**<separation step of subapertures along y*.*/
    int nht;        /**<number of height bins for each wfs pair*/ 
    int nsep;       /**<number of subaperture separations to use. equal to nhs usually*/
    int wfs0;       /**<first wfs in this pair*/
    int wfs1;       /**<second wfs in this pair*/
} CN2PAIR_T;
/**
   contains the data related to Cn2 Estimation.
 */
typedef struct CN2EST_T{
    int *wfscov;     /**<Whether this wfs participates in covariance computation.*/
    long nembed;      /**<size of array to embed the LGS gradients into*/
    long *embed;      /**<pointers to embed*/
    dcell *gxs;      /**<gradient x*/
    dcell *gys;      /**<gradient y*/
    dcell *cxs;      /**<curvature x*/
    dcell *cys;      /**<curvature y*/
    dcell *cc;       /**<cross-covariance for each pair*/
    int nwfspair;    /**<number of wfs pairs to use for cn2 estimation*/
    int ovs;         /**<Over sampling ratio in building the influence matrix*/
    int count;       /**<number of time steps we have accumulated the covariance*/
    struct CN2PAIR_T *pair; /**<information about each pair*/
    dcell *Pkn;      /**<Cn2 Estimation forward matrix*/
    dcell *iPkn;     /**<Cn2 Estimation matrix.*/
    dcell *wt;       /**<Estimated weighting of the layers*/
    dcell *ht;       /**<Estimated Height of the layers*/
    double hmax;     /**<maximum cn2 estimation when keepht!=2*/
    dmat *r0;        /**<Estimated r0*/
    /*the following are used to feed into reconstructor. */
    dmat *htrecon;   /**<layer heights for tomography*/
    dmat *os;        /**<over sampling factor of the layers in htrecon*/
    dmat *dx;        /**<sampling of each layer in reconstruction.*/
    dmat *dmht;      /**<For theta_2 printing*/
    dcell *wtrecon;   /**<layer weights for tomography*/
    spcell *wtconvert; /**<to convert wt from wt to wtrecon.*/
    double r0m;       /**<averaged r0 from all the pairs.>*/
    double l0;        /**<outer scale*/
} CN2EST_T;

CN2EST_T *cn2est_new(dmat *wfspair, dmat *wfstheta, loc_t *saloc, dmat *saa, const double saat, 
		     const double hs, dmat *htrecon, double hmax, int keepht, double l0);
void cn2est_est(CN2EST_T *cn2est, int verbose, int reset);
void cn2est_free(CN2EST_T *cn2est);
void cn2est_push(CN2EST_T *cn2est, dcell *gradol);
#endif