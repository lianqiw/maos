/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file cn2est.h
*/
/**
   contains the data related to Cn2 Estimation for each WFS pair.
*/
typedef struct cn2est_pair_t{
   int xstep;      /**<separation step of subapertures along x*/
   int ystep;      /**<separation step of subapertures along y*.*/
   int iht0;       /**<starting iht, may be negative or 0, for minimum height*/
   int iht1;       /**<ending iht (exclusive), for maximum height*/
   int nht;        /**<iht1-iht0+1*/
   int nsep;       /**<number of subaperture separations to use. equal to nhs usually*/
   int wfs0;       /**<first wfs in this pair*/
   int wfs1;       /**<second wfs in this pair*/
} cn2est_pair_t;
/**
   contains the data related to Cn2 Estimation.
 */
typedef struct cn2est_t{
   /*the following are constants*/
   struct cn2est_pair_t* pair; /**<information about each pair*/
   int* wfscov;     /**<Whether this wfs participates in covariance computation.*/
   long nembed;      /**<size of array to embed the LGS gradients into*/
   lmat* embed;      /**<pointers to embed*/
   lmat* mask;       /**<select subapertures that are both full and have neighbors to compute covariance*/

   dmat* overlapi;  /**<1./Number of overlapping subapertures for each separation*/
   int nsa;         /**<Number of subapertures*/
   int nwfs;        /**<number of wfs*/
   int nwfspair;    /**<number of wfs pairs to use for cn2 estimation*/
   int ovs;         /**<Over sampling ratio in building the influence matrix*/

   dcell* Pnk;      /**<Cn2 Estimation forward matrix*/
   dcell* iPnk;     /**<Cn2 Estimation matrix.*/

   dcell* ht;       /**<Estimated Height of the layers*/
   real hmax;     /**<maximum cn2 estimation when keepht!=2*/

   /*the following are used to feed into reconstructor. */
   dmat* htrecon;   /**<layer heights for tomography*/
   dmat* os;        /**<over sampling factor of the layers in htrecon*/
   dmat* dx;        /**<sampling of each layer in reconstruction.*/
   dmat* dmht;      /**<For theta_2 printing*/
   dspcell* wtconvert; /**<to convert wt from wt to wtrecon.*/
   real L0;        /**<outer scale*/

   /*the following are temporary data*/
   dcell* gxs;      /**<gradient x*/
   dcell* gys;      /**<gradient y*/
   ccell* curi;     /**<For FFT*/

   dcell* cov2;     /**<Covariance in 2d*/
   dcell* cov1;     /**<Cut of cov2 along wfs separation*/

   dcell* wt;       /**<Estimated weighting of the layers*/

   /*the following are runtime data for accumulation*/
   ccell* covc;     /**<Accumulation of FFT of Covariance in 2d*/
   int count;       /**<number of time steps we have accumulated the covariance*/

   /*the following are estimation results*/
   dmat* r0;        /**<Estimated r0*/
   real r0m;       /**<averaged r0 from all the pairs.>*/
   dcell* wtrecon;   /**<layer weights for tomography*/
} cn2est_t;

cn2est_t* cn2est_new(const dmat* wfspair,const dmat* wfstheta,const loc_t* saloc,
   const dmat* saa,const real saat,
   const dmat* hs,const dmat* htrecon,int keepht,real l0);
void cn2est_est(cn2est_t* cn2est,int verbose);
void cn2est_free(cn2est_t* cn2est);
void cn2est_push(cn2est_t* cn2est,const dcell* gradol);
cn2est_t* cn2est_all(const dmat* wfspair,dmat* wfstheta,const loc_t* saloc,
   const dmat* saa,const real saat,
   const dmat* hs,const dmat* htrecon,int keepht,real l0,dcell* grad);
void cn2est_reset(cn2est_t* cn2est);
#endif
