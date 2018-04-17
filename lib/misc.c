/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file misc.c
   
   Miscellaneous routines.
*/
#include "misc.h"

/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
out. pcalib2 part of bkgrnd2 is calibrated out.  */
void addnoise(dmat *A,              /**<[in/out]The pixel intensity array*/
	      rand_t* rstat,        /**<[in]The random stream*/
	      const double bkgrnd,  /**<[in]Real background in PDEs per pixel per frame*/
	      const double bkgrndc, /**<[in]Removed background in PDEs per pixel per frame*/
	      const dmat *bkgrnd2,  /**<[in]Real background in PDEs of each pixel per frame.*/
	      const dmat *bkgrnd2c, /**<[in]Removed background in PDEs of each pixel per frame.*/
	      const dmat *qe,       /**<[in]Pixel dependent Quantum Efficiency*/
	      const double rne,     /**<[in]Read out noise per pixel per read*/
	      double excess   /**<[in]Excess noise factor*/
	      ){
    long np=A->nx*A->ny;
    assert(!bkgrnd2 || bkgrnd2->nx*bkgrnd2->ny==np);
    assert(!bkgrnd2c || bkgrnd2c->nx*bkgrnd2c->ny==np);
    if(excess<1) excess=1;
    for(int ix=0; ix<np; ix++){
	double tot=A->p[ix]+bkgrnd+(bkgrnd2?bkgrnd2->p[ix]:0);
	double corr=bkgrndc+(bkgrnd2c?bkgrnd2c->p[ix]:0);
	double scale=1;
	if(qe){//the second qe factor is flat-field correction.
	    tot*=qe->p[ix];
	    scale=1./qe->p[ix];
	}
	A->p[ix]=(randp(rstat, tot*excess)+tot*(1.-excess)+rne*randn(rstat))/scale-corr;
    }
}
