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
	      const double rne      /**<[in]Read out noise per pixel per read*/
	      ){
    long np=A->nx*A->ny;
    if(bkgrnd2){
	assert(bkgrnd2->nx*bkgrnd2->ny==np);
	if(bkgrnd2c){
	    assert(bkgrnd2c->nx*bkgrnd2c->ny==np);
	    for(int ix=0; ix<np; ix++){
		A->p[ix]=randp(rstat,A->p[ix]+bkgrnd+bkgrnd2->p[ix])+rne*randn(rstat)-(bkgrndc+bkgrnd2c->p[ix]);
	    }
	}else{
	    for(int ix=0; ix<np; ix++){
		A->p[ix]=randp(rstat,A->p[ix]+bkgrnd+bkgrnd2->p[ix])+rne*randn(rstat)-bkgrndc ;
	    }
	}
    }else{
	if(bkgrnd2c){
	    assert(bkgrnd2c->nx*bkgrnd2c->ny==np);
	    for(int ix=0; ix<np; ix++){
		A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)-bkgrndc+rne*randn(rstat)-bkgrnd2c->p[ix];
	    }
	}else{
	    for(int ix=0; ix<np; ix++){
		A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)-bkgrndc+rne*randn(rstat);
	    }
	}
    }
}
