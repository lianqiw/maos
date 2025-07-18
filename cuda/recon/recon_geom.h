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
#ifndef AOS_CUDA_RECON_GEOM_H
#define AOS_CUDA_RECON_GEOM_H
#include "../math/cumath.h"
#include "../../maos/parms.h"
#include "../../maos/types.h"


/*Data for aperture bi-linear weighting, used in fitting*/
class w01_t{
	curmat W1;    /**< The aperture weighting, piston removal*/
	cusp   W0p;   /**< W0 for partial points*/
	cuimat W0f;   /**< index for fully illuminated points.*/
	Real   W0v;   /**< maximum Value of W0*/
	int     nxx;  /**< First dimension of grid*/
	mutable curmat pis;   /**< Temporary data*/
public:
	w01_t(const dsp* R_W0, const dmat* R_W1, int R_nxx);
	void apply(curcell& out, const curcell& in, stream_t& stream) const;
};

/**
   Contains geometry information that may be shared by
   tomo
   fit
   moao
*/

class curecon_geom:public nonCopyable{
public:
	int npsr=0, ndm=0;
	int delay=0, reconisim=0;
	Real dt=0;
	long* anx=0, *any=0;/*do not free*/
	long* anloc=0, *ngrad=0;/*do not free*/
	long* xnx=0, *xny=0;/*do not free*/
	cugrid_t pmap; /*Pupil map for tomo*/
	cugrid_t fmap; /*Pupil map for fit*/
	cugridcell amap;
	cugridcell xmap;/*Grid of xmap*/
	cugridcell xcmap;
	w01_t W01;    /**< The aperture weighting defined on floc*/
	
	curecon_geom(const parms_t* parms, const recon_t* recon);
	~curecon_geom(){}
};


#endif
