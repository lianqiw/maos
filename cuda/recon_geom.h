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
#ifndef AOS_CUDA_RECON_GEOM_H
#define AOS_CUDA_RECON_GEOM_H
#include "types.h"
/**
   Contains geometry information that may be shared by 
   tomo
   fit
   moao
*/
class curecon_geom{
public:
    int npsr, ndm;
    int delay, isim, isimr;
    curcell *cubic_cc;/*Cubic influence function of system DM*/
    cugrid_t *amap;/*Grid of amap*/
    cugrid_t *acmap;
    cugrid_t *xmap;/*Grid of xmap*/
    cugrid_t *xcmap;
    cugrid_t pmap; /*Pupil map for tomo*/
    cugrid_t fmap; /*Pupil map for fit*/
    W01_T *W01;    /**< The aperture weighting defined on floc*/
    long *xnx, *xny;/*do not free*/
    long *anx, *any;/*do not free*/
    long *anloc, *ngrad;/*do not free*/
    float *vx, *vy, dt; /*wind speed*/
    void init(const PARMS_T *parms, const RECON_T *recon);
    curecon_geom(const PARMS_T *parms=0, const RECON_T *recon=0):
	npsr(0),ndm(0),isim(0),isimr(0),cubic_cc(0),
	amap(0),acmap(0),xmap(0),xcmap(0),W01(0),
	xnx(0),xny(0),anx(0),any(0),
	anloc(0),ngrad(0),
	vx(0),vy(0),dt(0),delay(0){
	if(parms){
	    init(parms, recon);
	}
    }
    ~curecon_geom(){
	delete[] amap;
	delete[] acmap;
	delete[] xmap;
	delete[] xcmap;
	delete[] vx;
	delete[] vy;
	delete W01;
	delete cubic_cc;
    }
};
#endif
