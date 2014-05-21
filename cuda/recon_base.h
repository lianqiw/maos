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
#include "prop_wrap.h"

namespace cuda_recon{
struct dir_t{
    Real thetax;
    Real thetay;
    Real hs;
    int skip;
};
/*Data for aperture bi-linear weighting, used in fitting*/
class W01_T{
    curmat *W1;    /**< The aperture weighting, piston removal*/
    cusp   *W0p;   /**< W0 for partial points*/
    cumat<int>*W0f;   /**< index for fully illuminated points.*/
    Real   W0v;   /**< maximum Value of W0*/
    int     nxx;   /**< First dimension of grid*/
    curmat *pis;   /**< Temporary data*/
public:
    W01_T(const dsp *R_W0, const dmat *R_W1, int R_nxx);
    ~W01_T(){
	delete W1;
	delete W0p;
	delete W0f;
	delete pis;
    }
    void apply(Real *restrict out, const Real *in, int ndir, stream_t &stream);
};

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
    cugrid_t *xmap;/*Grid of xmap*/
    cugrid_t *xcmap;
    cugrid_t *amap;
    cugrid_t pmap; /*Pupil map for tomo*/
    cugrid_t fmap; /*Pupil map for fit*/
    W01_T *W01;    /**< The aperture weighting defined on floc*/
    long *xnx, *xny;/*do not free*/
    long *anx, *any;/*do not free*/
    long *anloc, *ngrad;/*do not free*/
    Real dt; 
    curecon_geom(const PARMS_T *parms=0, const RECON_T *recon=0);
    ~curecon_geom(){
	delete[] xmap;
	delete[] amap;
	delete[] xcmap;
	delete W01;
    }
};
class map_ray{
protected:
    PROP_WRAP_T *hdata;
    int nlayer, ndir;
public:
    map_ray():hdata(0),nlayer(0),ndir(0){};
    //from in to out
    void forward(Real **out, Real **in,  Real alpha, Real *wt, stream_t &stream){
	gpu_prop_grid_do<<<dim3(4,4,ndir==0?nlayer:ndir),dim3(PROP_WRAP_TX,4),0,stream>>>
	    (hdata, out, in, ndir, nlayer, alpha, wt, 'n');
    }
    //from out to in
    void backward(Real **out, Real **in, Real alpha, Real *wt, stream_t &stream){
	gpu_prop_grid_do<<<dim3(4,4,nlayer),dim3(PROP_WRAP_TX,4),0,stream>>>
	    (hdata, out, in, ndir, nlayer, alpha, wt, 't');
    }
    virtual ~map_ray(){
	PROP_WRAP_T *pcpu=new PROP_WRAP_T[ndir*nlayer];
	cudaMemcpy(&pcpu, hdata, sizeof(PROP_WRAP_T), cudaMemcpyDeviceToHost);
	for(int i=0; i<ndir*nlayer; i++){
	    if(pcpu[i].reverse){
		cudaFree(pcpu[i].reverse);
	    }
	}
	delete [] pcpu;
	cudaFree(hdata);
    }
};
/*Ray tracing from one/multiple layers to one/multiple directions*/
class map_l2d:public map_ray{
public:
    map_l2d(const cugrid_t &out, dir_t *dir, int _ndir, const cugrid_t *in, int _nlayer, Real dt=0);
};

/*Ray tracing from layer to layer, for caching.*/
class map_l2l:public map_ray{
public:
    map_l2l(const cugrid_t *out, const cugrid_t *in, int _nlayer);
};

}//namespace
#endif
