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
#ifndef AOS_CUDA_RECON_GEOM_H
#define AOS_CUDA_RECON_GEOM_H
#include "types.h"
#include "prop_wrap.h"

namespace cuda_recon{
struct dir_t{
    Real thetax; /**<Direction of guide star*/
    Real thetay; /**<Direction of guide star*/
    Real hs;     /**<Height of guide star*/
    Real misregx;/**<Misregistration*/
    Real misregy;/**<Misregistration*/
    Real delay;  /**<The delay used for tomo.predict*/
    int skip;    /**<Skip this direction*/
public:
    dir_t():thetax(0),thetay(0),hs(INFINITY),misregx(0),misregy(0),delay(0),skip(0){}
};
/*Data for aperture bi-linear weighting, used in fitting*/
class W01_T{
    curmat W1;    /**< The aperture weighting, piston removal*/
    cusp   W0p;   /**< W0 for partial points*/
    cumat<int>W0f;/**< index for fully illuminated points.*/
    Real   W0v;   /**< maximum Value of W0*/
    int     nxx;  /**< First dimension of grid*/
    curmat pis;   /**< Temporary data*/
public:
    W01_T(const dsp *R_W0, const dmat *R_W1, int R_nxx);
    void apply(Real *restrict out, const Real *in, int ndir, stream_t &stream);
};

/**
   Contains geometry information that may be shared by 
   tomo
   fit
   moao
*/

class curecon_geom:public nonCopyable{
public:
    int npsr, ndm;
    int delay, isim, reconisim;
    cugridcell xmap;/*Grid of xmap*/
    cugridcell xcmap;
    cugridcell amap;
    cugrid_t pmap; /*Pupil map for tomo*/
    cugrid_t fmap; /*Pupil map for fit*/
    W01_T W01;    /**< The aperture weighting defined on floc*/
    long *xnx, *xny;/*do not free*/
    long *anx, *any;/*do not free*/
    long *anloc, *ngrad;/*do not free*/
    Real dt; 
    curecon_geom(const PARMS_T *parms, const RECON_T *recon);
    ~curecon_geom(){ }
};
class map_ray:public nonCopyable{
protected:
    PROP_WRAP_T *hdata;
    int nlayer, ndir;
    void deinit(){
	if(hdata){
	    PROP_WRAP_T *pcpu=new PROP_WRAP_T[ndir*nlayer];
	    cudaMemcpy(pcpu, hdata, sizeof(PROP_WRAP_T), cudaMemcpyDeviceToHost);
	    for(int i=0; i<ndir*nlayer; i++){
		if(pcpu[i].reverse){
		    cudaFree(pcpu[i].reverse);
		}
	    }
	    delete [] pcpu;
	    cudaFree(hdata);
	    hdata=0;
	    nlayer=0;
	    ndir=0;
	}
    }
public:
    map_ray():hdata(0),nlayer(0),ndir(0){};
    void Init_l2d(const cugrid_t &out, const dir_t *dir, int _ndir, const cugridcell &in, Real dt=0);
    void Init_l2l(const cugridcell &out, const cugridcell &in);
    //from in to out
    void forward(Real **out, Real **in,  Real alpha, Real *wt, stream_t &stream){
	gpu_map2map_do<<<dim3(4,4,ndir==0?nlayer:ndir),dim3(PROP_WRAP_TX,4),0,stream>>>
	    (hdata, out, in, ndir, nlayer, alpha, wt, 'n');
    }
    //from out to in
    void backward(Real **out, Real **in, Real alpha, Real *wt, stream_t &stream){
	gpu_map2map_do<<<dim3(4,4,nlayer),dim3(PROP_WRAP_TX,4),0,stream>>>
	    (hdata, out, in, ndir, nlayer, alpha, wt, 't');
    }
    virtual ~map_ray(){
	deinit();
    }
    operator bool()const{
	return nlayer||ndir;
    }
};

}//namespace
#endif
