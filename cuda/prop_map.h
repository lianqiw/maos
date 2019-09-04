/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_PROP_WRAP_H
#define AOS_CUDA_PROP_WRAP_H
#include "common.h"
/*data to be used by kernel */
struct map2map_t{
    int offdir, offdirx, offdiry;
    int offps, offpsx, offpsy;
    const Real *cc;
    int nxdir,nydir;
    int nxps,nyps;
    Real dispx, dispy;
    Real xratio, yratio;
    int nx, ny;
    int isreverse;/*Indicate this is an reverse*/
    map2map_t *reverse;/*store pointer for the reversely prepared data.*/
    map2map_t(){ /*This is not good. zeros out already initialized childs.*/
	memset(this, 0, sizeof(*this));
    }
    void togpu(map2map_t *pgpu){
	if(reverse){
	    map2map_t *gpureverse;
	    DO(cudaMalloc(&gpureverse, sizeof(map2map_t)));
	    reverse->togpu(gpureverse);
	    delete reverse;
	    reverse=gpureverse;
	}
	DO(cudaMemcpy(pgpu, this, sizeof(map2map_t),cudaMemcpyHostToDevice));  
    }
};
#define WRAP_TX 32

typedef struct{
    int ix0[WRAP_TX];
    int iy0[WRAP_TX];
    int stepx;
    int stepy;
    int nn;
    int nx;
    int ny;
    int ndirx;
    int ndiry;
    int npsx;
    Real *pps;
    Real *pdir;
    Real dispx;
    Real dispy;
    Real xratio;
    Real yratio;//nohelp
    Real fx[9];//9 to avoid bank conflict
    Real fy[9];
}map2map_shared_t;

__global__ void map2map_do(map2map_t *data, Real *const*pdirs, Real *const*ppss, 
			   int ndir, int nps, Real alpha1, const Real *alpha2, char trans);

void map2map_prep(map2map_t*res, const cugrid_t &g_dir, const cugrid_t &gi,
		  Real dispx, Real dispy, const curmat &cc);

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

class map2map:public nonCopyable{
    map2map_t *hdata;
    int nlayer, ndir;
    void deinit(){
	if(hdata){
	    int np=(ndir==0?1:ndir)*nlayer;
	    map2map_t *pcpu=new map2map_t[np];
	    cudaMemcpy(pcpu, hdata, sizeof(map2map_t), cudaMemcpyDeviceToHost);
	    for(int i=0; i<np; i++){
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
    map2map():hdata(0),nlayer(0),ndir(0){};
    void init_l2d(const cugrid_t &out, const dir_t *dir, int _ndir, const cugridcell &in, Real dt=0);
    void init_l2l(const cugridcell &out, const cugridcell &in);
    //from in to out
    void forward(Real *const*out, const Real *const *in,  Real alpha, Real *wt, stream_t &stream){
	map2map_do<<<dim3(4,4,ndir==0?nlayer:ndir),dim3(WRAP_TX,4),0,stream>>>
	    (hdata, out, (Real*const*)in, ndir, nlayer, alpha, wt, 'n');
    }
    //from out to in
    void backward(const Real *const *out, Real *const*in, Real alpha, Real *wt, stream_t &stream){
	map2map_do<<<dim3(4,4,nlayer),dim3(WRAP_TX,4),0,stream>>>
	    (hdata, (Real*const*)out, in, ndir, nlayer, alpha, wt, 't');
    }
    virtual ~map2map(){
	deinit();
    }
    operator bool()const{
	return nlayer||ndir;
    }
};

void map_rot(curcell &out, const curcell &in, const curmat &wfsrot, int dir, stream_t &stream);
#endif
