/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_CUDA_DATA_H
#define AOS_CUDA_DATA_H
#include <map>
#include "../math/math.h"
#include "../recon/recon.h"
#include "perf.h"
#include "wfs.h"
#include "gpu_sim.h"
extern int NGPU;//Actual number of GPUs being used. 
extern int MAXGPU;
extern Array<int> GPUS;//Used GPUS in cuda (col 0) and nvidia-smi index (col 1). nrow may be more than NGPU
typedef Real AReal;
typedef Real GReal;

typedef struct{
    int ips;
    int nx0;
    int ny0;
    int offx;
    int offy;
    int isim;
    int isim_next;
    map_t *atm;
    Real ox;
    Real oy;
    Real *next_atm;
    pthread_t threads;
}atm_prep_t;
//Global data independent of GPU
class cuglobal_t{
public:
    //Field ordering is important in destructor calling.
    Array<int> evlgpu;
    Array<int> wfsgpu;
    dmat *atmscale; /**<Scaling of atmosphere due to r0 variation*/
    int recongpu;
    int atm_full;   /**<Indicate whether atm is loaded in full to GPU*/
	int atm_dtrat;  /**saves parms->atm.dtrat;*/
	int atm_interp; /**<saves parms->atm.interp */
	int atm_nps;    /**<saves parms->atm.nps;*/
    cuperf_g perf;
    Array<cuwfs_t>wfs;
    Array<atm_prep_t>atm_prep_data;
    curmat mvm;
    cuglobal_t():atmscale(0), recongpu(0), atm_full(0){
    }
    ~cuglobal_t(){
	    dfree(atmscale);
    }
};
//Per GPU data
class cudata_t{
public:
    int igpu;
    /**<for accphi */
    curmat reserve;   /**<Reserve some memory in GPU*/
    cumapcell atm;   /**<atmosphere: array of cumap_t */

    cumapcell dmreal;/**<DM: array of cumap_t */
    cumapcell dmproj;/**<DM: array of cumap_t */
	cumapcell opdr;  /**<Copied from curecon.opdr_map for evl.tomo */
    //int nps; /**<number of phase screens*/
    /*for perfevl */
    cuperf_t perf;
    stream_t perf_stream;/**<Default stream for perfevl. One per GPU. This allows sharing resource per GPU.*/
    /*for wfsgrad */
    Array<cupowfs_t>powfs;//This is created for each GPU.

    /*for recon */
    curecon_t *recon;
    /*for moao*/
    Array<cumapcell> dm_wfs;
    Array<cumapcell> dm_evl;
    /*for mvm*/
    curmat mvm_m;/*the control matrix*/
    NumArray<AReal, Gpu>mvm_a; /*contains act result from mvm_m*mvm_g*/
    NumArray<GReal, Gpu>mvm_g;/*the gradients copied from gpu*/
    stream_t mvm_stream;
    cudata_t():igpu(0),recon(0){
    }
    ~cudata_t(){
	    delete recon;
    }
};

#ifdef __APPLE__
extern pthread_key_t cudata_key;
static inline cudata_t* _cudata(){
    return (cudata_t*)pthread_getspecific(cudata_key);
}
#define cudata _cudata()
#else
extern __thread cudata_t *cudata;
#endif
extern cudata_t **cudata_all;/*use pointer array to avoid misuse. */
extern cuglobal_t *cuglobal;
void gpu_print_mem(const char *msg);
long gpu_get_free_mem(void);
/**
   switch to the next GPU and update the pointer.
*/
static inline void gpu_set(int igpu){
    extern Array<int> GPUS;
    if(igpu>=NGPU){
	    error("Invalid igpu=%d", igpu);
    }
    cudaSetDevice(GPUS[igpu]);
#ifdef __APPLE__
    pthread_setspecific(cudata_key, cudata_all[igpu]);
#else
    cudata=cudata_all[igpu];
#endif
    if(cudata->reserve){
		cudata->reserve.init();//release reserve
    }
}
/**
   returns next available GPU. Useful for assigning GPUs to particular wfs, evl, etc.
*/
static inline int gpu_next(float add=1){
    static float cur=-1;
    cur+=add;
    int ans=(int)ceil(cur);
    if(ans>=NGPU){
	    ans-=NGPU;
	    cur-=NGPU;
    }
    return ans;
}

#endif
