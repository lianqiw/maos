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
#ifndef AOS_CUDA_DATA_H
#define AOS_CUDA_DATA_H
#include <map>
#include "types.h"
extern int NGPU;
typedef Real ATYPE;
typedef Real GTYPE;
namespace cuda_recon{
class curecon_t;
class curecon_geom;
}
struct cuperf_t;
class cuwloc_t;
class cuwfs_t;

typedef struct cudata_t{ 
    int igpu;
    static int recongpu;
    static int *evlgpu;
    static int *wfsgpu;
    std::map<uint64_t, void*> *memhash;/*For reuse constant GPU memory*/
    std::map<void *, int> *memcount; /*Store count of reused memory*/
    std::map<void*, void*> *memcache;/*For reuse temp array for type conversion.*/
    pthread_mutex_t memmutex;
    /**<for accphi */
    void *reserve;   /**<Reserve some memory in GPU*/
    cumap_t *atm;   /**<atmosphere: array of cumap_t */
    static dmat *atmscale; /**<Scaling of atmosphere due to r0 variation*/
    cumap_t *dmreal;/**<DM: array of cumap_t */
    cumap_t *dmproj;/**<DM: array of cumap_t */
    int nps; /**<number of phase screens*/
    int ndm; /**<number of DM.*/
    /*for perfevl */
    cuperf_t *perf;
    /*for wfsgrad */
    cuwloc_t *powfs;
    static cuwfs_t *wfs;
    /*for recon */
    cuda_recon::curecon_t *recon;
    /*for moao*/
    cumap_t **dm_wfs;
    cumap_t **dm_evl;
    /*for mvm*/
    curmat *mvm_m;/*the control matrix*/
    ATYPE *mvm_a; /*contains act result from mvm_m*mvm_g*/
    ATYPE **mvm_a2;/*contains act copied from other gpus for sum*/
    GTYPE *mvm_g;/*the gradients copied from gpu*/
    stream_t *mvm_stream;
    cudata_t(){
	memset(this, 0, sizeof(cudata_t));
	memhash=new std::map<uint64_t, void*>;
	memcount=new std::map<void*, int>;
	memcache=new std::map<void*, void*>;
	pthread_mutex_init(&memmutex, 0);
    }
}cudata_t;
class cudata_t;
#ifdef __APPLE__
extern pthread_key_t cudata_key;
inline cudata_t* _cudata(){
    return (cudata_t*)pthread_getspecific(cudata_key);
}
#define cudata _cudata()
#else
extern __thread cudata_t *cudata;
#endif
extern cudata_t *cudata_all;/*use pointer array to avoid misuse. */

void gpu_print_mem(const char *msg);
long gpu_get_mem(void);
/**
   switch to the next GPU and update the pointer.
*/
inline void gpu_set(int igpu){
    extern int *GPUS;
    if(igpu>=NGPU){
	error("Invalid igpu=%d", igpu);
    }
    cudaSetDevice(GPUS[igpu]);
#ifdef __APPLE__
    pthread_setspecific(cudata_key, &cudata_all[igpu]);
#else
    cudata=cudata_all+igpu;
#endif
    if(cudata->reserve){
	cudaFree(cudata->reserve);
	cudata->reserve=NULL;
    }
}
/**
   returns next available GPU. Useful for assigning GPUs to particular wfs, evl, etc.
*/
inline int gpu_next(float add=1){
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
