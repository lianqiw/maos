/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_PERF_H
#define AOS_CUDA_PERF_H
#include "types.h"
/**
   cell array of nevl*1 is allocated using static data member to avoid duplicate.
 */
struct cuperf_t{
    static int nevl;
    static cuarray<int> nembed;
    static cuarray<int> psfsize;
    static cuarray<Real> wvls;
    static cuarray<cudaStream_t>   stream;
    static cuarray<cublasHandle_t> handle;
    static cuarray<cufftHandle>    plan;
    static curcell surf;
    static curcell opd;
    static curcell psfcl;
    static curcell psfcl_ngsr;
    static curcell opdcov;
    static curcell opdcov_ngsr;
    static curcell opdmean;
    static curcell opdmean_ngsr;
    static curcell cc_ol, cc_cl, coeff;
    static Real **ccb_ol, **ccb_cl;
    static pthread_mutex_t perfmutex;
    culoc_t locs;
    cuarray<cuarray<culoc_t> > locs_dm;
    curmat imcc;
    curmat amp;
    int    **embed;

    cuccell wvf;
    cuccell psfs;
    curcell psfol;
    curmat  opdcovol;
    curmat  opdmeanol;
    cuperf_t():embed(0){
    }
    ~cuperf_t();
};
#endif
