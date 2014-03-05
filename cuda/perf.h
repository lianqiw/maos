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
#ifndef AOS_CUDA_PERF_H
#define AOS_CUDA_PERF_H
#include "types.h"
/**
   cell array of nevl*1 is allocated using static data member to avoid duplicate.
 */
struct cuperf_t{
    static int *nembed;
    static int *psfsize;
    static float *wvls;    
    static cudaStream_t *stream;
    static cublasHandle_t *handle;
    static cufftHandle *plan;

    static curcell *surf;
    static curcell *opd;
    static curcell *psfcl;
    static curcell *psfcl_ngsr;
    static curcell *opdcov;
    static curcell *opdcov_ngsr;
    static curcell *opdmean;
    static curcell *opdmean_ngsr;
    static curcell *cc;
    static scell   *ccb;

    pthread_mutex_t mutex;
    culoc_t *locs;
    culoc_t ***locs_dm;
    float   *amp;
    int    **embed;

    cuccell *wvf;
    cuccell *psfs;
    curcell *psfol;
    curmat  *opdcovol;
    curmat  *opdmeanol;
    cuperf_t(){
	memset(this, 0, sizeof(*this));
    }
};
#endif
