/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_FDPCG_H
#define AOS_CUDA_FDPCG_H
#include "solve.h"
namespace cuda_recon{
typedef struct GPU_FDPCG_T{
    int nx;
    int ny;
    Real scale;
}GPU_FDPCG_T;
class curecon_geom;
class cufdpcg_t:public cucgpre_t{
    curecon_geom *grid;
    cumat<int> perm;   /**<permutation vector for fdpcg*/
    cuccell Mb;  /**<The main fdpcg block matrix*/
    cuarray<cufftHandle> fft;
    cuarray<cufftHandle> ffti;
    int      fftnc;/*total number of ffts*/
    int     *fftips;/*starting ips for each fft.*/
    cuccell xhat1, xhat2;
    int nb, bs, nby, nbz; 
    int scale;
    //int ixskip;//skip the ix due to redundant
    cumat<GPU_FDPCG_T> fddata;
    cufdpcg_t &operator=(const cufdpcg_t&in);
    cufdpcg_t(const cufdpcg_t&);
public:
    virtual ~cufdpcg_t(){
	free(fftips);
    }
    cufdpcg_t(FDPCG_T *fdpcg, curecon_geom *_grid);
    void update(FDPCG_T *fdpcg);
    void P(curcell &xout, const curcell &xin, stream_t &stream);
};
}//namespace
#endif
