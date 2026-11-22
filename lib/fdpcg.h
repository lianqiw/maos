/*
  Copyright 2009-2026 Lianqi Wang

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
#ifndef AOS_LIB_FDPCG_H
#define AOS_LIB_FDPCG_H
#include "loc.h"
/**
   \file lib/fdpcg.h
   Fourier domain preconditioner (FDPCG) for the tomography with conjugate
   gradients. The preconditioner approximates the tomography left hand side
   matrix in the Fourier domain with periodic boundary conditions, where it is
   block diagonal after grouping the spatial frequencies of all layers that
   couple through the subaperture sampling.
*/
/**
   contains data for Fourier Domain Preconditioner.
 */
typedef struct fdpcg_t{
    loccell *xloc;  /**<record recon->xloc*/
    csp *Minv;     /**<inverse of fourier domain preconditioner matrix M.*/
    ccell *Mbinv;  /**<block version of Minv. (in permuted order)*/
    lmat *perm;    /**<Permutation vector to get block diagonal matrix*/
    lmat *permhf;  /**<Permutation vector to be used when complex2real fft is used. Size is about half of perm.*/
    lmat *nx;      /**<Grid size along x of each layer*/
    lmat *ny;      /**<Grid size along y of each layer*/
    long nxtot;    /**<Total number of reconstructed phase points*/
    int square;    /**<Whether xloc is square*/
    int scale;     /**<Do we need to scale after fft.*/
    long nbx;      /**<Basic frequency range in x*/
    long nby;      /**<Basic frequency range in y nb=nbx*nby.*/
    long bs;       /**<Size of each diagonal block.*/
    /*Temperary Data used during execution*/
    /*Cannot save xhat here because we may be calling fdpcg from multiple
     * threads simultanously during assembly of control matrix*/
    //ccell *xhat, *xhat2, /**<FFT of each xout and multiplied. */
}fdpcg_t;
csp* fdpcg_saselect(long nx, long ny, real dx, loc_t* saloc, real* saa);
lmat* fdpcg_perm(const long* nx, const long* ny, const real* os, int bs, int nps, int shift, int half);
void fdpcg_g(cmat** gx, cmat** gy, long nx, long ny, real dx, real dsa);
csp* fdpcg_prop_trans(long nps, long nxp, long nyp, long* nx, long* ny, real dx, real* dispx, real* dispy);
csp* fdpcg_mhat(long nps, long* nx, long* ny, long nxp, long nyp, real dxp, loc_t* saloc, real* saa, const cmat* invpsd, long nwfs, const dmat* dispx, const dmat* dispy, const real* neai, const long* iwfs, int save);
fdpcg_t* fdpcg_new(csp* Mhat, long nps, const long* nx, const long* ny, const real* os, int half, loccell* xloc, int square, int save);
void fdpcg_apply(dcell** xout, const fdpcg_t* fdpcg, const dcell* xin, int nthread);
void fdpcg_free(fdpcg_t* fdpcg);
#endif
