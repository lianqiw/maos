/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_RECON_H
#define AOS_CUDA_RECON_H
#include "types.h"
typedef struct{
    curcell *gradin; /**< The grad to operator on*/
    curcell *neai;
    curcell *opdwfs;/**<Temporary*/
    curcell *grad;  /**<Temporary*/
    curcell *opdr;  /**<Reconstructed atm on xloc. Don't free to have warm
		       restart. Free with new seed*/
    curcell *dmfit; /**<Reconstructed DM command. Don't free to have warm
		       restart. Free with new seed.*/
    cudaStream_t     *wfsstream;
    cublasHandle_t   *wfshandle;
    cusparseHandle_t *wfssphandle;
    cudaEvent_t      *wfsevent;
    
    cudaStream_t     *psstream;

    cudaStream_t     *fitstream;
    cusparseHandle_t *fitsphandle;
    cublasHandle_t   *fithandle;

    cudaStream_t     cgstream;
    cublasHandle_t   cghandle;
    curcell *PTT;  /**< Global tip/tilt, DF removal */
    curcell *PDF;  /**< Differential focus removal */
    float *l2c;    /**< Laplacian */
    int   *zzi;    /**< Piston constraint coordinate*/
    float *zzv;    /**< Piston constraint value*/
    curmat *W1;    /**< The aperture weighting, piston removal*/
    cusp   *W0p;   /**< W0 for partial points*/
    int    *W0f;   /**< index for fully illuminated points.*/
    int     nW0f;  /**< Number of fully illuminated points.*/
    float   W0v;   /**< maximum Value of W0*/
    curcell *opdfit; /**<OPDs defined on ploc for fitting.*/
    curcell *opdfit2;/**<OPDs defined on ploc for fitting.*/
    cumuv_t FR;
    cumuv_t FL;
    int reconisim;
}curecon_t;


extern curecon_t *curecon;
void gpu_TomoR(curcell **xout, const void *A, curcell *grad, const float alpha);
void gpu_TomoL(curcell **xout, const float beta, const void *A, const curcell *xin, const float alpha);
void gpu_FitR(curcell **xout, const void *A, const curcell *xin, const float alpha);

#endif
