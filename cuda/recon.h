/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "pcg.h"

typedef struct W01_T{
    curmat *W1;    /**< The aperture weighting, piston removal*/
    cusp   *W0p;   /**< W0 for partial points*/
    int    *W0f;   /**< index for fully illuminated points.*/
    int     nW0f;  /**< Number of fully illuminated points.*/
    float   W0v;   /**< maximum Value of W0*/
    W01_T(){
	memset(this, 0, sizeof(W01_T));
    }
}W01_T;
typedef struct cumoao_t{
    curcell *fitNW;
    cuspcell *actslave;
    float *cubic_cc;   
    float dxa; /**<sampling of amap*/
    float oxa; /**<original of amap*/
    float oya; /**<original of amap*/
    int nxa;   /**<size of amap*/
    int nya;   /**<size of amap*/
    float dxf; /**<sampling of fmap*/
    float oxf; /**<original of fmap*/
    float oyf; /**<original of fmap*/
    int nxf;   /**<size of fmap*/
    int nyf;   /**<size of fmap*/
    /*aperture weighting*/
    W01_T *W01;
    curcell *rhs;/*right hand side.*/
    float *pis; /*a temporary variable*/
    curmat *xp; /*temporary array*/
    curmat *xp2;/*temporary array*/
    curmat *tmpNW;/*temporary array*/
    cumoao_t(){
	memset(this, 0, sizeof(*this));
    }
}cumoao_t;
typedef struct cufdpcg_t{
    int *perm;   /**<permutation vector for fdpcg*/
    long nxtot;  /**<total number of points*/
    cuccell *Mb;  /**<The main fdpcg block matrix*/
    cufftHandle *fft;
    cufftHandle *ffti;
    int      fftnc;/*total number of ffts*/
    int     *fftips;/*starting ips for each fft.*/
    cuccell *xhat1;
    int nby, nbz; 
    int scale;
    int half;
    ~cufdpcg_t(){
	cudaFree(perm);
	free(fftips);
	delete Mb;
	delete xhat1;
    }
    cufdpcg_t(){
	memset(this, 0, sizeof(*this));
    }
}cufdpcg_t;
typedef struct curecon_t{
    curcell *gradin; /**< The grad to operator on*/
    curcell *neai;
    curcell *opdwfs;/**<Temporary*/
    curcell *grad;  /**<Temporary*/
    curcell *opdr;  /**<Reconstructed atm on xloc. Don't free to have warm restart. Free with new seed*/
    curcell *opdr_vec; /**<Referencing opdr in vector form*/
    curcell *dmfit; /**<Reconstructed DM command. Don't free to have warm restart. Free with new seed.*/
    curcell *dmfit_vec;/**<Referencing dmfit in vector form.*/
    cufdpcg_t *fdpcg;
    stream_t    *psstream;
    stream_t    *dmstream;
    stream_t    *fitstream;
    stream_t    *cgstream;

    curcell *PTT;  /**< Global tip/tilt */
    curcell **PDF;  /**< Differential focus removal */
    curmat *ttf;

    float *l2c;    /**< Laplacian */
    int   *zzi;    /**< Piston constraint coordinate*/
    float *zzv;    /**< Piston constraint value*/
    W01_T *W01;    /**< The aperture weighting,*/
    curcell *opdfit; /**<OPDs defined on ploc for fitting.*/
    curcell *opdfit2;/**<OPDs defined on ploc for fitting.*/
    curmat *pis;     /**<contains result of W1'*opdfit*/
    float **cubic_cc;
    cumuv_t FR;
    cumuv_t FL;
    curcell *MVM;
    float (*floc)[2];/**<recon->floc*/
    int nfloc;       /**<recon->floc->nloc*/
    int reconisim;
    curcell *fitNW;/**< DM fitting low rank terms*/
    cuspcell *actslave;/**<DM fitting actuator slaving*/
    cumoao_t *moao;/**<moao configurations for GPU*/
    curcell **dm_wfs;/**<moao results for wfs for warm restart*/
    curcell **dm_evl;/**<moao results for evl for warm restart*/
    stream_t *moao_stream;
    //CBS Tomo
    cusp *RCl;/**<Tomography Cholesky factor*/
    int  *RCp;/**<Tomography Cholesky permutation vector*/
    curmat *RUp;
    curmat *RVp;
    curmat *RMI;//SVD
    
    cusp *FCl;
    int  *FCp;
    curmat *FUp;
    curmat *FVp;
    curmat *FMI;//SVD
    int disablelrt; /*1: disable t/t removal lrt in split tomo*/

    curcell *RFlgsx;
    curcell *RFngsx;
    CGTMP_T cgtmp_tomo;
    CGTMP_T cgtmp_fit;
    CGTMP_T cgtmp_moaowfs;
    CGTMP_T cgtmp_moaoevl;
    /*the following data reside in the gpu memory*/
    GPU_PROP_GRID_T *hxdata;
    GPU_PROP_GRID_T *hxtdata;
    GPU_GP_T *gpdata;
    GPU_FDPCG_T *fddata;
    curecon_t(){ /*This is not good. zeros out already initialized childs.*/
	memset(this, 0, sizeof(*this));
    }
}curecon_t;

__global__ void apply_W_do(float *restrict out, const float *restrict in, const int *W0f, 
			   float alpha, int nx, int n);
W01_T *gpu_get_W01(dsp *R_W0, dmat *R_W1);


void gpu_TomoR(curcell **xout, float beta, const void *A, const curcell *grad, float alpha, stream_t &stream);
void gpu_TomoRt(curcell **gout,float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_TomoL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_FitR (curcell **xout, float beta, const void *A, const curcell *xin, float alpha);
void gpu_FitRt(curcell **xout, float beta, const void *A, const curcell *xin, float alpha);
void gpu_FitL (curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_Tomo_fdprecond(curcell **xout, const void *A, const curcell *xin, cudaStream_t stream);

void cumuv(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha);
void cumuv_trans(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha);

void cuchol_solve(float *restrict out, cusp *Cl, int *Cp, const float *restrict in, 
		  cudaStream_t stream);

void gpu_tomo_test(SIM_T *simu);
void gpu_fit_test(SIM_T *simu);
#endif
