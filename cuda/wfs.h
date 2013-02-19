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
#ifndef AOS_CUDA_WFS_H
#define AOS_CUDA_WFS_H
#include <cusparse.h>
//#include <curand_kernel.h>
#include <cufft.h>

#define RAND_BLOCK 16
#define RAND_THREAD 32
typedef struct {
    fcomplex **nominal;/*array for each sa. */
    fcomplex **etf;
    int etfis1d;
}cudtf_t;

typedef struct cuwloc_t{
    float (*pts)[2];  /**<location of lower left OPD point in each sa*/
    float (*loc)[2];  /**<location of OPD points*/
    float dx;         /**<Sampling of OPD*/
    float dsa;        /**<Subaperture spacing*/
    int nloc;
    float (*saloc)[2];/**<Lower left corner of each sa. may be different by dx/2 than pts.*/
    int (*saptr)[2];  /**<pointer of subaperture in ploc*/
    int nsa;
    int nxsa;         /**<number of points in each subaperture in each dimension.*/
    cusp *GP;         /**<GP in col major*/
    float GPscale;    /**<Scale GP to fit in 2 byte int*/
    cumat<int> *GPp ; /**<GP for x/y grad in dense matrix format.*/
    int *embed;       /**<embed for field stop computation*/
    int nembed;       /**<embed for field stop computation*/
    curmat *fieldstop;/**<*mask for field stop computation*/
    struct cuwloc_t *llt;
}cuwloc_t;


typedef struct{
    cusparseHandle_t sphandle;
    cublasHandle_t handle;
    cudaStream_t stream;
    cuwloc_t *powfs;
    cusp *GS0t;         /**<For gtilt. is GS0t in col major */
    float (**imcc)[3];  /**<For ztilt.*/
    float  *neasim;     /**<The noise equivalent angles for each grad.*/
    float  *amp;        /**<Amplitude map*/
    cufftHandle plan1, plan2, plan3,plan_fs;   /**<FFTW plan if any*/
    cudtf_t *dtf;       /**<array for each wvl.*/
    float   *srot;      /**<angle to rotate PSF/OTF*/
    float  (**mtche)[2]; /**<matched filter gradient operator.*/
    float   *i0sum;     /**<sum of i0 for each subaperture.*/
    float  **bkgrnd2;   /**<background as an image*/
    float  **bkgrnd2c;  /**<calibration of background to subtract.*/
    float *cogcoeff;
    /*For LLT */
    float  *lltncpa;    /**<NCPA for llt*/
    float (**lltimcc)[3];
    float  *lltamp;
    int msa;            /**<Number of subapertures in each batch of FFT. <nsa to save memory in psf.*/
    cufftHandle lltplan_wvf, lltplan_otf;/**<FFTW plan for LLT*/
    curmat *opdadd;     /**<The ncpa and surface aberration.*/
    //curmat *gradphyoff; /**<The gradient offset for CoG*/
    curmat *gradoff;    /**<The gradient offset for ncpa_method=1.*/
    /*For random number of this wfs. */
    struct curandStateXORWOW *custat;
    int     custatb;/*allocated block */
    int     custatt;/*allocated thread */
    /*Run time data that changes */
    curmat *neareal;
    curmat *gradacc;    /**<For accumulating grads*/
    curcell *ints;       /**<For accumulating subaperture image.*/
    curcell *pistatout;  /**<For output pistatout*/
}cuwfs_t;

void gpu_wfsints(SIM_T *simu, float *phiout, curmat *gradref, int iwfs, int isim, cudaStream_t stream);

__global__ void cuztilt(float *restrict g, float *restrict opd, 
			const int nsa, const float dx, const int nx, float (**imcc)[3],
			const float (*orig)[2], const float*restrict amp, float alpha);
__global__ void cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
			    const fcomplex *restrict in, int ninx, int niny);
#endif
