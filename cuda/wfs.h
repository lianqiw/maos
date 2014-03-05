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
#include <cufft.h>

#define RAND_BLOCK 16
#define RAND_THREAD 32
struct cudtf_t{
    fcomplex **nominal;/*array for each sa. */
    fcomplex **etf;
    int etfis1d;
};
struct cullt_t{
    cupts_t *pts;
    culoc_t *loc;
};
typedef struct cuwloc_t{
    cupts_t *pts;   /**<location of lower left OPD point in each sa*/
    culoc_t *loc;  /**<location of OPD points. for reconstruction purpose only.*/
    culoc_t *saloc;/**<Lower left corner of each sa. may be different by dx/2 than pts.*/
    int *embed;       /**<embed for field stop computation*/
    int nembed;       /**<embed for field stop computation*/
    curmat *fieldstop;/**<*mask for field stop computation*/
    cullt_t *llt;
}cupowfs_t;
/**For matched filter update*/
class dither_t{
    int      imc;
    curcell *im0;
    curcell *imx;
    curcell *imy;
public:
    dither_t(int nsa, int pixpsax, int pixpsay);
    void reset();
    void acc(curcell *ints, float angle, cudaStream_t stream);
    void output(float a2m, int iwfs, int isim, cudaStream_t stream);
    ~dither_t(){
	delete im0;
	delete imx;
	delete imy;
    }
};
class cuwfs_t{
  public:
    stream_t *stream;
    cupowfs_t *powfs;
    culoc_t **loc_dm;  /**<Grid for ray tracing from DM to WFS*/
    culoc_t *loc_tel;  /**<Grid for ray tracing from Telescope to WFS*/
    cusp *GS0;         /**<For gtilt. is GS0t in col major */
    float (**imcc)[3];  /**<For ztilt.*/
    float  *neasim;     /**<The noise equivalent angles for each subaperture.*/
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
    curmat *lltncpa;    /**<NCPA for llt*/
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
    curmat *phiout;
    curmat *gradacc;    /**<For accumulating grads*/
    curmat *gradcalc;   /**<For outputing grads*/
    curmat *lltopd;
    float  *lltg;
    cucmat *lltwvf;
    cucmat *lltotfc;
    cucmat *wvf;
    cucmat *psf;
    cucmat *otf;
    curcell *ints;       /**<For accumulating subaperture image.*/
    curcell *pistatout;  /**<For output pistatout*/
    cuccell *wvfout;
    cucmat *psfout;
    cucmat *psfstat;
    /*For matched filter update*/
    dither_t *dither;
};

void gpu_wfsints(SIM_T *simu, float *phiout, curmat *gradref, int iwfs, int isim, cudaStream_t stream);

void cuztilt(float *restrict g, float *restrict opd, 
	     const int nsa, const float dx, const int nx, float (**imcc)[3],
	     const float (*orig)[2], const float*restrict amp, float alpha, cudaStream_t stream);
__global__ void cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
			    const fcomplex *restrict in, int ninx, int niny);
#endif
