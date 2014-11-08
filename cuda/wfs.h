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
#include <cufft.h>

#define RAND_BLOCK 8
#define RAND_THREAD 256
struct cudtf_t{
    cucmat *nominal;/*array for each sa. */
    cucmat *etf;
    cucmat *etf2;
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
    curcell *imb;
    curcell *imx;
    curcell *imy;
public:
    dither_t(int nsa, int pixpsax, int pixpsay);
    void acc(DITHER_T *dither, curcell *ints, Real angle, int nstat, cudaStream_t stream);
    ~dither_t(){
	delete imb;
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
    Real (**imcc)[3];  /**<For ztilt.*/
    Real  *neasim;     /**<The noise equivalent angles for each subaperture.*/
    Real  *amp;        /**<Amplitude map*/
    cufftHandle plan1, plan2, plan3,plan_fs;   /**<FFTW plan if any*/
    cudtf_t *dtf;       /**<array for each wvl.*/
    Real   *srot;      /**<angle to rotate PSF/OTF*/
    curmat *mtche;     /**<matched filter gradient operator.*/
    curmat *i0sum;     /**<sum of i0 for each subaperture.*/
    Real  **bkgrnd2;   /**<background as an image*/
    Real  **bkgrnd2c;  /**<calibration of background to subtract.*/
    Real *cogcoeff;
    /*For LLT */
    curmat *lltncpa;    /**<NCPA for llt*/
    Real (**lltimcc)[3];
    Real  *lltamp;
    int msa;            /**<Number of subapertures in each batch of FFT. <nsa to save memory in psf.*/
    cufftHandle lltplan_wvf, lltplan_otf;/**<FFTW plan for LLT*/
    curmat *opdadd;     /**<The ncpa and surface aberration.*/
    /*For random number of this wfs. */
    struct curandStateXORWOW *custat;
    int     custatb;/*allocated block */
    int     custatt;/*allocated thread */
    /*Run time data that changes */
    curmat *phiout;
    curmat *gradacc;    /**<For accumulating grads*/
    curmat *gradcalc;   /**<For outputing grads*/
    curmat *lltopd;
    Real  *lltg;
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

void gpu_wfsints(SIM_T *simu, Real *phiout, curmat *gradref, int iwfs, int isim, cudaStream_t stream);

void cuztilt(Real *restrict g, Real *restrict opd, 
	     const int nsa, const Real dx, const int nx, Real (**imcc)[3],
	     const Real (*orig)[2], const Real*restrict amp, Real alpha, cudaStream_t stream);
__global__ void cpcenter_do(Comp *restrict out, int noutx, int nouty,
			    const Comp *restrict in, int ninx, int niny);
#endif
