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
#ifndef AOS_CUDA_WFS_H
#define AOS_CUDA_WFS_H
#include <cufft.h>

#define RAND_BLOCK 8
#define RAND_THREAD 256
struct cudtf_t{
    cucmat nominal;/*array for each sa. */
    cucmat etf;
    cucmat etf2;
    int etfis1d;
};
struct cullt_t{
    cupts_t pts;
    culoc_t loc;
};
typedef struct cupowfs_t{
    cupts_t pts;   /**<location of lower left OPD point in each sa*/
    culoc_t loc;  /**<location of OPD points. for reconstruction purpose only.*/
    culoc_t saloc;
    cuarray<culoc_t> msaloc;/**<Mishaped saloc, for pywfs.*/
    int **embed;       /**<embed for field stop computation*/
    int *nembed;       /**<embed for field stop computation*/
    curcell fieldstop;/**<*mask for field stop computation*/
    cullt_t llt;
    /*For Pyramid*/
    struct PYWFS_T *pywfs;//points to powfs[ipowfs].pywfs
    cuccell pyramid;
    cucmat pynominal;
    curmat saa;
    curmat pyoff; //pywfs->gradoff
    cupowfs_t():embed(0),nembed(0),pywfs(0){};
}cupowfs_t;
/**For matched filter update*/
class dither_t{
    int      imc;
    curcell imb;
    curcell imx;
    curcell imy;
public:
    dither_t():imc(0){}
    dither_t(int nsa, int pixpsax, int pixpsay);
    void acc(DITHER_T *dither, curcell &ints, Real cs, Real ss, int nstat, cudaStream_t stream);
};
class cuwfs_t{//one for each WFS.
public:
    stream_t stream;
    cupowfs_t *powfs;
    cuarray<culoc_t> loc_dm;  /**<Grid for ray tracing from DM to WFS*/
    culoc_t loc_tel;  /**<Grid for ray tracing from Telescope to WFS*/
    cusp GS0;         /**<For gtilt. is GS0t in col major */
    curmat imcc;     /**<size is 9*nsa*/
    curmat neasim;     /**<The noise equivalent angles for each subaperture.*/
    curmat  amp;        /**<Amplitude map*/
    cufftHandle plan1, plan2, plan3,plan_fs;   /**<FFTW plan if any*/
    cuarray<cudtf_t> dtf;       /**<array for each wvl.*/
    curmat srot;      /**<angle to rotate PSF/OTF*/
    curmat mtche;     /**<matched filter gradient operator.*/
    curmat i0sum;     /**<sum of i0 for each subaperture.*/
    curmat bkgrnd2;   /**<background as an image*/
    curmat bkgrnd2c;  /**<calibration of background to subtract.*/
    curmat cogcoeff;
    /*For LLT */
    curmat lltncpa;    /**<NCPA for llt*/
    curmat lltimcc;      /**<size of 9x1*/
    curmat lltamp;
    int msa;            /**<Number of subapertures in each batch of FFT. <nsa to save memory in psf.*/
    cufftHandle lltplan_wvf, lltplan_otf;/**<FFTW plan for LLT*/
    curmat opdadd;     /**<The ncpa and surface aberration.*/
    /*For random number of this wfs. */
    struct curandStateXORWOW *custat;
    int     custatb;/*allocated block */
    int     custatt;/*allocated thread */
    /*Run time data that changes */
    curmat phiout;
    curmat gradacc;    /**<For accumulating grads*/
    curmat gradcalc;   /**<For outputing grads*/
    curmat lltopd;
    Real  *lltg; //pinned memory.
    cucmat lltwvf;
    cucmat lltotfc;
    cucmat wvf;
    cucmat psf;
    cucmat otf;
    curcell ints;       /**<For accumulating subaperture image.*/
    curcell pistatout;  /**<For output pistatout*/
    cuccell wvfout;
    cucmat psfout;
    cucmat psfstat;
    /*For matched filter update*/
    dither_t dither;
    /*For Pyramid WFS*/
    cuccell pywvf;//Original PSF from OPD
    cucmat pyotf; //Truncated OTF to be multiplied by pyramid and FFT.
    curmat pypsf; //stores psf during modulation.
    curmat isum;//stores sum of psf for all subapertures at each frame in GPU.
    cufftHandle plan_py;
    cuwfs_t():powfs(0),dtf(0),msa(0),custatb(0),custatt(0),lltg(0){
    }
};

void gpu_wfsints(SIM_T *simu, Real *phiout, curmat &gradref, int iwfs, int isim, cudaStream_t stream);

void cuztilt(Real *restrict g, Real *restrict opd, 
	     const int nsa, const Real dx, const int nx, Real *imcc,
	     const Real (*orig)[2], const Real*restrict amp, Real alpha, cudaStream_t stream);
__global__ void cpcenter_do(Comp *restrict out, int noutx, int nouty,
			    const Comp *restrict in, int ninx, int niny);
void pywfs_grad(curmat &grad, const curmat &ints, const curmat &saa, curmat &isum, const curmat &goff, Real gain, cudaStream_t stream);
void pywfs_ints(curmat &ints, curmat &phiout, cuwfs_t &cuwfs, Real siglev, cudaStream_t stream);
dsp *gpu_pywfs_mkg(const PARMS_T *parms, const POWFS_T *powfs, loc_t *aloc, int iwfs, int idm);
#endif
