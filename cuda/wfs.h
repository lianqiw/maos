#ifndef AOS_CUDA_WFS_H
#define AOS_CUDA_WFS_H
#include "cusparse.h"
#include "curand_kernel.h"
#include "cufft.h"
typedef curandState_t curandStat;

#define RAND_BLOCK 16
#define RAND_THREAD 32
typedef struct {
    fcomplex **nominal;//array for each sa.
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
    curmat *GPpx;     /**<GP for x grad in dense matrix format.*/
    curmat *GPpy;     /**<GP for y grad in dense matrix format.*/
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
    cufftHandle plan1, plan2;   /**<FFTW plan if any*/
    cudtf_t *dtf;       /**<array for each wvl.*/
    float   *srot;      /**<angle to rotate PSF/OTF*/
    float  (**mtche)[2]; /**<matched filter gradient operator.*/
    float  **bkgrnd2;   /**<background as an image*/
    float  **bkgrnd2c;  /**<calibration of background to subtract.*/
    //For LLT
    float  *lltncpa;    /**<NCPA for llt*/
    float (**lltimcc)[3];
    float  *lltamp;
    int msa;            /**<Number of subapertures in each batch of FFT. <nsa to save memory in psf.*/
    cufftHandle lltplan1, lltplan2;/**<FFTW plan for LLT*/
    curmat *opdadd;    /**<The ncpa and surface aberration.*/

    //Run time data that changes
    //For random number of this wfs.
    curandStat *custat;
    int     custatb;//allocated block
    int     custatt;//allocated thread
    float  *neareal;
    curmat *gradacc;    /**<For accumulating grads*/
    curmat *ints;       /**<For accumulating subaperture image.*/
    
}cuwfs_t;

void wfsints(SIM_T *simu, float *phiout, int iwfs, int isim, cudaStream_t stream);

__global__ void cuztilt(float *restrict g, float *restrict opd, 
			const int nsa, const float dx, const int nx, float (**imcc)[3],
			const float (*orig)[2], const float*restrict amp, float alpha);
__global__ void cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
			    const fcomplex *restrict in, int ninx, int niny);
#endif
