#ifndef AOS_CUDA_WFS_H
#define AOS_CUDA_WFS_H
typedef curandState_t curandStat;


#define RAND_BLOCK 16
#define RAND_THREAD 32
typedef struct {
    fcomplex **nominal;//array for each sa.
    //cusp_t **si;
    fcomplex **etf;
    int etfis1d;
}cudtf_t;

typedef struct cuwloc_t{
    float (*saorig)[2];
    float dx;
    float (*loc)[2];
    int nloc;
    int nsa;
    int nxsa;
    struct cuwloc_t *llt;
}cuwloc_t;


typedef struct{
    cusparseHandle_t sphandle;
    cudaStream_t stream;
    cusp_t *GS0;        /**<For gtilt*/
    float (**imcc)[3];  /**<For ztilt.*/
    float  *gradacc;    /**<For accumulating grads*/
    float  *ints;       /**<For accumulating subaperture image.*/
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
    cufftHandle lltplan1, lltplan2;/**<FFTW plan for LLT*/
    //For random number of this wfs.
    curandStat *custat;
    int     custatb;//allocated block
    int     custatt;//allocated thread
    float  *neareal;
}cuwfs_t;

extern cuwloc_t *cupowfs;
extern cuwfs_t *cuwfs;
extern cusparseMatDescr_t cuspdesc;
void wfsints(SIM_T *simu, float *phiout, int iwfs, int isim, cudaStream_t stream);

__global__ void cuztilt(float *restrict g, float *restrict opd, 
			const int nsa, const float dx, const int nx, float (**imcc)[3],
			const float (*orig)[2], const float*restrict amp, float alpha);
__global__ void cpcenter_do(fcomplex *restrict out, int noutx, int nouty,
			    const fcomplex *restrict in, int ninx, int niny);
#endif
