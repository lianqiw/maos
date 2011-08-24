#ifndef AOS_CUDA_RECON_H
#define AOS_CUDA_RECON_H
#include "curmat.h"
typedef struct{
    curcell *gradin; /**< The grad to operator on*/
    curcell *neai;
    curcell *opdwfs;/**<Temporary*/
    curcell *opdx;  /**<Temporary*/
    curcell *grad;  /**<Temporary*/
    cudaStream_t *wfsstream;
    cublasHandle_t *wfshandle;
    cusparseHandle_t *wfssphandle;
    cudaStream_t *psstream;
    curcell *PTT;  /**< Global tip/tilt, DF removal */
    curcell *PDF;  /**< Differential focus removal */
    float *l2c;    /**< Laplacian */
    int *zzi;      /**< Piston constraint coordinate*/
    float *zzv;    /**< Piston constraint value*/
}curecon_t;

extern curecon_t *curecon;


#endif
