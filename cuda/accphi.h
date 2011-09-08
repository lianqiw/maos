#ifndef AOS_CUDA_ACCPHI_H
#define AOS_CUDA_ACCPHI_H
#include "curmat.h"

void gpu_atm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax,const float thetay,
		 const float mispx, const float mispy, const float dtisim, const float atmalpha, cudaStream_t stream);
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, cumap_t *cudm,
		const float hs, const float thetax, const float thetay,
		const float mispx, const float mispy, const float dmalpha, cudaStream_t stream);

void prop_grid_match(curmat *out, float oxo, float oyo,
		     const curmat *in, float oxi, float oyi, float dxi,
		     float dispx, float dispy,
		     float alpha, cudaStream_t stream);
void gpu_prop_grid(curmat *out, float oxo, float oyo, float dxo,
		   curmat *in, float oxi, float oyi, float dxi,
		   float dispx, float dispy,
		   float alpha, char trans, cudaStream_t stream);
#endif

