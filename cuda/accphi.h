#ifndef AOS_CUDA_ACCPHI_H
#define AOS_CUDA_ACCPHI_H

void gpu_atm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax,const float thetay,
		 const float mispx, const float mispy, const float dtisim, const float atmalpha, cudaStream_t stream);
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax, const float thetay,
		const float mispx, const float mispy, const float dmalpha, cudaStream_t stream);
#endif

