#ifndef AOS_CUDA_PCG_H
#define AOS_CUDA_PCG_H
typedef void (*G_CGFUN)(curcell**, const void*, const curcell*, float);
typedef void (*G_PREFUN)(curcell**, const void*, const curcell*);

void gpu_pcg(curcell **px, 
	     G_CGFUN Amul, const void *A, 
	     G_PREFUN Mmul, const void *M, 
	     const curcell *b, int warm, int maxiter,
	     cudaStream_t stream);
#endif
