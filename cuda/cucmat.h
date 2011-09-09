#ifndef AOS_CUDA_CUCMAT_H
#define AOS_CUDA_CUCMAT_H

#include "utils.h"
#include "types.h"
cucmat *cucnew(int nx, int ny);
cucmat *cucnew(int nx, int ny, cudaStream_t stream);
void cucfree(cucmat *A);
cuccell* cuccellnew(int nx, int ny);
void curcellfree(curcell *A);
void cucwrite(const cucmat *A, const char *format, ...);

inline void cuczero(cucmat *A, cudaStream_t stream){
    if(A && A->p){
	DO(cudaMemsetAsync(A->p, 0, A->nx*A->ny*sizeof(fcomplex), stream));
    }
}
#endif
