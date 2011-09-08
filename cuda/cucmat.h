#ifndef AOS_CUDA_CUCMAT_H
#define AOS_CUDA_CUCMAT_H

#include "utils.h"
#include "types.h"
cucmat *cucnew(int nx, int ny);
void cucfree(cucmat *A);
cuccell* cuccellnew(int nx, int ny);
void curcellfree(curcell *A);
void cucwrite(const cucmat *A, const char *format, ...);

#endif
