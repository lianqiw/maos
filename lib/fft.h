#ifndef AOS_LIB_FFT_H
#define AOS_LIB_FFT_H
#include "type.h"
void cfft2plan(cmat *A, int dir);
void cfft2(cmat *A, int dir);
void cifft2(cmat *A, int dir);
void cfft2s(cmat *A, int dir);
void cfft2partialplan(cmat *A, int ncomp, int dir);
void cfft2partial(cmat *A, int ncomp, int dir);
void cfree_plan(cmat *A);
cmat *cffttreat(cmat *A);
#endif
