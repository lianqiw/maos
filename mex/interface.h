#ifdef MATLAB_MEX_FILE
#include <mex.h>
#include "../lib/type.h"
#include "../lib/loc.h"
mxArray *d2mx(const dmat*A);
mxArray *c2mx(const cmat *A);
mxArray *dsp2mx(const dsp *A);
dmat *mx2d(const mxArray *A);

dsp *mx2dsp(const mxArray *A);
LOC_T *mx2loc(const mxArray *A);
#endif
