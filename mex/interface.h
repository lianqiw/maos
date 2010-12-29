#ifdef MATLAB_MEX_FILE
#include <mex.h>
typedef struct dmat dmat;
typedef struct cmat cmat;
typedef struct dsp dsp;
typedef struct loc_t loc_t;
mxArray *d2mx(const dmat*A);
mxArray *c2mx(const cmat *A);
mxArray *dsp2mx(const dsp *A);
dmat *mx2d(const mxArray *A);

dsp *mx2dsp(const mxArray *A);
loc_t *mx2loc(const mxArray *A);
char *mx2str(const mxArray *A);
#endif
