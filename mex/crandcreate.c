#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
/*
  Create a random stream with input seed
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned int seed,nlem;
    rand_t *p;
    if(nrhs!=1 || nlhs!=1){
	mexErrMsgTxt("Usage: stat=crandcreate(seed)");
    }
    seed=(unsigned int)mxGetScalar(prhs[0]);
    nlem=(unsigned int)ceil((double)sizeof(rand_t)
			    /(double)sizeof(mxINT32_CLASS));
    plhs[0]=mxCreateNumericMatrix(nlem,1,mxINT32_CLASS,mxREAL);
    p=(rand_t*)mxGetPr(plhs[0]);
    fprintf(stderr,"Seed=%u\n",seed);
    seed_rand(p,seed);
}
