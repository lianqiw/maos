/*
  Create a random stream with input seed
*/
#include "random.h"
#include <mex.h>
#include <math.h>
#define USE_PTHREAD 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned int seed,nlem;
    struct_rand *p;
    if(nrhs!=1 || nlhs!=1) 
	mexErrMsgTxt("Accept one input one output\n");
    seed=(unsigned int)mxGetScalar(prhs[0]);
    nlem=(unsigned int)ceil((double)sizeof(struct_rand)
			    /(double)sizeof(mxINT32_CLASS));
    plhs[0]=mxCreateNumericMatrix(nlem,1,mxINT32_CLASS,mxREAL);
    p=(struct_rand*)mxGetPr(plhs[0]);
    fprintf(stderr,"Seed=%u\n",seed);
    seed_rand(p,seed);
}
