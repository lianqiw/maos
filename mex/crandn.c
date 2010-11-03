#include "random.h"
#include <mex.h>
#include <math.h>
#define USE_PTHREAD 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    rand_t *strand;
    double *p;
    int nelem=1,M;
    int i;
    if(nrhs<1) 
	mexErrMsgTxt("Except one or more input\n");
    strand=(rand_t*)mxGetPr(prhs[0]);
    if(nrhs==1){
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    }else if(nrhs>=2){
	M=mxGetNumberOfElements(prhs[1]);
	p=mxGetPr(prhs[1]);
	if(M==2){
	    plhs[0]=mxCreateDoubleMatrix((int)p[0],(int)p[1],mxREAL);
	}else if(M==1){
	    if(nrhs==2){
		plhs[0]=mxCreateDoubleMatrix((int)p[0],1,mxREAL);
	    }else if(nrhs==3){
		double *p2;
		p2=mxGetPr(prhs[2]);
		plhs[0]=mxCreateDoubleMatrix((int)p[0],(int)p2[0],mxREAL);
	    }else{
		mexErrMsgTxt("Unexpected input: accepts at most 3 inputs");
	    }
	}else{
	    mexErrMsgTxt("Unexpected input");
	}
    }
    nelem=mxGetNumberOfElements(plhs[0]);
    p=mxGetPr(plhs[0]);
    for(i=0; i<nelem; i++){
	p[i]=randn(strand);
    }
}
