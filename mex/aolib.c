#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE
#endif
#include "interface.h"

void dtrapz_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=2) mexErrMsgTxt("Expect 2 arguments\n");
    dmat* x=mx2d(prhs[0]);
    dmat* y=mx2d(prhs[1]);
    dmat* dtrapz_out=dtrapz(x,y);
    plhs[0]=any2mx(dtrapz_out);
    dfree(x);
    dfree(y);
}
void psdinterp1_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=2) mexErrMsgTxt("Expect 2 arguments\n");
    dmat* psdin=mx2d(prhs[0]);
    dmat* fnew=mx2d(prhs[1]);
    dmat* psdinterp1_out=psdinterp1(psdin,fnew);
    plhs[0]=any2mx(psdinterp1_out);
    dfree(psdin);
    dfree(fnew);
}
void print_usage(){
    printf("Usage:\n");
    printf("out=aolib('dtrapz',x,y)\n");
    printf("out=aolib('psdinterp1',psdin,fnew)\n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs<1){
	print_usage();
        return;
    }
    char *cmd=mxArrayToString(prhs[0]);
    if(!strcmp(cmd, "dtrapz")) dtrapz_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "psdinterp1")) psdinterp1_mex(nlhs, plhs, nrhs-1, prhs+1);
    else print_usage();
}

