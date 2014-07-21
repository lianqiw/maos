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
void KL_kolmogorov_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=3) mexErrMsgTxt("Expect 3 arguments\n");
    loc_t* loc=mx2loc(prhs[0]);
    double R=(double)mxGetScalar(prhs[1]);
    int nr=(int)mxGetScalar(prhs[2]);
    dmat* KL_kolmogorov_out=KL_kolmogorov(loc,R,nr);
    plhs[0]=any2mx(KL_kolmogorov_out);
}
void zernike_cov_kolmogorov_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=1) mexErrMsgTxt("Expect 1 arguments\n");
    int nr=(int)mxGetScalar(prhs[0]);
    dmat* zernike_cov_kolmogorov_out=zernike_cov_kolmogorov(nr);
    plhs[0]=any2mx(zernike_cov_kolmogorov_out);
}
void zernike_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=3) mexErrMsgTxt("Expect 3 arguments\n");
    loc_t* loc=mx2loc(prhs[0]);
    double R=(double)mxGetScalar(prhs[1]);
    int nr=(int)mxGetScalar(prhs[2]);
    dmat* zernike_out=zernike(loc,R,nr);
    plhs[0]=any2mx(zernike_out);
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
void diag_mod_cov_mex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs!=2) mexErrMsgTxt("Expect 2 arguments\n");
    dmat* mz=mx2d(prhs[0]);
    dmat* cov=mx2d(prhs[1]);
    dmat* diag_mod_cov_out=diag_mod_cov(mz,cov);
    plhs[0]=any2mx(diag_mod_cov_out);
    dfree(mz);
    dfree(cov);
}
void print_usage(){
    printf("Usage:\n");
    printf("out=aolib('dtrapz',x,y)\n");
    printf("out=aolib('KL_kolmogorov',loc,R,nr)\n");
    printf("out=aolib('zernike_cov_kolmogorov',nr)\n");
    printf("out=aolib('zernike',loc,R,nr)\n");
    printf("out=aolib('psdinterp1',psdin,fnew)\n");
    printf("out=aolib('diag_mod_cov',mz,cov)\n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if(nrhs<1){
	print_usage();
        return;
    }
    char *cmd=mxArrayToString(prhs[0]);
    if(!strcmp(cmd, "dtrapz")) dtrapz_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "KL_kolmogorov")) KL_kolmogorov_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "zernike_cov_kolmogorov")) zernike_cov_kolmogorov_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "zernike")) zernike_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "psdinterp1")) psdinterp1_mex(nlhs, plhs, nrhs-1, prhs+1);
    else if(!strcmp(cmd, "diag_mod_cov")) diag_mod_cov_mex(nlhs, plhs, nrhs-1, prhs+1);
    else print_usage();
}

