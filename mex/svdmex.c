#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_A,
	P_TOT,
    };
    enum{
	PL_U,
	PL_S,
	PL_V,
	PL_TOT,
    };
  if(nlhs!=PL_TOT || nrhs !=P_TOT){
      mexErrMsgTxt("Usage: [u,s,v]=svdmex(A)\n");
  }
#if USE_MKL
  omp_set_num_threads(&NCPU2)
#endif
  dmat *A=mx2d(prhs[P_A]);
  dmat *U=NULL, *Sdiag=NULL, *VT=NULL;
  dsvd(&U, &Sdiag, &VT, A);
  plhs[PL_U]=d2mx(U);
  plhs[PL_S]=d2mx(Sdiag);
  dmat *V=dtrans(VT);
  plhs[PL_V]=d2mx(V);
  dfree(U); dfree(Sdiag); dfree(VT); dfree(V); dfree(A); 
}
