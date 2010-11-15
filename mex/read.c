#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_FN_LEN 800
/*compile with 
  mex read.c -largeArrayDims*/
#include "io.h"
#include <complex.h>

static mxArray *readdata(file_t *fp){
    int magic;
    mxArray *out=NULL;
    if(fp->eof) return NULL;
    readfile(&magic, sizeof(int),1,fp);
    if(fp->eof){
	return NULL;
    }
    switch(magic){
    case MCC_ANY:
    case MCC_DBL:
    case MCC_CMP:
    case MC_CSP:
    case MC_SP:
    case MC_DBL:
    case MC_CMP:
    case MC_INT32:
    case MC_INT64:
	{
	    mwIndex nx,ny;
	    mwIndex ix;
	    readfile(&nx, sizeof(mwIndex), 1,fp);
	    readfile(&ny, sizeof(mwIndex), 1,fp);
	    if(fp->eof) return NULL;
	    out=mxCreateCellMatrix(nx,ny);
	    for(ix=0; ix<nx*ny; ix++){
		mxArray *tmp=readdata(fp);
		if(fp->eof){
		    break;
		}
		if(tmp){
		    mxSetCell(out, ix, tmp);
		}
	    }
	}
	break;
    case M_SP64:
    case M_SP32:
	{
	    mwIndex m,n,nzmax;
	    size_t size;
	    if(magic==M_SP32){
		size=4;
	    }else if(magic==M_SP64){
		size=8;
	    }else{
		error("Invalid magic\n");
	    }
	    uint64_t m2,n2,nzmax2;
	    readfile(&m2,sizeof(uint64_t),1,fp);
	    readfile(&n2,sizeof(uint64_t),1,fp);
	    if(m2!=0 && n2!=0){
		readfile(&nzmax2,sizeof(uint64_t),1,fp);
	    }else{
		nzmax2=0;
	    }
	    m=m2;
	    n=n2;
	    nzmax=nzmax2;
	    if(fp->eof) return NULL;
	    out=mxCreateSparse(m,n,nzmax,mxREAL);
	    if(m!=0 && n!=0){
		if(sizeof(mwIndex)==size){
		    readfile(mxGetJc(out), sizeof(mwIndex),n+1,fp);
		    readfile(mxGetIr(out), sizeof(mwIndex),nzmax, fp);
		}else{
		    long i;
		    mwIndex *Jc0=mxGetJc(out);
		    mwIndex *Ir0=mxGetIr(out);
		    void *Jc=malloc(size*(n+1));
		    void *Ir=malloc(size*nzmax);
		    readfile(Jc, size, n+1, fp);
		    readfile(Ir, size, nzmax, fp);
		    warning("Converting from %zu to %zu bytes\n",size,sizeof(mwIndex));
		    if(size==4){
			uint32_t* Jc2=Jc;
			uint32_t* Ir2=Ir;
			for(i=0; i<n+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else if(size==8){
			uint64_t* Jc2=Jc;
			uint64_t* Ir2=Ir;
			for(i=0; i<n+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else{
			mexErrMsgTxt("Invalid sparse format\n");
		    }
		}
		readfile(mxGetPr(out), sizeof(double), nzmax, fp);
	    }
	}
	break;
    case M_CSP64:
    case M_CSP32:/*complex sparse*/
	{
	    mwIndex m,n,nzmax;
	    size_t size;
	    if(magic==M_SP32){
		size=4;
	    }else if(magic==M_SP64){
		size=8;
	    }
	    uint64_t m2,n2,nzmax2;
	    readfile(&m2,sizeof(uint64_t),1,fp);
	    readfile(&n2,sizeof(uint64_t),1,fp);
	    if(m2!=0 && n2!=0){
		readfile(&nzmax2,sizeof(uint64_t),1,fp);
	    }else{
		nzmax2=0;
	    }
	    m=m2;
	    n=n2;
	    nzmax=nzmax2;
	    if(fp->eof) return NULL;
	    out=mxCreateSparse(m,n,nzmax,mxCOMPLEX);
	    if(m!=0 && n!=0){
		long i;
		if(sizeof(mwIndex)==size){
		    readfile(mxGetJc(out), sizeof(mwIndex),n+1,fp);
		    readfile(mxGetIr(out), sizeof(mwIndex),nzmax, fp);
		}else{
		    mwIndex *Jc0=mxGetJc(out);
		    mwIndex *Ir0=mxGetIr(out);
		    void *Jc=malloc(size*(n+1));
		    void *Ir=malloc(size*nzmax);
		    readfile(Jc, size, n+1, fp);
		    readfile(Ir, size, nzmax, fp);
		    if(size==4){
			uint32_t* Jc2=Jc;
			uint32_t* Ir2=Ir;
			for(i=0; i<n+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else if(size==8){
			uint64_t* Jc2=Jc;
			uint64_t* Ir2=Ir;
			for(i=0; i<n+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else{
			mexErrMsgTxt("Invalid sparse format\n");
		    }
		}
		dcomplex *tmp=malloc(sizeof(dcomplex)*nzmax);
		readfile(tmp, sizeof(dcomplex), nzmax, fp);
		double *Pr=mxGetPr(out);
		double *Pi=mxGetPi(out);
		for(i=0; i<nzmax; i++){
		    Pr[i]=creal(tmp[i]);
		    Pi[i]=cimag(tmp[i]);
		}
		free(tmp);
	    }
	}
	break;
    case M_DBL:/*double array*/
	{
	    long nx,ny;
	    readfile(&nx, sizeof(long),1,fp);
	    readfile(&ny, sizeof(long),1,fp);
	    if(fp->eof) return NULL;
	    out=mxCreateDoubleMatrix(nx,ny,mxREAL);
	    if(nx!=0 && ny!=0){
		readfile(mxGetPr(out), sizeof(double),nx*ny,fp);
	    }
	}
	break;
    case M_INT64:/*long array*/
	{
	    long nx,ny;
	    readfile(&nx, sizeof(long),1,fp);
	    readfile(&ny, sizeof(long),1,fp);
	    if(fp->eof) return NULL;
	    out=mxCreateNumericMatrix(nx,ny,mxINT64_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		/*Don't use sizeof(mxINT64_CLASS), it is just an integer, not a valid C type.*/
		readfile(mxGetPr(out), 8,nx*ny,fp);
	    }
	}
	break;
    case M_INT32:
	{
	    long nx,ny;
	    readfile(&nx, sizeof(long),1,fp);
	    readfile(&ny, sizeof(long),1,fp);
	    if(fp->eof) return NULL;
	    out=mxCreateNumericMatrix(nx,ny,mxINT32_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		readfile(mxGetPr(out), 4,nx*ny,fp);
	    }
	}
	break;
    case M_CMP:/*double complex array*/
	{
	    long nx,ny;
	    readfile(&nx, sizeof(long),1,fp);
	    readfile(&ny, sizeof(long),1,fp);
	    if(fp->eof) return NULL;
	    out=mxCreateDoubleMatrix(nx,ny,mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		dcomplex*tmp=malloc(sizeof(dcomplex)*nx*ny);
		readfile(tmp,sizeof(dcomplex),nx*ny,fp);
		double *Pr=mxGetPr(out);
		double *Pi=mxGetPi(out);
		long i;
		for(i=0; i<nx*ny; i++){
		    Pr[i]=creal(tmp[i]);
		    Pi[i]=cimag(tmp[i]);
		}
		free(tmp);
	    }
	}
	break;
    default:
	fprintf(stderr,"magic=%x\n",magic);
	warning("Unrecognized file. Please recompile the mex routines in the newest code\n");
	out=NULL;
    }
    return out;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    file_t *fp;
    char fn[MAX_FN_LEN];
    mxGetString(prhs[0],fn,MAX_FN_LEN+1);
    (void)nlhs;
    (void)nrhs;
  
    fp=openfile(fn,"r");
    plhs[0]=readdata(fp);
    test_eof(fp);
    closefile(fp);
}
