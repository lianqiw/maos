#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "io.h"

static mxArray *readdata(file_t *fp, mxArray **header){
    uint32_t magic;
    mxArray *out=NULL;
    uint64_t nx,ny;
    if(fp->eof) return NULL;
    char *header2=NULL;
    magic=read_magic(fp, &header2);

    if(header){
	if(header2)
	    *header=mxCreateString(header2);
	else
	    *header=mxCreateString("");
    }
    free(header2); header2=NULL;

    zfread(&nx, sizeof(uint64_t), 1,fp);
    zfread(&ny, sizeof(uint64_t), 1,fp);
    if(fp->eof) return NULL;
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
	    mwIndex ix;
	    if(fp->eof) return NULL;
	    out=mxCreateCellMatrix(nx,ny);
	    mxArray *header0=mxCreateCellMatrix(nx*ny+1,1);
	    for(ix=0; ix<nx*ny; ix++){
		mxArray *header3=NULL;
		mxArray *tmp=readdata(fp, &header3);
		if(fp->eof){
		    break;
		}
		if(tmp){
		    mxSetCell(out, ix, tmp);
		}
		if(header3){
		    mxSetCell(header0, ix, header3);
		}
	    }
	    if(header){
		mxSetCell(header0, nx*ny, *header);
		*header=header0;
	    }
	}
	break;
    case M_SP64:
    case M_SP32:
	{
	    size_t size;
	    if(magic==M_SP32){
		size=4;
	    }else if(magic==M_SP64){
		size=8;
	    }else{
		error("Invalid magic\n");
	    }
	    uint64_t nzmax;
	    if(nx!=0 && ny!=0){
		zfread(&nzmax,sizeof(uint64_t),1,fp);
	    }else{
		nzmax=0;
	    }
	    if(fp->eof) return NULL;
	    out=mxCreateSparse(nx,ny,nzmax,mxREAL);
	    if(nx!=0 && ny!=0 && nzmax!=0){
		if(sizeof(mwIndex)==size){/*Match*/
		    warning("Matching %ld byte\n", size);
		    zfread(mxGetJc(out), size,ny+1,fp);
		    zfread(mxGetIr(out), size,nzmax, fp);
		}else{
		    long i;
		    mwIndex *Jc0=mxGetJc(out);
		    mwIndex *Ir0=mxGetIr(out);
		    void *Jc=malloc(size*(ny+1));
		    void *Ir=malloc(size*nzmax);
		    zfread(Jc, size, ny+1, fp);
		    zfread(Ir, size, nzmax, fp);
		    if(size==4){
			warning("Converting 4 byte to 8 byte\n");
			uint32_t* Jc2=Jc;
			uint32_t* Ir2=Ir;
			for(i=0; i<ny+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else if(size==8){
			warning("Converting 8 byte to 4 byte\n");
			uint64_t* Jc2=Jc;
			uint64_t* Ir2=Ir;
			for(i=0; i<ny+1; i++){
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
		zfread(mxGetPr(out), sizeof(double), nzmax, fp);
	    }
	}
	break;
    case M_CSP64:
    case M_CSP32:/*complex sparse*/
	{
	    size_t size;
	    if(magic==M_CSP32){
		size=4;
	    }else if(magic==M_CSP64){
		size=8;
	    }
	    uint64_t nzmax;
	    if(nx!=0 && ny!=0){
		zfread(&nzmax,sizeof(uint64_t),1,fp);
	    }else{
		nzmax=0;
	    }
	    if(fp->eof) return NULL;
	    out=mxCreateSparse(nx,ny,nzmax,mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		long i;
		if(sizeof(mwIndex)==size){
		    zfread(mxGetJc(out), size,ny+1,fp);
		    zfread(mxGetIr(out), size,nzmax, fp);
		}else{
		    mwIndex *Jc0=mxGetJc(out);
		    mwIndex *Ir0=mxGetIr(out);
		    void *Jc=malloc(size*(ny+1));
		    void *Ir=malloc(size*nzmax);
		    zfread(Jc, size, ny+1, fp);
		    zfread(Ir, size, nzmax, fp);
		    if(size==4){
			uint32_t* Jc2=Jc;
			uint32_t* Ir2=Ir;
			for(i=0; i<ny+1; i++){
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
			for(i=0; i<ny+1; i++){
			    Jc0[i]=Jc2[i];
			}
			for(i=0; i<nzmax; i++){
			    Ir0[i]=Ir2[i];
			}
			free(Jc);
			free(Ir);
		    }else{
			info("size=%lu\n", size);
			mexErrMsgTxt("Invalid sparse format\n");
		    }
		}
		dcomplex *tmp=malloc(sizeof(dcomplex)*nzmax);
		zfread(tmp, sizeof(dcomplex), nzmax, fp);
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
	    out=mxCreateDoubleMatrix(nx,ny,mxREAL);
	    if(nx!=0 && ny!=0){
		zfread(mxGetPr(out), sizeof(double),nx*ny,fp);
	    }
	}
	break;
    case M_INT64:/*long array*/
	{
	    out=mxCreateNumericMatrix(nx,ny,mxINT64_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		/*Don't use sizeof(mxINT64_CLASS), it is just an integer, not a valid C type.*/
		zfread(mxGetPr(out), 8,nx*ny,fp);
	    }
	}
	break;
    case M_INT32:
	{
	    out=mxCreateNumericMatrix(nx,ny,mxINT32_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		zfread(mxGetPr(out), 4,nx*ny,fp);
	    }
	}
	break;
    case M_CMP:/*double complex array*/
	{
	    out=mxCreateDoubleMatrix(nx,ny,mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		dcomplex*tmp=malloc(sizeof(dcomplex)*nx*ny);
		zfread(tmp,sizeof(dcomplex),nx*ny,fp);
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
    case M_HEADER:
	break;
    default:
	fprintf(stderr,"magic=%x\n",magic);
	warning("Unrecognized file. Please recompile the mex routines in the newest code\n");
	out=NULL;
    }
    if(!out){
	out=mxCreateDoubleMatrix(0,0,mxREAL);
    }

    return out;
}
static char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    file_t *fp;
    if(nrhs!=1 || nlhs>2){
	mexErrMsgTxt("Usage:var=read('filename') or [var, header]=read('file name')\n");
    }
    char *fn=mx2str(prhs[0]);
    fp=openfile(fn,"rb");
    free(fn);
    if(nlhs==2){
	plhs[0]=readdata(fp, &plhs[1]);
    }else{
	plhs[0]=readdata(fp, NULL);
	test_eof(fp);
    }
    zfclose(fp);
}
