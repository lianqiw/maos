#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "io.h"
mxArray *INVALID=(mxArray*)1;
static mxArray *readdata(file_t *fp, mxArray **header, int start, int howmany){
    /*
      if start!=0 || howmany!=0 and data is cell, will only read cell from start to start+howmany
      if data is not cell and start==-1, will skip the data.
      Only the first call to readdata will possibly have howmany!=0
     */
    uint32_t magic;
    mxArray *out=NULL;
    long nx,ny;
    if(fp->eof) return NULL;
    header_t header2;
    if(read_header2(&header2, fp)){
	return NULL;
    }
    magic=header2.magic;
    if(header){
	if(header2.str)
	    *header=mxCreateString(header2.str);
	else
	    *header=mxCreateString("");
    }
    free(header2.str); header2.str=NULL;
    nx=header2.nx;
    ny=header2.ny;
    int start_save=start;
    if(iscell(magic)){
	if(iscell(magic) && start!=-1 && howmany==0){
	    if(start!=0){
		error("Invalid use\n");
	    }
	    howmany=nx*ny;
	}
    }else if(fp->isfits){
	if(howmany!=0){
	    /*first read of fits file. determine if we need it*/
	    if(start>0){
		start=-1;//skip this block.
	    }
	}
    }else{
	if((start!=0 && start!=-1) || howmany!=0){
	    error("invalid use");
	}
    }
    int iscell=0;
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
	    iscell=1;
	    mwIndex ix;
	    if(fp->eof) return NULL;
	    out=mxCreateCellMatrix(nx,ny);
	    mxArray *header0=mxCreateCellMatrix(nx*ny+1,1);
	    for(ix=0; ix<nx*ny; ix++){
		int start2=0;
		if(start==-1 || ix<start || ix>=start+howmany){
		    start2=-1;
		}
		mxArray *header3=NULL;
		mxArray *tmp=readdata(fp, &header3, start2, 0);
		if(fp->eof){
		    break;
		}
		if(tmp>0){
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
	    if(start==-1) error("Invalid use\n");
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
	    if(start==-1) error("Invalid use\n");
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
			error("Invalid sparse format\n");
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
	if(start==-1){
	    if(zfseek(fp, sizeof(double)*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
	    out=mxCreateDoubleMatrix(nx,ny,mxREAL);
	    if(nx!=0 && ny!=0){
		zfread(mxGetPr(out), sizeof(double),nx*ny,fp);
	    }
	}
	break;
    case M_FLT:/*float array. convert to double*/
	if(start==-1){
	    if(zfseek(fp, sizeof(float)*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
	    out=mxCreateDoubleMatrix(nx,ny,mxREAL);
	    if(nx!=0 && ny!=0){
		float *tmp=malloc(nx*ny*sizeof(float));
		zfread(tmp, sizeof(float),nx*ny, fp);
		double *p=mxGetPr(out);
		long i;
		for(i=0; i<nx*ny; i++){
		    p[i]=tmp[i];
		}
		free(tmp);
	    }
	}break;
    case M_INT64:/*long array*/
	if(start==-1){
	    if(zfseek(fp, 8*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
	    out=mxCreateNumericMatrix(nx,ny,mxINT64_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		/*Don't use sizeof(mxINT64_CLASS), it is just an integer, not a valid C type.*/
		zfread(mxGetPr(out), 8,nx*ny,fp);
	    }
	}
	break;
    case M_INT32:
	if(start==-1){
	    if(zfseek(fp, 4*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
	    out=mxCreateNumericMatrix(nx,ny,mxINT32_CLASS,mxREAL);
	    if(nx!=0 && ny!=0){
		zfread(mxGetPr(out), 4,nx*ny,fp);
	    }
	}
	break;
    case M_CMP:/*double complex array*/
	if(start==-1){
	    if(zfseek(fp, 16*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
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
    case M_ZMP:/*float complex array. convert to double*/
	if(start==-1){
	    if(zfseek(fp, 8*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=INVALID;
	}else{
	    out=mxCreateDoubleMatrix(nx,ny,mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		fcomplex*tmp=malloc(sizeof(fcomplex)*nx*ny);
		zfread(tmp,sizeof(fcomplex),nx*ny,fp);
		double *Pr=mxGetPr(out);
		double *Pi=mxGetPi(out);
		long i;
		for(i=0; i<nx*ny; i++){
		    Pr[i]=(double)crealf(tmp[i]);
		    Pi[i]=(double)cimagf(tmp[i]);
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
    start=start_save;
    if(!iscell && fp->isfits==1){/*fits file may contain extra extensions.*/
	fp->isfits++;
	int nx=0;
	mxArray **outarr=NULL;
	mxArray **headerarr=NULL;
	while(out){
	    nx++;
	    outarr=realloc(outarr, sizeof(mxArray*)*nx);
	    headerarr=realloc(headerarr, sizeof(mxArray*)*nx);
	    if(out==INVALID) out=NULL;
	    outarr[nx-1]=out;
	    if(header) headerarr[nx-1]=*header;
	    int start2=0;
	    if(howmany!=0){//selective reading.
		if(nx<start || nx+1>start+howmany){
		    start2=-1;//don't read next data.
		}
	    }
	    out=readdata(fp, header, start2, 0);
	}
	if(nx>1){/*set output.*/
	    out=mxCreateCellMatrix(nx, 1);
	    if(header) *header=mxCreateCellMatrix(nx, 1);
	    int i;
	    for(i=0; i<nx; i++){
		mxSetCell(out, i, outarr[i]);
		if(header) mxSetCell(*header, i, headerarr[i]);
	    }
	    free(outarr);
	    free(headerarr);
	}else{
	    out=outarr[0];
	    if(header) *header=headerarr[0];
	}
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
    char *fn=NULL;
    int start=0, howmany=0;
    switch(nrhs){
    case 3:
	howmany=(long)mxGetScalar(prhs[2]);//do not break
    case 2:
	start=(long)mxGetScalar(prhs[1]);//starting block to read. matlab index.
    case 1:
	fn=mx2str(prhs[0]);
	break;
    default:
	mexErrMsgTxt("Usage: [var, [header]]=read('filename' [,start] [,howmany]). [] means optional\n");
    }
    if(howmany>0){
	if(start>0){
	    start--;//convert to C index.
	}
	if(start<0){
	    start=0;
	}
    }
    fp=zfopen(fn,"rb");
    if(!fp){
	error("Unable to open file: %s\n", fn);
	return;
    }
    free(fn);
    switch(nlhs){
    case 2:
	plhs[0]=readdata(fp, &plhs[1], start, howmany); break;
    case 1:
	plhs[0]=readdata(fp, NULL, start, howmany); break;
    default:
	mexErrMsgTxt("Usage: [var, [header]]=read('filename' [,start] [,howmany]). [] means optional\n");
    }
    if(start==0 && howmany==0){
	zfeof(fp);
    }
    zfclose(fp);
}
