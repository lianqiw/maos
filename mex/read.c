/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
mxArray *SKIPPED=(mxArray*)(1);
static mxArray *readdata(file_t *fp, mxArray **header, int start, int howmany){
    /*
      if start!=0 || howmany!=0 and data is cell, will only read cell from start to start+howmany
      if data is not cell and start==-1, will skip the data.
      Only the first call to readdata will possibly have howmany!=0
     */
    if(fp->eof) return NULL;
    header_t header2;
    if(read_header2(&header2, fp)){
	return NULL;
    }
    uint32_t magic=header2.magic & 0xFFFF;
    if(magic==0){//end of file or empty file
	fp->eof=1;
	return NULL;
    }
    if(header){
	if(header2.str)
	    *header=mxCreateString(header2.str);
	else
	    *header=mxCreateString("");
    }
    free(header2.str); header2.str=NULL;
    long nx=header2.nx;
    long ny=header2.ny;
    int start_save=start;
    if(iscell(magic)){
	if(start!=-1 && howmany==0){
	    if(start!=0){
		error("Invalid use\n");
	    }
	    howmany=nx*ny-start;
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
    mxArray *out=NULL;
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
		if(tmp && tmp!=SKIPPED){
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
		size=0;
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
	    switch(magic){
	    case M_CSP32:
		size=4;break;
	    case M_CSP64:
		size=8;break;
	    default:
		size=0;
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
		    Pr[i]=tmp[i].x;
		    Pi[i]=tmp[i].y;
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
	    out=SKIPPED;
	}else{
	    out=mxCreateDoubleMatrix(nx,ny,mxREAL);
	    if(nx!=0 && ny!=0){
		zfread(mxGetPr(out), sizeof(double),nx*ny,fp);
	    }
	}
	break;
    case M_FLT:/*float array*/
	if(start==-1){
	    if(zfseek(fp, sizeof(float)*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=SKIPPED;
	}else{
	    mwSize nxy[2]={nx,ny};
	    out=mxCreateNumericArray(2, nxy,mxSINGLE_CLASS, mxREAL);
	    if(nx!=0 && ny!=0){
		zfread((float*)mxGetPr(out), sizeof(float),nx*ny, fp);
	    }
	}break;
    case M_INT64:/*long array*/
    case M_INT32:
    case M_INT16:
    case M_INT8:
	{
	    int byte=0;
	    mxClassID id;
	    switch(magic){
	    case M_INT64: byte=8; id=mxINT64_CLASS; break;
	    case M_INT32: byte=4; id=mxINT32_CLASS; break;
	    case M_INT16: byte=2; id=mxINT16_CLASS; break;
	    case M_INT8:  byte=1; id=mxINT8_CLASS; break;
	    default: id=0;
	    }
	    if(start==-1){
		if(zfseek(fp, byte*nx*ny, SEEK_CUR)){
		    error("Seek failed\n");
		}
		out=SKIPPED;
	    }else{
		out=mxCreateNumericMatrix(nx,ny,id,mxREAL);
		if(nx!=0 && ny!=0){
		    /*Don't use sizeof(mxINT64_CLASS), it is just an integer, 
		      not a valid C type.*/
		    zfread(mxGetPr(out), byte,nx*ny,fp);
		}
	    }
	}
	break;
    case M_CMP:/*double complex array*/
	if(start==-1){
	    if(zfseek(fp, 16*nx*ny, SEEK_CUR)){
		error("Seek failed\n");
	    }
	    out=SKIPPED;
	}else{
	    out=mxCreateDoubleMatrix(nx,ny,mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		dcomplex*tmp=malloc(sizeof(dcomplex)*nx*ny);
		zfread(tmp,sizeof(dcomplex),nx*ny,fp);
		double *Pr=mxGetPr(out);
		double *Pi=mxGetPi(out);
		long i;
		for(i=0; i<nx*ny; i++){
		    Pr[i]=tmp[i].x;
		    Pi[i]=tmp[i].y;
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
	    out=SKIPPED;
	}else{
	    mwSize nxy[2]={nx, ny};
	    out=mxCreateNumericArray(2, nxy,mxSINGLE_CLASS, mxCOMPLEX);
	    if(nx!=0 && ny!=0){
		fcomplex*tmp=malloc(sizeof(fcomplex)*nx*ny);
		zfread(tmp,sizeof(fcomplex),nx*ny,fp);
		float *Pr=(float*)mxGetPr(out);
		float *Pi=(float*)mxGetPi(out);
		long i;
		for(i=0; i<nx*ny; i++){
		    Pr[i]=tmp[i].x;
		    Pi[i]=tmp[i].y;
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
	int icell=0;
	mxArray **outarr=NULL;
	mxArray **headerarr=NULL;
	while(out){
	    icell++;
	    outarr=realloc(outarr, sizeof(mxArray*)*icell);
	    headerarr=realloc(headerarr, sizeof(mxArray*)*icell);
	    if(out==SKIPPED) out=NULL;
	    outarr[icell-1]=out;
	    if(header) headerarr[icell-1]=*header;
	    int start2=0;
	    if(howmany!=0){//selective reading.
		if(icell<start || icell+1>start+howmany){
		    start2=-1;//don't read next data.
		}
	    }
	    out=readdata(fp, header, start2, 0);
	}
	if(icell>1){/*set output.*/
	    out=mxCreateCellMatrix(icell, 1);
	    if(header) *header=mxCreateCellMatrix(icell, 1);
	    int i;
	    for(i=0; i<icell; i++){
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
void usage(){
    mexErrMsgTxt("Usage: [var, [header]]=read('filename' [,howmany] [, start]). [] means optional\n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    file_t *fp;
    char *fn=NULL;
    int start=0, howmany=0;
    switch(nrhs){
    case 3:
	start=(long)mxGetScalar(prhs[2]);//starting block to read. matlab index.
    case 2:
	howmany=(long)mxGetScalar(prhs[1]);//do not break
    case 1:
	fn=mx2str(prhs[0]);
	break;
    default:
	usage();
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
    case 0:
    case 1:
	plhs[0]=readdata(fp, NULL, start, howmany); break;
    case 2:
	plhs[0]=readdata(fp, &plhs[1], start, howmany); break;
    default:
	usage();
    }
    if(start==0 && howmany==0){
	int res=zfeof(fp);
	if(res){
	    warning("There is unread data: res=%d\n", res);
	}
    }
    zfclose(fp);
}
