/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



/*compile with 
  mex write.c -largeArrayDims
  usage write(filename, data);
*/
#include "io.h"
static char *mx2str(const mxArray *header){
    char *str=NULL;
    if(header && mxIsChar(header)){
	int nheader=mxGetM(header)*mxGetN(header)+1;
	str=(char*)malloc(nheader);
	mxGetString(header, str, nheader);
	//convert matlab \n (2 chars) to C \n (1 char with a space)
	const char *str2=str+strlen(str)-1;
	for(char *h=str; h<str2; h++){
	    if(h[0]=='\\' && h[1]=='n'){
		h[0]=' ';
		h[1]='\n';
	    }
	}
    }
    return str;
}
static void writedata(file_t *fp, const mxArray *arr, const mxArray *header){
    char *str=mx2str(header);
    header_t header2={0,0,0,0,str};
    //uint64_t m,n;
    if(arr){
	header2.ndim=mxGetNumberOfDimensions(arr);
	header2.dims=(mwSize*)mxGetDimensions(arr);
	if(header2.ndim==0){
	    header2.ntot=0;
	}else{
	    header2.ntot=1;
	    for(unsigned long id=0; id<header2.ndim; id++){
		header2.ntot*=header2.dims[id];
	    }
	}
	if(!header2.ntot){
	    arr=NULL;
	}
    }
    if(!arr){
	header2.ndim=0;
	header2.dims=0;
	header2.magic=M_DBL;
	write_header(&header2, fp);
	return;
    }
    if(mxIsCell(arr)){
	header2.magic=MCC_ANY;
	write_header(&header2, fp);
	for(size_t ix=0; ix<mxGetNumberOfElements(arr); ix++){
	    const mxArray *in=mxGetCell(arr, ix);
	    if(header && mxIsCell(header)){
		writedata(fp, in, mxGetCell(header, ix));
	    }else{
		writedata(fp, in, header);
	    }
	}
    }else{/*not cell.*/
	int ntot=0;
	int byte=0;
	int magic;
	if(mxIsComplex(arr)){
	    switch(mxGetClassID(arr)){
	    case mxDOUBLE_CLASS:
		byte=16;
		magic=M_CMP;
		break;
	    case mxSINGLE_CLASS:
		byte=8;
		magic=M_ZMP;
		break;
	    default:
		error("Please implement this class type\n");
	    }
	}else{
	    switch(mxGetClassID(arr)){
	    case mxDOUBLE_CLASS:
		byte=8;
		magic=M_DBL;
		break;
	    case mxSINGLE_CLASS:
		byte=4;
		magic=M_FLT;
		break;
	    case mxINT64_CLASS:
		byte=8;
		magic=M_INT64;
		break;
	    case mxINT32_CLASS:
		byte=4;
		magic=M_INT32;
		break;
	    case mxINT16_CLASS:
		byte=2;
		magic=M_INT16;
		break;
	    case mxINT8_CLASS:
		byte=1;
		magic=M_INT8;
		break;
	    default:
		error("Please implement this class type\n");
	    }
	}
	
	if(mxIsSparse(arr)){
	    if(header2.ndim!=2){
		zfclose(fp);
		error("Invalid format. Must have 2 dimensions.\n");
	    }
	    if(!fp->isfits){
		if(mxIsComplex(arr)){
		    if(sizeof(mwIndex)==4){
			header2.magic=M_CSP32;
		    }else{
			header2.magic=M_CSP64;
		    }
		}else{
		    if(sizeof(mwIndex)==4){
			header2.magic=M_SP32;
		    }else{
			header2.magic=M_SP64;
		    }
		}
		write_header(&header2, fp);
	    }
	    if(header2.ntot){
		long ny=header2.dims[1];
		const mwIndex *Jc=mxGetJc(arr);
		const mwIndex *Ir=mxGetIr(arr);
		const uint64_t nz=Jc[ny];
		if(!fp->isfits){
		    zfwrite(&nz, sizeof(uint64_t), 1, fp);
		}
		//Do not use getNz. It shows available slots, but not actual data.
		if(fp->isfits){
		    //write Jc header
		    header2.magic=(sizeof(mwIndex)==4?M_INT32:M_INT64);
		    header2.str="Jc";
		    header2.dims[0]=ny+1;
		    header2.dims[1]=1;
		    write_header(&header2, fp);
		}
		zfwrite(Jc, sizeof(mwIndex), ny+1, fp);
		if(fp->isfits){
		    //write Ir header
		    header2.magic=(sizeof(mwIndex)==4?M_INT32:M_INT64);
		    header2.str="Ir";
		    header2.dims[0]=nz;
		    header2.dims[1]=1;
		    write_header(&header2, fp);
		}
		zfwrite(Ir, sizeof(mwIndex), nz, fp);
		if(fp->isfits){
		    //Write P header
		    header2.magic=(mxIsComplex(arr)?M_CMP:M_DBL);
		    header2.str="P";
		    header2.dims[0]=nz;
		    header2.dims[1]=1;
		    write_header(&header2, fp);
		}
		ntot=nz;
	    }
	}else{
	    ntot=header2.ntot;
	    header2.magic=magic;
	    write_header(&header2, fp);
	}
#if !MX_HAS_INTERLEAVED_COMPLEX
	if(mxIsComplex(arr)){//convert from continuous to interleaved data.
	    void *Pr0=mxGetPr(arr);
	    void *Pi0=mxGetPi(arr);
	    void *tmp0=malloc(ntot*byte);
	    if(mxIsDouble(arr)){
		dcomplex *tmp=tmp0;
		double *Pr=Pr0;
		double *Pi=Pi0;
		for(size_t i=0; i<ntot; i++){
		    tmp[i].x=Pr[i];
		    tmp[i].y=Pi[i];
		}
	    }else{
		fcomplex *tmp=tmp0;
		float *Pr=Pr0;
		float *Pi=Pi0;
		for(size_t i=0; i<ntot; i++){
		    tmp[i].x=Pr[i];
		    tmp[i].y=Pi[i];
		}
	    }
	    
	    zfwrite(tmp0, byte, ntot, fp);
	    free(tmp0);
	}else
#endif
	    zfwrite(mxGetData(arr), byte, ntot, fp);
    }
    free(str);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    (void)nlhs;
    (void)plhs;
    file_t *fp;
    int ifn=1;
    const mxArray *header=NULL;
    if(nrhs==2){/*data and file*/
	ifn=1;
    }else if(nrhs==3){/*data, header, and file*/
	ifn=2;
	header=prhs[1];
    }else{
	mexErrMsgTxt("Usage: write(var,'file') or write(var, header, 'file')\n");
    }
    int nlen=mxGetM(prhs[ifn])*mxGetN(prhs[ifn])+1;
    char *fn=(char*)malloc(nlen);
    mxGetString(prhs[ifn],fn,nlen);
    fp=zfopen(fn,"wb");
    if(!fp){
	mexErrMsgTxt("Error writing file.\n");
	return;
    }
    free(fn);
    writedata(fp, prhs[0], header);
    zfclose(fp);
}
