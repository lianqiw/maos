#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*compile with 
mex write.c -largeArrayDims
usage write(filename, data);
*/
#include "io.h"

static void writedata(file_t *fp, int type, const mxArray *arr, const mxArray *header){
    char *str=NULL;
    if(header && mxIsChar(header)){
	int nheader=mxGetM(header)*mxGetN(header)+1;
	str=malloc(nheader);
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
    uint32_t magic;
    uint64_t m,n;
    if(!arr){
	m=0; n=0;
    }else{
	m=mxGetM(arr);
	n=mxGetN(arr);
    }
    if(!m || !n){
	arr=NULL;
    }
    if(arr && mxIsCell(arr)){
	magic=MCC_ANY;
	int issparse=0;
	mxArray *in;
	int type2;
	long ix;
	if(mxGetNumberOfElements(arr)==0){
	    in=NULL;
	    issparse=0;
	    type2=M_DBL;
	}else{
	    for(ix=0; ix<mxGetNumberOfElements(arr); ix++){
		in=mxGetCell(arr,ix);
		if(in){
		    if(mxIsSparse(in)){
			issparse=1;
			if(mxIsComplex(in)){
			    type2=MAT_CSP;
			}else{
			    type2=MAT_SP;
			}
		    }else{
			issparse=0;
			if(mxIsComplex(in)){
			    if(mxIsSingle(in)){
				type2=M_ZMP;
			    }else{
				type2=M_CMP;
			    }
			}else{
			    if(mxIsSingle(in)){
				type2=M_FLT;
			    }else{
				type2=M_DBL;
			    }
			}
		    }
		    break;
		}
	    }
	}
	if(!in){//all cell empty
	    issparse=0;
	    type2=M_DBL;
	}
	header_t header2={magic, m, n, str};
	write_header(&header2, fp);
	for(ix=0; ix<mxGetNumberOfElements(arr); ix++){
	    in=mxGetCell(arr, ix);
	    if(in && !mxIsEmpty(in) && mxIsSparse(in) !=issparse)
		error("can only save cell array of all sparse or all dense");
	    if(header && mxIsCell(header)){
		writedata(fp, type2, in, mxGetCell(header, ix));
	    }else{
		writedata(fp, type2, in, NULL);
	    }
	}
    }else{/*not cell.*/
	if(type == MAT_SP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		magic=M_SP32;
	    }else{
		magic=M_SP64;
	    }
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		mwIndex *Jc=mxGetJc(arr);
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite(mxGetPr(arr), sizeof(double), nzmax, fp);
	    }
	}else if(type == MAT_CSP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		magic=M_CSP32;
	    }else{
		magic=M_CSP64;
	    }
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		mwIndex *Jc=mxGetJc(arr);
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite_dcomplex(mxGetPr(arr),mxGetPi(arr),nzmax,fp);
	    }
	}else if(type == M_DBL || ((arr) && mxIsDouble(arr) && !mxIsComplex(arr))){
	    magic=M_DBL;
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		zfwrite(mxGetPr(arr), sizeof(double), m*n, fp);
	    }  
	}else if(type == M_FLT || ((arr) && mxIsSingle(arr) && !mxIsComplex(arr))){
	    magic=M_FLT;
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		zfwrite(mxGetPr(arr), sizeof(float), m*n, fp);
	    }
	}else if(type == M_CMP || ((arr)&& mxIsDouble(arr) && mxIsComplex(arr))){
	    magic=M_CMP;
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		zfwrite_dcomplex(mxGetPr(arr),mxGetPi(arr), m*n, fp);
	    }
	}else if(type == M_ZMP || ((arr)&& mxIsSingle(arr) && mxIsComplex(arr))){
	    magic=M_ZMP;
	    header_t header2={magic, m, n, str};
	    write_header(&header2, fp);
	    if(m!=0 && n!=0){
		zfwrite_fcomplex((float*)mxGetPr(arr),(float*)mxGetPi(arr), m*n, fp);
	    }
	}else{
	    error("Unrecognized data type");
	}
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
    char *fn=malloc(nlen);
    mxGetString(prhs[ifn],fn,nlen);
    fp=zfopen(fn,"wb");
    if(!fp){
	mexErrMsgTxt("Error writing file.\n");
	return;
    }
    free(fn);
    writedata(fp, 0, prhs[0], header);
    zfclose(fp);
}
