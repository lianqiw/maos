#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_FN_LEN 800
/*compile with 
mex write.c -largeArrayDims
usage write(filename, data);
*/
#include "io.h"

static void writedata(file_t *fp, int type, const mxArray *arr){
    int magic;
    long m,n;
    if(!arr){
	m=0; n=0;
    }else{
	m=mxGetM(arr);
	n=mxGetN(arr);
    }
    if(arr && mxIsCell(arr)){
	int issparse=0;
	mxArray *in;
	int type2;
	long ix;
	if(mxGetNumberOfElements(arr)==0){
	    in=NULL;
	    issparse=0;
	    magic=MC_DBL;
	    type2=M_DBL;
	}else{
	    in=mxGetCell(arr,0);
	    if(mxIsSparse(in)){
		issparse=1;
		if(mxIsComplex(in)){
		    magic=MC_CSP;
		    type2=M_CSP;
		}else{
		    magic=MC_SP;
		    type2=M_SP;
		}
	    }else{
		issparse=0;
		if(mxIsComplex(in)){
		    magic=MC_CMP;
		    type2=M_CMP;
		}else{
		    magic=MC_DBL;
		    type2=M_DBL;
		}
	    }
	}
	writefile(&magic, sizeof(int), 1, fp);
	writefile(&m, sizeof(long), 1, fp);
	writefile(&n, sizeof(long), 1, fp);
	for(ix=0; ix<mxGetNumberOfElements(arr); ix++){
	    in=mxGetCell(arr, ix);
	    if(in && mxIsSparse(in) !=issparse)
		error("can only save cell array of all sparse or all dense");
	    writedata(fp, type2, in);
	}
    }else if(type == M_SP || ((arr) && mxIsSparse(arr))){
	magic=M_SP;
	writefile(&magic, sizeof(int), 1, fp);
	writefile(&m, sizeof(long), 1, fp);
	writefile(&n, sizeof(long), 1, fp);
	if(m!=0 && n!=0){
	    mwIndex *Jc=mxGetJc(arr);
	    long nzmax=Jc[n];
	    writefile(&nzmax, sizeof(double), 1, fp);
	    writefile(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
	    writefile(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
	    writefile(mxGetPr(arr), sizeof(double), nzmax, fp);
	}
    }else if(type == M_CSP || ((arr) && mxIsSparse(arr))){
	magic=M_CSP;
	writefile(&magic, sizeof(int), 1, fp);
	writefile(&m, sizeof(long), 1, fp);
	writefile(&n, sizeof(long), 1, fp);
	if(m!=0 && n!=0){
	    mwIndex *Jc=mxGetJc(arr);
	    long nzmax=Jc[n];
	    writefile(&nzmax, sizeof(double), 1, fp);
	    writefile(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
	    writefile(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
	    writefile_complex(mxGetPr(arr),mxGetPi(arr),nzmax,fp);
	}
    }else if(type == M_DBL || ((arr)&& mxIsDouble(arr))){
	magic=M_DBL;
	writefile(&magic, sizeof(int), 1, fp);
	writefile(&m, sizeof(long), 1, fp);
	writefile(&n, sizeof(long), 1, fp);
	if(m!=0 && n!=0){
	    writefile(mxGetPr(arr), sizeof(double), m*n, fp);
	}  
    }else if(type == M_CMP || ((arr)&& mxIsDouble(arr))){
	magic=M_CMP;
	writefile(&magic, sizeof(int), 1, fp);
	writefile(&m, sizeof(long), 1, fp);
	writefile(&n, sizeof(long), 1, fp);
	if(m!=0 && n!=0){
	    writefile_complex(mxGetPr(arr),mxGetPi(arr), m*n, fp);
	}
    }else{
	mexErrMsgTxt("Unrecognized data type");
    }
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    file_t *fp;
    char fn[MAX_FN_LEN];
    mxGetString(prhs[1],fn,MAX_FN_LEN+1);
    /*if(fn[0]=='~'){
	const char *HOME=getenv("HOME");
	if(!HOME){
	    error("Unable to get user home\n");
	}
	int fnlen=strlen(fn);
	int homelen=strlen(HOME);
	if(fnlen+homelen>=MAX_FN_LEN){
	    error("MAX_FN_LEN is too small\n");
	}
	memmove(fn+homelen, fn+1, fnlen);
	memcpy(fn, HOME, homelen);
	}*/
    fp=openfile(fn,"w");
    (void)nlhs;
    (void)nrhs;
    (void)plhs;
    writedata(fp, 0, prhs[0]);
    closefile(fp);
}
