/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>

#include <fcntl.h>
#include <unistd.h>
#include "type.h"
#include "spbin.h"
#include "dsp.h"
#include "csp.h"
#ifdef USE_COMPLEX
#define T dcomplex
#define X(A) c##A
#define Y(A) c##A
#define M_SPT64 M_CSP64
#define M_SPT32 M_CSP32
#else
#define T double
#define X(A) d##A
#define Y(A) A
#define M_SPT64 M_SP64
#define M_SPT32 M_SP32
#endif
/**
   Function to write sparse matrix data into file pointed using a file
   pointer. Generally used by library developer.  We do not convert data during
   saving, but rather do the conversion during reading.

 */
void Y(spwritedata)(file_t *fp, const X(sp) *sp){
    uint32_t magic;
    if(sizeof(spint)==4)
	magic=M_SPT32;
    else if(sizeof(spint)==8)
	magic=M_SPT64;
    else
	error("Invalid");
    Y(spsort)((X(sp)*)sp);//sort the matrix to have the right order
    zfwrite(&magic, sizeof(uint32_t),1,fp);
    if(sp){
	uint64_t m,n,nzmax;
	m=sp->m;
	n=sp->n;
	nzmax=sp->p[n];//don't use sp->nzmax, which maybe larger than actual
	zfwritelarr(fp, 3, &m, &n, &nzmax);
	zfwrite(sp->p, sizeof(spint), n+1, fp);
	zfwrite(sp->i, sizeof(spint), nzmax, fp);
	zfwrite(sp->x ,sizeof(T),nzmax,fp);  
    }else{
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }
}
/**
   Function to read sparse matrix data from file pointer into memory. Used by
   library developer.
  */
X(sp) *Y(spreaddata)(file_t *fp, uint32_t magic){
    if(!magic){
	magic=read_magic(fp, NULL);
    }
    uint32_t size=0;
    if(magic==M_SPT64){
	size=8;
    }else if(magic==M_SPT32){
	size=4;
    }else{
	error("This is not a valid sparse matrix file\n");
    }
    uint64_t m,n,nzmax;
    zfreadlarr(fp, 2, &m, &n);
    X(sp) *out;
    if(m==0 || n==0){
	out=NULL;
    }else{
	zfread(&nzmax,sizeof(uint64_t),1,fp);
	out=Y(spnew)(m,n,nzmax);
	if(sizeof(spint)==size){//data match.
	    zfread(out->p, sizeof(spint), n+1, fp);
	    zfread(out->i, sizeof(spint), nzmax, fp);
	}else if(size==8){//convert uint64_t to uint32_t
	    assert(sizeof(spint)==4);
	    uint64_t *p=malloc(sizeof(uint64_t)*(n+1));
	    uint64_t *i=malloc(sizeof(uint64_t)*nzmax);
	    zfread(p, sizeof(uint64_t), n+1, fp);
	    zfread(i, sizeof(uint64_t), nzmax, fp);
	    for(unsigned long j=0; j<n+1; j++){
		out->p[j]=(spint)p[j];
	    }
	    for(unsigned long j=0; j<nzmax; j++){
		out->i[j]=(spint)i[j];
	    }
	    free(p);
	    free(i);
	}else if(size==4){//convert uint32_t to uint64_t
	    assert(sizeof(spint)==8);
	    uint32_t *p=malloc(sizeof(uint32_t)*(n+1));
	    uint32_t *i=malloc(sizeof(uint32_t)*nzmax);
	    zfread(p, sizeof(uint32_t), n+1, fp);
	    zfread(i, sizeof(uint32_t), nzmax, fp);
	    for(unsigned long j=0; j<n+1; j++){
		out->p[j]=(spint)p[j];
	    }
	    for(unsigned long j=0; j<nzmax; j++){
		out->i[j]=(spint)i[j];
	    }
	    free(p);
	    free(i);
	}
	zfread(out->x, sizeof(T),nzmax, fp);
    }
    return out;
}

/**
   User callable function to write sparse matrix into file. 

   Usage: spwrite(A,"A.bin.gz");
*/
void Y(spwrite)(const X(sp) *sp, const char *format,...){
    format2fn;
    //write the sparse matrix to file to later load from matlab
    file_t *fp=zfopen(fn,"wb");
    Y(spwritedata)(fp, sp);
    //don't worry about the warning of 0x401ee45 in valgrind. That is the IO 
    zfclose(fp);
}
/**
   User callable function to write cell array of sparse matrix into file. 

   Usage: spcellwrite(A,"A.bin.gz"); */
void Y(spcellwrite)(const Y(spcell) *spc, const char *format,...){
    format2fn;
    uint32_t magic=MCC_ANY;
    file_t *fp=zfopen(fn,"wb");
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    if(spc){
	uint64_t nx=spc->nx;
	uint64_t ny=spc->ny;
	zfwritelarr(fp, 2, &nx, &ny);
	for(unsigned long iy=0; iy<spc->ny; iy++){
	    for(unsigned long ix=0; ix<spc->nx; ix++){
		Y(spwritedata)(fp, spc->p[ix+iy*spc->nx]);
	    }
	}
    }else{
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }
    zfclose(fp);
}
/**
   User callable function to read sparse metrix from file. 

   Usage: A=spread("A.bin.gz");*/
X(sp)* Y(spread)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    X(sp) *out=Y(spreaddata)(fp, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}
/**
   User callable function to read cell array of sparse matrix
   from file. Usage: A=spcellread("A.bin.gz");
 */
Y(spcell) *Y(spcellread)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    uint32_t magic=read_magic(fp, NULL);
    if(!iscell(magic)){
	error("%s is not a sparse cell file. want %d, got %u\n", fn, MCC_ANY, magic);
    }
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    Y(spcell) *out;
    if(nx==0 || ny==0)
	out=NULL;
    else{
	out=Y(spcellnew)(nx,ny);
	for(unsigned long ix=0; ix<nx*ny; ix++){
	    out->p[ix]=Y(spreaddata)(fp, 0);
	}
    }
    zfeof(fp);
    zfclose(fp);
    return out;
}
