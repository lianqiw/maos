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
#ifdef DLONG
#define M_SPT M_SPT64
#else
#define M_SPT M_SPT32
#endif
/**
   Function to write sparse matrix data into file pointed using a file
   pointer. Generally used by library developer.  We do not convert data during
   saving, but rather do the conversion during reading.

 */
void Y(spwritedata)(file_t *fp, const X(sp) *sp){
    uint32_t magic=M_SPT;
    Y(spsort)((X(sp)*)sp);//sort the matrix to have the right order
    zfwrite(&magic, sizeof(uint32_t),1,fp);
    if(sp && sp->nzmax){
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
    uint64_t m,n,nzmax;
    zfreadlarr(fp, 2, &m, &n);
    X(sp) *out=NULL;
    if(m!=0 && n!=0){
	uint32_t magic2=0;
	switch(magic){
	case M_SPT64:
	    magic2=M_INT64;
	    break;
	case M_SPT32:
	    magic2=M_INT32;
	    break;
	default:
	    error("This is not a valid sparse matrix file. magic=%x\n", magic);
	}
	zfread(&nzmax,sizeof(uint64_t),1,fp);
	if(nzmax!=0){
	    out=Y(spnew)(m,n,nzmax);
	    readspintdata(fp, magic2, out->p, n+1);
	    readspintdata(fp, magic2, out->i, nzmax);
	    zfread(out->x, sizeof(T), nzmax, fp);
	}
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
