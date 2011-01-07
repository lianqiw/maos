/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <fcntl.h>
#include <unistd.h>
#include "type.h"
#include "matbin.h"
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"
#include "csp.h"
#ifdef USE_COMPLEX
#define T dcomplex
#define X(A) c##A
#define Y(A) c##A
#define M_T M_CMP
#define M_SPT64 M_CSP64
#define M_SPT32 M_CSP32
#else
#define T double
#define X(A) d##A
#define Y(A) A
#define M_T M_DBL
#define M_SPT64 M_SP64
#define M_SPT32 M_SP32
#endif
/**
   Contains routines to write/read dense/sparse matrix into/from file.
*/

/**
   Function to write dense matrix data into a file pointer. Generally used by
   library developer */
void X(writedata)(file_t *fp, const X(mat) *A){
    void *p=NULL;
    uint64_t nx=0, ny=0;
    if(A){
	p=A->p;
	nx=(uint64_t)A->nx;
	ny=(uint64_t)A->ny;
    }
    do_write(fp, 0, sizeof(T), M_T, p, nx, ny);
}
/**
   Function to write cell array of dense matrix data. into a file pointer
   Generally used by library developer
 */
void X(cellwritedata)(file_t *fp, const X(cell) *dc){
    uint32_t magic=MCC_ANY;
    zfwrite(&magic, sizeof(uint32_t), 1, fp);
    if(!dc){
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }else{
	uint64_t nx=dc->nx;
	uint64_t ny=dc->ny;
	zfwritelarr(fp, 2, &nx, &ny);
	for(unsigned long iy=0; iy<ny; iy++){
	    for(unsigned long ix=0; ix<nx; ix++){
		X(writedata)(fp, dc->p[ix+iy*nx]);
	    }
	}
    }
}
/**
   User callable function to write dense matrix into a file. Usage:
   dwrite(A,"A") for double matrix.
  */
void X(write)(const X(mat) *A, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    X(writedata)(fp, A);
    zfclose(fp);
}
/**
   User callable function to write cell array of dense matrix into a
   file. Usage: dcellwrite(A,"A.bin.gz") for double matrix cell. */
void X(cellwrite)(const X(cell) *dc, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    X(cellwritedata)(fp,dc);
    zfclose(fp);
}
/**
   Function to read dense matrix into memory from file pointer. Generally used
   by library developer.  */
X(mat) *X(readdata)(file_t *fp, uint32_t magic){
    if(!magic){
	magic=read_magic(fp, NULL);
    }
    if(magic!=M_T)
	error("This is not a X(mat) file\n");
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    X(mat) *out;
    if(nx==0 || ny==0)
	out=NULL;
    else{
	out=X(new)((long)nx,(long)ny);
	zfread(out->p,sizeof(T),nx*ny,fp);
    }
    return out;
}
/**
   Function to read dense matrix cell array into memory from file
   pointer. Generally used by library developer.  */
X(cell)* X(cellreaddata)(file_t *fp, uint32_t magic){
    if(!magic){
	magic=read_magic(fp, NULL);
    }
    if(!iscell(magic)){
	error("This is is not a X(mat) cell file. want %d, get %d\n",(int)MCC_ANY,(int)magic);
    }
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    X(cell) *out;
    if(nx==0 || ny==0)
	out=NULL;
    else{
	out=X(cellnew)((long)nx,(long)ny);
	for(unsigned long ix=0; ix<nx*ny; ix++){
	    out->p[ix]=X(readdata)(fp, 0);
	}
    }
    return out;
}
/**
   User callable function to read dense matrix into memory from file. Usage:
   A=dread("A.bin.gz"); for a dmat*/
X(mat)* X(read)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    X(mat) *out=X(readdata)(fp, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}
/**
   User callable function to read cell array of dense matrix into memory from
   file. Usage: A=dcellread("A.bin.gz"); for a double dcell.
 */
X(cell)* X(cellread)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    X(cell) *out=X(cellreaddata)(fp, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}
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
/**
   Open a file for write with mmmap.
 */
static int mmap_open(char *fn){
    char *fn2=procfn(fn,"w",0);
    if(strlen(fn2)>=7&&!strncmp(fn2+strlen(fn2)-7,".bin.gz",7)){
	error("new_mmap does not support gzip\n");
    }
    int fd=open(fn2,O_RDWR|O_CREAT,0600);
    if(fd==-1){
	perror("open");
	error("Error ");
    }
    free(fn2);
    return fd;
}
/**
   Initialize the header in the mmaped file.
 */
static inline void mmap_header(char *map, uint32_t magic, long nx, long ny){
    uint32_t *mm=(uint32_t*)map;
    mm[0]=magic;
    long *nn=(long*)(map+sizeof(uint32_t));
    nn[0]=nx;
    nn[1]=ny;
}
/**
   Create a new X(mat) matrix object, mmapped from file. be aware that the data
   is not 8-byte aligned. The file is truncated if already exists.*/
X(mat)* X(new_mmap)(long nx, long ny, const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn);
    size_t headersize=2*sizeof(long)+sizeof(uint32_t);
    size_t msize=nx*ny*sizeof(T)+headersize;
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    memset(map, 0, msize);
    close(fd);
    mmap_header(map, M_T, nx, ny);
    X(mat) *out=X(new_data)((T*)(map+headersize), nx, ny);
    out->type=MT_MMAP;
    return out;
}
/**
   Create a new X(cell) matrix cell object, mmapped from file. be aware that the data
   is not 8-byte aligned. The file is truncated if already exists.
*/
X(cell)* X(cellnew_mmap)(long nx, long ny, long *nnx, long *nny, const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn);
    long headersize=sizeof(long)*2+sizeof(uint32_t);
    long msize=headersize;
    for(long ix=0; ix<nx*ny; ix++){
	msize+=headersize+nnx[ix]*nny[ix]*sizeof(T);
    }
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    memset(map, 0, msize);
    close(fd);
    X(cell) *out=X(cellnew)(nx,ny);
    out->mmap=map;
    mmap_header(map, MCC_ANY, nx, ny);
    map+=headersize;//header of cell
    for(long ix=0; ix<nx*ny; ix++){
	mmap_header(map, M_T, nnx[ix], nny[ix]);
	out->p[ix]=X(new_data)((T*)(map+headersize), nnx[ix], nny[ix]);
	if(out->p[ix]) {
	    out->p[ix]->type=MT_MMAP;
	}
	map+=nnx[ix]*nny[ix]*sizeof(T)+headersize;
    }

    return out;
}
/**
   Create a new X(cell) matrix cell object, with identical blocks, mmapped from
   file. be aware that the data is not 8-byte aligned. The file is truncated if
   already exists.  */
X(cell)* X(cellnewsame_mmap)(long nx, long ny, long mx, long my, const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn);
    long headersize=sizeof(long)*2+sizeof(uint32_t);
    long msize=headersize;
    msize=headersize+nx*ny*(headersize+mx*my*sizeof(T));
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    memset(map, 0, msize);
    close(fd);
    X(cell) *out=X(cellnew)(nx,ny);
    out->mmap=map;
    mmap_header(map, MCC_ANY, nx, ny);
    map+=headersize;//header of cell
    for(long ix=0; ix<nx*ny; ix++){
	mmap_header(map, M_T, mx, my);
	out->p[ix]=X(new_data)((T*)(map+headersize), mx, my);
	if(out->p[ix]){
	    out->p[ix]->type=MT_MMAP;
	}
	map+=mx*my*sizeof(T)+headersize;
    }
    return out;
}


/**
   Scale a dcell array and save to file.
*/
void X(cellswrite)(X(cell) *A, double scale, const char *format, ...){
    format2fn;
    X(cell) *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    X(celladd)(&tmp, 0, A, scale);
    X(cellwrite)(tmp,"%s",fn);
    X(cellfree)(tmp);
}

/**
   Scale a dcell array and save to file.
*/
void X(swrite)(X(mat) *A, double scale, const char *format, ...){
    format2fn;
    X(mat) *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    X(add)(&tmp, 0, A, scale);
    X(write)(tmp,"%s",fn);
    X(free)(tmp);
}
