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
#include "matbin.h"
#include "dmat.h"
#include "cmat.h"
#include "defs.h"
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
	write_header(A->header, fp);
    }
    do_write(fp, 0, sizeof(T), M_T, p, nx, ny);
}
/**
   Function to write cell array of dense matrix data. into a file pointer
   Generally used by library developer
*/
void X(cellwritedata)(file_t *fp, const X(cell) *dc){
    if(dc) write_header(dc->header, fp);
    write_magic(MCC_ANY, fp);
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
   X(write)(A,"A") for double matrix.
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
    char *header=NULL;
    if(!magic){
	magic=read_magic(fp, &header);
    }
    if(magic!=M_T)
	error("%s is not a X(mat) file. magic=%x\n", zfname(fp), magic);
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    X(mat) *out;
    if(nx==0 || ny==0){
	out=NULL;
	free(header);
    }else{
	out=X(new)((long)nx,(long)ny);
	zfread(out->p,sizeof(T),nx*ny,fp);
	out->header=header;
    }
    return out;
}
/**
   Function to read dense matrix cell array into memory from file
   pointer. Generally used by library developer.  */
X(cell)* X(cellreaddata)(file_t *fp, uint32_t magic){
    char *header=NULL;
    if(!magic){
	magic=read_magic(fp, &header);
    }
    if(!iscell(magic)){
	error("This is is not a X(mat) cell file. want %d, get %d\n",(int)MCC_ANY,(int)magic);
    }
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    X(cell) *out;
    if(nx==0 || ny==0){
	out=NULL;
	free(header);
    }else{
	out=X(cellnew)((long)nx,(long)ny);
	out->header=header;
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
   file.

   Usage: A=dcellread("A.bin.gz"); for a double dcell.
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
   User callable function to read array of cell array of dense matrix from file. 

   Usage: A=dcellreadarr(&nx, &ny, filename);
*/
X(cell) **X(cellreadarr)(long *nxout, long *nyout, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn, "rb");
    uint32_t magic=read_magic(fp, NULL);
    if(!iscell(magic)){
	error("This is is not a cell file. want %d, get %d\n",(int)MCC_ANY,(int)magic);
    }
    uint64_t nx, ny;
    zfreadlarr(fp, 2, &nx, &ny); 
    X(cell) **out=calloc(nx*ny, sizeof(X(cell)*));
    for(long ic=0; ic<nx*ny; ic++){
	out[ic]=X(cellreaddata)(fp, 0);
    }
    *nxout=nx;
    *nyout=ny;
    zfeof(fp);
    zfclose(fp);
    return out;
}
/**
   User callable function to write array of cell array of dense matrix to file. 

   Usage: dcellwritearr(A, nx, ny, filename);
*/
void X(cellwritearr)(X(cell)**A, long nxin, long nyin, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    write_magic(MCC_ANY, fp);
    uint64_t nx, ny;
    nx=nxin;
    ny=nyin;
    zfwritelarr(fp, 2, &nx, &ny);
    for(long ic=0; ic<nxin*nyin; ic++){
	X(cellwritedata)(fp, A[ic]);
    }
    zfclose(fp);
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

/**
   Open a file for write with mmmap. We don't provide a access control here for
   generic usage of the function. Lock on a special dummy file for access
   control.
*/
static int mmap_open(char *fn, int rw){
    char *fn2=procfn(fn,rw?"w":"r",0);
    if(!fn2) return -1;
    if(fn2 && strlen(fn2)>=7&&!strncmp(fn2+strlen(fn2)-7,".bin.gz",7)){
	error("new_mmap does not support gzip\n");
    }
    int fd;
    if(rw){
	fd=open(fn2, O_RDWR|O_CREAT, 0600);
	if(fd!=-1 && ftruncate(fd, 0)){//truncate the file.
	    error("Unable to ftruncate file to 0 size\n");
	}
    }else{
	fd=open(fn2, O_RDONLY);
    }
    /*in read only mode, allow -1 to indicate failed. In write mode, fail.*/
    if(fd==-1 && rw){
	perror("open");
	error("Unable to create file %s\n", fn2);
    }
 
    free(fn2);
    return fd;
}
/**
   Initialize the header in the mmaped file.
*/
static inline void mmap_header_rw(char **p0, char **header0, uint32_t magic, long nx, long ny, const char *header){
    char *p=*p0;
    //Always have a header to align the data.
    if(header){
	uint64_t nlen=bytes_header(header)-24;
	((uint32_t*)p)[0]=(uint32_t)M_HEADER; p+=4;
	((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
	*header0=p;
	memcpy(p, header, strlen(header)+1); p+=nlen;
	((uint64_t*)p)[0]=(uint64_t)nlen; p+=8;
	((uint32_t*)p)[0]=(uint32_t)M_HEADER;p+=4;
    }else{
	*header0=NULL;
    }
    ((uint32_t*)p)[0]=(uint32_t)M_SKIP; p+=4;
    ((uint32_t*)p)[0]=(uint32_t)magic; p+=4;
    ((uint64_t*)p)[0]=(uint64_t)nx; p+=8;
    ((uint64_t*)p)[0]=(uint64_t)ny; p+=8;
    *p0=p;
}
/**
   Create a new X(mat) matrix object, mmapped from file. be aware that the data
   is not 8-byte aligned. The file is truncated if already exists in rw mode.*/
X(mat)* X(new_mmap)(long nx, long ny, char *header, const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn, 1);
    size_t metasize=3*8+bytes_header(header);
    size_t msize=nx*ny*sizeof(T)+metasize;
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize, (PROT_WRITE|PROT_READ), MAP_SHARED, fd, 0);
    if(!map){
	error("mmap failed\n");
    }
    char *map0=map;
    //memset(map, 0, msize);
    char *header0;
    mmap_header_rw(&map, &header0, M_T, nx, ny, header);
    X(mat) *out=X(new_data)(nx, ny, (T*)map);
    out->header=header0;
    out->mmap=mmap_new(fd, map0, msize);
    return out;
}
/**
   Create a new X(cell) matrix cell object, mmapped from file. be aware that the
   data is not 8-byte aligned. The file is truncated if already exists. We only
   add headers to individual dmat/cmat in the file, not in the cell.
*/
X(cell)* X(cellnew_mmap)(long nx, long ny, long *nnx, long *nny, char *header1, char **header2,
			 const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn, 1);
    long metasize=3*8;
    long msize=metasize+bytes_header(header1);
    for(long ix=0; ix<nx*ny; ix++){
	msize+=metasize+nnx[ix]*nny[ix]*sizeof(T)+bytes_header(header2?header2[ix]:NULL);
    }
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize, (PROT_WRITE|PROT_READ), MAP_SHARED, fd, 0);
    if(!map){
	error("mmap failed\n");
    }
    char *map0=map;
    //memset(map, 0, msize);
    char *header0;
    mmap_header_rw(&map, &header0, MCC_ANY, nx, ny, header1);
    X(cell) *out=X(cellnew)(nx,ny);
    out->mmap=mmap_new(fd, map0, msize);
    out->header=header0;
    for(long ix=0; ix<nx*ny; ix++){
	mmap_header_rw(&map, &header0, M_T, nnx[ix], nny[ix], header2?header2[ix]:NULL);
	out->p[ix]=X(new_data)(nnx[ix], nny[ix], (T*)map);
	map+=nnx[ix]*nny[ix]*sizeof(T);
	if(out->p[ix]) {
	    out->p[ix]->mmap=mmap_ref(out->mmap);//reference
	    out->p[ix]->header=header0;
	}
    }
    
    return out;
}
/**
   Create a new X(cell) matrix cell object, with identical blocks, mmapped from
   file. be aware that the data is not 8-byte aligned. The file is truncated if
   already exists.  */
X(cell)* X(cellnewsame_mmap)(long nx, long ny, long mx, long my, char *header,
			     const char *format, ...){
    if(!nx || !ny) return NULL;
    format2fn;
    int fd=mmap_open(fn, 1);
    long metasize=3*8;
    long msize=metasize;
    msize=metasize+bytes_header(header)+nx*ny*(metasize+mx*my*sizeof(T));
    if(ftruncate(fd, msize)){
	error("Error truncating file\n");
    }
    char *map=mmap(NULL, msize,(PROT_WRITE|PROT_READ), MAP_SHARED, fd, 0);
    if(!map){
	error("mmap failed\n");
    }
    char *map0=map;
    //memset(map, 0, msize);
    X(cell) *out=X(cellnew)(nx,ny);
    char *header0;
    mmap_header_rw(&map, &header0, MCC_ANY, nx, ny, header);
    out->mmap=mmap_new(fd, map0, msize);
    out->header=header0;
    for(long ix=0; ix<nx*ny; ix++){
	mmap_header_rw(&map, &header0, M_T, mx, my, NULL);
	out->p[ix]=X(new_data)(mx, my, (T*)map);
	map+=mx*my*sizeof(T);
	if(out->p[ix]){
	    out->p[ix]->mmap=mmap_ref(out->mmap);;
	    out->p[ix]->header=header0;
	}
    }
    return out;
}
/**
   Initialize the dimension from the header in the mmaped file.
*/
static inline void mmap_header_ro(char **p0, uint32_t *magic, long *nx, long *ny, char **header0){
    char *p=*p0;
    char *header=NULL;
    while(((uint32_t*)p)[0]==M_HEADER){
	p+=4;
	long nlen=((uint64_t*)p)[0];p+=8;
	header=p;
	p+=nlen;
	if(nlen == ((uint64_t*)p)[0]){
	    p+=8;
	    if(((uint32_t*)p)[0]==M_HEADER){
		p+=4;
	    }else{
		header=NULL;
	    }
	}else{
	    header=NULL;
	}
	if(!header){
	    error("Parse head failed\n");
	}
    }
    if(!header){
	p=*p0;
    }
    if(((uint32_t*)p)[0]==M_SKIP) p+=4;
    *magic=((uint32_t*)p)[0]; p+=4;
    *nx=((uint64_t*)p)[0]; p+=8;
    *ny=((uint64_t*)p)[0]; p+=8;
    *p0=p;
    if(header0) *header0=header;
}
static X(mat*) X(readdata_mmap)(char **map){
    long nx, ny;
    uint32_t magic;
    char *header;
    mmap_header_ro(map, &magic, &nx, &ny, &header);
    if(magic!=M_T){
	error("File has magic %d, we want %d\n", (int)magic, M_T);
    }
    X(mat)* out=X(new_data)(nx, ny, (T*)(*map));
    *map+=nx*ny*sizeof(T);
    if(out){
	out->header=header;
    }
    return out;
}
/**
   Map the file to memory in read only, shared mode.
*/
X(mat*) X(read_mmap)(const char *format, ...){
    format2fn;
    int fd=mmap_open(fn, 0);
    if(fd<0) return NULL;
    long msize=flen(fn);
    char *map=mmap(NULL, msize, PROT_READ, MAP_SHARED, fd, 0);
    if(!map){
	error("mmap failed\n");
    }
    char *map0=map;
    X(mat) *out=X(readdata_mmap)(&map);
    if(out){
	out->mmap=mmap_new(fd, map0,msize);
    }
    return out;
}

/**
   Map the file to memory in read only, shared mode.
*/
X(cell*) X(cellread_mmap)(const char *format, ...){
    format2fn;
    int fd=mmap_open(fn, 0);
    if(fd<0) return NULL;
    long msize=flen(fn);
    char *map=mmap(NULL, msize, PROT_READ, MAP_SHARED, fd, 0);
    if(!map){
	error("mmap failed\n");
    }
    char *map0=map;
    long nx, ny;
    uint32_t magic;
    char *header;
    mmap_header_ro(&map, &magic, &nx, &ny, &header);
    if(!iscell(magic)){
	error("We want a cell array, File has %d\n", (int)magic);
    }
    X(cell) *out=X(cellnew)(nx, ny);
    out->mmap=mmap_new(fd, map0, msize);
    out->header=header;
    for(long ix=0; ix<nx*ny; ix++){
	out->p[ix]=X(readdata_mmap)(&map);
	out->p[ix]->mmap=mmap_ref(out->mmap);
    }
    return out;
}
