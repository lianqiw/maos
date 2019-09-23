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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "type.h"
#include "mathdef.h"
#include "defs.h"


/**
   Create a new X(mat) matrix object, mmapped from file. The file is truncated if already exists in rw mode.*/
X(mat)* X(new_mmap)(long nx, long ny, const char *header, const char *format, ...){
    if(!nx || !ny) return NULL;
    if(disable_save){
	return X(new)(nx, ny);
    }
    format2fn;
    int fd=mmap_open(fn, 1);
    size_t metasize=3*8+bytes_header(header);//size of meta data. 
    size_t msize=nx*ny*sizeof(T)+metasize;//total size of file/memory
    if(ftruncate(fd, msize)){
	error("Error truncating file %s to size %zu\n", fn, msize);
    }
    char *map=(char*)mmap(NULL, msize, (PROT_WRITE|PROT_READ), MAP_SHARED, fd, 0);
    if(map==MAP_FAILED){
	perror("mmap");
	error("mmap failed\n");
    }
    char *map0=map;//save value of map
    /*memset(map, 0, msize); */
    const char *header0;
    mmap_header_rw(&map, &header0, M_T, nx, ny, header);
    mem_t *mem=mem_new(fd, map0, msize);
    X(mat) *out=X(new_do)(nx, ny, (T*)map, mem);
    mem_unref(mem);
    if(header) out->header=strdup(header);
    return out;
}
/**
   Create a new X(cell) matrix cell object, mmapped from file. The file is
   truncated if already exists. 
*/
X(cell)* X(cellnew_mmap)(long nx, long ny, long *nnx, long *nny, 
			 const char *header1, const char *header2[],
			 const char *format, ...){
    if(!nx || !ny) return NULL;
    if(disable_save){
	return X(cellnew3)(nx, ny, nnx, nny);
    }
    format2fn;
    int fd=mmap_open(fn, 1);
    long metasize=3*8;
    long msize=metasize+bytes_header(header1);
    for(long ix=0; ix<nx*ny; ix++){
	long mx=(long)nnx>0?nnx[ix]:(-(long)nnx);
	long my=nny?((long)nny>0?nny[ix]:(-(long)nny)):1;
	msize+=metasize+mx*my*sizeof(T)+bytes_header(header2?header2[ix]:NULL);
    }
    if(ftruncate(fd, msize)){
	error("Error truncating file %s to size %lu\n", fn, msize);
    }
    char *map=(char*)mmap(NULL, msize, (PROT_WRITE|PROT_READ), MAP_SHARED, fd, 0);
    if(map==MAP_FAILED){
	perror("mmap");
	error("mmap failed\n");
    }
    char *map0=map;
    /*memset(map, 0, msize); */
    const char *header0;
    mmap_header_rw(&map, &header0, MCC_ANY, nx, ny, header1);
    X(cell) *out=X(cellnew)(nx,ny);
    if(header0) out->header=strdup(header0);
    mem_t *mem=mem_new(fd, map0, msize);
    for(long ix=0; ix<nx*ny; ix++){
	long mx=(long)nnx>0?nnx[ix]:(-(long)nnx);
	long my=nny?((long)nny>0?nny[ix]:(-(long)nny)):1;
	mmap_header_rw(&map, &header0, M_T, mx, my, header2?header2[ix]:NULL);
	out->p[ix]=X(new_do)(mx, my, (T*)map, mem);
	memset(map, 0, mx*my*sizeof(T));//temporary
	map+=nnx[ix]*my*sizeof(T);
	if(header0) out->p[ix]->header=strdup(header0);
    }
    
    return out;
}
/**
   Create a new X(cell) matrix cell object, with identical blocks, mmapped from
   file. be aware that the data is not 8-byte aligned. The file is truncated if
   already exists.  */
X(cell)* X(cellnewsame_mmap)(long nx, long ny, long mx, long my, const char *header,
			     const char *format, ...){
    format2fn;
    return X(cellnew_mmap)(nx, ny, (long*)-mx, (long*)-my, header, NULL, "%s", fn);
}

static X(mat*) X(readdata_mmap)(char **map, mem_t *mem){
    long nx, ny;
    uint32_t magic;
    const char *header;
    mmap_header_ro(map, &magic, &nx, &ny, &header);
    if(magic!=M_T){
	error("File has magic %d, we want %d\n", (int)magic, M_T);
    }
    X(mat)* out=X(new_do)(nx, ny, (T*)(*map), mem);
    *map+=nx*ny*sizeof(T);
    if(out){
	out->header=strdup(header);
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
    char *map=(char*)mmap(NULL, msize, PROT_READ, MAP_SHARED, fd, 0);
    if(map==MAP_FAILED){	
	perror("mmap");
	error("mmap failed\n");
    }
    char *map0=map;
    X(mat) *out=X(readdata_mmap)(&map, mem_new(fd, map0,msize));
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
    char *map=(char*)mmap(NULL, msize, PROT_READ, MAP_SHARED, fd, 0);
    if(map==MAP_FAILED){
	perror("mmap");
	error("mmap failed\n");
    }
    char *map0=map;
    long nx, ny;
    uint32_t magic;
    const char *header;
    mmap_header_ro(&map, &magic, &nx, &ny, &header);
    if(!iscell(&magic)){
	error("We want a cell array, File has %x\n", (int)magic);
    }
    X(cell) *out=X(cellnew)(nx, ny);
    mem_t *mem=mem_new(fd, map0, msize);
    if(header) out->header=strdup(header);
    for(long ix=0; ix<nx*ny; ix++){
	out->p[ix]=X(readdata_mmap)(&map, mem);
    }
    return out;
}
