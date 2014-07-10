/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "mathdef.h"
#include "defs.h"
/**
   Contains routines to write/read dense/sparse matrix into/from file.
*/

/**
   Function to write dense matrix data into a file pointer. Generally used by
   library developer */
void X(writedata)(file_t *fp, const X(mat) *A){
    uint64_t nx=0, ny=0;
    if(A){
	nx=(uint64_t)A->nx;
	ny=(uint64_t)A->ny;
    }
    do_write(fp, 0, sizeof(T), M_T, A?A->header:NULL, A?A->p:NULL, nx, ny);
}
/**
   Function to write cell array of dense matrix data. into a file pointer
   Generally used by library developer
*/
void X(cellwritedata)(file_t *fp, const X(cell) *dc){
    uint64_t nx=0;
    uint64_t ny=0;
    if(dc){
	nx=dc->nx;
	ny=dc->ny;
    }
    header_t header={MCC_ANY, nx, ny, dc?dc->header:NULL};
    write_header(&header, fp);
    if(nx>0 && ny>0){
	for(long iy=0; iy<ny; iy++){
	    for(long ix=0; ix<nx; ix++){
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
   file. Usage: dcellwrite(A,"A") for double matrix cell. */
void X(cellwrite)(const X(cell) *dc, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    X(cellwritedata)(fp,dc);
    zfclose(fp);
}
/**
   Function to read dense matrix into memory from file pointer. Generally used
   by library developer.  */
X(mat) *X(readdata)(file_t *fp, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    uint64_t nx, ny;
    nx=header->nx;
    ny=header->ny;
    if(!nx || !ny) return NULL;
    X(mat) *out;
    out=X(new)((long)nx,(long)ny);
    out->header=header->str; header->str=NULL;
    if(nx>0 && ny>0){
	if(header->magic==M_T){
	    zfread(out->p,sizeof(T),nx*ny,fp);
	}else if(M_T==M_DBL && header->magic==M_FLT){
	    float *p=malloc(nx*ny*sizeof(float));
	    zfread(p, sizeof(float), nx*ny, fp);
	    for(int i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
	}else if(M_T==M_CMP && header->magic==M_ZMP){
	    fcomplex *p=malloc(nx*ny*sizeof(fcomplex));
	    zfread(p, sizeof(fcomplex), nx*ny, fp);
	    for(int i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
	}else{
	    error("%s is not a X(mat) file. magic=%x. We want %x\n", 
		  zfname(fp), header->magic, M_T);
	}
    }
    return out;
}
/**
   Function to read dense matrix cell array into memory from file
   pointer. Generally used by library developer.  */
X(cell)* X(cellreaddata)(file_t *fp, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    long nx,ny;
    header_t *headerc=check_cell(header, &nx, &ny);
    X(cell) *out;
    if(!headerc || !zfisfits(fp)){/*genuine cell array or bin data.*/
	out=X(cellnew)((long)nx,(long)ny);
	if(!headerc){
	    out->header=header->str; header->str=NULL;
	}
        for(long ix=0; ix<nx*ny; ix++){
	    out->p[ix]=X(readdata)(fp, headerc);
	}
    }else{/*fits format. need to read extensions*/
	if(ny!=1) error("Not expected\n");
	X(mat) **tmp=calloc(nx, sizeof(X(mat)*));
	tmp[nx-1]=X(readdata)(fp, headerc);
	while(!read_header2(headerc, fp)){
	    nx++;
	    tmp=realloc(tmp, nx*sizeof(X(mat)*));
	    tmp[nx-1]=X(readdata)(fp, headerc);
	}
	out=X(cellnew)(nx, 1);
	memcpy(out->p, tmp, sizeof(X(mat)*)*nx);
	free(tmp);
    }
    return out;
}
/**
   User callable function to read dense matrix into memory from file. Usage:
   A=dread("A"); for a dmat*/
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

   Usage: A=dcellread("A"); for a double dcell.
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
    header_t header;
    read_header(&header, fp);
    long nx, ny;
    header_t *headerc=check_cell(&header, &nx, &ny);
    if(!headerc){/*Is cell arr*/
	free(header.str); header.str=NULL;
    }
    X(cell) **out=calloc(nx*ny, sizeof(X(cell)*));
    if(!headerc || !zfisfits(fp)){/*read cellarr or read cell as single cell arr.*/
	for(long ic=0; ic<nx*ny; ic++){
	    out[ic]=X(cellreaddata)(fp, headerc);
	}
    }else{/*fits format. need to read extensions*/
	out[0]=X(cellreaddata)(fp, headerc);
	if(ny!=1) error("Not expected\n");
	while(!read_header2(headerc, fp)){
	    nx++;
	    out=realloc(out, nx*sizeof(X(cell)*));
	    out[nx-1]=X(cellreaddata)(fp, headerc);
	}
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
    header_t header={MCC_ANY, nxin, nyin, NULL};
    write_header(&header, fp);
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

#include "matmmap.c"
