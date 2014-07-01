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
/**
   This file contains routines for imat (matrix of integer numbers).
*/
/*
  Defines arrays of integers. We don't use the template mat.c because we only
  need very few functions for this data type and no numerical functions are
   required.
*/
#include <stdlib.h>
#include "imat.h"
#include "../sys/bin.h"
/**
   Allocate a new imat.
 */
imat* inew(long nx, long ny){
    imat *A=calloc(1, sizeof(imat));
    A->nx=nx;
    A->ny=ny;
    A->p=calloc(nx*ny, sizeof(long));
    return A;
}
void iresize(imat *A, long nx, long ny){
    if(!A) return;
    A->p=realloc(A->p, sizeof(long)*nx*ny);
    A->nx=nx;
    A->ny=ny;
}
/**
   Allocate a new icell.
 */
icell* icellnew(long nx, long ny){
    icell *A=calloc(1, sizeof(icell));
    A->nx=nx;
    A->ny=ny;
    A->p=calloc(nx*ny, sizeof(imat*));
    return A;
}
/**
   Free a imat.
 */
void ifree(imat *A){
    if(A) {
	free(A->p);
	free(A);
    }
}
/**
   Free a icell
*/
void icellfree(icell *A){
    if(A){
	for(long i=0; i<A->nx*A->ny; i++){
	    ifree(A->p[i]);
	}
	free(A);
    }
}
void iwritedata(file_t *fp, const imat *A){
    uint64_t nx=0, ny=0;
    if(A){
	nx=(uint64_t)A->nx;
	ny=(uint64_t)A->ny;
    }
    do_write(fp, 0, sizeof(long), M_INT64, NULL, A?A->p:NULL, nx, ny);
}
/**
   Function to write cell array of dense matrix data. into a file pointer
   Generally used by library developer
*/
void icellwritedata(file_t *fp, const icell *dc){
    uint64_t nx=0;
    uint64_t ny=0;
    if(dc){
	nx=dc->nx;
	ny=dc->ny;
    }
    header_t header={MCC_ANY, nx, ny, NULL};
    write_header(&header, fp);
    if(nx>0 && ny>0){
	for(unsigned long iy=0; iy<ny; iy++){
	    for(unsigned long ix=0; ix<nx; ix++){
		iwritedata(fp, dc->p[ix+iy*nx]);
	    }
	}
    }
}
/**
   User callable function to write dense matrix into a file. Usage:
   iwrite(A,"A") for double matrix.
*/
void iwrite(const imat *A, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    iwritedata(fp, A);
    zfclose(fp);
}

long isum(const imat *A){
    long res=0;
    for(long i=0; i<A->nx*A->ny; i++){
	res+=A->p[i];
    }
    return res;
}
