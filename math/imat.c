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
    A->id=M_LONG;
    A->nx=nx;
    A->ny=ny;
    A->p=calloc(nx*ny, sizeof(long));
    return A;
}
void ifree_do(imat *A, int keepdata){
    if(!A) return;
    if(!keepdata) free(A->p);
    free(A);
}
void iresize(imat *A, long nx, long ny){
    if(!A) return;
    A->p=realloc(A->p, sizeof(long)*nx*ny);
    if((nx*ny)>(A->nx*A->ny)){
	memset(A->p+A->nx*A->ny, 0, sizeof(long)*((nx*ny)-(A->nx*A->ny)));
    }
    A->nx=nx;
    A->ny=ny;
}

void iwritedata(file_t *fp, const imat *A){
    uint64_t nx=0, ny=0;
    if(A){
	nx=(uint64_t)A->nx;
	ny=(uint64_t)A->ny;
    }
    do_write(fp, 0, sizeof(long), M_LONG, NULL, A?A->p:NULL, nx, ny);
}
imat *ireaddata(file_t *fp, header_t *header){
    (void)fp; (void)header;
    error("Please do int conversion\n");
    return 0;
}

long isum(const imat *A){
    long res=0;
    for(long i=0; i<A->nx*A->ny; i++){
	res+=A->p[i];
    }
    return res;
}
