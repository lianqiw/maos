/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   Defines arrays of integers. We don't use the template mat.c because we only
   need very few functions for this data type and no numerical functions are
   required.
*/
#include <stdlib.h>
#include "imat.h"
imat* inew(long nx, long ny){
    imat *A=calloc(1, sizeof(imat));
    A->nx=nx;
    A->ny=ny;
    A->p=calloc(nx*ny, sizeof(long));
    return A;
}

icell* icellnew(long nx, long ny){
    icell *A=calloc(1, sizeof(icell));
    A->nx=nx;
    A->ny=ny;
    A->p=calloc(nx*ny, sizeof(imat*));
    return A;
}

void ifree(imat *A){
    if(A) {
	free(A->p);
	free(A);
    }
}

void icellfree(icell *A){
    if(A){
	for(long i=0; i<A->nx*A->ny; i++){
	    ifree(A->p[i]);
	}
	free(A);
    }
}
