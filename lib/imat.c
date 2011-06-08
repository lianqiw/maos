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
    A->p=calloc(nx*ny, sizeof(int));
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
