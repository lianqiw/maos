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
#include "cucmat.h"
#include "utils.h"
/**
   Create cucmat object.
*/
cucmat *cucnew(int nx, int ny){
    cucmat *out;
    out=(cucmat*)calloc(1, sizeof(cucmat));
    out->nref=(int*)calloc(1, sizeof(int));
    out->nref[0]=1;
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(fcomplex)));
    DO(cudaMemset(out->p, 0, nx*ny*sizeof(fcomplex)));
    out->nx=nx;
    out->ny=ny;
    return out;
}
/**
   Create cucmat object with existing vector. own=1: we own the pointer. own=0: we don't own the pointer.
*/
cucmat *cucnew(int nx, int ny, fcomplex *p, int own){
    cucmat *out;
    out=(cucmat*)calloc(1, sizeof(cucmat));
    if(own){
	out->nref=(int*)calloc(1, sizeof(int));
	out->nref[0]=1;
    }
    out->p=p;
    out->nx=nx;
    out->ny=ny;
    return out;
}
cucmat *cucnew(int nx, int ny, cudaStream_t stream){
    cucmat *out;
    out=(cucmat*)calloc(1, sizeof(cucmat));
    out->nref=(int*)calloc(1, sizeof(int));
    out->nref[0]=1;
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(fcomplex)));
    DO(cudaMemsetAsync(out->p, 0, nx*ny*sizeof(fcomplex),stream));
    out->nx=nx;
    out->ny=ny;
    return out;
}
cucmat *cucref(cucmat *A){
    if(!A) return NULL;
 
    cucmat *out=(cucmat*)calloc(1, sizeof(cucmat));
    memcpy(out, A, sizeof(cucmat));
    A->nref[0]++;
    return out;
}

void cucfree(cucmat *A){
    if(A){
	if(A->nref){
	    if(A->nref[0]==1){
		cudaFree(A->p);
		free(A->nref);
	    }else{
		A->nref[0]--;
		if(A->nref[0]<0){
		    error("Invalid nref=%d\n", A->nref[0]);
		}
	    }
	}
	free(A);
    }
}

cuccell* cuccellnew(int nx, int ny){
    cuccell *out=(cuccell*)calloc(1, sizeof(cuccell));
    out->p=(cucmat**)calloc(nx*ny, sizeof(void*));
    out->nx=nx;
    out->ny=ny;
    return out;
}
/** Allocate continuous memory for blocks of the same size*/
cuccell *cuccellnew(int nx, int ny, int mx, int my){
    cuccell *out=cuccellnew(nx, ny);
    out->m=cucnew(mx, my*nx*ny);
    for(int i=0; i<nx*ny; i++){
	out->p[i]=cucnew(mx, my, out->m->p+i*(mx*my), 0);
    }
    return out;
}
void cuccellfree(cuccell *A){
    if(!A) return;
    if(A->p){
	for(int i=0; i<A->nx*A->ny; i++){
	    cucfree(A->p[i]);
	}
	free(A->p);
    }
    cucfree(A->m);
    free(A);
}

void cucwritedata(const cucmat *A, file_t *fp){
    if(A && A->nx >0 && A->ny>0){
	cudaDeviceSynchronize();
	fcomplex *tmp=(fcomplex*)malloc(A->nx*A->ny*sizeof(fcomplex));
	cudaMemcpy(tmp, A->p, A->nx*A->ny*sizeof(fcomplex), cudaMemcpyDefault);
	cudaDeviceSynchronize();
	do_write(fp, 0, sizeof(fcomplex), M_ZMP, tmp, A->nx, A->ny);
	free(tmp);
    }else{
	do_write(fp, 0, sizeof(fcomplex), M_ZMP, NULL, 0, 0);
    }
}
void cucwrite(const cucmat *A, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    cucwritedata(A, fp);
    zfclose(fp);
}
