#include "cucmat.h"
#include "utils.h"
/**
   Createa cucmat object.
*/
cucmat *cucnew(int nx, int ny){
    cucmat *out;
    out=(cucmat*)calloc(1, sizeof(cucmat));
    out->ref=0;
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(fcomplex)));
    DO(cudaMemset(out->p, 0, nx*ny*sizeof(fcomplex)));
    out->nx=nx;
    out->ny=ny;
    return out;
}
void cucfree(cucmat *A){
    if(A){
	if(A->p){
	    cudaFree(A->p);
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



void cuccellfree(cuccell *A){
    if(!A) return;
    if(A->p){
	for(int i=0; i<A->nx*A->ny; i++){
	    cucfree(A->p[i]);
	}
	free(A->p);
    }
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
