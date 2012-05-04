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
   Createa curmat object.
*/
#ifndef AOS_CUDA_CUMAT_H
#define AOS_CUDA_CUMAT_H
/*
template <typename T>
inline cumat<T> *curef(cumat<T> *A){
    if(!A) return NULL;
    cumat<T> *out=(cumat<T>*)calloc(1, sizeof(cumat<T>));
    memcpy(out, A, sizeof(cumat<T>));
    A->nref[0]++;
    return out;
}
*/
template <typename T>
inline void cuzero(cumat<T> *A, cudaStream_t stream){
    if(A && A->p){
	DO(cudaMemsetAsync(A->p, 0, A->nx*A->ny*sizeof(T), stream));
    }
}

template <typename T>
inline void cuzero(cumat<T> *A){
    if(A && A->p){
	DO(cudaMemset(A->p, 0, A->nx*A->ny*sizeof(T)));
    }
}

/*
  The following are for cell.
*/

template <typename T>
inline cucell<T>* cucellnew(const int nx, const int ny){
    cucell<T> *out=(cucell<T>*)calloc(1, sizeof(cucell<T>));
    out->p=(cumat<T>**)calloc(nx*ny, sizeof(cumat<T>*));
    out->nx=nx;
    out->ny=ny;
    return out;
}

/** Allocate continuous memory for blocks of the same size*/
template <typename T>
inline cucell<T>* cucellnew(const int nx, const int ny, int mx, int my){
    cucell<T> *out=cucellnew<T>(nx, ny);
    out->m=new cumat<T>(mx*my*nx*ny, 1);
    for(int i=0; i<nx*ny; i++){
	out->p[i]=mx&&my?new cumat<T>(mx, my, out->m->p+i*(mx*my), 0):NULL;
    }
    return out;
}

/**
   Create a cellarray, mx and my should be nx*ny long.
*/
template <typename T, typename L>
inline cucell<T>* cucellnew(const int nx, const int ny, L *mx, L *my){
    long tot=0;
    for(int i=0; i<nx*ny; i++){
	tot+=mx[i]*(my?my[i]:1);
    }
    cucell<T> *out=cucellnew<T>(nx,ny);
    out->m=new cumat<T> (tot,1);
    tot=0;
    for(int i=0; i<nx*ny; i++){
	out->p[i]=mx[i]?new cumat<T>(mx[i],(my?my[i]:1),out->m->p+tot, 0):NULL;
	tot+=mx[i]*(my?my[i]:1);
    }
    return out;
}

template <typename T>
inline cucell<T>* cucellnew(const cucell<T> *in){
    cucell<T> *out;
    if(!in->m){
	out=cucellnew<T>(in->nx, in->ny);
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=new cumat<T>(in->p[i]->nx, in->p[i]->ny);
	}
    }else{
	int mx[in->nx*in->ny];
	int my[in->nx*in->ny];
	for(int i=0; i<in->nx*in->ny; i++){
	    mx[i]=in->p[i]->nx;
	    my[i]=in->p[i]->ny;
	}
	out=cucellnew<T,int>(in->nx, in->ny, (int*)mx, (int*)my);
    }
    return out;
}

template <typename T>
inline void cucellfree(cucell<T> *A){
    delete A;
    /*if(!A) return;
      if(A->p){
	for(int i=0; i<A->nx*A->ny; i++){
	    delete A->p[i];
	}
	free(A->p);
    }
    delete A->m;
    free(A);*/
}
template <typename T>
inline void cucellzero(cucell<T> *A, cudaStream_t stream){
    if(!A) return;
    if(A->m){
	cuzero(A->m, stream);
    }else{
	for(int i=0; i<A->nx*A->ny; i++){
	    cuzero(A->p[i], stream);
	}
    }
}

template <typename T, uint32_t magic>
    inline void cuwritedata(const cumat<T> *A, file_t *fp){
    if(A && A->nx>0 && A->ny>0){
	T *tmp=(T*)malloc(A->nx*A->ny*sizeof(T));
	cudaMemcpy(tmp, A->p, A->nx*A->ny*sizeof(T), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	do_write(fp, 0, sizeof(T), magic, NULL, tmp, A->nx, A->ny);
	free(tmp);
    }else{
	do_write(fp, 0, sizeof(T), magic, NULL, NULL, 0, 0);
    }
}
template <typename T, uint32_t magic>
    inline void cuwrite(const cumat<T> *A, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    cuwritedata<T, magic>(A, fp);
    zfclose(fp);
}

template <typename T, uint32_t magic>
    inline void cucellwrite(const cucell<T> *A, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    header_t header={MCC_ANY, A?A->nx:0, A?A->ny:0, NULL};
    write_header(&header, fp);
    if(A){
	for(int i=0; i<A->nx*A->ny; i++){
	    cuwritedata<T, magic>(A->p[i], fp);
	}
    }
    zfclose(fp);	
}

template <typename T>
inline void cucellcp(cucell<T> **out, const cucell<T>*in, cudaStream_t stream){
    if(!*out){
	*out=cucellnew<T>(in);
    }
    if(!in->m){
	for(int i=0; i<in->nx*in->ny; i++){
	    curcp(&(*out)->p[i], in->p[i], stream);
	}
    }else{
	if(!(*out)->m){
	    error("in is continuous, out is not\n");
	}
	curcp(&(*out)->m, in->m, stream);
    }
}
#endif
