/**
   Createa curmat object.
*/
#ifndef AOS_CUDA_CUMAT_H
#define AOS_CUDA_CUMAT_H

template <typename T>
inline cumat<T>* cunew(int nx, int ny){
    cumat<T> *out;
    out=(cumat<T>*)calloc(1, sizeof(cumat<T>));
    out->nref=(int*)calloc(1, sizeof(int));
    out->nref[0]=1;
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(T)));
    DO(cudaMemset(out->p, 0, nx*ny*sizeof(T)));
    out->nx=nx;
    out->ny=ny;
    return out;
}
template <typename T>
inline cumat<T> *cunew(int nx, int ny, T *p, int own){
    cumat<T> *out;
    out=(cumat<T>*)calloc(1, sizeof(cumat<T>));
    if(own){
	out->nref=(int*)calloc(1, sizeof(int));
	out->nref[0]=1;
    }
    out->p=p;
    out->nx=nx;
    out->ny=ny;
    return out;
}
template <typename T>
inline cumat<T> *cunew(int nx, int ny, cudaStream_t stream){
    cumat<T> *out;
    out=(cumat<T>*)calloc(1, sizeof(cumat<T>));
    out->nref=(int*)calloc(1, sizeof(int));
    out->nref[0]=1;
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(T)));
    DO(cudaMemsetAsync(out->p, 0, nx*ny*sizeof(T), stream));
    out->nx=nx;
    out->ny=ny;
    return out;
}
template <typename T>
inline cumat<T> *curef(cumat<T> *A){
    if(!A) return NULL;
    cumat<T> *out=(cumat<T>*)calloc(1, sizeof(cumat<T>));
    memcpy(out, A, sizeof(cumat<T>));
    A->nref[0]++;
    return out;
}

template <typename T>
inline void cufree(cumat<T> *A){
    if(A){
	if(A->nref){
	    if(A->nref[0]==1){
		cudaFree(A->p);
		free(A->nref);
		free(A->header);
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
template <typename T>
inline void cuzero(cumat<T> *A, cudaStream_t stream){
    if(A && A->p){
	DO(cudaMemsetAsync(A->p, 0, A->nx*A->ny*sizeof(T), stream));
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
    out->m=cunew<T>(mx*my*nx*ny, 1);
    for(int i=0; i<nx*ny; i++){
	out->p[i]=cunew<T>(mx, my, out->m->p+i*(mx*my), 0);
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
	tot+=mx[i]*my[i];
    }
    cucell<T> *out=cucellnew<T>(nx,ny);
    out->m=cunew<T> (tot,1);
    tot=0;
    for(int i=0; i<nx*ny; i++){
	out->p[i]=cunew<T>(mx[i],my[i],out->m->p+tot, 0);
	tot+=mx[i]*my[i];
    }
    return out;
}

template <typename T>
inline cucell<T>* cucellnew(const cucell<T> *in){
    cucell<T> *out;
    if(!in->m){
	out=cucellnew<T>(in->nx, in->ny);
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=cunew<T>(in->p[i]->nx, in->p[i]->ny);
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
    if(!A) return;
    if(A->p){
	for(int i=0; i<A->nx*A->ny; i++){
	    cufree(A->p[i]);
	}
	free(A->p);
    }
    cufree(A->m);
    free(A);
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
