/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

template <typename T, uint32_t magic>
    inline void cuwritedata(const cumat<T> *A, file_t *fp){
    if(A && A->nx>0 && A->ny>0){
	T *tmp=(T*)malloc(A->nx*A->ny*sizeof(T));
	cudaMemcpy(tmp, A->p, A->nx*A->ny*sizeof(T), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	writearr(fp, 0, sizeof(T), magic, NULL, tmp, A->nx, A->ny);
	free(tmp);
    }else{
	writearr(fp, 0, sizeof(T), magic, NULL, NULL, 0, 0);
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
    header_t header={MCC_ANY, A?(uint64_t)A->nx:0, A?(uint64_t)A->ny:0, NULL};
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
	*out=new cucell<T>(in);
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
