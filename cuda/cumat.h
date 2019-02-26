/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <typeinfo>
template <typename T>
static inline void cuwritedata(const Array<T, Gpu> &A, file_t *fp){
    uint32_t magic;
    if(typeid(T)==typeid(float)){
	magic=M_FLT;
    }else if(typeid(T)==typeid(double)){
	magic=M_DBL;
    }else if(typeid(T)==typeid(double2)){
	magic=M_CMP;
    }else if(typeid(T)==typeid(float2)){
	magic=M_ZMP;
    }else{
	error("Invalid type\n");
    }
    if(A && A.Nx()>0 && A.Ny()>0){
	T *tmp=(T*)malloc(A.Nx()*A.Ny()*sizeof(T));
	cudaMemcpy(tmp, A(), A.Nx()*A.Ny()*sizeof(T), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	writearr(fp, 0, sizeof(T), magic, NULL, tmp, A.Nx(), A.Ny());
	free(tmp);
    }else{
	writearr(fp, 0, sizeof(T), magic, NULL, NULL, 0, 0);
    }
}
template <typename T>
static inline void cuwrite(const Array<T, Gpu> &A, const char *format, ...)CHECK_ARG(2){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    cuwritedata<T>(A, fp);
    zfclose(fp);
}

template <typename T>
static inline void cuwrite(const cucell<T, Gpu> &A, const char *format, ...)CHECK_ARG(2){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    header_t header={MCC_ANY, A?(uint64_t)A.Nx():0, A?(uint64_t)A.Ny():0, NULL};
    write_header(&header, fp);
    if(A){
	for(int i=0; i<A.Nx()*A.Ny(); i++){
	    cuwritedata<T>(A[i], fp);
	}
    }
    zfclose(fp);	
}

template <typename T>
static inline void cucellcp(cucell<T, Gpu> &out, const cucell<T, Gpu>&in, cudaStream_t stream){
    if(!out){
	out=in.New();
    }
    if(!in.M()){
	for(int i=0; i<in.Nx()*in.Ny(); i++){
	    curcp(out[i], in[i], stream);
	}
    }else{
	if(!out.M()){
	    error("in is continuous, out is not\n");
	}
	curcp(out.M(), in.M(), stream);
    }
}
#endif
