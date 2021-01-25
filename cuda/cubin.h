/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "types.h"
template <typename T>
struct Magic{};
template <>
struct Magic<double>{
	enum{ magic=M_DBL };
};
template <>
struct Magic<float>{
	enum{ magic=M_FLT };
};
template <>
struct Magic<double2>{
	enum{ magic=M_CMP };
};
template <>
struct Magic<float2>{
	enum{ magic=M_ZMP };
};
template <>
struct Magic<int>{//always 32 bit int
	enum{ magic=M_INT32 };
};
template <>
struct Magic<long int>{//always 64 bit int. long is not.
	enum{ magic=M_INT64 };
};

//C++ does not support partial specialization of functions, we use overloading.

template <typename T>
static inline void cuwritedata_do(const Array<T, Gpu>& A, file_t* fp, uint32_t magic, const char* header){
	T* tmp=(T*)malloc(A.N()*sizeof(T));
	cudaMemcpy(tmp, A(), A.N()*sizeof(T), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	writearr(fp, 0, sizeof(T), magic, header, tmp, A.Nx(), A.Ny());
	free(tmp);
}
template <typename T>
static inline void cuwritedata_do(const Array<T, Cpu>& A, file_t* fp, uint32_t magic, const char* header){
	writearr(fp, 0, sizeof(T), magic, header, A(), A.Nx(), A.Ny());
}
template <typename T>
static inline void cuwritedata_do(const Array<T, Pinned>& A, file_t* fp, uint32_t magic, const char* header){
	writearr(fp, 0, sizeof(T), magic, header, A(), A.Nx(), A.Ny());
}

template <typename T, template<typename> class Dev>
static inline void cuwritedata(const Array<T, Dev>& A, file_t* fp){
	uint32_t magic=Magic<T>::magic;
	const char* header=A.header.length()?A.header.c_str():NULL;
	if(A.N()>0){
		cuwritedata_do(A, fp, magic, header);
	} else{
		writearr(fp, 0, sizeof(T), magic, header, NULL, 0, 0);
	}
}
template <typename T, template<typename> class Dev>
static inline void cuwrite(const Array<T, Dev>& A, const char* format, ...)CHECK_ARG(2){
	format2fn;
	file_t* fp=zfopen(fn, "wb");
	cuwritedata<T, Dev>(A, fp);
	zfclose(fp);
}
template <typename T, template<typename> class Dev>
static inline void cuwrite(const Cell<T, Dev>& A, const char* format, ...)CHECK_ARG(2){
	format2fn;
	file_t* fp=zfopen(fn, "wb");
	header_t header={MCC_ANY, A?(uint64_t)A.Nx():0, A?(uint64_t)A.Ny():0, NULL};
	write_header(&header, fp);
	if(A){
		for(int i=0; i<A.Nx()*A.Ny(); i++){
			cuwritedata<T, Dev>(A[i], fp);
		}
	}
	zfclose(fp);
}
template <typename T>
static inline void cucellcp(Cell<T, Gpu>& out, const Cell<T, Gpu>& in, cudaStream_t stream){
	if(!out){
		out=in.New();
	}
	if(!in.M()){
		for(int i=0; i<in.Nx()*in.Ny(); i++){
			curcp(out[i], in[i], stream);
		}
	} else{
		if(!out.M()){
			error("in is continuous, out is not\n");
		}
		curcp(out.M(), in.M(), stream);
	}
}

#endif
