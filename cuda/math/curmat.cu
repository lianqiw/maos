/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "curmat.h"
#include "utils.h"
#include "cublas.h"
#include "kernel.h"

/**
   out=out*beta+in*alpha;
*/
/*void Add(curmat& out, Real alpha, const curmat& in, Real beta, cudaStream_t stream){
	if(!in) return;
	if(!out||alpha==0){
		Copy(out, in, stream);
		if(Z(fabs)(beta-(Real)1)>EPS){
			scale_do<<<DIM(in.N(), 256), 0, stream>>>
				(out(), in.N(), beta);
		}
	} else{
		assert(in.N()==out.N());
		add_do<<<DIM(in.N(), 256), 0, stream>>>
			(out(), (Real *)NULL, alpha, in(), (Real *)NULL, beta, in.N());
	}
}*/
/**
   out=out*beta+abs2(in)*alpha;
*/
void curaddcabs2(curmat& out, Real alpha, const cucmat& in, Real beta, cudaStream_t stream){
	if(!out){
		out=curmat(in.Nx(), in.Ny());
	} else if(alpha==0){
		cuzero(out, stream);
	}
	addcabs2_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), alpha, in(), beta, in.N());
}
/**
   out=out+abs2(in)*alpha;
*/
void curaddcabs2(curmat &out, const cucmat &in, Real beta, cudaStream_t stream){
	if(!out){
		out=curmat(in.Nx(), in.Ny());
	}
	addcabs2_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), in(), beta, in.N());
}
void Scale(curmat &in, Real alpha, cudaStream_t stream){
	if(!in) return;
	if(alpha==0){
		in.Zero(stream);
	} else if(Z(fabs)(alpha-(Real)1)>EPS){
		scale_do<<<DIM(in.N(), 256), 0, stream>>>(in(), in.N(), alpha);
	}
}


/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmv(Real* c, Real alpha, const curmat& A, const Real* b, char trans, Real beta, stream_t& stream){
	CUBL(gemv)(stream.blas(), (trans=='t'||trans==1)?CUBLAS_OP_T:CUBLAS_OP_N, A.Nx(), A.Ny(), &beta, A(), A.Nx(), b, 1, &alpha, c, 1);
}
void curcellmm(curcell& C, Real alpha, const curcell& A, const curcell& B,
	const char trans[2], const Real beta, stream_t& stream){
	if(!A||!B) return;
	int ax, az;
	int nx, ny, nz;
	int bz, by;
	if(trans[0]=='n'||trans[0]=='N'){
		nx=A.Nx();
		ax=1; az=A.Nx();
		nz=A.Ny();
	} else{
		nx=A.Ny();
		az=1; ax=A.Nx();
		nz=A.Nx();
	}
	if(trans[1]=='n'||trans[0]=='N'){
		ny=B.Ny();
		bz=1; by=B.Nx();
		if(nz!=B.Nx()) error("mismatch\n");
	} else{
		ny=B.Nx();
		by=1; bz=B.Nx();
		if(nz!=B.Ny()) error("mismatch\n");
	}
	if(!C){
		C=curcell(nx, ny);
	} else{
		assert(C.Nx()==nx&&C.Ny()==ny);
		if(alpha==0){
			cuzero(C, stream);
		} else if(Z(fabs)(alpha-(Real)1)>EPS){
			curcellscale(C, alpha, stream);
		}
	}
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			for(int iz=0; iz<nz; iz++){
				if(A[ix*ax+iz*az]&&B[iz*bz+iy*by]){
					cugemm(C[ix+iy*nx], (Real)1., A[ix*ax+iz*az],
						B[iz*bz+iy*by], trans, beta, stream);
				}
			}
		}
	}
}


/*
  A=A*beta+B*alpha;
*/
void curcelladd(curcell& A, Real beta, const curcell& B, Real alpha, cudaStream_t stream){
	if(!B) return;
	if(!A){
		A=B.New();
	} else{
		assert(A.Nx()==B.Nx()&&A.Ny()==B.Ny());
	}
	if(A.M()&&B.M()){
		Add(A.M(), beta, B.M(), alpha, stream);
	} else{
		for(int i=0; i<B.Nx()*B.Ny(); i++){
			Add(A[i], beta, B[i], alpha, stream);
		}
	}
}
/*
void Add(curmat& A, Real beta, cudaStream_t stream){
	const int n=A.Nx()*A.Ny();
	add_do<<<DIM(n, 256), 0, stream>>>(A(), beta, n);
}*/
/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a+b*alpha*beta;
*/

/**
   out=out+in*alpha; beta, alpha lives in device memory.
*/
/*void Add(curmat& out, const curmat& in, Real* alpha, Real alpha2, cudaStream_t stream){
	if(!out){
		out=curmat(in.Nx(), in.Ny());
	}
	add_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), in(), alpha, alpha2, in.N());
}*/


/**
   A=A*beta+B*alpha; beta, alpha lives in device memory.
*/
void curcelladd(curcell& A, const curcell& B, Real* alpha, Real alpha2, cudaStream_t stream){
	if(!B) return;
	if(!A){
		A=B.New();
	} else{
		assert(A.Nx()==B.Nx()&&A.Ny()==B.Ny());
	}
	if(A.M()&&B.M()){
		Add(A.M(), B.M(), alpha, alpha2, stream);
	} else{
		for(int i=0; i<B.N(); i++){
			Add(A[i], B[i], alpha, alpha2, stream);
		}
	}
}

/**
   out=out*beta+in; beta, alpha lives in device memory.
*/
/*void Add(curmat& out, Real* alpha1, const curmat& in, cudaStream_t stream){
	if(!out){
		out=curmat(in.Nx(), in.Ny());
	}
	add_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), alpha1, 1.f, in(), in.N());
}*/

/**
   A=A*alpha1+B*alpha; alpha1, alpha lives in device memory.
*/
void curcelladd(curcell& A, Real* alpha1, const curcell& B, cudaStream_t stream){
	if(!B) return;
	if(!A){
		A=B.New();
	} else{
		assert(A.Nx()==B.Nx()&&A.Ny()==B.Ny());
	}
	if(A.M()&&B.M()){
		Add(A.M(), alpha1, B.M(), stream);
	} else{
		for(int i=0; i<B.Nx()*B.Ny(); i++){
			Add(A[i], alpha1, B[i], stream);
		}
	}
}



Real curinn(const curmat& a, const curmat& b, cudaStream_t stream){
	curmat res(1, 1);
	Real out;
	inn_wrap(res(), a(), b(), a.Nx()*a.Ny(), stream);
	DO(cudaMemcpyAsync(&out, res(), sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}

/**
   Sum all the elements in an array.
 */
void cursum2(Real* restrict res,/**<Result in GPU*/
	const curmat& a,   /**<Source in GPU*/
	cudaStream_t stream){
	DO(cudaMemsetAsync(res, 0, sizeof(Real), stream));
	sum_wrap(res, a(), a.Nx()*a.Ny(), stream);
}
/**
   Sum all the elements in an array, and return a value.
*/
Real cursum(const curmat& a, cudaStream_t stream){
	Real out;//result in CPU.
	curmat res(1, 1);
	sum_wrap(res, a(), a.Nx()*a.Ny(), stream);
	DO(cudaMemcpyAsync(&out, res(), sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}

/**
   Find the maximum value
*/
Real curmax(const curmat& a, cudaStream_t stream){
	Real out;
	curmat res(1, 1);
	max_wrap(res, a(), a.N(), stream);
	DO(cudaMemcpyAsync(&out, res(), sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}

/**
   Find the maximum value
*/
Real curmaxabs(const curmat& a, cudaStream_t stream){
	Real out;
	curmat res(1, 1);
	maxabs_wrap(res, a(), a.N(), stream);
	DO(cudaMemcpyAsync(&out, res(), sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}
/**
   Find the maximum value
*/
Real curcellmax(const curcell& a, cudaStream_t stream){
	int n=a.Nx()*a.Ny();
	Real out;
	curmat res(n+1, 1);
	for(int i=0; i<n; i++){
		int m=a[i].N();
		max_wrap(&res[i], a[i](), m, stream);
	}
	if(n>1){
		max_wrap(&res[n], res, n, stream);
	}
	DO(cudaMemcpyAsync(&out, &res[n>1?n:0], sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}
/**
   Find the maximum value
*/
Real curcellmaxabs(const curcell& a, cudaStream_t stream){
	int n=a.N();
	Real out;
	curmat res(n+1, 1);
	for(int i=0; i<n; i++){
		int m=a[i].N();
		maxabs_wrap(&res[i], a[i](), m, stream);
	}
	if(n>1){
		maxabs_wrap(&res[n], res, n, stream);
	}
	DO(cudaMemcpyAsync(&out, &res[n>1?n:0], sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	return out;
}
/**
   Scale elements
*/
void curcellscale(curcell& A, Real alpha, cudaStream_t stream){
	if(!A) return;
	if(A.M()){
		Scale(A.M(), alpha, stream);
	} else{
		for(int i=0; i<A.Nx()*A.Ny(); i++){
			Scale(A[i], alpha, stream);
		}
	}
}



void cucscale(cucmat& in, Real alpha, cudaStream_t stream){
	if(!in) return;
	if(alpha==0){
		cuzero(in, stream);
	} else if(Z(fabs)(alpha-1.f)>EPS){
		int n=in.N();
		scale_do<<<DIM(n, 256), 0, stream>>>(in(), n, alpha);
	}
}
