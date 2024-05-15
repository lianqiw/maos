/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
 * 
 * Wraps cuda cusparse, cusolver routines
 * */
#include "cublas.h"
#include "curmat.h"
#include "cucmat.h"
#include <cusolverDn.h>
int NULL_STREAM=0;
static cusparseMatDescr_t spdesc=NULL;
static __attribute((constructor)) void init(){
	DO(cusparseCreateMatDescr(&spdesc));
	cusparseSetMatType(spdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(spdesc, CUSPARSE_INDEX_BASE_ZERO);
	char *tmp=getenv("CUDA_LAUNCH_BLOCKING");
	if(tmp){
		int blocking=strtol(tmp, NULL, 10);
		if(blocking){
			warning("CUDA_LAUNCH_BLOCKING is enabled. Use only NULL stream\n");
			NULL_STREAM=1; //Use only default stream
		}
	}
}


static const char *scsrmv_err[]={
	"Success",
	"Not initialized",
	"Allocation failed",
	"Invalid value",
	"Archtecture mismatch",
	"Mapping error",
	"Execution failed",
	"Internal error",
	"Matrix type not supported"
};

/*
	y=A*x where A is sparse. x, y are vectors or matrices. (Slow for GS0.)
	Y = α op (A) ⋅ X + β Y or Y = α op (A) ⋅ X + β Y with β=1
	Converted cuSparse generic API for CUDA version 11 which removed old API.
*/

void cuspmul(Real *y, const cusp &A, const Real *x, int ncolvec, char trans, Real alpha, stream_t &stream){
	cusparseOperation_t opr;
	int istrans=(trans=='t'||trans==1);
	if(A.Type()==SP_CSC){
		istrans=!istrans;
	}
	if(istrans){
		opr=CUSPARSE_OPERATION_TRANSPOSE;
	} else{
		opr=CUSPARSE_OPERATION_NON_TRANSPOSE;
	}
	int ncol=0, nrow=0;
	switch(A.Type()){
	case SP_CSR:
		nrow=A.Nx(); ncol=A.Ny(); break;
	case SP_CSC:
		nrow=A.Ny(); ncol=A.Nx(); break;
	default:
		error("Invalid type: %d", A.Type());
	}
	int status;
	Real one=1.f;

#if CUDA_VERSION >= 10000
	/*
cusparseStatus_t
cusparseSpMV(cusparseHandle_t     handle,
					cusparseOperation_t  opA,
					const void*          alpha,
					cusparseSpMatDescr_t matA,
					cusparseDnVecDescr_t vecX,
					const void*          beta,
					cusparseDnVecDescr_t vecY,
					cudaDataType         computeType,
					cusparseSpMVAlg_t    alg,
					void*                externalBuffer)

cusparseStatus_t
cusparseSpMM(cusparseHandle_t     handle,
			cusparseOperation_t  opA,
			cusparseOperation_t  opB,
			const void*          alpha,
			cusparseSpMatDescr_t matA,
			cusparseDnMatDescr_t matB,
			const void*          beta,
			cusparseDnMatDescr_t matC,
			cudaDataType         computeType,
			cusparseSpMMAlg_t    alg,
			void*                externalBuffer)

	*/
	size_t bsize;
	void *buffer;
	size_t ny=istrans?ncol:nrow;
	size_t nx=istrans?nrow:ncol;
	if(ncolvec==1){
		cusparseDnVecDescr_t xv, yv;
		cusparseCreateDnVec(&xv, nx, (void *)x, CUDA_R);
		cusparseCreateDnVec(&yv, ny, (void *)y, CUDA_R);
#if CUDA_VERSION > 11010 
		cusparseSpMVAlg_t alg=CUSPARSE_SPMV_ALG_DEFAULT;
#else
		cusparseSpMVAlg_t alg=CUSPARSE_MV_ALG_DEFAULT;
#endif
		DO(status=cusparseSpMV_bufferSize(stream.sparse(), opr, &alpha, A.Desc(), xv, &one, yv, CUDA_R, alg, &bsize));
		DO(cudaMalloc(&buffer, bsize));
		DO(status=cusparseSpMV(stream.sparse(), opr, &alpha, A.Desc(), xv, &one, yv, CUDA_R, alg, buffer));
		DO(cudaFree(buffer));
		cusparseDestroyDnVec(yv);
		cusparseDestroyDnVec(xv);
	} else{
		cusparseDnMatDescr_t Bm, Cm;
		cusparseCreateDnMat(&Bm, nx, ncolvec, nx, (void *)x, CUDA_R, CUSPARSE_ORDER_COL);
		cusparseCreateDnMat(&Cm, ny, ncolvec, ny, (void *)y, CUDA_R, CUSPARSE_ORDER_COL);
		cusparseOperation_t opB=CUSPARSE_OPERATION_NON_TRANSPOSE;
#if CUDA_VERSION > 11010
		cusparseSpMMAlg_t alg=CUSPARSE_SPMM_ALG_DEFAULT;
#else
		cusparseSpMMAlg_t alg=CUSPARSE_MM_ALG_DEFAULT;
#endif
		DO(status=cusparseSpMM_bufferSize(stream.sparse(), opr, opB, &alpha,
										  A.Desc(), Bm, &one, Cm, CUDA_R, alg, &bsize));
		DO(cudaMalloc(&buffer, bsize));

		DO(status=cusparseSpMM(stream.sparse(), opr, opB, &alpha,
							   A.Desc(), Bm, &one, Cm, CUDA_R, alg, buffer));
		DO(cudaFree(buffer));
		cusparseDestroyDnMat(Bm);
		cusparseDestroyDnMat(Cm);
	}
#else		
	if(ncolvec==1){
		status=CUSP(csrmv)(stream.sparse(), opr,
			nrow, ncol, A.Nzmax(), &alpha, spdesc,
			A.Px(), A.Pp(), A.Pi(), x, &one, y);

	} else{
		int nleadx=istrans?nrow:ncol;
		int nleady=istrans?ncol:nrow;
		status=CUSP(csrmm)(stream.sparse(), opr,
			nrow, ncolvec, ncol, A.Nzmax(), &alpha, spdesc,
			A.Px(), A.Pp(), A.Pi(), x, nleadx, &one, y, nleady);
	}
#endif
	if(status!=0){
		error("cuspmul failed with status '%s'\n", scsrmv_err[status]);
	}
}

template<>
cublasStatus_t cublasGdgmm(cublasHandle_t handle, cublasSideMode_t mode,
						  int m, int n,
						  const float *A, int lda,
						  const float *x, int incx,
						  float *C, int ldc){
	return cublasSdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc);
}
cublasStatus_t cublasGdgmm(cublasHandle_t handle, cublasSideMode_t mode,
						  int m, int n,
						  const double *A, int lda,
						  const double *x, int incx,
						  double *C, int ldc){
	return cublasDdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc);
}
template<>
cublasStatus_t cublasGgemm(cublasHandle_t handle,
						   cublasOperation_t transa, cublasOperation_t transb,
						   int m, int n, int k,
						   const float *alpha,
						   const float *A, int lda,
						   const float *B, int ldb,
						   const float *beta,
						   float *C, int ldc){
	return cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
template<>
cublasStatus_t cublasGgemm(cublasHandle_t handle,
						   cublasOperation_t transa, cublasOperation_t transb,
						   int m, int n, int k,
						   const double *alpha,
						   const double *A, int lda,
						   const double *B, int ldb,
						   const double *beta,
						   double *C, int ldc){
	return cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
