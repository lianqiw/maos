/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#include <unistd.h>
#include "mathdef.h"
#include "defs.h" //Defines T, X, etc 
#include "blas.h" //defines blas headers
/*
  Group operations related to BLAS/LAPACK to this file.
*/
//#if !defined(COMP_SINGLE) && !defined(COMP_COMPLEX)
//If an external declaration is marked weakandthat symbol does not exist during linking(possibly dynamic) the address of the symbol will evaluate to NULL.
void (*X(svd_ext))(X(mat) **U_, XR(mat) **S_, X(mat) **Vt_, const X(mat) *A_)=NULL;
void (*X(svd_pow_ext))(X(mat) *A_, real power, real thres)=NULL;
void (*X(gemm_ext))(X(mat) **out, const real beta, const X(mat) *A, const X(mat) *B, const char trans[2], const real alpha)=NULL;
//#endif
/**
   compute matrix product using blas dgemm.;
   C=beta*C+ alpha *trans(A)*trans(B); if C exist.
*/
void X(mm)(X(mat)** C0, const T beta, const X(mat)* A, const X(mat)* B,
	const char trans[2], const T alpha){
	if(!A||!B||A->nx==0||B->nx==0) return;
#if !defined(COMP_SINGLE) && !defined(COMP_COMPLEX)
	if(dgemm_ext&&NX(A)>500&&NY(A)>500&&NY(B)>100){
		dbg_once("Using external routine for dmm of size %ldx%ldx%ld\n", NX(A), NY(A), NY(B));
		dgemm_ext(C0,beta,A,B,trans,alpha);
		return;
	}
#endif
	ptrdiff_t m, n, k, lda, ldb, ldc, k2;
	if(trans[0]=='T'||trans[0]=='t'||trans[0]=='C'||trans[0]=='c'){
		m=A->ny; k=A->nx;
	} else{
		m=A->nx; k=A->ny;
	}
	if(trans[1]=='T'||trans[1]=='t'||trans[0]=='C'||trans[0]=='c'){
		n=B->nx;
		k2=B->ny;
	} else{
		n=B->ny;
		k2=B->nx;
	}
	if(k!=k2) error("dmm: Matrix doesn't match: A: %tdx%td, B: %tdx%td, trans=%s\n",
		m, k, k2, n, trans);
	X(init)(C0, m, n);
	X(mat)* C=*C0;
	lda=A->nx;
	ldb=B->nx;
	ldc=C->nx;
	Z(gemm)(&trans[0], &trans[1], &m, &n, &k, &alpha,
		P(A), &lda, P(B), &ldb, &beta, P(C), &ldc);
}


/**
   inplace invert a small square SPD matrix using lapack dposv_, usually
   (A'*w*A).  by solving Ax=I; copy x to A.  dposv_ modifies A also. be
   careful
*/
void X(invspd_inplace)(X(mat)* A){
	if(!A) return;
	if(A->nx!=A->ny) error("Must be a square matrix");
	ptrdiff_t info=0, N=A->nx;
	const char uplo='U';
	/* B is identity matrix */
	T* B=mycalloc(N*N, T);
	for(long i=0;i<N;i++){
		B[i+i*N]=1;
	}
	Z(posv)(&uplo, &N, &N, P(A), &N, B, &N, &info);
	if(info!=0){
		writebin(A, "posv");
		error("posv_ failed, info=%td. data saved to posv.\n", info);
	}
	memcpy(P(A), B, sizeof(T)*N*N);
	free(B);
}

/**
   out of place version of dinvspd_inplace
*/
X(mat)* X(invspd)(const X(mat)* A){
	if(!A) return NULL;
	X(mat)* out=NULL;
	X(cp)(&out, A);
	X(invspd_inplace)(out);
	return out;
}

/**
   inplace invert a general square matrix using lapack dgesv_
*/
void X(inv_inplace)(X(mat)* A){
	if(!A) return;
	if(A->nx!=A->ny) error("Must be a square matrix, but is %ldx%ld\n", A->nx, A->ny);
	ptrdiff_t info=0, N=A->nx;
	T* B=mycalloc(N*N, T);
	for(int i=0;i<N;i++){
		B[i+i*N]=1;
	}
	ptrdiff_t* ipiv=mycalloc(N, ptrdiff_t);
	Z(gesv)(&N, &N, P(A), &N, ipiv, B, &N, &info);
	if(info!=0){
		writebin(A, "gesv");
		error("dgesv_ failed, info=%td. data saved to posv.\n", info);
	}
	memcpy(P(A), B, sizeof(T)*N*N);
	free(B);
	free(ipiv);
}

/**
   out of place version of dinv
*/
X(mat)* X(inv)(const X(mat)* A){
	if(!A) return NULL;
	X(mat)* out=NULL;
	X(cp)(&out, A);
	X(inv_inplace)(out);
	return out;
}

/**
   compute the pseudo inverse of matrix A with weigthing of full matrix W or
   sparse matrix weighting Wsp.  For full matrix, wt can be either W or diag (W)
   for diagonal weighting.  B=inv(A'*W*A+tikcr)*A'*W;
   thres is the threshold to truncate eigenvalues.
*/
X(mat)* X(pinv2)(const X(mat)* A, const_anyarray W_, R thres){
	const cell* W=W_.c;
	if(!A) return NULL;
	X(mat)* AtW=NULL;
	/*Compute AtW=A'*W */
	if(W){
		if(ismat(W)){//dense
			const X(mat)* wt=(X(mat)*)W;
			if(wt->ny==wt->nx){
				X(mm)(&AtW, 0, A, wt, "tn", 1);
			} else if(wt->ny==1&&wt->nx==A->nx){
				AtW=X(new)(A->ny, A->nx);
				for(long iy=0; iy<A->ny; iy++){
					for(long ix=0;ix<A->nx; ix++){
						P(AtW, iy, ix)=P(A, ix, iy)*P(wt, ix);
					}
				}
			} else{
				error("Invalid format\n");
			}
		} else if(issp(W)){//sparse
			const X(sp)* Wsp=X(sp_cast)(W);
			X(mulsp)(&AtW, A, Wsp, "tn", 1);
		} else{
			error("Invalid data type :%u\n", W->id);
		}
	} else{
		AtW=X(trans)(A);
	}
	/*Compute cc=A'*W*A */
	X(mat)* cc=NULL;
	X(mm) (&cc, 0, AtW, A, "nn", 1);
	/*Compute inv of cc */
	/*X(invspd_inplace)(cc); */
	if(X(isnan(cc))){
		writebin(cc, "cc_isnan");
		writebin(A, "A_isnan");
		writebin(AtW, "AtW_isnan");
		writecell(W, "W_isnan");
	}
	X(svd_pow)(cc, -1, thres);/*invert the matrix using SVD. safe with small eigen values. */
	X(mat)* out=NULL;
	/*Compute (A'*W*A)*A'*W */
	X(mm) (&out, 0, cc, AtW, "nn", 1);
	X(free)(AtW);
	X(free)(cc);
	return out;
}

/**
   A convenient wrapper.
 */
X(mat)* X(pinv)(const X(mat)* A, const_anyarray W){
	return X(pinv2)(A, W, 1e-14);
}
/**
   computes out=out*alpha+exp(A*beta) using scaling and squaring method.
   Larger scaling is more accurate but slower. Tested against matlab expm

   First, scaling the matrix by 2^-n times so it is smaller than 1/threshold.
   Then compute exp(A*beta*2^-n) using Tylor expansion.
   Then compute the final answer by taking exponential of 2^n.
*/
void X(expm)(X(mat)** out, R alpha, const X(mat)* A, R beta){
	const int accuracy=10;//How many terms in Taylor expansion to evaluate
	const R threshold=10;
	X(mat)* m_small=0;
	X(mat)* m_exp1=0, * m_power=0, * m_power1=0;
	//first determine the scaling needed
	int scaling=0;
	{
		R norm=sqrt(X(norm)(A));
		R max2=X(maxabs)(A);
		if(norm<max2){
			norm=max2;
		}
		scaling=(int)ceil(log2(fabs(norm*beta*threshold)));
		if(scaling<0) scaling=0;
	}
	X(add)(&m_small, 0, A, beta*exp2(-scaling));
	X(mat)* result=X(new)(A->nx, A->ny);
	X(addI)(result, 1);
	X(cp)(&m_power, m_small);
	//Compute the exponential using Taylor expansion.
	R factorial_i=1.0;
	for(int i=1; i<accuracy; i++){
		factorial_i*=i;
		//m_exp += M_power/factorial(i)
		X(add)(&result, 1., m_power, 1./(factorial_i));
		//m_power *= m_small;
		X(mm)(&m_power1, 0, m_power, m_small, "nn", 1);
		X(cp)(&m_power, m_power1);
	}
	//squaring step
	for(int i=0; i<scaling; i++){
		X(cp)(&m_exp1, result);
		X(mm)(&result, 0, m_exp1, m_exp1, "nn", 1);
	}
	X(add)(out, alpha, result, 1);
	if(X(isnan)(*out)){
		static int count=-1; count++;
		dbg("scaling=%d. exp2=%g\n", scaling, exp2(-scaling));
		writebin(A, "error_expm_%d_%f", count, beta);
		error("expm returns NaN\n");
	}
	X(free)(m_small);
	X(free)(m_exp1);
	X(free)(m_power);
	X(free)(m_power1);
	X(free)(result);
}
/**
   compute inv(dmcc(A, wt))
*/
X(mat)* X(imcc)(const X(mat)* A, const X(mat)* wt){
	X(mat)* mcc=X(mcc)(A, wt);
	X(invspd_inplace)(mcc);
	return mcc;
}

/**
   Compute the cholesky decomposition of a symmetric semi-definit dense matrix.
*/
X(mat)* X(chol)(const X(mat)* A){
	if(!A) return NULL;
	if(A->nx!=A->ny) error("dchol requires square matrix\n");
	X(mat)* B=X(dup)(A);
	if(A->nx==1){
		P(B,0)=sqrt(P(B,0));
		return B;
	}
	ptrdiff_t n=B->nx;
	ptrdiff_t info=0;//some take 4 byte, some take 8 byte in 64 bit machine.
	Z(potrf)("L", &n, P(B), &n, &info);
	if(info){
		writebin(A, "error_chol_A");
		if(A->nx+A->ny<20){
			X(show)(A, "chol_A");
		}
		if(info<0){
			error("The %td-th parameter has an illegal value\n", -info);
		} else{
			error("The leading minor of order %td is not posite definite\n", info);
		}
	} else{/*Zero out the upper diagonal. For some reason they are not zero. */
		X(mat)* Bp=B;
		for(long iy=0; iy<A->ny; iy++){
			for(long ix=0; ix<iy; ix++){
				P(Bp, ix, iy)=0;
			}
		}
	}
	return B;
}

DEF_ENV_FLAG(SVD_SDD_THRES, 0);
/**
   Compute SVD of a general matrix A.
   A=U*diag(S)*V';
   diag(S) is returned.
*/
void X(svd)(X(mat)** U, XR(mat)** Sdiag, X(mat)** VT, const X(mat)* A){
	if(isempty(A)){
		if(U) *U=0;
		if(VT) *VT=0;
		if(Sdiag) *Sdiag=0;
		return;
	}
//#if !defined(COMP_SINGLE) && !defined(COMP_COMPLEX)
	if(NX(A)>500&&X(svd_ext)){
		dbg_once("Using external routine for dsvd of size %ldx%ld\n", NX(A), NY(A));
		X(svd_ext)(U, Sdiag, VT, A);
		return;
	} 
//#endif
	int fd=-1;
	if(A->nx>2048&&!OMP_IN_PARALLEL){
	//Prevent multiple processes class gesvd simultaneously.
		char fnlock[PATH_MAX];
		snprintf(fnlock, PATH_MAX, "%s/%s", TEMP, "svd");
		fd=lock_file(fnlock, 1);
	}
	ptrdiff_t M=(ptrdiff_t)A->nx;
	ptrdiff_t N=(ptrdiff_t)A->ny;
	/*if((Sdiag&&*Sdiag)||(U&&*U)||(VT&&*VT)){
	warning("Sdiag,U,VT should all be NULL. discard their value\n");
	}*/
	X(mat)* tmp=X(dup)(A);
	ptrdiff_t nsvd=M<N?M:N;
	XR(mat)* s=XR(new)(nsvd, 1);
	X(mat)* u=X(new)(M, nsvd);
	X(mat)* vt=X(new)(nsvd, N);
	ptrdiff_t lwork=-1;
	T work0[1];
	ptrdiff_t info=0;
	
	char jobuv='S';
	if(SVD_SDD_THRES&&M>SVD_SDD_THRES){
		ptrdiff_t nmax=M<N?N:M;
		ptrdiff_t* iwork=mymalloc(nmax*8, ptrdiff_t);
		info("(dgesdd)");
#ifdef COMP_COMPLEX
		R* rwork=0;
		Z(gesdd)(&jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work0, &lwork, rwork, iwork, &info);
#else
		Z(gesdd)(&jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work0, &lwork, iwork, &info);
#endif
		lwork=(ptrdiff_t)REAL(work0[0]);
		T* work1=mymalloc(lwork, T);
#ifdef COMP_COMPLEX
		rwork=mymalloc(nsvd*MAX(5*nsvd+7, 2*nmax+2*nsvd+1), R);
		Z(gesdd)(&jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work1, &lwork, rwork, iwork, &info);
		free(rwork);
#else
		Z(gesdd)(&jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work1, &lwork, iwork, &info);
#endif
		free(work1);
		free(iwork);
	} else{
#ifdef COMP_COMPLEX
		R* rwork=0;
		Z(gesvd)(&jobuv, &jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work0, &lwork, rwork, &info);
#else
		Z(gesvd)(&jobuv, &jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work0, &lwork, &info);
#endif
		lwork=(ptrdiff_t)REAL(work0[0]);
		T* work1=mymalloc(lwork, T);
#ifdef COMP_COMPLEX
		rwork=mymalloc(nsvd*5, R);
		Z(gesvd)(&jobuv, &jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work1, &lwork, rwork, &info);
		free(rwork);
#else
		Z(gesvd)(&jobuv, &jobuv, &M, &N, P(tmp), &M, P(s), P(u), &M, P(vt), &nsvd, work1, &lwork, &info);
#endif
		free(work1);
	}
	if(info){
		writebin(A, "A_svd_failed");
		if(info<0){
			error("The %td-th argument has an illegal value\n", info);
		} else{
			error("svd: dbdsqr doesn't converge. info is %td\n", info);
		}
	}
	if(Sdiag) *Sdiag=s; else XR(free)(s);
	if(U) *U=u; else X(free)(u);
	if(VT) *VT=vt; else X(free)(vt);
	X(free)(tmp);
	if(fd>=0){
		close(fd);
	}
}
/**
   Compute SVD with caching.
 */
void X(svd_cache)(X(mat)** U, XR(mat)** Sdiag, X(mat)** VT, const X(mat)* A){
	char fnsvd[PATH_MAX];
	int do_cache=0;
	if(A->nx>512){
		//Cache the result
		uint32_t key=0;
		key=hashlittle(P(A), A->nx*A->ny*sizeof(T), key);
		snprintf(fnsvd, sizeof fnsvd, "%s/.aos/cache/svd", HOME);
		if(!exist(fnsvd)){
			mymkdir("%s", fnsvd);
		}
		long avail=available_space(fnsvd);
		long need=A->nx*A->ny*sizeof(T)*3;
		if(avail>need){
			do_cache=1;
			snprintf(fnsvd, sizeof fnsvd, "%s/.aos/cache/svd/svd_%ld_%u.bin", HOME, A->nx, key);
		}
	}
	if(!do_cache){
		X(svd)(U, Sdiag, VT, A);
	} else{
		cell* in=0;
		CACHE_FILE(in, fnsvd, ({in=readbin("%s", fnsvd);}), 
				({X(svd)(U, Sdiag, VT, A);
					in=cellnew(3, 1);
					P(in,0)=(cell *)*U;
					P(in,1)=(cell *)*Sdiag;
					P(in,2)=(cell *)*VT;}),
				({writecell(in, "%s", fnsvd);}));
		*U=X(mat_cast)(P(in,0)); P(in,0)=0;
		*Sdiag=XR(mat_cast)(P(in,1)); P(in,1)=0;
		*VT=X(mat_cast)(P(in,2)); P(in,2)=0;
		cellfree(in);
	}
}

/**
   computes pow(A,power) in place using svd.
   positive thres: Drop eigenvalues that are smaller than thres * max eigen value
   negative thres: Drop eigenvalues that are smaller than thres * previous eigen value (sorted descendantly).
*/
void X(svd_pow)(X(mat)* A, /**<[in/out] The matrix*/
	R power, /**<[in] power of eigen values. usually -1 for inverse*/ 
	R thres /**<[in] SVD inverse threshold*/
	){
	if(isempty(A)) return;
#if !defined(COMP_SINGLE) && !defined(COMP_COMPLEX)
	if(NX(A)>500&&dsvd_pow_ext){
		dbg_once("Using external routine for dsvd_pow of size %ldx%ld\n", NX(A), NY(A));
		dsvd_pow_ext(A, power, thres);
		return;
	}
#endif
	XR(mat)* Sdiag=NULL;
	X(mat)* U=NULL;
	X(mat)* VT=NULL;
	X(svd)(&U, &Sdiag, &VT, A);

	/*eigen values below the threshold will not be used. the first is the biggest. */
	R maxeig=fabs(P(Sdiag,0));
	R thres0=fabs(thres)*maxeig;
	for(long i=0; i<Sdiag->nx; i++){
		if(fabs(P(Sdiag,i))>thres0){/*only do with  */
			if(thres<0){/*compare adjacent eigenvalues*/
				thres0=P(Sdiag,i)*(-thres);
			}
			P(Sdiag,i)=pow(P(Sdiag,i), power);
		} else{
			for(int j=i; j<Sdiag->nx; j++){
				P(Sdiag,j)=0;
			}
			//long skipped=Sdiag->nx-i;
			break;
		}
	}
	for(long iy=0; iy<VT->ny; iy++){
		T* p=PCOL(VT, iy);
		for(long ix=0; ix<VT->nx; ix++){
			p[ix]*=P(Sdiag,ix);
		}
	}
#ifdef COMP_COMPLEX
	X(mm)(&A, 0, VT, U, "cc", 1);
#else
	X(mm)(&A, 0, VT, U, "tt", 1);
#endif
	X(free)(U);
	X(free)(VT);
	XR(free)(Sdiag);
}

/**
   Inplace Invert a SPD matrix. It is treated as a block matrix
*/
X(cell)* X(cellinvspd)(X(cell)* A){
	X(mat)* Ab=X(cell2m)(A);
	if(!Ab) return NULL;
	X(invspd_inplace)(Ab);
	X(cell)* B=NULL;
	X(2cell)(&B, Ab, A);
	X(celldropzero)(B, 0);
	X(free)(Ab);
	return B;
}

/**
   Inplace Invert a matrix. It is treated as a block matrix.
*/
X(cell)* X(cellinv)(X(cell)* A){
	X(mat)* Ab=X(cell2m)(A);
	X(inv_inplace)(Ab);
	X(cell)* B=NULL;
	X(2cell)(&B, Ab, A);
	X(celldropzero)(B, 0);
	X(free)(Ab);
	return B;
}
/**
   invert each component of the dcell. Each cell is treated
   as an individual matrix.
*/
X(cell)* X(cellinvspd_each)(X(cell)* A){
	X(cell)* out=NULL;
	X(cellcp)(&out, A);
	for(int i=0; i<out->nx*out->ny; i++){
		X(invspd_inplace)(P(out,i));
	}
	return out;
}

/**
   compute the pseudo inverse of block matrix A.  A is n*p cell, wt n*n cell or
   sparse cell.  \f$B=inv(A'*W*A)*A'*W\f$  */
X(cell)* X(cellpinv2)(const X(cell)* A, /**<[in] The matrix to pseudo invert*/
	const_anyarray W_,    /**<[in] The weighting matrix. dense or sparse*/
	R thres          /**<[in] SVD inverse threshold*/
	){
	if(!A) return NULL;
	X(cell)* wA=NULL;
	const cell *W=W_.c;
	if(W){
		if(NY(W)==NX(A)){
			X(cellmm)(&wA, W, A, "nn", 1);
		}else if(NX(A)==NX(W) && NY(A)==NY(W)){
			wA=X(celldup)(A);
			X(cellcwm)(wA, (X(cell)*)W);
		}else{
			warning("W is invalid, ignore it\n");
		}
	}
	if(!wA){
		wA=X(cellref)(A);
	}

	X(cell)* ata=NULL;
	X(cellmm)(&ata, wA, A, "tn", 1);
	X(cellsvd_pow)(ata, -1, thres);
	X(cell)* out=NULL;
	X(cellmm)(&out, ata, wA, "nt", 1);
	X(cellfree)(wA);
	X(cellfree)(ata);
	return out;
}
X(cell)* X(cellpinv)(const X(cell)* A, /**<[in] The matrix to pseudo invert*/
	const X(spcell)* W /**<[in] The weighting matrix. */){
	return X(cellpinv2)(A, W?W:NULL, 1e-14);
}

/**
   compute the power of a block matrix using svd method. First convert it do
   X(mat), do the power, and convert back to block matrix.
*/
void X(cellsvd_pow)(X(cell)* A, R power, R thres){
	X(mat)* Ac=X(cell2m)(A);
	X(svd_pow)(Ac, power, thres);
	X(2cell)(&A, Ac, NULL);
	X(celldropzero)(A, 0);
	X(free)(Ac);
}
