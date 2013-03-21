/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#define MAT_VERBOSE 0
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>

#include "../sys/sys.h"
#include "random.h"

#include "mathmisc.h"
#include "dsp.h"
#include "ssp.h"
#include "csp.h"
#include "dmat.h"
#include "smat.h"
#include "cmat.h"
#include "fft.h"
#include "matbin.h"
#include "loc.h"
#include "defs.h"/*Defines T, X, etc */

/*
  Group operations related to BLAS/LAPACK to this file.
*/
/**
   compute matrix product using blas dgemm with beta=1;
   C=beta*C+ alpha *trans(A)*trans(B); if C exist.
*/
void X(mm)(X(mat)**C0, const X(mat) *A, const X(mat) *B,   
	   const char trans[2], const T alpha){
    if(!A || !B) return;
    int m,n,k,lda,ldb,ldc,k2;
    if (trans[0]=='T' || trans[0]=='t' || trans[0]=='C' || trans[0]=='c'){
	m=A->ny; k=A->nx;
    }else{
	m=A->nx; k=A->ny;
    }
    if (trans[1]=='T' || trans[1]=='t'|| trans[0]=='C' || trans[0]=='c'){
	n=B->nx;
	k2=B->ny;
    }else{
	n=B->ny;
	k2=B->nx;
    }
    if(k!=k2) error("dmm: Matrix doesn't match: A: %dx%d, B: %dx%d\n",
		    m, k, k2, n);
    if(!*C0){
	*C0=X(new)(m,n); 
    }else if(m!=(*C0)->nx || n!=(*C0)->ny){
	error("dmm: Matrix doesn't match: C: %dx%d, C0: %ldx%ld\n", 
	      m, n, (*C0)->nx, (*C0)->ny);
    }
    X(mat) *C=*C0;
    lda=A->nx;
    ldb=B->nx;
    ldc=C->nx;
    const T beta=1;
    Z(gemm)(&trans[0], &trans[1], &m,&n,&k,&alpha, 
	    A->p, &lda, B->p, &ldb, &beta, C->p,&ldc);
}


/**
   inplace invert a small square SPD matrix using lapack dposv_, usually
   (A'*w*A).  by solving Ax=I; copy x to A.  dposv_ modifies A also. be
   careful
*/
void X(invspd_inplace)(X(mat) *A){
    if(!A) return;
    if(A->nx!=A->ny) error("Must be a square matrix");
    int info=0, N=A->nx;
    const char uplo='U';
    /* B is identity matrix */
    T *B=calloc(N*N,sizeof(T));
    for(long i=0;i<N;i++)
	B[i+i*N]=1;
    Z(posv)(&uplo, &N, &N, A->p, &N, B, &N, &info);
    if(info!=0){
	X(write)(A,"posv");
	error("posv_ failed, info=%d. data saved to posv.\n",info);
    }
    memcpy(A->p, B, sizeof(T)*N*N);
    free(B);
}

/**
   out of place version of dinvspd_inplace
*/
X(mat)* X(invspd)(const X(mat) *A){
    if(!A) return NULL;
    X(mat) *out=NULL;
    X(cp)(&out, A);
    X(invspd_inplace)(out);
    return out;
}



/**
   inplace invert a general square matrix using lapack dgesv_
*/
void X(inv_inplace)(X(mat)*A){
    if(!A) return;
    if(A->nx!=A->ny) error("Must be a square matrix");
    int info=0, N=A->nx;
    T *B=calloc(N*N,sizeof(T));
    for(int i=0;i<N;i++){
	B[i+i*N]=1;
    }
    int *ipiv=calloc(N, sizeof(int));
    Z(gesv)(&N, &N, A->p, &N, ipiv, B, &N, &info);
    if(info!=0){
	X(write)(A,"gesv");
	error("dgesv_ failed, info=%d. data saved to posv.\n",info);
    }
    memcpy(A->p, B, sizeof(T)*N*N);
    free(B);
    free(ipiv);
}

/**
   out of place version of dinv
*/
X(mat)* X(inv)(const X(mat) *A){
    if(!A) return NULL;
    X(mat) *out=NULL;
    X(cp)(&out, A);
    X(inv_inplace)(out);
    return out;
}

/**
   compute the pseudo inverse of matrix A with weigthing of full matrix W or
   sparse matrix weighting Wsp.  For full matrix, wt can be either W or diag (W)
   for diagonal weighting.  B=inv(A'*W*A)*A'*W; */
X(mat) *X(pinv)(const X(mat) *A, const X(mat) *wt, const X(sp) *Wsp){
    if(!A) return NULL;
    X(mat) *AtW=NULL;
    /*Compute AtW=A'*W */
    if(wt){
	if(Wsp){
	    error("Both wt and Wsp are supplied. Not supported\n");
	}
	if(wt->ny==wt->nx){
	    X(mm)(&AtW, A, wt, "tn", 1);
	}else if(wt->ny==1){
	    AtW=X(new)(A->ny,A->nx);
	    PMAT(A,pA);
	    PMAT(AtW,pAtW);
	    T *w=wt->p;
	    for(long iy=0; iy<A->ny; iy++){
		for(long ix=0;ix<A->nx; ix++){
		    pAtW[ix][iy]=pA[iy][ix]*w[ix];
		}
	    }
	}else{
	    error("Invalid format\n");
	}
    }else{
	if(Wsp){
	    X(mat)*At = X(trans)(A);
	    X(mulsp)(&AtW, At, Wsp, 1);
	    X(free)(At);
	}else{
	    AtW=X(trans)(A);
	}
    }
    /*Compute cc=A'*W*A */
    X(mat) *cc=NULL;
    X(mm) (&cc, AtW, A, "nn", 1);
    /*Compute inv of cc */
    /*X(invspd_inplace)(cc); */
    if(X(isnan(cc))){
	X(write)(cc,"cc_isnan");
	X(write)(A,"A_isnan");
	X(write)(wt,"wt_isnan");
	Y(spwrite)(Wsp, "Wsp_isnan");
    }
    X(svd_pow)(cc,-1,0,1e-14);/*invert the matrix using SVD. safe with small eigen values. */
    X(mat) *out=NULL;
    /*Compute (A'*W*A)*A'*W */
    X(mm) (&out, cc, AtW, "nn", 1);
    X(free)(AtW);
    X(free)(cc);
    return out;
}

/**
   compute inv(dmcc(A, wt))
*/
X(mat) *X(imcc)(const X(mat) *A, const X(mat) *wt){
    X(mat) *mcc=X(mcc)(A,wt);
    X(invspd_inplace)(mcc);
    return mcc;
}


/**
   Apply tikhonov regularization to A
*/
void X(tikcr)(X(mat) *A, T thres){
    XR(mat) *S=NULL;
    X(svd)(NULL,&S,NULL,A);
    T val=S->p[0]*thres;
    XR(free)(S);
    X(addI)(A,val);
}

/**
   Compute the cholesky decomposition of a symmetric semi-definit dense matrix.
*/
X(mat)* X(chol)(const X(mat) *A){
    if(!A) return NULL;
    if(A->nx!=A->ny) error("dchol requires square matrix\n");
    X(mat) *B = X(dup)(A);
    int uplo='L', n=B->nx, info;
    Z(potrf)(&uplo, &n, B->p, &n, &info);
    if(info){
	if(info<0){
	    error("The %d-th parameter has an illegal value\n", -info);
	}else{
	    error("The leading minor of order %d is not posite denifite\n", info);
	}
    }else{/*Zero out the upper diagonal. For some reason they are not zero. */
	PDMAT(B, Bp);
	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<iy; ix++){
		Bp[iy][ix]=0;
	    }
	}
    }
    return B;
}




/**
   Compute SVD of a general matrix A. 
   A=U*diag(S)*V';
   diag(S) is returned.
*/
void X(svd)(X(mat) **U, XR(mat) **Sdiag, X(mat) **VT, const X(mat) *A){
    char jobuv='S';
    int M=(int)A->nx;
    int N=(int)A->ny;
    /*if((Sdiag&&*Sdiag)||(U&&*U)||(VT&&*VT)){
	warning("Sdiag,U,VT should all be NULL. discard their value\n");
	}*/
    X(mat) *tmp=X(dup)(A);
    int nsvd=M<N?M:N;
    XR(mat) *s=XR(new)(nsvd,1);
    X(mat) *u=X(new)(M,nsvd);
    X(mat) *vt=X(new)(nsvd,N);
    int lwork=-1;
    T work0[1];
    int info=0;
#ifdef USE_COMPLEX
    R *rwork=malloc(nsvd*5*sizeof(R));
    Z(gesvd)(&jobuv,&jobuv,&M,&N,tmp->p,&M,s->p,u->p,&M,vt->p,&nsvd,work0,&lwork,rwork,&info);
#else
    Z(gesvd)(&jobuv,&jobuv,&M,&N,tmp->p,&M,s->p,u->p,&M,vt->p,&nsvd,work0,&lwork,&info);
#endif
    lwork=(int)(work0[0]);
    T *work1=malloc(sizeof(T)*lwork);
#ifdef USE_COMPLEX
    Z(gesvd)(&jobuv,&jobuv,&M,&N,tmp->p,&M,s->p,u->p,&M,vt->p,&nsvd,work1,&lwork,rwork,&info);
#else
    Z(gesvd)(&jobuv,&jobuv,&M,&N,tmp->p,&M,s->p,u->p,&M,vt->p,&nsvd,work1,&lwork,&info);
#endif
    free(work1);
    if(info){
	X(write)(A,"A_svd_failed");
	if(info<0){
	    error("The %d-th argument has an illegal value\n",info);
	}else{
	    error("svd: dbdsqr doesn't converge. info is %d\n",info);
	}
    }
    if(Sdiag) *Sdiag=s; else XR(free)(s);
    if(U) *U=u; else X(free)(u);
    if(VT) *VT=vt; else X(free)(vt);
    X(free)(tmp);
#ifdef USE_COMPLEX
    free(rwork);
#endif
}

/**
   Compute the eigen values and, optionally, eigen vectors of a real symmetric
matrix. Notice that here Sdiag is in ascending order, which is different from
X(svd).  */
void X(evd)(X(mat) **U, XR(mat) **Sdiag,const X(mat) *A){
    assert(A->nx==A->ny && A->nx>0);
    *Sdiag=XR(new)(A->nx,1);
    char jobz=U?'V':'N';
    char uplo='U';
    int lda=A->nx;
    T worksize;
    int lwork=-1;
    int info;
    X(mat) *atmp=X(dup)(A);
#ifdef USE_COMPLEX
    R *rwork=malloc((3*A->nx-2)*sizeof(R));
    Z(heev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, &worksize, &lwork,rwork, &info);
    lwork=(int)worksize;
    T *work=malloc(sizeof(T)*lwork);
    Z(heev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, work, &lwork,rwork, &info);
    free(rwork);
#else
    Z(syev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, &worksize, &lwork, &info);
    lwork=(int)worksize;
    R *work=malloc(sizeof(R)*lwork);
    Z(syev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, work, &lwork, &info);
#endif
    if(info){
	X(write)(A,"A_evd_failed");
	if(info<0)
	    error("The %dth argument had an illegal value\n", -info);
	else
	    error("The %dth number of elements did not converge to zero\n", info);
    }
    if(U) *U=atmp; else X(free)(atmp);
    free(work);
}

/**
   computes pow(A,power) in place using svd or evd. if issym==1, use evd, otherwise use svd
   positive thres: Drop eigenvalues that are smaller than thres * max eigen value
   negative thres: Drop eigenvalues that are smaller than thres * previous eigen value (sorted descendantly).
*/
void X(svd_pow)(X(mat) *A, double power, int issym, double thres){
    int use_evd=0;
    if(A->nx!=A->ny){
	warning("dpow is only good for square arrays.\n");
    }else if(issym){
	use_evd=1;
    }
    XR(mat) *Sdiag=NULL;
    X(mat) *U=NULL;
    X(mat) *VT=NULL;
    double maxeig;
    if(use_evd){
	X(evd)(&U, &Sdiag, A);
	/*eigen values below the threshold will not be used. the last is the biggest. */
	maxeig=FABS(Sdiag->p[Sdiag->nx-1]);
	VT=X(trans)(U);
    }else{
	X(svd)(&U, &Sdiag, &VT, A);
	/*eigen values below the threshold will not be used. the first is the biggest. */
	maxeig=FABS(Sdiag->p[0]);
    }
    double thres0=fabs(thres)*maxeig;
    for(long i=0; i<Sdiag->nx; i++){
	if(FABS(Sdiag->p[i])>thres0){/*only do with  */
	    if(thres<0){/*compare adjacent eigenvalues*/
		thres0=Sdiag->p[i]*(-thres);
	    }
	    Sdiag->p[i]=pow(Sdiag->p[i],power);
	}else{
	    for(int j=i; j<Sdiag->nx; j++){
		Sdiag->p[j]=0;
	    }
	    //long skipped=Sdiag->nx-i;
	    break;
	}
    }
    for(long iy=0; iy <VT->ny; iy++){
	T *p=VT->p+iy*VT->nx;
	for (long ix=0; ix<VT->nx; ix++){
	    p[ix]*=Sdiag->p[ix];
	}
    }
    X(zero)(A);
#ifdef USE_COMPLEX
    X(mm)(&A,VT,U,"cc",1);
#else
    X(mm)(&A,VT,U,"tt",1);
#endif
    X(free)(U);
    X(free)(VT);
    XR(free)(Sdiag);
}



/**
   Compute A*B and add to C0.
   C=C+trans(A)*trans(B)*alpha
       
   2009-11-09: There was initially a beta parameter
   It was implemented wrongly for beta!=1 because for every 
   call to dmm, the already accumulated ones are scaled.
   removed beta.
*/
void X(cellmm)(X(cell) **C0, const X(cell) *A, const X(cell) *B, 
	       const char trans[2], const double alpha){
    if(!A || !B) return;
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    if(trans[0]=='n'||trans[0]=='N'){
	nx=A->nx; 
	ax=1; az=A->nx;
	nz=A->ny;
    }else{ 
	nx=A->ny;
	az=1; ax=A->nx;
	nz=A->nx;
    }
    if(trans[1]=='n'||trans[0]=='N'){
	ny=B->ny; 
	bz=1; by=B->nx;
	if(nz!=B->nx) error("mismatch\n");
    }else{
	ny=B->nx;
	by=1; bz=B->nx;
	if(nz!=B->ny) error("mismatch\n");
    }
    if(!*C0){
	*C0=X(cellnew)(nx,ny);
    }
    X(cell) *C=*C0;
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    for(int iz=0; iz<nz; iz++){
		if(A->p[ix*ax+iz*az]&&B->p[iz*bz+iy*by]){
		    X(mm)(&C->p[ix+iy*nx],A->p[ix*ax+iz*az], 
			  B->p[iz*bz+iy*by],trans,alpha);
		}
	    }
	}
    }
}



/**
   Inplace Invert a SPD matrix. It is treated as a block matrix
*/
X(cell)* X(cellinvspd)(X(cell) *A){
    X(mat) *Ab=X(cell2m)(A);
    if(!Ab) return NULL;
    X(invspd_inplace)(Ab);
    X(cell) *B=NULL;
    X(2cell)(&B, Ab, A);
    X(celldropzero)(B,0);
    X(free)(Ab);
    return B;
}

/**
   Inplace Invert a matrix. It is treated as a block matrix.
*/
X(cell)* X(cellinv)(X(cell) *A){
    X(mat) *Ab=X(cell2m)(A);
    X(inv_inplace)(Ab);
    X(cell) *B=NULL;
    X(2cell)(&B, Ab, A);
    X(celldropzero)(B,0);
    X(free)(Ab);
    return B;
}
/**
   invert each component of the dcell. Each cell is treated
   as an individual matrix.
*/
X(cell)* X(cellinvspd_each)(X(cell) *A){
    X(cell) *out=NULL;
    X(cellcp)(&out,A);
    for(int i=0; i<out->nx*out->ny; i++){
	X(invspd_inplace)(out->p[i]);
    }
    return out;
}

/**
   compute the pseudo inverse of block matrix A.  A is n*p cell, wt n*n cell or
   sparse cell.  \f$B=inv(A'*W*A)*A'*W\f$  */
X(cell)* X(cellpinv)(const X(cell) *A,    /**<[in] The matrix to pseudo invert*/
		     const X(cell) *wt,   /**<[in] Use a dense matrix for weighting*/
		     const Y(spcell) *Wsp /**<[in] Use a sparse matrix for weighting*/
		     ){
    if(!A) return NULL;
    X(cell) *wA=NULL;
    if(wt){
	if(Wsp){
	    error("Both wt and Wsp are specified.\n");
	}
	assert(wt->nx==A->nx && wt->ny==A->nx);
	X(cellmm)(&wA, wt, A, "nn",1);
    }else{
	if(Wsp){
	    assert(Wsp->nx==A->nx && Wsp->ny==A->nx);
	    Y(spcellmulmat)(&wA,Wsp,A,1);
	}else{
	    wA=X(cellref)(A);
	}
    }

    X(cell) *ata=NULL;
    X(cellmm)(&ata,wA,A,"tn",1);
    X(cell) *iata=X(cellsvd_pow)(ata, -1, 0, 1e-14);
    X(cellfree)(ata);
    X(cell) *out=NULL;
    X(cellmm)(&out, iata, wA, "nt",1);
    X(cellfree)(wA);
    X(cellfree)(iata);
    return out;
}


/**
   compute the power of a block matrix using svd method. First convert it do
   X(mat), do the power, and convert back to block matrix.
*/
X(cell) *X(cellsvd_pow)(X(cell) *A, double power, int issym, double thres){
    X(mat) *Ac=X(cell2m)(A);
    X(svd_pow)(Ac, power, issym, thres);
    X(cell)*B=NULL;
    X(2cell)(&B, Ac, A);
    X(celldropzero)(B,0);
    X(free)(Ac);
    return B;
}


/**
   Apply tickholov regularization of relative thres to cell array by
   converting it to mat
*/
void X(celltikcr)(X(cell) *A, double thres){
    X(mat) *Ab=X(cell2m)(A);
    X(tikcr)(Ab,thres);
    X(2cell)(&A,Ab,NULL);
    X(free)(Ab);
}

