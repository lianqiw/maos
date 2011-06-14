/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
#include "random.h"
#include "misc.h"
#include "mathmisc.h"
#include "dsp.h"
#include "csp.h"
#include "dmat.h"
#include "cmat.h"
#include "fft.h"
#include "matbin.h"
#include "loc.h"
#if MAT_VERBOSE == 1
#define matinfo(A...) {fprintf(stderr, A);}
#else
#define matinfo(A...)
#endif
#ifndef MAT_TYPE
#define MAT_TYPE
#ifndef USE_COMPLEX
#define X(A) d##A
#define Y(A) A
#define Z(A) d##A##_
#define T double
#define M_T M_DBL
#define M_TT M_DMAT
#define REAL(A) (A)
#define ABS(A) fabs(A)
#define SQRT(A) sqrt(A)
#define RANDU(A) randu(A)
#define RANDN(A) randn(A)
#define PRINT(A) printf("%10.3e",A);
#define CONJ(x) (x)
#define dot_do dotdbl
#else
#define X(A) c##A
#define Y(A) c##A
#define Z(A) z##A##_ //blas/lapack convention
#define T dcomplex
#define M_T M_CMP
#define M_TT M_CMAT
#define REAL(A) creal(A)
#define ABS(A) cabs(A)
#define SQRT(A) csqrt(A)
#define RANDU(A) (randu(A)+I*randu(A))
#define RANDN(A) (randn(A)+I*randn(A))
#define PRINT(A) printf("(%10.3e %10.3eI)",creal(A),cimag(A));
#define CONJ(x) conj(x)
#define dot_do dotcmp
#endif
#endif
#include "dmat.h"
#include "cmat.h"
#define PMAT(A,pp) T (*restrict pp)[(A)->nx]=(void *)(A)->p
#define PCELL(M,P) X(mat)* (*restrict P)[(M)->nx]=(void*)(M)->p

vtbl X(mat_vtbl)={M_TT,
		  (vtbl_write)X(write),
		  (vtbl_writedata)X(writedata),
		  (vtbl_read)X(read),
		  (vtbl_readdata)X(readdata)};
/**
   Work horse function that creats the matrix object. if p is
   NULL, memory is allocated. Allocation for X(mat) objects
   should always be created directly or indirectly through
   X(new_do) so that changing data structure is straight
   forward.
*/
static inline X(mat) *X(new_do)(long nx, long ny, T *p, int ref){
    if(nx==0 || ny==0) return NULL;
    X(mat) *out=calloc(1, sizeof(X(mat)));
    out->vtbl=&X(mat_vtbl);
    out->nx=nx;
    out->ny=ny;
    if(ref){//the data does not belong to us.
	if(!p && nx*ny!=0){
	    error("When ref is 1, p must not be NULL.\n");
	}
	out->p=p;
    }else{
	if(!p){
	    p=calloc((nx*ny), sizeof(T));
	}
	out->p=p;
	out->nref=calloc(1, sizeof(T));
	out->nref[0]=1;
    }
    return out;
}
/**
   Creat a X(mat) object to reference an already existing
   vector.  free the X(mat) object won't free the existing
   vector.
*/
X(mat) *X(new_ref)(long nx, long ny, T *p){
    return X(new_do)(nx,ny,p,1);
}

/**
   Creat a X(mat) object with already allocated memory
   chunk. the memory is freed when the memory is freed.
*/
X(mat) *X(new_data)(long nx, long ny, T *p){
    return X(new_do)(nx,ny,p,0);
}

/**
   Create a new T matrix object. initialized all to zero.
*/
X(mat) *X(new)(long nx, long ny){
    return X(new_do)(nx,ny,NULL,0);
}

/**
   Free the X(mat), but keep the data.
*/
void X(free_keepdata)(X(mat) *A){
    X(free_do)(A,1);
}
/**
   free a X(mat) object. if keepdata!=0, will not free A->p.
*/
void X(free_do)(X(mat) *A, int keepdata){
    if(!A) return;
    if(A->nref){
	if(A->nref[0]==1){
#ifdef USE_COMPLEX
	    cfree_plan(A);
#endif
	    if(A->header){
		double count=search_header_num(A->header, "count");
		if(count>0){
		    info("count=%g, scaling the data\n", count);
		    X(scale)(A, 1./count);
		}
	    }
	    if(!keepdata && A->p){
		if(A->mmap){//data is mmap'ed.
		    mmap_unref(A->mmap);
		}else{
		    free(A->p);
		    free(A->header);
		}
	    }
	    free(A->nref);
	}else if(A->nref[0]>1){
	    A->nref[0]--;
	}else{
	    error("The ref is less than 1. unlikely!!!:%ld\n",A->nref[0]);
	}
    }else{//data does not belong to us.
	free(A->header);
    }
    free(A);
}

/**
   Resize a matrix by adding or removing columns or rows. Data is kept
   whever possible.
*/
void X(resize)(X(mat) *A, long nx, long ny){
    if(A->nref[0]>1){
	error("Resizing a referenced vector\n");
    }
    if(nx==0 || ny==0){
	warning("resizing to (%ld, %ld). Free the vector.\n", nx,ny);
	X(free)(A);
    }
    if(A->nx==nx || A->ny==1){
	A->p=realloc(A->p, sizeof(T)*nx*ny);
	if(nx*ny>A->nx*A->ny){
	    memset(A->p+A->nx*A->ny, 0, (nx*ny-A->nx*A->ny)*sizeof(T));
	}
    }else{
	warning("column vector length is not preserved!!!");
	T *p=calloc(nx*ny,sizeof(T));
	long minx=A->nx<nx?A->nx:nx;
	long miny=A->ny<ny?A->ny:ny;
	for(long iy=0; iy<miny; iy++){
	    memcpy(p+iy*nx, A->p+iy*A->nx, sizeof(T)*minx);
	}
	free(A->p);
	A->p=p;
    }
    A->nx=nx;
    A->ny=ny;
}

/**
   creat a X(mat) reference an existing X(mat). Use the reference carefully.
*/
X(mat) *X(ref)(X(mat) *in){
    if(!in) return NULL;
    X(mat) *out=calloc(1, sizeof(X(mat)));
    if(!out){
	error("Allocation failed\n");
    }
    memcpy(out,in,sizeof(X(mat)));
    out->nref[0]++;
    return out;
}

/**
   create an new X(mat) reference another with different shape.
*/
X(mat) *X(ref_reshape)(X(mat) *in, long nx, long ny){
    X(mat) *out=X(ref)(in);
    if(in->nx*in->ny!=nx*ny){
	error("Must not change number of elements\n");
    }
    out->nx=nx;
    out->ny=ny;
    return out;
}

/**
   creat a new X(mat) referencing columns in existing
   X(mat). reference counted. not used
*/
X(mat)* X(refcols)(X(mat) *in, long icol, long ncol){
    X(mat) *out=calloc(1, sizeof(X(mat)));
    out->nx=in->nx;
    out->ny=ncol;
    out->p=in->p+icol*out->nx;
    return out;
}

/**
   Create a new sub matrix of nx*ny starting from(sx,sy)
*/
X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny){
    if(nx<=0){
	nx=in->nx-sx;
    }
    if(ny<=0){
	ny=in->ny-sy;
    }
    X(mat)*out=X(new)(nx, ny);
    if(sx+nx>in->nx || sy+ny>in->ny){
	error("Invalid parameter range\n");
    }
    PMAT(in, pin);
    PMAT(out, pout);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    pout[iy][ix]=pin[iy+sy][ix+sx];
	}
    }
    return out;
}
/**
 * Check for Nan in elements
 */
int X(isnan)(const X(mat)*A){
    for(long i=0; i<A->nx*A->ny; i++){
	if(isnan(A->p[i])){
	    return 1;
	}
    }
    return 0;
}
/**
   concatenate two matrixes into 1 along dimension "dim"
*/
X(mat)* X(cat)(const X(mat) *in1, const X(mat) *in2, int dim){
    if(!in2){
	if(in1){
	    return X(dup)(in1);
	}else{
	    return NULL;
	}
    }else if(!in1){
	return X(dup)(in2);
    }

    X(mat) *out=NULL;

    if(dim==1){
	//along x.
	if(in1->ny!=in2->ny){
	    error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld)\n", 
		  in1->nx, in1->ny, in2->nx, in2->ny);
	}
	out=X(new)(in1->nx+in2->nx, in1->ny);
	PMAT(out,pout);
	PMAT(in1,pin1);
	PMAT(in2,pin2);
	for(long iy=0; iy<in1->ny; iy++){
	    memcpy(pout[iy],pin1[iy], in1->nx*sizeof(T));
	    memcpy(pout[iy]+in1->nx, pin2[iy], in2->nx*sizeof(T));
	}
    }else if(dim==2){
	//along y.
	if(in1->nx!=in2->nx){
	    error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld)\n", 
		  in1->nx, in1->ny, in2->nx, in2->ny);
	}
	out=X(new)(in1->nx, in1->ny+in2->ny);
	memcpy(out->p, in1->p, in1->nx*in1->ny*sizeof(T));
	memcpy(out->p+in1->nx*in1->ny,in2->p,in2->nx*in2->ny*sizeof(T));
    }else{
	error("Invalid dim\n");
    }
    return out;
}

/**
   free a X(mat) array.
*/
void X(arrfree)(X(mat) **As, int n){
    int i;
    if(!As) return;
    for(i=0; i<n; i++){
	X(free)(As[i]);
    }
    free(As);
}

/**
   find the maximum value of a X(mat) object
*/
double X(max)(const X(mat) *A){
    double max,min;
#ifdef USE_COMPLEX
    maxmincmp(A->p,A->nx*A->ny,&max,&min,NULL);
#else
    maxmindbl(A->p,A->nx*A->ny,&max,&min);
#endif
    return max;
}

/**
   find the minimum value of a X(mat) object
*/
double X(min)(const X(mat) *A){
    double max,min;
#ifdef USE_COMPLEX
    maxmincmp(A->p,A->nx*A->ny,&max,&min,NULL);
#else
    maxmindbl(A->p,A->nx*A->ny,&max,&min);
#endif
    return min;
}

/**
   duplicate a X(mat) array
*/
X(mat) *X(dup)(const X(mat) *in){
    X(mat) *out=NULL;
    X(cp)(&out, in);
    return out;
}

/**
   copy the values from one X(mat) to another.
*/
void X(cp)(X(mat) **out0, const X(mat) *in){
    if(in && in->nx!=0 && in->ny!=0){
	if(!*out0) 
	    *out0=X(new)(in->nx, in->ny);
	else{
	    assert(in->nx==(*out0)->nx && in->ny == (*out0)->ny);
	}
	X(mat) *out=*out0;
	memcpy(out->p, in->p, in->nx*in->ny*sizeof(T));
    }else{
	X(zero)(*out0);
    }
}

/**
   transpose a X(mat) object
*/
X(mat) *X(trans)(const X(mat) *A){
    if(!A) return NULL;
    X(mat) *B=X(new)(A->ny,A->nx);
    if(A->nx==1 || A->ny==1){
	memcpy(B->p, A->p, A->nx*A->ny*sizeof(T));
    }else{
	PMAT(B,Bp);
	PMAT(A,Ap);
	for(int ix=0; ix<A->nx; ix++){
	    for(int iy=0; iy<A->ny; iy++){
		Bp[ix][iy]=Ap[iy][ix];
	    }
	}
    }
    return B;
}

/**
   set values of each element in a X(mat) to val.
*/
void X(set)(X(mat) *A, const T val){
    if(A){
	for(long i=0; i<A->nx*A->ny; i++){
	    A->p[i]=val;
	}
    }
}

/**
   compute the norm2 of A
*/
double X(norm2)(const X(mat)*A){
    double out=0;
    for(int i=0; i<A->nx; i++){
	out+=(double)(A->p[i]*CONJ(A->p[i]));
    }
    return out;
}

/**
   Fill A with random uniform numbers between [0, 1]*max
*/
void X(randu)(X(mat) *A, const T max, rand_t *rstat){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]=RANDU(rstat)*max;
    }
}
/**
   Fill A with random normal distribution numbers with
   standard deviation of sigma.
*/
void X(randn)(X(mat) *A, const T sigma, rand_t *rstat){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]=RANDN(rstat)*sigma;
    }
}

/**
   display a X(mat) matrix.
*/
void X(show)(const X(mat) *A, const char *format, ...){
    if(!A) return;
    format2fn;
    info2("Displaying content of %s:\n",fn);
    PMAT(A,p);
    long colmax=10;
    long icol,i,j;
    for(icol=0; icol<ceil((double)A->ny/(double)colmax); icol++){
	int ncol=(icol+1)*colmax;
	if(ncol>A->ny) ncol=A->ny;
	printf("Cols %ld to %d\n", icol, ncol-1);
	for(j=0; j<A->nx; j++){
	    for(i=icol*colmax; i<ncol; i++){
		PRINT(p[i][j]);
	    }
	    printf("\n");
	}
    }
}

/**
   scale each element of A by w
*/
void X(scale)(X(mat) *A, T w){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]*=w;
    }
}

/**
   create sum of all the elements in A.
*/
T X(sum)(const X(mat) *A){
    T v=0;
    if(A){
	T *restrict p=A->p;
	/**
	   Loops like this will only be vectorized with -ffast-math because
	   different orders of accumulation give different results for floating
	   point numbers.
	*/
	for(int i=0; i<A->nx*A->ny; i++){
	    v+=p[i];
	}
    }
    return v;
}

/**
   compute B=bc*B+ac*A
   behavior changed on 2009-11-02. if A is NULL, don't do anything.
*/
void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac){
    if(A){
	if(!*B0){
	    *B0=X(new)(A->nx, A->ny); 
	    bc=0;//no bother to accumulate.
	}
	X(mat) *B=*B0;
	assert(A->nx==B->nx && A->ny == B->ny);
	if(fabs(bc)>EPS){
	    for(int i=0; i<A->nx*A->ny; i++){
		B->p[i]=B->p[i]*bc+A->p[i]*ac;
	    }
	}else{//just assign
	    for(int i=0; i<A->nx*A->ny; i++){
		B->p[i]=A->p[i]*ac;
	    }
	}
    }
}

/**
   compute the inner product of A and B. (inner product)
*/

T X(inn)(const X(mat)*A, const X(mat) *B){
    if(!A || !B) return 0;
    assert(A->nx==B->nx && A->ny==1 && B->ny==1);
    T out=0;
    for(int i=0; i<A->nx; i++){
	out+=A->p[i]*B->p[i];
    }
    if(isnan(out)){
	error("NaN found\n");
    }
    return out;
}

/**
   compute weighted dot product a'*(w*b)
*/
T X(wdot)(const T *a, const X(mat) *w, const T *b){
    PMAT(w,pw);
    T res=0;
    for(int j=0; j<w->ny; j++){
	for(int i=0; i<w->nx; i++){
	    res+=pw[j][i]*a[i]*b[j];
	}
    }
    return res;
}

/**
   special version of dwdot for just 2 element vectors.
*/
T X(wdot2)(const T *a, const X(mat) *w, const T *b){
    assert(w->nx==2 && w->ny==2);
    PMAT(w,W);
    T res;
    res=a[0]*(W[0][0]*b[0]+W[1][0]*b[1])
	+a[1]*(W[0][1]*b[0]+W[1][1]*b[1]);
    return res;
}

/**
   special version of dwdot for just 3 element  vectors.
*/
T X(wdot3)(const T *a, const X(mat) *w, const T *b){
    assert(w->nx==3 && w->ny==3);
    PMAT(w,W);
    T res;
    res=a[0]*(W[0][0]*b[0]+W[1][0]*b[1]+W[2][0]*b[2])
	+a[1]*(W[0][1]*b[0]+W[1][1]*b[1]+W[2][1]*b[2])
	+a[2]*(W[0][2]*b[0]+W[1][2]*b[1]+W[2][2]*b[2]);
    return res;
}

/**
   Compute component wise multiply B=B.*A
*/
void X(cwm)(X(mat) *B, const X(mat) *A){
    assert(A->nx==B->nx && A->ny==B->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	B->p[i]*=A->p[i];
    }
}

/**
   multiply a X(mat) matrix with a vector and accumulate to y:
   y+=A*x*alpha
*/
void X(mulvec)(T *restrict y, const X(mat) * restrict A,
	       const T *restrict x, const T alpha){
    assert(y && x && A);
    PMAT(A,Ap);
    if(ABS(alpha-1)>1.e-15){
	for(int iy=0; iy<A->ny; iy++){
	    T tmp=x[iy]*alpha;
	    for(int ix=0; ix<A->nx; ix++){
		y[ix]+=Ap[iy][ix]*tmp;
	    }
	}
    }else{
	for(int iy=0; iy<A->ny; iy++){
	    for(int ix=0; ix<A->nx; ix++){
		y[ix]+=Ap[iy][ix]*x[iy];
	    }
	}
    }
}

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
    if(k!=k2) error("dmm: Matrix doesn't match\n");
    if(!*C0){
	*C0=X(new)(m,n); 
    }else if(m!=(*C0)->nx || n!=(*C0)->ny){
	error("dmm: Matrix doesn't match\n");
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
    // B is identity matrix
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
    for(int i=0;i<N;i++)
	B[i+i*N]=1;
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
    }else{//Zero out the upper diagonal. For some reason they are not zero.
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
   compute (A'*W*A); where diag(W)=wt
*/
X(mat) *X(mcc)(const X(mat) *A, const X(mat) *wt){
    assert(A->nx==wt->nx && wt->ny==1);
    int nmod=A->ny;
    int nsa2=A->nx;
    X(mat) *ata=X(new)(nmod, nmod);;
    PMAT(ata,ATA);
    PMAT(A,Ap);
    for(int imod=0; imod<nmod; imod++){
	for(int jmod=imod; jmod<nmod; jmod++){
	    T tmp=0;
	    for(long ik=0; ik<nsa2; ik++){
		tmp+=Ap[imod][ik]*Ap[jmod][ik]*wt->p[ik];
	    }
	    ATA[imod][jmod]=tmp;
	    if(imod!=jmod)
		ATA[jmod][imod]=ATA[imod][jmod];
	}
    }
    return ata;
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
   compute (A*W*A'); where diag(W)=wt
*/
X(mat) *X(tmcc)(const X(mat) *A, const X(mat) *wt){
    assert(A->ny==wt->nx && wt->ny==1);
    int nmod=A->nx;
    int nsa2=A->ny;
    X(mat) *ata=X(new)(nmod, nmod);;
    PMAT(ata,ATA);
    PMAT(A,Ap);
    for(int imod=0; imod<nmod; imod++){
	for(int jmod=imod; jmod<nmod; jmod++){
	    T tmp=0;
	    for(int k=0; k<nsa2; k++){
		tmp+=Ap[k][imod]*Ap[k][jmod]*wt->p[k];
	    }
	    ATA[imod][jmod]=tmp;
	    if(imod!=jmod)
		ATA[jmod][imod]=ATA[imod][jmod];
	}
    }
    return ata;
}

/**
   compute the pseudo inverse of matrix A with weigthing of full matrix W or
   sparse matrix weighting Wsp.  For full matrix, wt can be either W or diag (W)
   for diagonal weighting.  B=inv(A'*W*A)*A'*W; */
X(mat) *X(pinv)(const X(mat) *A, const X(mat) *wt, const X(sp) *Wsp){
    if(!A) return NULL;
    X(mat) *AtW=NULL;
    //Compute AtW=A'*W
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
    //Compute cc=A'*W*A
    X(mat) *cc=NULL;
    X(mm) (&cc, AtW, A, "nn", 1);
    //Compute inv of cc
    //X(invspd_inplace)(cc);
    if(X(isnan(cc))){
	X(write)(cc,"cc_isnan");
	X(write)(A,"A_isnan");
	X(write)(wt,"wt_isnan");
	Y(spwrite)(Wsp, "Wsp_isnan");
    }
    X(svd_pow)(cc,-1,0,1e-14);//invert the matrix using SVD. safe with small eigen values.
    X(mat) *out=NULL;
    //Compute (A'*W*A)*A'*W
    X(mm) (&out, cc, AtW, "nn", 1);
    X(free)(AtW);
    X(free)(cc);
    return out;
}

/**
   compute the relative difference betwee two vectors.
   ||A-B||/||A|| using norm2. for debugging purpose.
*/
T X(diff)(const X(mat) *A, const X(mat) *B){
    X(mat) *C=NULL;
    X(cp)(&C,A);
    X(add)(&C,1,B,-1);
    T d=sqrt(X(norm2)(C)*2/(X(norm2)(C)+X(norm2)(B)));
    return isnan(d)?0:d;
}

/**
   a new gray pixel map generation based on bilinear influence
   functions used in mkw.  creates slightly larger map.  add
   an filled circle.  cx,cy,r are in unit of pixels
*/
void X(circle)(X(mat) *A, double cx, double cy, double r, T val){
    int nres=100;
    const double res=1./(double)(nres);
    const double res1=1./(double)(nres);
    const double res2=res1*res1*4.;
    double resm=(double)(nres-1)/2.;
    double r2=r*r;
    double r2l=(r-1.5)*(r-1.5);
    double r2u=(r+2.5)*(r+2.5);
    PMAT(A,As);
    for(int iy=0; iy<A->ny; iy++){
	double r2y=(iy-cy)*(iy-cy);
	for(int ix=0; ix<A->nx; ix++){
	    double r2r=(ix-cx)*(ix-cx)+r2y;
	    if(r2r<r2l) 
	    	As[iy][ix]+=val;
	    else if(r2r<r2u){
		double tot=0.;
		for(int jy=0; jy<nres; jy++){
		    double iiy=iy+(jy-resm)*2*res;
		    double rr2y=(iiy-cy)*(iiy-cy);
		    double wty=1.-fabs(iy-iiy);
		    for(int jx=0; jx<nres; jx++){
			double iix=ix+(jx-resm)*2*res;
			double rr2r=(iix-cx)*(iix-cx)+rr2y;
			double wtx=1.-fabs(ix-iix);
			if(rr2r<r2)
			    tot+=res2*wty*wtx;
		    }
		}
		As[iy][ix]+=tot*val;
	    }
	}
    }
}

/**
   similar to X(circle). but don't actually compute the
   weights. just test the corners;
*/
void X(circle_symbolic)(X(mat) *A, double cx, double cy, double r){
    double r2=r*r;
    double r2u=(r+2.5)*(r+2.5);
    PMAT(A,As);
    for(int iy=0; iy<A->ny; iy++){
	double r2y=(iy-cy)*(iy-cy);
	for(int ix=0; ix<A->nx; ix++){
	    double r2r=(ix-cx)*(ix-cx)+r2y;
	    if(r2r<r2) 
	    	As[iy][ix]=1;
	    else if(r2r<r2u){
		for(int jy=-1; jy<2; jy++){
		    double iiy=iy+jy;
		    double rr2y=(iiy-cy)*(iiy-cy);
		    for(int jx=-1; jx<2; jx++){
			double iix=ix+jx;
			double rr2r=pow(iix-cx,2)+rr2y;
			if(rr2r<r2){
			    As[iy][ix]=1;
			    goto next;
			}
		    }
		}
	    }
	next:
	    continue;
	}
    }
}

/**
   shift frequency components by n/2
*/
void X(fftshift)(X(mat) *A){
    size_t i;
    const size_t nx=A->nx;
    const size_t ny=A->ny;
    assert ((nx&1)==0);
    const size_t nx2=nx/2;
    const size_t ny2=ny/2;
    const size_t nx2d=nx2*sizeof(T);
    T *tmp=(T*)malloc(nx2d);
    T *data=A->p;
    if(ny==1){
	memcpy(tmp,data,nx2d);
	memcpy(data,data+nx2,nx2d);
	memcpy(data+nx2,tmp,nx2d);
    }else{
	assert((ny&1)==0);
	for(i=0; i<ny2; i++){
	    memcpy(tmp,data+i*nx,nx2d);
	    memcpy(data+i*nx,data+(i+ny2)*nx+nx2, nx2d);
	    memcpy(data+(i+ny2)*nx+nx2,tmp, nx2d);
	    memcpy(tmp,data+i*nx+nx2,nx2d);
	    memcpy(data+i*nx+nx2,data+(i+ny2)*nx, nx2d); 
	    memcpy(data+(i+ny2)*nx,tmp, nx2d);
	}
    }
    
    free(tmp);
}

/**
   reorder B
   and embed/crop into center of A 
   \verbatim
   4 * * 3
   * * * *
   * * * *
   2 * * 1
   \endverbatim
   to
   \verbatim
   1 2 
   3 4
   \endverbatim
*/
void X(cpcorner2center)(X(mat) *A, const X(mat)*B){
    const size_t nx=A->nx;
    const size_t ny=A->ny;
    T *Ap=A->p;
    memset(Ap, 0, sizeof(T)*nx*ny);
    const size_t ninx=B->nx;
    const size_t niny=B->ny;
    const T * Bp=B->p;
    assert((nx&1)==0 && (ny&1)==0 && (ninx&1)==0 && (niny&1)==0);

    const int ny2=(ny<niny)?ny/2:niny/2;
    const int nx2=(nx<ninx)?nx/2:ninx/2;
    const int xskip=nx/2-nx2;
    const int yskip=ny/2-ny2;
    T* Ap0=Ap+yskip*nx+xskip;
    const int nx2d=nx2*sizeof(T);
    for(int i=0; i<ny2; i++){
	memcpy(Ap0+i*nx, Bp+(niny-ny2+i)*ninx+(ninx-nx2),nx2d); 
	memcpy(Ap0+i*nx+nx2, Bp+(niny-ny2+i)*ninx, nx2d); 
	memcpy(Ap0+(i+ny2)*nx, Bp+i*ninx+(ninx-nx2), nx2d); 
	memcpy(Ap0+(i+ny2)*nx+nx2, Bp+i*ninx, nx2d); 
    }
}

/**
   cyclic shift A by nx and ny to B.
   \verbatim
   4   3     1   2 
      
   2   1 to  3   4
   \endverbatim
*/
void X(shift)(X(mat) **B0, const X(mat) *A, int sx, int sy){
    if(!*B0){
	*B0=X(new)(A->nx, A->ny);
    }
    X(mat) *B=*B0;

    const int nx=A->nx; 
    const int ny=A->ny;
    sx=sx%nx; if(sx<0) sx+=nx;
    sy=sy%ny; if(sy<0) sy+=ny;
    if(sx!=0 || sy!=0){
	int dy=ny-sy;
	int dx=nx-sx;
	for(int iy=0; iy<sy; iy++){
	    memcpy(B->p+iy*nx, A->p+(dy+iy)*nx+dx, sx*sizeof(T));//3
	    memcpy(B->p+iy*nx+sx, A->p+(dy+iy)*nx, dx*sizeof(T));//4
	}
	for(int iy=sy; iy<ny; iy++){
	    memcpy(B->p+iy*nx, A->p+(iy-sy)*nx+dx, sx*sizeof(T));//1
	    memcpy(B->p+iy*nx+sx, A->p+(iy-sy)*nx, dx*sizeof(T));//2
	}
    }else{
	memcpy(B->p, A->p, sizeof(T)*nx*ny);
    }
}

/**
   rotate the vectors CCW.
   same as rotate coordinate theta CW.
   A(:,1) is x, A(:,2) is y.
*/
void X(rotvec)(X(mat) *A, const double theta){
    if(A->ny!=2) error("Wrong dimension\n");
    const double ctheta=cos(theta);
    const double stheta=sin(theta);
    PMAT(A,Ap);
    for(int i=0; i<A->nx; i++){
	T tmp=Ap[0][i]*ctheta-Ap[1][i]*stheta;
	Ap[1][i]=Ap[0][i]*stheta+Ap[1][i]*ctheta;
	Ap[0][i]=tmp;
    }
}

/**
   rotate a 2x2 covariance matrix A by theta CCW
   (coordinate rotate -theta CCW) or from ra to xy
   coordinate.  R*A*R';
*/
void X(rotvecnn)(X(mat) **B0, const X(mat) *A, double theta){
    assert(A->nx==2 && A->ny==2);
    if(!*B0) 
	*B0=X(new)(2,2);
    X(mat) *B=*B0;
    assert(B->nx==2 && B->ny==2);
    const T ctheta=cos(theta);
    const T stheta=sin(theta);
    T tmp[2][2];
    PMAT(A,Ap);
    PMAT(B,Bp);
    //first apply left R
    tmp[0][0]=ctheta*Ap[0][0]-stheta*Ap[0][1];
    tmp[1][0]=ctheta*Ap[1][0]-stheta*Ap[1][1];
    tmp[0][1]=stheta*Ap[0][0]+ctheta*Ap[0][1];
    tmp[1][1]=stheta*Ap[1][0]+ctheta*Ap[1][1];
    //then apply right R'
    
    Bp[0][0]=ctheta*tmp[0][0]-stheta*tmp[1][0];
    Bp[1][0]=stheta*tmp[0][0]+ctheta*tmp[1][0];
    Bp[0][1]=ctheta*tmp[0][1]-stheta*tmp[1][1];
    Bp[1][1]=stheta*tmp[0][1]+ctheta*tmp[1][1];

}

/**
   T matrix vector multiply optimized for just
   three values.  y=A*x;
*/
void X(mulvec3)(T *y, const X(mat) *A, const T *x){
    assert(A->nx==3 && A->ny==3);
    PMAT(A,Ap);
    //calculate y=A*x for 3.
    y[0]=Ap[0][0]*x[0]+Ap[1][0]*x[1]+Ap[2][0]*x[2];
    y[1]=Ap[0][1]*x[0]+Ap[1][1]*x[1]+Ap[2][1]*x[2];
    y[2]=Ap[0][2]*x[0]+Ap[1][2]*x[1]+Ap[2][2]*x[2];
}


/**
   Compute thresholded center of gravity. The threshold
   is absolute value. bkgrnd is removed from i0 when
   computing cog.  offset is the offset of the reference
   point (cog=0) from the physical center.
   all length are given in terms of pixel.
*/
void X(cog)(double *grad,const X(mat) *i0,double offsetx,
	    double offsety, double thres, double bkgrnd){
    double sum=0,sumx=0,sumy=0;
    double iI;
    PMAT(i0,pi0);
    for(int iy=0; iy<i0->ny; iy++){
	for(int ix=0; ix<i0->nx; ix++){
	    iI=REAL(pi0[iy][ix])-bkgrnd;
	    if(iI>thres){
		sum+=iI;
		sumx+=iI*ix;
		sumy+=iI*iy;
	    }
	}
    }
    if(fabs(sum)>0){
	grad[0]=sumx/sum-((double)(i0->ny-1)*0.5+offsetx);
	grad[1]=sumy/sum-((double)(i0->nx-1)*0.5+offsety);
    }else{
	grad[0]=0;
	grad[1]=0;
    }
}

/**
   Shift the image in A to center on physical
   center+[offsetx,offsety] using cog and fft.
*/
void X(shift2center)(X(mat) *A, double offsetx, double offsety){
    double grad[2];
    double Amax=X(max)(A);
    X(cog)(grad,A,offsetx,offsety,Amax*0.1,Amax*0.2);
    if(fabs(grad[0])>0.1 || fabs(grad[1])>0.1){
	//info("Before shift, residual grad is %g %g\n",grad[0],grad[1]);
	cmat *B=cnew(A->nx,A->ny);
	cfft2plan(B,-1);
	cfft2plan(B,1);
#ifdef USE_COMPLEX
	ccp(&B,A);
#else
	ccpd(&B,A);
#endif
	double scale=1./(A->nx*A->ny);
	cfftshift(B);
	cfft2(B,-1);
	ctilt(B,-grad[0],-grad[1],0);
	cfft2(B,1);
	cfftshift(B);
	cscale(B,scale);
#ifdef USE_COMPLEX
	ccp(&A,B);
#else
	creal2d(&A,0,B,1);
#endif
	X(cog)(grad,A,offsetx,offsety,Amax*0.1,Amax*0.2);
	//info("After shift, residual grad is %g %g\n",grad[0],grad[1]);
	cfree(B);
    }
}

/**
   Limit numbers in A to within [min, max]. used for DM clipping.
*/
int X(clip)(X(mat) *A, double min, double max){
    if(!A) return 0;
    if(isinf(min)==-1 && isinf(max)==1) return 0;
    if(max<=min){
	error("upper light should be larger than lower limit\n");
    }
    T *restrict Ap=A->p;
    int nclip=0;
    for(long i=0; i<A->nx *A->ny; i++){
	if(REAL(Ap[i])>max) {
	    Ap[i]=max;
	    nclip++;
	}else if(REAL(Ap[i])<min) {
	    Ap[i]=min;
	    nclip++;
	}
    }
    return nclip;
}

/**
   OrthNormalize column vector in Mod, with weighting from vector amp.
   <Mod|wt|Mod> is equal to sum(wt).
   2010-07-21: Bug found: The result is not orthonormal. cause: nonvalid is not initialized to 0.
*/
void X(gramschmidt)(X(mat) *Mod, double *amp){
    const int nmod=Mod->ny;
    const long nx=Mod->nx;
    T wtsum=(T)nx;
    if(amp){
	wtsum=dblsum(amp, nx);
    }
    int nonvalid[nmod];
    memset(nonvalid, 0, sizeof(int)*nmod);
    PMAT(Mod,pMod);
    for(int imod=0; imod<nmod; imod++){
	if(imod>0){//orthogonalize
	    T cross;
	    //compute dot product.
	    for(int jmod=0; jmod<imod; jmod++){
		if(nonvalid[jmod]) continue;
		cross=-dot_do(pMod[imod],pMod[jmod],amp,nx)/wtsum;
		for(long ix=0; ix<nx; ix++){
		    pMod[imod][ix]+=cross*pMod[jmod][ix];
		}
	    }
	}
	
	T norm=SQRT(dot_do(pMod[imod],pMod[imod],amp,nx)/wtsum);
	if(ABS(norm)>1.e-15){
	    norm=1./norm;
	    for(long ix=0; ix<nx; ix++){
		pMod[imod][ix]*=norm;
	    }
	}else{
	    nonvalid[imod]=1;
	    warning("Column %d is not independent on other columns\n",imod);
	}
    }
}

/**
   A=A*B, where diag(B)=s
*/
void X(muldiag)(X(mat) *A, X(mat) *s){
    assert(A->ny==s->nx && s->ny==1);
    PMAT(A,pA);
    const T *ps=s->p;
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    pA[iy][ix]*=ps[iy];
	}
    }
}

/**
   Raise all elements to power power
*/
void X(cwpow)(X(mat)*A, double power){
    if(!A) return;
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=Y(pow)(A->p[i],power);
    }
}

/**
   Compute SVD of a general matrix A. 
   A=U*diag(S)*V';
   diag(S) is returned.
*/
void X(svd)(X(mat) **U, dmat **Sdiag, X(mat) **VT, const X(mat) *A){
    char jobuv='S';
    int M=(int)A->nx;
    int N=(int)A->ny;
    /*if((Sdiag&&*Sdiag)||(U&&*U)||(VT&&*VT)){
	warning("Sdiag,U,VT should all be NULL. discard their value\n");
	}*/
    X(mat) *tmp=X(dup)(A);
    int nsvd=M<N?M:N;
    dmat *s=dnew(nsvd,1);
    X(mat) *u=X(new)(M,nsvd);
    X(mat) *vt=X(new)(nsvd,N);
    int lwork=-1;
    T work0[1];
    int info=0;
#ifdef USE_COMPLEX
    warning("Not tested\n");
    double *rwork=malloc(nsvd*5*sizeof(double));
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
    if(Sdiag) *Sdiag=s; else dfree(s);
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
void X(evd)(X(mat) **U, dmat **Sdiag,const X(mat) *A){
    assert(A->nx==A->ny && A->nx>0);
    *Sdiag=dnew(A->nx,1);
    char jobz=U?'V':'N';
    char uplo='U';
    int lda=A->nx;
    T worksize;
    int lwork=-1;
    int info;
    X(mat) *atmp=X(dup)(A);
#ifdef USE_COMPLEX
    double *rwork=malloc((3*A->nx-2)*sizeof(double));
    Z(heev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, &worksize, &lwork,rwork, &info);
    lwork=(int)worksize;
    T *work=malloc(sizeof(double)*lwork);
    Z(heev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, work, &lwork,rwork, &info);
    free(rwork);
#else
    Z(syev)(&jobz, &uplo, &lda, atmp->p, &lda, (*Sdiag)->p, &worksize, &lwork, &info);
    lwork=(int)worksize;
    T *work=malloc(sizeof(double)*lwork);
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
*/
void X(svd_pow)(X(mat) *A, double power, int issym, double thres){
    int use_evd=0;
    if(A->nx!=A->ny){
	warning("dpow is only good for square arrays.\n");
    }else if(issym){
	use_evd=1;
    }
    dmat *Sdiag=NULL;
    X(mat) *U=NULL;
    X(mat) *VT=NULL;
    double maxeig;
    if(use_evd){
	X(evd)(&U, &Sdiag, A);
	//eigen values below the threshold will not be used. the last is the biggest.
	maxeig=fabs(Sdiag->p[Sdiag->nx-1]);
	VT=X(trans)(U);
    }else{
	X(svd)(&U, &Sdiag, &VT, A);
	//eigen values below the threshold will not be used. the first is the biggest.
	maxeig=fabs(Sdiag->p[0]);
    }
    thres*=maxeig;
    long skipped=0;
    double mineig=INFINITY;
    for(long i=0; i<Sdiag->nx; i++){
	if(fabs(Sdiag->p[i])>thres){//only do with 
	    Sdiag->p[i]=pow(Sdiag->p[i],power);
	}else{
	    if(fabs(Sdiag->p[i])<mineig) 
		mineig=fabs(Sdiag->p[i]);
	    Sdiag->p[i]=0;
	    skipped++;
	}
    }
    for(long iy=0; iy <VT->ny; iy++){
	T *p=VT->p+iy*VT->nx;
	for (long ix=0; ix<VT->nx; ix++){
	    p[ix]*=Sdiag->p[ix];
	}
    }
    X(zero)(A);
    X(mm)(&A,U,VT,"nn",1);
    
    X(free)(U);
    X(free)(VT);
    dfree(Sdiag);
}

/**
   add val to diagonal values of A.
*/
void X(addI)(X(mat) *A, T val){
    if(!A) return;
    if(A->nx!=A->ny)
	warning("daddI: A is not square\n");
    long M=A->nx<A->ny?A->nx:A->ny;
    for(long i=0; i<M; i++){
	A->p[i+i*A->nx]+=val;
    } 
}

/**
   Apply tikhonov regularization to A
*/
void X(tikcr)(X(mat) *A, T thres){
    dmat *S=NULL;
    X(svd)(NULL,&S,NULL,A);
    T val=S->p[0]*thres;
    dfree(S);
    X(addI)(A,val);
}

/**
   y=y+alpha*x*A;
   implemented by transposing x,y index in sptmulmat implementation
   TESTED OK.
*/
void X(mulsp)(X(mat) **yout, const X(mat) *x,const X(sp) *A, const T alpha){
    if(A&&x){
	long icol, ix;
	if(!*yout){
	    *yout=X(new)(x->nx, A->n);
	}
	X(mat) *y=*yout;
	assert(x->nx==y->nx);
	if(x->nx==1){
	    Y(sptmulvec)(y->p, A, x->p, alpha);
	}else{
	    int jcol;
	    PMAT(y,Y); PMAT(x,X);
	    if(ABS(alpha-1.)<1.e-100){
		for(icol=0; icol<A->n; icol++){
		    for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(jcol=0; jcol<y->nx; jcol++){
			    Y[icol][jcol]+=A->x[ix]*X[A->i[ix]][jcol];
			}
		    }
		}
	    }else{
		for(icol=0; icol<A->n; icol++){
		    for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(jcol=0; jcol<y->nx; jcol++){
			    Y[icol][jcol]+=alpha*A->x[ix]*X[A->i[ix]][jcol];
			}
		    }
		}
	    }
	}
    }
}

/**
   Create log spaced vector.
*/
X(mat)* X(logspace)(double emin, double emax, long n){
    X(mat)* out=X(new)(n,1);
    double esep=(emax-emin)/(n-1);
    for(long i=0; i<n; i++){
	double ex=emin+esep*i;
	out->p[i]=pow(10, ex);
    }
    return out;
}

/**
   Create linearly spaced vector.
*/
X(mat)* X(linspace)(double min, double dx, long n){
    X(mat)* out=X(new)(n,1);
    for(long i=0; i<n; i++){
	out->p[i]=min+dx*i;
    }
    return out;
}

/**
   Interpolate using linear interp. xin is the coordinate of yin. xnew is the
   coordinate of the output.
*/
X(mat)* X(interp1)(dmat *xin, dmat *yin, dmat *xnew){
    long nmax=xin->nx;
    long nmax1=nmax-1;
    double xminl=(xin->p[0]);
    double xmaxl=(xin->p[nmax-1]);
    double xsep=(xmaxl-xminl)/(double)(nmax1);
    double xsep1=1./xsep;
    if(fabs(xsep+xminl-xin->p[1])>1.e-3){
	error("Xin is not linearly spaced\n");
    }
    if(xin->ny!=1 || xnew->ny!=1){
	error("Either xin or xnew is in wrong format\n");
    }

    X(mat) *ynew=X(new)(xnew->nx, xnew->ny);
    PMAT(yin, pyin);
    PMAT(ynew, pynew);
    for(long iy=0; iy<ynew->ny; iy++){
	for(long ix=0; ix<ynew->nx; ix++){
	    double xx=((xnew->p[ix])-xminl)*xsep1;
	    if(xx<0 || xx>nmax1){
		pynew[iy][ix]=0;
	    }else{
		long xxm=ifloor(xx);
		double xxw=xx-xxm;
		pynew[iy][ix]=xxw*pyin[iy][xxm+1]+(1.-xxw)*pyin[iy][xxm];
	    }
	}
    }
    return ynew;
}

/**
   Interpolate using log(xin) and log(xnew)
   xin is the coordinate of yin. xnew is the coordinate of the output.
*/
X(mat)* X(interp1log)(dmat *xin, dmat *yin, dmat *xnew){
    long nmax=xin->nx;
    long nmax1=nmax-1;
    double xminl=log10(xin->p[0]);
    double xmaxl=log10(xin->p[nmax-1]);
    double xsep=(xmaxl-xminl)/(double)(nmax1);
    double xsep1=1./xsep;
    if(fabs(xsep+xminl-log10(xin->p[1]))>1.e-3){
	info("xsep=%g, nmax1=%ld, xminl=%g, xmaxl=%g\n",xsep, nmax1, xminl, xmaxl);
	error("Xin is not logrithmicly spaced\n");
    }
    if(xin->ny!=1 || xnew->ny!=1){
	error("Either xin or xnew is in wrong format\n");
    }

    X(mat) *ynew=X(new)(xnew->nx, xnew->ny);
    PMAT(yin, pyin);
    PMAT(ynew, pynew);
    for(long iy=0; iy<ynew->ny; iy++){
	for(long ix=0; ix<ynew->nx; ix++){
	    double xx=(log10(xnew->p[ix])-xminl)*xsep1;
	    if(xx<0 || xx>nmax1){
		pynew[iy][ix]=0;
	    }else{
		long xxm=ifloor(xx);
		double xxw=xx-xxm;
		pynew[iy][ix]=xxw*pyin[iy][xxm+1]+(1.-xxw)*pyin[iy][xxm];
	    }
	}
    }
    return ynew;
}

#ifndef USE_COMPLEX
/**
   embed a ninx*niny matrix in into A with optional rotation by -theta CCW
   (coordinate rotate theta CCW) around the fft center. Used to rotate the PSF
   from x-y to radial-azimuthal coordinate in radial format CCD.  
   \todo{
   merge this definition with cembed in cmat.c
   }
*/
void X(embed)(X(mat) *restrict A, X(mat) *restrict B, const double theta){
    
    const long ninx=B->nx;
    const long niny=B->ny;
    const long noutx=A->nx;
    const long nouty=A->ny;
    memset(A->p, 0, sizeof(T)*noutx*nouty);
    if(fabs(theta)<1.e-10){//no rotation.
	const long skipx=(noutx-ninx)/2;
	const long skipy=(nouty-niny)/2;
	long ixstart=0, ixend=ninx;
	long iystart=0, iyend=niny;
	if(skipx<0){
	    ixstart=-skipx;
	    ixend=ninx+skipx;
	}
	if(skipy<0){
	    iystart=-skipy;
	    iyend=niny+skipy;
	}
	PMAT(A, pA);
	PMAT(B, pB);
	for(long iy=iystart; iy<iyend; iy++){
	    T *outi=&pA[skipy+iy][skipx+ixstart];
	    T *ini =&pB[iy][ixstart];
	    memcpy(outi, ini, sizeof(T)*(ixend-ixstart));
	}
    }else{
	PMAT(A, outs);
	PMAT(B, ins);
	const double ctheta=cos(theta);
	const double stheta=sin(theta);
	double x2,y2;
	double x,y;
	long ninx2=ninx/2;
	long noutx2=noutx/2;
	long niny2=niny/2;
	long nouty2=nouty/2;
	long ix2, iy2;
	for(long iy=0; iy<nouty; iy++){ 
	    y=(double)(iy-nouty2); 
	    for(long ix=0; ix<noutx; ix++){ 
		x=(double)(ix-noutx2); 
		x2=x*ctheta+y*stheta+ninx2; 
		y2=-x*stheta+y*ctheta+niny2; 
		if(x2>0 && x2<ninx-1 && y2>0 && y2<niny-1){ 
		    ix2=ifloor(x2); 
		    iy2=ifloor(y2); 
		    x2=x2-ix2; 
		    y2=y2-iy2; 
		    outs[iy][ix] =
			ins[iy2][ix2]*((1.-x2)*(1.-y2))
			+ins[iy2][ix2+1]*(x2*(1.-y2))
			+ins[iy2+1][ix2]*((1-x2)*y2)
			+ins[iy2+1][ix2+1]*(x2*y2); 
		} 
	    } 
	} 
    }
}

/**
   For each entry in A, call repeatly to collect its histogram, centered at
   center, spaced by spacing, for n bins in total. center if at bin n/2.  */
void X(histfill)(dmat **out, const X(mat)* A,
		 double center, double spacing, int n){
    if(!A || !A->p) return;
    int nn=A->nx*A->ny;
    if(!*out){
	*out=dnew(n,nn);
    }
    PDMAT(*out,Op);
    const T *restrict Ap=A->p;
    const double spacingi=1./spacing;
    const int noff=n/2;
    const int n1=n-1;
    for(long i=0; i<A->nx*A->ny; i++){
	int ind=(int)round(REAL(Ap[i]-center)*spacingi)+noff;
	if(ind<0) ind=0;
	if(ind>n1) ind=n1;
	Op[i][ind]++;
    }
}
#endif

/**
   1D Cubic spline interpolation preparation.
   if x has only 1 column: x is the coordinate. y is the function value. 
   if x has two columns: first column is the coordinate, y is null.

   It is upto the user to make sure that the coordinate is increasingly ordered
   and evenly spaced .

   If the values of a function \f$f(x)\f$ and its derivative are know at x=0,
   and x=1 (normalized coordinate), then the function can be interpolated on the
   interval [0,1] using a third degree polynomial. This is called cubic
   interpolation. The formula of this polynomial can be easily derived.

   A third degree polynomial and its derivative:
   \f[
   f(x)=ax^3+bx^2+cx+d
   \f]
   \f[
   f(x)=3ax^3+2bx+c
   \f]
   The coefficients can be derived from the value and derivatives:
   \f{eqnarray*}{
   a&=&2f(0)-2f(1)+f^\prime (0)+f^\prime(1)\\
   b&=&-3f(0)+3f(1)-2f^\prime(0)-f^\prime(0)\\
   c&=&f^\prime(0)\\
   d&=&f(0)\\
   \f}
   the derivatives can be computed as
   \f{eqnarray*}{
   f^\prime(0)&=&\frac{f(1)-f(-1)}{2}\\
   f^\prime(1)&=&\frac{f(2)-f(0)}{2}\\
   \f}
   so we have the formula
   \f{eqnarray*}{
   a&=&-0.5 f(-1) + 1.5 f(0) - 1.5 f(1) + 0.5 f(2)\\
   b&=&     f(-1) - 2.5 f(0) + 2   f(1) - 0.5 f(2)\\
   c&=&-0.5 f(-1)            + 0.5 f(1)           \\
   d&=&                 f(0)                      \\
   \f}

   for the boundary pints, replace 
   \f[f^\prime(0)=(f(1)-f(-1))/2\f] by
   \f[f^\prime(0)=(f(1)-f(0))\f]
   Otehr type of boundaries are handled in the same way.

   see http://www.paulinternet.nl/?page=bicubicx */
X(mat) *X(spline_prep)(X(mat) *x, X(mat) *y){
    T *px,*py;
    const long nx=x->nx;
    px=x->p;
    switch(x->ny){
    case 1:
	py=y->p;
	break;
    case 2:
	assert(y==NULL);
	py=x->p+nx;
	break;
    default:
	py=NULL;
	error("Invalid input\n");
    }
    X(mat) *coeff=X(new)(4,nx);
    T xsep=(px[nx-1]-px[0])/(nx-1);
    double thres=ABS(xsep)*1.e-5;
  
    PMAT(coeff,pc);
    T ypriv,ynext;
    for(long ix=0; ix<nx-1; ix++){
	if(ABS(px[ix+1]-px[ix]-xsep)>thres){
	    error("The coordinate is not evenly spaced\n");
	}
	if(UNLIKELY(ix==0)){
	    ypriv=2*py[ix]-py[ix+1];
	}else{
	    ypriv=py[ix-1];
	}
	if(UNLIKELY(ix==nx-2)){
	    ynext=2*py[ix+1]-py[ix];
	}else{
	    ynext=py[ix+2];
	}
	pc[ix][0]=-0.5*ypriv+1.5*py[ix]-1.5*py[ix+1]+0.5*ynext;//a
	pc[ix][1]=     ypriv-2.5*py[ix]+2.0*py[ix+1]-0.5*ynext;//b
	pc[ix][2]=-0.5*ypriv           +0.5*py[ix+1];//c
	pc[ix][3]=               py[ix] ;//d
	/*
	  For any point within this bin, with normalized coordinate t (0<t<1);
	  y(t)=a*pow(t,3)+b*pow(t,2)+c*t+d;
	*/
    }
    return coeff;
}
/**
   Evluate the cubic spline represented by nx5 matrix coeff, at location array xnew.
 */
X(mat)* X(spline_eval)(X(mat) *coeff, X(mat)* x, X(mat) *xnew){
    assert(coeff->nx==4);
    const long nx=coeff->ny;
    PMAT(coeff,pc);
    T xmin=x->p[0];
    T xsep1=(double)(nx-1)/(x->p[nx-1]-xmin);
    X(mat) *out=X(new)(xnew->nx, xnew->ny);
    for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
	double xn=REAL((xnew->p[ix]-xmin)*xsep1);
	long xnf=floor(xn);
	if(xnf<0) xnf=0;
	if(xnf>nx-2) xnf=nx-2;
	xn=xn-xnf;
	T xn2=xn*xn;
	T xn3=xn2*xn;
	out->p[ix]=pc[xnf][0]*xn3+pc[xnf][1]*xn2+pc[xnf][2]*xn+pc[xnf][3];
    }
    return out;
}
/**
   Do 1D cubic spline all at once by calling X(spline_prep) and X(spline_evald)
*/
X(mat)* X(spline)(X(mat) *x,X(mat) *y,X(mat) *xnew){
    X(mat) *coeff=X(spline_prep)(x,y);
    X(mat) *out=X(spline_eval)(coeff,x,xnew);
    X(free)(coeff);
    return out;
}
/**
   2D cubic spline interpolation preparation. x is the x coordinate vector of
 the 2-d grid. y is the y coordinate vector of the 2-d grid. z is defined on the
 2-d grid.  It is upto the user to make sure that the coordinate is increasingly
 ordered and evenly spaced .

 The boundaries are handled in the same way is X(spline). i.e. replace 
 \f[f^\prime(0)=(f(1)-f(-1))/2\f] by
 \f[f^\prime(0)=(f(1)-f(0))\f]
 Otehr type of boundaries are handled in the same way.
*/

X(cell)* X(bspline_prep)(X(mat)*x, X(mat)*y, X(mat) *z){
    const long nx=x->nx;
    const long ny=y->nx;
    assert(x->ny==1 && y->ny ==1 && z->nx==nx && z->ny==ny);
    X(cell)*coeff=X(cellnew)(nx,ny);
    PCELL(coeff,pc);
  
    PMAT(z,p);
    T p00,p01,p02,p03,p10,p11,p12,p13,p20,p21,p22,p23,p30,p31,p32,p33;
    for(long iy=0; iy<ny-1; iy++){
	for(long ix=0; ix<nx-1; ix++){
	    if(iy==0){
		if(ix==0){
		    p00=2*(2*p[iy][ix]-p[iy][ix+1])-(2*p[iy+1][ix]-p[iy+1][ix+1]);//from a
		}else{
		    p00=2*p[iy][ix-1]-p[iy+1][ix-1];//from b
		}
		p01=2*p[iy][ix]-p[iy+1][ix];
		p02=2*p[iy][ix+1]-p[iy+1][ix+1];
		if(ix==nx-2){
		    p03=2*(p[iy][ix+1]*2-p[iy][ix])-(p[iy+1][ix+1]*2-p[iy+1][ix]);//from n
		}else{
		    p03=2*p[iy][ix+2]-p[iy+1][ix+2];//from m
		}
	    }else{
		if(ix==0){
		    p00=2*p[iy-1][ix]-p[iy-1][ix+1];//a from b
		}else{
		    p00=p[iy-1][ix-1];//b
		}
		p01=p[iy-1][ix];
		p02=p[iy-1][ix+1];
		if(ix==nx-2){
		    p03=p[iy-1][ix+1]*2-p[iy-1][ix];//n from m
		}else{
		    p03=p[iy-1][ix+2];//m
		}
	    }
	    if(ix==0){
		p10=p[iy][ix]*2-p[iy][ix+1];//from c
	    }else{
		p10=p[iy][ix-1];//c
	    }
	    p11=p[iy][ix];
	    p12=p[iy][ix+1];
	    if(ix==nx-2){
		p13=p[iy][ix+1]*2-p[iy][ix];//from d
	    }else{
		p13=p[iy][ix+2];//d
	    }
	    if(ix==0){
		p20=p[iy+1][ix]*2-p[iy+1][ix+1];//from e
	    }else{
		p20=p[iy+1][ix-1];//e
	    }
	    p21=p[iy+1][ix];
	    p22=p[iy+1][ix+1];
	    if(ix==nx-2){
		p23=p[iy+1][ix+1]*2-p[iy+1][ix];//from f
	    }else{
		p23=p[iy+1][ix+2];//f
	    }
	    if(iy==ny-2){
		if(ix==0){
		    p30=2*(p[iy+1][ix]*2-p[iy+1][ix+1])-(p[iy][ix]*2-p[iy][ix+1]);//from h
		}else{
		    p30=2*p[iy+1][ix-1]-p[iy][ix-1];//from g
		}
		p31=2*p[iy+1][ix]-p[iy][ix];
		p32=2*p[iy+1][ix+1]-p[iy][ix+1];
		if(ix==nx-2){
		    p33=2*(2*p[iy+1][ix+1]-p[iy+1][ix])-(2*p[iy][ix+1]-p[iy][ix]);//from j
		}else{
		    p33=2*p[iy+1][ix+2]-p[iy][ix+2];//from i
		}
	    }else{
		if(ix==0){
		    p30=p[iy+2][ix]*2-p[iy+2][ix+1];//h from g
		}else{
		    p30=p[iy+2][ix-1];//g
		}
		p31=p[iy+2][ix];
		p32=p[iy+2][ix+1];
		if(ix==nx-2){
		    p33=2*p[iy+2][ix+1]-p[iy+2][ix];//j from i
		}else{
		    p33=p[iy+2][ix+2];//i
		}
	    }
	    pc[iy][ix] = X(new)(4,4);
	    PMAT(pc[iy][ix],ppc);
	    ppc[0][0] = p11;
	    ppc[0][1] = -.5*p10 + .5*p12;
	    ppc[0][2] = p10 - 2.5*p11 + 2*p12 - .5*p13;
	    ppc[0][3] = -.5*p10 + 1.5*p11 - 1.5*p12 + .5*p13;
	    ppc[1][0] = -.5*p01 + .5*p21;
	    ppc[1][1] = .25*p00 - .25*p02 - .25*p20 + .25*p22;
	    ppc[1][2] = -.5*p00 + 1.25*p01 - p02 + .25*p03 + .5*p20 - 1.25*p21 + p22 - .25*p23;
	    ppc[1][3] = .25*p00 - .75*p01 + .75*p02 - .25*p03 - .25*p20 + .75*p21 - .75*p22 + .25*p23;
	    ppc[2][0] = p01 - 2.5*p11 + 2*p21 - .5*p31;
	    ppc[2][1] = -.5*p00 + .5*p02 + 1.25*p10 - 1.25*p12 - p20 + p22 + .25*p30 - .25*p32;
	    ppc[2][2] = p00 - 2.5*p01 + 2*p02 - .5*p03 - 2.5*p10 + 6.25*p11 - 5*p12 + 1.25*p13 + 2*p20 - 5*p21 + 4*p22 - p23 - .5*p30 + 1.25*p31 - p32 + .25*p33;
	    ppc[2][3] = -.5*p00 + 1.5*p01 - 1.5*p02 + .5*p03 + 1.25*p10 - 3.75*p11 + 3.75*p12 - 1.25*p13 - p20 + 3*p21 - 3*p22 + p23 + .25*p30 - .75*p31 + .75*p32 - .25*p33;
	    ppc[3][0] = -.5*p01 + 1.5*p11 - 1.5*p21 + .5*p31;
	    ppc[3][1] = .25*p00 - .25*p02 - .75*p10 + .75*p12 + .75*p20 - .75*p22 - .25*p30 + .25*p32;
	    ppc[3][2] = -.5*p00 + 1.25*p01 - p02 + .25*p03 + 1.5*p10 - 3.75*p11 + 3*p12 - .75*p13 - 1.5*p20 + 3.75*p21 - 3*p22 + .75*p23 + .5*p30 - 1.25*p31 + p32 - .25*p33;
	    ppc[3][3] = .25*p00 - .75*p01 + .75*p02 - .25*p03 - .75*p10 + 2.25*p11 - 2.25*p12 + .75*p13 + .75*p20 - 2.25*p21 + 2.25*p22 - .75*p23 - .25*p30 + .75*p31 - .75*p32 + .25*p33;

	}
    }
    return coeff;
}

/**
   Evaluate 2D cubic spline at location defined 2-d arrays by xnew, ynew
*/
X(mat) *X(bspline_eval)(X(cell)*coeff, X(mat) *x, X(mat) *y, X(mat) *xnew, X(mat) *ynew){
    const long nx=x->nx;
    const long ny=y->nx;
    T xmin=x->p[0];
    T ymin=y->p[0];
    T xsep1=(double)(nx-1)/(x->p[nx-1]-xmin);
    T ysep1=(double)(ny-1)/(y->p[ny-1]-ymin);
    assert(xnew->nx == ynew->nx && xnew->ny == ynew->ny);
    X(mat)*znew=X(new)(xnew->nx, xnew->ny);
    PCELL(coeff,pc);
    for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
	double xm=REAL((xnew->p[ix]-xmin)*xsep1);
	long xmf=floor(xm);
	if(xmf<0) xmf=0;
	if(xmf>nx-2) xmf=nx-2;
	xm=xm-xmf;

	double ym=REAL((ynew->p[ix]-ymin)*ysep1);
	long ymf=floor(ym);
	if(ymf<0) ymf=0;
	if(ymf>ny-2) ymf=ny-2;
	ym=ym-ymf;
	
	T xm2=xm *xm;
	T xm3=xm2*xm;
	T ym2=ym *ym;
	T ym3=ym2*ym;
	PMAT(pc[ymf][xmf],ppc);
	znew->p[ix]= ppc[0][0] + ppc[0][1] * xm + ppc[0][2] * xm2 + ppc[0][3] * xm3 +
	    ppc[1][0] * ym + ppc[1][1] * ym * xm + ppc[1][2] * ym * xm2 + ppc[1][3] * ym * xm3 +
	    ppc[2][0] * ym2 + ppc[2][1] * ym2 * xm + ppc[2][2] * ym2 * xm2 + ppc[2][3] * ym2 * xm3 +
	    ppc[3][0] * ym3 + ppc[3][1] * ym3 * xm + ppc[3][2] * ym3 * xm2 + ppc[3][3] * ym3 * xm3;

    }
    return znew;
}
/**
   Do a component wise log10 on each element of A.
*/
void X(cwlog10)(X(mat) *A){
    double ratio=1./log(10);
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=Y(log)(A->p[i])*ratio;
    }
}
/**
   Embeding an OPD defined on loc to another array. *out is dmat or cmat depends
   on iscomplex. Do the embeding using locstat to have best speed.
   reverse = 0 : from oin to out: out=out*alpha+in*beta
   reverse = 1 : from out to oin: in=in*beta+out*alpha
*/
void X(embed_locstat)(X(mat) **restrict out, double alpha,
		      loc_t *restrict loc, 
		      double *restrict oin, double beta, int reverse){
    loc_create_stat(loc);
    locstat_t *restrict locstat=loc->stat;
    if(!*out){
	if(reverse == 0){
	    *out=X(new)(locstat->nrow, locstat->ncol);
	}else{
	    error("For reverse embedding the array needs to be non-empty\n");
	}
    }else{
	if((*out)->nx < locstat->nrow || (*out)->ny < locstat->ncol){
	    error("Preallocated array %ldx%ld is too small, we need %ldx%ld\n",
		  (*out)->nx, (*out)->ny, locstat->nrow, locstat->ncol);
	}
    }
    PMAT(*out, p);
    double dx1=1./locstat->dx;
    long xoff0=((*out)->nx - locstat->nrow +1)/2;
    long yoff0=((*out)->ny - locstat->ncol +1)/2;

    for(long icol=0; icol<locstat->ncol; icol++){
	long xoff=(long)round((locstat->cols[icol].xstart-locstat->xmin)*dx1);
	long yoff=(long)round((locstat->cols[icol].ystart-locstat->ymin)*dx1);
	long pos1=locstat->cols[icol].pos;
	long pos2=locstat->cols[icol+1].pos;
	T *restrict dest=&p[yoff+yoff0][xoff+xoff0];
	if(!reverse){
	    if(oin){
		const double *restrict oin2=oin+pos1;
		if(fabs(alpha)>EPS){
		    for(long ix=0; ix<pos2-pos1; ix++){
			dest[ix]=dest[ix]*alpha+oin2[ix]*beta;
		    }
		}else{
		    for(long ix=0; ix<pos2-pos1; ix++){
			dest[ix]=oin2[ix]*beta;
		    }
		}
	    }else{
		if(fabs(alpha)>EPS){
		    for(long ix=0; ix<pos2-pos1; ix++){
			dest[ix]=dest[ix]*alpha+beta;
		    }
		}else{
		    for(long ix=0; ix<pos2-pos1; ix++){
			dest[ix]=beta;
		    }
		}
	    }
	}else{
	    double *restrict oin2=oin+pos1;
	    if(fabs(beta)>EPS){
		for(long ix=0; ix<pos2-pos1; ix++){
		    oin2[ix]=oin2[ix]*beta+alpha*REAL(dest[ix]);
		}
	    }else{
		for(long ix=0; ix<pos2-pos1; ix++){
		    oin2[ix]=alpha*REAL(dest[ix]);
		}
	    }
	}
    }
}
/**
   Calculate number of pixels having values larger than or equal to half of
maximum. Useful to compute fwhm. */
long X(fwhm)(X(mat) *A){
    if(!A) return 0;
    double hm=0.5*X(max)(A);
    long fwhm=0;
    for(long ix=0; ix<A->nx*A->ny; ix++){
	if(ABS(A->p[ix])>=hm){
	    fwhm++;
	}
    }
    return fwhm;
}
