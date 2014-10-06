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

/*
  Notice that this file compiles in two modes. 
  1) Compiled alone for T
  2) Compiled by embedding into cX(sp) for complex types. So 

  For any nonstatic function here abc, put a line in sparsedef.sh
  #define abc cabc
*/

#include <math.h>
#include <search.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../sys/sys.h" 

#include "mathmisc.h"
#include "type.h"
#include "mathdef.h"
#include "suitesparse.h"
#include "defs.h"
#include "suitesparse.c"

/**
   Create a nx*ny X(sp) matrix with memory for nmax max
   elements allocated.
*/
X(sp)* X(spnew)(long nx, long ny, long nzmax){
    X(sp) *sp;
    if(nx<0) nx=0;
    if(ny<0) ny=0;
    sp = calloc(1, sizeof(X(sp)));
    sp->id=M_SPT;
    if(nzmax>0){
	sp->p=malloc((ny+1)*sizeof(spint));
	sp->i=malloc(nzmax*sizeof(spint));
	sp->x=malloc(nzmax*sizeof(T));
    }
    sp->m=nx;
    sp->n=ny;
    sp->nzmax=nzmax;
    sp->nz=-1;
    sp->nref=calloc(1,sizeof(int));
    sp->nref[0]=1;
    return sp;
}

/**
   reference a sparse object.
*/
X(sp) *X(spref)(X(sp) *A){
    if(!A) return NULL;
    X(sp) *out = calloc(1, sizeof(X(sp)));
    if(!A->nref){
	A->nref=calloc(1, sizeof(int));
	A->nref[0]=1;
    }
    memcpy(out,A,sizeof(X(sp)));
    atomicadd(out->nref,1);
    return out;
}

/**
   copy a X(sp) matrix to another.
*/
X(sp) *X(spdup)(const X(sp) *A){
    if(!A) return NULL;
    long nmax=A->p[A->n];
    X(sp) *out;
    out=X(spnew)(A->m, A->n, nmax);
    memcpy(out->p, A->p, sizeof(spint)*(A->n+1));
    memcpy(out->i, A->i, sizeof(spint)*nmax);
    memcpy(out->x, A->x, sizeof(T)*nmax);
    return out;
}

/**
   Create a new X(sp) matrix of the same size as A.
*/
X(sp) *X(spnew2)(const X(sp) *A){
    return X(spnew)(A->m, A->n, A->p[A->n]);
}
/**
   Create a new X(sp) matrix and fill in uniform random
   numbers with filling factor of 'fill'
*/
X(sp)* X(spnewrandu)(int nx, int ny, const T mean, 
		     R fill,rand_t *rstat){
    if(fill>1) fill=1.;
    if(fill<0) fill=0.;
    const long nzmax=nx*ny;
    long nz1=nx*ny*fill*4;
    if(nz1>nzmax) nz1=nzmax;
    X(sp) *A=X(spnew)(nx,ny,nz1);
    spint *pp=A->p;
    spint *pi=A->i;
    T *px=A->x;
    long count=0;
    R thres=1.-fill;
    for(int icol=0; icol<A->n; icol++){
	pp[icol]=count;
	for(int irow=0; irow<A->m; irow++){
	    if(randu(rstat)>thres){
		pi[count]=irow;
		px[count]=RANDU(rstat)*mean;
		count++;
		if(count>nz1){
		    /*check out of bound; */
		    nz1=nz1*2; if(nz1>nzmax) nz1=nzmax;
		    X(spsetnzmax)(A,nz1);
		    /*the pointers may change */
		    pp=A->p;
		    pi=A->i;
		    px=A->x;
		}
	    }
	}
    }
    pp[A->n]=count;
    X(spsetnzmax)(A,count);
    return A;
}

/**
   resize a X(sp) matrix
*/
void X(spsetnzmax)(X(sp) *sp, long nzmax){
    if(sp->nzmax!=nzmax){
	sp->i=realloc(sp->i, sizeof(spint)*nzmax);
	sp->x=realloc(sp->x, sizeof(T)*nzmax);
	sp->nzmax=nzmax;
    }
}
/**
 * free a X(sp) matrix*/
void X(spfree_do)(X(sp) *sp){
    if(!sp) return;
    if(!sp->nref || sp->nref[0]<=1){
	if(sp->nref[0]!=1){
	    warning("nref should nevre be less than 1\n");
	}
	free(sp->x);
	free(sp->p);
	free(sp->i);
	if(sp->nref){
	    free(sp->nref);
	}else{
	    warning("X(sp) was corrected incorrectly\n");
	}
    }else{
	atomicadd(sp->nref, -1);
    }
    free(sp);
}

/**
 * Display a X(sp) array*/
void X(spdisp)(const X(sp) *sp){
    long ic,ir;
    long imax;
    if(sp->nzmax==0){
	info("X(spdisp): All zeros\n");
    }else{
	info("X(spdisp):\n");
	for(ic=0; ic<sp->n; ic++){
	    imax=-1;
	    for(ir=sp->p[ic];ir<sp->p[ic+1];ir++){ 
#ifdef USE_COMPLEX
		printf("(%ld,%ld)=(%g,%g)\n", 
		       (long)sp->i[ir], (long)ic, creal(sp->x[ir]),cimag(sp->x[ir]));
#else		
		printf("(%ld,%ld)=%g\n", (long)sp->i[ir], (long)ic, sp->x[ir]);
#endif
		if(sp->i[ir]>imax){
		    imax=sp->i[ir];
		}else{
		    warning("Wrong order");
		}
	    }
	}
    }
}
/**
 * Check a X(sp) array for wrong orders. Return 1 if is lower triangle, 2 if
 * upper triangle, and 3 if diagonal.*/
int X(spcheck)(const X(sp) *sp){
    int not_lower=0;
    int not_upper=0;
    if(sp){
	long ic,ir;
	long imax;
	for(ic=0; ic<sp->n; ic++){
	    imax=-1;
	    if(sp->p[ic+1]<sp->p[ic]){
		error("p in column %ld is smaller than %ld\n",ic+1,ic);
	    }
	    for(ir=sp->p[ic];ir<sp->p[ic+1];ir++){ 
		if(sp->i[ir]>imax){
		    imax=sp->i[ir];
		}else{
		    warning("Wrong order at column %ld",ic);
		}
		if(sp->i[ir]<ic) not_lower=1;
		if(sp->i[ir]>ic) not_upper=1;
	    }
	    if(imax>=sp->m){
		error("imax=%ld exceeds column size at column %ld\n",imax,ic);
	    }
	}
	if(sp->p[sp->n]!=sp->nzmax){
	    warning("real nzmax is %ld, allocated is %ld\n",(long)sp->p[sp->n],sp->nzmax);
	}
    }
    return (not_lower?0:1) | (not_upper?0:2);
}
/**
 * inplace scale X(sp) matrix elements.*/
void X(spscale)(X(sp) *A, const T beta){
    if(A){
	if(A->nref[0]>1){
	    error("spscale on referenced dsp\n");
	}
	for(long i=0; i<A->p[A->ny]; i++){
	    A->x[i]*=beta;
	}
    }
}
/**
 * inplace scale a X(dspcell) object*/
void X(spcellscale)(X(spcell) *A, const T beta){
    for(int i=0; i<A->nx*A->ny; i++){
	X(spscale)(A->p[i],beta);
    }
}
/**
 * Create a new sparse matrix with diagonal elements set to vec*alpha*/
X(sp)* X(spnewdiag)(long N, T *vec, T alpha){
    X(sp) *out=X(spnew)(N,N,N);
    spint *pp=out->p;
    spint *pi=out->i;
    T *px=out->x;
    long count=0;
    if(vec){
	for(long icol=0; icol<out->n; icol++){
	    pp[icol]=count;
	    pi[count]=icol;
	    px[count]=vec[icol]*alpha;
	    count++;
	}
    }else{
	for(long icol=0; icol<out->n; icol++){
	    pp[icol]=count;
	    pi[count]=icol;
	    px[count]=alpha;
	    count++;
	}
    }
    pp[out->n]=count;
    return out;
}
/**
   Extract diagonal element of A and return
*/
X(mat) *X(spdiag)(const X(sp) *A){
    if(A->m!=A->n){
	error("Only work for square matrix\n");
    }
    X(mat) *out=X(new)(A->m,1);
    for(long icol=0; icol<A->n; icol++){
	for(long irow=A->p[icol]; irow<A->p[icol+1]; irow++){
	    long row=A->i[irow];
	    if(row==icol){
		out->p[icol]=A->x[irow];
	    }
	}
    }
    return out;
}
/**
   Multiply a X(sp) matrix inplace with a diagonal weighting matrix whose
   diagonal values are stored in w.
   W_ii=w_i; W_ij=0 if i!=j
   A=A*W*alpha; 
   W is a diagonal X(sp) matrix. diag(W) is w
   multiply w[i] to all numbers in column[i] 
*/
void X(spmuldiag)(X(sp) *restrict A, const T* w, T alpha){
    if(A && w){
	for(long icol=0; icol<A->n; icol++){
	    const T wi=w[icol]*alpha;
	    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		A->x[ix]*=wi;
	    }
	}
    }
}
typedef struct sp_thread_t {
    const X(sp) *A;
    const T *x;
    T **ytmp;
    T *y;
    T alpha;
    long nthread;
}sp_thread_t;

static void X(spmulvec_thread_do_mul)(thread_t *info){
    sp_thread_t *data=info->data;
    const X(sp) *A = data->A;
    const T *restrict x = data->x;
    T *restrict y = data->ytmp[info->ithread]=calloc(A->m, sizeof(T));
    for(long icol=info->start; icol<info->end; icol++){
    	for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
	    y[A->i[ix]]+=A->x[ix]*x[icol];
	}
    }
}
static void X(spmulvec_thread_do_acc)(thread_t *info){
    sp_thread_t *data=info->data;
    T *y    = data->y;
    T **ytmp= data->ytmp;
    
    T alpha = data->alpha;
    long nc = data->nthread;
   
    for(long ix=info->start; ix<info->end; ix++){
	T tmp=0;
	for(long ic=0; ic<nc; ic++){
	    tmp+=ytmp[ic][ix];
	}
	y[ix]+=alpha*tmp;
    }
    
}
/**
   Multiply a sparse with a vector using multithread. Speed up is not signicant
because need to allocate new memory.  */
void X(spmulvec_thread)(T *restrict y, const X(sp) *A, 
			const T * restrict x, T alpha, int nthread){
    if(!A || !x) return; 	
    assert(y);
  
    sp_thread_t data;
    data.A=A;
    data.y=y;
    data.x=x;
    data.alpha=alpha;
    data.ytmp=calloc(nthread,sizeof(T*));
    data.nthread=nthread;
    thread_t mul[nthread];
    thread_t acc[nthread];
    thread_prep(mul, 0, A->n, nthread, X(spmulvec_thread_do_mul), &data);
    CALL_THREAD(mul, 1);
    thread_prep(acc, 0, A->m, nthread, X(spmulvec_thread_do_acc), &data);
    CALL_THREAD(acc, 1);
    for(int ithread=0; ithread<nthread; ithread++){
	free(data.ytmp[ithread]);
    }
    free(data.ytmp);
   
}
/**
 * sparse matrix multiply with a vector*/
void X(spmulvec)(T *restrict y, const X(sp) *A, 
		 const T * restrict x, T alpha){
    if(A && x){
	long icol, ix;
	assert(y);

	for(icol=0; icol<A->n; icol++){
	    for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		y[A->i[ix]]+=alpha*A->x[ix]*x[icol];
	    }
	    
	}
    }
}
static void X(sptmulvec_thread_do)(thread_t *info){
    const long icol_min=info->start;
    const long icol_max=info->end;
    /*It takes 0.0079 s to long the latest threads with thread_pool!!! too bad. */
    sp_thread_t *data=info->data;
    const X(sp) *A = data->A;
    const T *restrict x = data->x;
    T *restrict y = data->y;
    const T alpha=data->alpha;
   
    for(long icol=icol_min; icol<icol_max; icol++){
	for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
	    y[icol]+=alpha*CONJ(A->x[ix])*x[A->i[ix]];
	}
    }
    
    return;
}
/**
 * Threaded version of sptmulvec*/
void X(sptmulvec_thread)(T *restrict y, const X(sp) *A, 
			 const T * restrict x,const T alpha){
    if(!A || !x) return;
    assert(y);
    sp_thread_t data;
    data.A=A;
    data.y=y;
    data.x=x;
    data.alpha=alpha;
    thread_t mul[A->n];
    thread_prep(mul, 0, A->n, A->n, X(sptmulvec_thread_do), &data);
    CALL_THREAD(mul, 1);
  
}
/**
 * Multiply transpose of a sparse matrix with a vector*/
void X(sptmulvec)(T *restrict y, const X(sp) *A, 
		  const T * restrict x,const T alpha){
    if(A && x){
	/*y=y+alpha*A'*x; */
	assert(y);
	long icol, ix;

	for(icol=0; icol<A->n; icol++){
	    for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		y[icol]+=alpha*CONJ(A->x[ix])*x[A->i[ix]];
	    }
	   
	}
    }
}
/**
 * Multiply a sparse matrix with the real part of a complex vector*/
void X(spmulcreal)(T *restrict y, const X(sp) *A, 
		   const RI * restrict x, T alpha){
    /*y=y+alpha*A*creal(x); */
    if(A && x){
	long icol, ix;

	for(icol=0; icol<A->n; icol++){
	    for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		y[A->i[ix]]+=alpha*A->x[ix]*creal(x[icol]);
	    }
	}

    }
}

/**
   y=y+alpha*x*A;
   implemented by transposing x,y index in sptmulmat implementation
   TESTED OK.
*/
void X(mulsp)(X(mat) **yout, const X(mat) *x,const X(sp) *A, const T alpha){
    if(A&&x){
	long icol, ix;
	X(init)(yout, x->nx, A->n);
	X(mat) *y=*yout;
	assert(x->nx==y->nx && x->ny==A->m);
	if(x->nx==1){
	    X(sptmulvec)(y->p, A, x->p, alpha);
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
 * Multiply a sparse matrix X(sp) with a dense matrix X(mat): y+=alpha*A*x;
 */
void X(spmulmat)(X(mat) **yout, const X(sp) *A, const X(mat) *x, 
		 const T alpha){
    if(A&&x){
	/* y=y+alpha*A*x; */
	long icol, ix;
	if(!*yout){
	    *yout=X(new)(A->m, x->ny); 
	}else{
	    assert((*yout)->nx==A->m);
	}
	X(mat) *y=*yout;
	assert(x->ny==y->ny && A->n==x->nx);
	if(x->ny==1){
	    X(spmulvec)(y->p, A, x->p,  alpha);
	}else{
	    int jcol;
	    PMAT(y,Y);
	    PMAT(x,X);
	 
	    for(icol=0; icol<A->n; icol++){
		for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    for(jcol=0; jcol<y->ny; jcol++){
			Y[jcol][A->i[ix]]+=alpha*A->x[ix]*X[jcol][icol];
		    }
		}
	    }
	}
    }
}
/**
   y=y+alpha*A'*x;
*/
void X(sptmulmat)(X(mat) **yout, const X(sp) *A, const X(mat) *x, const T alpha){
    if(A&&x){
	long icol, ix;
	if(!*yout){
	    *yout=X(new)(A->n, x->ny);
	}else{
	    assert((*yout)->nx==A->n);
	}
	X(mat) *y=*yout;
	assert(x->ny==y->ny && A->m==x->nx);
	if(x->ny==1){
	    X(sptmulvec)(y->p, A, x->p, alpha);
	}else{
	    int jcol;
	    PMAT(x,X);
	    PMAT(y,Y);
	  
	    for(icol=0; icol<A->n; icol++){
		for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    for(jcol=0; jcol<y->ny; jcol++){
			Y[jcol][icol]+=alpha*CONJ(A->x[ix])*X[jcol][A->i[ix]];
		    }
		}
	    }
	}
    }
}
/**
 * Multiply two matrices with weighting by sparse matrix.  return y'*(A*x)*/
T X(spwdinn)(const X(mat) *y, const X(sp) *A, const X(mat) *x){
    /*X(sp) weighted ddot. */
    /*computes y'*(A*x). x,y are vectors */
    T res=0;
    if(x && y){
	if(A){
	    assert(x->ny==1 && y->ny==1 && A->m==y->nx && A->n==x->nx);
	    for(int icol=0; icol<A->n; icol++){
		for(int ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    res+=y->p[A->i[ix]]*A->x[ix]*x->p[icol];
		}
	    }
	}else{
	    res=X(inn)(x,y);
	}
    }
    return res;
}
/**
 * Multiply two cell arrays with weighting by sparse matrix*/
T X(spcellwdinn)(const X(cell) *y, const X(spcell) *A, const X(cell) *x){
    /*computes y'*(A*x) */
    T res=0;
    if(x && y){
	if(A){
	    assert(x->ny==1 && y->ny==1 && A->nx==y->nx && A->ny==x->nx);
	    /*PSPCELL(A,Ap); */
	    X(sp) *(*Ap)[A->nx]=(X(sp) *(*)[A->nx])A->p;
	    for(int iy=0; iy<A->ny; iy++){
		for(int ix=0; ix<A->nx; ix++){
		    res+=X(spwdinn)(y->p[ix], Ap[iy][ix], x->p[iy]);
		}
	    }
	}else{
	    res = X(cellinn)(x,y);
	}
    }
    return res;
}
/**
   Multiply a cell with a sparse cell.

   \f$C0+=A*B*alpha\f$.
*/

void X(cellmulsp)(X(cell) **C0, const X(cell) *A, const X(spcell) *B, R alpha){
    if(!A || !B) return;
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    const char trans[2]="nn";
    if(trans[0]=='n'||trans[0]=='N'){
	nx=A->nx; 
	ax=1; az=A->nx;
	nz=A->ny;
    }else{ 
	nx=A->ny;
	az=1; ax=A->nx;
	nz=A->nx;
    }
    if(trans[1]=='n'||trans[1]=='N'){
	ny=B->ny; 
	bz=1; by=B->nx;
	if(nz!=B->nx) error("mismatch\n");
    }else{
	ny=B->nx;
	by=1; bz=B->nx;
	if(nz!=B->ny) error("mismatch\n");
    }
    if(!*C0){
	*C0=cellnew(nx,ny);
    }
    X(cell) *C=*C0;
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    for(int iz=0; iz<nz; iz++){
		if(A->p[ix*ax+iz*az] && B->p[iz*bz+iy*by]){
		    X(mulsp)(&C->p[ix+iy*nx],A->p[ix*ax+iz*az], 
			     B->p[iz*bz+iy*by],alpha);
		}
	    }
	}
    }
}

#define SPCELLMULMAT_DO						\
    for(int iz=0; iz<nz; iz++){					\
	if(trans){						\
	    X(sptmulmat)(&C->p[ix+iy*nx], A->p[iz+ix*A->nx],	\
			 B->p[iz+iy*B->nx], alpha);		\
	}else{							\
	    X(spmulmat)(&C->p[ix+iy*nx], A->p[ix+iz*A->nx],	\
			B->p[iz+iy*B->nx], alpha);		\
	}							\
    }
static void X(spcellmulmat2)(X(cell) **C0, const X(spcell)*A, const X(cell)*B, 
		      const T alpha, const int trans){
    if(!A || !B) return;
    int nx,ny,nz;
    if(trans){
	nx=A->ny; 
	nz=A->nx;
    }else{
	nx=A->nx; 
	nz=A->ny;
    }
    ny=B->ny; 
    assert(nz==B->nx);
    
    if(!*C0){
	*C0=cellnew(nx,ny);
    }
    X(cell) *C=*C0;
    assert(C->nx==nx && C->ny==ny);
    
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    SPCELLMULMAT_DO;
	}
    }
}
/**
 * Multiply a dspcell with a dense cell: C=C+A*B*alpha*/
void X(spcellmulmat)(X(cell) **C, const X(spcell)*A, const X(cell)*B, const T alpha){
    return X(spcellmulmat2)(C,A,B,alpha,0);
}
/**
 * C=C+A'*B*alpha*/
void X(sptcellmulmat)(X(cell) **C, const X(spcell)*A, const X(cell)*B, const T alpha){
    return X(spcellmulmat2)(C,A,B,alpha,1);
}
typedef struct{
    int nx;
    int nz;
    int trans;
    const X(spcell) *A;
    const X(cell) *B;
    X(cell) *C;
    T alpha;
}dspcell_thread_t;

static void X(spcellmulmat_thread_do)(thread_t *info){
    dspcell_thread_t *data=info->data;
    const X(spcell) *A=data->A;
    const X(cell) *B=data->B;
    X(cell) *C=data->C;
    const T alpha=data->alpha;
    const int nx=data->nx;
    const int nz=data->nz;
    const int trans=data->trans;
    for(int ic=info->start; ic<info->end; ic++){
	int iy=ic/nx;
	int ix=ic-nx*iy;
	SPCELLMULMAT_DO;
    }
}
/**
   A common routine that can handle dspcellmulmat_thread and sptcellmulmat_thread
 */
static void X(spcellmulmat_thread2)(X(cell) **C0, const X(spcell)*A, 
				    const X(cell)*B, const T alpha, 
				    const int trans){
    /*C=C+A*B*alpha */
    if(!A || !B) return;
    int nx,ny,nz;
    if(trans){
	nx=A->ny; 
	nz=A->nx;
    }else{
	nx=A->nx; 
	nz=A->ny;
    }
    ny=B->ny; 
    assert(nz==B->nx);
    
    if(!*C0){
	*C0=cellnew(nx,ny);
    }
    X(cell) *C=*C0;
    if(C->nx!=nx || C->ny!=ny) error("Mismatch\n");
    int ntot=nx*ny;
    dspcell_thread_t data;
    data.nx=nx;
    data.nz=nz;
    data.trans=trans;
    data.A=A;
    data.B=B;
    data.C=C;
    data.alpha=alpha;
    thread_t mul[ntot];
    thread_prep(mul, 0, ntot, ntot, X(spcellmulmat_thread_do), &data);
    CALL_THREAD(mul, 1);
}
/**
 * threaded version of X(spcellmulmat)*/
void X(spcellmulmat_thread)(X(cell) **C, const X(spcell)*A, 
			    const X(cell)*B, const T alpha){
    X(spcellmulmat_thread2)(C,A,B,alpha,0);
}
/**
 * threaded version of X(sptcellmulmat*/
void X(sptcellmulmat_thread)(X(cell) **C, const X(spcell)*A, 
			     const X(cell)*B, const T alpha){
    X(spcellmulmat_thread2)(C,A,B,alpha,1);
}
/**
   Data for use in dspcellmulmat_each
 */
typedef struct{
    X(cell) *xout;
    X(spcell) *A;
    X(cell) *xin; 
    T alpha;
    int trans;
}EACH_T;
/**
   Multiply each cell in dspcell with each cell in dcell.
*/
static void X(spcellmulmat_each_do)(thread_t *info){
    EACH_T *data=info->data;
    if(data->trans){
	for(int ic=info->start; ic<info->end; ic++){
	    X(sptmulmat)(&data->xout->p[ic], data->A->p[ic], data->xin->p[ic], data->alpha);
	}
    }else{
	for(int ic=info->start; ic<info->end; ic++){
	    X(spmulmat)(&data->xout->p[ic], data->A->p[ic], data->xin->p[ic], data->alpha);
	}
    }
}
/**
   Threaded multiply each cell in dspcell with each cell in dcell.
*/
void X(spcellmulmat_each)(X(cell) **xout, X(spcell) *A, X(cell) *xin, 
			  T alpha, int trans){
    if(!*xout){
	*xout=cellnew(xin->nx, xin->ny);
    }
    assert(xin->ny==1);
    EACH_T data;
    data.xout=*xout;
    data.A=A;
    data.xin=xin;
    data.alpha=alpha;
    data.trans=trans;
    int tot=data.xout->nx;
    thread_t info[tot];
    thread_prep(info, 0, tot, tot, X(spcellmulmat_each_do), &data);
    CALL_THREAD(info,1);
}

/**
 * Convert sparse matrix into dense matrix and add to output:
 * out0=out0+full(A)*alpha*/
void X(spfull)(X(mat) **out0, const X(sp) *A, const T alpha){
    if(!A)
	return;
    /**
       add A*f to dense matrix located in p;
    */
    long nx=A->m;
    long icol,ix,irow;
    if(!*out0){
	*out0=X(new)(A->m, A->n);
    }
    X(mat) *out=*out0;
    assert(out->nx==A->m && out->ny==A->n);
    PMAT(out,pp);
    for(icol=0; icol<A->n; icol++){
	for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
	    irow=A->i[ix];
	    if(irow>=nx)
		error("invalid row:%ld, %ld",irow,nx);
	    pp[icol][irow]+=alpha*A->x[ix];
	}
    }
}
X(sp) *X(2sp)(X(mat)*A){
    if(!A) return 0;
    X(sp) *out=X(spnew)(A->nx, A->ny, A->nx*A->ny);
    out->p[0]=0;
    for(long icol=0; icol<A->ny; icol++){
	out->p[icol+1]=(icol+1)*A->nx;
	for(long irow=0; irow<A->nx; irow++){
	    long ix=out->p[icol]+irow;
	    out->i[ix]=irow;
	}
    }
    memcpy(out->x, A->p, sizeof(T)*A->nx*A->ny);
    return out;
}
/** 
 * Convert the transpose of a sparse matrix into dense matrix and add to output:
 * out0=out0+full(A')*alpha;*/
void X(sptfull)(X(mat) **out0, const X(sp) *A, const T alpha){
    if(!A) return;
    /**
       add A*f to dense matrix located in p;
    */
    long nx=A->m;
    long icol,ix,irow;
    if(!*out0){
	*out0=X(new)(A->n, A->m);
    }
    X(mat) *out=*out0;
    assert(out->nx==A->n && out->ny==A->m);
    PMAT(out,pp);
    for(icol=0; icol<A->n; icol++){
	for(ix=A->p[icol]; ix<A->p[icol+1]; ix++){
	    irow=A->i[ix];
	    if(irow>=nx)
		error("invalid row:%ld, %ld",irow,nx);
	    pp[irow][icol]+=alpha*A->x[ix];
	}
    }
}
/**
 * Convert sparse cell to dense matrix cell: out0=out0+full(A)*alpha*/
void X(spcellfull)(X(cell) **out0, const X(spcell) *A, const T alpha){
    if(!A) return;
    X(cell) *out=*out0;
    if(!out){
	out=*out0=cellnew(A->nx, A->ny);
    }else{
	assert(out->nx==A->nx && out->ny==A->ny);
    }
    PSPCELL(A,pA);
    PCELL(out,pout);
    for(int iy=0; iy<A->ny; iy++){
	for(int ix=0; ix<A->nx; ix++){
	    X(spfull)(&pout[iy][ix], pA[iy][ix], alpha);
	}
    }
}
/**
 * Convert transpose of sparse cell to dense matrix cell: out0=out0+full(A')*alpha*/
void X(sptcellfull)(X(cell) **out0, const X(spcell) *A, const T alpha){
    if(!A) return;
    X(cell) *out=*out0;
    if(!out){
	out=*out0=cellnew(A->ny, A->nx);
    }else{
	assert(out->nx==A->ny && out->ny==A->nx);
    }
    PSPCELL(A,pA);
    PCELL(out, pout);
    for(int iy=0; iy<A->ny; iy++){
	for(int ix=0; ix<A->nx; ix++){
	    X(spfull)(&pout[ix][iy], pA[iy][ix], alpha);
	}
    } 
}
/**
 * Added two sparse matrices: return A*a+B*b*/
X(sp) *X(spadd2)(X(sp) *A,X(sp)*B,T a,T b){
    X(sp) *C=X(ss_add)(A,B,a,b);
    X(ss_dropzeros)(C);
    return C;
}
/**
 * Add a sparse matrix to another: A0=A0+B*/
void X(spadd)(X(sp) **A0, const X(sp) *B){
    /*add B to A. */
    if(B){
	if(!*A0) 
	    *A0=X(spdup)(B);
	else{
	    if((*A0)->m!=B->m || (*A0)->n!=B->n) {
		error("X(sp) matrix mismatch: (%ldx%ld) vs (%ldx%ld\n",
		      (*A0)->m, (*A0)->n, B->m, B->n);
	    }
	    X(sp) *res=X(ss_add)(*A0,B,1.,1.);
	    X(ss_dropzeros)(res);
	    X(spfree)(*A0);
	    *A0=res; res=0;
	}
    }
}
/**
 * Add a sparse cell to another: A0=A0+B */
void X(spcelladd)(X(spcell) **A0, const X(spcell) *B){
    if(B){
	if(!*A0){
	    *A0=cellnew(B->nx, B->ny);
	}
	for(int i=0; i<B->nx*B->ny; i++){
	    X(spadd)(&((*A0)->p[i]), B->p[i]);
	}
    }
}
/**
   Add alpha times identity to a sparse matrix
*/
void X(spaddI)(X(sp) **A, R alpha){
    assert((*A)->m==(*A)->n);
    X(sp) *B=X(spnewdiag)((*A)->m,NULL,alpha);
    X(spadd)(A,B);
    X(spfree)(B);
}
/**
   Add alpha times identity to sparse array.
 */
void X(spcelladdI)(X(spcell) *A, R alpha){
    assert(A->nx==A->ny);
    for(int ii=0; ii<A->ny; ii++){
	X(spaddI)(&A->p[ii+ii*A->nx],alpha);
    }
}
/**
 * Transpose a sparse array*/
X(sp) *X(sptrans)(const X(sp) *A){
    if(!A) return NULL;
    X(sp) *res=X(ss_transpose)(A,1);
    X(ss_dropzeros)(res);
    return res;
}
/**
 * Multiply two sparse arrays: return A*B*/
X(sp) *X(spmulsp)(const X(sp) *A, const X(sp) *B){		
    /*return C=(A*B) */
    if(!A || !B) return NULL;
    X(sp) *C=X(ss_multiply)(A, B);
    X(ss_dropzeros)(C);
    return C;
}
/**
 * Multiply the transpose of a sparse with another: return A'*B*/
X(sp) *X(sptmulsp)(const X(sp) *A, const X(sp) *B){
    /*return A'*B; */
    /*fixme : may need to improve this so that tranpose of A is not necessary. */
    X(sp) *At=X(sptrans)(A);
    X(sp) *C=X(spmulsp)(At, B);
    X(spfree)(At);
    X(ss_dropzeros)(C);
    return C;
}
/**
 * Multiply two sparse arrays and add to the third: C0=C0+A*B*scale*/
void X(spmulsp2)(X(sp) **C0, const X(sp) *A, const X(sp) *B, 
		 const T scale){
    /*return C=C+ alpha*(A*B) */
    if(!A || !B) return;
    X(sp) *res=X(ss_multiply)(A, B);
    if(ABS(scale-1.)>EPS){
	X(spscale)(res, scale);
    }
    if(!*C0) 
	*C0=res;
    else{
	X(spadd)(C0, res);
	X(spfree)(res);
    }
    X(ss_dropzeros)(*C0);
}
/**
 * Multiply two sparse cell*/
X(spcell) *X(spcellmulspcell)(const X(spcell) *A, 
			      const X(spcell) *B, 
			      const T scale){
    /*return C=A*B; */
    if(!A || !B) return NULL;
    if(A->ny!=B->nx) error("mismatch\n");
    X(spcell) *C=cellnew(A->nx, B->ny);
    PSPCELL(A,Ap);
    PSPCELL(B,Bp);
    PSPCELL(C,Cp);
    /*
    X(sp) *(*Ap)[A->nx] = (X(sp) *(*)[A->nx]) A->p;
    X(sp) *(*Bp)[B->nx] = (X(sp) *(*)[B->nx]) B->p;
    X(sp) *(*Cp)[C->nx] = (X(sp) *(*)[C->nx]) C->p;
    */
    for(int iy=0; iy<B->ny; iy++){
	for(int ix=0; ix<A->nx; ix++){
	    Cp[iy][ix]=NULL;
	    for(int iz=0; iz<A->ny; iz++){
		X(spmulsp2)(&Cp[iy][ix],Ap[iz][ix],Bp[iy][iz],scale);
	    }
	}
    }
    return C;
}

/**
 * Transpose a sparse cell*/
X(spcell) *X(spcelltrans)(const X(spcell) *spc){
    if(!spc) return NULL;
    long nx,ny;
    nx=spc->nx;
    ny=spc->ny;
    X(spcell) *spct=cellnew(ny,nx);
    
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    spct->p[iy+ix*ny]=X(sptrans)(spc->p[ix+iy*nx]);
	}
    }
    return spct;
}
/**
 * Free a sparse cell data*/
void X(spcellfree_do)(X(spcell) *spc){
    if(!spc || !spc->p) return;
    for(int ix=0; ix<spc->nx*spc->ny; ix++){
	X(spfree)(spc->p[ix]);
    }
    free(spc->p);
    free(spc);
}
/**
 * Concatenate two sparse array along dim dimension*/
X(sp) *X(spcat)(const X(sp) *A, const X(sp) *B, int dim){
    X(sp) *C=NULL;
    if(dim==0){
	error("Not implemented\n");
	/*
	  |A|
	  |B|
	*/
    }else if(dim==1){
	/*|AB|*/
	if(A->m != B->m){
	    error("X(sp) matrix doesn't match\n");
	}
	const long nzmax=A->p[A->n]+B->p[B->n];
	C=X(spnew)(A->m, A->n+B->n, nzmax);
	memcpy(C->p, A->p, A->n*sizeof(spint));
	memcpy(C->i, A->i, A->p[A->n]*sizeof(spint));
	memcpy(C->x, A->x, A->p[A->n]*sizeof(T));
	memcpy(C->i+A->p[A->n], B->i, B->p[B->n]*sizeof(spint));
	memcpy(C->x+A->p[A->n], B->x, B->p[B->n]*sizeof(T));
	const long Anzmax=A->p[A->n];
	for(long i=0; i<B->n+1; i++){
	    C->p[i+A->n]=Anzmax+B->p[i];
	}
    }else{
	error("Wrong dimension\n");
    }
    return C;
}
/**
 * Concatenate a dspcell to sparse array*/
X(sp) *X(spcell2sp)(const X(spcell) *A){
    /*convert X(spcell) to sparse. */
    if(A->nx*A->ny==1){/*There is a single cell */
	return X(spref)(A->p[0]);
    }
    PSPCELL(A,Ap);
    long nx=0,ny=0,nzmax=0;
    long nnx[A->nx];
    long nny[A->ny];
    for(long ix=0; ix<A->nx; ix++){
	for(long iy=0; iy<A->ny; iy++){
	    if(Ap[iy][ix]) {
		nnx[ix]=Ap[iy][ix]->m;
		nx+=Ap[iy][ix]->m;
		break;
	    }
	}
    }
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    if(Ap[iy][ix]) {
		nny[iy]=Ap[iy][ix]->n;;
		ny+=Ap[iy][ix]->n;
		break;
	    }
	}
    }
    for(long i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    nzmax+=A->p[i]->p[A->p[i]->n];
	}
    }
    X(sp) *out=X(spnew)(nx,ny,nzmax);
    long count=0;
    long jcol=0;
    for(long iy=0; iy<A->ny; iy++){
	for(long icol=0; icol<nny[iy]; icol++){
	    out->p[jcol+icol]=count;
	    long kr=0;
	    for(long ix=0; ix<A->nx; ix++){
		if(Ap[iy][ix]){
		    for(long ir=Ap[iy][ix]->p[icol]; 
			ir<Ap[iy][ix]->p[icol+1]; ir++){
			out->x[count]=Ap[iy][ix]->x[ir];
			out->i[count]=Ap[iy][ix]->i[ir]+kr;
			count++;
		    }
		}
		kr+=nnx[ix];
	    }
	}
	jcol+=nny[iy];
    }
    out->p[ny]=count;
    if(count>nzmax){
	error("X(spcell2sp) gets Wrong results. count=%ld, nzmax=%ld\n",count,nzmax);
    }
    /*nzmax maybe smaller than A->p[A->n]  */
    /*because nzmax simply show the slots available. */
    return out;
}

/**
 * Sum elements of sparse array along dimension dim*/
X(mat) *X(spsum)(const X(sp) *A, int dim){
    /*Sum X(sp) matrix along col or row to form a vector */
    X(mat) *v=NULL;
    T *p;
    switch(dim){
    case 1:/*sum along col */
	v=X(new)(1,A->n);
	p=v->p;
	for(int icol=0; icol<A->n; icol++){
	    for(int irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		p[icol]+=A->x[irow];
	    }
	}
	break;
    case 2:/*sum along row */
	v=X(new)(A->m,1);
	p=v->p;
	for(int icol=0; icol<A->n; icol++){
	    for(int irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		p[A->i[irow]]+=A->x[irow];
	    }
	}
	break;
    default:
	error("Invalid\n");
    }
    return v;
}
/**
 * Sum abs of elements of sparse array along dimension dim*/
X(mat) *X(spsumabs)(const X(sp) *A, int col){
    X(mat) *v=NULL;
    T *p;
    switch(col){
    case 1:/*sum along col */
	v=X(new)(1,A->n);
	p=v->p;
	for(int icol=0; icol<A->n; icol++){
	    for(int irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		p[icol]+=ABS(A->x[irow]);
	    }
	}
	break;
    case 2:/*sum along row */
	v=X(new)(A->m,1);
	p=v->p;
	for(int icol=0; icol<A->n; icol++){
	    for(int irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		p[A->i[irow]]+=ABS(A->x[irow]);
	    }
	}
	break;
    default:
	error("Invalid\n");
    }
    return v;
}
/**
   Clean up a sparse array by dropping zeros*/
void X(spclean)(X(sp) *A){
    X(ss_dropzeros)(A);
}
/**
 * Multiply a dspcell with vectors*/
void X(spcellmulvec)(T *restrict yc, const X(spcell) *Ac, 
		     const T * restrict xc, T alpha){
    /*y=y+alpha*A*creal(x); Ac X(sp) cell. xc is vector. */
    if(Ac && xc){
	const T *restrict x=xc;
	for(long icy=0; icy<Ac->ny; icy++){
	    T *y=yc;
	    for(long icx=0; icx<Ac->nx; icx++){
		const X(sp) *A=Ac->p[icx+icy*Ac->nx];
		X(spmulvec)(y,A,x,alpha);
		y+=A->m;
	    }
	    x+=Ac->p[icy*Ac->nx]->n;
	}
    }
}
/**
   Drop elements that are EPS times the largest value.
*/
void X(spdroptol)(X(sp) *A, R thres){
    if(thres<EPS) thres=EPS;
    R maxv;
    X(maxmin)(A->x,A->nzmax,&maxv,NULL);
    X(ss_droptol)(A, maxv*thres);
}
/**
   Drop elements that are EPS times the largest value.
*/
void X(spcelldroptol)(X(spcell) *A, R thres){
    for(int i=0; i<A->nx*A->ny; i++){
	X(spdroptol)(A->p[i], thres);
    }
}

typedef struct spelem{
    int i;
    T x;
}spelem;
static int spelemcmp(const spelem *A, const spelem *B){
    return A->i-B->i;
}
/**
   Make sure the elements are sorted correctly. Does not change the location of
data. can be done without harm. */
void X(spsort)(X(sp) *A){
    if(!A || A->n==0 || A->m==0) return;
    long nelem_max=A->nzmax/A->n*2;
    spelem *col=malloc(nelem_max*sizeof(spelem));
    for(long i=0; i<A->n; i++){
	long nelem=(A->p[i+1]-A->p[i]);
	if(nelem==0) continue;
	if(nelem>nelem_max){
	    nelem_max=nelem;
	    col=realloc(col, nelem_max*sizeof(spelem));
	}
	for(long j=0; j<nelem; j++){
	    col[j].i=A->i[A->p[i]+j];
	    col[j].x=A->x[A->p[i]+j];
	}
	qsort(col, nelem, sizeof(spelem), (int(*)(const void*,const void*))spelemcmp);
	for(long j=0; j<nelem; j++){
	    A->i[A->p[i]+j]=col[j].i;
	    A->x[A->p[i]+j]=col[j].x;
	}
    }
    free(col);

}
/**
   Make sure the elements are sorted correctly.
*/
void X(spcellsort)(X(spcell) *A){
    for(int i=0; i<A->nx*A->ny; i++){
	X(spsort)(A->p[i]);
    }
}
/**
   symmetricize a X(sp) matrix and drop values below a
   threshold. 
*/
void X(spsym)(X(sp) *A){
    X(sp) *B=X(sptrans)(A);
    X(spadd)(&A,B);
    X(spscale)(A,0.5);
    X(spfree)(B);
    X(spdroptol)(A,EPS);
    X(spsort)(A);/*This is important to make chol work. */
}
/**
   symmetricize a X(sp) cell and drop values below a
   threshold. 
*/
void X(spcellsym)(X(spcell) *A){
    X(spcell) *B=X(spcelltrans)(A);
    X(spcelladd)(&A,B);
    X(spcellfree)(B);
    X(spcellscale)(A,0.5);
    X(spcelldroptol)(A,EPS);
    X(spcellsort)(A);
}

/**
   Create a X(sp) convolution operator C with
   C(i,j)=A(i-j);
   A must be very X(sp) with only a view non-zero value otherwise C will be too full.
*/
X(sp) *X(spconvolvop)(X(mat) *A){
    /*First collect statistics on A. */
    long nini=10;
    T *vals=calloc(nini, sizeof(T));
    long *sepx=calloc(nini, sizeof(spint));
    long *sepy=calloc(nini, sizeof(spint));
    long count=0;
    const long nx=A->nx;
    const long ny=A->ny;
    const long nn=nx*ny;
    PMAT(A,PA);
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    if(ABS(PA[iy][ix])>0){
		vals[count]=PA[iy][ix];
		sepx[count]=ix;
		sepy[count]=iy;
		count++;
	    }
	    if(count>=nini){
		nini*=2;
		vals=realloc(vals, sizeof(T)*nini);
		sepx=realloc(sepx, sizeof(spint)*nini);
		sepy=realloc(sepy, sizeof(spint)*nini);
	    }
	}
    }
    if(count>10){
	warning("Number of coupled points %ld is too large\n",count);
    }
    long nsep=count;
    X(sp) *out=X(spnew)(nn,nn,nn*count);
    spint *pp=out->p;
    spint *pi=out->i;
    T *px=out->x;
    count=0;
    long icol=0;
    for(long iiy=0; iiy<A->ny; iiy++){
	for(long iix=0; iix<A->nx; iix++){
	    pp[icol]=count;
	    icol++;
	    for(int irow=0; irow<nsep; irow++){
		long jix=(iix+sepx[irow])%nx;
		long jiy=(iiy+sepy[irow])%ny;
		pi[count]=jix+jiy*nx;
		px[count]=vals[irow];
		count++;
	    }
	}
    }
    pp[nn]=count;
    free(vals);
    free(sepx);
    free(sepy);
    X(spsort)(out);
    return out;
}
/**
   Permute rows of X(sp) matrix A;
   let pi be the inverse order of p. pi=invperm(p);
   if reverse is 0
   B(:,i)=A(p,i); or B(pi,i)=A(:,i)
   if reverse is 1
   B(p,i)=A(:,i); or B(:,i)=A(pi,i);
*/
static X(sp) *X(sppermcol)(const X(sp) *A, int reverse, long *p){
    if(!p) return X(spdup)(A);
    X(sp) *B=X(spnew)(A->m,A->n,A->nzmax);
    long *perm;
    if(reverse){
	perm=p;
    }else{
	perm=invperm(p,A->m);
    }

    for(long icol=0; icol<A->n; icol++){
	B->p[icol]=A->p[icol];
	for(long irow=A->p[icol]; irow<A->p[icol+1]; irow++){
	    long row=A->i[irow];
	    B->i[irow]=perm[row];
	    B->x[irow]=A->x[irow];
	}
    }
    B->p[A->n]=A->p[A->n];
    if(!reverse){
	free(perm);
    }
    return B;
}
/**
   Permute rows and columns of X(sp) matrix A;
*/
X(sp) *X(spperm)(X(sp) *A, int reverse, long *pcol, long *prow){
    X(sp) *out;
    if(pcol){
	out=X(sppermcol)(A,reverse,pcol);
    }else{
	out=X(spref)(A);
    }
    if(prow){
	X(sp) *Ap=X(sptrans)(out);
	X(sp) *App=X(sppermcol)(Ap,reverse,prow);
	X(spfree)(Ap);
	X(spfree)(out);
	out=X(sptrans)(App);
	X(spfree)(App);
    }
    return out;
}

/**
   Invert a SPD X(sp) matrix that is block diagonal with
   block sizes of bs.
*/
#ifndef USE_SINGLE
X(sp) *X(spinvbdiag)(const X(sp) *A, long bs){
    if(A->m!=A->n){
	error("Must be a square matrix\n");
    }
    long nb=A->m/bs;
    X(sp) *B=X(spnew)(A->m, A->n, nb*bs*bs);
    X(mat) *bk=X(new)(bs,bs);
    PMAT(bk,pbk);
    for(long ib=0;ib<nb; ib++){
	long is=ib*bs;/*starting col */
	X(zero)(bk);

	for(long icol=is; icol<is+bs; icol++){
	    for(long irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		long row=A->i[irow];
		long ind=row-is;
		if(ind<0 || ind>=bs){
		    info("solving block %ld\n",ib);
		    error("The array is not block diagonal matrix or not calculated properly.\n");
		}
		pbk[icol-is][ind]=A->x[irow];
	    }
	}
	X(inv_inplace)(bk);
	for(long icol=is; icol<is+bs; icol++){
	    B->p[icol]=icol*bs;
	    for(long irow=0; irow<bs; irow++){
		B->i[B->p[icol]+irow]=irow+is;
		B->x[B->p[icol]+irow]=pbk[icol-is][irow];
	    }
	}
    }
    B->p[A->n]=nb*bs*bs;
    X(free)(bk);
    return B;
}
#endif
/**
   Extrat the diagonal blocks of size bs into cell arrays.
*/
X(cell) *X(spblockextract)(const X(sp) *A, long bs){
    if(A->m!=A->n){
	error("Must be a square matrix\n");
    }
    long nb=A->m/bs;
    X(cell) *out=cellnew(nb,1);
    for(long ib=0;ib<nb; ib++){
	long is=ib*bs;/*starting col */
	out->p[ib]=X(new)(bs,bs);
	PMAT(out->p[ib],pbk);
	for(long icol=is; icol<is+bs; icol++){
	    for(long irow=A->p[icol]; irow<A->p[icol+1]; irow++){
		long row=A->i[irow];
		long ind=row-is;
		if(ind<0 || ind>=bs){
		    info("solving block %ld\n",ib);
		    error("The array is not block diagonal matrix or not calculated property\n");
		}
		pbk[icol-is][ind]=A->x[irow];
	    }
	}
    }
    return out;
}
