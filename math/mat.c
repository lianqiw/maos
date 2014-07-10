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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include "../sys/sys.h"
#include "random.h"
#include "mathdef.h"
#include "defs.h"/*Defines T, X, etc */
/**
   Work horse function that creats the matrix object. if p is
   NULL, memory is allocated. Allocation for X(mat) objects
   should always be created directly or indirectly through
   X(new_do) so that changing data structure is straight
   forward.
*/
static inline X(mat) *X(new_do)(long nx, long ny, T *p, int ref){
    X(mat) *out=calloc(1, sizeof(X(mat)));
    out->id=M_T;
    out->nx=nx;
    out->ny=ny;
    if(ref){/*the data does not belong to us. */
	if(!p){
	    error("When ref is 1, p must not be NULL.\n");
	}
	out->p=p;
    }else{
	if(!p){
	    p=calloc((nx*ny), sizeof(T));
	}
	out->p=p;
	out->nref=calloc(1, sizeof(int));
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
   chunk. the memory is freed when the X(mat) is freed.
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
   check the size of matrix if exist. Otherwise create it
*/
void X(init)(X(mat)**A, long nx, long ny){
    if(*A){
	assert((*A)->nx == nx && (*A)->ny==ny);
    }else{
	*A=X(new)(nx, ny);
    }
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
    int free_extra=0;
    if(A->nref){
	int nref=atomicadd(A->nref, -1);
	if(!nref){
	    if(A->header){
		long count=search_header_num(A->header, "count");
		if(!isnan(count) && count>0){
		    error("deprecated: count=%ld, scaling the data\n", count);
		    X(scale)(A, (R)(1./count));
		}
	    }
	    if(!keepdata && A->p){
		if(A->mmap){/*data is mmap'ed. */
		    mmap_unref(A->mmap);
		}else{
		    free_extra=1;
		    free(A->p);
		}
	    }
	    free(A->nref);
	}else if(nref<0){
	    error("The ref is less than 0. something wrong!!!:%d\n",A->nref[0]);
	}
    }else{
	free_extra=1;
    }
    if(free_extra){
	X(fft_free_plan)(A->fft);
	free(A->header);
    }
    free(A);
}

/**
   Resize a matrix by adding or removing columns or rows. Data is kept
   whever possible.
*/
void X(resize)(X(mat) *A, long nx, long ny){
    if(!A->nref || A->nref[0]>1){
	error("Resizing a referenced vector\n");
    }
    if(A->nx==nx || A->ny==1){
	A->p=realloc(A->p, sizeof(T)*nx*ny);
	if(nx*ny>A->nx*A->ny){
	    memset(A->p+A->nx*A->ny, 0, (nx*ny-A->nx*A->ny)*sizeof(T));
	}
    }else{/*copy memory to preserve data*/
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
X(mat) *X(ref)(const X(mat) *in){
    if(!in) return NULL;
    X(mat) *out=calloc(1, sizeof(X(mat)));
    memcpy(out,in,sizeof(X(mat)));
    if(!in->nref){
	extern quitfun_t quitfun;
	if(quitfun==&default_quitfun){
	    warning_once("Referencing non-referrable data. This may cause error.\n");
	}
    }else{
	atomicadd(in->nref, 1);
    }
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
	/*along x. */
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
	/*along y. */
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
   Compute max, min and sum. Has to handle NAN nicely */
void X(maxmin)(const T *restrict p, long N, R *max, R *min){
    R a,b;
    long i;
    a=-INFINITY;
    b=INFINITY;
    for(i=0; i<N; i++){
	R tmp=MAG(p[i]);
	if(!isnan(tmp)){
	    if(tmp>a) a=tmp;
	    if(tmp<b) b=tmp;
	}
    }
    if(max)*max=a; 
    if(min)*min=b; 
}

/**
   find the maximum value of a X(mat) object
*/
R X(max)(const X(mat) *A){
    R max,min;
    X(maxmin)(A->p, A->nx*A->ny, &max, &min);
    return max;
}

/**
   find the minimum value of a X(mat) object
*/
R X(min)(const X(mat) *A){
    R max,min;
    X(maxmin)(A->p, A->nx*A->ny, &max, &min);
    return min;
}
/**
   find the maximum of abs of a X(mat) object
*/
R X(maxabs)(const X(mat) *A){
    R max,min;
    X(maxmin)(A->p, A->nx*A->ny, &max, &min);
    max=FABS(max);
    min=FABS(min);
    return max>min?max:min;
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
    if(in && in->nx && in->ny){
	if(!*out0) 
	    *out0=X(new)(in->nx, in->ny);
	else{
	    assert(in->nx==(*out0)->nx && in->ny == (*out0)->ny);
	}
	X(mat) *out=*out0;
	if(out->p!=in->p){
	    memcpy(out->p, in->p, in->nx*in->ny*sizeof(T));
	}
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
R X(norm2)(const X(mat)*A){
    R out=0;
    for(long i=0; i<A->nx*A->ny; i++){
	out+=(R)(A->p[i]*CONJ(A->p[i]));
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
    int colmax=10;
    int iset,i,j;
    int nset=(A->ny+colmax-1)/colmax;
    for(iset=0; iset<nset; iset++){
	int ncol=(iset+1)*colmax;
	if(ncol>A->ny) ncol=A->ny;
	printf("Cols %d to %d\n", iset, ncol-1);
	for(j=0; j<A->nx; j++){
	    for(i=iset*colmax; i<ncol; i++){
		PRINT(p[i][j]);
	    }
	    printf("\n");
	}
    }
}

/**
   scale each element of A by w
*/
void X(scale)(X(mat) *A, R w){
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
	    if(!isnan(p[i])){
		v+=p[i];
	    }
	}
    }
    return v;
}
/**
   trace
*/
T X(trace)(const X(mat)*A){
    T trace=0;
    long n=MIN(A->nx, A->ny);
    for(long i=0; i<n; i++){
	trace+=A->p[i*(1+A->nx)];
    }
    return trace;
}
/**
   compute B=bc*B+ac*A
   behavior changed on 2009-11-02. if A is NULL, don't do anything.
*/
void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac){
    if(A){
	if(!*B0){
	    *B0=X(new)(A->nx, A->ny); 
	    bc=0;/*no bother to accumulate. */
	}
	X(mat) *B=*B0;
	if(A->nx!=B->nx || A->ny != B->ny){
	    error("A is %ldx%ld, B is %ldx%ld. They should match\n",
		  A->nx, A->ny, B->nx, B->ny);
	}
	if(FABS(bc)>EPS){
	    for(int i=0; i<A->nx*A->ny; i++){
		B->p[i]=B->p[i]*bc+A->p[i]*ac;
	    }
	}else{/*just assign */
	    for(int i=0; i<A->nx*A->ny; i++){
		B->p[i]=A->p[i]*ac;
	    }
	}
    }
}
/**
   Add a scalar to matrix
*/
void X(adds)(X(mat*)A, const T ac){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]+=ac;
    }
}
/**
   compute the inner product of A and B. (inner product)
*/

T X(inn)(const X(mat)*A, const X(mat) *B){
    if(!A || !B) return 0;
    assert(A->nx==B->nx && A->ny==B->ny);
    T out=0;
    for(int i=0; i<A->nx*A->ny; i++){
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
   component-wise multiply of three matrices.
   A=A.*B.*C
*/
void X(cwm3)(X(mat) *restrict A, const X(mat) *restrict B, const X(mat) *restrict C){
    if(!B){
	X(cwm)(A,C);
    }else if (!C){
	X(cwm)(A,B);
    }else{
	assert(A && B && C);
	assert(A->nx==B->nx && A->nx==C->nx&&A->ny==B->ny && A->ny==C->ny);
	/*component-wise multiply A=A.*B */
	const long ntot=A->nx*A->ny;
	for(long i=0; i<ntot; i++){
	    A->p[i]=A->p[i]*B->p[i]*C->p[i];
	}
    }
}

/**
   Component-wise multiply each column of A with B
   A(:,i)=A(:,i).*B;
 */
void X(cwmcol)(X(mat) *restrict A, const X(mat) *restrict B){
    if (!B) return;
    assert(A->nx==B->nx && B->ny==1);
    T (*As)[A->nx]=(T(*)[A->nx])A->p;
    T *B1=B->p;
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    As[iy][ix]*=B1[ix];
	}
    }
}

/**
   component-wise multiply of columns of A with combination of B1 and B2:
   A(:,i)=A(:,i)*(B1*wt1+B2*wt2);
 */
void X(cwmcol2)(X(mat) *restrict A, 
	      const T *restrict B1, const R wt1,
	      const T *restrict B2, const R wt2){
    assert(A && A->p); 
    assert(B1);
    T (*As)[A->nx]=(T(*)[A->nx])A->p;
    if(B2){
	for(long ix=0; ix<A->nx; ix++){
	    T junk=B1[ix]*wt1+B2[ix]*wt2;
	    for(long iy=0; iy<A->ny; iy++){
		As[iy][ix]*=junk;
	    }
	}
    }else{
	for(long ix=0; ix<A->nx; ix++){
	    T junk=B1[ix]*wt1;
	    for(long iy=0; iy<A->ny; iy++){
		As[iy][ix]*=junk;
	    }
	}
    }
}

/**
   component wise multiply of 2d complex matrix A,W and 1d vector B.
   A(:,i)=A(:,i).*W(:,i).*B;
*/
void X(cwm3col)(X(mat) *restrict A,const X(mat) *restrict W,const X(mat) *restrict B){

    if(!W){
	X(cwmcol)(A,B);
    }else{
	assert(A->nx==B->nx&& A->nx==W->nx&&A->ny==W->ny&&B->ny==1);
	T (*As)[A->nx]=(T(*)[A->nx])A->p;
	T (*Ws)[W->nx]=(T(*)[W->nx])W->p;
	T *B1=B->p;

	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<A->nx; ix++){
		As[iy][ix]=As[iy][ix]*Ws[iy][ix]*B1[ix];
	    }
	}
    }

}
/**
   Component wise multiply each row of A with B.
   A(i,:)=A(i,:)*B
*/
void X(cwmrow)(X(mat) *restrict A, const X(mat) *restrict B){
    if(!A || !B) return;
    T (*As)[A->nx]=(T(*)[A->nx])A->p;
    T *B1=B->p;
    assert(A->ny==B->nx && B->ny==1);
    for(long iy=0; iy<A->ny; iy++){
	T junk=B1[iy];
	for(long ix=0; ix<A->nx; ix++){
	    As[iy][ix]*=junk;
	}
    }
}

/**
   component-wise multiply of rows of A with combination of B1 and B2:
   A(i,:)=A(i,:).*(B1*wt1+B2*wt2);
 */
void X(cwmrow2)(X(mat) *restrict A, 
	     const T *restrict B1, const R wt1,
	     const T *restrict B2, const R wt2){
    assert(A && A->p); 
    assert(B1);
    T (*As)[A->nx]=(T(*)[A->nx])A->p;
    if(B2){
	for(long iy=0; iy<A->ny; iy++){
	    T junk=B1[iy]*wt1+B2[iy]*wt2;
	    for(long ix=0; ix<A->nx; ix++){
		As[iy][ix]*=junk;
	    }
	}
    }else{
	for(long iy=0; iy<A->ny; iy++){
	    T junk=B1[iy]*wt1;
	    for(long ix=0; ix<A->nx; ix++){
		As[iy][ix]*=junk;
	    }
	}
    }
}

/**
   Component wise division B=B./A. 0/0 is replace by value;
 */
void X(cwdiv)(X(mat) *B, const X(mat) *A, T value){
    assert(A->nx==B->nx && A->ny==B->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	B->p[i]/=A->p[i];
	if(isnan(REAL(B->p[i]))) B->p[i]=value;
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
   compute the relative difference betwee two vectors.
   sqrt(||A-B||/||A||) using norm2. for debugging purpose.
*/
T X(diff)(const X(mat) *A, const X(mat) *B){
    X(mat) *C=NULL;
    X(cp)(&C,A);
    X(add)(&C,1,B,-1);
    T d=SQRT(X(norm2)(C)*2/(X(norm2)(C)+X(norm2)(B)));
    X(free)(C);
    return isnan(d)?0:d;
}

/**
   a new gray pixel map generation based on bilinear influence functions used in
   mkw.  creates slightly larger map.  add an filled circle.  cx,cy,r are in
   unit of meter, dx,dy specify the unit of pixel along x or y dimension.
*/
void X(circle)(X(mat) *A, R cx, R cy, R dx, R dy, R r, T val){
    int nres=100;
    const R res=(R)(1./nres);
    const R res1=(R)(1./nres);
    const R res2=(R)(res1*res1*4.);
    R resm=(R)((nres-1)*0.5);
    R r2=r*r;
    R r2l=(r-1.5)*(r-1.5);
    R r2u=(r+2.5)*(r+2.5);
    PMAT(A,As);
    for(int iy=0; iy<A->ny; iy++){
	R r2y=(iy*dy-cy)*(iy*dy-cy);
	for(int ix=0; ix<A->nx; ix++){
	    R r2r=(ix*dx-cx)*(ix*dx-cx)+r2y;
	    T val2=0;
	    if(r2r<r2l) {
		val2=val;
	    }else if(r2r<r2u){
		T tot=0.;
		for(int jy=0; jy<nres; jy++){
		    R iiy=iy+(jy-resm)*2*res;
		    R rr2y=(iiy*dy-cy)*(iiy*dy-cy);
		    R wty=1.-FABS(iy-iiy);
		    for(int jx=0; jx<nres; jx++){
			R iix=ix+(jx-resm)*2*res;
			R rr2r=(iix*dx-cx)*(iix*dx-cx)+rr2y;
			R wtx=1.-FABS(ix-iix);
			if(rr2r<r2){
			    tot+=res2*wty*wtx;
			}
		    }
		}
		val2=tot*val;
	    }
	    As[iy][ix]+=val2;
	}
    }
}
/**
   Similar to X(circle), by multiply instead of add.
*/
void X(circle_mul)(X(mat) *A, R cx, R cy, R dx, R dy, R r, T val){
    int nres=100;
    const R res=1./(R)(nres);
    const R res1=1./(R)(nres);
    const R res2=res1*res1*4.;
    R resm=(R)(nres-1)/2.;
    R r2=r*r;
    R r2l=(r-1.5)*(r-1.5);
    R r2u=(r+2.5)*(r+2.5);
    PMAT(A,As);
    for(int iy=0; iy<A->ny; iy++){
	R r2y=(iy*dy-cy)*(iy*dy-cy);
	for(int ix=0; ix<A->nx; ix++){
	    R r2r=(ix*dx-cx)*(ix*dx-cx)+r2y;
	    R val2=0;
	    if(r2r<r2l) {
		val2=val;
	    }else if(r2r<r2u){
		R tot=0.;
		for(int jy=0; jy<nres; jy++){
		    R iiy=iy+(jy-resm)*2*res;
		    R rr2y=(iiy*dy-cy)*(iiy*dy-cy);
		    R wty=1.-FABS(iy-iiy);
		    for(int jx=0; jx<nres; jx++){
			R iix=ix+(jx-resm)*2*res;
			R rr2r=(iix*dx-cx)*(iix*dx-cx)+rr2y;
			R wtx=1.-FABS(ix-iix);
			if(rr2r<r2)
			    tot+=res2*wty*wtx;
		    }
		}
		val2=tot*val;
	    }
	    As[iy][ix]*=val2;
	}
    }
}

/**
   Unlike X(circle), we don't use bilinear influence function, but use pure gray
   pixel instead. Also the values are black/white, no gray.
*/
void X(circle_symbolic)(X(mat) *A, R cx, R cy, R dx, R dy, R r){
    R r2=r*r;
    R r2l=(r-1.5)*(r-1.5);
    R r2u=(r+2.5)*(r+2.5);
    PMAT(A,As);
    for(int iy=0; iy<A->ny; iy++){
	R r2y=(iy*dy-cy)*(iy*dy-cy);
	for(int ix=0; ix<A->nx; ix++){
	    R r2r=(ix*dx-cx)*(ix*dx-cx)+r2y;
	    if(r2r<r2l){
	    	As[iy][ix]=1;
	    }else if(r2r>r2u){
		continue;//do not set to 0.
	    }else{
		for(R jy=-0.5; jy<1; jy++){
		    R iiy=iy+jy;
		    R rr2y=(iiy*dy-cy)*(iiy*dy-cy);
		    for(R jx=-0.5; jx<1; jx++){
			R iix=ix+jx;
			R rr2r=(iix*dx-cx)*(iix*dx-cx)+rr2y;
			if(rr2r<r2){
			    As[iy][ix]=1;
			    continue;
			}
		    }
		}
	    }
	}
    }
}


/**
   shift frequency components by n/2
*/
void X(fftshift)(X(mat) *A){
    long i;
    const long nx=A->nx;
    const long ny=A->ny;
    if((nx&1)==1){
	warning("nx=%ld is not multiple of 2\n", nx);
    }
    const long nx2=nx/2;
    const long ny2=ny/2;
    const long nx2d=nx2*sizeof(T);
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
    const long nx=A->nx;
    const long ny=A->ny;
    T *Ap=A->p;
    const long ninx=B->nx;
    const long niny=B->ny;
    if(nx>ninx || ny>niny){
	memset(Ap, 0, sizeof(T)*nx*ny);
    }
    const T * Bp=B->p;
    assert((nx&1)==0 && (ny&1)==0 && (ninx&1)==0 && (niny&1)==0);

    const int ny2=MIN(ny,niny)>>1;
    const int nx2=MIN(nx,ninx)>>1;
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
	    memcpy(B->p+iy*nx, A->p+(dy+iy)*nx+dx, sx*sizeof(T));/*3 */
	    memcpy(B->p+iy*nx+sx, A->p+(dy+iy)*nx, dx*sizeof(T));/*4 */
	}
	for(int iy=sy; iy<ny; iy++){
	    memcpy(B->p+iy*nx, A->p+(iy-sy)*nx+dx, sx*sizeof(T));/*1 */
	    memcpy(B->p+iy*nx+sx, A->p+(iy-sy)*nx, dx*sizeof(T));/*2 */
	}
    }else{
	memcpy(B->p, A->p, sizeof(T)*nx*ny);
    }
}

/**
   rotate the column vectors CCW.
   same as rotate coordinate theta CW.
   A(:,1) is x, A(:,2) is y.
*/
void X(rotvec)(X(mat) *A, const R theta){
    if(A->ny!=2) error("Wrong dimension\n");
    const R ctheta=cos(theta);
    const R stheta=sin(theta);
    PMAT(A,Ap);
    for(int i=0; i<A->nx; i++){
	T tmp=Ap[0][i]*ctheta-Ap[1][i]*stheta;
	Ap[1][i]=Ap[0][i]*stheta+Ap[1][i]*ctheta;
	Ap[0][i]=tmp;
    }
}

/**
   rotate the row vectors CCW.
   same as rotate coordinate theta CW.
   A(:,1) is x, A(:,2) is y.
*/
void X(rotvect)(X(mat) *A, const R theta){
    if(A->nx!=2) error("Wrong dimension\n");
    const R ctheta=cos(theta);
    const R stheta=sin(theta);
    PMAT(A,Ap);
    for(int i=0; i<A->ny; i++){
	T tmp=Ap[i][0]*ctheta-Ap[i][1]*stheta;
	Ap[i][1]=Ap[i][0]*stheta+Ap[i][1]*ctheta;
	Ap[i][0]=tmp;
    }
}

/**
   rotate a 2x2 covariance matrix A by theta CCW
   (coordinate rotate -theta CCW) or from ra to xy
   coordinate.  R*A*R';
*/
void X(rotvecnn)(X(mat) **B0, const X(mat) *A, R theta){
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
    /*first apply left R */
    tmp[0][0]=ctheta*Ap[0][0]-stheta*Ap[0][1];
    tmp[1][0]=ctheta*Ap[1][0]-stheta*Ap[1][1];
    tmp[0][1]=stheta*Ap[0][0]+ctheta*Ap[0][1];
    tmp[1][1]=stheta*Ap[1][0]+ctheta*Ap[1][1];
    /*then apply right R' */
    
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
    /*calculate y=A*x for 3. */
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
void X(cog)(R *grad,const X(mat) *im,R offsetx,
	    R offsety, R thres, R bkgrnd){
    R sum=0,sumx=0,sumy=0;
    R iI;
    PMAT(im,pim);
    for(int iy=0; iy<im->ny; iy++){
	for(int ix=0; ix<im->nx; ix++){
	    iI=REAL(pim[iy][ix])-bkgrnd;
	    if(iI>thres){
		sum+=iI;
		sumx+=iI*ix;
		sumy+=iI*iy;
	    }
	}
    }
    if(FABS(sum)>0){
	grad[0]=sumx/sum-((R)(im->nx-1)*0.5+offsetx);
	grad[1]=sumy/sum-((R)(im->ny-1)*0.5+offsety);
    }else{
	grad[0]=0;
	grad[1]=0;
    }
}

/**
   Shift the image in A to center on physical
   center+[offsetx,offsety] using cog and fft.
*/
#ifndef USE_SINGLE
void X(shift2center)(X(mat) *A, R offsetx, R offsety){
    R grad[2];
    R Amax=X(max)(A);
    X(cog)(grad,A,offsetx,offsety,Amax*0.1,Amax*0.2);
    if(FABS(grad[0])>0.1 || FABS(grad[1])>0.1){
	/*info("Before shift, residual grad is %g %g\n",grad[0],grad[1]); */
	XC(mat) *B=XC(new)(A->nx,A->ny);
	XC(fft2plan)(B,-1);
	XC(fft2plan)(B,1);
#ifdef USE_COMPLEX
	XC(cp)(&B,A);
#else
	XC(cpd)(&B,A);
#endif
	R scale=1./(A->nx*A->ny);
	XC(fftshift)(B);
	XC(fft2)(B,-1);
	XC(tilt)(B,-grad[0],-grad[1],0);
	XC(fft2)(B,1);
	XC(fftshift)(B);
	XC(scale)(B,scale);
#ifdef USE_COMPLEX
	XC(cp)(&A,B);
#else
	XC(real2d)(&A,0,B,1);
#endif
	X(cog)(grad,A,offsetx,offsety,Amax*0.1,Amax*0.2);
	/*info("After shift, residual grad is %g %g\n",grad[0],grad[1]); */
	XC(free)(B);
    }
}
#endif
/**
   Limit numbers in A to within [min, max]. used for DM clipping.
*/
int X(clip)(X(mat) *A, R min, R max){
    if(!A) return 0;
    if(!isfinite(min)==-1 && !isfinite(max)==1) return 0;
    if(max<=min){
	error("upper light should be larger than lower limit\n");
    }
    T *restrict Ap=A->p;
    int nclip=0;
    for(long i=0; i<A->nx *A->ny; i++){
	R Ar=REAL(Ap[i]);
	if(Ar>max) {
	    Ap[i]=max;
	    nclip++;
	}else if(Ar<min) {
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
void X(gramschmidt)(X(mat) *Mod, R *amp){
    const int nmod=Mod->ny;
    const long nx=Mod->nx;
    R wtsum=(R)nx;
    if(amp){
#ifdef USE_SINGLE
	wtsum=fltsum(amp, nx);
#else
	wtsum=dblsum(amp, nx);
#endif
    }
    int nonvalid[nmod];
    memset(nonvalid, 0, sizeof(int)*nmod);
    PMAT(Mod,pMod);
    for(int imod=0; imod<nmod; imod++){
	if(imod>0){/*orthogonalize */
	    R cross;
	    /*compute dot product. */
	    for(int jmod=0; jmod<imod; jmod++){
		if(nonvalid[jmod]) continue;
		cross=-DOT(pMod[imod],pMod[jmod],amp,nx)/wtsum;
		for(long ix=0; ix<nx; ix++){
		    pMod[imod][ix]+=cross*pMod[jmod][ix];
		}
	    }
	}
	
	R norm=SQRT(DOT(pMod[imod],pMod[imod],amp,nx)/wtsum);
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
void X(muldiag)(X(mat) *A, const X(mat) *s){
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
   A=B*A*B, where diag(B)=s
*/
void X(muldiag2)(X(mat) *A, const X(mat) *s){
    assert(A->ny==s->nx && s->ny==1);
    PMAT(A,pA);
    const T *ps=s->p;
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    pA[iy][ix]*=ps[iy]*ps[ix];
	}
    }
}
/**
   Raise all elements to power power
*/
void X(cwpow)(X(mat)*A, R power){
    if(!A) return;
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=POW(A->p[i],power);
    }
}

/**
   compute exponential of all elements after scaling by alpha
*/
void X(cwexp)(X(mat)*A, R alpha){
    if(!A) return;
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=EXP(A->p[i]*alpha);
    }
}

/**
   Raise all elements above thres*maxabs(A) to pow power. Set others to zero.
*/
void X(cwpow_thres)(X(mat) *A, R power, R thres){
    thres*=X(maxabs)(A);
    for(long i=0; i<A->nx*A->ny; i++){
	if(ABS(A->p[i])>thres){
	    A->p[i]=POW(A->p[i], power);
	}else{
	    A->p[i]=0;
	}
    }
}
/**
   Compute polynomial functional value inplace
 */
void X(polyval)(X(mat) *A, XR(mat)*p){
    if(p->nx==1 && p->ny>1){
	p->nx=p->ny;
	p->ny=1;
    }
    int np=p->nx;
    for(long i=0; i<A->nx*A->ny; i++){
	T tmp=0;
	for(long ip=0; ip<np; ip++){
	    tmp+=p->p[ip]*POW(A->p[i], np-ip-1);
	}
	A->p[i]=tmp;
    }
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
	assert(x->nx==y->nx && x->ny==A->m);
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
X(mat)* X(logspace)(R emin, R emax, long n){
    X(mat)* out=X(new)(n,1);
    R esep=(emax-emin)/(n-1);
    for(long i=0; i<n; i++){
	R ex=emin+esep*i;
	out->p[i]=pow(10, ex);
    }
    return out;
}

/**
   Create linearly spaced vector.
*/
X(mat)* X(linspace)(R min, R dx, long n){
    X(mat)* out=X(new)(n,1);
    for(long i=0; i<n; i++){
	out->p[i]=min+dx*i;
    }
    return out;
}
/**
   Check whether xin is linearly spaced
*/
static int X(islinear)(const X(mat)*xin){
    long nmax=xin->nx;
    long nmax1=nmax-1;
    R xminl=(xin->p[0]);
    R xmaxl=(xin->p[nmax-1]);
    R xsep=(xmaxl-xminl)/(R)(nmax1);
    if(FABS(xsep+xminl-xin->p[1])>xsep*1.e-3){
	return 0;
    }else{
	return 1;
    }
}
/**
   Check whether xin is logrithmically spaced
 */
static int X(islog)(const X(mat)*xin){
    long nmax=xin->nx;
    long nmax1=nmax-1;
    R xminl=log10(xin->p[0]);
    R x1=log10(xin->p[1]);
    R xmaxl=log10(xin->p[nmax1]);
    R xsep=(xmaxl-xminl)/(R)(nmax1);
    if(!isfinite(xsep) || FABS(xsep+xminl-x1)>xsep*1.e-3){
	return 0;
    }else{
	return 1;
    }
}
/**
   Interpolate using linear interp. xin is the coordinate of yin. xnew is the
   coordinate of the output.
*/
X(mat)* X(interp1linear)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew){
    if(!X(islinear)(xin)){
	error("xin is not linearly spaced\n");
    }
    if(xin->ny!=1 || xnew->ny!=1){
	error("Either xin or xnew is in wrong format\n");
    }
    long nmax=xin->nx;
    long nmax1=nmax-1;
    R xminl=(xin->p[0]);
    R xmaxl=(xin->p[nmax-1]);
    R xsep=(xmaxl-xminl)/(R)(nmax1);
    R xsep1=1./xsep;
    X(mat) *ynew=X(new)(xnew->nx, xnew->ny);
    PMAT(yin, pyin);
    PMAT(ynew, pynew);
    for(long iy=0; iy<ynew->ny; iy++){
	for(long ix=0; ix<ynew->nx; ix++){
	    R xx=((xnew->p[ix])-xminl)*xsep1;
	    long xxm=ifloor(xx);
	    if(xxm<0){
		pynew[iy][ix]=pyin[iy][0];
	    }else if(xxm>=nmax1){
		pynew[iy][ix]=pyin[iy][nmax1];
	    }else{
		R xxw=xx-xxm;
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
X(mat)* X(interp1log)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew){
    if(!X(islog)(xin)){
	error("xin is not logrithmically spaced\n");
    }
    if(xin->ny!=1 || xnew->ny!=1){
	error("Either xin or xnew is in wrong format\n");
    }
    long nmax=xin->nx;
    long nmax1=nmax-1;
    R xminl=log10(xin->p[0]);
    R xmaxl=log10(xin->p[nmax-1]);
    R xsep=(xmaxl-xminl)/(R)(nmax1);
    R xsep1=1./xsep;
    X(mat) *ynew=X(new)(xnew->nx, xnew->ny);
    PMAT(yin, pyin);
    PMAT(ynew, pynew);
    for(long iy=0; iy<ynew->ny; iy++){
	for(long ix=0; ix<ynew->nx; ix++){
	    R xx=(log10(xnew->p[ix])-xminl)*xsep1;
	    long xxm=ifloor(xx);
	    if(xxm<0){
		pynew[iy][ix]=pyin[iy][0];
	    }else if(xxm>=nmax1){
		pynew[iy][ix]=pyin[iy][nmax1];
	    }else{
		R xxw=xx-xxm;
		pynew[iy][ix]=xxw*pyin[iy][xxm+1]+(1.-xxw)*pyin[iy][xxm];
	    }
	}
    }
    return ynew;
}
#ifndef USE_COMPLEX
X(mat)* X(interp1)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew){
    int free_xy=0;
    X(mat)*ynew=NULL;
    if(!yin){
	yin=X(new_ref)(xin->nx, 1, xin->p+xin->nx);
	xin=X(new_ref)(xin->nx, 1, xin->p);
	free_xy=1;
    }
    if(X(islinear)(xin)){
	ynew=X(interp1linear)(xin, yin, xnew);
    }else if(X(islog)(xin)){
        ynew=X(interp1log)(xin, yin, xnew);
    }else{//arbitrary spacing
	if(xin->ny!=1 || xnew->ny!=1){
	    error("Either xin or xnew is in wrong format\n");
	}
	ynew=X(new)(xnew->nx, xnew->ny); 
	PMAT(yin, pyin);
	PMAT(ynew, pynew);
	int curpos=0;
	for(long ix=0; ix<ynew->nx; ix++){
	    for(curpos=0; curpos<xin->nx-2; curpos++){
		if(xnew->p[ix]>xin->p[curpos] && xnew->p[ix]<xin->p[curpos+1]){
		    break;
		}
	    }
	    
	    R xx=((xnew->p[ix])-xin->p[curpos])/(xin->p[curpos+1]-xin->p[curpos]);
	    for(long iy=0; iy<ynew->ny; iy++){
		pynew[iy][ix]=xx*pyin[iy][curpos+1]+(1.-xx)*pyin[iy][curpos];
	    }
	}
    }
    if(free_xy){
	X(mat) *tmp=(X(mat)*)xin; X(free)(tmp);
	tmp=(X(mat)*)yin; X(free)(tmp);
    }
    return ynew;
}

#endif
#ifndef USE_COMPLEX
/**
   embed a ninx*niny matrix in into A with optional rotation by -theta CCW
   (coordinate rotate theta CCW) around the fft center. Used to rotate the PSF
   from x-y to radial-azimuthal coordinate in radial format CCD. A may be bigger or smaller than B.
   \todo{
   merge this definition with cembed in cmat.c
   }
*/

/*blend B into center of A with width of overlap. The center
  (size is B->nx-overlap, B->ny-overlap) of A is replaced by
  center of B . The overlapping area is blended*/
void X(blend)(X(mat) *restrict A, X(mat) *restrict B, int overlap){
    const long ninx=B->nx;
    const long niny=B->ny;
    const long noutx=A->nx;
    const long nouty=A->ny;
    const long skipx=(noutx-ninx)/2;
    const long skipy=(nouty-niny)/2;
    long ixstart=0, ixlen=ninx;
    long iystart=0, iylen=niny;

    if(skipx<0){
	ixstart=-skipx;
	ixlen+=2*skipx;
    }
    if(skipy<0){
	iystart=-skipy;
	iylen+=2*skipy;
    }
    PMAT(A, pA);
    PMAT(B, pB);
    R wty, wtx;
    for(long iy=0; iy<iylen; iy++){
	T *outi=&pA[iystart+skipy+iy][ixstart+skipx];
	T *ini =&pB[iystart+iy][ixstart];
	if(iy<overlap){
	    wty=(R)iy/(R)(overlap-1);
	}else if(iylen-iy-1<overlap){
	    wty=(R)(iylen-iy-1)/(R)(overlap-1);
	}else{
	    wty=1;
	}
	for(long ix=0; ix<ixlen; ix++){
	    if(ix<overlap){
		wtx=(R)ix/(R)(overlap-1);
	    }else if(ixlen-ix-1<overlap){
		wtx=(R)(ixlen-ix-1)/(R)(overlap-1);
	    }else{
		wtx=1;
	    }
	    outi[ix]=outi[ix]*(1-wtx*wty)+ini[ix]*wtx*wty;
	}
    }
}
/**
   For each entry in A, call repeatly to collect its histogram, centered at
   center, spaced by spacing, for n bins in total. center if at bin n/2.  */
void X(histfill)(X(mat) **out, const X(mat)* A,
		 R center, R spacing, int n){
    if(!A || !A->p) return;
    int nn=A->nx*A->ny;
    if(!*out){
	*out=X(new)(n,nn);
    }
    PMAT(*out,Op);
    const T *restrict Ap=A->p;
    const R spacingi=1./spacing;
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

   see www.paulinternet.nl/?page=bicubicx */
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
    R thres=ABS(xsep)*1.e-5;
  
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
	pc[ix][0]=-0.5*ypriv+1.5*py[ix]-1.5*py[ix+1]+0.5*ynext;/*a */
	pc[ix][1]=     ypriv-2.5*py[ix]+2.0*py[ix+1]-0.5*ynext;/*b */
	pc[ix][2]=-0.5*ypriv           +0.5*py[ix+1];/*c */
	pc[ix][3]=               py[ix] ;/*d */
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
    T xsep1=(T)(nx-1)/(x->p[nx-1]-xmin);
    X(mat) *out=X(new)(xnew->nx, xnew->ny);
    for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
	R xn=REAL((xnew->p[ix]-xmin)*xsep1);
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
   Do a component wise log10 on each element of A.
*/
void X(cwlog10)(X(mat) *A){
    R ratio=1./log(10);
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=LOG(A->p[i])*ratio;
    }
}
/**
   Do a component wise log10 on each element of A.
*/
void X(cwlog)(X(mat) *A){
    for(long i=0; i<A->nx*A->ny; i++){
	A->p[i]=LOG(A->p[i]);
    }
}

void X(embed)(X(mat) *restrict A, const X(mat) *restrict B, const R theta){
  
    const long ninx=B->nx;
    const long niny=B->ny;
    const long noutx=A->nx;
    const long nouty=A->ny;
    memset(A->p, 0, sizeof(T)*noutx*nouty);
    if(FABS(theta)<1.e-10){/*no rotation. */
	const long skipx=(noutx-ninx-1)/2;//-1 to handle odd case
	const long skipy=(nouty-niny-1)/2;
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
	const R ctheta=cos(theta);
	const R stheta=sin(theta);
	R x2,y2;
	R x,y;
	long ninx2=ninx/2;
	long noutx2=noutx/2;
	long niny2=niny/2;
	long nouty2=nouty/2;
	long ix2, iy2;
	for(long iy=0; iy<nouty; iy++){ 
	    y=(R)(iy-nouty2); 
	    for(long ix=0; ix<noutx; ix++){ 
		x=(R)(ix-noutx2); 
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
   Calculate number of pixels having values larger than or equal to half of
maximum. Useful to compute fwhm. */
long X(fwhm)(X(mat) *A){
    if(!A) return 0;
    R hm=0.5*X(max)(A);
    long fwhm=0;
    for(long ix=0; ix<A->nx*A->ny; ix++){
	if(ABS(A->p[ix])>=hm){
	    fwhm++;
	}
    }
    return fwhm;
}
#ifndef USE_COMPLEX
static int sort_ascend(const T*A, const T*B){
    if ((*A)>(*B)) return 1; 
    else return -1;
}
static int sort_descend(const T*A, const T*B){
    if ((*A)>(*B)) return -1; 
    else return 1;
}
/*
  Sort all columns of A, in ascending order if ascend is non zero, otherwise in descending order.
 */
void X(sort)(X(mat) *A, int ascend){
    for(int i=0; i<A->ny; i++){
	if(ascend){
	    qsort(A->p+i*A->nx, A->nx, sizeof(R), 
		  (int(*)(const void*,const void*))sort_ascend);
	}else{
	    qsort(A->p+i*A->nx, A->nx, sizeof(R), 
		  (int(*)(const void*,const void*))sort_descend);
	}
    }
}

typedef struct{
    X(mat) *enc; /**<Output*/
    X(mat) *dvec;/**<Radius wanted*/
    X(mat) *phat; /**<processed image.*/
    int type;  
}ENC_T;

static void X(enc_thread)(thread_t *pdata){
    ENC_T *data=pdata->data;
    const X(mat) *dvec=data->dvec;
    X(mat) *enc=data->enc;
    PMAT(data->phat, ppsf);
    int type=data->type;
    const R *restrict dr=dvec->p;
    const long ncomp2=data->phat->nx;
    const long ncomp=ncomp2/2;
    const R dk=1./ncomp2;
    const R pi2=2*M_PI;
    if(type==0){
	X(mat) *ksinc=X(new)(dvec->nx, ncomp2);
	PMAT(ksinc, pks);
	/*Cache the data. */
	for(long iy=0; iy<ncomp2; iy++){
	    R ky=(iy<ncomp?iy:iy-ncomp2)*dk;
	    for(long ir=pdata->start; ir<pdata->end; ir++){
		pks[iy][ir]=sinc(ky*dr[ir])*dr[ir];
	    }
	}
	for(long iy=0; iy<ncomp2; iy++){
	    for(long ix=0; ix<ncomp2; ix++){
		for(long ir=pdata->start; ir<pdata->end; ir++){
		    R s=pks[iy][ir]*pks[ix][ir];
		    enc->p[ir]+=s*ppsf[iy][ix];
		}
	    }
	}
    }else{
	for(long iy=0; iy<ncomp2; iy++){
	    R ky=(iy<ncomp?iy:iy-ncomp2)*dk;
	    for(long ix=0; ix<ncomp2; ix++){
		R kx=(ix<ncomp?ix:ix-ncomp2)*dk;
		switch(type){
		case -1: {/*azimuthal average. dr is radius */
		    R k=sqrt(kx*kx+ky*ky);
		    for(long ir=pdata->start; ir<pdata->end; ir++){
			R s=j0(k*pi2*dr[ir]);
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 0:
		    break;
		case 1: {/*Encircled energy. dr is diameter */
		    R k=sqrt(kx*kx+ky*ky);
		    for(long ir=pdata->start; ir<pdata->end; ir++){
			const R r=dr[ir]*0.5;
			const R tmp=k*pi2*r;
			R s=j1(tmp)*r/k;
			if(!ix && !iy) s=pi2*r*r;/*special case. */
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 2:/*Enstripped energe in a slit. */
		    error("To implement: Do FFT only along 1-d\n");
		    break;
		default:
		    error("Not implemented\n");
		}
	    }
	}
    }
}
/**
   Compute the enclosed energy or azimuthal average of a.
*/
X(mat) *X(enc)(X(mat) *psf, /**<The input array*/
	       X(mat) *dvec,/**<The diameter for enclosed energy, or radius for azimuthal average*/
	       int type,  /**<The type. -1: azimuthal average, 0: within a square, 1: within a circle, 2: within a slit*/
	       int nthread
    ){
    R rmax=ceil(X(max)(dvec))+1;
    long ncomp;
    ncomp=nextfftsize(rmax*2);//avoid wrapping
    long ncomp_max=psf->nx>psf->ny?psf->nx:psf->ny;
    X(mat) *psfc;
    if(ncomp_max > ncomp){
	psfc=X(new)(ncomp, ncomp);
	X(embed)(psfc, psf, 0);
    }else{
	ncomp=ncomp_max;
	psfc=X(ref)(psf);
    }
    long ncomp2=ncomp*2;
    XC(mat) *psf2=XC(new)(ncomp2, ncomp2);
    XC(fft2plan)(psf2, -1);
    XC(embedd)(psf2, psfc, 0);
    X(free)(psfc);
    XC(fftshift)(psf2);
    XC(fft2)(psf2, -1);
    X(mat) *phat=NULL;
    XC(real2d)(&phat, 0, psf2, 1);
    X(scale)(phat, pow((R)ncomp2,-2));
    XC(free)(psf2);
    X(mat) *enc=X(new)(dvec->nx, 1);    
    ENC_T data={enc, dvec, phat, type};
    thread_t info[nthread];
    thread_prep(info, 0, dvec->nx, nthread, X(enc_thread), &data);
    CALL_THREAD(info, 0);
    X(free)(phat);
    return enc;
}

#endif

X(cell) *X(cellnew)(long nx, long ny){
    return (X(cell*)) cellnew(nx, ny);
}
/**
   create an new X(cell) similar to A in shape.
   When a cell is empty, it is created with a (0,0) array and cannot be overriden.
*/
X(cell) *X(cellnew2)(const X(cell) *A){
    X(cell) *out=X(cellnew)(A->nx, A->ny);
    long tot=0;
    for(long i=0; i<A->nx*A->ny; i++){
	if(!isempty(A->p[i])){
	    tot+=A->p[i]->nx*A->p[i]->ny;
	}
    }
    out->m=X(new)(tot,1);
    tot=0;
    for(int i=0; i<A->nx*A->ny; i++){
	if(!isempty(A->p[i])){
	    out->p[i]=X(new_ref)(A->p[i]->nx, A->p[i]->ny, out->m->p+tot);
	    tot+=A->p[i]->nx*A->p[i]->ny;
	}else{
	    out->p[i]=X(new)(0,0);//place holder to avoid been overriden.
	}
    }
    return out;
}
/**
   Create an new X(cell) with X(mat) specified. Each block is stored continuously in memory.
 */
/*X(cell) *X(cellnew3s)(long nx, long ny, long *nnx, long *nny){
    X(cell) *out=X(cellnew)(nx,ny);
    for(long i=0; i<nx*ny; i++){
	out->p[i]=X(new)(nnx[i], (nny?nny[i]:1));
    }
    return out;
    }*/
/**
   Create an new X(cell) with X(mat) specified. Each block is stored continuously in memory.
 */
X(cell) *X(cellnew3)(long nx, long ny, long *nnx, long *nny){
    X(cell) *out=X(cellnew)(nx,ny);
    long tot=0;
    for(long i=0; i<nx*ny; i++){
	tot+=nnx[i]*(nny?nny[i]:1);
    }
    if(!tot) return NULL;
    out->m=X(new)(tot,1);
    tot=0;
    for(long i=0; i<nx*ny; i++){
	out->p[i]=X(new_ref)(nnx[i], (nny?nny[i]:1), out->m->p+tot);
	tot+=nnx[i]*(nny?nny[i]:1);
    }
    return out;
}

/**
   Free a X(cell) object.
*/
void X(cellfree_do)(X(cell) *dc){
    if(!dc) return;
    if(dc->p){
	if(dc->header){
	    R count=search_header_num(dc->header, "count");
	    if(!isnan(count) && count>0){
		info("count=%g, scaling the data\n", count);
		X(cellscale)(dc, 1./count);
	    }
	}
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    X(free)(dc->p[ix]);
	}
	if(dc->mmap){
	    mmap_unref(dc->mmap);
	}else{
	    free(dc->header);
	}
	free(dc->p);dc->p=0;
    }
    if(dc->m) X(free)(dc->m);
    if(dc->fft) X(fft_free_plan)(dc->fft);
    free(dc);
}
/**
   check the size of matrix if exist. Otherwise create it
*/
void X(cellinit)(X(cell)**A, long nx, long ny){
    if(*A){
	assert((*A)->nx == nx && (*A)->ny==ny);
    }else{
	*A=X(cellnew)(nx, ny);
    }
}
/**
   creat a X(cell) reference an existing X(cell) by referencing the
   elements.
*/
X(cell) *X(cellref)(const X(cell) *in){
    if(!in) return NULL;
    X(cell) *out=X(cellnew)(in->nx, in->ny);
    if(in->m){
	out->m=X(ref)(in->m);
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=X(new_ref)(in->p[i]->nx, in->p[i]->ny, in->p[i]->p);
	}
    }else{
	for(int i=0; i<in->nx*in->ny; i++){
	    out->p[i]=X(ref)(in->p[i]);
	}
    }
    return out;
}

/**
   duplicate a X(cell) object.
*/
X(cell) *X(celldup)(const X(cell) *in){
    X(cell) *out=NULL;
    X(cellcp)(&out, in);
    return out;
}

/**
   copy the values from one X(cell) to another.
*/
void X(cellcp)(X(cell)** out0, const X(cell) *in){
    if(in){
	if(!*out0){
	    *out0=X(cellnew2)(in);
	}else{
	    assert((*out0)->nx==in->nx && (*out0)->ny==in->ny);
	}
	X(cell)* out=*out0;
	for(int i=0; i<in->nx*in->ny; i++){
	    X(cp)(&out->p[i], in->p[i]);
	}
    }else{/*if input is empty, zero the output. do not free. */
	X(cellzero)(*out0);
    }
}

/**
   setting all elements of a X(cell) to zero.
*/
void X(cellzero)(X(cell) *dc){
    if(dc){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    X(zero)(dc->p[ix]);
	}
    }
}

/**
   setting all elements of a X(cell) to alpha.
*/
void X(cellset)(X(cell)*dc, T alpha){
    if(dc){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    X(set)(dc->p[ix],alpha);
	}
    }
}

/**
   transpose a X(cell) object
*/
X(cell) *X(celltrans)(const X(cell) *A){
    if(!A) return NULL;
    X(cell) *B=X(cellnew)(A->ny, A->nx);
    X(mat)* (*Bp)[B->nx]=(void*)B->p;
    X(mat)* (*Ap)[A->nx]=(void*)A->p;
    for(int ix=0; ix<A->nx; ix++){
	for(int iy=0; iy<A->ny; iy++){
	    Bp[ix][iy]=X(trans)(Ap[iy][ix]);
	}
    }
    return B;
}

/**
   compute norm2.
*/
R X(cellnorm2)(const X(cell) *A){
    R out=0;
    for(int i=0; i<A->nx*A->ny; i++){
	out+=X(norm2)(A->p[i]);
    }
    return out;
}

/**
   scale each element of A.
*/
void X(cellscale)(X(cell) *A, R w){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	X(scale)(A->p[i],w);
    }
}
/**
   Obtain the dimensions.
*/
static void X(celldim)(const X(cell) *A, long *nx, long *ny,
		       long **nxs, long **nys){
    X(mat) *(*Ap)[A->nx] = (X(mat) *(*)[A->nx])A->p;
    *nxs=calloc(A->nx, sizeof(long));
    *nys=calloc(A->ny, sizeof(long));
    *nx=0;
    *ny=0;
    for(long ix=0; ix<A->nx; ix++){
	for(long iy=0; iy<A->ny; iy++){
	    if(!isempty(Ap[iy][ix])){
		*nx+=Ap[iy][ix]->nx;
		(*nxs)[ix]=Ap[iy][ix]->nx;
		break;
	    }
	}
    }
    for(long iy=0; iy<A->ny; iy++){
	for(long ix=0; ix<A->nx; ix++){
	    if(!isempty(Ap[iy][ix])){
		*ny+=Ap[iy][ix]->ny;
		(*nys)[iy]=Ap[iy][ix]->ny;
		break;
	    }
	}
    }
}

/**
   reduce nx*ny cell matrix to 1*ny if dim=1 and nx*1 if dim=2
*/
X(cell) *X(cellreduce)(const X(cell)*A, int dim){
    if(!A) return NULL;
    X(cell)* out=NULL;
    long nx, ny, *nxs, *nys;
    X(celldim)(A, &nx, &ny, &nxs, &nys);
    if(nx==0 || ny==0) return NULL;
    PCELL(A,pA);
    if(dim==1){
	out=X(cellnew)(1, A->ny);
	for(long iy=0; iy<A->ny; iy++){
	    if(nys[iy]==0) continue;
	    out->p[iy]=X(new)(nx,nys[iy]);
	    for(long icol=0; icol<nys[iy]; icol++){
		long kr=0;
		for(long ix=0; ix<A->nx; ix++){
		    if(!isempty(pA[iy][ix])){
			memcpy(out->p[iy]->p+icol*nx+kr,
			       pA[iy][ix]->p+icol*nxs[ix],
			       nxs[ix]*sizeof(T));
		    }			
		    kr+=nxs[ix];
		}
	    }
	}
    }else if(dim==2){
	out=X(cellnew)(A->nx,1);
	for(long ix=0; ix<A->nx; ix++){
	    if(nxs[ix]==0) continue;
	    out->p[ix]=X(new)(nxs[ix],ny);
	    long kr=0;
	    for(long iy=0; iy<A->ny; iy++){
		if(!isempty(pA[iy][ix])){
		    memcpy(out->p[ix]->p+kr*nxs[ix],
			   pA[iy][ix]->p,
			   nxs[ix]*nys[iy]*sizeof(T));
		}
		kr+=nys[iy];
	    }
	}
    }else{
	error("Invalid dim=%d\n",dim);
    }
    free(nxs);
    free(nys);
    return out;
}

/**
   concatenate two cell matrices along dimenstion 'dim'
*/
X(cell) *X(cellcat)(const X(cell) *A, const X(cell) *B, int dim){
    if(!A){
	if(!B){
	    return NULL;
	}else{
	    return X(celldup)(B);
	}
    }else if(!B){
	return X(celldup)(A);
    }

    X(cell) *out=NULL;
    PCELL(A,pA);
    PCELL(B,pB);

    if(dim==1){
	/*along x. */
	if(A->ny!=B->ny){
	    error("Mismatch: A is (%ld, %ld), B is (%ld, %ld)\n",
		  A->nx, A->ny, B->nx, B->ny);
	}
	out=X(cellnew)(A->nx+B->nx, A->ny);
	PCELL(out,pout);
	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<A->nx; ix++){
		pout[iy][ix]=X(dup)(pA[iy][ix]);
	    }
	    for(long ix=0; ix<B->nx; ix++){
		pout[iy][ix+A->nx]=X(dup)(pB[iy][ix]);
	    }
	}
    }else if(dim==2){
	/*along y. */
	if(A->nx!=B->nx){
	    error("Mismatch. A is (%ld, %ld), B is (%ld, %ld)\n", 
		  A->nx, A->ny, B->nx, B->ny);
	}
	out=X(cellnew)(A->nx, A->ny+B->ny);
	PCELL(out,pout);
	for(long iy=0; iy<A->ny; iy++){
	    for(long ix=0; ix<A->nx; ix++){
		pout[iy][ix]=X(dup)(pA[iy][ix]);
	    }
	}
	for(long iy=0; iy<B->ny; iy++){
	    for(long ix=0; ix<B->nx; ix++){
		pout[iy+A->ny][ix]=X(dup)(pB[iy][ix]);
	    }
	}
    }else{
	error("Invalid dim\n");
    }
    return out;
}

/**
   concatenate coresponding elements of each X(cell). They must
   have the same shape.
*/
X(cell) *X(cellcat_each)(const X(cell) *A, const X(cell) *B, int dim){
    if(!A){
	if(!B){
	    return NULL;
	}else{
	    return X(celldup)(B);
	}
    }else if(!B){
	return X(celldup)(A);
    }
    if(A->nx!=B->nx || A->ny!=B->ny){
	error("Mismatch: (%ld %ld), (%ld %ld)\n",A->nx, A->ny, B->nx, B->ny);
    }
    X(cell) *out=X(cellnew)(A->nx, A->ny);
    for(long ix=0; ix<A->nx*A->ny; ix++){
	out->p[ix]=X(cat)(A->p[ix], B->p[ix], dim);
    }
    return out;
}

/**
   drop empty rows or columns. (size of *A0 is changed.
*/
void X(celldropempty)(X(cell) **A0, int dim){
    X(cell) *A=*A0;
    if(!A) return;
    PCELL(A,pA);
    if(dim==1){
	/*drop rows */
	int keep[A->nx];
	int ndrop=0;
	for(int ix=0; ix<A->nx; ix++){
	    keep[ix]=0;
	    for(int iy=0; iy<A->ny; iy++){
		if(!isempty(pA[iy][ix])){
		    keep[ix]=1;
		    break;
		}
	    }
	    if(keep[ix]==0){
		ndrop++;
	    }
	}
	if(ndrop){
	    if(ndrop==A->nx){
		X(cellfree)(A);
		*A0=NULL;
	    }else{
		X(cell) *B=calloc(1, sizeof(X(cell)*));
		B->p=calloc((A->nx-ndrop)*A->ny, sizeof(X(mat)*));
		PCELL(B,pB);
		int count=0;
		for(int ix=0; ix<A->nx; ix++){
		    if(keep[ix]){
			if(count!=ix){
			    for(int iy=0; iy<A->ny; iy++){
				pB[iy][count]=pA[iy][ix];
			    }
			}
			count++;
		    }else{
			warning("row %d dropped\n", ix);
		    }
		}
		free(A->p); free(A);
		A=B;
	    }
	}
    }else if(dim==2){
	/*drop cols */
	int count=0;
	for(int iy=0; iy<A->ny; iy++){
	    int keep=0;
	    for(int ix=0; ix<A->nx; ix++){
		if(!isempty(pA[iy][ix])){
		    keep=1;
		    break;
		}
	    }
	    if(keep){
		for(int ix=0; ix<A->nx; ix++){
		    pA[count][ix]=pA[iy][ix];
		}
		count++;
	    }else{
		/*warning("Col %d dropped\n", iy); */
	    }
	}
	A->ny=count;
	if(count==0){
	    X(cellfree)(A);
	    *A0=NULL;
	}else{
	    A->p=realloc(A->p,sizeof(X(mat)*)*A->ny*A->nx);
	}
    }else{
	error("Invalid dim: %d\n",dim);
    }
 
}

/**
   add one to another.  B=B*bc+A*ac
*/
void X(celladd)(X(cell) **B0, R bc, const X(cell) *A,const R ac){
    if(A){
	X(cell) *B=*B0;
	if(!B){
	    B=*B0=X(cellnew2)(A); bc=0;
	}else{
	    assert(A->nx==B->nx && A->ny == B->ny);
	}
	for(int i=0; i<A->nx*A->ny; i++){
	    X(add)(&B->p[i], bc, A->p[i], ac);
	}
    }
}

/**
   Compute the inner produce of two dcell.
*/
T X(cellinn)(const X(cell)*A, const X(cell)*B){
    if(!A || !B) return 0; 
    if(A->nx!=B->nx || A->ny!=1 || B->ny!=1) error("mismatch\n");
    T out=0;
    for(int i=0; i<A->nx; i++){
	out+=X(inn)(A->p[i], B->p[i]);
    }
    return out;
}

/**
   Component wise multiply of two dcell
   B=A.*B*alpha
*/
void X(cellcwm)(X(cell) *B, const X(cell) *A){
    if(A->nx!=B->nx || A->ny !=B->ny) error("mismatch\n");
    for(int i=0; i<A->nx*A->ny; i++){
	X(cwm)(B->p[i], A->p[i]);
    }
}

/**
   Convert a block matrix to a matrix.
*/
X(mat) *X(cell2m)(const X(cell) *A){
    if(A->nx*A->ny==1){
	return X(dup)(A->p[0]);
    }
    PCELL(A, Ap);
    //X(mat) *(*Ap)[A->nx] = (X(mat) *(*)[A->nx])A->p;
    long nx,ny,*nxs,*nys;
    X(celldim)(A,&nx,&ny,&nxs,&nys);
    X(mat) *out=X(new)(nx,ny);
    long jcol=0;
    for(long iy=0; iy<A->ny; iy++){
	for(long icol=0; icol<nys[iy]; icol++){
	    long kr=0;
	    for(long ix=0; ix<A->nx; ix++){
		if(!isempty(Ap[iy][ix])){
		    memcpy(out->p+((icol+jcol)*nx+kr),
			   Ap[iy][ix]->p+icol*nxs[ix],
			   nxs[ix]*sizeof(T));
		}
		kr+=nxs[ix];
	    }
	}
	jcol+=nys[iy];
    }
    free(nxs);
    free(nys);
    return out;
}

/**
   convert a vector to cell using dimensions specified in dims. Reference the vector
*/
X(cell)* X(2cellref)(const X(mat) *A, long*dims, long ndim){
    long nx=0;
    for(int ix=0; ix<ndim; ix++){
	nx+=dims[ix];
    }
    if(A->ny!=1){
	error("Use d2cell2 instead for non vectors\n");
    }
    if(nx!=A->nx ){
	error("Shape doesn't agree. nx=%ld, nx=%ld\n", nx,A->nx);
    }
    long kr=0;
    X(cell) *B=X(cellnew)(ndim,1);
    B->m=X(ref)(A);
    for(long ix=0; ix<ndim; ix++){
	B->p[ix]=X(new_ref)(dims[ix],1,A->p+kr);
	kr+=dims[ix];
    }
    return B;
}

/**
   make A a cell array using shape information from ref if *B is NULL
*/
void X(2cell)(X(cell) **B, const X(mat) *A, const X(cell) *ref){
    long nx,ny,*nxs,*nys;
    if(*B) ref=*B;/*use B as reference. */
    X(celldim)(ref, &nx, &ny, &nxs, &nys);
    if(nx!=A->nx || ny!=A->ny){
	error("Shape doesn't agree. Reference is %ldx%ld but input is %ldx%ld\n",
	      nx,ny,A->nx,A->ny);
    }
    if(!*B){
	*B=X(cellnew2)(ref);
	PCELL((*B),Bp);
	for(long iy=0; iy<ref->ny; iy++){
	    for(long ix=0; ix<ref->nx; ix++){
		Bp[iy][ix]=X(new)(nxs[ix],nys[iy]);
	    }
	}
    }
    PCELL((*B),Bp);
    long jcol=0;
    for(long iy=0; iy<ref->ny; iy++){
	for(long icol=0; icol<nys[iy]; icol++){
	    long kr=0;
	    for(long ix=0; ix<ref->nx; ix++){
		if(nxs[ix]>0){
		    memcpy(Bp[iy][ix]->p+icol*Bp[iy][ix]->nx,
			   A->p+((icol+jcol)*nx+kr),
			   nxs[ix]*sizeof(T));
		    kr+=nxs[ix];
		}
	    }
	}
	jcol+=nys[iy];
    }
    free(nxs);
    free(nys);
}

/**
   drop empty blocks (zero). Size of B is not modified.
*/
void X(celldropzero)(X(cell) *B, R thres){
    PCELL(B,Bp);
    for(long iy=0; iy<B->ny; iy++){
	for(long ix=0; ix<B->nx; ix++){
	    X(mat) *tmp=Bp[iy][ix];
	    if(!tmp) continue;
	    int hasnonzero=0;
	    for(int ixy=0; ixy<tmp->nx*tmp->ny; ixy++){
		if(ABS(tmp->p[ixy])>thres){
		    hasnonzero=1;
		    break;
		}
	    }
	    if(!hasnonzero){
		X(free)(Bp[iy][ix]);
		Bp[iy][ix]=NULL;
		/*warning("Dropped block (%ld, %ld)\n", ix, iy); */
	    }
	}
    }
}

/**
   compute ||A-B||/||A||
   use mean.
*/
R X(celldiff)(const X(cell) *A, const X(cell) *B){
    X(cell) *C=NULL;
    X(cellcp)(&C,A);
    X(celladd)(&C,1,B,-1);
    R d=sqrt(X(cellnorm2)(C)*2/(X(cellnorm2)(C)+X(cellnorm2)(B)));
    return isnan(d)?0:d;
}

/**
   clip a X(cell) array to max at 'max', min at 'min'
*/
int X(cellclip)(X(cell) *Ac, R min, R max){
    if(!Ac || !Ac->p) return 0;
    if(!isfinite(min)==-1 && !isfinite(max)==1) return 0;
    int nclip=0;
    for(long i=0; i<Ac->nx*Ac->ny; i++){
	nclip+=X(clip)(Ac->p[i],min,max);
    }
    return nclip;
}


/**
   Multiply a cell with a sparse cell.

  \f$C0+=A*B*alpha\f$.
*/
#ifndef USE_SINGLE
void X(cellmulsp)(X(cell) **C0, const X(cell) *A, const Y(spcell) *B, R alpha){
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
		if(A->p[ix*ax+iz*az] && B->p[iz*bz+iy*by]){
		    X(mulsp)(&C->p[ix+iy*nx],A->p[ix*ax+iz*az], 
			     B->p[iz*bz+iy*by],alpha);
		}
	    }
	}
    }
}
#endif
/**
   add a to diagonal elements of A;
*/
void X(celladdI)(X(cell) *A, R a){
    if(A->nx!=A->ny) 
	error("A must be symmetric\n");
    for(int ib=0; ib<A->nx; ib++){
	X(addI)(A->p[ib+ib*A->nx], a);
    }
}

/**
   raise each cell in the cell array to power of power.
*/
void X(cellcwpow)(X(cell)*A, R power){
    for(long ib=0; ib<A->nx*A->ny; ib++){
	X(cwpow)(A->p[ib],power);
    }
}
/**
   Create a new sub cell matrix of nx*ny starting from(sx,sy)
*/
X(cell) *X(cellsub)(const X(cell) *in, long sx, long nx, long sy, long ny){
    if(nx<=0){
	nx=in->nx-sx;
    }
    if(ny<=0){
	ny=in->ny-sy;
    }
    X(cell)*out=X(cellnew)(nx, ny);
    if(sx+nx>in->nx || sy+ny>in->ny){
	error("Invalid parameter range\n");
    }
    PCELL(in, pin);
    PCELL(out, pout);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    pout[iy][ix]=pin[iy+sy][ix+sx];
	}
    }
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
		    p00=2*(2*p[iy][ix]-p[iy][ix+1])-(2*p[iy+1][ix]-p[iy+1][ix+1]);/*from a */
		}else{
		    p00=2*p[iy][ix-1]-p[iy+1][ix-1];/*from b */
		}
		p01=2*p[iy][ix]-p[iy+1][ix];
		p02=2*p[iy][ix+1]-p[iy+1][ix+1];
		if(ix==nx-2){
		    p03=2*(p[iy][ix+1]*2-p[iy][ix])-(p[iy+1][ix+1]*2-p[iy+1][ix]);/*from n */
		}else{
		    p03=2*p[iy][ix+2]-p[iy+1][ix+2];/*from m */
		}
	    }else{
		if(ix==0){
		    p00=2*p[iy-1][ix]-p[iy-1][ix+1];/*a from b */
		}else{
		    p00=p[iy-1][ix-1];/*b */
		}
		p01=p[iy-1][ix];
		p02=p[iy-1][ix+1];
		if(ix==nx-2){
		    p03=p[iy-1][ix+1]*2-p[iy-1][ix];/*n from m */
		}else{
		    p03=p[iy-1][ix+2];/*m */
		}
	    }
	    if(ix==0){
		p10=p[iy][ix]*2-p[iy][ix+1];/*from c */
	    }else{
		p10=p[iy][ix-1];/*c */
	    }
	    p11=p[iy][ix];
	    p12=p[iy][ix+1];
	    if(ix==nx-2){
		p13=p[iy][ix+1]*2-p[iy][ix];/*from d */
	    }else{
		p13=p[iy][ix+2];/*d */
	    }
	    if(ix==0){
		p20=p[iy+1][ix]*2-p[iy+1][ix+1];/*from e */
	    }else{
		p20=p[iy+1][ix-1];/*e */
	    }
	    p21=p[iy+1][ix];
	    p22=p[iy+1][ix+1];
	    if(ix==nx-2){
		p23=p[iy+1][ix+1]*2-p[iy+1][ix];/*from f */
	    }else{
		p23=p[iy+1][ix+2];/*f */
	    }
	    if(iy==ny-2){
		if(ix==0){
		    p30=2*(p[iy+1][ix]*2-p[iy+1][ix+1])-(p[iy][ix]*2-p[iy][ix+1]);/*from h */
		}else{
		    p30=2*p[iy+1][ix-1]-p[iy][ix-1];/*from g */
		}
		p31=2*p[iy+1][ix]-p[iy][ix];
		p32=2*p[iy+1][ix+1]-p[iy][ix+1];
		if(ix==nx-2){
		    p33=2*(2*p[iy+1][ix+1]-p[iy+1][ix])-(2*p[iy][ix+1]-p[iy][ix]);/*from j */
		}else{
		    p33=2*p[iy+1][ix+2]-p[iy][ix+2];/*from i */
		}
	    }else{
		if(ix==0){
		    p30=p[iy+2][ix]*2-p[iy+2][ix+1];/*h from g */
		}else{
		    p30=p[iy+2][ix-1];/*g */
		}
		p31=p[iy+2][ix];
		p32=p[iy+2][ix+1];
		if(ix==nx-2){
		    p33=2*p[iy+2][ix+1]-p[iy+2][ix];/*j from i */
		}else{
		    p33=p[iy+2][ix+2];/*i */
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
    T xsep1=(R)(nx-1)/(x->p[nx-1]-xmin);
    T ysep1=(R)(ny-1)/(y->p[ny-1]-ymin);
    assert(xnew->nx == ynew->nx && xnew->ny == ynew->ny);
    X(mat)*zz=X(new)(xnew->nx, xnew->ny);
    PCELL(coeff,pc);
    for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
	R xm=REAL((xnew->p[ix]-xmin)*xsep1);
	long xmf=floor(xm);
	if(xmf<0) xmf=0;
	if(xmf>nx-2) xmf=nx-2;
	xm=xm-xmf;

	R ym=REAL((ynew->p[ix]-ymin)*ysep1);
	long ymf=floor(ym);
	if(ymf<0) ymf=0;
	if(ymf>ny-2) ymf=ny-2;
	ym=ym-ymf;
	
	T xm2=xm *xm;
	T xm3=xm2*xm;
	T ym2=ym *ym;
	T ym3=ym2*ym;
	PMAT(pc[ymf][xmf],ppc);
	zz->p[ix]= ppc[0][0] + ppc[0][1] * xm + ppc[0][2] * xm2 + ppc[0][3] * xm3 +
	    ppc[1][0] * ym + ppc[1][1] * ym * xm + ppc[1][2] * ym * xm2 + ppc[1][3] * ym * xm3 +
	    ppc[2][0] * ym2 + ppc[2][1] * ym2 * xm + ppc[2][2] * ym2 * xm2 + ppc[2][3] * ym2 * xm3 +
	    ppc[3][0] * ym3 + ppc[3][1] * ym3 * xm + ppc[3][2] * ym3 * xm2 + ppc[3][3] * ym3 * xm3;

    }
    return zz;
}
