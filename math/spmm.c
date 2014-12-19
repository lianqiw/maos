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
  Handling multiplication between sparse and dense matrix/vector.
*/
typedef struct sp_thread_t {
    const X(sp) *A;
    const T *x;
    T **ytmp;
    T *y;
    T alpha;
    long nthread;
}sp_thread_t;

/**
 * Multiply a sparse matrix with the real part of a complex vector: y=y+A*creal(x)*alpha
 */
void X(spmulcreal)(T *restrict y, const X(sp) *A, const RI * restrict x, T alpha){
    if(A && x){
	for(long icol=0; icol<A->n; icol++){
	    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		y[A->i[ix]]+=alpha*A->x[ix]*creal(x[icol]);
	    }
	}

    }
}

/**
 * sparse matrix multiply with a vector: y=y+op(A)*x*alpha
 op(A)=A  if trans=='n'
 op(A)=A' if trans=='t'
 op(A)=CONJ(A') if trans=='c'
*/
void X(spmulvec)(T *restrict y, const X(sp) *A, const T * restrict x, char trans, T alpha){
    if(A && x){
	assert(y);
	if(trans=='n'){
	    for(long icol=0; icol<A->n; icol++){
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    y[A->i[ix]]+=alpha*A->x[ix]*x[icol];
		}
	    }
	}else if(trans=='t'){
	    OMPTASK_FOR(icol, 0, A->n){
		T tmp=0;
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    tmp+=alpha*(A->x[ix])*x[A->i[ix]];
		}
		y[icol]+=tmp;
	    }
	    OMPTASK_END;
	}else if(trans=='c'){
	    OMPTASK_FOR(icol, 0, A->n){
		T tmp=0;
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    tmp+=alpha*CONJ(A->x[ix])*x[A->i[ix]];
		}
		y[icol]+=tmp;
	    }
	    OMPTASK_END;
	}else{
	    error("Invalid trans=%c\n", trans);
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
	X(init)(yout, x->nx, A->n);
	X(mat) *y=*yout;
	assert(x->nx==y->nx && x->ny==A->m);
	if(x->nx==1){
	    X(spmulvec)(y->p, A, x->p, 't', alpha);
	}else{
	    PMAT(y,Y); PMAT(x,X);
	    if(ABS(alpha-1)<1.e-100){
		for(long icol=0; icol<A->n; icol++){
		    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(long jcol=0; jcol<y->nx; jcol++){
			    Y[icol][jcol]+=A->x[ix]*X[A->i[ix]][jcol];
			}
		    }
		}
	    }else{
		for(long icol=0; icol<A->n; icol++){
		    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(long jcol=0; jcol<y->nx; jcol++){
			    Y[icol][jcol]+=alpha*A->x[ix]*X[A->i[ix]][jcol];
			}
		    }
		}
	    }
	}
    }
}

/**
   sparse matrix multiply with dens matrix X(mat): y=y+op(A)*x*alpha
   op(A)=A  if trans=='n'
   op(A)=A' if trans=='t'
   op(A)=CONJ(A') if trans=='c'
*/
void X(spmm)(X(mat) **yout, const X(sp) *A, const X(mat) *x, char trans, const T alpha){
    if(A&&x){
	/* y=y+alpha*A*x; */
	long xnx, ynx;
	if(trans=='n'){
	    ynx=A->nx;
	    xnx=A->ny;
	}else{
	    ynx=A->ny;
	    xnx=A->nx;
	}
	X(init)(yout, ynx, x->ny);
	X(mat) *y=*yout;
	if(xnx!=x->nx){
	    error("mismatch: xnx=%ld, x->nx=%ld\n", xnx, x->nx);
	}
	if(x->ny==1){
	    X(spmulvec)(y->p, A, x->p, trans, alpha);
	}else{
	    PMAT(y,Y);
	    PMAT(x,X);
	    if(trans=='n'){
		for(long icol=0; icol<A->n; icol++){
		    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(long jcol=0; jcol<y->ny; jcol++){
			    Y[jcol][A->i[ix]]+=alpha*A->x[ix]*X[jcol][icol];
			}
		    }
		}
	    }else if(trans=='t'){
		OMPTASK_FOR(icol, 0, A->n){
		    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(long jcol=0; jcol<y->ny; jcol++){
			    Y[jcol][icol]+=alpha*(A->x[ix])*X[jcol][A->i[ix]];
			}
		    }
		}
		OMPTASK_END;
	    }else if(trans=='c'){
		OMPTASK_FOR(icol, 0, A->n){
		    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
			for(long jcol=0; jcol<y->ny; jcol++){
			    Y[jcol][icol]+=alpha*CONJ(A->x[ix])*X[jcol][A->i[ix]];
			}
		    }
		}
		OMPTASK_END;
	    }else{
		error("Invalid trans=%c\n", trans);
	    }
	}
    }
}

/**
   Multiply a cell with a sparse cell.

   \f$C0+=A*B*alpha\f$. where C0, and A are dense.
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
    if(ABS(scale-1)>EPS){
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
   Compute A*B and add to C0.
   C0=C0+op(A)*op(B)*alpha;
   op(A)=A  if trans[0]=='n'
   op(A)=A' if trans[0]=='t'
   op(A)=CONJ(A') if trans[0]=='c'

   op(B)=B  if trans[1]=='n'
   op(B)=B' if trans[1]=='t'
   op(B)=CONJ(B') if trans[1]=='c'
   
   A may be dense or sparse matrix.
*/
/*
  2009-11-09: There was initially a beta parameter It was implemented wrongly
  for beta!=1 because for every call to dmm, the already accumulated ones are
  scaled.  removed beta.
*/
void X(cellmm)(X(cell) **C0, const void *A_, const X(cell) *B, const char trans[2], const R alpha){
    const cell *A=(cell*)A_;
    if(!A || !B) return;
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    if(trans[0]=='n'){
	nx=A->nx; 
	ax=1; az=A->nx;
	nz=A->ny;
    }else if(trans[0]=='c' || trans[0]=='t'){ 
	nx=A->ny;
	az=1; ax=A->nx;
	nz=A->nx;
    }else{
	error("Invalid trans[0]=%c\n", trans[0]);
	nx=0; az=0; ax=0; nz=0;
    }
    if(trans[1]=='n'){
	ny=B->ny; 
	bz=1; by=B->nx;
	if(nz!=B->nx) error("mismatch\n");
    }else if(trans[1]=='c' || trans[1]=='t'){
	ny=B->nx;
	by=1; bz=B->nx;
	if(nz!=B->ny) error("mismatch\n");
    }else{
	error("Invalid trans[1]=%c\n", trans[1]);
	ny=0; by=0; bz=0;
    }
    X(cellinit)(C0, nx, ny);
    X(cell) *C=*C0;
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
#if _OPENMP >= 200805
#pragma omp task firstprivate(ix,iy) if(nx*ny>1)
#endif
	    for(int iz=0; iz<nz; iz++){
		if(A->p[ix*ax+iz*az] && B->p[iz*bz+iy*by]){
		    switch(A->p[ix*ax+iz*az]->id){
		    case M_T://dense A
			X(mm)(C->p+ix+iy*nx,1.,X(mat_cast)(A->p[ix*ax+iz*az]), B->p[iz*bz+iy*by],trans,alpha);
			break;
		    case M_SPT:
			if(trans[1]!='n'){
			    error("Not implemented\n");
			}
			X(spmm)(C->p+ix+iy*nx,(X(sp*))A->p[ix*ax+iz*az], B->p[iz*bz+iy*by],trans[0],alpha);
			break;
		    }
		}
	    }
	}
    }
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
}
