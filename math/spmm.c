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
typedef struct mm_t{
    int nx, ny, nz;
    int ax, az;
    int bz, by;
}mm_t;
/*
  determine parameters to perform for loop based on trans
*/
static mm_t parse_trans(const cell *A, const cell *B, const char trans[2]){
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    if(trans[0]=='n' || trans[0]=='N'){
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
    if(trans[1]=='n' || trans[1]=='N'){
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
    mm_t out={nx, ny, nz, ax, az, bz, by};
    return out;
}
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
/*
void X(mulsp)(X(mat) **yout, const X(mat) *x,const X(sp) *A, const T alpha){
    if(A&&x){
	X(init)(yout, x->nx, A->n);
	X(mat) *y=*yout;
	assert(x->nx==y->nx && x->ny==A->m);
	if(x->nx==1){
	    X(spmulvec)(y->p, A, x->p, 't', alpha);
	}else{
	    PMAT(y,Y); PMAT(x,X);
	    if(ABS(alpha-(T)1.)<1.e-100){
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
*/
/**
   Unified sparse matrix multiply with dense matrix to handle sparse*dense or dense*sparse with flexibility.
   op(y)=op(y)+op(A)*op(x)*alpha
   op(A)=A  if trans=='n'
   op(A)=A' if trans=='t'
   op(A)=CONJ(A') if trans=='c'
*/
static void X(spmm_do)(X(mat) **yout, const X(sp) *A, const X(mat) *x, const char trans[2], const int transy, const T alpha){
    if(!A || !x) return;
    mm_t D=parse_trans((cell*)A, (cell*)x, trans);
    if(transy){
	X(init)(yout, D.ny, D.nx);
    }else{
	X(init)(yout, D.nx, D.ny);
    }
    X(mat) *y=*yout;
    if(x->ny==1 && trans[1]=='n' && transy==0){
	X(spmulvec)(y->p, A, x->p, trans[0], alpha);
    }else{
	PMAT(y,Y);
	PMAT(x,X);

#define no_conj(A) (A)
#define do_conj(A) CONJ(A)
#define no_trans(A,i,j) A[j][i]
#define do_trans(A,i,j) A[i][j]
#define LOOP_NORMA(py, yny,  conjA, px, conjx)				\
	for(long icol=0; icol<A->n; icol++){				\
	    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){		\
		for(long jcol=0; jcol<yny; jcol++){			\
		    py(Y, A->i[ix], jcol)+=alpha*conjA(A->x[ix])*conjx(px(X, icol, jcol)); \
		}							\
	    }								\
	}

#define LOOP_TRANSA(py, yny, conjA, px, conjx)				\
	{								\
	    OMPTASK_FOR(icol, 0, A->n){					\
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){		\
		    for(long jcol=0; jcol<yny; jcol++){			\
			py(Y, icol, jcol)+=alpha*conjA(A->x[ix])*conjx(px(X, A->i[ix], jcol)); \
		    }							\
		}							\
	    }								\
	    OMPTASK_END;						\
	}
	if(transy){
	    if(trans[0]=='n'){
		switch(trans[1]){
		case 'n':
		    LOOP_NORMA(do_trans, y->ny, no_conj, no_trans, no_conj);break;
		case 'N'://conjuate, but not transpose
		    LOOP_NORMA(do_trans, y->nx, no_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_NORMA(do_trans, y->nx, no_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_NORMA(do_trans, y->nx, no_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='N'){//conjuate, but not transpose
		switch(trans[1]){
		case 'n':
		    LOOP_NORMA(do_trans, y->nx, do_conj, no_trans, no_conj);break;
		case 'N'://conjuate, but not transpose
		    LOOP_NORMA(do_trans, y->nx, do_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_NORMA(do_trans, y->nx, do_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_NORMA(do_trans, y->nx, do_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='t'){
		switch(trans[1]){
		case 'n':
		    LOOP_TRANSA(do_trans, y->nx, no_conj, no_trans, no_conj);break;
		case 'N':
		    LOOP_TRANSA(do_trans, y->nx, no_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_TRANSA(do_trans, y->nx, no_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_TRANSA(do_trans, y->nx, no_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='c'){
		switch(trans[1]){
		case 'n':
		    LOOP_TRANSA(do_trans, y->nx, do_conj, no_trans, no_conj);break;
		case 'N':
		    LOOP_TRANSA(do_trans, y->nx, do_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_TRANSA(do_trans, y->nx, do_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_TRANSA(do_trans, y->nx, do_conj, do_trans, do_conj);break;
		}
	    }else{
		error("Invalid trans=%s\n", trans);
	    }
	}else{
	    if(trans[0]=='n'){
		switch(trans[1]){
		case 'n':
		    LOOP_NORMA(no_trans, y->ny, no_conj, no_trans, no_conj);break;
		case 'N'://conjuate, but not transpose
		    LOOP_NORMA(no_trans, y->ny, no_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_NORMA(no_trans, y->ny, no_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_NORMA(no_trans, y->ny, no_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='N'){//conjuate, but not transpose
		switch(trans[1]){
		case 'n':
		    LOOP_NORMA(no_trans, y->ny, do_conj, no_trans, no_conj);break;
		case 'N'://conjuate, but not transpose
		    LOOP_NORMA(no_trans, y->ny, do_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_NORMA(no_trans, y->ny, do_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_NORMA(no_trans, y->ny, do_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='t'){
		switch(trans[1]){
		case 'n':
		    LOOP_TRANSA(no_trans, y->ny, no_conj, no_trans, no_conj);break;
		case 'N':
		    LOOP_TRANSA(no_trans, y->ny, no_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_TRANSA(no_trans, y->ny, no_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_TRANSA(no_trans, y->ny, no_conj, do_trans, do_conj);break;
		}
	    }else if(trans[0]=='c'){
		switch(trans[1]){
		case 'n':
		    LOOP_TRANSA(no_trans, y->ny, do_conj, no_trans, no_conj);break;
		case 'N':
		    LOOP_TRANSA(no_trans, y->ny, do_conj, no_trans, do_conj);break;
		case 't':
		    LOOP_TRANSA(no_trans, y->ny, do_conj, do_trans, no_conj);break;
		case 'c':
		    LOOP_TRANSA(no_trans, y->ny, do_conj, do_trans, do_conj);break;
		}
	    }else{
		error("Invalid trans=%s\n", trans);
	    }
	}
    }
}
/**
   sparse matrix multiply with dense matrix
*/
void X(spmm)(X(mat) **yout, const X(sp) *A, const X(mat) *x, const char trans[2], const T alpha){
    X(spmm_do)(yout, A, x, trans, 0, alpha);
}
static char reverse_trans(char in){
    switch(in){
    case 'n':
	return 't';
    case 't':
	return 'n';
    case 'c':
	return 'N';
    default:
	error("Invalid\n");
	return 'a';
    }
}
/**
   dense matrix multiply with sparse matrix, implemented by transpose all: 
   y=x*A-->y'=A'*x'; conjugation is handled carefully.
 */
void X(mulsp)(X(mat) **yout, const X(mat) *x, const X(sp) *A, const char trans[2], const T alpha){
    char trans2[2];
    trans2[1]=reverse_trans(trans[0]);
    trans2[0]=reverse_trans(trans[1]);
    X(spmm_do)(yout, A, x, trans2, 1, alpha);
}
/**
   Multiply a cell with a sparse cell.

   \f$C0+=A*B*alpha\f$. where C0, and A are dense.
*/
/*
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
    }*/

/**
 * Multiply two sparse arrays and add to the third: C0=C0+op(A)*op(B)*scale*/
void X(spmulsp2)(X(sp) **C0, const X(sp) *A, const X(sp) *B, const char trans[2], const T scale){
    /*return C=C+ alpha*(A*B) */
    if(!A || !B) return;
    X(sp)* At=0;
    X(sp)* Bt=0;
    int transa=0, transb=0;
    if(trans[0]=='t'||trans[0]=='c'){
	At=X(sptrans)(A);
	if(trans[0]=='c'){
	    X(spconj)(At);
	}
    }
    if(trans[1]=='t'||trans[1]=='t'){
	Bt=X(sptrans)(B);
	if(trans[1]=='c'){
	    X(spconj)(Bt);
	}
    }
    X(sp) *res=X(ss_multiply)(At?At:A, Bt?Bt:B);
    X(spsort)(res);
    if(ABS(scale-(T)1.)>EPS){
	X(spscale)(res, scale);
    }
    if(!*C0) 
	*C0=res;
    else{
	X(spadd)(C0, res);
	X(spfree)(res);
    }
    X(spfree)(At);
    X(spfree)(Bt);
    X(ss_dropzeros)(*C0);
}
X(sp) *X(spmulsp)(const X(sp) *A, const X(sp) *B, const char trans[2]){
    X(sp) *res=0;
    X(spmulsp2)(&res, A, B, trans, 1);
    X(ss_dropzeros)(res);
    return res;
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
		X(spmulsp2)(&Cp[iy][ix],Ap[iz][ix],Bp[iy][iz],"nn",scale);
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
void X(cellmm)(void *C0_, const void *A_, const void *B_, const char trans[2], const R alpha){
    if(!A_ || !B_) return;
    const cell *A=cell_cast(A_);
    const cell *B=cell_cast(B_);
    cell **C0=(cell**)C0_;
    mm_t D=parse_trans(A, B, trans);
    cellinit(C0, D.nx, D.ny);
    cell *C=*C0;
    for(int iy=0; iy<D.ny; iy++){
	for(int ix=0; ix<D.nx; ix++){
#if _OPENMP >= 200805
#pragma omp task firstprivate(ix,iy) if(D.nx*D.ny>1)
#endif
	    for(int iz=0; iz<D.nz; iz++){
		if(A->p[ix*D.ax+iz*D.az] && B->p[iz*D.bz+iy*D.by]){
		    switch(B->p[iz*D.bz+iy*D.by]->id){
		    case M_T: //dense B
			{
			    X(mat)*Bi_d=X(mat_cast)(B->p[iz*D.bz+iy*D.by]);
			    X(mat)** Ci_d=(X(mat)**)(C->p+ix+iy*D.nx);
			    switch(A->p[ix*D.ax+iz*D.az]->id){
			    case M_T://dense A * dense B
				{
				    X(mat)*Ai_d=X(mat_cast)(A->p[ix*D.ax+iz*D.az]);
				    X(mm)(Ci_d,1.,Ai_d, Bi_d,trans,alpha);
				}
				break;
			    case M_SPT://sparse A * dense B
				{
				    X(sp)*Ai_sp=X(sp_cast)(A->p[ix*D.ax+iz*D.az]);
				    X(spmm)(Ci_d, Ai_sp, Bi_d,trans,alpha);
				}
				break;
			    }
			}
			break;
		    case M_SPT: //sparse B
			{
			    X(sp)*Bi_sp=X(sp_cast)(B->p[iz*D.bz+iy*D.by]);
			    switch(A->p[ix*D.ax+iz*D.az]->id){
			    case M_T://dense A * sparse B. dense C
				{
				    X(mat)**Ci_d=(X(mat)**)(C->p+ix+iy*D.nx);
				    X(mat)* Ai_d=X(mat_cast)(A->p[ix*D.ax+iz*D.az]);
				    X(mulsp)(Ci_d, Ai_d, Bi_sp,trans,alpha);
				}
				break;
			    case M_SPT://sparse A * sparse B. sparse or dense C.
				{
				    X(sp)*Ai_sp=X(sp_cast)(A->p[ix*D.ax+iz*D.az]);
				    if(!C->p[ix+iy*D.nx] || C->p[ix+iy*D.nx]->id==M_SPT){//C is sparse
					X(sp)** Ci_sp=(X(sp)**)(C->p+ix+iy*D.nx);
					X(spmulsp2)(Ci_sp,Ai_sp, Bi_sp,trans,alpha);
				    }else if(C->p[ix+iy*D.nx]->id==M_T){//C is dense
					X(mat)** Ci_d=(X(mat)**)(C->p+ix+iy*D.nx);
					X(mat)* Afull=0; X(spfull)(&Afull, Ai_sp, 1);
					X(mat)* Bfull=0; X(spfull)(&Bfull, Bi_sp, 1);
					X(mm)(Ci_d,1., Afull, Bfull, trans, alpha);
					X(free)(Afull);
					X(free)(Bfull);
				    }else{
					error("Invalid C id=%ld\n", C->p[ix+iy*D.nx]->id);
				    }
				}
				break;
			    }
			}
			break;
		    default:
			error("Invalid B id=%ld\n", B->p[ix*D.ax+iz*D.az]->id);
		    }
		}
	    }
	}
    }
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
}
void *X(cellmm2)(const void *A_, const void *B_, const char trans[2]){
    void *res=0;
    X(cellmm)(&res, A_, B_, trans, 1);
    return res;
}

void X(celladdI)(void *A_, T alpha){
    if(!A_) return;
    cell *A=cell_cast(A_);
    assert(A->nx==A->ny);
    for(int ii=0; ii<A->ny; ii++){
	if(!IND(A, ii, ii)){
	    continue;
	}else if(IND(A, ii, ii)->id==M_T){
	    X(addI)(X(mat_cast)(IND(A,ii,ii)), alpha);
	}else if(IND(A, ii, ii)->id==M_SPT){
	    X(spaddI)(X(sp_cast)(IND(A,ii,ii)),alpha);
	}else{
	    error("Invalid id=%ld", IND(A,ii,ii)->id);
	}
    }
}
/**
   Takes parameters of X(mat), X(sp), X(cell), X(spcell)
 */
void X(celladd)(void *A_, T ac, const void *B_, T bc){
    if(!A_ || !B_) return;
    cell *B=(cell*)(B_);
    cell **pA=(cell**)A_;
    if(B->id==MCC_ANY){//cell
	cellinit(pA, B->nx, B->ny);
	cell *A=*pA;
	for(int i=0; i<B->nx*B->ny; i++){
	    X(celladd)(A->p+i, ac, B->p[i], bc);
	}
    }else{//non cell
	if(!*pA || (*pA)->id==M_T){//A is dense
	    if(B->id==M_T){
		X(add)((X(mat)**)pA, ac, (X(mat)*)B, bc);
	    }else{
		if(ac!=1){
		    X(scale)((X(mat*))*pA, ac);
		}
		X(spfull)((X(mat)**)(pA), (X(sp)*)B, bc);
	    }
	}else if((*pA)->id==M_SPT){
	    if(B->id==M_SPT){
		X(sp)* tmp=X(spadd2)((X(sp)*)(*pA), ac,
				     (X(sp)*)B, bc);
		X(spmove)((X(sp)*)(*pA), (X(sp)*)B);
		X(spfree)(tmp);
	    }else{
		error("Adding dense to sparse matrix makes a dense\n");
	    }
	}else{
	    error("Invalid operand\n");
	}
    }
}

