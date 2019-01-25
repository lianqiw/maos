/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	if(nz!=B->nx) error("mismatch: trans=%s, A is %ldx%ld, B is %ldx%ld\n", trans, A->nx, A->ny, B->nx, B->ny);
    }else if(trans[1]=='c' || trans[1]=='t'){
	ny=B->nx;
	by=1; bz=B->nx;
	if(nz!=B->ny) error("mismatch: trans=%s, A is %ldx%ld, B is %ldx%ld\n", trans, A->nx, A->ny, B->nx, B->ny);
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
void X(spmulcreal)(T *restrict y, const X(sp) *A, const RI * restrict x, R alpha){
    if(A && x){
	for(long icol=0; icol<A->ny; icol++){
	    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		y[A->i[ix]]+=alpha*(A->x[ix]*creal(x[icol]));
	    }
	}

    }
}

/**
 * sparse matrix multiply with a vector: y=y+op(A)*x*alpha
 op(A)=A  if trans=='n'
 op(A)=A' if trans=='t'
 op(A)=conj(A') if trans=='c'
*/
void X(spmulvec)(T *restrict y, const X(sp) *A, const T * restrict x, char trans, T alpha){
    if(A && x){
	assert(y);
	if(trans=='n'){
	    for(long icol=0; icol<A->ny; icol++){
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    y[A->i[ix]]+=alpha*A->x[ix]*x[icol];
		}
	    }
	}else if(trans=='t'){
	    OMPTASK_FOR(icol, 0, A->ny){
		T tmp=0;
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    tmp+=alpha*(A->x[ix])*x[A->i[ix]];
		}
		y[icol]+=tmp;
	    }
	    OMPTASK_END;
	}else if(trans=='c'){
	    OMPTASK_FOR(icol, 0, A->ny){
		T tmp=0;
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){
		    tmp+=alpha*conj(A->x[ix])*x[A->i[ix]];
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
   Unified sparse matrix multiply with dense matrix to handle sparse*dense or dense*sparse with flexibility.
   op(y)=op(y)+op(A)*op(x)*alpha
   op(A)=A  if trans=='n'
   op(A)=A' if trans=='t'
   op(A)=conj(A') if trans=='c'
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

#define no_conj(A) (A)
#define do_conj(A) conj(A)
#define no_trans(A,i,j) P(A,i,j)
#define do_trans(A,i,j) P(A,j,i)
#define LOOP_NORMA(py, yny,  conjA, px, conjx)				\
	for(long icol=0; icol<A->ny; icol++){				\
	    for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){		\
		for(long jcol=0; jcol<yny; jcol++){			\
		    py(y, A->i[ix], jcol)+=alpha*conjA(A->x[ix])*conjx(px(x, icol, jcol)); \
		}							\
	    }								\
	}

#define LOOP_TRANSA(py, yny, conjA, px, conjx)				\
	{								\
	    OMPTASK_FOR(icol, 0, A->ny){					\
		for(long ix=A->p[icol]; ix<A->p[icol+1]; ix++){		\
		    for(long jcol=0; jcol<yny; jcol++){			\
			py(y, icol, jcol)+=alpha*conjA(A->x[ix])*conjx(px(x, A->i[ix], jcol)); \
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
 * Multiply two sparse arrays and add to the third: C0=C0+op(A)*op(B)*scale*/
void X(spmulsp2)(X(sp) **C0, const X(sp) *A, const X(sp) *B, const char trans[2], const T scale){
    /*return C=C+ alpha*(A*B) */
    if(!A || !B) return;
    X(sp)* At=0;
    X(sp)* Bt=0;
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
    if(fabs(scale-(T)1.)>EPS){
	X(spscale)(res, scale);
    }
    if(!*C0) 
	*C0=res;
    else{
	X(spadd)(C0,1,res,1);
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
   Compute A*B and add to C0.
   C0=C0+op(A)*op(B)*alpha;
   op(A)=A  if trans[0]=='n'
   op(A)=A' if trans[0]=='t'
   op(A)=conj(A') if trans[0]=='c'

   op(B)=B  if trans[1]=='n'
   op(B)=B' if trans[1]=='t'
   op(B)=conj(B') if trans[1]=='c'
   
   A may be dense or sparse matrix.
*/
/*
  2009-11-09: There was initially a beta parameter It was implemented wrongly
  for beta!=1 because for every call to dmm, the already accumulated ones are
  scaled.  removed beta.
*/
void X(cellmm)(void *C0_, const void *A_, const void *B_, const char trans[2], const R alpha){
    if(!A_ || !B_) return;
    const cell *A=(const cell*)(A_);
    const cell *B=(const cell*)(B_);
    cell **C0=(cell**)C0_;
    if(iscell(A) && iscell(B)){
	//multiplication of cells.
	mm_t D=parse_trans(A, B, trans);
	cellinit(C0, D.nx, D.ny);
	cell *C=*C0;
	for(int iy=0; iy<D.ny; iy++){
	    for(int ix=0; ix<D.nx; ix++){
/*#if _OPENMP >= 200805
#pragma omp task firstprivate(ix,iy) if(D.nx*D.ny>1 && omp_in_parallel())
#endif*/
		for(int iz=0; iz<D.nz; iz++){
		    X(cellmm)(C->p+ix+iy*D.nx, A->p[ix*D.ax+iz*D.az], B->p[iz*D.bz+iy*D.by], trans, alpha);
		}
/*#if _OPENMP >= 200805
#pragma omp taskwait
#endif*/
	    }
	}
	
    }else{
	//multiplication between dense or sparse 
	if(ismat(B)){//dense B
	    if(ismat(A)){ //dense A
		X(mm)((X(mat)**)C0, 1, (X(mat)*)A, (X(mat)*)B, trans, alpha);
	    }else if(issp(A)){//sparse A
		X(spmm)((X(mat)**)C0, (X(sp)*)A, (X(mat)*)B, trans, alpha);
	    }else{
		error("Invalid A type, id=%u.\n", A->id);
	    }
	}else if(issp(B)){//sparse B
	    if(ismat(A)){//dense A
		X(mulsp)((X(mat)**)C0, (X(mat)*)A, (X(sp)*)B, trans, alpha);
	    }else if(issp(A)){//sparse A
		if(!*C0 || issp(*C0)){//result is sparse
		    X(spmulsp2)((X(sp)**)C0, (X(sp)*)A, (X(sp)*)B, trans, alpha);
		}else if(ismat(*C0)){//result is dense
		    X(mat)*Afull=0; X(spfull)(&Afull, (X(sp)*)A, 'n',1);
		    X(mat)*Bfull=0; X(spfull)(&Bfull, (X(sp)*)B, 'n',1);
		    X(mm)((X(mat)**)C0, 1, Afull, Bfull, trans, alpha);
		    X(free)(Afull);
		    X(free)(Bfull);
		}else{
		    error("Invalid C type, id=%u. \n", (*C0)->id);
		}
	    }else{
		error("Invalid A type, id=%u.\n", A->id);
	    }
	}else{
	    error("Invalid B type, id=%u.\n", B->id);
	}
    }
}

cell *X(cellmm2)(const void *A_, const void *B_, const char trans[2]){
    cell *res=0;
    X(cellmm)(&res, A_, B_, trans, 1);
    return res;
}

void X(celladdI)(void *A_, T alpha){
    if(!A_) return;
    cell *A=cell_cast(A_);
    assert(A->nx==A->ny);
    for(int ii=0; ii<A->ny; ii++){
	if(!P(A, ii, ii)){
	    continue;
	}else if(ismat(P(A, ii, ii))){
	    X(addI)(X(mat_cast)(P(A,ii,ii)), alpha);
	}else if(issp(P(A, ii, ii))){
	    X(spaddI)(X(sp_cast)(P(A,ii,ii)),alpha);
	}else{
	    error("Invalid id=%u", P(A,ii,ii)->id);
	}
    }
}
/**
   Takes parameters of X(mat), X(sp), X(cell), X(spcell): A=A*ac+B*bc;
 */
void X(celladd)(void *A_, R ac, const void *B_, R bc){
    if(!A_ || !B_ || bc==0) return;
    cell *B=(cell*)(B_);
    cell **pA=(cell**)A_;
    if(iscell(B)){//cell
	cellinit2(pA, B);
	cell *A=*pA;
	for(int i=0; i<B->nx*B->ny; i++){
	    X(celladd)(A->p+i, ac, B->p[i], bc);
	}
    }else{//non cell
	if(!*pA || ismat(*pA)){//A is dense
	    if(ismat(B)){//Add dense to dense
		X(add)((X(mat)**)pA, ac, (X(mat)*)B, bc);
	    }else{//add sparse to dense
		if(ac!=1){
		    X(scale)((X(mat*))*pA, ac);
		}
		X(spfull)((X(mat)**)(pA), (X(sp)*)B, 'n',bc);
	    }
	}else if(issp(*pA)){
	    if(issp(B)){//add sparse to sparse
		X(sp)* tmp=X(spadd2)((X(sp)*)(*pA), ac, (X(sp)*)B, bc);
		X(spmove)((X(sp)*)(*pA), (X(sp)*)tmp);
		free(tmp);
	    }else{
		error("Adding dense to sparse matrix is not allowed.\n");
	    }
	}else{
	    error("Invalid operand\n");
	}
    }
}
/**
   Takes parameters of X(mat), X(sp), X(cell), X(spcell): Copy B to A;
 */
void X(cellcp)(void *A_, const void *B_){
    if(!B_){
	X(cellscale)(*((cell**)A_), 0);
    }else{
	X(celladd)(A_, 0, B_, 1);
    }
}

/**
   scale each element of A.
*/
void X(cellscale)(void *A_, R w){
    if(!A_) return;
    cell *A=(cell*)A_;
    if(iscell(A_)){
	for(int i=0; i<A->nx*A->ny; i++){
	    X(cellscale)(A->p[i],w);
	}
    }else{
	if(ismat(A)){
	    X(scale)((X(mat*))A, w);
	}else if(issp(A)){
	    X(spscale)((X(sp*))A, w);
	}else{
	    error("Invalid type: id=%u\n", A->id);
	}
    }
}


/**
   setting all elements of a X(cell) to zero.
*/
void X(cellzero)(void *dc){
    X(cellscale)(dc, 0);
}

/**
   Convert X(cell) or X(spcell) to X(mat)
*/
X(mat)* X(cell2m)(const void *A_){
    if(!A_) return 0;
    cell *A=cell_cast(A_);
    long nx,ny,*nxs,*nys;
    celldim(A,&nx,&ny,&nxs,&nys);
    if(!nx || !ny) return 0;
    X(mat) *out=X(new)(nx,ny);
    long jcol=0;
    for(long iy=0; iy<A->ny; iy++){
	for(long icol=0; icol<nys[iy]; icol++){
	    long kr=0;
	    for(long ix=0; ix<A->nx; ix++){
		if(!isempty(P(A,ix,iy))){
		    T *pout=out->p+((icol+jcol)*nx+kr);
		    if(ismat(P(A,ix,iy))){
			memcpy(pout, ((X(mat*))P(A,ix,iy))->p+icol*nxs[ix], nxs[ix]*sizeof(T));
		    }else if(issp(P(A, ix, iy))){
			//convert sparse col to full
			X(sp*)Asp=(X(sp*))P(A,ix,iy);
			for(long j=Asp->p[icol]; j<Asp->p[icol+1]; j++){
			    pout[Asp->i[j]]=Asp->x[j];
			}			
		    }
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
