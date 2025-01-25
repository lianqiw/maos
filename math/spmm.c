/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "mathdef.h"
#include "defs.h"

typedef struct mm_t{
	int nx, ny, nz;
	int ax, az;
	int bz, by;
}mm_t;
/*
  determine parameters to perform for loop based on trans
*/
static mm_t parse_trans(const_anyarray A_, const_anyarray B_, const char trans[2]){
	const cell *A=A_.c;
	const cell *B=B_.c;
	int ax, az;
	int nx, ny, nz;
	int bz, by;
	if(trans[0]=='n'||trans[0]=='N'){
		nx=A->nx;
		ax=1; az=A->nx;
		nz=A->ny;
	} else if(trans[0]=='c'||trans[0]=='t'){
		nx=A->ny;
		az=1; ax=A->nx;
		nz=A->nx;
	} else{
		error("Invalid trans[0]=%c\n", trans[0]);
		nx=0; az=0; ax=0; nz=0;
	}
	if(trans[1]=='n'||trans[1]=='N'){
		ny=B->ny;
		bz=1; by=B->nx;
		if(nz!=B->nx) error("mismatch: trans=%s, A is %ldx%ld, B is %ldx%ld\n", trans, A->nx, A->ny, B->nx, B->ny);
	} else if(trans[1]=='c'||trans[1]=='t'){
		ny=B->nx;
		by=1; bz=B->nx;
		if(nz!=B->ny) error("mismatch: trans=%s, A is %ldx%ld, B is %ldx%ld\n", trans, A->nx, A->ny, B->nx, B->ny);
	} else{
		error("Invalid trans[1]=%c\n", trans[1]);
		ny=0; by=0; bz=0;
	}
	mm_t out={nx, ny, nz, ax, az, bz, by};
	return out;
}
/*
  Handling multiplication between sparse and dense matrix/vector.
*/
typedef struct sp_thread_t{
	const X(sp)* A;
	const T* x;
	T** ytmp;
	T* y;
	T alpha;
	long nthread;
}sp_thread_t;

/**
 * Multiply a sparse matrix with the real part of a complex vector: y=y+A*creal(x)*alpha
 */
void X(spmulcreal)(T* restrict y, const X(sp)* A, const RC* restrict x, R alpha){
	if(A&&x&&y){
		for(long icol=0; icol<A->ny; icol++){
			for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
				y[A->pi[ix]]+=alpha*(A->px[ix]*creal(x[icol]));
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
void X(spmv)(X(mat)* y, const X(sp)* A, const X(mat)* restrict x, char trans, T alpha){
	if(A&&x&&y){
		assert(y);
		if(trans=='n'){
			if(PN(x)!=NY(A)||PN(y)!=NX(A)){
				error("Invalid dimension in y+=op(A)*x*alpha with op=%c: y is %ldx%ld, A is %ldx%ld, x is %ldx%ld", 
				trans, NX(y),NY(y),NX(A),NY(A),NX(x),NY(x));
			}
			for(long icol=0; icol<A->ny; icol++){
				for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					P(y,A->pi[ix])+=alpha*A->px[ix]*P(x,icol);
				}
			}
		} else if(trans=='t'){
			if(PN(x)!=NX(A)||PN(y)!=NY(A)){
				error("Invalid dimension in y+=op(A)*x*alpha with op=%c: y is %ldx%ld, A is %ldx%ld, x is %ldx%ld",
				trans, NX(y), NY(y), NX(A), NY(A), NX(x), NY(x));
			}
OMP_FOR(4)
			for(long icol=0L; icol<A->ny; icol++){			
				T tmp=0;
				for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					tmp+=alpha*(A->px[ix])*P(x,A->pi[ix]);
				}
				P(y,icol)+=tmp;
			}
		} else if(trans=='c'){
			if(PN(x)!=NX(A)||PN(y)!=NY(A)){
				error("Invalid dimension in y+=op(A)*x*alpha with op=%c: y is %ldx%ld, A is %ldx%ld, x is %ldx%ld",
				trans, NX(y), NY(y), NX(A), NY(A), NX(x), NY(x));
			}
OMP_FOR(4)			
			for(long icol=0L; icol<A->ny; icol++){
				T tmp=0;
				for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					tmp+=alpha*conj(A->px[ix])*P(x,A->pi[ix]);
				}
				P(y,icol)+=tmp;
			}
		} else{
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
static void X(spmm_do)(X(mat)** yout, const X(sp)* A, const X(mat)* x, const char trans[2], const int transy, const T alpha){
	if(!A||!x) return;
	mm_t D=parse_trans(A, x, trans);
	X(init)(yout, transy?D.ny:D.nx, transy?D.nx:D.ny);	
	X(mat)* y=*yout;
	if(x->ny==1&&trans[1]=='n'&&transy==0){
		X(spmv)(y, A, x, trans[0], alpha);
	} else{

#define no_conj(A) (A)
#define do_conj(A) conj(A)
#define no_trans(A,i,j) P(A,i,j)
#define do_trans(A,i,j) P(A,j,i)
#define LOOP_NORMA(ppy, yny,  conjA, ppx, conjx)				\
	for(long icol=0; icol<A->ny; icol++){				\
	    for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){		\
		for(long jcol=0; jcol<yny; jcol++){			\
		    ppy(y, A->pi[ix], jcol)+=alpha*conjA(A->px[ix])*conjx(ppx(x, icol, jcol)); \
		}							\
	    }								\
	}

#define LOOP_TRANSA(ppy, yny, conjA, ppx, conjx)				\
	{								\
		OMP_FOR(4)\
		for(long icol=0; icol<A->ny; icol++){\
			for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){		\
				for(long jcol=0; jcol<yny; jcol++){			\
				ppy(y, icol, jcol)+=alpha*conjA(A->px[ix])*conjx(ppx(x, A->pi[ix], jcol)); \
				}							\
			}							\
	    }								\
	}
		if(transy){
			if(trans[0]=='n'){
				switch(trans[1]){
				case 'n':
					LOOP_NORMA(do_trans, y->nx, no_conj, no_trans, no_conj);break;
				case 'N'://conjuate, but not transpose
					LOOP_NORMA(do_trans, y->nx, no_conj, no_trans, do_conj);break;
				case 't':
					LOOP_NORMA(do_trans, y->nx, no_conj, do_trans, no_conj);break;
				case 'c':
					LOOP_NORMA(do_trans, y->nx, no_conj, do_trans, do_conj);break;
				}
			} else if(trans[0]=='N'){//conjuate, but not transpose
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
			} else if(trans[0]=='t'){
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
			} else if(trans[0]=='c'){
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
			} else{
				error("Invalid trans=%s\n", trans);
			}
		} else{
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
			} else if(trans[0]=='N'){//conjuate, but not transpose
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
			} else if(trans[0]=='t'){
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
			} else if(trans[0]=='c'){
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
			} else{
				error("Invalid trans=%s\n", trans);
			}
		}
	}
}
/**
   sparse matrix multiply with dense matrix
*/
void X(spmm)(X(mat)** yout, const X(sp)* A, const X(mat)* x, const char trans[2], const T alpha){
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
void X(mulsp)(X(mat)** yout, const X(mat)* x, const X(sp)* A, const char trans[2], const T alpha){
	char trans2[2];
	trans2[1]=reverse_trans(trans[0]);
	trans2[0]=reverse_trans(trans[1]);
	X(spmm_do)(yout, A, x, trans2, 1, alpha);
}

/**
   Multiply two sparse arrays and add to the third: C0=C0+op(A)*op(B)*scale
*/
void X(spmulsp2)(X(sp)** C0, const X(sp)* A, const X(sp)* B, const char trans[2], const T scale){
	/*return C=C+ alpha*(A*B) */
	if(!A||!B) return;
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
	X(sp)* res=X(ss_multiply)(At?At:A, Bt?Bt:B);
	X(spsort)(res);
	if(fabs(scale-(T)1.)>EPS){
		X(spscale)(res, scale);
	}
	if(!*C0)
		*C0=res;
	else{
		X(spadd)(C0, 1, res, 1);
		X(spfree)(res);
	}
	X(spfree)(At);
	X(spfree)(Bt);
	X(ss_dropzeros)(*C0);
}
/**
   Multiply two sparse arrays and return the result op(A)*op(B)
*/
X(sp)* X(spmulsp)(const X(sp)* A, const X(sp)* B, const char trans[2]){
	X(sp)* res=0;
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

   2009-11-09: There was initially a beta parameter It was implemented wrongly
   for beta!=1 because for every call to dmm, the already accumulated ones are
   scaled.  removed beta.
*/
void X(cellmm)(panyarray C0_, const_anyarray A_, const_anyarray B_, const char trans[2], const R alpha){
	cell **C0=C0_.c;
	const cell *A=A_.c;
	const cell *B=B_.c;
	if(!A||!B) return;
	if(iscell(A)&&iscell(B)){//multiplication of cells.
		mm_t D=parse_trans(A, B, trans);
		cellinit(C0, D.nx, D.ny);
		cell* C=*C0;
OMP_FOR_COLLAPSE(2, 8)
		for(int iy=0; iy<D.ny; iy++){
			for(int ix=0; ix<D.nx; ix++){
				for(int iz=0; iz<D.nz; iz++){
					X(cellmm)(P(C)+ix+iy*D.nx, A->p[ix*D.ax+iz*D.az], B->p[iz*D.bz+iy*D.by], trans, alpha);
				}
			}
		}
	} else{//multiplication between dense or sparse 
		if(ismat(B)){//dense B
			if(ismat(A)){ //dense A
				X(mm)((X(mat)**)C0, 1, (X(mat)*)A, (X(mat)*)B, trans, alpha);
			} else if(issp(A)){//sparse A
				X(spmm)((X(mat)**)C0, (X(sp)*)A, (X(mat)*)B, trans, alpha);
			} else{
				error("Invalid A type, id=%u.\n", A->id);
			}
		} else if(issp(B)){//sparse B
			if(ismat(A)){//dense A
				X(mulsp)((X(mat)**)C0, (X(mat)*)A, (X(sp)*)B, trans, alpha);
			} else if(issp(A)){//sparse A
				if(!*C0||issp(*C0)){//result is sparse
					X(spmulsp2)((X(sp)**)C0, (X(sp)*)A, (X(sp)*)B, trans, alpha);
				} else if(ismat(*C0)){//result is dense
					X(mat)* Afull=0; X(spfull)(&Afull, (X(sp)*)A, 'n', 1);
					X(mat)* Bfull=0; X(spfull)(&Bfull, (X(sp)*)B, 'n', 1);
					X(mm)((X(mat)**)C0, 1, Afull, Bfull, trans, alpha);
					X(free)(Afull);
					X(free)(Bfull);
				} else{
					error("Invalid C type, id=%u. \n", (*C0)->id);
				}
			} else{
				error("Invalid A type, id=%u.\n", A->id);
			}
		} else{
			error("Invalid B type, id=%u.\n", B->id);
		}
	}
}
void X(mm_cell)(X(mat)** C0, const cell* A, const X(mat)* B, const char trans[2], const R alpha){
	if(!A||!B||!C0) return;
	X(cellmm)((cell**)C0, A, B, trans, alpha);
}
/*void X(cellmm)(X(cell)** C0, const X(cell)* A, const X(cell)* B, const char trans[2], const R alpha){
	if(!A || !B || !C0) return;
	X(cellmm)((cell**)C0, A, B, trans, alpha);
}*/
/**
   a different interface for multiplying cells.
 */
X(cell)* X(cellmm2)(const X(cell)* A, const X(cell)* B, const char trans[2]){
	X(cell)* res=0;
	X(cellmm)(&res, A, B, trans, 1);
	return res;
}

/**
   Add alpha to diagnonal elements of A_. A_ can be matrix or sparse matrix.
 */
void X(celladdI)(anyarray A_, T alpha){
	cell* A=A_.c;
	if(!A) return;
	assert(A->nx==A->ny);
OMP_FOR(8)
	for(int ii=0; ii<A->ny; ii++){
		if(!P(A, ii, ii)){
		} else if(ismat(P(A, ii, ii))){
			X(addI)(X(mat_cast)(P(A, ii, ii)), alpha);
		} else if(issp(P(A, ii, ii))){
			X(spaddI)(X(sp_cast)(P(A, ii, ii)), alpha);
		} else{
			error("Invalid id=%u", P(A, ii, ii)->id);
		}
	}
}

/**
   Calculate A_=A_*ac+B_*bc;

   Takes parameters of matrix, sparse matrix, or cell array of them.
 */
void X(celladd)(panyarray pA_, R ac, const_anyarray B_, R bc){
	const cell* B=B_.c;
	cell** pA=pA_.c;
	if(!pA){
		error("pA is null\n");
		return; 
	}
	if(!B){
		if(ac==1){
			return;//no operation
		}
		X(cellscale)(*pA, ac);
	}
	
	if(iscell(B)){//cell
		cellinit2(pA, B);
		cell* A=*pA;
OMP_FOR(8)
		for(int i=0; i<B->nx*B->ny; i++){
			X(celladd)(&P(A,i), ac, P(B,i), bc);
		}
	} else{//non cell
		if(!*pA||ismat(*pA)){//A is dense
			if(ismat(B)){//Add dense to dense
				X(add)((X(mat)**)pA, ac, (X(mat)*)B, bc);
			} else{//add sparse to dense
				if(ac!=1){
					X(scale)((X(mat*))*pA, ac);
				}
				X(spfull)((X(mat)**)(pA), (X(sp)*)B, 'n', bc);
			}
		} else if(issp(*pA)){
			if(issp(B)){//add sparse to sparse
				X(sp)* tmp=X(spadd2)((X(sp)*)(*pA), ac, (X(sp)*)B, bc);
				X(spmove)((X(sp)*)(*pA), (X(sp)*)tmp);
				free(tmp);
			} else{
				error("Adding dense to sparse matrix is not supported.\n");
			}
		} else{
			error("Invalid operand\n");
		}
	}
}
/**
   Copy B to A.

   Takes parameters of matrix, sparse matrix, or cell array of them.
 */
void X(cellcp)(X(cell)** A_, const X(cell)* B_){
	if(!B_){
		if(A_ && *A_) X(cellzero)(*A_);
	} else{
		X(celladd)(A_, 0, B_, 1);
	}
}

/**
   scale each element of A.
*/
void X(cellscale)(anyarray A_, R w){
	cell* A=A_.c;
	if(!A) {
		return;
	}else if(iscell(A)){
OMP_FOR(8)
		for(int i=0; i<A->nx*A->ny; i++){
			X(cellscale)(P(A,i), w);
		}
	} else{
		if(ismat(A)){
			X(scale)((X(mat*))A, w);
		} else if(issp(A)){
			X(spscale)((X(sp*))A, w);
		} else if(A){
			error("Invalid type: id=%u\n", A->id);
		}
	}
}

/**
   Setting all elements of a cell to zero.
*/
void X(cellzero)(X(cell)* dc){
	if(dc) X(cellscale)(dc, 0);
}

/**
   Convert dense or sparse matrix cell to matrix.
*/
X(mat)* X(cell2m)(const_anyarray A_){
	const cell* A=A_.c;
	if(!A) return 0;
	long nx=0, ny=0, *nxs=NULL, *nys=NULL;
	celldim(A, &nx, &ny, &nxs, &nys);
	
	X(mat)* out=X(new)(nx, ny);
	long jcol=0;
	for(long iy=0; iy<A->ny; iy++){
		for(long icol=0; icol<nys[iy]; icol++){
			long kr=0;
			for(long ix=0; ix<A->nx; ix++){
				if(!isempty(P(A, ix, iy))){
					T* pout=&P(out, kr, icol+jcol);
					if(ismat(P(A, ix, iy))){
						memcpy(pout, PCOL(((X(mat*))P(A, ix, iy)), icol), nxs[ix]*sizeof(T));
					} else if(issp(P(A, ix, iy))){//convert sparse col to full
						X(sp*)Asp=(X(sp*))P(A, ix, iy);
						for(long j=Asp->pp[icol]; j<Asp->pp[icol+1]; j++){
							pout[Asp->pi[j]]=Asp->px[j];
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

/**
 * Rotate a or multiple 2d array ccw by angle
 * */
void X(maprot)(anyarray A_, real theta){
	const cell* A=A_.c; 
	if(!A|| !theta) return;
	if(iscell(A)){
		for(int i=0; i<PN(A); i++){
			X(maprot)(P(A,i),theta);
		}
	}else if(ismat(A)){
		X(mat)*Am=X(mat_cast)(A);
		X(mat)*B=X(dup)(Am);
		X(zero)(Am);
		X(embed)(Am, B, theta);
		X(free)(B);
	}else{
		error("Invalid usage\n");
	}
}
