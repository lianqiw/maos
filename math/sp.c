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




#include <search.h>
#include "../sys/sys.h" 
#include "mathmisc.h"
#include "type.h"
#include "mathdef.h"
#include "defs.h"
#include "suitesparse.c"
#define assert_sp(A) assert(!A || A->id==M_SPT)
/**
   Create a nx*ny sparse matrix with memory for nmax max
   elements allocated.
*/
X(sp)* X(spnew)(long nx, long ny, long nzmax){
	if(nx<0) nx=0;
	if(ny<0) ny=0;
	X(sp)* sp=mycalloc(1, X(sp));
	sp->id=M_SPT;
	sp->pp=mymalloc((ny+1), spint);
	if(nzmax>0){
		sp->pi=mymalloc(nzmax, spint);
		sp->px=mymalloc(nzmax, T);
	}
	sp->nx=nx;
	sp->ny=ny;
	sp->nzmax=nzmax;
	sp->nref=mycalloc(1, unsigned int);
	sp->nref[0]=1;	
	return sp;
}
static void X(spfree_content)(X(sp)* sp){
	if(!sp) return;
	assert(issp(sp));
	if(sp->fp){
		writedata_by_id(sp->fp, sp, 0, 0);
		zfclose(sp->fp);
		sp->fp=NULL;
	}
	if(sp->nref && !atomic_sub_fetch(sp->nref, 1)){
		free(sp->px);
		free(sp->pp);
		free(sp->pi);
		free(sp->nref);
	}
}
/**
 * free a sparse matrix*/
void X(spfree_do)(X(sp)* sp){
	if(sp){
		X(spfree_content)(sp);
		free(sp);
	}
}

/**
   reference a sparse object.
*/
X(sp)* X(spref)(X(sp)* A){
	if(!A) return NULL;
	assert_sp(A);
	X(sp)* out=mycalloc(1, X(sp));
	if(!A->nref){
		warning_once("Referencing not owned data. This may cause error.\n");
	} else{
		atomic_add_fetch(A->nref,1);
	}
	memcpy(out, A, sizeof(X(sp)));
	return out;
}
/**
   move the matrix from res to A without change value of A itself.
*/
void X(spmove)(X(sp)* A, X(sp)* res){
	if(res&&A){
		assert(issp(res)&&issp(A));
		X(spfree_content)(A);
		memcpy(A, res, sizeof(X(sp)));
		memset(res, 0, sizeof(X(sp)));
	} else{
		if(!(!res&&!A)){
			error("Trying to move an NULL matrix\n");
		}
	}
}

/**
   copy a sparse matrix to another.
*/
X(sp)* X(spdup)(const X(sp)* A){
	if(!A) return NULL;
	assert_sp(A);
	long nmax=A->pp[A->ny];
	X(sp)* out;
	out=X(spnew)(A->nx, A->ny, nmax);
	memcpy(out->pp, A->pp, sizeof(spint)*(A->ny+1));
	if(nmax){
		memcpy(out->pi, A->pi, sizeof(spint)*nmax);
		memcpy(out->px, A->px, sizeof(T)*nmax);
	}
	return out;
}

/**
   Create a new sparse matrix of the same size as A.
*/
X(sp)* X(spnew2)(const X(sp)* A){
	if(!A) return NULL;
	assert_sp(A);
	return X(spnew)(A->nx, A->ny, A->pp[A->ny]);
}
/**
   Create a new sparse matrix and fill in uniform random
   numbers with filling factor of 'fill'
*/
X(sp)* X(spnewrandu)(long nx, long ny, const T mean, R fill, rand_t* rstat){
	if(fill>1) fill=1.;
	if(fill<0) fill=0.;
	const long nzmax=nx*ny;
	long nz1=nx*ny*fill*4;
	if(nz1>nzmax) nz1=nzmax;
	X(sp)* A=X(spnew)(nx, ny, nz1);
	spint* pp=A->pp;
	spint* pi=A->pi;
	T* px=A->px;
	long count=0;
	R thres=1.-fill;
	for(int icol=0; icol<A->ny; icol++){
		pp[icol]=count;
		for(int irow=0; irow<A->nx; irow++){
			if(randu(rstat)>thres){
				pi[count]=irow;
				px[count]=RANDU(rstat)*mean;
				count++;
				if(count>nz1){
					/*check out of bound; */
					nz1=nz1*2; if(nz1>nzmax) nz1=nzmax;
					X(spsetnzmax)(A, nz1);
					/*the pointers may change */
					pp=A->pp;
					pi=A->pi;
					px=A->px;
				}
			}
		}
	}
	pp[A->ny]=count;
	X(spsetnzmax)(A, count);
	return A;
}
/**
   Cast a pointer to sparse array.
 */
X(sp)* X(sp_cast)(const cell* A){
	return (A&&A->nx&&A->ny&&issp(A))?(X(sp)*)A:NULL;
}
/**
   resize a sparse matrix
*/
void X(spsetnzmax)(X(sp)* sp, long nzmax){
	if(!sp) return;
	if(!issp(sp)){
		warning("Not a sparse matrix\n");
		return;
	}
	if(sp->nzmax!=nzmax){
		sp->pi=myrealloc(sp->pi, nzmax, spint);
		sp->px=myrealloc(sp->px, nzmax, T);
		sp->nzmax=nzmax;
	}
}

/**
 * Display a sparse array*/
void X(spdisp)(const X(sp)* sp){
	long ic, ir;
	long imax;
	assert(issp(sp));
	if(sp->nzmax==0){
		dbg("X(spdisp): All zeros\n");
	} else{
		dbg("X(spdisp):\n");
		for(ic=0; ic<sp->ny; ic++){
			imax=-1;
			for(ir=sp->pp[ic];ir<sp->pp[ic+1];ir++){
#ifdef COMP_COMPLEX
				printf("(%ld,%ld)=(%g,%g)\n",
					(long)sp->pi[ir], (long)ic, REAL(sp->px[ir]), IMAG(sp->px[ir]));
#else		
				printf("(%ld,%ld)=%g\n", (long)sp->pi[ir], (long)ic, sp->px[ir]);
#endif
				if(sp->pi[ir]>imax){
					imax=sp->pi[ir];
				} else{
					warning("Wrong order");
				}
			}
		}
	}
}
/**
 * Check a sparse array for wrong orders. Return 1 if is lower triangle, 2 if
 * upper triangle, and 3 if diagonal.*/
int X(spcheck)(const X(sp)* sp){
	assert(issp(sp));
	int not_lower=0;
	int not_upper=0;
	if(sp){
		long ic, ir;
		long imax;
		for(ic=0; ic<sp->ny; ic++){
			imax=-1;
			if(sp->pp[ic+1]<sp->pp[ic]){
				error("p in column %ld is smaller than %ld\n", ic+1, ic);
			}
			for(ir=sp->pp[ic];ir<sp->pp[ic+1];ir++){
				if(sp->pi[ir]>imax){
					imax=sp->pi[ir];
				} else{
					warning("Wrong order at column %ld", ic);
				}
				if(sp->pi[ir]<ic) not_lower=1;
				if(sp->pi[ir]>ic) not_upper=1;
			}
			if(imax>=sp->nx){
				error("imax=%ld exceeds column size at column %ld\n", imax, ic);
			}
		}
		if(sp->pp[sp->ny]!=sp->nzmax){
			warning("real nzmax is %ld, allocated is %ld\n", (long)sp->pp[sp->ny], sp->nzmax);
		}
	}
	return (not_lower?0:1)|(not_upper?0:2);
}
/**
 * inplace scale sparse matrix elements.*/
void X(spscale)(X(sp)* A, const T beta){
	if(!A) return;
	assert_sp(A);
	if(A->nref[0]>1){
		warning("spscale on referenced dsp\n");
	}
	for(long i=0; i<A->pp[A->ny]; i++){
		A->px[i]*=beta;
	}
}

/**
 * inplace scale sparse matrix elements.*/
void X(spscalex)(X(sp)* A, const X(mat)* xs){
	if(!A) return;
	assert_sp(A);
	if(A->nx!=xs->nx){
		error("sparse matrix[%ldx%ld] does not match vector [%ldx%ld]\n",
			A->nx, A->ny, xs->nx, xs->ny);
	}
	for(long iy=0; iy<A->ny; iy++){
		for(long i=A->pp[iy]; i<A->pp[iy+1]; i++){
			long ix=A->pi[i];
			A->px[i]*=P(xs,ix);
		}
	}
}

/**
 * inplace scale sparse matrix elements.*/
void X(spscaley)(X(sp)* A, const X(mat)* ys){
	if(!A) return;
	assert_sp(A);
	if(A->ny!=ys->nx){
		error("sparse matrix[%ldx%ld] does not match vector [%ldx%ld]\n",
			A->nx, A->ny, ys->nx, ys->ny);
	}
	for(long iy=0; iy<A->ny; iy++){
		for(long i=A->pp[iy]; i<A->pp[iy+1]; i++){
			A->px[i]*=P(ys,iy);
		}
	}
}
/**
   Reference a spcell array.
 */
X(spcell)* X(spcellref)(const X(spcell)* A){
	if(!A) return NULL;
	X(spcell)* out=(X(spcell*))cellnew(A->nx, A->ny);
	for(long i=0; i<A->nx*A->ny; i++){
		P(out, i)=X(spref)(P(A, i));
	}
	return out;
}
/**
   cast a cell object to X(spcell) after checking
 */
X(spcell)* X(spcell_cast)(const cell* A){
	if(!A) return 0;
	int err=0;
	for(int i=0; i<A->nx*A->ny; i++){
		if(P(A,i)&&!issp(P(A,i))){
			warning("A[%d] is not sparse\n", i);
			err++;
		}
	}
	if(err){
		print_backtrace();
	}
	return (X(spcell)*)A;
}
/**
 * inplace scale a X(dspcell) object*/
void X(spcellscale)(X(spcell)* A, const T beta){
	for(int i=0; i<A->nx*A->ny; i++){
		X(spscale)(P(A,i), beta);
	}
}
/**
 * Create a new sparse matrix with diagonal elements set to vec*alpha*/
X(sp)* X(spnewdiag)(long N, T* vec, T alpha){
	X(sp)* out=X(spnew)(N, N, N);
	spint* pp=out->pp;
	spint* pi=out->pi;
	T* px=out->px;
	long count=0;
	if(vec){
		for(long icol=0; icol<out->ny; icol++){
			pp[icol]=count;
			pi[count]=icol;
			px[count]=vec[icol]*alpha;
			count++;
		}
	} else{
		for(long icol=0; icol<out->ny; icol++){
			pp[icol]=count;
			pi[count]=icol;
			px[count]=alpha;
			count++;
		}
	}
	pp[out->ny]=count;
	return out;
}
/**
   Extract diagonal element of A and return
*/
X(mat)* X(spdiag)(const X(sp)* A){
	if(!A) return NULL;
	if(A->nx!=A->ny){
		error("Only implemented for square matrix\n");
	}
	assert_sp(A);
	X(mat)* out=X(new)(A->nx, 1);
	for(long icol=0; icol<A->ny; icol++){
		for(long irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
			long row=A->pi[irow];
			if(row==icol){
				P(out,icol)=A->px[irow];
			}
		}
	}
	return out;
}
/**
   Multiply a sparse matrix inplace with a diagonal weighting matrix whose
   diagonal values are stored in w.
   W_ii=w_i; W_ij=0 if i!=j
   A=A*W*alpha;
   W is a diagonal sparse matrix. diag(W) is w
   multiply w[i] to all numbers in column[i]
*/
void X(spmuldiag)(X(sp)* restrict A, const T* w, T alpha){
	if(A&&w){
		assert_sp(A);
		for(long icol=0; icol<A->ny; icol++){
			const T wi=w[icol]*alpha;
			for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
				A->px[ix]*=wi;
			}
		}
	}
}

/**
 * Multiply two vectors with weighting by sparse matrix.  return y'*(A*x)*/
T X(spwdot)(const X(mat)* y, const X(sp)* A, const X(mat)* x){
	/*X(sp) weighted ddot. */
	/*computes y'*(A*x). x,y are vectors */
	T res=0;
	if(x&&y){
		if(A){
			if(!(issp(A) &&x->ny==1&&y->ny==1&&A->nx==y->nx&&A->ny==x->nx)){
				error("Parameters mismatch\n");
			}
			for(int icol=0; icol<A->ny; icol++){
				for(int ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					res+=P(y, A->pi[ix])*A->px[ix]*P(x, icol);
				}
			}
		} else{
			res=X(dot)(x, y);
		}
	}
	return res;
}
/**
 * Multiply two cell arrays with weighting by sparse matrix*/
T X(spcellwdot)(const X(cell)* y, const X(spcell)* A, const X(cell)* x){
	/*computes y'*(A*x) */
	T res=0;
	if(x&&y){
		if(A){
			if(!(issp(A)&&x->ny==1&&y->ny==1&&A->nx==y->nx&&A->ny==x->nx)){
				error("Parameters mismatch\n");
			}
			/*PSX(cell)* Ap=A; */
			for(int iy=0; iy<A->ny; iy++){
				for(int ix=0; ix<A->nx; ix++){
					res+=X(spwdot)(P(y, ix), P(A, ix, iy), P(x, iy));
				}
			}
		} else{
			res=X(celldot)(x, y);
		}
	}
	return res;
}

/**
 * Convert sparse matrix into dense matrix and add to output:
 * out0=out0+full(A)*alpha*/
void X(spfull)(X(mat)** out0, const X(sp)* A, const char trans, const T alpha){
	if(!A) return;
		/**
		   add A*f to dense matrix located in p;
		*/
	if(trans=='n'){
		if(issp(A)){//sparse
			long nx=A->nx;
			long icol, ix, irow;
			if(!*out0){
				*out0=X(new)(A->nx, A->ny);
			}
			X(mat)* out=*out0;
			assert(out->nx==A->nx&&out->ny==A->ny);
			for(icol=0; icol<A->ny; icol++){
				for(ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					irow=A->pi[ix];
					if(irow>=nx)
						error("invalid row:%ld, %ld", irow, nx);
					P(out, irow, icol)+=alpha*A->px[ix];
				}
			}
		} else{
			X(add)(out0, 1, (X(mat)*)A, alpha);
		}
	} else if(trans=='t'){
		if(issp(A)){
			long nx=A->nx;
			long icol, ix, irow;
			if(!*out0){
				*out0=X(new)(A->ny, A->nx);
			}
			X(mat)* out=*out0;
			assert(out->nx==A->ny&&out->ny==A->nx);
			for(icol=0; icol<A->ny; icol++){
				for(ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
					irow=A->pi[ix];
					if(irow>=nx)
						error("invalid row:%ld, %ld", irow, nx);
					P(out, icol, irow)+=alpha*A->px[ix];
				}
			}
		} else{
			X(mat)* tmp=X(trans)((X(mat)*)A);
			X(add)(out0, 1, tmp, alpha);
			X(free)(tmp);
		}
	} else{
		error("trans=%c is invalid\n", trans);
	}
}
/**
   Convert dense matrix A to sparse matrix with threshold.
 */
X(sp)* X(2sp)(X(mat)* A, R thres){
	if(!A) return 0;
	assert(ismat(A));
	X(sp)* out=X(spnew)(A->nx, A->ny, A->nx*A->ny);
	long count=0;
	for(long icol=0; icol<A->ny; icol++){
		out->pp[icol]=count;
		for(long irow=0; irow<A->nx; irow++){
			T val=P(A, irow, icol);
			if(fabs(val)>thres){
				out->pi[count]=irow;
				out->px[count]=val;
				count++;
			}
		}
	}
	out->pp[A->ny]=count;
	X(spsetnzmax)(out, count);
	return out;
}

/**
 * Added two sparse matrices: return A*a+B*b*/
X(sp)* X(spadd2)(const X(sp)* A, T a, const X(sp)* B, T b){
	assert(issp(A)&&issp(B));
	if(A->nx!=B->nx||A->ny!=B->ny){
		error("X(sp) matrix mismatch: (%ldx%ld) vs (%ldx%ld\n",
			A->nx, A->ny, B->nx, B->ny);
	}
	X(sp)* C=X(ss_add)(A, B, a, b);
	if(!C){
		warning("ss_add returned null\n");
	}
	X(ss_dropzeros)(C);
	return C;
}
/**
 * Add a sparse matrix to another: A0=A0+B*/
void X(spadd)(X(sp)** A0, T alpha, const X(sp)* B, T beta){
	/*add B to A. */
	if(B){
		if(!*A0){
			*A0=X(spdup)(B);
			if(beta!=(T)1){
				X(spscale)(*A0, beta);
			}
		} else{
			X(sp)* res=X(spadd2)(*A0, alpha, B, beta);
			X(spmove)(*A0, res);
			free(res);
		}
	}
}

/**
   Add alpha times identity to a sparse matrix.
   If X(sp) is not symmetric, only add diagonal to first nm*nm block for nm=min(nx,ny)
*/
void X(spaddI)(X(sp)* A, T alpha){
	X(spsort)(A);//make sure it is sorted correctly.
	//First check for missing diagonal elements
	long nm=MIN(A->nx, A->ny);
	long missing=0;
	for(long icol=0; icol<nm; icol++){
		int found=0;
		for(long ix=A->pp[icol]; ix<A->pp[icol+1]; ix++){
			if(A->pi[ix]==icol){
				found=1;
				break;
			}
		}
		if(!found) missing++;
	}
	long nzmax=A->pp[A->ny];
	if(missing){//expanding storage
		A->px=myrealloc(A->px, (nzmax+missing), T);
		A->pi=myrealloc(A->pi, (nzmax+missing), spint);
	}
	missing=0;
	for(long icol=0; icol<A->ny; icol++){
		A->pp[icol]+=missing;
		int found=0;
		long ix=0;
		if(icol<nm){
			for(ix=A->pp[icol]; ix<A->pp[icol+1]+missing; ix++){
				if(A->pi[ix]==icol){
					found=1;
					A->px[ix]+=alpha;
					break;
				} else if(A->pi[ix]>icol){ //insertion place
					break;
				}
			}
			if(!found){
				memmove(A->px+ix+1, A->px+ix, sizeof(T)*(nzmax+missing-ix));
				memmove(A->pi+ix+1, A->pi+ix, sizeof(spint)*(nzmax+missing-ix));
				A->pi[ix]=icol;
				A->px[ix]=alpha;
				missing++;
			}
		}
	}
	A->pp[A->ny]+=missing;
	A->nzmax=A->pp[A->ny];
}
/**
 * Transpose a sparse array*/
X(sp)* X(sptrans)(const X(sp)* A){
	if(!A) return NULL;
	assert_sp(A);
	X(sp)* res=X(ss_transpose)(A, 1);
	X(ss_dropzeros)(res);
	return res;
}
/**
   Take conjugation elementwise.
 */
void X(spconj)(X(sp)* A){
#ifdef COMP_COMPLEX
	const long nzmax=A->pp[A->ny];
	for(long i=0; i<nzmax; i++){
		A->px[i]=conj(A->px[i]);
	}
#else
	(void)A;
#endif
}
/**
   Expand sparse matrix cell to dense matrix cell.
 */
void X(spcellfull)(X(cell)** out0, const X(spcell)* A, const char trans, const T alpha){
	if(!A) return;
	if(!*out0){
		if(trans=='n'){
			*out0=X(cellnew)(A->nx, A->ny);
		} else if(trans=='t'){
			*out0=X(cellnew)(A->ny, A->nx);
		} else{
			error("trans=%c is invalid\n", trans);
		}
	}
	for(int iy=0; iy<A->ny; iy++){
		for(int ix=0; ix<A->nx; ix++){
			if(trans=='n'){
				X(spfull)(&P((*out0), ix, iy), P(A, ix, iy), trans, alpha);
			} else{
				X(spfull)(&P((*out0), iy, ix), P(A, ix, iy), trans, alpha);
			}
		}
	}
}
/**
 * Transpose a sparse cell*/
X(spcell)* X(spcelltrans)(const X(spcell)* spc){
	if(!spc) return NULL;
	long nx, ny;
	nx=spc->nx;
	ny=spc->ny;
	X(spcell)* spct=(X(spcell*))cellnew(ny, nx);

	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(spct, iy, ix)=X(sptrans)(P(spc, ix, iy));
		}
	}
	return spct;
}

/**
 * Concatenate two sparse array along dim dimension*/
X(sp)* X(spcat)(const X(sp)* A, const X(sp)* B, int dim){
	X(sp)* C=NULL;
	if(!B){
		C=X(spdup)(A);
	} else if(!A){
		C=X(spdup)(B);
	} else if(dim==0){
		error("Not implemented\n");
		/*
		  |A|
		  |B|
		*/
	} else if(dim==1){
	/*|AB|*/
		assert(issp(A)&&issp(B));
		if(A->nx!=B->nx){
			error("X(sp) matrix doesn't match\n");
		}
		const long nzmax=A->pp[A->ny]+B->pp[B->ny];
		if(nzmax==0){
			return 0;
		}
		C=X(spnew)(A->nx, A->ny+B->ny, nzmax);
		memcpy(C->pp, A->pp, A->ny*sizeof(spint));
		memcpy(C->pi, A->pi, A->pp[A->ny]*sizeof(spint));
		memcpy(C->px, A->px, A->pp[A->ny]*sizeof(T));
		memcpy(C->pi+A->pp[A->ny], B->pi, B->pp[B->ny]*sizeof(spint));
		memcpy(C->px+A->pp[A->ny], B->px, B->pp[B->ny]*sizeof(T));
		const long Anzmax=A->pp[A->ny];
		for(long i=0; i<B->ny+1; i++){
			C->pp[i+A->ny]=Anzmax+B->pp[i];
		}
	} else{
		error("Wrong dimension\n");
	}
	return C;
}
/**
 * Concatenate a dspcell to sparse array*/
X(sp)* X(spcell2sp)(const X(spcell)* A){
	/*convert X(spcell) to sparse. */
	if(A->nx*A->ny==1){/*There is a single cell */
		return X(spref)(P(A, 0));
	}
	long nx=0, ny=0, nzmax=0;
	long nnx[A->nx];//number of rows of each row block
	long nny[A->ny];//number of cols for each column block
	for(long ix=0; ix<A->nx; ix++){
		nnx[ix]=0;
		for(long iy=0; iy<A->ny; iy++){
			if(P(A, ix, iy)){
				if(!issp(P(A, ix, iy))){
					error("Not every cell is sparse\n");
				}
				if(!nnx[ix]){
					nnx[ix]=P(A, ix, iy)->nx;
					nx+=P(A, ix, iy)->nx;
				} else if(nnx[ix]!=P(A, ix, iy)->nx){
					error("block (%ld, %ld) has wrong number of rows. expect %ld got %ld\n",
					ix, iy, nnx[ix], P(A, ix, iy)->nx);
				}
			}
		}
	}
	for(long iy=0; iy<A->ny; iy++){
		nny[iy]=0;
		for(long ix=0; ix<A->nx; ix++){
			if(P(A, ix, iy)){
				if(!nny[iy]){
					nny[iy]=P(A, ix, iy)->ny;
					ny+=P(A, ix, iy)->ny;
				} else if(nny[iy]!=P(A, ix, iy)->ny){
					error("block (%ld, %ld) has wrong number of cols. expect %ld got %ld\n",
					ix, iy, nny[iy], P(A, ix, iy)->ny);
				}
			}
		}
	}
	for(long i=0; i<A->nx*A->ny; i++){
		if(P(A,i)){
			nzmax+=P(A,i)->pp[P(A,i)->ny];
		}
	}
	X(sp)* out=X(spnew)(nx, ny, nzmax);
	long count=0;
	long jcol=0;
	for(long iy=0; iy<A->ny; iy++){
		for(long icol=0; icol<nny[iy]; icol++){
			out->pp[jcol+icol]=count;
			long kr=0;
			for(long ix=0; ix<A->nx; ix++){
				if(P(A, ix, iy)){
					for(long ir=P(A, ix, iy)->pp[icol];
						ir<P(A, ix, iy)->pp[icol+1]; ir++){
						out->px[count]=P(A, ix, iy)->px[ir];
						out->pi[count]=P(A, ix, iy)->pi[ir]+kr;
						count++;
					}
				}
				kr+=nnx[ix];
			}
		}
		jcol+=nny[iy];
	}
	out->pp[ny]=count;
	if(count!=nzmax){
		error("X(spcell2sp) gets Wrong results. count=%ld, nzmax=%ld\n", count, nzmax);
	}
	/*nzmax maybe smaller than A->pp[A->n]  */
	/*because nzmax simply show the slots available. */
	return out;
}

/**
 * Sum elements of sparse array along dimension dim*/
X(mat)* X(spsum)(const X(sp)* A, int dim){
	if(!A) return 0;
	/*Sum sparse matrix along col or row to form a vector */
	X(mat)* v=NULL;
	assert_sp(A);
	T* p;
	switch(dim){
	case 1:/*sum along col */
		v=X(new)(1, A->ny);
		p=P(v);
		for(int icol=0; icol<A->ny; icol++){
			for(int irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				p[icol]+=A->px[irow];
			}
		}
		break;
	case 2:/*sum along row */
		v=X(new)(A->nx, 1);
		p=P(v);
		for(int icol=0; icol<A->ny; icol++){
			for(int irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				p[A->pi[irow]]+=A->px[irow];
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
X(mat)* X(spsumabs)(const X(sp)* A, int col){
	if(!A) return 0;
	assert_sp(A);
	X(mat)* v=NULL;
	T* p;
	switch(col){
	case 1:/*sum along col */
		v=X(new)(1, A->ny);
		p=P(v);
		for(int icol=0; icol<A->ny; icol++){
			for(int irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				p[icol]+=fabs(A->px[irow]);
			}
		}
		break;
	case 2:/*sum along row */
		v=X(new)(A->nx, 1);
		p=P(v);
		for(int icol=0; icol<A->ny; icol++){
			for(int irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				p[A->pi[irow]]+=fabs(A->px[irow]);
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
void X(spclean)(X(sp)* A){
	if(!A) return;
	assert_sp(A);
	X(ss_dropzeros)(A);
}

/**
   Drop elements that are EPS times the largest value.
*/
void X(spdroptol)(X(sp)* A, R thres){
	if(!A) return;
	assert_sp(A);
	if(thres<EPS) thres=EPS;
	R maxv;
	X(vecmaxmin)(A->px, A->nzmax, &maxv, NULL);
	X(ss_droptol)(A, maxv*thres);
}
/**
   Drop elements that are EPS times the largest value.
*/
void X(spcelldroptol)(X(spcell)* A, R thres){
	for(int i=0; i<A->nx*A->ny; i++){
		X(spdroptol)(P(A,i), thres);
	}
}
/**
 * Drop empty columns
*/
lmat *X(spdropemptycol)(X(sp)*A){
	if (!A) return NULL;
	lmat *pmap=lnew(A->ny, 1);
	long jcol=0;
	for(long icol=0; icol<A->ny; icol++){
		if(icol!=jcol){//there are skipped cols
			A->pp[jcol]=A->pp[icol];
		}
		P(pmap, jcol)=icol;
		if(A->pp[icol+1]>A->pp[icol]){//valid col
			jcol++;
		}
	}
	A->pp[jcol]=A->pp[A->ny];
	A->ny=jcol;
	lresize(pmap, jcol, 1);
	//no need to modify A->pp vector
	return pmap;
}
typedef struct spelem{
	int pi;
	T px;
}spelem;
static int spelemcmp(const spelem* A, const spelem* B){
	return A->pi-B->pi;
}
/**
   Make sure the elements are sorted correctly. Does not change the location of
data. can be done without harm. */
void X(spsort)(X(sp)* A){
	if(!A||A->ny==0||A->nx==0) return;
	assert_sp(A);
	long nelem_max=A->nzmax/A->ny*2;
	spelem* col=mymalloc(nelem_max, spelem);
	for(long i=0; i<A->ny; i++){
		long nelem=(A->pp[i+1]-A->pp[i]);
		if(nelem==0) continue;
		if(nelem>nelem_max){
			nelem_max=nelem;
			col=myrealloc(col, nelem_max, spelem);
		}
		for(long j=0; j<nelem; j++){
			col[j].pi=A->pi[A->pp[i]+j];
			col[j].px=A->px[A->pp[i]+j];
		}
		qsort(col, nelem, sizeof(spelem), (int(*)(const void*, const void*))spelemcmp);
		for(long j=0; j<nelem; j++){
			A->pi[A->pp[i]+j]=col[j].pi;
			A->px[A->pp[i]+j]=col[j].px;
		}
	}
	free(col);

}
/**
   Make sure the elements are sorted correctly.
*/
void X(spcellsort)(X(spcell)* A){
	for(int i=0; i<A->nx*A->ny; i++){
		X(spsort)(P(A,i));
	}
}
/**
   symmetricize a sparse matrix and drop values below a
   threshold.
*/
void X(spsym)(X(sp)** A){
	X(sp)* B=X(sptrans)(*A);
	X(spadd)(A, (T)1, B, (T)1);
	X(spscale)(*A, 0.5);
	X(spfree)(B);
	X(spdroptol)(*A, EPS);
	X(spsort)(*A);/*This is important to make chol work. */
}
/**
   symmetricize a sparse cell and drop values below a
   threshold.
*/
void X(spcellsym)(X(spcell)** A){
	X(spcell)* B=X(spcelltrans)(*A);
	X(spcelladd)(A, 1, B, 1);
	X(spcellfree)(B);
	X(spcellscale)(*A, 0.5);
	X(spcelldroptol)(*A, EPS);
	X(spcellsort)(*A);
}

/**
   Create a sparse convolution operator C with
   C(i,j)=A(i-j);
   A must be very sparse with only a view non-zero value otherwise C will be too full.
*/
X(sp)* X(spconvolvop)(X(mat)* A){
	/*First collect statistics on A. */
	long nini=10;
	T* vals=mycalloc(nini, T);
	spint* sepx=mycalloc(nini, spint);
	spint* sepy=mycalloc(nini, spint);
	long count=0;
	const long nx=A->nx;
	const long ny=A->ny;
	const long nn=nx*ny;
	X(mat)* PA=A;
	for(long iy=0; iy<A->ny; iy++){
		for(long ix=0; ix<A->nx; ix++){
			if(fabs(P(PA, ix, iy))>0){
				vals[count]=P(PA, ix, iy);
				sepx[count]=ix;
				sepy[count]=iy;
				count++;
			}
			if(count>=nini){
				nini*=2;
				vals=myrealloc(vals, nini, T);
				sepx=myrealloc(sepx, nini, spint);
				sepy=myrealloc(sepy, nini, spint);
			}
		}
	}
	if(count>10){
		warning("Number of coupled points %ld is too large\n", count);
	}
	long nsep=count;
	X(sp)* out=X(spnew)(nn, nn, nn*count);
	spint* pp=out->pp;
	spint* pi=out->pi;
	T* px=out->px;
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
   Permute rows of sparse matrix A;
   let pi be the inverse order of p. pi=invperm(p);
   if reverse is 0
   B(:,i)=A(p,i); or B(pi,i)=A(:,i)
   if reverse is 1
   B(p,i)=A(:,i); or B(:,i)=A(pi,i);
*/
static X(sp)* X(sppermcol)(const X(sp)* A, int reverse, long* p){
	if(!p) return X(spdup)(A);
	X(sp)* B=X(spnew)(A->nx, A->ny, A->nzmax);
	long* perm;
	if(reverse){
		perm=p;
	} else{
		perm=invperm(p, A->nx);
	}

	for(long icol=0; icol<A->ny; icol++){
		B->pp[icol]=A->pp[icol];
		for(long irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
			long row=A->pi[irow];
			B->pi[irow]=perm[row];
			B->px[irow]=A->px[irow];
		}
	}
	B->pp[A->ny]=A->pp[A->ny];
	if(!reverse){
		free(perm);
	}
	return B;
}
/**
   Permute rows and columns of sparse matrix A;
*/
X(sp)* X(spperm)(X(sp)* A, int reverse, long* pcol, long* prow){
	if(!A) return 0;
	assert_sp(A);
	X(sp)* out;
	if(pcol){
		out=X(sppermcol)(A, reverse, pcol);
	} else{
		out=X(spref)(A);
	}
	if(prow){
		X(sp)* Ap=X(sptrans)(out);
		X(sp)* App=X(sppermcol)(Ap, reverse, prow);
		X(spfree)(Ap);
		X(spfree)(out);
		out=X(sptrans)(App);
		X(spfree)(App);
	}
	return out;
}

/**
   Invert a SPD sparse matrix that is block diagonal with
   block sizes of bs.
*/
X(sp)* X(spinvbdiag)(const X(sp)* A, long bs){
	if(!A) return 0;
	assert_sp(A);
	if(A->nx!=A->ny){
		error("Must be a square matrix\n");
	}
	long nb=A->nx/bs;
	X(sp)* B=X(spnew)(A->nx, A->ny, nb*bs*bs);
	X(mat)* bk=X(new)(bs, bs);
	X(mat)* pbk=bk;
	for(long ib=0;ib<nb; ib++){
		long is=ib*bs;/*starting col */
		X(zero)(bk);

		for(long icol=is; icol<is+bs; icol++){
			for(long irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				long row=A->pi[irow];
				long ind=row-is;
				if(ind<0||ind>=bs){
					dbg("solving block %ld\n", ib);
					error("The array is not block diagonal matrix or not calculated properly.\n");
				}
				P(pbk, ind, icol-is)=A->px[irow];
			}
		}
		X(inv_inplace)(bk);
		for(long icol=is; icol<is+bs; icol++){
			B->pp[icol]=icol*bs;
			for(long irow=0; irow<bs; irow++){
				B->pi[B->pp[icol]+irow]=irow+is;
				B->px[B->pp[icol]+irow]=P(pbk, irow, icol-is);
			}
		}
	}
	B->pp[A->ny]=nb*bs*bs;
	X(free)(bk);
	return B;
}
/**
   Extrat the diagonal blocks of size bs into cell arrays.
*/
X(cell)* X(spblockextract)(const X(sp)* A, long bs){
	if(!A) return 0;
	assert_sp(A);
	if(A->nx!=A->ny){
		error("Must be a square matrix\n");
	}
	long nb=A->nx/bs;
	X(cell)* out=X(cellnew_same)(nb, 1, bs, bs);
	for(long ib=0;ib<nb; ib++){
		long is=ib*bs;/*starting col */
		//P(out,ib)=X(new)(bs, bs);
		X(mat)* pbk=P(out,ib);
		for(long icol=is; icol<is+bs; icol++){
			for(long irow=A->pp[icol]; irow<A->pp[icol+1]; irow++){
				long row=A->pi[irow];
				long ind=row-is;
				if(ind<0||ind>=bs){
					dbg("solving block %ld\n", ib);
					error("The array is not block diagonal matrix or not calculated property\n");
				}
				P(pbk, ind, icol-is)=A->px[irow];
			}
		}
	}
	return out;
}
