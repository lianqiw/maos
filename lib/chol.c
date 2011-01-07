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
#include "cholmod.h"
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"
#include "matbin.h"
#include "chol.h"
#if defined(DLONG)
#define MOD(A) cholmod_l_##A
#define ITYPE CHOLMOD_LONG
#else
#define MOD(A) cholmod_##A
#define ITYPE CHOLMOD_INT
#endif
struct spchol{
    cholmod_factor *L;
    cholmod_common *c;
};
/**
\file chol.c
Wraps the CHOLESKY Library to provide a simple interface.*/

/**
   convert out dsp spase type to cholmod_sparse type.
 */
static cholmod_sparse *sp2chol(dsp *A){
    cholmod_sparse *B=calloc(1, sizeof(cholmod_sparse));
    B->nrow=A->m;
    B->ncol=A->n;
    B->nzmax=A->nzmax;
    B->p=A->p;//do not duplicate.
    B->i=A->i;
    B->x=A->x;
    B->stype=1;//assume data is symmetric. upper triangular part is used.
    B->itype=ITYPE;
    B->xtype=CHOLMOD_REAL;
    B->dtype=CHOLMOD_DOUBLE;
    B->sorted=0;
    B->packed=1;
    B->nz=NULL;
    return B;
}
/**
   Convert our dmat type to cholmod_dense type.
 */
static cholmod_dense* d2chol(const dmat *A){
    cholmod_dense* B=calloc(1, sizeof(cholmod_dense));
    B->nrow=A->nx;
    B->ncol=A->ny;
    B->nzmax=A->nx*A->ny;
    B->d=A->nx;
    B->x=A->p;
    B->z=NULL;
    B->xtype=CHOLMOD_REAL;
    B->dtype=CHOLMOD_DOUBLE;
    return B;
}
/**
   Convert spchol to sparse matrix.
   keep=1: A is kept intact, otherwise destroyed. The matrix is in lower left side.

   This routine works only if the factor is using simplicity factor. not
   supernodal. so, in chol_factorize, we have set c.final_super=0, so that the
   final result is always in simplicity factor.  */
static dsp* chol_sp(spchol *A, int keep){
    if(!A) return NULL;

    cholmod_factor *L;
    if(keep){
	L=MOD(copy_factor)(A->L, A->c);
    }else{
	L=A->L;
    }
    cholmod_sparse *B=MOD(factor_to_sparse)(L, A->c);
    dsp *out=spnew(B->nrow, B->ncol, 0);
    out->p=B->p;
    out->i=B->i;
    out->x=B->x;
    out->nzmax=B->nzmax;
    free(B);
    MOD(free_factor)(&L, A->c);
    if(!keep){
	MOD(finish)(A->c);
	free(A);
    }
    return out;
}
/**
   Factorize a sparse array into LL'.
*/
spchol* chol_factorize(dsp *A_in){
    if(!A_in) return NULL;
    TIC;tic;
    spchol *out=calloc(1, sizeof(spchol));
    out->c=calloc(1, sizeof(cholmod_common));
    MOD(start)(out->c);
    cholmod_sparse *A=sp2chol(A_in);
    out->c->status=CHOLMOD_OK;
    out->c->final_super=0;//we want a simple result
    {
	//Try AMD ordering only. SLOW
	/*
	  out->c.nmethods=1;
	  out->c.method[0].ordering=CHOLMOD_AMD;
	  out->c.postorder=1;
	  out->c.supernodal=CHOLMOD_SIMPLICIAL;//force simplicial only.
	*/
    }
    info2("analyzing...");
    out->L=MOD(analyze)(A,out->c);
    info2("factoring...");
    if(!out->L) error("Analyze failed\n");
    MOD(factorize)(A,out->L, out->c);
    free(A);
    toc2("done.");
    return out;
}
/**
   Solve A*x=Y where the cholesky factor of A is stored in A.
*/
void chol_solve(dmat **x, spchol *A, const dmat *y){
    //solve A*x=Y;
    cholmod_dense *y2=d2chol(y);//share pointer.
    if(A->L->xtype==0) error("A->L is pattern only!\n");
    cholmod_dense *x2=MOD(solve)(CHOLMOD_A,A->L,y2,A->c);
    if(!x2) error("chol_solve failed\n");
    if(x2->z){
	error("why is this?\n");
    }
    if(x2->nzmax!=x2->nrow*x2->ncol || x2->d!=x2->nrow){
	error("Fix here\n");
    }
    free(y2);
    if(!*x){
	*x=dnew_data(x2->x,x2->nrow,x2->ncol);//takes over the owner of x2->x.
    }else{
	if((*x)->nx!=x2->nrow || (*x)->ny!=x2->ncol){
	    error("Matrix mismatch\n");
	}
	memcpy((*x)->p,x2->x,sizeof(double)*((*x)->nx)*((*x)->ny));
	free(x2->x);
    }
    free(x2);//x2->x is kept.
}
/**
   Free cholesky factor.*/
void chol_free_do(spchol *A){
    if(A){
	MOD(free_factor)(&A->L, A->c);
	MOD(finish)(A->c);
	free(A->c);
	free(A);
    }
}
/**
   Save cholesky factor: the lower left side and permutation vector.*/
void chol_save(spchol *A, const char *format,...){
    format2fn;
    char *fn2=malloc(strlen(fn)+10);
    dsp *Chol=chol_sp(A, 1);
    spwrite(Chol,"%s_C.bin",fn);
    spfree(Chol);
    writeint64(A->L->Perm,A->L->n,1,"%s_P.bin.gz",fn);//fixme: what happens in 32bit machine?
    free(fn2);
}
/**
   Convert the internal data type cholesky factor into the lower left diagonal
   and permutation vector
 */
void chol_convert(dsp **Cs, long **Cp, spchol *A, int keep){
    /*
       Convert the internal format spchol to a simple sparse and a reordering vector.
       Keep=1: A is kept.
       Keep=0: A is destroyed.
    */
    if(!A){
	*Cs=NULL;
	*Cp=NULL;
    }
    *Cp=malloc(sizeof(long)*A->L->n);//fixme: is this right in 32 bit machine?
    memcpy(*Cp, A->L->Perm, sizeof(long)*A->L->n);
    *Cs=chol_sp(A, keep);
}
/**
   forward permutation.
*/
static inline void chol_perm_f(dmat **out, long *perm, const dmat *in){
    if(!*out){
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx == in->nx && (*out)->ny == in->ny);
    }
    PDMAT(in,pin);
    PDMAT(*out,pout);
    for(int icy=0; icy<in->ny; icy++){
	for(int icx=0; icx<in->nx; icx++){
	    pout[icy][icx]=pin[icy][perm[icx]];
	}
    }
}
/**
   backward permutation.
*/
static inline void chol_perm_b(dmat **out, long *perm, const dmat *in){
    if(!*out){
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx == in->nx && (*out)->ny == in->ny);
    }
    PDMAT(in,pin);
    PDMAT(*out,pout);
    for(int icy=0; icy<in->ny; icy++){
	for(int icx=0; icx<in->nx; icx++){
	    pout[icy][perm[icx]]=pin[icy][icx];
	}
    }
}
/*
  2010-08-09: The following two routines are a bit slower than the cholmod solver.
*/
/**
   Solve A*x=Y where the A=LL' and L is stored in A.
   
   Solve the cholesky backsubstitution when it's expressed in sparse matrix and
   a permutation vector. Notice only the lower left side of the sparse matrix is
   stored.
   
   The original matrix B=A*A';
   solve B*x=y or A*(A'*x)=y 
   first solve A\\y
   then solve A'\\(A\\y)
   
*/
void chol_solve_lower(dmat **x, dsp *A, long *perm, const dmat *y){
    assert(A && A->m==A->n && A->m==y->nx);
    if(!*x){
	*x=dnew(y->nx,y->ny);
    }else{
	assert((*x)->nx==y->nx && (*x)->ny==y->ny);
    }
    dmat *y2=NULL;
    chol_perm_f(&y2, perm, y);
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    if(y2->ny==1){
	//Solve L\y
	double *py=y2->p;
	
	for(long icol=0; icol<A->n; icol++){
	    //assert(Ai[Ap[icol]]==icol);//lower triangular matrix.
	    py[icol]/=Ax[Ap[icol]];
	    double val=-py[icol];
	    for(long irow=Ap[icol]+1; irow<Ap[icol+1]; irow++){
		py[Ai[irow]]+=val*Ax[irow];//update in place.
	    }
	}

	//Solve L'\y;
	for(long icol=A->n-1; icol>-1; icol--){
	    double sum=0;
	    //We do in reverse order to increase memory reuse. 1.5xFaster than forward order.
	    for(long irow=Ap[icol+1]-1; irow>Ap[icol]; irow--){
		sum+=Ax[irow]*py[Ai[irow]];
	    }
	    py[icol]=(py[icol]-sum)/Ax[Ap[icol]];
	}
    }else{
	//Solve L\y
	PDMAT(y2, py);
	for(long icol=0; icol<A->n; icol++){
	    double AxI=1./Ax[Ap[icol]];
	    for(long iy=0; iy<y2->ny; iy++){
		py[iy][icol]*=AxI;
		double val=-py[iy][icol];
		for(long irow=Ap[icol]+1; irow<Ap[icol+1]; irow++){
		    py[iy][Ai[irow]]+=val*Ax[irow];//update in place.
		}
	    }
	}

	//Solve L'\y;
	for(long icol=A->n-1; icol>-1; icol--){
	    double AxI=1./Ax[Ap[icol]];
	    for(long iy=0; iy<y2->ny; iy++){
		double sum=0;
		//We do in reverse order to increase memory reuse. 1.5xFaster than forward order.
		for(long irow=Ap[icol+1]-1; irow>Ap[icol]; irow--){
		    sum+=Ax[irow]*py[iy][Ai[irow]];
		}
		py[iy][icol]=(py[iy][icol]-sum)*AxI;
	    }
	}
    }
    chol_perm_b(x, perm, y2);
    dfree(y2);
}
/**
   Solve A*x=Y where the A=U'U and U is stored in A.

   Solve the cholesky backsubstitution when it's expressed in sparse matrix and
   a permutation vector. Notice only the lower left side of the sparse matrix is
   stored.
   
   The original matrix B=A'*A;
*/
void chol_solve_upper(dmat **x, dsp *A, long *perm, const dmat *y){
    /*
  
    */
    assert(A && A->m==A->n && A->m==y->nx);
    if(!*x){
	*x=dnew(y->nx,y->ny);
    }else{
	assert((*x)->nx==y->nx && (*x)->ny==y->ny);
    }
    dmat *y2=NULL;
    chol_perm_f(&y2, perm, y);
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    if(y2->ny==1){
	//Solve R'\y
	double *py=y2->p;

	for(long icol=0; icol<A->m; icol++){
	    double sum=0;
	    for(long irow=Ap[icol]; irow<Ap[icol+1]-1; irow++){
		sum+=Ax[irow]*py[Ai[irow]];
	    }
	    //assert(Ai[Ap[icol+1]-1]==icol);//confirm upper right triangular
	    py[icol]=(py[icol]-sum)/Ax[Ap[icol+1]-1];
	}

	//Solve R\y
	for(long icol=A->m-1; icol>-1; icol--){
	    py[icol]/=Ax[Ap[icol+1]-1];
	    double val=-py[icol];
	    //We do in reverse order to increase memory reuse. 1.5xFaster than forward order.
	    for(long irow=Ap[icol+1]-2; irow>Ap[icol]-1; irow--){
		py[Ai[irow]]+=val*Ax[irow];
	    }
	}
    }else{
	PDMAT(y2,py);
	//Solve R'\y
	for(long icol=0; icol<A->m; icol++){
	    for(long iy=0; iy<y2->ny; iy++){
		double sum=0;
		for(long irow=Ap[icol]; irow<Ap[icol+1]-1; irow++){
		    sum+=Ax[irow]*py[iy][Ai[irow]];
		}
		//assert(Ai[Ap[icol+1]-1]==icol);//confirm upper right triangular
		py[iy][icol]=(py[iy][icol]-sum)/Ax[Ap[icol+1]-1];
	    }
	}
	
	//Solve R\y
	for(long icol=A->m-1; icol>-1; icol--){
	    for(long iy=0; iy<y2->ny; iy++){
		py[iy][icol]/=Ax[Ap[icol+1]-1];
		double val=-py[iy][icol];
		//We do in reverse order to increase memory reuse. 1.5xFaster than forward order.
		for(long irow=Ap[icol+1]-2; irow>Ap[icol]-1; irow--){
		    py[iy][Ai[irow]]+=val*Ax[irow];
		}
	    }
	}
    }
    chol_perm_b(x, perm, y2);
    dfree(y2);
}

