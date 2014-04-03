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
#include <unistd.h>
#include "../sys/sys.h"
#include "cholmod.h"
#include "mathdef.h"
#include "chol.h"
#if defined(DLONG)
#define MOD(A) cholmod_l_##A
#define ITYPE CHOLMOD_LONG
#else
#define MOD(A) cholmod_##A
#define ITYPE CHOLMOD_INT
#endif
#define CHOL_SIMPLE 1 //essential for chol_convert to work.

/**
   Convert our dsp spase type to cholmod_sparse type. Data is shared. 
*/
static cholmod_sparse *sp2chol(const dsp *A){
    cholmod_sparse *B=calloc(1, sizeof(cholmod_sparse));
    B->nrow=A->m;
    B->ncol=A->n;
    B->nzmax=A->nzmax;
    B->p=A->p;/*do not duplicate. */
    B->i=A->i;
    B->x=A->x;
    B->stype=1;/*assume data is symmetric. upper triangular part is used. */
    B->itype=ITYPE;
    B->xtype=CHOLMOD_REAL;
    B->dtype=CHOLMOD_DOUBLE;
    B->sorted=0;
    B->packed=1;
    B->nz=NULL;
    return B;
}
/**
   Convert cholmod_sparse to dsp. Data is shared
*/
static dsp *chol2sp(const cholmod_sparse *B){
    dsp *A;
    A=spnew(B->nrow, B->ncol, 0);
    A->p=B->p;
    A->i=B->i;
    A->x=B->x;
    A->nzmax=B->nzmax;
    return A;
}
/**
   Convert our dmat type to cholmod_dense type.
*/
static cholmod_dense* d2chol(const dmat *A, int start, int end){
    cholmod_dense* B=calloc(1, sizeof(cholmod_dense));
    if(end==0) end=A->ny;
    B->nrow=A->nx;
    B->ncol=end-start;
    B->nzmax=B->nrow*B->ncol;
    B->d=A->nx;
    B->x=A->p+start*A->nx;
    B->z=NULL;
    B->xtype=CHOLMOD_REAL;
    B->dtype=CHOLMOD_DOUBLE;
    return B;
}

/**
   Factorize a sparse array into LL' with reordering.
*/
spchol* chol_factorize(dsp *A_in){
    if(!A_in) return NULL;
    spchol *out=calloc(1, sizeof(spchol));
    out->c=calloc(1, sizeof(cholmod_common));
    MOD(start)(out->c);
    cholmod_sparse *A=sp2chol(A_in);
    out->c->status=CHOLMOD_OK;
#if CHOL_SIMPLE
	out->c->final_super=0;/*we want a simplicity result. */
	out->c->final_ll=1;   /*Leave in LL instead of LDL format. */
	out->c->final_asis=0; /*do the conversion as shown above. */
#endif    
    /*Try AMD ordering only: SLOW */
    /*
      out->c.nmethods=1;
      out->c.method[0].ordering=CHOLMOD_AMD;
      out->c.postorder=1;
      out->c.supernodal=CHOLMOD_SIMPLICIAL; force simplicial only. 
    */
    
    out->L=MOD(analyze)(A,out->c);
    if(!out->L) {
	info("\nCholmod error:");
	switch(out->c->status){
	case CHOLMOD_OK:
	    info2("Success\n");break;
	case CHOLMOD_NOT_INSTALLED:
	    info2("Method not installed\n"); break;
	case CHOLMOD_OUT_OF_MEMORY:
	    info2("Out of memory\n");break;
	case CHOLMOD_TOO_LARGE:
	    info2("Integer overflow occured\n"); break;
	case CHOLMOD_INVALID:
	    info2("Invalid input\n"); break;
	case CHOLMOD_NOT_POSDEF:
	    info2("Warning: Matrix not positive definite\n"); break;
	case CHOLMOD_DSMALL:
	    info2("Warning: D for LDL' or diag(L) for LL' has tiny absolute value\n"); break;
	default:
	    info2("Unknown error\n");
	}
	warning("Common->status=%d\n", out->c->status);
	error("Analyze failed\n");
    }
    MOD(factorize)(A,out->L, out->c);
#if CHOLMOD_SIMPLE
    if(!out->c->final_asis){
	/*Our solver is much slower than the simplicity solver, or the supernodal solver. */
	warning2("Converted to our format.");
	cholmod_factor *L=out->L;
	out->Cp=L->Perm; L->Perm=NULL;
	dsp *C=out->Cl=spnew(L->n, L->n, 0);
	C->p=L->p;L->p=NULL;
	C->i=L->i;L->i=NULL;
	C->x=L->x;L->x=NULL;
	C->nzmax=L->nzmax;
	MOD(free_factor)(&out->L, out->c);
	MOD(finish)(out->c);
	free(out->c);
	out->c=NULL;
	out->L=NULL;
    }
#endif
    free(A);/*just free our reference.*/
    return out;
}

/**
   Convert the internal data type cholesky factor into the lower left diagonal
   and permutation vector.  The internal data A->L is freed if keep is 0.

   This routine works only if the factor is using simplicity factor. not
   supernodal. so, in chol_factorize, we have set c.final_super=0, so that the
   final result is always in simplicity factor.
*/
void chol_convert(spchol *A, int keep){
    if(!A || !A->L || A->Cp || A->Cl) return;
#if ! CHOL_SIMPLE
    error("chol_convert only work with CHOL_SIMPLE=1\n");
#endif
    cholmod_factor *L=A->L;
    A->Cp=malloc(sizeof(spint)*A->L->n);
    memcpy(A->Cp, A->L->Perm, sizeof(spint)*A->L->n);
    if(keep){
	L=MOD(copy_factor)(A->L, A->c);
    }else{
	L=A->L;
    }
    cholmod_sparse *B=MOD(factor_to_sparse)(L, A->c);
    A->Cl=chol2sp(B);
    free(B);
    MOD(free_factor)(&L, A->c);
    if(!keep){
	MOD(finish)(A->c);
	free(A->c);
	A->c=NULL;
	A->L=NULL;
    }
}
/**
   Save cholesky factor and permutation vector to file.*/
void chol_save(spchol *A, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    dsp *C=A->Cl?A->Cl:A->Cu;
    long nc=0;
    if(C){/*Save our easy to use format. */
	header_t header={MCC_ANY, 2, 1, NULL};
	write_header(&header, fp);
	spwritedata(fp, C);
	do_write(fp, 0, sizeof(spint), M_SPINT, "Cp", A->Cp, C->m, 1);
    }else if(A->L){/*Save native cholmod format. */
	cholmod_factor *L=A->L;
	char str[1024];
	snprintf(str,1024,
		 "n=%zd\n"
		 "minor=%zd\n"
		 "nzmax=%zd\n"
		 "nsuper=%zd\n"
		 "ssize=%zd\n"
		 "xsize=%zd\n"
		 "maxcsize=%zd\n"
		 "maxesize=%zd\n"
		 "ordering=%d\n"
		 "is_ll=%d\n"
		 "is_super=%d\n"
		 "is_monotonic=%d\n"
		 "itype=%d\n"
		 "xtype=%d\n"
		 "dtype=%d\n"
		 ,L->n, L->minor,L->nzmax,L->nsuper,L->ssize,L->xsize,L->maxcsize,L->maxesize,
		 L->ordering,L->is_ll,L->is_super,L->is_monotonic,L->itype,L->xtype,L->dtype);
	nc=L->is_super?7:8;
	header_t header={MCC_ANY, nc, 1, str};
	write_header(&header, fp);
	do_write(fp, 0, sizeof(spint), M_SPINT, "Perm", L->Perm, L->n, 1);
	do_write(fp, 0, sizeof(spint), M_SPINT, "ColCount", L->ColCount, L->n, 1);
	if(L->is_super==0){/*Simplicity */
	    do_write(fp, 0, sizeof(spint), M_SPINT, "p", L->p, L->n+1, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "i", L->i, L->nzmax, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "nz", L->nz, L->n, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "next", L->next, L->n+2, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "prev", L->prev, L->n+2, 1);
	}else{
	    do_write(fp, 0, sizeof(spint), M_SPINT, "super", L->super, L->nsuper+1, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "pi", L->pi, L->nsuper+1, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "px", L->px, L->nsuper+1, 1);
	    do_write(fp, 0, sizeof(spint), M_SPINT, "s", L->s, L->ssize, 1);
	}
	do_write(fp, 0, sizeof(double),M_DBL, "x", L->x, L->is_super?L->xsize:L->nzmax, 1);
    }
    zfclose(fp);
}

/**
   Read cholesky factor form file. In the matlab convention, the upper triangle
   is saved. */
spchol *chol_read(const char *format, ...){
    format2fn;
    spchol *A=calloc(1, sizeof(spchol));
    file_t *fp=zfopen(fn, "rb");
    header_t header;
    read_header(&header, fp);
    if(!iscell(header.magic)){
	error("%s does not contain cell array\n", fn);
    }
    long ncx=header.nx;
    long ncy=header.ny;;
    if(ncx*ncy==2){/*Contains Cl(Cu) and Perm */
	dsp *C=spreaddata(fp, 0);
	int type=spcheck(C);
	if(type & 1){
	    A->Cl=C;
	}else if(type & 2){
	    A->Cu=C;
	}else{
	    error("C is not upper or lower\n");
	}
	long nxp, nyp;
	A->Cp=readspint(fp, &nxp, &nyp);
	if(nxp*nyp!=C->m){
	    error("Cp is in wrong format\n");
	}
    }else{/*Native cholmod format exists. Read it. */
	if(!header.str){
	    error("File %s does not contain a cholmod_factor\n", fn);
	}
#define READ_SIZE_T(A) L->A=(size_t)search_header_num(header.str,#A)
#define READ_INT(A) L->A=(int)search_header_num(header.str,#A)
	cholmod_factor *L=A->L=calloc(1, sizeof(cholmod_factor));
	READ_SIZE_T(n);
	READ_SIZE_T(minor);
	READ_SIZE_T(nzmax);
	READ_SIZE_T(nsuper);
	READ_SIZE_T(ssize);
	READ_SIZE_T(xsize);
	READ_SIZE_T(maxcsize);
	READ_SIZE_T(maxesize);
	READ_INT(ordering);
	READ_INT(is_ll);
	READ_INT(is_super);
	READ_INT(is_monotonic);
	READ_INT(itype);
	READ_INT(xtype);
	READ_INT(dtype);
#undef READ_SIZE_T
#undef READ_INT
	long nx, ny;
	header_t header2;
#define READSPINT(A,N) L->A=readspint(fp, &nx, &ny);			\
	if((N)!=nx*ny) error("%s has wrong length: wanted %ld, got %ld\n", #A, (long)(N), nx*ny);
#define READDBL(A,N) read_header(&header2,fp);				\
	if(header2.magic!=M_DBL) error("Invalid magic: wanted %u, got %u\n", M_DBL, header2.magic); \
	nx=header2.nx; ny=header2.ny;					\
	if(nx*ny!=(N)) error("%s has wrong length: wanted %ld, got %ld\n", #A, (long)(N), nx*ny); \
	L->A=malloc(sizeof(double)*nx*ny);				\
	zfread(L->A, sizeof(double), nx*ny, fp);
	
	READSPINT(Perm, L->n);
	READSPINT(ColCount, L->n);
	if(L->is_super==0){/*Simplicity */
	    info2("Reading simplicity cholmod_factor\n");
	    READSPINT(p, L->n+1);
	    READSPINT(i, L->nzmax);
	    READSPINT(nz, L->n);
	    READSPINT(next, L->n+2);
	    READSPINT(prev, L->n+2);
	}else{
	    info2("Reading supernodal cholmod_factor\n");
	    READSPINT(super, L->nsuper+1);
	    READSPINT(pi, L->nsuper+1);
	    READSPINT(px, L->nsuper+1);
	    READSPINT(s, L->ssize);
	}
	READDBL(x, L->is_super?L->xsize:L->nzmax);
	A->c=calloc(1, sizeof(cholmod_common));
	MOD(start)(A->c);
#undef READSPINT
#undef READDDBL
    }
    zfclose(fp);
    return A;
}

typedef struct{
    dmat *x;
    spchol *A;
    const dmat *y;
}CHOLSOLVE_T;
/*
  Solve section of the columns in multi-threaded way.
*/
static void chol_solve_each(thread_t *info){
    CHOLSOLVE_T *data=info->data;
    dmat *x=data->x;
    spchol *A=data->A;
    cholmod_dense *y2=d2chol(data->y, info->start, info->end);
    cholmod_dense *x2=MOD(solve)(CHOLMOD_A, A->L, y2, A->c);
    memcpy(x->p+info->start*x->nx, x2->x, x->nx*(info->end-info->start)*sizeof(double));
    free(y2);
    free(x2->x);
    free(x2);
}
/**
   Solve A*x=y where the cholesky factor of A is stored in A.
*/
void chol_solve(dmat **x, spchol *A, dmat *y){
    if(!A->L){
	if(A->Cl){
	    chol_solve_lower(x, A, y);
	}else if(A->Cu){
	    chol_solve_upper(x, A, y);
	}else{
	    error("There is no cholesky factor\n");
	}
    }else{
	assert(A->L->xtype!=0);/* error("A->L is pattern only!\n"); */
	if(y->ny==1){
	    cholmod_dense *y2=d2chol(y, 0, y->ny);/*share pointer. */
	    cholmod_dense *x2=NULL;
	    x2=MOD(solve)(CHOLMOD_A,A->L,y2,A->c);
	    if(!x2) error("chol_solve failed\n");
	    if(x2->z) error("why is this?\n");
	    if(x2->nzmax!=x2->nrow*x2->ncol || x2->d!=x2->nrow){
		error("Fix here\n");
	    }
	    if(!*x){
		*x=dnew_data(x2->nrow,x2->ncol, x2->x);/*takes over the owner of x2->x. */
	    }else{
		if((*x)->nx!=x2->nrow || (*x)->ny!=x2->ncol){
		    error("Matrix mismatch\n");
		}
		memcpy((*x)->p,x2->x,sizeof(double)*((*x)->nx)*((*x)->ny));
		free(x2->x);
	    }
	    free(y2);/*don't do dfree */
	    free(x2);/*don't do dfree */
	}else{
	    if(!*x) *x=dnew(y->nx, y->ny);
	    CHOLSOLVE_T data={*x,A,y};
	    thread_t info[NTHREAD];
	    thread_prep(info, 0, y->ny, NTHREAD, chol_solve_each, &data);
	    CALL_THREAD(info, NTHREAD, 1);
	}
    } 
}
/**
   Solve A*x=y where Y is sparse. Not good idea to use dsp ** as chol_solve
   because the result may have different nzmax.
*/
dsp *chol_spsolve(spchol *A, const dsp *y){
    assert(A->L->xtype!=0);/* error("A->L is pattern only!\n"); */
    cholmod_sparse *y2=sp2chol(y);
    cholmod_sparse *x2=MOD(spsolve)(CHOLMOD_A,A->L,y2,A->c);
    if(!x2) error("chol_solve failed\n");
    if(x2->z) error("why is this?\n");
    dsp *x=spnew(x2->nrow, x2->ncol, 0);
    x->p=x2->p;
    x->i=x2->i;
    x->x=x2->x;
    x->nzmax=x2->nzmax;
    free(y2);/*don't do spfree */
    free(x2);/*don't do spfree */
    return x;
}

/**
   forward permutation.
*/
static inline void chol_perm_f(dmat **out, spint *perm, const dmat *in){
    if(!*out){
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx == in->nx && (*out)->ny == in->ny);
    }
    PDMAT(in,pin);
    PDMAT(*out,pout);
    if(*out==in){/*Do each column in place. */
	double *tmp=malloc(in->nx*sizeof(double));
	for(int icy=0; icy<in->ny; icy++){
	    for(int icx=0; icx<in->nx; icx++){
		tmp[icx]=pin[icy][perm[icx]];
	    }
	    memcpy(pout[icy], tmp, sizeof(double)*in->nx);
	}
	free(tmp);
    }else{
	for(int icy=0; icy<in->ny; icy++){
	    for(int icx=0; icx<in->nx; icx++){
		pout[icy][icx]=pin[icy][perm[icx]];
	    }
	}
    }
}
/**
   backward permutation.
*/
static inline void chol_perm_b(dmat **out, spint *perm, const dmat *in){
    if(!*out){
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx == in->nx && (*out)->ny == in->ny);
    }
    PDMAT(in,pin);
    PDMAT(*out,pout);
    if(*out==in){/*Do each column in place. */
	double *tmp=malloc(in->nx*sizeof(double));
	for(int icy=0; icy<in->ny; icy++){
	    for(int icx=0; icx<in->nx; icx++){
		tmp[perm[icx]]=pin[icy][icx];
	    }
	    memcpy(pout[icy], tmp, sizeof(double)*in->nx);
	}
	free(tmp);
    }else{
	for(int icy=0; icy<in->ny; icy++){
	    for(int icx=0; icx<in->nx; icx++){
		pout[icy][perm[icx]]=pin[icy][icx];
	    }
	}
    }
}
typedef struct{
    spchol *C;
    dmat *y2;
}CHOL_LOWER_T;
/*
  The performance of this is poorer than chol_solve_each. done in place
*/
static void chol_solve_lower_each(thread_t *info){
    CHOL_LOWER_T *data=info->data;
    dsp *A=data->C->Cl;
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    dmat *y2=data->y2;
    info2("Lower solving %ld x %ld, %ld\n", y2->nx, info->start, info->end);
    /*Solve L\y */
    PDMAT(y2, py);
    for(long icol=0; icol<A->n; icol++){
	double AxI=1./Ax[Ap[icol]];
	for(long iy=info->start; iy<info->end; iy++){
	    py[iy][icol]*=AxI;
	    double val=-py[iy][icol];
	    for(long irow=Ap[icol]+1; irow<Ap[icol+1]; irow++){
		py[iy][Ai[irow]]+=val*Ax[irow];/*update in place. */
	    }
	}
    }

    /*Solve L'\y; */
    for(long icol=A->n-1; icol>-1; icol--){
	double AxI=1./Ax[Ap[icol]];
	for(long iy=info->start; iy<info->end; iy++){
	    double sum=0;
	    /*We do in reverse order to increase memory reuse. 1.5xFaster than forward order. */
	    for(long irow=Ap[icol+1]-1; irow>Ap[icol]; irow--){
		sum+=Ax[irow]*py[iy][Ai[irow]];
	    }
	    py[iy][icol]=(py[iy][icol]-sum)*AxI;
	}
    }
}
/*
  The performance of this is poorer than chol_solve_each. done in place
*/
static void chol_solve_upper_each(thread_t *info){
    CHOL_LOWER_T *data=info->data;
    dsp *A=data->C->Cu;
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    dmat *y2=data->y2;
    info2("Upper solving %ld x (%ld to %ld)\n", y2->nx, info->start, info->end);
    /*Solve L\y */
    PDMAT(y2, py);
    /*Solve R'\y */
    for(long icol=0; icol<A->m; icol++){
	double AxI=1./Ax[Ap[icol+1]-1];
	for(long iy=info->start; iy<info->end; iy++){
	    double sum=0;
	    for(long irow=Ap[icol]; irow<Ap[icol+1]-1; irow++){
		sum+=Ax[irow]*py[iy][Ai[irow]];
	    }
	    /*assert(Ai[Ap[icol+1]-1]==icol);//confirm upper right triangular */
	    py[iy][icol]=(py[iy][icol]-sum)*AxI;
	}
    }
	
    /*Solve R\y */
    for(long icol=A->m-1; icol>-1; icol--){
	double AxI=1./Ax[Ap[icol+1]-1];
	for(long iy=info->start; iy<info->end; iy++){
	    py[iy][icol]*=AxI;
	    double val=-py[iy][icol];
	    /*We do in reverse order to increase memory reuse. 1.5xFaster than forward order. */
	    for(long irow=Ap[icol+1]-2; irow>Ap[icol]-1; irow--){
		py[iy][Ai[irow]]+=val*Ax[irow];
	    }
	}
    }
}
/*
  2010-08-09: The following two routines are a bit slower than the cholmod
  solver because the later uses supernodal paritition.  */
/**
   Solve A*x=Y where the A=LL' and L is stored in A.
   
   Solve the cholesky backsubstitution when it's expressed in sparse matrix and
   a permutation vector. Notice only the lower left side of the sparse matrix is
   stored.
   
   The original matrix B=A*A';
   solve B*x=y or A*(A'*x)=y 
   first solve A\\y
   then solve A'\\(A\\y)

   This is slower than chol_solve in cholmod. I observed high cpu usage burst
   during chol_solve, which is probably due to call of blas which is
   multi-threaded. I don't call blas here.  */
void chol_solve_lower(dmat **x, spchol *C, dmat *y){
    if(!C->Cl){
	info2("Converting spchol\n");
	chol_convert(C, 1);
    }
    dsp *A=C->Cl;
    spint *perm=C->Cp;
    assert(A && A->m==A->n && A->m==y->nx);
    if(!*x){
	*x=dnew(y->nx,y->ny);
    }else{
	assert((*x)->nx==y->nx && (*x)->ny==y->ny);
    }
    dmat *y2=NULL;
    if(*x==y) y2=y;/*do inplace. */
    chol_perm_f(&y2, perm, y);
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    if(y2->ny==1){
	/*Solve L\y */
	double *py=y2->p;
	
	for(long icol=0; icol<A->n; icol++){
	    /*assert(Ai[Ap[icol]]==icol);//lower triangular matrix. */
	    py[icol]/=Ax[Ap[icol]];
	    double val=-py[icol];
	    for(long irow=Ap[icol]+1; irow<Ap[icol+1]; irow++){
		py[Ai[irow]]+=val*Ax[irow];/*update in place. */
	    }
	}

	/*Solve L'\y; */
	for(long icol=A->n-1; icol>-1; icol--){
	    double sum=0;
	    /*We do in reverse order to increase memory reuse. 1.5xFaster than forward order.*/
	    for(long irow=Ap[icol+1]-1; irow>Ap[icol]; irow--){
		sum+=Ax[irow]*py[Ai[irow]];
	    }
	    py[icol]=(py[icol]-sum)/Ax[Ap[icol]];
	}
    }else{/*When there are multiple columns. */
	CHOL_LOWER_T data={C,y2};
	thread_t info[NTHREAD];
	thread_prep(info, 0, y2->ny, NTHREAD, chol_solve_lower_each, &data);
	CALL_THREAD(info, NTHREAD, 1);
    }
    chol_perm_b(x, perm, y2);
    if(*x!=y2) dfree(y2);
}
/**
   Solve A*x=Y where the A=U'U and U is stored in A.

   Solve the cholesky backsubstitution when it's expressed in sparse matrix and
   a permutation vector. Notice only the lower left side of the sparse matrix is
   stored.
   
   The original matrix B=A'*A;
*/
void chol_solve_upper(dmat **x, spchol *C, dmat *y){
    if(!C->Cu){
	if(!C->Cl){
	    info2("Converting spchol\n");
	    chol_convert(C, 1);
	}
	info2("Transposing C->Cl");
	C->Cu=sptrans(C->Cl);
    }
    dsp *A=C->Cu;
    spint *perm=C->Cp;
    assert(A && A->m==A->n && A->m==y->nx);
    if(!*x){
	*x=dnew(y->nx,y->ny);
    }else{
	assert((*x)->nx==y->nx && (*x)->ny==y->ny);
    }
    dmat *y2=NULL;
    if(*x==y) y2=y;/*do inplace. */
    chol_perm_f(&y2, perm, y);
    double *Ax=A->x;
    spint *Ap=A->p;
    spint *Ai=A->i;
    if(y2->ny==1){
	/*Solve R'\y */
	double *py=y2->p;

	for(long icol=0; icol<A->m; icol++){
	    double sum=0;
	    for(long irow=Ap[icol]; irow<Ap[icol+1]-1; irow++){
		sum+=Ax[irow]*py[Ai[irow]];
	    }
	    py[icol]=(py[icol]-sum)/Ax[Ap[icol+1]-1];
	}

	/*Solve R\y */
	for(long icol=A->m-1; icol>-1; icol--){
	    py[icol]/=Ax[Ap[icol+1]-1];
	    double val=-py[icol];
	    /*We do in reverse order to increase memory reuse. 1.5xFaster than
	      forward order.*/
	    for(long irow=Ap[icol+1]-2; irow>Ap[icol]-1; irow--){
		py[Ai[irow]]+=val*Ax[irow];
	    }
	}
    }else{
	CHOL_LOWER_T data={C,y2};
	thread_t info[y2->ny];
	thread_prep(info, 0, y2->ny, y2->ny, chol_solve_upper_each, &data);
	CALL_THREAD(info, y2->ny, 1);
    }
    chol_perm_b(x, perm, y2);
    if(*x!=y2) dfree(y2);
}

/**
   Free cholesky factor.*/
void chol_free_do(spchol *A){
    if(A){
	if(A->L){
	    MOD(free_factor)(&A->L, A->c);
	    MOD(finish)(A->c);
	    free(A->c);
	}
	if(A->Cp) free(A->Cp);
	if(A->Cl) spfree(A->Cl);
	if(A->Cu) spfree(A->Cu);
	free(A);
    }
}
