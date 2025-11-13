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
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#if HAVE_SUITESPARSE_CHOLMOD_H
#include <suitesparse/cholmod.h>
#else
#include <cholmod.h>
#endif

#include "../math/mathdef.h"
#ifdef I
#undef I
#endif


#include "chol.h"
#if defined(DLONG)
#define MOD(A) cholmod_l_##A
#define CHOL_ITYPE CHOLMOD_LONG //CHOL_ITYPE is the same as spint.
#else
#define MOD(A) cholmod_##A
#define CHOL_ITYPE CHOLMOD_INT
#endif

/**
   2019-10-25: As of this writing, the cholmod library does not support single precision.
*/
#define AS_DOUBLE 1 //must always be true
#if !defined(COMP_SINGLE) || AS_DOUBLE
#define CHOL_DTYPE CHOLMOD_DOUBLE
typedef double chol_real;
#else
#define CHOL_DTYPE CHOLMOD_SINGLE
typedef float chol_real;
#endif

#define CHOL_SIMPLE 1 //essential for chol_convert to work.

#define DO_CONVERT(p1, p2, t1, t2, size)	\
    if(sizeof(t1)!=sizeof(t2)){			\
		p1=mymalloc(size, t1);			\
		for(long i=0; i<(long)size; i++){	\
		    ((t1*)p1)[i]=((t2*)p2)[i];		\
		}					\
    }else{					\
		p1=(t1*)p2;					\
    }
#define DO_COPY(p1, p2, t1, t2, size)		\
    if(sizeof(t1)!=sizeof(t2)){			\
		for(long i=0; i<(long)size; i++){	\
		    ((t1*)p1)[i]=((t2*)p2)[i];		\
		}					\
    }else{					\
		memcpy(p1, p2, sizeof(t1)*size);	\
    }
/**
   Convert our dsp spase type to cholmod_sparse type. Data is shared.
*/
static cholmod_sparse* dsp2chol(const dsp* A){
	cholmod_sparse* B=mycalloc(1, cholmod_sparse);
	B->nrow=A->nx;
	B->ncol=A->ny;
	B->nzmax=A->pp[A->ny];
	B->p=A->pp;/*do not duplicate. */
	B->i=A->pi;
	DO_CONVERT(B->x, A->px, chol_real, real, B->nzmax);
	B->dtype=CHOL_DTYPE;
	B->itype=CHOL_ITYPE;
	B->xtype=CHOLMOD_REAL;
	B->stype=1;/*assume data is symmetric. upper triangular part is used. */
	B->sorted=1;
	B->packed=1;
	B->nz=NULL;
	return B;
}

/**
   Convert cholmod_sparse to dsp. Data is shared
*/
static dsp* chol2sp(cholmod_sparse**pB){
	if(!pB || !*pB) return NULL;
	cholmod_sparse *B=*pB;
	dsp* A=dspnew(B->nrow, B->ncol, 0);
	A->pp=(spint*)B->p;
	A->nzmax=A->pp[A->ny];
	A->pi=(spint*)B->i;
	DO_CONVERT(A->px, B->x, real, chol_real, A->nzmax);
	if(A->px!=B->x){
		free(B->x);
	}
	free(B);//content is moved to A
	*pB=NULL;
	return A;
}
/**
   Convert our dmat type to cholmod_dense type.
*/
static cholmod_dense* d2chol(const dmat* A, int start, int end){
	cholmod_dense* B=mycalloc(1, cholmod_dense);
	if(end==0) end=A->ny;
	B->nrow=A->nx;
	B->ncol=end-start;
	B->nzmax=B->nrow*B->ncol;
	B->d=A->nx;
	B->z=NULL;
	B->xtype=CHOLMOD_REAL;
	DO_CONVERT(B->x, PCOL(A, start), chol_real, real, B->nzmax);
	B->dtype=CHOL_DTYPE;
	return B;
}

/**
   Factorize a sparse array into LL' with reordering.
*/
spchol* chol_factorize(const dsp* A_in){
	if(!A_in) return NULL;
	spchol* out=mycalloc(1, spchol);
	out->c=mycalloc(1, cholmod_common);
	MOD(start)(out->c);//This sets dtype to CHOLMOD_DOUBLE
	out->c->status=CHOLMOD_OK;
	out->c->itype=CHOL_ITYPE;
#if CHOL_SIMPLE
	out->c->final_super=0;/*we want a simplicity result. */
	out->c->final_ll=1;   /*Leave in LL instead of LDL format. */
	out->c->final_asis=0; /*0: do the conversion. */
	out->c->final_pack=1; /*In CSC format*/
#endif    
	/*Try AMD ordering only: SLOW */
	/*
	  out->c.nmethods=1;
	  out->c.method[0].ordering=CHOLMOD_AMD;
	  out->c.postorder=1;
	  out->c.supernodal=CHOLMOD_SIMPLICIAL; force simplicial only.
	*/
	int success=0;
	cholmod_sparse* A=dsp2chol(A_in);//borrows content.
	out->L=MOD(analyze)(A, out->c);
	if(out->L){
		success=MOD(factorize)(A, out->L, out->c);
	}
	if(!success){
		dbg("\nCholmod error:");
		switch(out->c->status){
		case CHOLMOD_OK:
			info("Success\n");break;
		case CHOLMOD_NOT_INSTALLED:
			info("Method not installed\n"); break;
		case CHOLMOD_OUT_OF_MEMORY:
			info("Out of memory\n");break;
		case CHOLMOD_TOO_LARGE:
			info("Integer overflow occured\n"); break;
		case CHOLMOD_INVALID:
			info("Invalid input\n"); break;
		case CHOLMOD_NOT_POSDEF:
			info("Warning: Matrix not positive definite\n"); break;
		case CHOLMOD_DSMALL:
			info("Warning: D for LDL' or diag(L) for LL' has tiny absolute value\n"); break;
		default:
			info("Unknown error\n");
		}
		warning("Common->status=%d\n", out->c->status);
		error("Cholmod factorize failed\n");
	}
#if CHOL_SIMPLE 
	if(!out->c->final_asis && sizeof(real)!=sizeof(chol_real) && 0){
		/*Our solver is much slower than the cholmod simplicity solver, or the
	 	* supernodal solver. Keep original data to use cholmod solver*/
	 	info("Converting chol result to sparse\n");
		chol_convert(out, 0);
		out->Cu=dsptrans(out->Cl); dspfree(out->Cl);
	}
#endif
	if(A->x!=A_in->px){
		free(A->x);
	}
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
void chol_convert(spchol* A, int keep){
	if(!A||!A->L||A->Cp||A->Cl) return;
#if ! CHOL_SIMPLE
	error("chol_convert only work with CHOL_SIMPLE=1\n");
#endif
	cholmod_factor* L=A->L;
	if(keep){
		A->Cp=(spint*)A->L->Perm;
		L=MOD(copy_factor)(A->L, A->c);
	} else{
		A->Cp=mymalloc(A->L->n, spint);
		memcpy(A->Cp, A->L->Perm, sizeof(spint)*A->L->n);
		L=A->L;
	}
	cholmod_sparse* B=MOD(factor_to_sparse)(L, A->c);
	A->Cl=chol2sp(&B); //moves content of B.
	MOD(free_factor)(&L, A->c);
	if(!keep){
		MOD(finish)(A->c);
		free(A->c);
		A->c=NULL;
		A->L=NULL;
	}
}
/**
   Alternative interface
*/
dsp* chol_factorize2(lmat** Cp, const dsp* A_in){
	spchol* C=chol_factorize(A_in);
	if(!C->Cl){
		chol_convert(C, 1);
	}
	if(Cp){
		*Cp=lnew(NX(A_in),1);
		DO_COPY(P(*Cp), C->Cp, long, spint, NX(A_in));
	}
	dsp* out=C->Cl; C->Cl=0;
	chol_free_do(C);
	return out;
}
/**
   Save cholesky factor and permutation vector to file.*/
void chol_save(spchol* A, const char* format, ...){
	format2fn;
	file_t* fp=zfopen(fn, "wb");
	dsp* C=A->Cl?A->Cl:A->Cu;
	long nc=0;
	if(C){/*Save our easy to use format. */
		header_t header={MCC_ANY, 2, 1, NULL};
		write_header(&header, fp);
		dspwritedata(fp, C);
		writearr(fp, 0, sizeof(spint), M_SPINT, "Cp", A->Cp, C->nx, 1);
	} else if(A->L){/*Save native cholmod format. */
		cholmod_factor* L=A->L;
		char str[1024];
		snprintf(str, 1024,
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
			, L->n, L->minor, L->nzmax, L->nsuper, L->ssize, L->xsize, L->maxcsize, L->maxesize,
			L->ordering, L->is_ll, L->is_super, L->is_monotonic, L->itype, L->xtype, L->dtype);
		nc=L->is_super?7:8;
		header_t header={MCC_ANY, (uint64_t)nc, 1, str};
		write_header(&header, fp);
		writearr(fp, 0, sizeof(spint), M_SPINT, "Perm", L->Perm, L->n, 1);
		writearr(fp, 0, sizeof(spint), M_SPINT, "ColCount", L->ColCount, L->n, 1);
		if(L->is_super==0){/*Simplicity */
			writearr(fp, 0, sizeof(spint), M_SPINT, "p", P(L), L->n+1, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "i", L->i, L->nzmax, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "nz", L->nz, L->n, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "next", L->next, L->n+2, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "prev", L->prev, L->n+2, 1);
		} else{
			writearr(fp, 0, sizeof(spint), M_SPINT, "super", L->super, L->nsuper+1, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "pi", L->pi, L->nsuper+1, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "px", L->px, L->nsuper+1, 1);
			writearr(fp, 0, sizeof(spint), M_SPINT, "s", L->s, L->ssize, 1);
		}
		if(L->dtype==CHOLMOD_DOUBLE){
			writearr(fp, 0, sizeof(double), M_DBL, "x", L->x, L->is_super?L->xsize:L->nzmax, 1);
		} else{
			writearr(fp, 0, sizeof(float), M_FLT, "x", L->x, L->is_super?L->xsize:L->nzmax, 1);
		}
	}
	zfclose(fp);
}

/**
   Read cholesky factor form file. In the matlab convention, the upper triangle
   is saved. */
spchol* chol_read(const char* format, ...){
	format2fn;
	spchol* A=mycalloc(1, spchol);
	file_t* fp=zfopen(fn, "rb");
	header_t header={0,0,0,0};
	read_header(&header, fp);
	if(!iscell(&header.id)){
		error("%s does not contain cell array\n", fn);
	}
	long ncx=header.nx;
	long ncy=header.ny;;
	if(ncx*ncy==2){/*Contains Cl(Cu) and Perm */
		dsp* C=dspreaddata(fp, 0);
		int type=dspcheck(C);
		if(type&1){
			A->Cl=C;
		} else if(type&2){
			A->Cu=C;
		} else{
			error("C is not upper or lower\n");
		}
		long nxp, nyp;
		A->Cp=readspint(fp, &nxp, &nyp);
		if(nxp*nyp!=C->nx){
			error("Cp is in wrong format\n");
		}
	} else{/*Native cholmod format exists. Read it. */
		if(!header.str){
			error("File %s does not contain a cholmod_factor\n", fn);
		}
#define READ_SIZE_T(A) L->A=(size_t)search_keyword_num(header.str,#A)
#define READ_INT(A) L->A=(int)search_keyword_num(header.str,#A)
		cholmod_factor* L=A->L=mycalloc(1, cholmod_factor);
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
		header_t header2={0,0,0,0};
#define READSPINT(A,N) L->A=readspint(fp, &nx, &ny);			\
	if((long)(N)!=nx*ny) error("%s has wrong length: wanted %ld, got %ld\n", #A, (long)(N), nx*ny);
#define READDBL(A,N) read_header(&header2,fp);				\
	nx=header2.nx; ny=header2.ny;					\
	if(nx*ny!=(long)(N)) error("%s has wrong length: wanted %ld, got %ld\n", #A, (long)(N), nx*ny); \
	L->A=mymalloc(nx*ny,real);					\
	readvec(L->A, M_REAL, header2.id, sizeof(real), nx*ny, fp);	
		READSPINT(Perm, L->n);
		READSPINT(ColCount, L->n);
		if(L->is_super==0){/*Simplicity */
			info("Reading simplicity cholmod_factor\n");
			READSPINT(p, L->n+1);
			READSPINT(i, L->nzmax);
			READSPINT(nz, L->n);
			READSPINT(next, L->n+2);
			READSPINT(prev, L->n+2);
		} else{
			info("Reading supernodal cholmod_factor\n");
			READSPINT(super, L->nsuper+1);
			READSPINT(pi, L->nsuper+1);
			READSPINT(px, L->nsuper+1);
			READSPINT(s, L->ssize);
		}
		READDBL(x, L->is_super?L->xsize:L->nzmax);
		A->c=mycalloc(1, cholmod_common);
		MOD(start)(A->c);
#undef READSPINT
#undef READDDBL
	}
	zfclose(fp);
	return A;
}

typedef struct{
	dmat* x;
	spchol* A;
	const dmat* y;
}CHOLSOLVE_T;
static void chol_solve_cholmod(dmat** x, const spchol* A, const dmat* y, long start, long end){
	cholmod_dense* y2=d2chol(y, start, end);
	cholmod_dense* x2=MOD(solve)(CHOLMOD_A, A->L, y2, A->c);
	if(y2->x!=PCOL(y, start)){
		free(y2->x);//data is converted. else: referenced
	}
	if(!x2||x2->z||x2->nzmax!=x2->nrow*x2->ncol||x2->d!=x2->nrow||start+x2->ncol!=(size_t)end){
		error("chol_solve failed or returns unexpected answer.\n");
	}

	if(!*x&&start==0&&sizeof(real)==sizeof(chol_real)){
		*x=dnew_do(x2->nrow, x2->ncol, (real*)x2->x, mem_new(x2->x));
	} else{
		if(!(*x)){
			*x=dnew(x2->nrow, end);
		} else if((*x)->nx!=(long)x2->nrow||(*x)->ny<(long)end){
			error("Matrix mismatch\n");
		}
		DO_COPY(PCOL(*x, start), x2->x, real, chol_real, (x2->nrow*x2->ncol));
		free_default(x2->x);
	}
	free(y2);/*keep data */
	free_default(x2);/*keep data */
}

/**
   Solve A*x=y where Y is sparse. Not good idea to use dsp ** as chol_solve
   because the result may have different nzmax.
*/
dsp* chol_spsolve(spchol* A, const dsp* y){
	assert(A->L->xtype!=0);/* error("A->L is pattern only!\n"); */
	cholmod_sparse* y2=dsp2chol(y);
	cholmod_sparse* x2=MOD(spsolve)(CHOLMOD_A, A->L, y2, A->c);
	if(!x2||x2->z) error("chol_solve failed or returns unexpected answer.\n");
	if(y2->x!=y->px) free(y2->x);
	free(y2);/*don't do spfree */
	dsp *x=chol2sp(&x2); //moves content
	return x;
}
/**
   forward permutation. In place permutation does not work.
*/
static inline void chol_perm(dmat** out, spint* perm, const dmat* in,
	long start_out, long start_in, long ncol, int direction){
	if(!*out){
		*out=dnew(in->nx, ncol+start_out);
	} else if(*out==in){
		error("Cannot do in place permutation");
	} else{
		assert((*out)->nx==in->nx&&(*out)->ny>=ncol+start_out);
	}
	for(long icy=0; icy<ncol; icy++){
		real* pout=PCOL(*out, icy+start_out);
		real* pin=PCOL(in, icy+start_in);
		if(direction==1){//forward
			for(long icx=0; icx<in->nx; icx++){
				pout[icx]=pin[perm[icx]];
			}
		} else{
			for(long icx=0; icx<in->nx; icx++){
				pout[perm[icx]]=pin[icx];
			}
		}
	}
}

/*
  The performance of this is poorer than chol_solve_each. done in place
*/
static void chol_solve_lower(const dsp* A, dmat* y2, long start, long end){
	real* Ax=A->px;
	spint* Ap=A->pp;
	spint* Ai=A->pi;
	/*Solve L\y */
	for(long icol=0; icol<A->ny; icol++){
		real AxI=1./Ax[Ap[icol]];
		for(long iy=start; iy<end; iy++){
			P(y2, icol, iy)*=AxI;
			real val=-P(y2, icol, iy);
			for(long irow=Ap[icol]+1; irow<Ap[icol+1]; irow++){
				P(y2, Ai[irow], iy)+=val*Ax[irow];/*update in place. */
			}
		}
	}

	/*Solve L'\y; */
	for(long icol=A->ny-1; icol>-1; icol--){
		real AxI=1./Ax[Ap[icol]];
		for(long iy=start; iy<end; iy++){
			real sum=0;
			/*We do in reverse order to increase memory reuse. 1.5xFaster than forward order. */
			for(long irow=Ap[icol+1]-1; irow>Ap[icol]; irow--){
				sum+=Ax[irow]*P(y2, Ai[irow], iy);
			}
			P(y2, icol, iy)=(P(y2, icol, iy)-sum)*AxI;
		}
	}
}
/*
  The performance of this is poorer than chol_solve_each. done in place
*/
static void chol_solve_upper(const dsp* A, dmat* y2, long start, long end){
	real* Ax=A->px;
	spint* Ap=A->pp;
	spint* Ai=A->pi;
	/*Solve R'\y */
	for(long icol=0; icol<A->nx; icol++){
		real AxI=1./Ax[Ap[icol+1]-1];
		for(long iy=start; iy<end; iy++){
			real sum=0;
			for(long irow=Ap[icol]; irow<Ap[icol+1]-1; irow++){
				sum+=Ax[irow]*P(y2, Ai[irow], iy);
			}
			/*assert(Ai[Ap[icol+1]-1]==icol);//confirm upper right triangular */
			P(y2, icol, iy)=(P(y2, icol, iy)-sum)*AxI;
		}
	}

	/*Solve R\y */
	for(long icol=A->nx-1; icol>-1; icol--){
		real AxI=1./Ax[Ap[icol+1]-1];
		for(long iy=start; iy<end; iy++){
			P(y2, icol, iy)*=AxI;
			real val=-P(y2, icol, iy);
			/*We do in reverse order to increase memory reuse. 1.5xFaster than forward order. */
			for(long irow=Ap[icol+1]-2; irow>Ap[icol]-1; irow--){
				P(y2, Ai[irow], iy)+=val*Ax[irow];
			}
		}
	}
}
//Solve from column start to end. *x and y can be the same.
static void chol_solve_each(dmat** x, spchol* A, const dmat* y, long start, long end){
	if(end==0) end=y->ny;
	if(A->L){
		chol_solve_cholmod(x, A, y, start, end);
	} else{
		dmat* y2=0;
		chol_perm(&y2, A->Cp, y, 0, start, end-start, 1);
		if(A->Cl){
			chol_solve_lower(A->Cl, y2, 0, end-start);
		} else if(A->Cu){
			chol_solve_upper(A->Cu, y2, 0, end-start);
		} else{
			error("There is no valid Cholesky factor\n");
		}
		chol_perm(x, A->Cp, y2, start, 0, end-start, -1);
		dfree(y2);
	}
}
/*
  Solve section of the columns in multi-threaded way.
*/
static void* chol_solve_thread(thread_t* info){
	CHOLSOLVE_T* data=(CHOLSOLVE_T*)info->data;
	chol_solve_each(&data->x, data->A, data->y, info->start, info->end);
	return NULL;
}
/**
   Solve A*x=y where the cholesky factor of A is stored in A.
*/
void chol_solve(dmat** x, spchol* A, dmat* y){
	if(y->ny==1){
		chol_solve_each(x, A, y, 0, 0);
	} else{//Multi-threaded call
		if(!*x) *x=dnew(y->nx, y->ny);
		CHOLSOLVE_T data={*x,A,y};
		thread_t *tdata=thread_prep(0, y->ny, NTHREAD, chol_solve_thread, &data);
		CALL_THREAD(tdata, 1);
		free(tdata);
	}
	if(0){
		static int count=-1; count++;
		writebin(*x, "chol_solve_x_%d", count);
		writebin(y, "chol_solve_y_%d", count);
	}
}
/**
   Free cholesky factor.*/
void chol_free_do(spchol* A){
	if(A){
		if(A->L){
			MOD(free_factor)(&A->L, A->c);
			MOD(finish)(A->c);
			free(A->c);
		} else{
			free(A->Cp);
		}
		if(A->Cl) dspfree(A->Cl);
		if(A->Cu) dspfree(A->Cu);
		free(A);
	}
}
