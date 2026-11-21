/*
  Copyright 2009-2026 Lianqi Wang
  
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
  Adapted from CSparse

  CSparse: a Concise Sparse Matrix package.
  Copyright (c) 2006-2017, Timothy A. Davis.
  http://www.suitesparse.com

  --------------------------------------------------------------------------------

  CSparse is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  CSparse is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this Module; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

/*
  Sort the columns in ascending order before exporting to matlab.
*/


#include <limits.h>
#include "mathdef.h"
#include "defs.h"
/*Obtain from SuiteSparse ss_multiply. Gather together */
#define SS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SS_FLIP(i) (-(i)-2)
#define SS_UNFLIP(i) (((i) < 0) ? SS_FLIP(i) : (i))
#define SS_MARKED(w,j) (w [j] < 0)
#define SS_MARK(w,j) { w [j] = SS_FLIP (w [j]) ; }

/**
   wrapper for malloc. Allocates at least one element.

   @param[in] n Number of elements.
   @param[in] size Size of each element in bytes.
   @return Newly allocated memory; the caller frees it.
*/
static void* ss_malloc(spint n, size_t size){
	return (malloc(SS_MAX(n, 1)*size));
}

/**
   wrapper for calloc. Allocates at least one element.

   @param[in] n Number of elements.
   @param[in] size Size of each element in bytes.
   @return Newly allocated zeroed memory; the caller frees it.
*/
static void* ss_calloc(spint n, size_t size){
	return (calloc(SS_MAX(n, 1), size));
}
/**
   wrapper for realloc. Allocates at least one element.

   @param[in] p Memory to reallocate.
   @param[in] n Number of elements.
   @param[in] size Size of each element in bytes.
   @param[out] ok Set to 1 on success and 0 on failure.
   @return The reallocated memory, or the original p on failure.
*/
static void* ss_realloc(void* p, spint n, size_t size, spint* ok){
	void* pnew;
	pnew=realloc(p, SS_MAX(n, 1)*size); /* realloc the block */
	*ok=(pnew!=NULL);                  /* realloc fails if pnew is NULL */
	return (pnew?pnew:p);             /* return original p if failure */
}

/**
   change the max # of entries sparse matrix

   @param[in,out] A Sparse matrix to resize.
   @param[in] nzmax New maximum number of entries. If <=0, use the actual number of entries.
   @return 1 on success, 0 on failure or if A is NULL.
*/
static spint ss_sprealloc(X(sp)* A, spint nzmax){
	spint ok, oki, okj=1, okx=1;
	if(!A) return (0);
	if(nzmax<=0) nzmax=(A->pp[A->ny]);
	A->pi=(spint*)ss_realloc(A->pi, nzmax, sizeof(spint), &oki);
	if(A->px) A->px=(T*)ss_realloc(A->px, nzmax, sizeof(T), &okx);
	ok=(oki&&okj&&okx);
	if(ok) A->nzmax=nzmax;
	return (ok);
}

/**
   free workspace and return a sparse matrix result

   @param[in] C The result sparse matrix. Freed if ok is 0.
   @param[in] w Workspace to free.
   @param[in] x Workspace to free.
   @param[in] ok Whether the operation succeeded.
   @return C if ok, otherwise NULL.
*/
static X(sp)* ss_done(X(sp)* C, void* w, void* x, spint ok){
	free(w);                       /* free workspace */
	free(x);
	return (ok?C:(X(spfree_do)(C), (X(sp*))NULL));   /* return result if OK, else free it */
}

/**
   x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse. New
   row indices are appended to the pattern of C.

   @param[in] A Sparse matrix.
   @param[in] j Column of A.
   @param[in] beta Scaling factor.
   @param[in,out] w Workspace marking rows already in the current column of C.
   @param[in,out] x Dense accumulator. May be NULL to compute the pattern only.
   @param[in] mark Mark value of the current column.
   @param[in,out] C Result matrix whose row indices are updated.
   @param[in] nz Current number of entries in C.
   @return The updated number of entries in C, or -1 if inputs are invalid.
*/
static spint ss_scatter(const X(sp)* A, spint j, T beta, spint* w, T* x, spint mark, X(sp)* C, spint nz){
	spint i, p, * Ap, * Ai, * Ci;
	T* Ax;
	if(!(A)||!w||!(C)) return (-1);     /* check inputs */
	Ap=A->pp; Ai=A->pi; Ax=A->px; Ci=C->pi;
	for(p=Ap[j]; p<Ap[j+1]; p++){
		i=Ai[p];                            /* A(i,j) is nonzero */
		if(w[i]<mark){
			w[i]=mark;                      /* i is new entry in column j */
			Ci[nz++]=i;                     /* add i to pattern of C(:,j) */
			if(x) x[i]=beta*Ax[p];      /* x(i) = beta*A(i,j) */
		} else if(x) x[i]+=beta*Ax[p];    /* i exists in C(:,j) already */
	}
	return (nz);
}

/**
   Sparse matrix multiplication C = A*B.

   @param[in] A Left sparse matrix.
   @param[in] B Right sparse matrix.
   @return Newly allocated sparse matrix C; the caller frees it.
*/
X(sp)* X(ss_multiply) (const X(sp)* A, const X(sp)* B){
	spint p, nz=0, anz, * Cp, * Ci, * Bp, bnz, * w, values, * Bi;
	long m, n, j;
	T* x, * Bx, * Cx;
	X(sp)* C;
	if(!(A)||!(B)){
		error("A or B is not available\n");
		return (NULL);      /* check inputs */
	}
	if(A->ny!=B->nx){
		error("Matrix mismatch: A->ny=%ld != B->nx=%ld\n", A->ny, B->nx);
		return (NULL);
	}
	m=A->nx; anz=A->pp[A->ny];
	n=B->ny; Bp=B->pp; Bi=B->pi; Bx=B->px; bnz=Bp[n];
	w=(spint*)ss_calloc(m, sizeof(spint));                    /* get workspace */
	values=(A->px!=NULL)&&(Bx!=NULL);
	x=values?(T*)ss_malloc(m, sizeof(T)):NULL; /* get workspace */
	C=X(spnew)(m, n, anz+bnz);        /* allocate result */
	if(!C||!w||(values&&!x)){
		error("Out of memory\n");
		return (ss_done(C, w, x, 0));
	}
	Cp=C->pp;
	for(j=0; j<n; j++){
		if(nz+m>C->nzmax&&!ss_sprealloc(C, 2*(C->nzmax)+m)){
			error("Out of memory. \n");
			return (ss_done(C, w, x, 0));             /* out of memory */
		}
		Ci=C->pi; Cx=C->px;         /* C->pi and C->px may be reallocated */
		Cp[j]=nz;                   /* column j of C starts here */
		for(p=Bp[j]; p<Bp[j+1]; p++){
			nz=ss_scatter(A, Bi[p], Bx?Bx[p]:1, w, x, j+1, C, nz);
		}
		if(values) for(p=Cp[j]; p<nz; p++) Cx[p]=x[Ci[p]];
	}
	Cp[n]=nz;                       /* finalize the last column of C */
	ss_sprealloc(C, 0);               /* remove extra space from C */
	return (ss_done(C, w, x, 1));     /* success; free workspace, return C */
}

/**
   C = alpha*A + beta*B

   @param[in] A First sparse matrix.
   @param[in] B Second sparse matrix, same size as A.
   @param[in] alpha Scaling of A.
   @param[in] beta Scaling of B.
   @return Newly allocated sparse matrix C; the caller frees it.
*/
X(sp)* X(ss_add) (const X(sp)* A, const X(sp)* B, T alpha, T beta){
	spint p, j, nz=0, anz, * Cp, * Ci, * Bp, bnz, * w, values;
	long m, n; //to avoid overflow when multying
	T* x, * Bx, * Cx;
	X(sp)* C;
	if(!(A)||!(B)){
		error("A or B is not available\n");
		return (NULL);         /* check inputs */
	}
	if(A->nx!=B->nx||A->ny!=B->ny){
		error("Matrix mismatch\n");
		return (NULL);
	}
	m=A->nx; anz=A->pp[A->ny];
	n=B->ny; Bp=B->pp; Bx=B->px; bnz=Bp[n];
	w=(spint*)ss_calloc(m, sizeof(spint));                       /* get workspace */
	values=(A->px!=NULL)&&(Bx!=NULL);
	x=values?(T*)ss_malloc(m, sizeof(T)):NULL;    /* get workspace */
	spint cnz=anz+bnz;
	if(cnz>m*n){
		cnz=m*n;
	}
	C=X(spnew)(m, n, cnz);           /* allocate result*/
	if(!C||!w||(values&&!x)){
		error("Out of memory\n");
		return (ss_done(C, w, x, 0));
	}
	Cp=C->pp; Ci=C->pi; Cx=C->px;
	for(j=0; j<n; j++){
		Cp[j]=nz;                   /* column j of C starts here */
		nz=ss_scatter(A, j, alpha, w, x, j+1, C, nz);   /* alpha*A(:,j)*/
		nz=ss_scatter(B, j, beta, w, x, j+1, C, nz);    /* beta*B(:,j) */
		if(values) for(p=Cp[j]; p<nz; p++) Cx[p]=x[Ci[p]];
	}
	Cp[n]=nz;                       /* finalize the last column of C */
	ss_sprealloc(C, 0);               /* remove extra space from C */
	return (ss_done(C, w, x, 1));     /* success; free workspace, return C */
}

/**
   drop entries for which fkeep(A(i,j)) is false.

   @param[in,out] A Sparse matrix, modified in place.
   @param[in] fkeep Function returning nonzero for entries to keep.
   @param[in] other Extra data passed to fkeep.
   @return Number of remaining entries if OK, else -1.
*/
static spint ss_fkeep(X(sp)* A, spint(*fkeep) (spint, spint, T, void*), void* other){
	spint p, nz=0, * Ap, * Ai;
	long j, n;
	T* Ax;
	if(!(A)||!fkeep) return (-1);    /* check inputs */
	n=A->ny; Ap=A->pp; Ai=A->pi; Ax=A->px;
	for(j=0; j<n; j++){
		p=Ap[j];                        /* get current location of col j */
		Ap[j]=nz;                       /* record new location of col j */
		for(; p<Ap[j+1]; p++){
			if(fkeep(Ai[p], j, Ax?Ax[p]:1, other)){
				if(Ax) Ax[nz]=Ax[p];  /* keep A(i,j) */
				Ai[nz++]=Ai[p];
			}
		}
	}
	Ap[n]=nz;                           /* finalize A */
	ss_sprealloc(A, 0);                   /* remove extra space from A */
	return (nz);
}
/**
   return true if aij is not zero.

   @param[in] i Row index. Unused.
   @param[in] j Column index. Unused.
   @param[in] aij Value of the entry.
   @param[in] other Unused.
   @return 1 if abs(aij)>1e-50, 0 otherwise.
*/
static spint ss_nonzero(spint i, spint j, T aij, void* other){
	(void)i;
	(void)j;
	(void)other;
	return ABS(aij)>1e-50;
}
/**
   drop zeros in the sparse matrix.

   @param[in,out] A Sparse matrix, modified in place.
   @return Number of remaining entries, or -1 if A is NULL.
*/
spint X(ss_dropzeros) (X(sp)* A){
	return (ss_fkeep(A, &ss_nonzero, NULL));  /* keep all nonzero entries */
}
/**
   whether value is above threshold.

   @param[in] i Row index. Unused.
   @param[in] j Column index. Unused.
   @param[in] aij Value of the entry.
   @param[in] tol Pointer to the threshold (real).
   @return 1 if abs(aij)>*tol, 0 otherwise.
*/
static spint ss_tol(spint i, spint j, T aij, void* tol){
	(void)i;
	(void)j;
	return (ABS(aij)>*((real*)tol));
}
/**
   drop values whose magnitude is not above tol.

   @param[in,out] A Sparse matrix, modified in place.
   @param[in] tol Threshold.
   @return Number of remaining entries, or -1 if A is NULL.
*/
spint X(ss_droptol) (X(sp)* A, real tol){
	return (ss_fkeep(A, &ss_tol, &tol));    /* keep all large entries */
}
/**
   p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c

   @param[out] p Cumulative sum, of length n+1.
   @param[in,out] c Counts of length n. Replaced by p[0..n-1] on return.
   @param[in] n Number of elements in c.
   @return Sum of c[0..n-1] (as real to avoid overflow), or -1 if inputs are NULL.
*/
static real ss_cumsum(spint* p, spint* c, spint n){
	spint i, nz=0;
	real nz2=0;
	if(!p||!c) return (-1);     /* check inputs */
	for(i=0; i<n; i++){
		p[i]=nz;
		nz+=c[i];
		nz2+=c[i];              /* also in real to avoid spint overflow */
		c[i]=p[i];             /* also copy p[0..n-1] back into c[0..n-1]*/
	}
	p[n]=nz;
	return (nz2);                  /* return sum (c [0..n-1]) */
}

/**
   Return the transpose of sparse matrix A as a new matrix. If values>0, the
   values are conjugated (for complex); otherwise they are copied as is.

   @param[in] A Sparse matrix.
   @param[in] values Whether to conjugate the values.
   @return Newly allocated transpose; the caller frees it.
*/
X(sp)* X(ss_transpose) (const X(sp)* A, spint values){
	spint p, q, * Cp, * Ci, * Ap, * Ai, * w;
	long n, m, j;
	T* Cx, * Ax;
	X(sp)* C;
	if(!(A)) return (NULL);    /* check inputs */
	m=A->nx; n=A->ny; Ap=A->pp; Ai=A->pi; Ax=A->px;
	C=X(spnew)(n, m, Ap[n]);       /* allocate result */
	w=(spint*)ss_calloc(m, sizeof(spint));                      /* get workspace */
	if(!C||!w) return (ss_done(C, w, NULL, 0));       /* out of memory */
	Cp=C->pp; Ci=C->pi; Cx=C->px;
	for(p=0; p<Ap[n]; p++) w[Ai[p]]++;          /* row counts */
	ss_cumsum(Cp, w, m);                                 /* row pointers */
	for(j=0; j<n; j++){
		for(p=Ap[j]; p<Ap[j+1]; p++){
			Ci[q=w[Ai[p]]++]=j; /* place A(i,j) as entry C(j,i) */
			if(Cx) Cx[q]=(values>0)?CONJ(Ax[p]):Ax[p];
		}
	}
	return (ss_done(C, w, NULL, 1));  /* success; free w and return C */
}
