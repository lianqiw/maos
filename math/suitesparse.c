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
 wrapper for malloc */
static void* ss_malloc(spint n, size_t size){
	return (malloc(SS_MAX(n, 1)*size));
}

/**
 wrapper for calloc */
static void* ss_calloc(spint n, size_t size){
	return (calloc(SS_MAX(n, 1), size));
}
/**
 wrapper for realloc */
static void* ss_realloc(void* p, spint n, size_t size, spint* ok){
	void* pnew;
	pnew=realloc(p, SS_MAX(n, 1)*size); /* realloc the block */
	*ok=(pnew!=NULL);                  /* realloc fails if pnew is NULL */
	return (pnew?pnew:p);             /* return original p if failure */
}

/**
 change the max # of entries sparse matrix */
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
 free workspace and return a sparse matrix result */
static X(sp)* ss_done(X(sp)* C, void* w, void* x, spint ok){
	free(w);                       /* free workspace */
	free(x);
	return (ok?C:(X(spfree_do)(C), (X(sp*))NULL));   /* return result if OK, else free it */
}

/**
 x = x + beta * A(:,j), where x is a dense vector and A(:,j) is dsp */
static spint ss_scatter(const X(sp)* A, spint j, T beta, spint* w, T* x, spint mark,
	X(sp)* C, spint nz){
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
 C = A*B */
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
 C = alpha*A + beta*B */
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
 drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
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
   return true of is not zero.
*/
static spint ss_nonzero(spint i, spint j, T aij, void* other){
	(void)i;
	(void)j;
	(void)other;
	return ABS(aij)>1e-50;
}
/**
   drop zeros in the sparse matrix.
 */
spint X(ss_dropzeros) (X(sp)* A){
	return (ss_fkeep(A, &ss_nonzero, NULL));  /* keep all nonzero entries */
}
/**
   whether value is below threshold.
 */
static spint ss_tol(spint i, spint j, T aij, void* tol){
	(void)i;
	(void)j;
	return (ABS(aij)>*((real*)tol));
}
/**
   drop values below threashold of tol.
*/
spint X(ss_droptol) (X(sp)* A, real tol){
	return (ss_fkeep(A, &ss_tol, &tol));    /* keep all large entries */
}
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
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
