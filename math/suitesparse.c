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

/*Taken from suitesparse package from UFL. */

/*my defines. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../sys/sys.h"
#include "type.h"
/*Obtain from SuiteSparse ss_multiply. Gather together */
typedef spint SS_INT;

#ifdef USE_COMPLEX
#ifdef USE_SINGLE
#define cs zsp
#define SS_ENTRY fcomplex
#define SS_REAL(x) crealf(x)
#define SS_IMAG(x) cimagf(x)
#define SS_CONJ(x) conjf(x)
#define SS_ABS(x) cabsf(x)
#define SS_MAX(a,b) ((cabsf(a) > cabsf(b)) ? (a) : (b))
#define SS_MIN(a,b) ((cabsf(a) < cabsf(b)) ? (a) : (b))
#else
#define cs csp
#define SS_ENTRY dcomplex
#define SS_REAL(x) creal(x)
#define SS_IMAG(x) cimag(x)
#define SS_CONJ(x) conj(x)
#define SS_ABS(x) cabs(x)
#define SS_MAX(a,b) ((cabs(a) > cabs(b)) ? (a) : (b))
#define SS_MIN(a,b) ((cabs(a) < cabs(b)) ? (a) : (b))
#endif
#else
#ifdef USE_SINGLE
#define cs ssp
#define SS_ENTRY float
#define SS_REAL(x) (x)
#define SS_IMAG(x) (0.)
#define SS_CONJ(x) (x)
#define SS_ABS(x) fabsf(x)
#define SS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#else
#define cs dsp
#define SS_ENTRY double
#define SS_REAL(x) (x)
#define SS_IMAG(x) (0.)
#define SS_CONJ(x) (x)
#define SS_ABS(x) fabs(x)
#define SS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#endif

#define SS_FLIP(i) (-(i)-2)
#define SS_UNFLIP(i) (((i) < 0) ? SS_FLIP(i) : (i))
#define SS_MARKED(w,j) (w [j] < 0)
#define SS_MARK(w,j) { w [j] = SS_FLIP (w [j]) ; }
#define SS_CSC(A) (A && (A->nz == -1))
#define SS_TRIPLET(A) (A && (A->nz >= 0))

/**
 wrapper for free */
static void *ss_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of ss_free */
}

/**
 free a sparse matrix */
static cs *ss_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    ss_free (A->p) ;
    ss_free (A->i) ;
    ss_free (A->x) ;
    return (ss_free (A)) ;      /* free the cs struct and return NULL */
}

/**
 wrapper for malloc */
static void *ss_malloc (SS_INT n, size_t size)
{
    return (malloc (SS_MAX (n,1) * size)) ;
}

/**
 wrapper for calloc */
static void *ss_calloc (SS_INT n, size_t size)
{
    return (calloc (SS_MAX (n,1), size)) ;
}
/**
 wrapper for realloc */
static void *ss_realloc (void *p, SS_INT n, size_t size, SS_INT *ok)
{
    void *pnew ;
    pnew = realloc (p, SS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

/**
 allocate a sparse matrix (triplet form or compressed-column form) */
static cs *ss_spalloc (SS_INT m, SS_INT n, SS_INT nzmax, SS_INT values, SS_INT triplet)
{
    cs *A = ss_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = SS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = ss_malloc (triplet ? nzmax : n+1, sizeof (SS_INT)) ;
    A->i = ss_malloc (nzmax, sizeof (SS_INT)) ;
    A->x = values ? ss_malloc (nzmax, sizeof (SS_ENTRY)) : NULL ;
    A->nref=calloc(1,sizeof(long));
    A->nref[0]=1;
    return ((!A->p || !A->i || (values && !A->x)) ? ss_spfree (A) : A) ;
}

/**
 change the max # of entries sparse matrix */
static SS_INT ss_sprealloc (cs *A, SS_INT nzmax)
{
    SS_INT ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (SS_CSC (A)) ? (A->p [A->n]) : A->nz ;
    A->i = ss_realloc (A->i, nzmax, sizeof (SS_INT), &oki) ;
    if (SS_TRIPLET (A)) A->p = ss_realloc (A->p, nzmax, sizeof (SS_INT), &okj) ;
    if (A->x) A->x = ss_realloc (A->x, nzmax, sizeof (SS_ENTRY), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/**
 free workspace and return a sparse matrix result */
static cs *ss_done (cs *C, void *w, void *x, SS_INT ok)
{
    ss_free (w) ;                       /* free workspace */
    ss_free (x) ;
    return (ok ? C : ss_spfree (C)) ;   /* return result if OK, else free it */
}

/**
 x = x + beta * A(:,j), where x is a dense vector and A(:,j) is dsp */
static SS_INT ss_scatter (const cs *A, SS_INT j, SS_ENTRY beta, SS_INT *w, SS_ENTRY *x, SS_INT mark,
    cs *C, SS_INT nz)
{
    SS_INT i, p, *Ap, *Ai, *Ci ;
    SS_ENTRY *Ax ;
    if (!SS_CSC (A) || !w || !SS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

/**
 C = A*B */
cs* Y(ss_multiply) (const cs *A, const cs *B)
{
    SS_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    SS_ENTRY *x, *Bx, *Cx ;
    cs *C ;
    if (!SS_CSC (A) || !SS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    w = ss_calloc (m, sizeof (SS_INT)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? ss_malloc (m, sizeof (SS_ENTRY)) : NULL ; /* get workspace */
    C = ss_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    if (!C || !w || (values && !x)) return (ss_done (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !ss_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (ss_done (C, w, x, 0)) ;             /* out of memory */
        } 
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = ss_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    ss_sprealloc (C, 0) ;               /* remove extra space from C */
    return (ss_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 C = alpha*A + beta*B */
cs* Y(ss_add) (const cs *A, const cs *B, SS_ENTRY alpha, SS_ENTRY beta)
{
    SS_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    SS_ENTRY *x, *Bx, *Cx ;
    cs *C ;
    if (!SS_CSC (A) || !SS_CSC (B)) return (NULL) ;         /* check inputs */
    if (A->m != B->m || A->n != B->n) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = ss_calloc (m, sizeof (SS_INT)) ;                       /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? ss_malloc (m, sizeof (SS_ENTRY)) : NULL ;    /* get workspace */
    C = ss_spalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    if (!C || !w || (values && !x)) return (ss_done (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = ss_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = ss_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    ss_sprealloc (C, 0) ;               /* remove extra space from C */
    return (ss_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
static SS_INT ss_fkeep (cs *A, SS_INT (*fkeep) (SS_INT, SS_INT, SS_ENTRY, void *), void *other)
{
    SS_INT j, p, nz = 0, n, *Ap, *Ai ;
    SS_ENTRY *Ax ;
    if (!SS_CSC (A) || !fkeep) return (-1) ;    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    Ap [n] = nz ;                           /* finalize A */
    ss_sprealloc (A, 0) ;                   /* remove extra space from A */
    return (nz) ;
}
/**
   return true of is not zero.
*/
static SS_INT ss_nonzero (SS_INT i, SS_INT j, SS_ENTRY aij, void *other)
{
    (void) i;
    (void) j;
    (void) other;
    return SS_ABS(aij)>1e-50 ;
}
/**
   drop zeros in the sparse matrix.
 */
SS_INT Y(ss_dropzeros) (cs *A)
{
    return (ss_fkeep (A, &ss_nonzero, NULL)) ;  /* keep all nonzero entries */
}
/**
   whether value is below threshold.
 */
static SS_INT ss_tol (SS_INT i, SS_INT j, SS_ENTRY aij, void *tol)
{
    (void) i;
    (void) j;
    return (SS_ABS (aij) > *((double *) tol)) ;
}
/**
   drop values below threashold of tol.
*/ 
SS_INT Y(ss_droptol) (cs *A, double tol)
{
    return (ss_fkeep (A, &ss_tol, &tol)) ;    /* keep all large entries */
}
static double ss_cumsum (SS_INT *p, SS_INT *c, SS_INT n)
{
    SS_INT i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid SS_INT overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

cs* Y(ss_transpose) (const cs *A, SS_INT values)
{
    SS_INT p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    SS_ENTRY *Cx, *Ax ;
    cs *C ;
    if (!SS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = ss_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    w = ss_calloc (m, sizeof (SS_INT)) ;                      /* get workspace */
    if (!C || !w) return (ss_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    ss_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = (values > 0) ? SS_CONJ (Ax [p]) : Ax [p] ;
        }
    }
    return (ss_done (C, w, NULL, 1)) ;  /* success; free w and return C */
}
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
