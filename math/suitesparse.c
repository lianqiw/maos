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
/*Obtain from SuiteSparse cs_multiply. Gather together */
typedef spint CS_INT;

#ifdef USE_COMPLEX
#ifdef USE_SINGLE
#define cs zsp
#define CS_ENTRY fcomplex
#define CS_REAL(x) crealf(x)
#define CS_IMAG(x) cimagf(x)
#define CS_CONJ(x) conjf(x)
#define CS_ABS(x) cabsf(x)
#define CS_MAX(a,b) ((cabsf(a) > cabsf(b)) ? (a) : (b))
#define CS_MIN(a,b) ((cabsf(a) < cabsf(b)) ? (a) : (b))
#else
#define cs csp
#define CS_ENTRY dcomplex
#define CS_REAL(x) creal(x)
#define CS_IMAG(x) cimag(x)
#define CS_CONJ(x) conj(x)
#define CS_ABS(x) cabs(x)
#define CS_MAX(a,b) ((cabs(a) > cabs(b)) ? (a) : (b))
#define CS_MIN(a,b) ((cabs(a) < cabs(b)) ? (a) : (b))
#endif
#else
#ifdef USE_SINGLE
#define cs ssp
#define CS_ENTRY float
#define CS_REAL(x) (x)
#define CS_IMAG(x) (0.)
#define CS_CONJ(x) (x)
#define CS_ABS(x) fabsf(x)
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#else
#define cs dsp
#define CS_ENTRY double
#define CS_REAL(x) (x)
#define CS_IMAG(x) (0.)
#define CS_CONJ(x) (x)
#define CS_ABS(x) fabs(x)
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#endif

#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

/**
 wrapper for free */
static void *cs_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/**
 free a sparse matrix */
static cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    cs_free (A->p) ;
    cs_free (A->i) ;
    cs_free (A->x) ;
    return (cs_free (A)) ;      /* free the cs struct and return NULL */
}

/**
 wrapper for malloc */
static void *cs_malloc (CS_INT n, size_t size)
{
    return (malloc (CS_MAX (n,1) * size)) ;
}

/**
 wrapper for calloc */
static void *cs_calloc (CS_INT n, size_t size)
{
    return (calloc (CS_MAX (n,1), size)) ;
}
/**
 wrapper for realloc */
static void *cs_realloc (void *p, CS_INT n, size_t size, CS_INT *ok)
{
    void *pnew ;
    pnew = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

/**
 allocate a sparse matrix (triplet form or compressed-column form) */
static cs *cs_spalloc (CS_INT m, CS_INT n, CS_INT nzmax, CS_INT values, CS_INT triplet)
{
    cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = cs_malloc (triplet ? nzmax : n+1, sizeof (CS_INT)) ;
    A->i = cs_malloc (nzmax, sizeof (CS_INT)) ;
    A->x = values ? cs_malloc (nzmax, sizeof (CS_ENTRY)) : NULL ;
    A->nref=calloc(1,sizeof(long));
    A->nref[0]=1;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/**
 change the max # of entries sparse matrix */
static CS_INT cs_sprealloc (cs *A, CS_INT nzmax)
{
    CS_INT ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;
    A->i = cs_realloc (A->i, nzmax, sizeof (CS_INT), &oki) ;
    if (CS_TRIPLET (A)) A->p = cs_realloc (A->p, nzmax, sizeof (CS_INT), &okj) ;
    if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (CS_ENTRY), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/**
 free workspace and return a sparse matrix result */
static cs *cs_done (cs *C, void *w, void *x, CS_INT ok)
{
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? C : cs_spfree (C)) ;   /* return result if OK, else free it */
}

/**
 x = x + beta * A(:,j), where x is a dense vector and A(:,j) is dsp */
static CS_INT cs_scatter (const cs *A, CS_INT j, CS_ENTRY beta, CS_INT *w, CS_ENTRY *x, CS_INT mark,
    cs *C, CS_INT nz)
{
    CS_INT i, p, *Ap, *Ai, *Ci ;
    CS_ENTRY *Ax ;
    if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
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
cs* Y(cs_multiply) (const cs *A, const cs *B)
{
    CS_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    CS_ENTRY *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    w = cs_calloc (m, sizeof (CS_INT)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? cs_malloc (m, sizeof (CS_ENTRY)) : NULL ; /* get workspace */
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (cs_done (C, w, x, 0)) ;             /* out of memory */
        } 
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 C = alpha*A + beta*B */
cs* Y(cs_add) (const cs *A, const cs *B, CS_ENTRY alpha, CS_ENTRY beta)
{
    CS_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    CS_ENTRY *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;         /* check inputs */
    if (A->m != B->m || A->n != B->n) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = cs_calloc (m, sizeof (CS_INT)) ;                       /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? cs_malloc (m, sizeof (CS_ENTRY)) : NULL ;    /* get workspace */
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/**
 drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
static CS_INT cs_fkeep (cs *A, CS_INT (*fkeep) (CS_INT, CS_INT, CS_ENTRY, void *), void *other)
{
    CS_INT j, p, nz = 0, n, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    if (!CS_CSC (A) || !fkeep) return (-1) ;    /* check inputs */
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
    cs_sprealloc (A, 0) ;                   /* remove extra space from A */
    return (nz) ;
}
/**
   return true of is not zero.
*/
static CS_INT cs_nonzero (CS_INT i, CS_INT j, CS_ENTRY aij, void *other)
{
    (void) i;
    (void) j;
    (void) other;
    return CS_ABS(aij)>1e-50 ;
}
/**
   drop zeros in the sparse matrix.
 */
CS_INT Y(cs_dropzeros) (cs *A)
{
    return (cs_fkeep (A, &cs_nonzero, NULL)) ;  /* keep all nonzero entries */
}
/**
   whether value is below threshold.
 */
static CS_INT cs_tol (CS_INT i, CS_INT j, CS_ENTRY aij, void *tol)
{
    (void) i;
    (void) j;
    return (CS_ABS (aij) > *((double *) tol)) ;
}
/**
   drop values below threashold of tol.
*/ 
CS_INT Y(cs_droptol) (cs *A, double tol)
{
    return (cs_fkeep (A, &cs_tol, &tol)) ;    /* keep all large entries */
}
static double cs_cumsum (CS_INT *p, CS_INT *c, CS_INT n)
{
    CS_INT i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid CS_INT overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

cs* Y(cs_transpose) (const cs *A, CS_INT values)
{
    CS_INT p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    CS_ENTRY *Cx, *Ax ;
    cs *C ;
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = cs_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    w = cs_calloc (m, sizeof (CS_INT)) ;                      /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cs_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = (values > 0) ? CS_CONJ (Ax [p]) : Ax [p] ;
        }
    }
    return (cs_done (C, w, NULL, 1)) ;  /* success; free w and return C */
}
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
