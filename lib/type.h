/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_TYPE_H
#define AOS_TYPE_H
#include "common.h"
#include <fftw3.h>
#include "misc.h"
#if defined(DLONG)
typedef long spint; //we use 32 bit in sparse index.
#else
typedef int spint;
#endif
/**
   \file lib/type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.  */
enum{
    MT_NORMAL=0,
    MT_REF,
    MT_MMAP
};
/**
   An arrays of 1-d plans that are used to do 2-d FFTs only over specified region.
 */
typedef struct PLAN1D_T{
    int ncomp;         /**< For a NxN array, only convert center ncomp*ncomp to Fourier space. */
    fftw_plan plan[3]; /**< Array of plans for 1-d FFT */
}PLAN1D_T;

/*
  Separate definition of struct with typedef. Put definition of struct in a
private file and typedef in a public interface to hide the struct from the user
*/

/**
   a double matrix object contains 2-d array of double numbers
 */
typedef struct dmat{
    double *p;  /**< the pointer to allocated memory. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    long *nref; /**< reference count */
    long type;  /**< specify whether this is allocated or in shared memory. */
}dmat;

/**
   a double complex matrix object contains 2-d arrays of
   double complex numbers. */
typedef struct cmat{
    dcomplex *p; /**< the pointer to allocated memory. */
    long nx;     /**< number of rows */
    long ny;     /**< number of columns */
    long *nref;  /**< reference count */
    long type;   /**< specify whether this is allocated or in shared memory. */
    fftw_plan  plan[3];  /**< Stores array of FFT plans for forward and backward FFT. */
    PLAN1D_T *plan1d[3];/**< Stores array of 1-D FFT plans. See PLAN1D_T */
}cmat;

/**
   an 2-d block matrix of cmat.
 */
typedef struct ccell{
    cmat **p;   /**< Contains an array of pointers to cmat. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    void *mmap; /**< not NULL if mmaped. unmap the mmaped memory*/
}ccell;
/**
   an 2-d block matrix of dmat.
 */
typedef struct dcell{
    dmat **p;   /**< Contains an array of pointers to dmat. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    void *mmap; /**< not NULL if mmaped. unmap the mmaped memory*/
}dcell;

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;
/**
   a sparse array of double numbers stored in
   compressed column format, i.e. MATLAB format */
typedef struct dsp{
    long nzmax ;/**< maximum number of entries */
    long m ;	/**< number of rows */
    long n ;	/**< number of columns */
    spint *p ;   /**< column pointers (size n+1) or col indlces (size nzmax) when nz!=-1 */
    spint *i ;   /**< row indices, size nzmax */
    double *x ;	/**< numerical values, size nzmax */
    long nz ;   /**< number of entries in triplet matrix, -1 for compressed-col */
    long *nref; /**< reference counting like dmat */
}dsp;
/**
   an 2-d array of sparse.
 */
typedef struct spcell{
    dsp **p;    /**< Contains an array of pointers to dsp. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
}spcell;
/**
   a sparse array of double complex numbers stored in
   compressed column format */
typedef struct csp{
    long nzmax ;/**< maximum number of entries */
    long m ;	/**< number of rows */
    long n ;	/**< number of columns */
    spint *p ;   /**< column pointers (size n+1) or col indlces (size nzmax) */
    spint *i ;   /**< row indices, size nzmax */
    dcomplex *x;/**< numerical values, size nzmax */
    long nz ;   /**< # of entries in triplet matrix, -1 for compressed-col */
    long *nref; /**< reference counting like cmat*/
}csp;
/**
   an 2-d array of csp.
 */
typedef struct cspcell{
    csp **p;    /**< Contains an array of pointers to csp. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
}cspcell;

#define PDMAT(M,P)   double (*restrict P)[(M)->nx]=(double(*)[(M)->nx])(M)->p
#define PCMAT(M,P) dcomplex (*restrict P)[(M)->nx]=(dcomplex(*)[(M)->nx])(M)->p
#define PDCELL(M,P)   dmat* (*restrict P)[(M)->nx]=(dmat*(*)[(M)->nx])(M)->p
#define PCCELL(M,P)   cmat* (*restrict P)[(M)->nx]=(cmat*(*)[(M)->nx])(M)->p
#define PSPCELL(M,P)   dsp* (*restrict P)[(M)->nx]=(dsp*(*)[(M)->nx])(M)->p
#define PDSPCELL(M,P)  dsp* (*restrict P)[(M)->nx]=(dsp*(*)[(M)->nx])(M)->p
#define PCSPCELL(M,P)  csp* (*restrict P)[(M)->nx]=(csp*(*)[(M)->nx])(M)->p

#define dfree(A) ({dfree_do(A,0);A=NULL;})
#define cfree(A) ({cfree_do(A,0);A=NULL;})

#define spfree(A)      {spfree_do(A); A=NULL;}
#define spcellfree(A)  {spcellfree_do(A); A=NULL;}
#define cspfree(A)     {cspfree_do(A); A=NULL;}
#define cspcellfree(A) {cspcellfree_do(A); A=NULL;}

#define dcp2(A,B)  memcpy(A->p,B->p,sizeof(double)*A->nx*A->ny);
#define dcellfree(A) ({dcellfree_do(A);A=NULL;})
#define ccellfree(A) ({ccellfree_do(A);A=NULL;})
#define dcellfreearr(A,n) ({for(int in=0; in<n; in++){dcellfree(A[in]);};free(A);})

#define cabs2(A) (pow(creal(A),2)+pow(cimag(A),2))

#define AOS_CMAT(A) c##A
#define AOS_CSPARSE(A) c##A

#define AOS_DMAT(A) d##A
#define AOS_SPARSE(A) A
#endif
