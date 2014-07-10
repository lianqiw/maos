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

#ifndef AOS_LIB_MATH_H
#define AOS_LIB_MATH_H
#include "../sys/sys.h"
#include "type.h"
#include "mat.h"
#include "sp.h"
#include "fft.h"
#include "matbin.h"
#include "spbin.h"
#include "cell.h"
#include "chol.h"
#include "imat.h"
#include "cellarr.h"
#include "mathmisc.h"
#include "loc.h"
#include "locbin.h"
#include "random.h"

#define cabs2f(A)     (powf(crealf(A),2)+powf(cimagf(A),2))
#define cabs2(A)     (pow(creal(A),2)+pow(cimag(A),2))

#ifdef __clang__
/*clang doesnot accept the gcc version in C++ mode*/
#define PALL(T,A,pp) typedef T pp##_ptr[(A)->nx]; pp##_ptr *pp=(pp##_ptr*)(A)->p
#else 
/*clang version caused <anonymous> must be uninitailized error in cuda code.*/
#define PALL(T,A,pp) T (*pp)[((cell*)A)->nx]=(T(*)[((cell*)A)->nx])((cell*)A)->p
#endif
#define PDMAT(M,P)   PALL(double,M,P)
#define PDCELL(M,P)  PALL(dmat*,M,P)
#define dfree(A)     ({dfree_do((A),0);(A)=NULL;})
#define dcp2(A,B)    memcpy(A->p,B->p,sizeof(double)*A->nx*A->ny)
#define dcellfree(A) ({dcellfree_do(A);A=NULL;})
#define dcellfreearr(A,n) ({for(int in=0; A&&in<n; in++){dcellfree(A[in]);};free(A);A=NULL;})
#define dzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(double))
#define dhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(double), key)

#define PSMAT(M,P)   PALL(float,M,P)
#define PSCELL(M,P)  PALL(smat*,M,P)
#define sfree(A)     ({sfree_do((A),0);(A)=NULL;})
#define scp2(A,B)    memcpy(A->p,B->p,sizeof(float)*A->nx*A->ny)
#define scellfree(A) ({scellfree_do(A);A=NULL;})
#define scellfreearr(A,n) ({for(int in=0; A&&in<n; in++){scellfree(A[in]);};free(A);A=NULL;})
#define szero(A) if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(float))
#define shash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(float), key)

#define PCMAT(M,P)   PALL(dcomplex,M,P)
#define PCCELL(M,P)  PALL(cmat*,M,P)
#define cfree(A)     ({cfree_do(A,0);A=NULL;})
#define ccellfree(A) ({ccellfree_do(A);A=NULL;})
#define ccellfreearr(A,n) ({for(int in=0; A&&in<n; in++){ccellfree(A[in]);};free(A);A=NULL;})
#define czero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(dcomplex))
#define chash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(dcomplex), key)

#define PZMAT(M,P)   PALL(fcomplex,M,P)
#define PZCELL(M,P)  PALL(zmat*,M,P) 
#define zfree(A)     ({zfree_do(A,0);A=NULL;})
#define zcellfree(A) ({zcellfree_do(A);A=NULL;})
#define zzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(fcomplex))
#define zhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(fcomplex), key)

#define AOS_CMAT(A) c##A
#define AOS_CSP(A)  c##A
#define AOS_DMAT(A) d##A
#define AOS_DSP(A)  A
#define AOS_SSP(A)  s##A
#define AOS_SMAT(A) s##A
#define AOS_ZMAT(A) z##A
#define AOS_ZSP(A)  z##A

//#ifdef USE_SINGLE
//Single
AOS_MAT_DEF(AOS_SMAT,AOS_SMAT,AOS_SSP,float,float)
AOS_MAT_DEF(AOS_ZMAT,AOS_SMAT,AOS_ZSP,fcomplex,float)
AOS_CMAT_DEF(AOS_ZMAT,AOS_SMAT,AOS_ZSP,fcomplex,float)

AOS_MATBIN_DEF(AOS_SMAT,AOS_SSP,float)
AOS_MATBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

//AOS_CELL_DEF(AOS_SMAT,AOS_SSP,float,float)
//AOS_CELL_DEF(AOS_ZMAT,AOS_ZSP,fcomplex,float)

AOS_SP_DEF(AOS_SMAT,AOS_SSP,float,float,fcomplex)
AOS_SP_DEF(AOS_ZMAT,AOS_ZSP,fcomplex,float,fcomplex)

AOS_SPBIN_DEF(AOS_SMAT,AOS_SSP,float)
AOS_SPBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

AOS_FFT_DEF(AOS_SMAT)
AOS_FFT_DEF(AOS_ZMAT)
//#else
//Double
AOS_MAT_DEF(AOS_DMAT,AOS_DMAT,AOS_DSP,double,double)
AOS_MAT_DEF(AOS_CMAT,AOS_DMAT,AOS_CSP,dcomplex,double)
AOS_MATBIN_DEF(AOS_DMAT,AOS_DSP,double)
AOS_MATBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)
AOS_CMAT_DEF(AOS_CMAT,AOS_DMAT,AOS_CSP,dcomplex,double)

//AOS_CELL_DEF(AOS_DMAT,AOS_DSP,double,double)
//AOS_CELL_DEF(AOS_CMAT,AOS_CSP,dcomplex,double)

AOS_SP_DEF(AOS_DMAT,AOS_DSP,double,double,dcomplex)
AOS_SP_DEF(AOS_CMAT,AOS_CSP,dcomplex,double,dcomplex)

AOS_SPBIN_DEF(AOS_DMAT,AOS_DSP,double)
AOS_SPBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)

AOS_FFT_DEF(AOS_DMAT)
AOS_FFT_DEF(AOS_CMAT)

//#endif

//#define dread(A...) readbin((READFUN)dreaddata, A)
#define mapread(A...) readbin((READFUN)mapreaddata, A)
#define locread(A...) readbin((READFUN)locreaddata, A)
#define mapwrite(out, A...) writebin((WRITEFUN)mapwritedata, out, A)
#define locwrite(out, A...) writebin((WRITEFUN)locwritedata, out, A)

#define mapcellread(A...) cellread(M_MAP64, A)
#define loccellread(A...) cellread(M_LOC64, A)
#define dccellread(A...) cellread(M_DBL, A)
#define dcccellread(A...) cellread(M_DBL, A)
#define cccellread(A...) cellread(M_CMP, A)
#define ccccellread(A...) cellread(M_CMP, A)

#endif
