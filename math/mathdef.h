/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "matmath.h"
#include "sp.h"
#include "fft.h"
#include "matbin.h"
#include "spbin.h"
#include "cell.h"
#include "chol.h"
#include "zfarr.h"
#include "mathmisc.h"
#include "locbin.h"
#include "random.h"

#define AOS_LMAT(A) l##A
#define AOS_CMAT(A) c##A
#define AOS_DMAT(A) d##A
#define AOS_SMAT(A) s##A
#define AOS_ZMAT(A) z##A
#define isempty(A) (!(A) || !(A)->nx || !(A)->ny)
//#ifdef USE_SINGLE
//Single
AOS_MAT_DEF(AOS_SMAT,float)
AOS_MAT_DEF(AOS_ZMAT,fcomplex)

AOS_MATMATH_DEF(AOS_SMAT,AOS_SMAT,float,float)
AOS_MATMATH_DEF(AOS_ZMAT,AOS_SMAT,fcomplex,float)

AOS_CMATMATH_DEF(AOS_ZMAT,AOS_SMAT,fcomplex,float)

AOS_MATBIN_DEF(AOS_SMAT,float)
AOS_MATBIN_DEF(AOS_ZMAT,fcomplex)


AOS_SP_DEF(AOS_SMAT,float,float,fcomplex)
AOS_SP_DEF(AOS_ZMAT,fcomplex,float,fcomplex)

AOS_SPBIN_DEF(AOS_SMAT,float)
AOS_SPBIN_DEF(AOS_ZMAT,fcomplex)

AOS_FFT_DEF(AOS_SMAT)
AOS_FFT_DEF(AOS_ZMAT)
//#else
//Double
AOS_MAT_DEF(AOS_DMAT,double)
AOS_MAT_DEF(AOS_CMAT,dcomplex)

AOS_MATMATH_DEF(AOS_DMAT,AOS_DMAT,double,double)
AOS_MATMATH_DEF(AOS_CMAT,AOS_DMAT,dcomplex,double)

AOS_CMATMATH_DEF(AOS_CMAT,AOS_DMAT,dcomplex,double)

AOS_MATBIN_DEF(AOS_DMAT,double)
AOS_MATBIN_DEF(AOS_CMAT,dcomplex)

AOS_SP_DEF(AOS_DMAT,double,double,dcomplex)
AOS_SP_DEF(AOS_CMAT,dcomplex,double,dcomplex)

AOS_SPBIN_DEF(AOS_DMAT,double)
AOS_SPBIN_DEF(AOS_CMAT,dcomplex)

AOS_FFT_DEF(AOS_DMAT)
AOS_FFT_DEF(AOS_CMAT)

//#endif
AOS_MAT_DEF(AOS_LMAT, long)
AOS_MATBIN_DEF(AOS_LMAT,long)

#define cabs2f(A)     (powf(crealf(A),2)+powf(cimagf(A),2))
#define cabs2(A)     (pow(creal(A),2)+pow(cimag(A),2))

#if __cplusplus >= 201103L
//Work around C++ restrictions
#define PALL(T,A,pp) const long pp##nx=(A)->nx; auto pp=(T(*)[pp##nx])(A)->p
#else 
//C99
#define PALL(T,A,pp) const long pp##nx=(A)->nx; T (*pp)[pp##nx]=(T(*)[pp##nx])(A)->p
#endif
#define PDMAT(M,P)   PALL(double,M,P)
#define PDCELL(M,P)  PALL(dmat*,M,P)
#define dfree(A)     ({dfree_do((A),0);(A)=NULL;})
#define dcp2(A,B)    memcpy(A->p,B->p,sizeof(double)*A->nx*A->ny)
#define dcellfree(A) ({cellfree_do(A);A=NULL;})
#define dcellfreearr(A,n) ({for(int in=0; A&&in<n; in++){dcellfree(A[in]);};free(A);A=NULL;})
#define dzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(double))
#define dhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(double), key)

#define PSMAT(M,P)   PALL(float,M,P)
#define PSCELL(M,P)  PALL(smat*,M,P)
#define sfree(A)     ({sfree_do((A),0);(A)=NULL;})
#define scp2(A,B)    memcpy(A->p,B->p,sizeof(float)*A->nx*A->ny)
#define scellfree(A) ({cellfree_do(A);A=NULL;})
#define scellfreearr(A,n) ({for(int in=0; A&&in<n; in++){scellfree(A[in]);};free(A);A=NULL;})
#define szero(A) if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(float))
#define shash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(float), key)

#define PCMAT(M,P)   PALL(dcomplex,M,P)
#define PCCELL(M,P)  PALL(cmat*,M,P)
#define cfree(A)     ({cfree_do(A,0);A=NULL;})
#define ccellfree(A) ({cellfree_do(A);A=NULL;})
#define ccellfreearr(A,n) ({for(int in=0; A&&in<n; in++){ccellfree(A[in]);};free(A);A=NULL;})
#define czero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(dcomplex))
#define chash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(dcomplex), key)

#define PZMAT(M,P)   PALL(fcomplex,M,P)
#define PZCELL(M,P)  PALL(zmat*,M,P) 
#define zfree(A)     ({zfree_do(A,0);A=NULL;})
#define zcellfree(A) ({cellfree_do(A);A=NULL;})
#define zzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(fcomplex))
#define zhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(fcomplex), key)

#define PLMAT(M,P)   PALL(long,M,P)
#define PLCELL(M,P)  PALL(lmat*,M,P) 
#define lfree(A)     ({lfree_do(A,0);A=NULL;})
#define lcellfree(A) ({cellfree_do(A);A=NULL;})
#define lzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(long))
#define lhash(A,key) hashlittle((A)->p, (A)->nx*(A)->ny*sizeof(long), key)

#define cellfree(A) ({cellfree_do(A); A=0;})

#define mapwrite(out, A...) write_by_id((void*)out, M_MAP64, A)
#define mapread(A...)    (map_t*)read_by_id(M_MAP64, 0, A)
#define mapcellread(A...) (mapcell*)read_by_id(M_MAP64, 1, A)

#define locwrite(out, A...) write_by_id((void*)out, M_LOC64, A)
#define locread(A...)    (loc_t*)read_by_id(M_LOC64, 0, A)
#define loccellread(A...) (loccell*)read_by_id(M_LOC64, 1, A)
/** Read needs type checking, so don't use readbin*/
#define dread(A...)    dmat_cast(read_by_id(M_DBL, 0, A))
#define dcellnew (dcell*)cellnew
#define dcellreaddata(fp, header) dcell_cast(readdata_by_id(fp, M_DBL, 1, header))
#define dcellread(A...) (dcell*)read_by_id(M_DBL, 1, A)
#define dccellread(A...) (dccell*)read_by_id(M_DBL, 2, A)
#define dcccellread(A...) (dcccell*)read_by_id(M_DBL, 3, A)

#define sread(A...)    smat_cast(read_by_id(M_FLT, 0, A))
#define scellnew (scell*)cellnew
#define scellreaddata(fp, header) scell_cast(readdata_by_id(fp, M_FLT, 1, header))
#define scellread(A...) (scell*)read_by_id(M_FLT, 1, A)
#define sccellread(A...) (sccell*)read_by_id(M_FLT, 2, A)
#define scccellread(A...) (scccell*)read_by_id(M_FLT, 3, A)

#define cread(A...)    cmat_cast(read_by_id(M_CMP, 0, A))
#define ccellnew (ccell*)cellnew
#define ccellreaddata(fp, header) ccell_cast(readdata_by_id(fp, M_CMP, 1, header))
#define ccellread(A...) (ccell*)read_by_id(M_CMP, 1, A)
#define cccellread(A...) (cccell*)read_by_id(M_CMP, 2, A)
#define ccccellread(A...) (ccccell*)read_by_id(M_CMP, 3, A)

#define zread(A...)    zmat_cast(read_by_id(M_ZMP, 0, A))
#define zcellnew (zcell*)cellnew
#define zcellreaddata(fp, header) zcell_cast(readdata_by_id(fp, M_ZMP, 1, header))
#define zcellread(A...) (zcell*)read_by_id(M_ZMP, 1, A)
#define zccellread(A...) (zccell*)read_by_id(M_ZMP, 2, A)
#define zcccellread(A...) (zcccell*)read_by_id(M_ZMP, 3, A)

#define lread(A...) lmat_cast(read_by_id(M_LONG, 0, A))
#define lcellnew (lcell*)cellnew
#define lcellreaddata(fp, header) lcell_cast(readdata_by_id(fp, M_LONG, 1, header))
#define lcellread(A...) (lcell*)read_by_id(M_LONG, 1, A)
#define lccellread(A...) (lccell*)read_by_id(M_LONG, 2, A)
#define lcccellread(A...) (lcccell*)read_by_id(M_LONG, 3, A)

#define dspread(A...) dsp_cast(read_by_id(M_DSP, 0, A))
#define dspcellread(A...) dspcell_cast(read_by_id(M_DSP, 1, A))
#define sspread(A...) ssp_cast(read_by_id(M_SSP, 0, A))
#define sspcellread(A...) sspcell_cast(read_by_id(M_SSP, 1, A))
#define cspread(A...) csp_cast(read_by_id(M_CSP, 0, A))
#define cspcellread(A...) cspcell_cast(read_by_id(M_CSP, 1, A))
#define zspread(A...) zsp_cast(read_by_id(M_ZSP, 0, A))
#define zspcellread(A...) zspcell_cast(read_by_id(M_ZSP, 1, A))
#endif
