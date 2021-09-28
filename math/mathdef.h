/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
   \file mathdef.h

   Includes macro definitions for all types.
*/
#ifndef AOS_LIB_MATH_H
#define AOS_LIB_MATH_H
#include "../sys/sys.h"
#include "type.h"
#include "dmath.h"
#include "cmath.h"
#include "smath.h"
#include "zmath.h"
#include "lmath.h"

#include "cell.h"
#include "chol.h"
#include "zfarr.h"
#include "mathmisc.h"
#include "loc.h"
#include "map.h"
#include "random.h"
/*
#define dnew(nx, ny) (funtrace_set,dnew_do(nx,ny,NULL,0));funtrace_unset
#define snew(nx, ny) (funtrace_set,snew_do(nx,ny,NULL,0));funtrace_unset
#define znew(nx, ny) (funtrace_set,znew_do(nx,ny,NULL,0));funtrace_unset
#define cnew(nx, ny) (funtrace_set,cnew_do(nx,ny,NULL,0));funtrace_unset
#define lnew(nx, ny) (funtrace_set,lnew_do(nx,ny,NULL,0));funtrace_unset
*/
#define isempty(A) (!(A) || !(A)->nx || !(A)->ny)

#define abs2(A)      ((A)*(A))
#define cabs2f(A)    (abs2(crealf(A))+abs2(cimagf(A)))
#define cabs2(A)     (abs2(creal(A))+abs2(cimag(A)))
/*!free a dmat and zero the pointer.*/
#define dfree(A)     if(A){dfree_do(A);A=NULL;}
#define dcellfree(A) if(A){cellfree_do((A)->base);A=NULL;}
#define dcellfreearr(A,n) if((A)&&(n)>0){for(int in=0; in<n; in++){dcellfree(A[in]);};free(A);A=NULL;}
/*!free a smat and zero the pointer.*/
#define sfree(A)     if(A){sfree_do(A);A=NULL;}
#define scellfree(A) if(A){cellfree_do((A)->base);A=NULL;}
#define scellfreearr(A,n) if((A)&&(n)>0){for(int in=0; A&&in<n; in++){scellfree(A[in]);};free(A);A=NULL;}
/*!free a cmat and zero the pointer.*/
#define cfree(A)     if(A){cfree_do(A);A=NULL;}
#define ccellfree(A) if(A){cellfree_do((A)->base);A=NULL;}
#define ccellfreearr(A,n) if((A)&&(n)>0){for(int in=0; A&&in<n; in++){ccellfree(A[in]);};free(A);A=NULL;}
/*!free a zmat and zero the pointer.*/
#define zfree(A)     if(A){zfree_do(A);A=NULL;}
#define zcellfree(A) if(A){cellfree_do((A)->base);A=NULL;}
/*!free a lmat and zero the pointer.*/
#define lfree(A)     if(A){lfree_do(A);A=NULL;}
#define lcellfree(A) if(A){cellfree_do((A)->base);A=NULL;}
/*!free a dsp and zero the pointer*/
#define dspfree(A)      if(A){dspfree_do(A); A=NULL;}
#define dspcellfree(A)  if(A){cellfree_do((A)->base); A=NULL;}
/*!free a ssp and zero the pointer*/
#define sspfree(A)      if(A){sspfree_do(A); A=NULL;}
#define sspcellfree(A)  if(A){cellfree_do((A)->base); A=NULL;}
/*!free a ssp and zero the pointer*/
#define cspfree(A)     if(A){cspfree_do(A); A=NULL;}
#define cspcellfree(A) if(A){cellfree_do((A)->base); A=NULL;}
/*!free a zsp and zero the pointer*/
#define zspfree(A)     if(A){zspfree_do(A); A=NULL;}
#define zspcellfree(A) if(A){cellfree_do((A)->base); A=NULL;}

#define mapwrite(out, A...) write_by_id(out?out->base:NULL, M_MAP, A)
#define mapread(A...)    (map_t*)read_by_id(M_MAP, 0, A)

#define mapcellread(A...) (mapcell*)read_by_id(M_MAP, 1, A)
#define mapcellnew (mapcell*)cellnew
#define mapccellnew (mapccell*)cellnew

#define rmapread(A...)    (rmap_t*)read_by_id(M_RECTMAP, 0, A)    
#define rmapwrite(out, A...)   write_by_id(out?out->base:NULL, M_RECTMAP, A)
#define rmapcellnew  (rmapcell*)cellnew
#define rmapccellnew (rmapccell*)cellnew

#define locwrite(out, A...) write_by_id(out?out->base:NULL, M_LOC, A)
#define locread(A...)    (loc_t*)read_by_id(M_LOC, 0, A)
#define loccellread(A...) (loccell*)read_by_id(M_LOC, 1, A)
#define loccellnew (loccell*)cellnew
#define locccellnew (locccell*)cellnew
/* Read needs type checking, so don't use readbin*/
/**read a dmat*/
#define dread(A...)    dmat_cast(read_by_id(M_REAL, 0, A))
#define dcellnew (dcell*)cellnew
#define dccellnew (dccell*)cellnew
#define dcellreaddata(fp, header) dcell_cast(readdata_by_id(fp, M_REAL, 1, header))
#define dcellread(A...) (dcell*)read_by_id(M_REAL, 1, A)
#define dccellread(A...) (dccell*)read_by_id(M_REAL, 2, A)
#define dcccellread(A...) (dcccell*)read_by_id(M_REAL, 3, A)
/**read a smat*/
#define sread(A...)    smat_cast(read_by_id(M_FLT, 0, A))
#define scellnew (scell*)cellnew
#define sccellnew (sccell*)cellnew
#define scellreaddata(fp, header) scell_cast(readdata_by_id(fp, M_FLT, 1, header))
#define scellread(A...) (scell*)read_by_id(M_FLT, 1, A)
#define sccellread(A...) (sccell*)read_by_id(M_FLT, 2, A)
#define scccellread(A...) (scccell*)read_by_id(M_FLT, 3, A)
/**read a cmat*/
#define cread(A...)    cmat_cast(read_by_id(M_COMP, 0, A))
#define ccellnew (ccell*)cellnew
#define cccellnew (cccell*)cellnew
#define ccellreaddata(fp, header) ccell_cast(readdata_by_id(fp, M_COMP, 1, header))
#define ccellread(A...) (ccell*)read_by_id(M_COMP, 1, A)
#define cccellread(A...) (cccell*)read_by_id(M_COMP, 2, A)
#define ccccellread(A...) (ccccell*)read_by_id(M_COMP, 3, A)
/**read a cmat*/
#define zread(A...)    zmat_cast(read_by_id(M_ZMP, 0, A))
#define zcellnew (zcell*)cellnew
#define zccellnew (zccell*)cellnew
#define zcellreaddata(fp, header) zcell_cast(readdata_by_id(fp, M_ZMP, 1, header))
#define zcellread(A...) (zcell*)read_by_id(M_ZMP, 1, A)
#define zccellread(A...) (zccell*)read_by_id(M_ZMP, 2, A)
#define zcccellread(A...) (zcccell*)read_by_id(M_ZMP, 3, A)
/**read a lmat*/
#define lread(A...) lmat_cast(read_by_id(M_LONG, 0, A))
#define lcellnew (lcell*)cellnew
#define lccellnew (lccell*)cellnew
#define lcellreaddata(fp, header) lcell_cast(readdata_by_id(fp, M_LONG, 1, header))
#define lcellread(A...) (lcell*)read_by_id(M_LONG, 1, A)
#define lccellread(A...) (lccell*)read_by_id(M_LONG, 2, A)
#define lcccellread(A...) (lcccell*)read_by_id(M_LONG, 3, A)

/**read a dsp*/
#define dspread(A...) dsp_cast(read_by_id(M_DSP, 0, A))
#define dspcellread(A...) dspcell_cast(read_by_id(M_DSP, 1, A))
#define sspread(A...) ssp_cast(read_by_id(M_SSP, 0, A))
#define sspcellread(A...) sspcell_cast(read_by_id(M_SSP, 1, A))
/**read a zsp*/
#define cspread(A...) csp_cast(read_by_id(M_CSP, 0, A))
#define cspcellread(A...) cspcell_cast(read_by_id(M_CSP, 1, A))
#define zspread(A...) zsp_cast(read_by_id(M_ZSP, 0, A))
#define zspcellread(A...) zspcell_cast(read_by_id(M_ZSP, 1, A))
#define dspcellnew (dspcell*)cellnew
#define dspccellnew (dspccell*)cellnew
#define cspcellnew (cspcell*)cellnew
#define cspccellnew (cspccell*)cellnew

#define writebin_delayed(A,format...) \
{\
  char fnstore[PATH_MAX];\
  snprintf(fnstore, PATH_MAX, format);\
  const char *fn=fnstore;\
  if((fnstore[0]=='-')||disable_save) fn=NULL;\
  if(A && fn){\
    A->fn=strdup(fn);\
  }\
}
#endif
