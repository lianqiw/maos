#ifndef AOS_LIB_ZMAT_H
#define AOS_LIB_ZMAT_H
/**
   \file dmat.h Contains the mathematically functions regarding to dmat and dcell object
*/
#include "mat.h"
#include "cell.h"
#include "matbin.h"
AOS_MAT_DEF   (AOS_ZMAT,AOS_SMAT,AOS_SSP,fcomplex,float)
AOS_CELL_DEF  (AOS_ZMAT,AOS_ZSP,fcomplex)
AOS_MATBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

//The following are only for dmat.
//dmat *denc(dmat *A, dmat *dvec, int type, int nthread);
#endif
