#ifndef AOS_LIB_SMAT_H
#define AOS_LIB_SMAT_H
/**
   \file dmat.h Contains the mathematically functions regarding to dmat and dcell object
*/
#include "mat.h"
#include "cell.h"
#include "matbin.h"
AOS_MAT_DEF   (AOS_SMAT,AOS_SMAT,AOS_SSP,float,float)
AOS_CELL_DEF  (AOS_SMAT,AOS_SSP,float)
AOS_MATBIN_DEF(AOS_SMAT,AOS_SSP,float)

//The following are only for dmat.
//dmat *denc(dmat *A, dmat *dvec, int type, int nthread);
#endif
