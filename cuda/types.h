/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
#include "common.h"
#include <cuComplex.h>
#define fcomplex cuFloatComplex
#define dcomplex cuDoubleComplex
template <typename T>
struct cumat{
    T *p;
    int nx;
    int ny;
    int *nref;
    char *header;
    cumat(int nxi, int nyi, T *pi=NULL, int own=1, cudaStream_t stream=0)
	: p(pi), nx(nxi), ny(nyi), nref(NULL), header(NULL){
	if(!p && nxi!=0 && nyi!=0){
	    DO(cudaMalloc(&p, nxi*nyi*sizeof(T)));
	    DO(cudaMemset(p, 0, nxi*nyi*sizeof(T)));
	}
	if(own){
	    nref=new int[1];
	    nref[0]=1;
	}
    }
    cumat(int nxi, int nyi, cudaStream_t stream){
	cumat<T>(nxi, nyi, NULL, 1, stream);
    }
    ~cumat(){
	if(nref){
	    nref[0]--;
	    if(nref[0]==0){
		cudaFree(p);
		free(nref);
		if(header) free(header);
	    }else if(nref[0]<0){
		error("Invalid nref=%d\n", nref[0]);
	    }
	}
    }
    cumat<T>* ref(){
	cumat<T>*res=new cumat<T>(nx, ny, p, 0);
	res->nref=nref;
	nref[0]++;
	return res;
    }
};

template <typename T>
struct cucell{
    cumat<T> **p;
    int nx;
    int ny;
    cumat<T> *m; /*contains the continuous data*/
};
typedef struct cumat<float>    curmat;
typedef struct cumat<fcomplex> cucmat;
typedef struct cucell<float>  curcell;
typedef struct cucell<fcomplex>  cuccell;

typedef struct{
    int *p;
    int *i;
    float *x;
    int nx;
    int ny;
    int nzmax;
}cusp;
typedef struct{
    cusp **p;
    int nx;
    int ny;
}cuspcell;
typedef struct{
    cuspcell *Mt;
    curcell *U;
    curcell *V;
}cumuv_t;

typedef struct{
    float (*loc)[2];/*in device. */
    float dx;
    int nloc;
}culoc_t;

typedef struct cumap_t:curmat{
    float ox, oy;
    float dx;
    float ht;
    float vx, vy;
    float *cubic_cc; /*coefficients for cubic influence function. */
    cumap_t(int nxi, int nyi, float oxi=0, float oyi=0, float dxi=0, float hti=0, float vxi=0, float vyi=0):
	curmat(nxi, nyi),ox(oxi),oy(oyi),dx(dxi),ht(hti),vx(vxi),vy(vyi),cubic_cc(NULL){};

}cumap_t;

#endif
