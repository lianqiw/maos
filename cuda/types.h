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
#include "kernel.h"
template <typename T>
struct cumat{
    T *p;
    int nx;
    int ny;
    int *nref;
    char *header;
    cumat(int nxi, int nyi, T *pi=NULL, int own=1)
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
    ~cumat(){
	if(nref){
	    nref[0]--;
	    if(nref[0]==0){
		cudaFree(p);
		delete nref;
		if(header) free(header);
	    }else if(nref[0]<0){
		error("Invalid nref=%d\n", nref[0]);
	    }
	}
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	if(this && p){
	    if(stream==(cudaStream_t)-1){
		DO(cudaMemset(p, 0, nx*ny*sizeof(T)));
	    }else{
		DO(cudaMemsetAsync(p, 0, nx*ny*sizeof(T), stream));
	    }
	}
    }
    cumat<T>* ref(){
	if(nref) nref[0]++;
	cumat<T>* res=new cumat<T>(nx, ny, p, 0);
	res->nref=nref;
	return res;
    }
};

template <typename T>
struct cucell{
    cumat<T> **p;
    int nx;
    int ny;
    cumat<T> *m; /*contains the continuous data*/
    void init(const int nxi, const int nyi){
	nx=nxi;
	ny=nyi;
	p=(cumat<T>**)calloc(nx*ny, sizeof(void*));
	m=NULL;
    }
    cucell(const int nx, const int ny){
	init(nx, ny);
    }
    void init(const int nx, const int ny, int mx, int my){
	init(nx,ny);
	m=new cumat<T>(mx*my*nx*ny,1);
	for(int i=0; i<nx*ny; i++){
	    p[i]=mx&&my?new cumat<T>(mx, my, m->p+i*(mx*my), 0):NULL;
	}
    }
    cucell(const int nx, const int ny, int mx, int my){
	init(nx, ny, mx, my);
    }
    template <typename L>
    void init(const int nx, const int ny, L *mx, L *my){
	init(nx,ny);
	long tot=0;
	for(int i=0; i<nx*ny; i++){
	    tot+=mx[i]*(my?my[i]:1);
	}
	m=new cumat<T> (tot,1);
	tot=0;
	for(int i=0; i<nx*ny; i++){
	    p[i]=mx[i]?new cumat<T>(mx[i],(my?my[i]:1),m->p+tot, 0):NULL;
	    tot+=mx[i]*(my?my[i]:1);
	}
    }
    template <typename L>
    cucell(const int nx, const int ny, L *mx, L *my){
	init(nx, ny, mx, my);
    }
    cucell(const cucell<T>*in){
	if(!in->m){
	    init(in->nx, in->ny);
	    for(int i=0; i<in->nx*in->ny; i++){
		p[i]=new cumat<T>(in->p[i]->nx, in->p[i]->ny);
	    }
	}else{
	    int mx[in->nx*in->ny];
	    int my[in->nx*in->ny];
	    for(int i=0; i<in->nx*in->ny; i++){
		mx[i]=in->p[i]->nx;
		my[i]=in->p[i]->ny;
	    }
	    init<int>(in->nx, in->ny, (int*)mx, (int*)my);
	}
    }
    ~cucell(){
	for(int i=0; i<nx*ny; i++){
	    delete p[i];
	}
	delete m;
	free(p);
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	if(!this) return;
	if(m){
	    m->zero();
	}else{
	    for(int i=0; i<nx*ny; i++){
		p[i]->zero(stream);
	    }
	}
    }
};
typedef struct cumat<float>    curmat;
typedef struct cumat<fcomplex> cucmat;
typedef struct cucell<float>  curcell;
typedef struct cucell<fcomplex>  cuccell;

typedef struct cusp{
    int *p;
    int *i;
    float *x;
    int nx;
    int ny;
    int nzmax;
    ~cusp(){
	cudaFree(p);
	cudaFree(i);
	cudaFree(x);
    }
}cusp;
typedef struct cuspcell{
    cusp **p;
    int nx;
    int ny;
    ~cuspcell(){
	for(int i=0; i<nx*ny; i++){
	    delete p[i];
	}
	free(p);
    }
}cuspcell;
typedef struct{
    cuspcell *Mt;
    curcell *U;
    curcell *V;
    stream_t *fitstream;
    stream_t *dmstream;
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
    cumap_t(int nxi, int nyi, float *p=NULL, int own=1, float oxi=0, float oyi=0, float dxi=0, float hti=0, float vxi=0, float vyi=0):
	curmat(nxi, nyi,p,own),ox(oxi),oy(oyi),dx(dxi),ht(hti),vx(vxi),vy(vyi),cubic_cc(NULL){};

}cumap_t;

typedef struct mulock{
    int lock;
    pthread_mutex_t mutex;
    mulock(int dolock=1):lock(dolock){
	if(lock){
	    pthread_mutex_init(&mutex, NULL);
	    LOCK(mutex);
	}
    }
    ~mulock(){
	if(lock){
	    UNLOCK(mutex);
	}
    }
}mulock;

template <typename T>
void initzero(cumat<T> **A, int nx, int ny){
    /*zero array if exist, otherwise allocate and zero*/
    if(*A){
	(*A)->zero();
    }else{
	*A=new cumat<T>(nx,ny);
    }
}

#endif
