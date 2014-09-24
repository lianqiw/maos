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
#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
#include "common.h"
#include "kernel.h"
template <typename T>
class cumat{
 public:
    T *p;
    long nx;
    long ny;
    int *nref;
    char *header;
    
    /*Constructors */
    cumat(long nxi, long nyi, T *pi=NULL, int own=1)
	:nx(nxi),ny(nyi),p(pi),nref(NULL),header(NULL){
	if(!p && nxi >0 && nyi >0){
	    DO(cudaMalloc(&p, nx*ny*sizeof(T)));
	    DO(cudaMemset(p, 0, nx*ny*sizeof(T)));
	    own=1;
	}
	if(own){
	    nref=new int;
	    nref[0]=1;
	}
    }
    virtual ~cumat(){
	if(nref && !atomicadd(nref, -1)){
	    DO(cudaFree(p));
	    delete nref;
	    if(header) free(header);
	}
	p=NULL;
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
    cumat* ref(int vector=0){
	if(!this) return NULL;
	if(nref) atomicadd(nref, 1);
	cumat* res;
	if(vector){
	    res=new cumat<T>(nx*ny, 1, p, 0);
	}else{
	    res=new cumat<T>(nx, ny, p, 0);
	}
	res->nref=nref;
	return res;
    }
    cumat<T>*trans(stream_t &stream);
    operator T*(){
	return p;
    }
    T&operator ()(int ix, int iy){
	return p[ix+nx*iy];
    }
};

template <typename T>
class cucell{
 public:
    cumat<T> **p;
    long nx;
    long ny;
    cumat<T> *m; /*contains the continuous data*/
    T **pm; /*contains the data pointer in each cell in gpu.*/
    T **pm_cpu;/*contains the data pointer in each cell in cpu.*/
    void p2pm(cudaStream_t stream=(cudaStream_t)-1){
	if(!p) error("p must not be null\n");
	if(!pm_cpu){
	    pm_cpu=(T**)malloc4async(nx*ny*sizeof(T*));
	}
	for(long i=0; i<nx*ny; i++){
	    pm_cpu[i]=p[i]?p[i]->p:NULL;
	}
	if(!pm){
	    DO(cudaMalloc(&pm, sizeof(T*)*nx*ny));
	}
	if(stream==(cudaStream_t)-1){
	    cudaMemcpy(pm, pm_cpu, sizeof(T*)*nx*ny,cudaMemcpyHostToDevice);
	}else{
	    cudaMemcpyAsync(pm, pm_cpu, sizeof(T*)*nx*ny,cudaMemcpyHostToDevice, stream);
	}
    }
    void init(const long nxi, const long nyi){
	nx=nxi;
	ny=nyi;
	p=(cumat<T>**)calloc(nx*ny, sizeof(void*));
	m=NULL;
	pm=NULL;
	pm_cpu=NULL;
    }
    cucell(const long nxi, const long nyi){
	init(nxi, nyi);
    }
    void init(const long nxi, const long nyi, long mx, long my, T *pin=NULL){
	init(nxi, nyi);
	m=new cumat<T>(mx*my*nxi*nyi,1,pin,pin?0:1);
	for(int i=0; i<nxi*nyi; i++){
	    p[i]=(mx&&my)?new cumat<T>(mx, my, m->p+i*(mx*my), 0):NULL;
	}
	p2pm();
    }
    cucell(const long _nx, const long _ny, long mx, long my, T *pin=NULL){
	init(_nx, _ny, mx, my, pin);
    }
    template <typename L>
    void init(const long _nx, const long _ny, L *mx, L *my, T *pin=NULL){
	init(_nx,_ny);
	long tot=0;
	for(long i=0; i<_nx*_ny; i++){
	    tot+=mx[i]*(my?my[i]:1);
	}
	m=new cumat<T> (tot,1,pin,pin?0:1);
	tot=0;
	for(long i=0; i<_nx*_ny; i++){
	    p[i]=mx[i]?new cumat<T>(mx[i],(my?my[i]:1),m->p+tot, 0):NULL;
	    tot+=mx[i]*(my?my[i]:1);
	}
	p2pm();
    }
    template <typename L>
    cucell(const long _nx, const long _ny, L *mx, L *my, T *pin=NULL){
	init(_nx, _ny, mx, my, pin);
    }
    cucell(const cucell<T>*in){
	if(!in->m){
	    init(in->nx, in->ny);
	    for(long i=0; i<in->nx*in->ny; i++){
		p[i]=new cumat<T>(in->p[i]->nx, in->p[i]->ny);
	    }
	}else{
	    long mx[in->nx*in->ny];
	    long my[in->nx*in->ny];
	    for(long i=0; i<in->nx*in->ny; i++){
		mx[i]=in->p[i]->nx;
		my[i]=in->p[i]->ny;
	    }
	    init<long>(in->nx, in->ny, mx, my);
	}
    }
    void replace(T *pnew, int free_original, cudaStream_t stream){
	/*replace the data with a new set. free original data if free_original
	  is set. we don't own the pnew.*/
	if(free_original){
	    if(m){
		cumat<T> *m2=new cumat<T>(m->nx, m->ny, pnew, 0);
		delete m;
		m=m2;
	    }else{
		for(long i=0; i<nx*ny; i++){
		    if(p[i]){
			cumat<T> *p2=new cumat<T>(p[i]->nx, p[i]->ny, (T*)-1, 0);
			delete p[i];
			p[i]=p2;
		    }
		}
	    }
	}
	if(m){
	    if(m->nref){
		free(m->nref);
		m->nref=0;
	    }
	    m->p=pnew;
	}
	for(long i=0; i<nx*ny; i++){
	    if(p[i]){
		p[i]->p=pnew;
		pnew+=p[i]->nx*p[i]->ny;
	    }
	}
	p2pm(stream);
    }
    ~cucell(){
	for(long i=0; i<nx*ny; i++){
	    delete p[i];
	}
	delete m;
	free(p);
	cudaFree(pm);
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	if(!this) return;
	if(m){
	    m->zero();
	}else{
	    for(long i=0; i<nx*ny; i++){
		p[i]->zero(stream);
	    }
	}
    }
};
typedef class cumat<Real>    curmat;
typedef class cumat<Comp> cucmat;
typedef class cucell<Real>  curcell;
typedef class cucell<Comp>  cuccell;
enum TYPE_SP{
    SP_CSC,
    SP_CSR,
};
class cusp{
 public:
    int *p;
    int *i;
    Real *x;
    int nx;
    int ny;
    int nzmax;
    int *nref;
    enum TYPE_SP type;
    //cusp():p(NULL),i(NULL),x(NULL),nx(0),ny(0),nzmax(0),type(SP_CSC){}
    cusp(const dsp *in=0, int tocsr=1);
    ~cusp(){
	if(nref && !atomicadd(nref, -1)){
	    cudaFree(p);
	    cudaFree(i);
	    cudaFree(x);
	    delete nref;
	}
    }
    void toCSR(){
	if(type==SP_CSC){
	    trans();
	    type=SP_CSR;
	}
    }
    void toCSC(){
	if(type==SP_CSR){
	    trans();
	    type=SP_CSC;
	}
    }
    void trans();/*covnert to CSR mode by transpose*/
    cusp* ref(void){
	if(!this) return NULL;
	if(nref) atomicadd(nref, 1);
	cusp* res=(cusp*)malloc(sizeof(cusp));
	memcpy(res, this, sizeof(*this));
	return res;
    }
};
class cuspcell{
 public:
    cusp **p;
    int nx;
    int ny;
    cuspcell(int _nx, int _ny):nx(_nx),ny(_ny){
	p=(cusp**)calloc(nx*ny, sizeof(void*));
    }
    ~cuspcell(){
	for(long i=0; i<nx*ny; i++){
	    delete p[i];
	}
	free(p);
    }
};

void cp2gpu(Real (* restrict *dest)[2], const loc_t *src);
class culoc_t{
public:
    Real (*p)[2];/*in device. */
    Real dx;
    Real dy;
    long nloc;
    culoc_t(loc_t *in=0):p(0),dx(0),dy(0),nloc(0){
	if(!in) return;
	dx=in->dx;
	dy=in->dy;
	nloc=in->nloc;
	cp2gpu(&p, in);
    }
    ~culoc_t(){
	if(p) cudaFree(p);
    }
};
class cupts_t:public culoc_t{
public:
    Real dxsa;
    int nxsa;
    cupts_t(pts_t *in=0):culoc_t((loc_t*)in),dxsa(0),nxsa(0){
	if(!in) return;
	dxsa=in->dx;
	nxsa=in->nx;
    }
};
curmat* gpu_dmcubic_cc(Real iac);
template <typename T>
void gpu2gpu(cumat<T> **dest, cumat<T> *source){
    if(!*dest){
	*dest=new cumat<T>(source->nx, source->ny);
    }else{
	assert((*dest)->nx*(*dest)->ny==source->nx*source->ny);
    }
    cudaMemcpy((*dest)->p, source->p, sizeof(T)*source->nx*source->ny, cudaMemcpyDeviceToDevice);
}
/**
   Specifies the grid.
*/
class cugrid_t{
public:
    long  nx, ny;
    Real ox, oy;
    Real dx, dy;
    Real ht;
    Real vx, vy;
    curmat *cubic_cc; /*coefficients for cubic influence function. */
    void init(const map_t *in){
	nx=in->nx;
	ny=in->ny;
	ox=in->ox;
	oy=in->oy;
	dx=in->dx;
	dy=in->dy;
	ht=in->h;
	vx=in->vx;
	vy=in->vy;
	if(cubic_cc){
	    delete cubic_cc;
	}
	if(in->cubic){
	    cubic_cc=gpu_dmcubic_cc(in->iac);
	}else{
	    cubic_cc=0;
	}
    }
    void zero(){
	nx=ny=0;
	ox=oy=dx=dy=ht=vx=vy=0;
	cubic_cc=0;
    }
    void init(const cugrid_t *in){
	if(in){
	    nx=in->nx;
	    ny=in->ny;
	    ox=in->ox;
	    oy=in->oy;
	    dx=in->dx;
	    dy=in->dy;
	    ht=in->ht;
	    vx=in->vx;
	    vy=in->vy;
	    if(in->cubic_cc){
		gpu2gpu(&cubic_cc, in->cubic_cc);
	    }
	}else{
	    zero();
	}
    }
    cugrid_t(cugrid_t *in){
	init(in);
    }
    cugrid_t(const map_t *in){
	if(in){
	    init(in);
	}else{
	    zero();
	}
    }
    cugrid_t(long nxi=0, long nyi=0,Real oxi=0, Real oyi=0, Real dxi=0, Real dyi=0, Real hti=0, Real vxi=0, Real vyi=0,curmat *_cubic_cc=0):nx(nxi),ny(nyi),ox(oxi),oy(oyi),dx(dxi),dy(dyi),ht(hti),vx(vxi),vy(vyi),cubic_cc(_cubic_cc){}
    ~cugrid_t(){
	delete cubic_cc;
    }
    cugrid_t scale(Real sc)const{
	return cugrid_t(nx,ny,ox*sc,oy*sc,dx*sc,dy*sc,ht,vx,vy,cubic_cc);
    }
    cugrid_t operator *(Real sc)const{
	return cugrid_t(nx,ny,ox*sc,oy*sc,dx*sc,dy*sc,ht,vx,vy,cubic_cc);
    }
};
class cumap_t:public cugrid_t{
public:
    curmat *p;
    /*Init the data, p*/
    cumap_t(cugrid_t *grid=0):cugrid_t(grid),p(0){
	if(nx>0 && ny>0) {
	    p=new curmat(nx,ny);
	}
    }
    cumap_t(map_t *map):cugrid_t(map){
	if(nx>0 && ny>0) {
	    p=new curmat(nx,ny);
	}
    }
    operator curmat*()const{
	return p;
    }
};


template <typename T>
void initzero(cumat<T> **A, long nx, long ny){
    /*zero array if exist, otherwise allocate and zero*/
    if(*A){
	(*A)->zero();
    }else{
	*A=new cumat<T>(nx,ny);
    }
}
#endif
