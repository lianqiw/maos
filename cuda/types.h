/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
class nonCopyable{
private:
    nonCopyable& operator=(nonCopyable&);
    nonCopyable(const nonCopyable&);
protected:
    nonCopyable(){}
};
//Standard CPU memory
class Cpu{
public:
    static void *calloc(size_t size){
	return ::calloc(size, 1);
    }
    static void free(void *p){
	::free(p);
    }
    static void zero(void *p, size_t size, cudaStream_t stream){
	(void) stream;
	if(p) memset(p, 0, size);
    }
};
//Pinned (page locked) CPU memory
class Pinned{
public:
    static void *calloc(size_t size){
	void *p=0;
	cudaMallocHost(&p, size);
	memset(p, 0, size);
	return p;
    }
    static void free(void *p){
	cudaFreeHost(p);
    }
    static void zero(void *p, size_t size, cudaStream_t stream){
	(void) stream;
	if(p) memset(p, 0, size);
    }
};

//Gpu memory
class Gpu{
public:
    static void *calloc(size_t size){
	void *p=0;
	DO(cudaMalloc(&p, size));
	DO(cudaMemset(p, 0, size));
	return p;
    }
    static void free(void *p){
	DO(cudaFree(p));
    }
    static void zero(void *p, size_t size, cudaStream_t stream){
	if(p){
	    if(stream==(cudaStream_t)-1){
		DO(cudaMemset(p, 0, size));
	    }else{
		DO(cudaMemsetAsync(p, 0, size, stream));
	    }
	}
    }
};
template <typename T, class Dev=Cpu >
class Array{
private:
    T *p;
    long nx;
    long ny;
    int *nref;
public:
    char *header;
    operator T*(){
	return p;
    }
    operator const T*()const{
	return p;
    }
    T* operator()(){
	return p;
    }
    const T* operator()() const{
	return p;
    }
    T *Col(int icol){
	return p+nx*icol;
    }
    const T *Col(int icol)const{
	return p+nx*icol;
    }
 
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long N()const{
	return nx*ny;
    }
    operator bool()const{
	return (nx && ny);
    }
    
    T& operator()(int i){
	return p[i];
    }
    const T& operator()(int i)const{
	return p[i];
    }
    T&operator ()(int ix, int iy){
	return p[ix+nx*iy];
    }
    const T&operator ()(int ix, int iy)const{
	return p[ix+nx*iy];
    }
    T *operator+(int off){
	return p+off;
    }
    const T*operator+(int off)const{
	return p+off;
    }
    void init(long nxi, long nyi){
	deinit();
	nx=nxi;
	ny=nyi;
	p=(T*)Dev::calloc(nx*ny*sizeof(T));
	nref=new int;
	nref[0]=1;
    }
    //Constructors 
    Array(long nxi=0, long nyi=1, T *pi=NULL, int own=1)
	:nx(nxi),ny(nyi),p(pi),nref(NULL),header(NULL){
	if(nxi >0 && nyi >0){
	    if(!p){
		p=(T*)Dev::calloc(nx*ny*sizeof(T));
		own=1;
	    }
	    if(own){
		nref=new int;
		nref[0]=1;
	    }
	}
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	Dev::zero(p, nx*ny*sizeof(T), stream);
    }
    void deinit(){
	if(nref && !atomicadd(nref, -1)){
	    Dev::free(p);
	    delete nref;
	    if(header) free(header);
	}
	p=NULL;
    }
    ~Array(){
	deinit();
    }
    Array(const Array &in):p(in.p),nx(in.nx),ny(in.ny),nref(in.nref),header(in.header){
	if(nref) nref[0]++;
    }
    Array &operator=(const Array &in){
	if(&p!=&in.p){//prevent self assignment
	    deinit();
	    p=in.p;
	    nx=in.nx;
	    ny=in.ny;
	    nref=in.nref;
	    header=in.header;
	    if(nref) nref[0]++;
	}
	return *this;
    }
    bool operator==(const Array&in){
	return p==in.p;
    }
    Array Vector(){
	Array tmp=*this;
	tmp.nx=tmp.nx*tmp.ny;
	tmp.ny=1;
	return tmp;
    }
    Array trans(stream_t &stream);
};
template <typename T, class Dev>
class cucell{
private:
    typedef Array<T,Dev> Tmat;
    Tmat *p;
    long nx;
    long ny;
    Tmat m; /*contains the continuous data*/
    int *nref;
public:
    //T **pm; /*contains the data pointer in each cell in gpu.*/
    //T **pm_cpu;/*contains the data pointer in each cell in cpu.*/
    Array<T*,Pinned>pm_cpu;
    Array<T*,Gpu>pm;
    Tmat &M(){
	return m;
    }
    const Tmat&M() const{
	return m;
    }
    Tmat* operator()(){
	return p;
    }
    const Tmat* operator()() const{
	return p;
    }
    
    Tmat& operator()(int i){
	return p[i];
    }
    const Tmat& operator()(int i)const{
	return p[i];
    }
    Tmat& operator[](int i){
	return p[i];
    }
    const Tmat& operator[](int i)const{
	return p[i];
    }
    Tmat& operator()(int ix, int iy){
	return p[ix+iy*nx];
    }
    const Tmat& operator()(int ix, int iy)const{
	return p[ix+iy*nx];
    }
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long N()const{
	return nx*ny;
    }
    Tmat *operator+(int off){
	return p+off;
    }
    const Tmat*operator+(int off)const{
	return p+off;
    }
    bool operator==(const cucell&in){
	return p==in.p;
    }
    void p2pm(cudaStream_t stream=(cudaStream_t)-1){
	if(!nref || nref[0]>1){
	    error("cannot act on referenced array\n");
	}
	if(!p) error("p must not be null\n");
	if(!pm_cpu){
	    //pm_cpu=(T**)malloc4async(nx*ny*sizeof(T*));
	    pm_cpu.init(nx, ny);
	}
	for(long i=0; i<nx*ny; i++){
	    pm_cpu[i]=p[i]();
	}
	if(!pm){
	    pm.init(nx,ny);
	}
	if(stream==(cudaStream_t)-1){
	    cudaMemcpy(pm(), pm_cpu(), sizeof(T*)*nx*ny,cudaMemcpyHostToDevice);
	}else{
	    cudaMemcpyAsync(pm(), pm_cpu(), sizeof(T*)*nx*ny,cudaMemcpyHostToDevice, stream);
	}
    }
    void init(long nxi, long nyi){
	nx=nxi;
	ny=nyi;
	if(nx&&ny) {
	    p=new Tmat[nx*ny];
	    nref=new int;
	    nref[0]=1;
	}else{
	    p=0;
	    nref=0;
	}
	//pm=NULL;
	//pm_cpu=NULL;
    }
    cucell(long nxi=0, long nyi=1){
	init(nxi, nyi);
    }
    cucell(long nxi, long nyi, long mx, long my, T *pin=NULL){
	init(nxi, nyi);
	if(mx && my){
	    m=Tmat(mx*my*nxi*nyi,1,pin,pin?0:1);
	    for(int i=0; i<nxi*nyi; i++){
		p[i]=Tmat(mx, my, m()+i*(mx*my), 0);
	    }
	    p2pm();
	}
    }
    template <typename L>
    cucell(const long _nx, const long _ny, L *mx, L *my, T *pin=NULL){
	init(_nx,_ny);
	long tot=0;
	for(long i=0; i<_nx*_ny; i++){
	    tot+=mx[i]*(my?my[i]:1);
	}
	m=Tmat(tot,1,pin,pin?0:1);
	tot=0;
	for(long i=0; i<_nx*_ny; i++){
	    if(mx[i]){
		p[i]=Tmat(mx[i],(my?my[i]:1),m()+tot, 0);
		tot+=mx[i]*(my?my[i]:1);
	    }
	}
	p2pm();
    }
    cucell(const cucell&in):p(in.p),nx(in.nx),ny(in.ny),m(in.m),nref(in.nref),pm(in.pm),pm_cpu(in.pm_cpu){
	if(nref){
	    nref[0]++;
	}
    }
    cucell& operator=(const cucell&in){
	if(&p!=&in.p){
	    reset();
	    p=in.p;
	    nx=in.nx;
	    ny=in.ny;
	    m=in.m;
	    nref=in.nref;
	    if(nref) nref[0]++;
	    pm=in.pm;
	    pm_cpu=in.pm_cpu;
	}
	return *this;
    }
    ~cucell(){
	reset();
    }
    void reset(){
	if(nref && !atomicadd(nref, -1)){
	    delete [] p;
	    pm_cpu.deinit();
	    pm.deinit();
	    //free4async(pm_cpu);
	    //cudaFree(pm);
	    delete nref;
	}
	p=0;
	nref=0;
	nx=0;
	ny=0;
    }
    //Create a similar cell from this cell.
    cucell New()const{
	if(!m){
	    return cucell(nx, ny);
	}else{
	    long mx[nx*ny];
	    long my[nx*ny];
	    for(long i=0; i<nx*ny; i++){
		mx[i]=p[i].Nx();
		my[i]=p[i].Ny();
	    }
	    return cucell(nx, ny, mx, my);
	}
    }
 
    void replace(T *pnew, cudaStream_t stream){
	/*replace the data with a new set. free original data if free_original
	  is set. we don't own the pnew.*/
	if(m){
	    m=Tmat(m.Nx(), m.Ny(), pnew, 0);//replace m
	}

	for(long i=0; i<nx*ny; i++){
	    p[i]=Tmat(p[i].Nx(), p[i].Ny(), pnew, 0);//don't own the data pnew
	    pnew+=p[i].N();
	}
	p2pm(stream);
    }
 
    operator bool()const{
	return (nx && ny);
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	if(m){
	    m.zero(stream);
	}else{
	    for(long i=0; i<nx*ny; i++){
		p[i].zero(stream);
	    }
	}
    }
};
typedef class Array<int, Gpu>   cuimat;
typedef class Array<Real, Gpu>   curmat;
typedef class Array<Comp, Gpu>   cucmat;
typedef class cucell<Real, Gpu>  curcell;
typedef class cucell<Comp, Gpu>  cuccell;
enum TYPE_SP{
    SP_CSC,//compressed sparse column major. 
    SP_CSR,//compressed sparse row major. 
};
class cusp{
    int *p;
    int *i;
    Real *x;
    int nx;
    int ny;
    int nzmax;
    int *nref;
    enum TYPE_SP type;
 public:
    enum TYPE_SP Type()const{
	return type;
    }
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long Nzmax()const{
	return nzmax;
    }
    int *Pp(){
	return p;
    }
    const int *Pp() const{
	return p;
    }
    int *Pi(){
	return i;
    }
    const int *Pi() const{
	return i;
    }

    Real *Px(){
	return x;
    }
    const Real *Px() const{
	return x;
    }
    cusp():p(0),i(0),x(0),nx(0),ny(0),nzmax(0),nref(0),type(SP_CSC){}
    cusp(const dsp *in, int tocsr);
    cusp(const cusp&in):p(in.p),i(in.i),x(in.x),nx(in.nx),ny(in.ny),nzmax(in.nzmax),nref(in.nref),type(in.type){
	if(nref) nref[0]++;
    }
    cusp &operator=(const cusp&in){
	if(&p!=&in.p){
	    reset();
	    p=in.p;
	    i=in.i;
	    x=in.x;
	    nx=in.nx;
	    ny=in.ny;
	    nzmax=in.nzmax;
	    nref=in.nref;
	    if(nref) nref[0]++;
	    type=in.type;
	}
	return *this;
    }
    ~cusp(){
	reset();
    }
    void reset(){
	if(nref && !atomicadd(nref, -1)){
	    cudaFree(p);
	    cudaFree(i);
	    cudaFree(x);
	    delete nref;
	}
	p=0; i=0; x=0; nref=0;
    }
 
    void trans();/*covnert to CSR mode by transpose*/
    /*cusp* ref(void){
	if(nref) atomicadd(nref, 1);
	cusp* res=(cusp*)malloc(sizeof(cusp));
	memcpy(res, this, sizeof(*this));
	return res;
	}*/
    operator bool()const{
	  return nx && ny;
    }
};
/*
template <typename T>
class Array{
    T *p;
    int nx;
    int ny;
    int *nref;
 public:
    Array(int _nx=0, int _ny=1):p(0),nx(_nx),ny(_nx?_ny:0),nref(0){
	if(nx && ny){
	    p=new T[nx*ny];
	    nref=new int;
	    nref[0]=1;
	}
    }
    Array(const Array&in):p(in.p),nx(in.nx),ny(in.ny),nref(in.nref){
	if(nref) nref[0]++;
    }
    Array &operator=(const Array&in){
	if(&p!=&in.p){
	    p=in.p;
	    nx=in.nx;
	    ny=in.ny;
	    nref=in.nref;
	    if(nref){
		nref[0]++;
	    }
	}
	return *this;
    }
    ~Array(){
	if(nref && !atomicadd(nref, -1)){
	    delete []p;
	    delete nref;
	}
    }
    operator bool()const {
	return nx && ny;
    }
    T * operator()(){
	return p;
    }
    const T *operator()() const{
	return p;
    }
    T& operator()(int i){
	return p[i];
    }
    const T& operator()(int i)const{
	return p[i];
    }
    T& operator[](int i){
	return p[i];
    }
    const T& operator[](int i)const{
	return p[i];
    }
    T& operator()(int ix, int iy){
	return p[ix+iy*nx];
    }
    const T& operator()(int ix, int iy)const{
	return p[ix+iy*nx];
    }
    long Nx()const{
	return nx;
    }
    long Ny()const{
	return ny;
    }
    long N()const{
	return nx*ny;
    }
    T *operator+(int off){
	return p+off;
    }
    const T*operator+(int off)const{
	return p+off;
    }
    };*/
//template <typename T> using Array=Array<T,Cpu>;
typedef Array<cusp> cuspcell;
typedef Array<curcell> curccell;
typedef Array<curccell> curcccell;
void cp2gpu(curmat &dest, const loc_t *src);
class culoc_t{
private:
    curmat p;/*in device. */
    Real dx;
    Real dy;
public:
    Real2* operator()(){
	return (Real2*)(p());
    }
    const Real2*operator()()const{
	return (Real2*)(p());
    }
    long Nloc()const{
	return p.Ny();
    }
    Real Dx()const{
	return dx;
    }
    Real Dy()const{
	return dy;
    }
    //No need custom copy assignment operator or copy constructor.
    culoc_t(const loc_t *in=0):dx(0),dy(0){
	if(in){
	    dx=in->dx;
	    dy=in->dy;
	    cp2gpu(p, in);
	}
    }
};
class cupts_t:public culoc_t{
    Real dxsa;
    long nxsa;
public:
    Real Dxsa(){
	return dxsa;
    }
    long Nxsa(){
	return nxsa;
    }
    cupts_t(pts_t *in=0):culoc_t((loc_t*)in),dxsa(0),nxsa(0){
	if(!in) return;
	dxsa=in->dx;
	nxsa=in->nx;
    }
};
curmat gpu_dmcubic_cc(Real iac);
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
    Real iac;
    curmat cubic_cc; /*coefficients for cubic influence function. */
    //use default copy assignment operator and copy constructor

    cugrid_t &operator=(const map_t *in){
	if(in){
	    nx=in->nx;
	    ny=in->ny;
	    ox=in->ox;
	    oy=in->oy;
	    dx=in->dx;
	    dy=in->dy;
	    ht=in->h;
	    vx=in->vx;
	    vy=in->vy;
	    if(fabs(iac-in->iac)>fabs(iac+in->iac)*1e-5){
		iac=in->iac;
		cubic_cc=gpu_dmcubic_cc(in->iac);
	    }
	}
	return *this;
    }
    cugrid_t():nx(0),ny(0),ox(0),oy(0),dx(0),dy(0),ht(0),vx(0),vy(0),iac(0){}
    cugrid_t Scale(Real sc)const{
	cugrid_t tmp(*this);
	tmp.ox*=sc;
	tmp.oy*=sc;
	tmp.dx*=sc;
	tmp.dy*=sc;
	return tmp;
    }
    operator bool(){
	return nx&&ny;
    }
};
class cumap_t:public cugrid_t{
public:
    curmat p;
    /*Init the data, p*/
    cumap_t(){}
    cumap_t(cugrid_t &grid):cugrid_t(grid){
	if(nx>0 && ny>0) {
	    p=curmat(nx,ny);
	}
    }
    cumap_t& operator=(map_t*in){
	if(in){
	    cugrid_t::operator=(in);
	    if(p.Nx()!=nx || p.Ny()!=ny){
		p=curmat(nx, ny);
	    }
	}
	return *this;
    }
    operator const curmat&()const{
	return p;
    }
    operator curmat&(){
	return p;
    }
    Real *operator()(){
	return p();
    }
    const Real *operator()()const{
	return p();
    }
};

typedef Array<cumap_t> cumapcell;
typedef Array<cugrid_t> cugridcell;
template <typename T, class Dev>
void initzero(Array<T, Dev> &A, long nx, long ny){
    /*zero array if exist, otherwise allocate and zero*/
    if(A){
	A.zero();
    }else{
	A=Array<T, Dev>(nx,ny);
    }
}

#define cuzero(A,B...) (A).zero(B)
#endif
