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

template <typename T>
class cumat{
private:
    T *p;
    long nx;
    long ny;
    int *nref;
public:
    char *header;
 
    T *P(){
	return (T*)p;
    }
    const T *P() const{
	return (const T*)p;
    }
    operator T*(){
	return P();
    }
    operator const T*()const{
	return P();
    }
    T *Col(int icol){
	return P()+nx*icol;
    }
    const T *Col(int icol)const{
	return P()+nx*icol;
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
    T&operator ()(int ix, int iy){
	return P()[ix+nx*iy];
    }
    const T&operator ()(int ix, int iy)const{
	return P()[ix+nx*iy];
    }
    T *operator+(int off){
	return P()+off;
    }
    const T*operator+(int off)const{
	return P()+off;
    }
    //Constructors 
    cumat(long nxi=0, long nyi=0, T *pi=NULL, int own=1)
	:nx(nxi),ny(nyi),p(pi),nref(NULL),header(NULL){
	if(nxi >0 && nyi >0){
	    if(!p){
		DO(cudaMalloc(&p, nx*ny*sizeof(T)));
		DO(cudaMemset(p, 0, nx*ny*sizeof(T)));
		//info("%p allocated with size %ld\n", p, nx*ny*sizeof(T));
		own=1;
	    }
	    if(own){
		nref=new int;
		nref[0]=1;
	    }
	}
    }
    void reset(){
	if(nref && !atomicadd(nref, -1)){
	    //warning("%p free with size %ld\n", p, nx*ny*sizeof(T));
	    DO(cudaFree(p));
	    delete nref;
	    if(header) free(header);
	}
	p=NULL;
    }
    ~cumat(){
	reset();
    }
    cumat(const cumat &in):p(in.p),nx(in.nx),ny(in.ny),nref(in.nref),header(in.header){
	if(nref) nref[0]++;
    }
    cumat &operator=(const cumat &in){
	if(&p!=&in.p){//prevent self assignment
	    reset();
	    p=in.p;
	    nx=in.nx;
	    ny=in.ny;
	    nref=in.nref;
	    header=in.header;
	    if(nref) nref[0]++;
	}
	return *this;
    }
    bool operator==(const cumat&in){
	return P()==in.P();
    }
    void zero(cudaStream_t stream=(cudaStream_t)-1){
	if(p){
	    if(stream==(cudaStream_t)-1){
		DO(cudaMemset(p, 0, nx*ny*sizeof(T)));
	    }else{
		DO(cudaMemsetAsync(p, 0, nx*ny*sizeof(T), stream));
	    }
	}
    }
    cumat Vector(){
	cumat tmp=*this;
	tmp.nx=tmp.nx*tmp.ny;
	tmp.ny=1;
	return tmp;
    }
    cumat<T>trans(stream_t &stream);
};

template <typename T>
class cucell{
private:
    cumat<T> *p;
    long nx;
    long ny;
    cumat<T> m; /*contains the continuous data*/
    int *nref;
public:
    T **pm; /*contains the data pointer in each cell in gpu.*/
    T **pm_cpu;/*contains the data pointer in each cell in cpu.*/

    cumat<T> *P(){
	return p;
    }
    const cumat<T> *P() const{
	return p;
    }
    cumat<T> &M(){
	return m;
    }
    const cumat<T>&M() const{
	return m;
    }
    cumat<T>& operator()(int i){
	return p[i];
    }
    const cumat<T>& operator()(int i)const{
	return p[i];
    }
    cumat<T>& operator[](int i){
	return p[i];
    }
    const cumat<T>& operator[](int i)const{
	return p[i];
    }
    cumat<T>& operator()(int ix, int iy){
	return p[ix+iy*nx];
    }
    const cumat<T>& operator()(int ix, int iy)const{
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
    cumat<T> *operator+(int off){
	return p+off;
    }
    const cumat<T>*operator+(int off)const{
	return p+off;
    }
    bool operator==(const cucell&in){
	return P()==in.P();
    }
    void p2pm(cudaStream_t stream=(cudaStream_t)-1){
	if(!nref || nref[0]>1){
	    error("cannot act on referenced array\n");
	}
	if(!p) error("p must not be null\n");
	if(!pm_cpu){
	    pm_cpu=(T**)malloc4async(nx*ny*sizeof(T*));
	}
	for(long i=0; i<nx*ny; i++){
	    pm_cpu[i]=p[i].P();
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
    void init(long nxi, long nyi){
	nx=nxi;
	ny=nyi;
	if(nx&&ny) {
	    p=new cumat<T>[nx*ny];
	    nref=new int;
	    nref[0]=1;
	}else{
	    p=0;
	    nref=0;
	}
	pm=NULL;
	pm_cpu=NULL;
    }
    cucell(long nxi=0, long nyi=1){
	init(nxi, nyi);
    }
    cucell(long nxi, long nyi, long mx, long my, T *pin=NULL){
	init(nxi, nyi);
	if(mx && my){
	    m=cumat<T>(mx*my*nxi*nyi,1,pin,pin?0:1);
	    for(int i=0; i<nxi*nyi; i++){
		p[i]=cumat<T>(mx, my, m.P()+i*(mx*my), 0);
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
	m=cumat<T>(tot,1,pin,pin?0:1);
	tot=0;
	for(long i=0; i<_nx*_ny; i++){
	    if(mx[i]){
		p[i]=cumat<T>(mx[i],(my?my[i]:1),m.P()+tot, 0);
		tot+=mx[i]*(my?my[i]:1);
	    }
	}
	p2pm();
    }
    cucell(const cucell<T>&in):p(in.p),nx(in.nx),ny(in.ny),m(in.m),nref(in.nref),pm(in.pm),pm_cpu(in.pm_cpu){
	if(nref){
	    nref[0]++;
	}
    }
    cucell& operator=(const cucell<T>&in){
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
	    free4async(pm_cpu);
	    cudaFree(pm);
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
	    m=cumat<T>(m.Nx(), m.Ny(), pnew, 0);//replace m
	}

	for(long i=0; i<nx*ny; i++){
	    p[i]=cumat<T>(p[i].Nx(), p[i].Ny(), pnew, 0);//don't own the data pnew
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
typedef class cumat<Real>   curmat;
typedef class cumat<Comp>   cucmat;
typedef class cucell<Real>  curcell;
typedef class cucell<Comp>  cuccell;
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
template <typename T>
class cuarray{
    T *p;
    int nx;
    int ny;
    int *nref;
 public:
    cuarray(int _nx=0, int _ny=1):p(0),nx(_nx),ny(_nx?_ny:0),nref(0){
	if(nx && ny){
	    p=new T[nx*ny];
	    nref=new int;
	    nref[0]=1;
	}
    }
    cuarray(const cuarray&in):p(in.p),nx(in.nx),ny(in.ny),nref(in.nref){
	if(nref) nref[0]++;
    }
    cuarray &operator=(const cuarray&in){
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
    ~cuarray(){
	if(nref && !atomicadd(nref, -1)){
	    delete []p;
	    delete nref;
	}
    }
    operator bool()const {
	return nx && ny;
    }
    T *P(){
	return p;
    }
    const T *P() const{
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
    operator T*(){
	return P();
    }
    operator const T*()const{
	return P();
    }
};
typedef cuarray<cusp> cuspcell;
typedef cuarray<curcell> curccell;
typedef cuarray<curccell> curcccell;
void cp2gpu(curmat &dest, const loc_t *src);
class culoc_t{
private:
    curmat p;/*in device. */
    Real dx;
    Real dy;
public:
    operator Real2*(){
	return P();
    }
    operator const Real2*()const{
	return P();
    }
    Real2* P(){
	return (Real2*)(p.P());
    }
    const Real2* P()const{
	return (const Real2*)(p.P());
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
    /*
    cugrid_t(const cugrid_t &in):nx(in.nx), ny(in.ny), ox(in.ox), oy(in.oy), dx(in.dx), dy(in.dy),
				 ht(in.ht), vx(in.vx), vy(in.vy),
				 iac(in.iac),cubic_cc(gpu_dmcubic_cc(in.iac)){
    }
    cugrid_t &operator=(const cugrid_t &in){
	if(&nx!=&in.nx){
	    nx=in.nx;
	    ny=in.ny;
	    ox=in.ox;
	    oy=in.oy;
	    dx=in.dx;
	    dy=in.dy;
	    ht=in.ht;
	    vx=in.vx;
	    vy=in.vy;
	    iac=in.iac;
	    cubic_cc=in.cubic_cc;
	}
	return *this;
	}*/
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
    /*cugrid_t(long nxi=0, long nyi=0,Real oxi=0, Real oyi=0, Real dxi=0, Real dyi=0, Real hti=0, Real vxi=0, Real vyi=0):nx(nxi),ny(nyi),ox(oxi),oy(oyi),dx(dxi),dy(dyi),ht(hti),vx(vxi),vy(vyi),cubic_cc(_cubic_cc){
    }
    cugrid_t scale(Real sc)const{
	return cugrid_t(nx,ny,ox*sc,oy*sc,dx*sc,dy*sc,ht,vx,vy,cubic_cc);
    }
    cugrid_t operator *(Real sc)const{
	return cugrid_t(nx,ny,ox*sc,oy*sc,dx*sc,dy*sc,ht,vx,vy,cubic_cc);
	}*/
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
    /*cumap_t(map_t *map=0):cugrid_t(map){
	if(nx>0 && ny>0) {
	    p=curmat(nx,ny);
	}
	}*/
    operator const curmat&()const{
	return p;
    }
    operator curmat&(){
	return p;
    }
    Real *P(){
	return p.P();
    }
    const Real *P()const{
	return p.P();
    }
};

typedef cuarray<cumap_t> cumapcell;
typedef cuarray<cugrid_t> cugridcell;
template <typename T>
void initzero(cumat<T> &A, long nx, long ny){
    /*zero array if exist, otherwise allocate and zero*/
    if(A){
	A.zero();
    }else{
	A=cumat<T>(nx,ny);
    }
}

#define cuzero(A,B...) (A).zero(B)
#endif
