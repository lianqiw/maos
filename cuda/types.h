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
#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
#include <typeinfo>
#include <string>
using std::string;
#include "common.h"
#include "kernel.h"
class nonCopyable{
private:
	nonCopyable& operator=(nonCopyable&);
	nonCopyable(const nonCopyable&);
protected:
	nonCopyable(){}
};
/**
   Cpu, Pinned, and GPU template classes are created to allow different memory
   allocation.
*/
//Standard CPU memory. 
template<typename T>
class Cpu{
	T val;
public:
	static void Zero(T* p, size_t size, cudaStream_t stream=0){
		(void)stream;
		if(p) memset(p, 0, size*sizeof(T));
	}
	void* operator new[](size_t size){
		return ::calloc(size, 1);
	}
	void operator delete[](void* p){
		::free(p);
	}
	static void Copy(T* pout, const T* pin, size_t size, cudaStream_t stream=0){
		(void) stream;
		memcpy(pout, pin, size*sizeof(T));
		//dbg("Cpu::Copy %p to %p for %lux%lu bytes\n", pin, pout, size, sizeof(T));
	}
	/*T* cpuPointer(T* p, size_t size){
		return p;
	}*/
};

//Pinned (page locked) CPU memory
template<typename T>
class Pinned:public Cpu<T>{
public:
	void* operator new[](size_t size){
		void* p=0;
		cudaMallocHost(&p, size);
		memset(p, 0, size);
		return p;
	}
		void operator delete[](void* p){
		cudaFreeHost(p);
	}
};
extern int cuda_dedup;
//Gpu memory
template<typename T>
class Gpu{
	T val;
public:
	void* operator new[](size_t size){
		void* p=0;
		DO(cudaMalloc(&p, size));
		DO(cudaMemset(p, 0, size));
		return p;
	}
	void operator delete[](void* p){
		DO(cudaFree(p));
	}
	static void Zero(T* p, size_t size, cudaStream_t stream=0){
		if(p){
			DO(cudaMemsetAsync(p, 0, size*sizeof(T), stream));
		}
	}
	static void Copy(T* pout, const T* pin, size_t size, cudaStream_t stream=0){
		DO(cudaMemcpyAsync(pout, pin, size*sizeof(T), D2D, stream));
	}
	/*T* cpuPointer(T* p, size_t size){
		T* p2=mymalloc(size, T);
		DO(cudaMemcpy(p2, p, size*sizeof(T), cudaMemcpyDeviceToHost));
		return p2;
	}
	void cpuPointerFree(T* p){
		free(p);
	}*/
};

template <typename T, template<typename> class Dev=Cpu >
class Array;
/**
   Functions to Set memory to zero. Array of fundamental types are treated
   differently from array of Arrays.
*/
template <typename T, template<typename> class Dev >
void Zero(Dev<T>* p, long n, cudaStream_t stream=0){
	Dev<T>::Zero((T*)p, n, stream);
	if(sizeof(T)>8){
		warning("call dev zero for %s. This will lead to code error.\n", typeid(T).name());
	}
}

//partially specialize for array of array
template <typename T, template<typename> class Dev >
void Zero(Cpu<Array<T, Dev> >* p_, long n, cudaStream_t stream=0){
	Array<T, Dev>* p=(Array<T, Dev>*)p_;
	for(long i=0; i<n; i++){
		Zero((Dev<T>*)p[i](), p[i].N(), stream);
	}
}
/**
   RefP is a reference counting for pointers.
*/
template <typename T, template<typename> class Dev>
class RefP{
protected:
	T* p;  //The memory address being used.
private:
	Dev<T>* p0; //The allocated memory address. p>=p0;
	int* nref;  //The reference counter. Delete p0 only if nref is valid and has value of 1.
	void init_(long n=0, T* pin=0, int own=1){
	//Initialize from uninitialized state.
		if(n>0){
			p0=(Dev<T>*)pin;
			if(!p0){
				p0=new Dev<T>[n];
				own=1;
			}
			if(own){
				nref=mymalloc(1, int);
				nref[0]=1;
			}
			p=(T*)p0;
		}
	}
public:
	int NRef() const{//return number of reference
		return nref?nref[0]:0;
	}
	void deinit(){
		if(nref&&!atomicadd(nref, -1)){
			delete[] p0;
			myfree(nref);
		}
		nref=0;
		p0=0;
		p=0;
	}

	//Constructors and related
	RefP():p(0), p0(0), nref(0){}
	explicit RefP(long n, T* pin=0, int own=1):p(0), p0(0), nref(0){
		init_(n, pin, own);
	}
	//Create a new pointer with offset for p.
	RefP(const RefP& pin, long offset=0):p(pin.p+offset), p0(pin.p0), nref(pin.nref){
		if(nref) atomicadd(nref, 1);
	}
	/*
	  template<template <typename> class Dev2>
	  RefP &operator=(const RefP<T, Dev2> &in);
	*/

	RefP& operator=(const RefP& in){
		if(this!=&in){
			deinit();//20190927: destroy old data (fix)
			p=in.p;
			p0=in.p0;
			nref=in.nref;
			if(nref) atomicadd(nref, 1);
		}
		return *this;
	}

	void init(long n=0){
		deinit();
		init_(n);
	}
	//Destructors and related
	~RefP(){
		deinit();
	}

	//Access operators
	//() operator
	T* operator()(){
		return (T*)p;
	}
	const T* operator()() const{
		return (T*)p;
	}
	//conversion operator
	operator T* (){
		return (T*)p;
	}
	operator const T* ()const{
		return (T*)p;
	}
	T& operator()(int i){
		return p[i];
	}
	const T& operator()(int i)const{
		return p[i];
	}
	bool operator==(const RefP& in){
		return p==in.p;
	}
	T* operator+(int i){
		return p+i;
	}
	const T* operator+(int i)const{
		return p+i;
	}
};

/**
   Generic array of basic types and classes.
   Need to call new on p to handle cases when T is a class.
*/
template <typename T, template<typename> class Dev>
class Array:private RefP<T, Dev>{
	typedef RefP<T, Dev> Parent;
protected:
	long nx;
	long ny;
public:
	string header;
	using Parent::deinit;
	using Parent::operator+;
	using Parent::operator T*;
	using Parent::operator const T*;
	using Parent::operator();
	using Parent::p;
	using Parent::NRef;
	~Array(){
		nx=0;
		ny=0;
	}
	T& operator ()(int ix, int iy){
		return p[ix+nx*iy];
	}
	const T& operator ()(int ix, int iy)const{
		return p[ix+nx*iy];
	}

	T* Col(int icol){
		return p+nx*icol;
	}
	const T* Col(int icol)const{
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
		return (nx&&ny);
	}

	void init(long nxi=0, long nyi=1){
		nx=nxi;
		ny=nyi;
		Parent::init(nxi*nyi);
	}
	//Constructors 
	Array(long nxi=0, long nyi=1, T* pi=NULL, int own=1)
		:Parent(nxi* nyi, pi, own), nx(nxi), ny(nyi){}

		//Create a reference with offset.
		//Array(long nxi,long nyi,const Parent& pi,long offset=0):Parent(pi,offset),nx(nxi),ny(nyi){
		//}
		//Create a reference with offset.
	Array(long nxi, long nyi, const Array& pi, long offset=0):Parent(pi, offset), nx(nxi), ny(nyi){}
	//Use default destructor

	//Need to handle both basic types and classes. Use template function.
	//Cannot partially specialize single member function.
	void Zero(cudaStream_t stream=0){
		::Zero((Dev<T>*)p, nx*ny, stream);
	}

	Array(const Array& in):Parent(in), nx(in.nx), ny(in.ny), header(in.header){}
	Array& operator=(const Array& in){
		if(this!=&in){
			Parent::operator=(in);
			nx=in.nx;
			ny=in.ny;
			header=in.header;
		}
		return *this;
	}

	Array Vector()const{
		Array tmp=*this;
		tmp.nx=tmp.nx*tmp.ny;
		tmp.ny=1;
		return tmp;
	}
	Array trans(stream_t& stream);

	//Copy the data. Zero data if input data is empty. Reallocate if mismatch.
	void Copy(const Array& in, cudaStream_t stream=0){
		if(this!=&in){
			if(!in){
				this->Zero(stream);
			} else {
				if(this->N()!=in.N()){
					if(NRef()>1){
						dbg("Copying into referenced data that has different size.\n");
					}
					if(this->N()!=0){
						dbg("Reinitiate data during copying, old size is %ldx%ld, new size is %ldx%ld\n", Nx(), Ny(), in.Nx(), in.Ny());
					}
					deinit();
					init(in.Nx(), in.Ny());
				}
				Dev<T>::Copy((*this)(), in(), in.N(), stream);
			}
		}
	}
};

/**
   Cell is a special Array that can stores multiple Arrays of data in continuous
   region.
*/

template <typename T, template<typename> class Dev=Cpu >
class Cell:public Array<Array<T, Dev>, Cpu>{
private:
	typedef Array<T, Dev> TMat;
	typedef Array<Array<T, Dev>, Cpu> Parent;
	TMat m; /*contains the continuous data*/
protected:
	using Parent::nx;
	using Parent::ny;
public:
	Array<T*, Pinned>pm_cpu;/*contains the data pointer in each cell in gpu.*/
	Array<T*, Gpu>pm;/*contains the data pointer in each cell in cpu.*/
	using Parent::operator();
	using Parent::p;
	TMat& M(){
		return m;
	}
	const TMat& M() const{
		return m;
	}

	void p2pm(cudaStream_t stream=0){
		if(nx&&ny){
			pm_cpu.init(nx, ny);
			pm.init(nx, ny);
			for(long i=0; i<nx*ny; i++){
				pm_cpu[i]=p[i]();
			}
			DO(cudaMemcpyAsync(pm(), pm_cpu(), sizeof(T*)*nx*ny, cudaMemcpyHostToDevice, stream));
		}
	}

	Cell(long nxi=0, long nyi=1):Parent(nxi, nyi){
		p2pm();
	}
	Cell(long nxi, long nyi, long mx, long my, T* pin=NULL):Parent(nxi, nyi){
		if(mx&&my){
			m=TMat(mx*my*nxi*nyi, 1, pin, 0);
			for(int i=0; i<nxi*nyi; i++){
				p[i]=TMat(mx, my, m, i*(mx*my));
			}
			p2pm();
		}
	}
	template <typename L>
	Cell(long nxi, long nyi, L* mx, L* my, T* pin=NULL):Parent(nxi, nyi){
		long tot=0;
		for(long i=0; i<nxi*nyi; i++){
			tot+=mx[i]*(my?my[i]:1);
		}
		m=TMat(tot, 1, pin, 0);
		tot=0;
		for(long i=0; i<nxi*nyi; i++){
			if(mx[i]){
				p[i]=TMat(mx[i], (my?my[i]:1), m, tot);
				tot+=p[i].N();
			}
		}
		p2pm();
	}

	Cell(const Cell& in):Parent(in), m(in.m), pm_cpu(in.pm_cpu), pm(in.pm){}

	Cell& operator=(const Cell& in){
		if(this!=&in){
			Parent::operator=(in);
			m=in.m;
			pm_cpu=in.pm_cpu;
			pm=in.pm;
		}
		return *this;
	}
	//Use default copy constructor and copy assignment operator
	void init(long nxi=0, long nyi=1){
		Parent::init(nxi*nyi);
		pm_cpu.init();
		pm.init();
	}
	//Create a similar cell from this cell.
	Cell New()const{
		if(!m){
			return Cell(nx, ny);
		} else{
			long mx[nx*ny];
			long my[nx*ny];
			for(long i=0; i<nx*ny; i++){
				mx[i]=p[i].Nx();
				my[i]=p[i].Ny();
			}
			return Cell(nx, ny, mx, my);
		}
	}

	void replace(T* pnew, cudaStream_t stream){
	/*replace the data with a new set. we don't own the pnew.*/
		if(m){
			m=TMat(m.Nx(), m.Ny(), pnew, 0);//replace m
		}

		for(long i=0; i<nx*ny; i++){
			p[i]=TMat(p[i].Nx(), p[i].Ny(), pnew, 0);//don't own the data pnew
			pnew+=p[i].N();
		}
		p2pm(stream);
	}
};
typedef class Array<int, Gpu>   cuimat;
typedef class Array<Real, Gpu>   curmat;
typedef class Array<Comp, Gpu>   cucmat;
typedef class Cell<Real, Gpu>  curcell;
typedef class Cell<Comp, Gpu>  cuccell;
enum TYPE_SP{
	SP_CSC,//compressed sparse column major. 
	SP_CSR,//compressed sparse row major. 
};
class cusp{
	int* p;
	int* i;
	Real* x;
	int nx;
	int ny;
	int nzmax;
	int* nref;
	enum TYPE_SP type;
#if __CUDACC_VER_MAJOR__ >= 10	
	cusparseSpMatDescr_t desc;
#else
	void* desc;
#endif
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
	int* Pp(){
		return p;
	}
	const int* Pp() const{
		return p;
	}
	int* Pi(){
		return i;
	}
	const int* Pi() const{
		return i;
	}

	Real* Px(){
		return x;
	}
	const Real* Px() const{
		return x;
	}
#if __CUDACC_VER_MAJOR__ >= 10	
	cusparseSpMatDescr_t Desc() const{
		return desc;
	}
#endif
	cusp():p(0), i(0), x(0), nx(0), ny(0), nzmax(0), nref(0), type(SP_CSC), desc(0){}
	cusp(const dsp* in, int tocsr=0, int transp=0);
	cusp(const cusp& in):p(in.p), i(in.i), x(in.x), nx(in.nx), ny(in.ny), nzmax(in.nzmax), nref(in.nref), type(in.type), desc(in.desc){
		if(nref) nref[0]++;
	}
	cusp& operator=(const cusp& in){
		if(this!=&in){
			deinit();
			p=in.p;
			i=in.i;
			x=in.x;
			nx=in.nx;
			ny=in.ny;
			nzmax=in.nzmax;
			nref=in.nref;
			if(nref) nref[0]++;
			type=in.type;
			desc=in.desc;
		}
		return *this;
	}
	~cusp(){
		deinit();
	}
	void deinit(){
		if(nref&&!atomicadd(nref, -1)){
#if __CUDACC_VER_MAJOR__ >= 10	
			cusparseDestroySpMat(desc);
#endif
			cudaFree(p);
			cudaFree(i);
			cudaFree(x);
			myfree(nref);
		}
		p=0; i=0; x=0; nref=0;
	}

	void trans();/*covnert to CSR mode by transpose*/

	operator bool()const{
		return nx&&ny;
	}
};
typedef Array<cusp> cuspcell;
typedef Array<curcell> curccell;
typedef Array<curccell> curcccell;
void cp2gpu(curmat& dest, const loc_t* src);
class culoc_t{
private:
	curmat p;/*in device. */
	Real dx;
	Real dy;
public:
	real xmax, xmin, ymax, ymin;
	Real2* operator()(){
		return (Real2*)(p());
	}
	const Real2* operator()()const{
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
	culoc_t(const loc_t* in=0):dx(0), dy(0){
		if(in){
			dx=in->dx;
			dy=in->dy;
			cp2gpu(p, in);
			dmaxmin(in->locx, in->nloc, &xmax, &xmin);
			dmaxmin(in->locy, in->nloc, &ymax, &ymin);
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
	cupts_t(pts_t* in=0):culoc_t((loc_t*)in), dxsa(0), nxsa(0){
		if(in){
			dxsa=in->dx;
			nxsa=in->nx;
		}
	}
};
/**
   Specifies the grid.
*/
class cugrid_t{
public:
	long nx, ny;
	Real ox, oy;
	Real dx, dy;
	Real ht;
	Real vx, vy;
	curmat cubic_cc; /*coefficients for cubic influence function. */
	//use default copy assignment operator and copy constructor

	curmat iac2cc(Real iac){
		if(iac){
			Real cc[5];
			Real cubicn=1.f/(1.f+2.f*iac);
			cc[0]=1.f*cubicn;
			cc[1]=(4.f*iac-2.5f)*cubicn;
			cc[2]=(1.5f-3.f*iac)*cubicn;
			cc[3]=(2.f*iac-0.5f)*cubicn;
			cc[4]=(0.5f-iac)*cubicn;
			curmat res(5, 1);
			DO(cudaMemcpy(res, cc, 5*sizeof(Real), cudaMemcpyHostToDevice));
			return res;
		} else{
			return curmat();
		}
	}
	cugrid_t& operator=(const map_t* in){
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
			cubic_cc=iac2cc(in->iac);
		}
		return *this;
	}
	cugrid_t():nx(0), ny(0), ox(0), oy(0), dx(0), dy(0), ht(0), vx(0), vy(0){}
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
	cumap_t(cugrid_t& grid):cugrid_t(grid){
		if(nx>0&&ny>0){
			p=curmat(nx, ny);
		}
	}
	cumap_t& operator=(map_t* in){
		if(in){
			cugrid_t::operator=(in);
			if(p.Nx()!=nx||p.Ny()!=ny){
				p=curmat(nx, ny);
			}
		}
		return *this;
	}
	operator const curmat& ()const{
		return p;
	}
	operator curmat& (){
		return p;
	}
	Real* operator()(){
		return p();
	}
	const Real* operator()()const{
		return p();
	}
};

typedef Array<cumap_t> cumapcell;
typedef Array<cugrid_t> cugridcell;
template <typename T, template<typename> class Dev >
void initzero(Array<T, Dev>& A, long nx, long ny){
	/*zero array if exist, otherwise allocate and zero*/
	if(A){
		A.Zero();
	} else{
		A.init(nx, ny);
	}
}

#define cuzero(A,B...) (A).Zero(B)
#define cucp(A,B...) (A).Copy(B)
#endif

