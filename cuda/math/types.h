/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	void *operator new[](size_t size){
		//return ::calloc(size, 1);
		return (void*) mycalloc(size, char);
	}
	void operator delete[](void* p){
		//::free(p);
		myfree(p);
	}
	static void Zero(T *p, size_t size, cudaStream_t stream=0){
		(void)stream;
		if(p) memset(p, 0, size*sizeof(T));
	}
	static void Copy(T *pout, const T *pin, size_t size, cudaStream_t stream=0){
		(void)stream;
		memcpy(pout, pin, size*sizeof(T));
		//dbg("Cpu::Copy %p to %p for %lux%lu bytes\n", pin, pout, size, sizeof(T));
	}
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
	using Cpu<T>::Zero;
	using Cpu<T>::Copy;
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
};

template <typename T, template<typename> class Dev=Cpu >
class Array;

/**
   RefP is a reference counting for pointers.
*/
template <typename T, template<typename> class Dev>
class RefP{
protected:
	T* p;  //The memory address being used.
private:
	Dev<T>* p0; //The allocated memory address. p>=p0;
	unsigned int* nref;  //The reference counter. Delete p0 only if nref is valid and has value of 1.
	void init_(long n=0, T* pin=0, int own=1){
	//Initialize from uninitialized state.
		if(n>0){
			p0=(Dev<T>*)pin;
			if(!p0){
				p0=new Dev<T>[n];
				own=1;
			}
			if(own){
				nref=mymalloc(1, unsigned int);
				nref[0]=1;
			}
			p=(T*)p0;
		}
	}
public:
	typedef T type;
	int NRef() const{//return number of reference
		return nref?nref[0]:0;
	}
	void deinit(){
		if(nref&&!atomic_sub_fetch(nref, 1)){
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
		if(nref) atomic_add_fetch(nref, 1);
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
			if(nref) atomic_add_fetch(nref, 1);
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
   
   Does not defines operations that are different between Numerical and Cell array.
*/
template <typename T, template<typename> class Dev>
class Array:private RefP<T, Dev>{
	typedef RefP<T, Dev> Parent;
protected:
	long nx;
	long ny;
public:
	string keywords;
	using Parent::deinit;
	using Parent::operator+;
	using Parent::operator T*;
	using Parent::operator const T*;
	using Parent::operator();
	using Parent::p;
	using Parent::NRef;
	typedef T type;
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
	T& R(int ix, int iy=0){//relaxed indexing
		return p[(nx==1?0:ix)+nx*(ny==1?0:iy)];
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
	Array():Parent(0,NULL,1),nx(0),ny(0){}
	Array(long nxi, long nyi=1, T* pi=NULL, int own=1)
		:Parent(nxi*nyi, pi, own), nx(nxi), ny(nyi){}

	//Create a reference with offset.
	Array(long nxi, long nyi, const Array& pi, long offset=0):Parent(pi, offset), nx(nxi), ny(nyi){}
	//Use default destructor

	Array(const Array& in):Parent(in), nx(in.nx), ny(in.ny), keywords(in.keywords){}
	Array& operator=(const Array& in){
		if(this!=&in){
			Parent::operator=(in);
			nx=in.nx;
			ny=in.ny;
			keywords=in.keywords;
		}
		return *this;
	}

	Array Vector()const{
		Array tmp=*this;
		tmp.nx=tmp.nx*tmp.ny;
		tmp.ny=1;
		return tmp;
	}
};
///Used for array in cpu of any other array.
template <typename T>
class CellArray:public Array<T, Cpu>{
public:
	typedef Array<T, Cpu> Parent;
	using Parent::deinit;
	using Parent::operator+;
	using Parent::operator T *;
	using Parent::operator const T *;
	using Parent::operator();
	using Parent::p;
	using Parent::NRef;
	using Parent::Col;
	using Parent::Nx;
	using Parent::Ny;
	using Parent::N;
	using Parent::operator bool;
	using Parent::init;
	CellArray Vector()const{
		CellArray tmp(this->N(), 1);
		for(long i=0; i<this->N(); i++){
			tmp[i]=(*this)[i].Vector();
		}
		return tmp;
	}
	CellArray trans(stream_t &stream);
	//using Parent::Copy;
	CellArray(long nxi, long nyi=1, T *pi=NULL, int own=1):Parent(nxi, nyi, pi, own){};
	CellArray(long nxi, long nyi, const CellArray &pi, long offset=0):Parent(nxi, nyi, pi, offset){};
	CellArray():Parent(){};

	void Zero(cudaStream_t stream=0){
		for(long i=0; i<N(); i++){
			(*this)(i).Zero(stream);
		}
	}

	//Copy the data. Zero data if input data is empty. Reallocate if mismatch.
	void Copy(const CellArray &in, cudaStream_t stream=0){
		if(this!=&in){
			if(!in){
				this->Zero(stream);
				warning("Copying from empty input, use Zero().\n");
			} else{
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
				for(long i=0; i<N(); i++){
					(*this)(i).Copy(in(i), stream);
				}
			}
		}
	}
};

///Only for array of numerical data.
template <typename T, template<typename> class Dev>
class NumArray:public Array<T, Dev>{
#if __cpp_static_assert >= 200410L
	static_assert(sizeof(T)<=16, "NumArray can only be used with numerical data");
#endif	
public:
	typedef Array<T, Dev> Parent;
	using Parent::deinit;
	using Parent::operator+;
	using Parent::operator T *;
	using Parent::operator const T *;
	using Parent::operator();
	using Parent::p;
	using Parent::NRef;
	using Parent::Col;
	using Parent::Nx;
	using Parent::Ny;
	using Parent::N;
	using Parent::operator bool;
	using Parent::init;
#if CUDA_VERSION>10000
	cudaDataType dtype();
	mutable cusparseDnVecDescr_t vdesc=0;//for cuspmul
	mutable cusparseDnMatDescr_t mdesc=0;//for cuspmul
#endif
	NumArray trans(stream_t &stream);
	//using Parent::Copy;
	NumArray():Parent(){};
	NumArray(long nxi, long nyi=1, T *pi=NULL, int own=1):Parent(nxi, nyi, pi, own){};
	NumArray(long nxi, long nyi, const NumArray& pi, long offset=0):Parent(nxi, nyi, pi, offset){};
	NumArray(const dmat *A);//C Array wrapper. Only specilize for matching type
	NumArray(const cmat *A);
	NumArray(const smat *A);
	NumArray(const zmat *A);
#if CUDA_VERSION>10000
	NumArray(const NumArray &in):Parent(in){}//don't copy vdesc, mdesc
	NumArray &operator=(const NumArray &in){
		if(this!=&in){
			Parent::operator=(in);
			vdesc=NULL;//don't copy vdesc, mdesc
			mdesc=NULL;
		}
		return *this;
	}
	~NumArray(){
		if(vdesc) {DO(cusparseDestroyDnVec(vdesc)); vdesc=NULL;}
		if(mdesc) {DO(cusparseDestroyDnMat(mdesc)); mdesc=NULL;}
	}
#endif
	NumArray Vector()const{
		NumArray tmp=*this;
		tmp.nx=tmp.nx*tmp.ny;
		tmp.ny=1;
		return tmp;
	}
	//Implementations that are different from regular Array.
	void Zero(cudaStream_t stream=0){
		Dev<T>::Zero(p, N(), stream);
	}

	//Copy the data. Zero data if input data is empty. Reallocate if mismatch.
	void Copy(const NumArray &in, cudaStream_t stream=0){
		if(this!=&in){
			if(!in){
				this->Zero(stream);
				warning("Copying from empty input, ignored.\n");
			} else{
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
   Cell is a special Array of array of numerical data that can store multiple Arrays of numerical data in continuous region.
*/

template <typename T, template<typename> class Dev=Cpu >
class Cell:public CellArray<NumArray<T, Dev> >{
private:
	typedef NumArray<T, Dev> TMat;
	typedef CellArray<NumArray<T, Dev> > Parent;
	TMat m; /*contains the continuous data*/
	Array<T *, Pinned>pm_cpu;/*contains the data pointer in each cell in gpu.*/
protected:
	using Parent::nx;
	using Parent::ny;
public:
	Array<T*, Gpu>pm;/*contains the data pointer in each cell in cpu.*/
	using Parent::operator();
	using Parent::p;
	using Parent::Copy;
	using Parent::Zero;
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
	Cell Vector()const{
		Cell tmp(this->N(), 1);
		for(long i=0; i<this->N(); i++){
			tmp[i]=(*this)[i].Vector();
		}
		tmp.pm=this->pm;
		tmp.pm_cpu=this->pm_cpu;
		tmp.m=this->m;
		return tmp;
	}
};
typedef class NumArray<int, Gpu>   cuimat;
typedef class NumArray<Real, Gpu>   curmat;
typedef class NumArray<Comp, Gpu>   cucmat;
typedef class NumArray<real, Gpu>   cudmat;//same dtype as cpu dmat
typedef class NumArray<comp, Gpu>   cuzmat;//same dtype as cpu cmat
typedef class Cell<Real, Gpu>  curcell;
typedef class Cell<Comp, Gpu>  cuccell;
typedef class NumArray<real, Cpu>   crmat;//equivalent to rmat in c
typedef class NumArray<comp, Cpu>   ccmat;//equivalent to cmat in c

template <typename T>
class dtype{
	int operator()();
};

enum TYPE_SP{
	SP_CSC,//compressed sparse column major. 
	SP_CSR,//compressed sparse row major. 
};
class cusp;
/**
 * @brief The content of cusp is moved into a separate struct so that reference can be done correctly.
 * 
 */
struct cusp_ref{
	friend class cusp;
	friend void cuspmul(curmat &y, const cusp &A, const curmat &x, long ncolvec, char trans, Real alpha, stream_t &stream);
	Spint *p=NULL;
	Spint *i=NULL;
	Real *x=NULL;
	long nx=0;
	long ny=0;
	long nzmax=0;
	mutable void *bspmv=NULL; //buffer for spmv cuda-version 12.4 and 12.5 crash if a different buffer is used for spmv with the same matrix.
	mutable void *bspmm=NULL; //buffer for spmm
	mutable long nbspmm=0; //number of columns for the bspmm
	enum TYPE_SP type=SP_CSC;
	unsigned int count=1;
#if CUDA_VERSION>=10000
	cusparseSpMatDescr_t desc=NULL;
#else
	void *desc=NULL;
#endif
	~cusp_ref(){
		if(count){
			warning("Deleted while count is %d which should be 0\n", count);
		}
#if CUDA_VERSION >= 10000
		if(desc) cusparseDestroySpMat(desc);
#endif			
		if(bspmv) cudaFree(bspmv);
		if(bspmm) cudaFree(bspmm);

		cudaFree(p);
		cudaFree(i);
		cudaFree(x);
	}
};
class cusp{
	cusp_ref *ref=NULL;
public:
	friend void cuspmul(curmat &y, const cusp &A, const curmat &x, long ncolvec, char trans, Real alpha, stream_t &stream);
	//The following operations assume the caller have checked that the cusp is not empty.
	enum TYPE_SP Type()const{
		return ref->type;
	}
	long Nx()const{
		return ref->nx;
	}
	long Ny()const{
		return ref->ny;
	}
	long Nzmax()const{
		return ref->nzmax;
	}
	Spint* Pp(){
		return ref->p;
	}
	const Spint* Pp() const{
		return ref->p;
	}
	Spint* Pi(){
		return ref->i;
	}
	const Spint* Pi() const{
		return ref->i;
	}

	Real* Px(){
		return ref->x;
	}
	const Real* Px() const{
		return ref->x;
	}
	cusp(){}
	cusp(const dsp* in, int tocsr=0, int transp=0);
	cusp(const cusp& in):ref(in.ref){
		if(ref) ref->count++;
	}
	cusp& operator=(const cusp& in){
		if(this!=&in && this->ref!=in.ref){
			deinit();
			ref=in.ref;
			if(ref) ref->count++;
		}
		return *this;
	}
	~cusp(){
		deinit();
	}
	void deinit(){
		if(ref && !atomic_sub_fetch(&ref->count, 1)){
			delete ref;
		}
	}

	void trans();/*covnert to CSR mode by transpose*/

	operator bool()const{
		return ref && ref->nx&&ref->ny;
	}
};
typedef CellArray<cusp> cuspcell;
typedef CellArray<curcell> curccell;
typedef CellArray<curccell> curcccell;
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
			nxsa=in->nxsa;
		}
	}
};
curmat iac2cc(Real iac);

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
	Real dratio;/**<demagnification ratio (>1 for shrinking beam)*/
	curmat cubic_cc; /*coefficients for cubic influence function. */
	//use default copy assignment operator and copy constructor

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
			dratio=in->dratio;
			cubic_cc=iac2cc(in->iac);
		}
		return *this;
	}
	cugrid_t():nx(0), ny(0), ox(0), oy(0), dx(0), dy(0), ht(0), vx(0), vy(0),dratio(1){}
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
	void Zero(cudaStream_t stream=0){
		p.Zero(stream);
	}
	void Copy(const cumap_t& in, cudaStream_t stream=0){
		cugrid_t::operator=(in);
		p.Copy(in.p, stream);
	}
};

typedef CellArray<cumap_t> cumapcell;
typedef CellArray<cugrid_t> cugridcell;

template <typename T>
void Copy(CellArray<T> &out, const CellArray<T>& in, cudaStream_t stream=0){
	out.Copy(in, stream);
}
//Copy between same device
template <typename T, template<typename> class Dev>
void Copy(NumArray<T, Dev> &out, const NumArray<T, Dev> &in, cudaStream_t stream=0){
	out.Copy(in, stream);
}
//Copy between different device, same or different type
template <typename T, template<typename> class Dev, typename T2, template<typename> class Dev2>
void Copy(NumArray<T, Dev> &out, const NumArray<T2, Dev2> &in, cudaStream_t stream=0){
	cp2gpu(out(), in(), in.Nx(), in.Ny(), stream);
}
static inline void Copy(cumap_t&out, const cumap_t &in, cudaStream_t stream=0){
	out.cugrid_t::operator=(in);
	out.Copy(in, stream);
}

template <typename T, template<typename> class Dev >
void initzero(NumArray<T, Dev>& A, long nx, long ny){
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

