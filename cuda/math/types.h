/*
  Copyright 2009-2026 Lianqi Wang
  
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
	/**
	   Protected default constructor; copying is disabled.
	*/
	nonCopyable(){}
};
template <typename T>
struct Magic{
	enum{ magic=0};
};
template <>
struct Magic<double>{
	enum{ magic=M_DBL };
};
template <>
struct Magic<float>{
	enum{ magic=M_FLT };
};
template <>
struct Magic<double2>{
	enum{ magic=M_CMP };
};
template <>
struct Magic<float2>{
	enum{ magic=M_ZMP };
};
template <>
struct Magic<int>{//always 32 bit int
	enum{ magic=M_INT32 };
};
template <>
struct Magic<long int>{//always 64 bit int. long is not.
	enum{ magic=M_INT64 };
};

#define ASSERT_NUMERIC\
	static_assert(Magic<T>::magic!=0)

/**
   Cpu, Pinned, and GPU template classes are created to allow different memory
   allocation.
*/
//Standard CPU memory. 
template<typename T>
class Cpu{
	T val;
public:
	/**
	   Allocate zero initialized host memory.

	   @param[in] size Number of bytes.
	   @return Allocated memory. Freed with delete[].
	*/
	void *operator new[](size_t size){
		return (void *)mycalloc(size, char);
	}
	/**
	   Free host memory allocated by new[].

	   @param[in] p Memory to free.
	*/
	void operator delete[](void* p){
		myfree(p);
	}
	/**
	   Zero host memory after synchronizing stream if it is not 0.

	   @param[out] p Memory to zero. Nothing is done if NULL.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to synchronize first.
	*/
	static void Zero(T *p, size_t size, cudaStream_t stream=0){
		if(stream) CUDA_SYNC_STREAM;
		if(p) memset(p, 0, size*sizeof(T));
	}
	/**
	   Copy host memory after synchronizing stream if it is not 0.

	   @param[out] pout Destination.
	   @param[in] pin Source.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to synchronize first.
	*/
	static void Copy(T *pout, const T *pin, size_t size, cudaStream_t stream=0){
		if(stream) CUDA_SYNC_STREAM;
		memcpy(pout, pin, size*sizeof(T));
		//dbg("Cpu::Copy %p to %p for %lux%lu bytes\n", pin, pout, size, sizeof(T));
	}
	/**
	   Copy host memory with type conversion after synchronizing stream if it is
	   not 0.

	   @param[out] pout Destination.
	   @param[in] pin Source.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to synchronize first.
	*/
	template<typename S>
	static void Copy(T *pout, const S *pin, size_t size, cudaStream_t stream=0){
		if(stream) CUDA_SYNC_STREAM;
		for(size_t i=0; i<size; i++){
			pout[i]=pin[i];
		}
	}
	/**
	   Write host array to file.

	   @param[in] p Array of nx*ny.
	   @param[in] nx Number of rows.
	   @param[in] ny Number of columns.
	   @param[in] keywords Header keywords. May be NULL.
	   @param[in,out] fp File to write to.
	   @param[in] stream Unused.
	*/
	static void Write(T *p, long nx, long ny, const char *keywords, file_t *fp, cudaStream_t stream=0){
		writearr(fp, 0, sizeof(T), Magic<T>::magic, keywords, p, nx, ny);
	}
};
template<typename T>
class CpuObj{
	T val;
public:
	/**
	   Allocate zero initialized host memory for objects.

	   @param[in] size Number of bytes.
	   @return Allocated memory. Freed with delete[].
	*/
	void *operator new[](size_t size){
		return (void *)mycalloc(size, char);
	}
	/**
	   Free host memory allocated by new[].

	   @param[in] p Memory to free.
	*/
	void operator delete[](void *p){
		myfree(p);
	}
	/**
	   Copy objects by assignment after synchronizing stream if it is not 0.

	   @param[out] pout Destination.
	   @param[in] pin Source.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to synchronize first.
	*/
	static void Copy(T *pout, const T *pin, size_t size, cudaStream_t stream=0){
		if(stream) CUDA_SYNC_STREAM;
		for(size_t i=0; i<size; i++){
			pout[i]=pin[i];
		}
		//dbg("Cpu::Copy %p to %p for %lux%lu bytes\n", pin, pout, size, sizeof(T));
	}
};
//Pinned (page locked) CPU memory
template<typename T>
class Pinned:public Cpu<T>{
public:
	/**
	   Allocate zero initialized page locked host memory.

	   @param[in] size Number of bytes.
	   @return Allocated memory. Freed with delete[].
	*/
	void* operator new[](size_t size){
		void* p=0;
		cudaMallocHost(&p, size);
		memset(p, 0, size);
		return p;
	}
	/**
	   Free page locked host memory.

	   @param[in] p Memory to free.
	*/
	void operator delete[](void* p){
		cudaFreeHost(p);
	}
	using Cpu<T>::Zero;
	using Cpu<T>::Copy;
	using Cpu<T>::Write;
};
extern int cuda_dedup;
//Gpu memory
template<typename T>
class Gpu{
	T val;
public:
	/**
	   Allocate zero initialized device memory.

	   @param[in] size Number of bytes.
	   @return Allocated memory. Freed with delete[].
	*/
	void* operator new[](size_t size){
		void* p=0;
		DO(cudaMalloc(&p, size));
		DO(cudaMemset(p, 0, size));
		return p;
	}
	/**
	   Free device memory.

	   @param[in] p Memory to free.
	*/
	void operator delete[](void* p){
		DO(cudaFree(p));
	}
	/**
	   Zero device memory asynchronously.

	   @param[out] p Memory to zero. Nothing is done if NULL.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to run in.
	*/
	static void Zero(T* p, size_t size, cudaStream_t stream=0){
		if(p){
			DO(cudaMemsetAsync(p, 0, size*sizeof(T), stream));
		}
	}
	/**
	   Copy device memory asynchronously.

	   @param[out] pout Destination.
	   @param[in] pin Source.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to run in.
	*/
	static void Copy(T* pout, const T* pin, long size, cudaStream_t stream=0){
		ASSERT_NUMERIC;
		DO(cudaMemcpyAsync(pout, pin, size*sizeof(T), D2D, stream));
	}
	/**
	   Copy device memory with type conversion. See copy_do().

	   @param[out] pout Destination.
	   @param[in] pin Source.
	   @param[in] size Number of elements.
	   @param[in] stream Stream to run in.
	*/
	template<typename S>
	static void Copy(T *pout, const S *pin, long size, cudaStream_t stream=0){
		ASSERT_NUMERIC;
		copy_do<<<DIM(size, 256), 0, stream>>>(pout, pin, size);
	}
	/**
	   Copy device array to host and write it to file. Synchronizes stream.

	   @param[in] p Device array of nx*ny.
	   @param[in] nx Number of rows.
	   @param[in] ny Number of columns.
	   @param[in] keywords Header keywords. May be NULL.
	   @param[in,out] fp File to write to.
	   @param[in] stream Stream to copy in.
	*/
	static void Write(T *p, long nx, long ny, const char *keywords, file_t *fp, cudaStream_t stream=0){
		ASSERT_NUMERIC;
		T *tmp=mymalloc(nx*ny, T);
		DO(cudaMemcpyAsync(tmp, p, nx*ny*sizeof(T), D2H, stream));
		CUDA_SYNC_STREAM;
		writearr(fp, 0, sizeof(T), Magic<T>::magic, keywords, tmp, nx, ny);
		myfree(tmp);
	}
	/**
	   Copy device array to host and write it to a file named by format.

	   @param[in] p Device array of nx*ny.
	   @param[in] nx Number of rows.
	   @param[in] ny Number of columns.
	   @param[in] keywords Header keywords. May be NULL.
	   @param[in] format printf style format of the file name, followed by arguments.
	*/
	static void Write(T *p, long nx, long ny, const char *keywords, const char *format, ...){
		ASSERT_NUMERIC;
		format2fn;
		file_t *fp=zfopen(fn, "wb");
		Write(p, nx, ny, keywords, fp);
		zfclose(fp);
	}
};

//template <typename T, template<typename> class Dev=Cpu >
//class Array;

/**
   RefP is a reference counting for pointers.
*/
template <typename T, template<typename> class Dev=Cpu>
class RefP{
protected:
	T* p=0;  //The memory address being used.
private:
	Dev<T>* p0=0; //The allocated memory address. p>=p0;
	long n0=0; //the allocated memory length.
	unsigned int* nref=0;  //The reference counter. Delete p0 only if nref is valid and has value of 1.
public:
	/**
	   Initialize the pointer. Old data is released unless it already has the
	   same size and pin is NULL.

	   @param[in] n Number of elements. Nothing is allocated if 0.
	   @param[in] pin Existing memory to use. Allocated (zeroed) if NULL.
	   @param[in] own Whether to free pin when the last reference is released.
	   Always owned if allocated here.
	*/
	void init(long n=0, T* pin=0, int own=1){
		if(p){
			if(n0==n&&p==(T *)p0&&!pin){//same data size already exists. no action.
				//dbg("Same size data already exists, no action.\n");
				return;
			}else{
				deinit();//destroy old data
			}
		}
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
			n0=n;
		}
	}

	typedef T type;
	/**
	   @return Number of references, or 0 if the memory is not owned.
	*/
	int NRef() const{
		return nref?nref[0]:0;
	}
	/**
	   Release the reference and free the memory if this is the last owner.
	*/
	void deinit(){
		if(nref&&!atomic_sub_fetch(nref, 1)){
			delete[] p0;
			myfree(nref);
		}
		nref=0;
		p0=0;
		p=0;
		n0=0;
	}

	//Constructors and related
	/**
	   Create an empty pointer.
	*/
	RefP(){}
	/**
	   Create and initialize. See init().

	   @param[in] n Number of elements.
	   @param[in] pin Existing memory to use. Allocated if NULL.
	   @param[in] own Whether to free pin when the last reference is released.
	*/
	explicit RefP(long n, T* pin=0, int own=1){
		init(n, pin, own);
	}
	/**
	   Create a new reference to existing data with offset for p.

	   @param[in] pin Pointer to reference.
	   @param[in] offset Offset in number of elements.
	*/
	RefP(const RefP& pin, long offset=0):p(pin.p+offset), p0(pin.p0), nref(pin.nref){
		if(nref) atomic_add_fetch(nref, 1);
	}
	/*
	  template<template <typename> class Dev2>
	  RefP &operator=(const RefP<T, Dev2> &in);
	*/

	/**
	   Release the current data and reference the data of in.

	   @param[in] in Pointer to reference.
	   @return Reference to this.
	*/
	RefP& operator=(const RefP& in){
		if(this!=&in){
			deinit();//20190927: destroy old data (fix)
			p=in.p;
			p0=in.p0;
			nref=in.nref;
			n0=in.n0;
			if(nref) atomic_add_fetch(nref, 1);
		}
		return *this;
	}

	//Destructors and related
	/**
	   Release the reference. See deinit().
	*/
	virtual ~RefP(){
		deinit();
	}

	//Access operators
	//() operator
	/**
	   @return The data pointer.
	*/
	T* operator()(){
		return (T*)p;
	}
	/**
	   @return The data pointer.
	*/
	const T* operator()() const{
		return (T*)p;
	}
	//conversion operator
	/**
	   @return The data pointer.
	*/
	operator T* (){
		return (T*)p;
	}
	/**
	   @return The data pointer.
	*/
	operator const T* ()const{
		return (T*)p;
	}
	/**
	   @param[in] i Index.
	   @return Reference to the i-th element.
	*/
	T& operator()(int i){
		return p[i];
	}
	/**
	   @param[in] i Index.
	   @return Reference to the i-th element.
	*/
	const T& operator()(int i)const{
		return p[i];
	}
	/**
	   @param[in] in Pointer to compare with.
	   @return Whether both point to the same address.
	*/
	bool operator==(const RefP& in){
		return p==in.p;
	}
	/*T* operator+(int i){
		return p+i;
	}
	const T* operator+(int i)const{
		return p+i;
	}*/
};
	
/**
   Generic array of basic types and classes.
   
   Does not defines operations that are different between Numerical and NumCell array.
*/
template <typename T, template<typename> class Dev=Cpu>
class Array:private RefP<T, Dev>{
	typedef RefP<T, Dev> Parent;
protected:
	long nx=0;
	long ny=0;
	using Parent::p;
public:
	string keywords;
	using Parent::deinit;
	//using Parent::operator+;
	using Parent::operator T*;
	using Parent::operator const T*;
	using Parent::operator();
	using Parent::NRef;
	typedef T type;
	/**
	   Destructor. The data is released by RefP.
	*/
	~Array(){
		nx=0;
		ny=0;
	}
	/**
	   @param[in] ix Row index.
	   @param[in] iy Column index.
	   @return Reference to element (ix, iy).
	*/
	T& operator ()(int ix, int iy){
		return p[ix+nx*iy];
	}
	/**
	   @param[in] ix Row index.
	   @param[in] iy Column index.
	   @return Reference to element (ix, iy).
	*/
	const T& operator ()(int ix, int iy)const{
		return p[ix+nx*iy];
	}
	/**
	   Relaxed indexing: an index along a dimension of size 1 is ignored.

	   @param[in] ix Row index.
	   @param[in] iy Column index.
	   @return Reference to the element.
	*/
	T& R(int ix, int iy=0){
		return p[(nx==1?0:ix)+nx*(ny==1?0:iy)];
	}
	/**
	   @param[in] icol Column index.
	   @return Pointer to the start of column icol.
	*/
	T* Col(int icol){
		return p+nx*icol;
	}
	/**
	   @param[in] icol Column index.
	   @return Pointer to the start of column icol.
	*/
	const T* Col(int icol)const{
		return p+nx*icol;
	}

	/**
	   @return Number of rows.
	*/
	long Nx()const{
		return nx;
	}
	/**
	   @return Number of columns.
	*/
	long Ny()const{
		return ny;
	}
	/**
	   @return Total number of elements.
	*/
	long N()const{
		return nx*ny;
	}
	/**
	   @return Whether the array is not empty.
	*/
	explicit operator bool()const{
		return (nx&&ny);
	}

	/**
	   Reallocate the array. Existing data of the same size is kept.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	*/
	void init(long nxi=0, long nyi=1){
		nx=nxi;
		ny=nyi;
		Parent::init(nxi*nyi);
	}
	//Constructors 
	/**
	   Create an empty array.
	*/
	Array(){}
	/**
	   Create an array.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Existing memory to use. Allocated (zeroed) if NULL.
	   @param[in] own Whether to free pi when the last reference is released.
	*/
	Array(long nxi, long nyi=1, T* pi=NULL, int own=1)
		:Parent(nxi*nyi, pi, own), nx(nxi), ny(nyi){}

	/**
	   Create a reference to part of an existing array with offset.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Array to reference.
	   @param[in] offset Offset in number of elements.
	*/
	Array(long nxi, long nyi, const Array& pi, long offset=0):Parent(pi, offset), nx(nxi), ny(nyi){}
	//Use default destructor

	/**
	   Create a reference to an existing array.

	   @param[in] in Array to reference.
	*/
	Array(const Array& in):Parent(in), nx(in.nx), ny(in.ny), keywords(in.keywords){}
	/**
	   Release the current data and reference in.

	   @param[in] in Array to reference.
	   @return Reference to this.
	*/
	Array& operator=(const Array& in){
		if(this!=&in){
			Parent::operator=(in);
			nx=in.nx;
			ny=in.ny;
			keywords=in.keywords;
		}
		return *this;
	}

	/**
	   @return A reference to the same data reshaped as a column vector.
	*/
	Array Vector()const{
		Array tmp=*this;
		tmp.nx=tmp.nx*tmp.ny;
		tmp.ny=1;
		return tmp;
	}
};
///Used for array in cpu of any other array.
template <typename T>
class CellArray:public Array<T, CpuObj>{
	typedef Array<T, CpuObj> Parent;
protected:
	using Parent::nx;
	using Parent::ny;
	using Parent::p;
public:
	using Parent::keywords;
	using Parent::deinit;
	//using Parent::operator+;
	using Parent::operator T *;
	using Parent::operator const T *;
	using Parent::operator();
	using Parent::NRef;
	using Parent::Col;
	using Parent::Nx;
	using Parent::Ny;
	using Parent::N;
	using Parent::operator bool;
	using Parent::init;
	/**
	   @return A new cell array (column vector) whose blocks are reshaped as column vectors.
	*/
	CellArray Vector()const{
		CellArray tmp(this->N(), 1);
		for(long i=0; i<this->N(); i++){
			tmp[i]=(*this)[i].Vector();
		}
		return tmp;
	}
	/**
	   Transpose the cell array.

	   @param[in] stream Stream to run in.
	   @return The transposed cell array.
	*/
	CellArray trans(stream_t &stream);
	//using Parent::Copy;
	/**
	   Create a cell array.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Existing memory to use. Allocated if NULL.
	   @param[in] own Whether to free pi when the last reference is released.
	*/
	CellArray(long nxi, long nyi=1, T *pi=NULL, int own=1):Parent(nxi, nyi, pi, own){};
	/**
	   Create a reference to part of an existing cell array with offset.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Cell array to reference.
	   @param[in] offset Offset in number of elements.
	*/
	CellArray(long nxi, long nyi, const CellArray &pi, long offset=0):Parent(nxi, nyi, pi, offset){};
	/**
	   Create an empty cell array.
	*/
	CellArray():Parent(){};

	/**
	   Zero each block.

	   @param[in] stream Stream to run in.
	*/
	void Zero(cudaStream_t stream=0){
		for(long i=0; i<N(); i++){
			(*this)(i).Zero(stream);
		}
	}

	//Copy the data. Zero data if input data is empty. Reallocate if mismatch.
	/*virtual void Copy(const CellArray &in, cudaStream_t stream=0){
		if(this!=&in){
			if(!in){
				dbg("Copying from empty input, no action is taken.\n");
			} else{
				if(this->N()!=in.N()){
					if(NRef()>1){
						dbg("Copying into referenced data that has different size.\n");
					}
					if(this->N()!=0){
						dbg("Reinitiate data during copying, old size is %ldx%ld, new size is %ldx%ld\n", Nx(), Ny(), in.Nx(), in.Ny());
					}
					//deinit();
					init(in.Nx(), in.Ny());
				}
				for(long i=0; i<N(); i++){
					(*this)(i).Copy(in(i), stream);
				}
			}
		}
	}*/
	/**
	   Write the cell array to file.

	   @param[in] stream Stream to copy in.
	   @param[in,out] fp File to write to.
	*/
	void Write(cudaStream_t stream, file_t *fp) const{
		header_t header={MCC_ANY, (uint64_t)Nx(), (uint64_t)Ny(), NULL};
		write_header(&header, fp);
		for(int i=0; i<N(); i++){
			p[i].Write(stream, fp);
		}
	}
	/**
	   Write the cell array to a file named by format.

	   @param[in] stream Stream to copy in.
	   @param[in] format printf style format of the file name, followed by arguments.
	*/
	void Write(cudaStream_t stream, const char *format, ...) const{
		format2fn;
		file_t *fp=zfopen(fn, "wb");
		Write(stream, fp);
		zfclose(fp);
	}
};

///Only for array of numerical data.
template <typename T, template<typename> class Dev=Cpu>
class NumArray:public Array<T, Dev>{
#if __cpp_static_assert >= 200410L
	static_assert(sizeof(T)<=16, "NumArray can only be used with numerical data");
#endif	
	typedef Array<T, Dev> Parent;
protected:
	using Parent::nx;
	using Parent::ny;
	using Parent::p;
public:
	using Parent::keywords;
	using Parent::deinit;
	//using Parent::operator+;
	using Parent::operator T *;
	using Parent::operator const T *;
	using Parent::operator();
	using Parent::NRef;
	using Parent::Col;
	using Parent::Nx;
	using Parent::Ny;
	using Parent::N;
	using Parent::operator bool;
	using Parent::init;
#if CUDA_VERSION>10000
	/**
	   @return CUDA data type of the elements. Specialized in types.cu.
	*/
	cudaDataType dtype();
	mutable cusparseDnVecDescr_t vdesc=0;//for cuspmul
	mutable cusparseDnMatDescr_t mdesc=0;//for cuspmul
#endif
	/**
	   Transpose the array.

	   @param[in] stream Stream to run in.
	   @return A newly allocated transposed array.
	*/
	NumArray trans(stream_t &stream){
		NumArray B(ny,nx);
		transpose<<<dim3(16, 16), dim3(16, 16), 0, stream>>>
			(B(), p, nx, ny);
		return B;
	}
	//using Parent::Copy;
	/**
	   Create an empty array.
	*/
	NumArray():Parent(){};
	/**
	   Create an array.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Existing memory to use. Allocated (zeroed) if NULL.
	   @param[in] own Whether to free pi when the last reference is released.
	*/
	NumArray(long nxi, long nyi=1, T *pi=NULL, int own=1):Parent(nxi, nyi, pi, own){};
	/**
	   Create a reference to part of an existing array with offset.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] pi Array to reference.
	   @param[in] offset Offset in number of elements.
	*/
	NumArray(long nxi, long nyi, const NumArray& pi, long offset=0):Parent(nxi, nyi, pi, offset){};
	/**
	   C Array wrapper by weak reference. Only specialized for matching type in
	   types.cu.

	   @param[in] A CPU matrix.
	*/
	NumArray(const dmat *A);
	/**
	   C Array wrapper. See above.

	   @param[in] A CPU matrix.
	*/
	NumArray(const cmat *A);
	/**
	   C Array wrapper. See above.

	   @param[in] A CPU matrix.
	*/
	NumArray(const smat *A);
	/**
	   C Array wrapper. See above.

	   @param[in] A CPU matrix.
	*/
	NumArray(const zmat *A);
#if CUDA_VERSION>10000
	/**
	   Create a reference to an existing array. The cusparse descriptors are not
	   copied.

	   @param[in] in Array to reference.
	*/
	NumArray(const NumArray &in):Parent(in){}
	/**
	   Reference in. The cusparse descriptors are not copied.

	   @param[in] in Array to reference.
	   @return Reference to this.
	*/
	NumArray &operator=(const NumArray &in){
		if(this!=&in){
			Parent::operator=(in);
			vdesc=NULL;//don't copy vdesc, mdesc
			mdesc=NULL;
		}
		return *this;
	}
	/**
	   Destroy the cusparse descriptors.
	*/
	~NumArray(){
		if(vdesc) {DO(cusparseDestroyDnVec(vdesc)); vdesc=NULL;}
		if(mdesc) {DO(cusparseDestroyDnMat(mdesc)); mdesc=NULL;}
	}
#endif
	/**
	   @return A reference to the same data reshaped as a column vector.
	*/
	NumArray Vector()const{
		NumArray tmp=*this;
		tmp.nx=tmp.nx*tmp.ny;
		tmp.ny=1;
		return tmp;
	}
	//Implementations that are different from regular Array.
	/**
	   Zero the data.

	   @param[in] stream Stream to run in.
	*/
	void Zero(cudaStream_t stream=0){
		Dev<T>::Zero(p, N(), stream);
	}
	/**
	   Set every element to alpha.

	   @param[in] alpha Value to set.
	   @param[in] stream Stream to run in.
	*/
	void Set(T alpha, cudaStream_t stream=0){
		if(N()){
			set_do<<<DIM(N(), 256), 0, stream>>>(p(), alpha, N());
		}
	}
	/**
	   Copy the data. Zero data if input data is empty. Reallocate if mismatch.

	   @param[in] in Array to copy from.
	   @param[in] stream Stream to run in.
	*/
	template<typename S>
	void Copy(const NumArray<S, Dev> &in, cudaStream_t stream=0){
		if((void*)this!=(void*)&in){
			if(!in){
				dbg("Copying from empty input, no action.\n");
				this->Zero(stream);
			} else{
				if(this->N()!=in.N()){
					if(NRef()>1){
						dbg("Copying into referenced data that has different size.\n");
					}
					if(this->N()!=0){
						dbg("Reinitiate data during copying, old size is %ldx%ld, new size is %ldx%ld\n", Nx(), Ny(), in.Nx(), in.Ny());
					}
					//edinit()
					init(in.Nx(), in.Ny());
				}
				Dev<T>::Copy((*this)(), in(), in.N(), stream);
			}
		}
	}
	/**
	   Write the array to file.

	   @param[in] stream Stream to copy in.
	   @param[in,out] fp File to write to.
	*/
	void Write(cudaStream_t stream, file_t *fp) const{
		const char *keys=keywords.length()?keywords.c_str():NULL;
		Dev<T>::Write(p, Nx(), Ny(), keys, fp, stream);
	}
	/**
	   Write the array to a file named by format.

	   @param[in] stream Stream to copy in.
	   @param[in] format printf style format of the file name, followed by arguments.
	*/
	void Write(cudaStream_t stream, const char *format, ...) const{
		format2fn;
		file_t *fp=zfopen(fn, "wb");
		Write(stream, fp);
		zfclose(fp);
	}
};

/**
   NumCell is a special Array of array of numerical data that can store multiple Arrays of numerical data in continuous region.
*/

template <typename T, template<typename> class Dev=Cpu>
class NumCell:public CellArray<NumArray<T, Dev> >{
private:
	typedef NumArray<T, Dev> TMat;
	typedef CellArray<NumArray<T, Dev> > Parent;
	TMat m; /*contains the continuous data*/
	Array<T *, Pinned>pm_cpu;/*contains the data pointer in each cell in gpu.*/
protected:
	using Parent::p;
	using Parent::nx;
	using Parent::ny;
public:
	Array<T*, Gpu>pm;/*contains the data pointer in each cell in cpu.*/
	using Parent::operator();
	using Parent::N;
	using Parent::Nx;
	using Parent::Ny;
	using Parent::NRef;
	//using Parent::Copy;
	using Parent::Zero;
	using Parent::Write;
	/**
	   @return Reference to the continuous data array. Empty if blocks are not continuous.
	*/
	TMat& M(){
		return m;
	}
	/**
	   @return Reference to the continuous data array. Empty if blocks are not continuous.
	*/
	const TMat& M() const{
		return m;
	}

	/**
	   Update the device array pm of block data pointers.

	   @param[in] stream Stream to copy in.
	*/
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

	/**
	   Create a cell array with empty blocks.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	*/
	NumCell(long nxi=0, long nyi=1):Parent(nxi, nyi){
		p2pm();
	}
	/**
	   Create a cell array with blocks of the same size in continuous memory.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] mx Number of rows of each block.
	   @param[in] my Number of columns of each block.
	   @param[in] pin Existing memory to use (not owned). Allocated if NULL.
	*/
	NumCell(long nxi, long nyi, long mx, long my, T* pin=NULL):Parent(nxi, nyi){
		if(mx&&my){
			m=TMat(mx*my*nxi*nyi, 1, pin, 0);
			for(int i=0; i<nxi*nyi; i++){
				p[i]=TMat(mx, my, m, i*(mx*my));
			}
			p2pm();
		}
	}
	/**
	   Create a cell array with blocks of different sizes in continuous memory.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	   @param[in] mx Number of rows of each block. Blocks with 0 are left empty.
	   @param[in] my Number of columns of each block. 1 if NULL.
	   @param[in] pin Existing memory to use (not owned). Allocated if NULL.
	*/
	template <typename L>
	NumCell(long nxi, long nyi, L* mx, L* my, T* pin=NULL):Parent(nxi, nyi){
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

	/**
	   Create a reference to an existing cell array.

	   @param[in] in Cell array to reference.
	*/
	NumCell(const NumCell& in):Parent(in), m(in.m), pm_cpu(in.pm_cpu), pm(in.pm){}

	/**
	   Reference in.

	   @param[in] in Cell array to reference.
	   @return Reference to this.
	*/
	NumCell& operator=(const NumCell& in){
		if(this!=&in){
			Parent::operator=(in);
			m=in.m;
			pm_cpu=in.pm_cpu;
			pm=in.pm;
		}
		return *this;
	}
	//Use default copy constructor and copy assignment operator
	/**
	   Reallocate the cell array with empty blocks.

	   @param[in] nxi Number of rows.
	   @param[in] nyi Number of columns.
	*/
	void init(long nxi=0, long nyi=1){
		Parent::init(nxi, nyi);
		m=TMat();//blocks are no longer in continuous memory
		p2pm();
	}
	/**
	   @return A newly allocated cell array with the same layout (in continuous memory if this one is).
	*/
	NumCell New()const{
		if(!m){
			return NumCell(nx, ny);
		} else{
			NumArray<long>mx(nx*ny);
			NumArray<long>my(nx*ny);
			for(long i=0; i<nx*ny; i++){
				mx[i]=p[i].Nx();
				my[i]=p[i].Ny();
			}
			return NumCell(nx, ny, mx(), my());
		}
	}

	/**
	   Replace the data with a new set. We don't own pnew.

	   @param[in] pnew New data with the same total size, laid out block after block.
	   @param[in] stream Stream to update pm in.
	*/
	void replace(T* pnew, cudaStream_t stream){
		if(m){
			m=TMat(m.Nx(), m.Ny(), pnew, 0);//replace m
		}

		for(long i=0; i<nx*ny; i++){
			p[i]=TMat(p[i].Nx(), p[i].Ny(), pnew, 0);//don't own the data pnew
			pnew+=p[i].N();
		}
		p2pm(stream);
	}
	/**
	   @return A new cell array (column vector) whose blocks are reshaped as column vectors, sharing the data.
	*/
	NumCell Vector()const{
		NumCell tmp(this->N(), 1);
		for(long i=0; i<this->N(); i++){
			tmp[i]=(*this)[i].Vector();
		}
		tmp.pm=this->pm;
		tmp.pm_cpu=this->pm_cpu;
		tmp.m=this->m;
		return tmp;
	}
	/**
	   Copy the data. Zero data if input data is empty. Reallocate if size
	   mismatch. The continuous data M() is copied if both have it with the
	   same size, otherwise the data is copied block by block.

	   @param[in] in Cell array to copy from.
	   @param[in] stream Stream to run in.
	*/
	void Copy(const NumCell &in, cudaStream_t stream=0){
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
					*this=in.New();
				}
				if(m&&in.m&&m.N()==in.m.N()){
					m.Copy(in.m, stream);
				} else{
					for(long i=0; i<N(); i++){
						p[i].Copy(in.p[i], stream);
					}
					p2pm(stream);
				}
				//dbg("NumCell: Copy, m is %ldx%ld\n", m.Nx(), m.Ny());
			}
		}
	}
};
typedef class NumArray<int, Gpu>   cuimat;
typedef class NumArray<Real, Gpu>   curmat;
typedef class NumArray<Comp, Gpu>   cucmat;
typedef class NumArray<real, Gpu>   cudmat;//same dtype as cpu dmat
typedef class NumArray<comp, Gpu>   cuzmat;//same dtype as cpu cmat
typedef class NumCell<Real, Gpu>  curcell;
typedef class NumCell<Comp, Gpu>  cuccell;
typedef class NumArray<real>   crmat;//equivalent to rmat in c
typedef class NumArray<comp>   ccmat;//equivalent to cmat in c

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
	/**
	   Free the device data, buffers and descriptor.
	*/
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
	/**
	   @return Storage format (SP_CSC or SP_CSR).
	*/
	enum TYPE_SP Type()const{
		return ref->type;
	}
	/**
	   @return Number of rows.
	*/
	long Nx()const{
		return ref->nx;
	}
	/**
	   @return Number of columns.
	*/
	long Ny()const{
		return ref->ny;
	}
	/**
	   @return Number of nonzeros.
	*/
	long Nzmax()const{
		return ref->nzmax;
	}
	/**
	   @return Device pointer to the column (CSC) or row (CSR) pointers.
	*/
	Spint* Pp(){
		return ref->p;
	}
	/**
	   @return Device pointer to the column (CSC) or row (CSR) pointers.
	*/
	const Spint* Pp() const{
		return ref->p;
	}
	/**
	   @return Device pointer to the row (CSC) or column (CSR) indices.
	*/
	Spint* Pi(){
		return ref->i;
	}
	/**
	   @return Device pointer to the row (CSC) or column (CSR) indices.
	*/
	const Spint* Pi() const{
		return ref->i;
	}

	/**
	   @return Device pointer to the values.
	*/
	Real* Px(){
		return ref->x;
	}
	/**
	   @return Device pointer to the values.
	*/
	const Real* Px() const{
		return ref->x;
	}
	/**
	   Create an empty sparse matrix.
	*/
	cusp(){}
	/** Copy a CPU sparse matrix to GPU. See utils.cu. */
	cusp(const dsp* src_csc, int tocsr=0, int transp=0);
	/**
	   Create a reference to an existing sparse matrix.

	   @param[in] in Sparse matrix to reference.
	*/
	cusp(const cusp& in):ref(in.ref){
		if(ref) ref->count++;
	}
	/**
	   Release the current data and reference in.

	   @param[in] in Sparse matrix to reference.
	   @return Reference to this.
	*/
	cusp& operator=(const cusp& in){
		if(this!=&in && this->ref!=in.ref){
			deinit();
			ref=in.ref;
			if(ref) ref->count++;
		}
		return *this;
	}
	/**
	   Release the reference. See deinit().
	*/
	~cusp(){
		deinit();
	}
	/**
	   Release the reference and free the data if this is the last one.
	*/
	void deinit(){
		if(ref && !atomic_sub_fetch(&ref->count, 1)){
			delete ref;
		}
	}

	/**
	   Convert to CSR mode by transpose.
	*/
	void trans();

	/**
	   @return Whether the sparse matrix is not empty.
	*/
	explicit operator bool()const{
		return ref && ref->nx&&ref->ny;
	}
	/**
	   Not supported. Errors out.

	   @param[in] in Unused.
	   @param[in] stream Unused.
	*/
	void Copy(const cusp& in, cudaStream_t stream=0){
		error("Copy is not supported\n");
	}
};
typedef CellArray<cusp> cuspcell;
typedef CellArray<curcell> curccell;
typedef CellArray<curccell> curcccell;
class culoc_t{
private:
	curmat p;/*2*n. transposed from CPU DMAT(loc)*/
	Real dx;
	Real dy;
public:
	real xmax, xmin, ymax, ymin;
	/**
	   @return Device pointer to the (x, y) coordinates.
	*/
	Real2* operator()(){
		return (Real2*)(p());
	}
	/**
	   @return Device pointer to the (x, y) coordinates.
	*/
	const Real2* operator()()const{
		return (Real2*)(p());
	}
	/**
	   @return Number of points.
	*/
	long Nloc()const{
		return p.Ny();
	}
	/**
	   @return Sampling along x.
	*/
	Real Dx()const{
		return dx;
	}
	/**
	   @return Sampling along y.
	*/
	Real Dy()const{
		return dy;
	}
	//No need custom copy assignment operator or copy constructor.
	/** Copy loc_t to GPU. See types.cu. */
	culoc_t(const loc_t* in=0);
};
class cupts_t:public culoc_t{
	Real dxsa;
	long nxsa;
public:
	/**
	   @return Subaperture size.
	*/
	Real Dxsa(){
		return dxsa;
	}
	/**
	   @return Number of points along each dimension of a subaperture.
	*/
	long Nxsa(){
		return nxsa;
	}
	/**
	   Copy pts_t to GPU.

	   @param[in] in CPU pts. An empty object is created if NULL.
	*/
	cupts_t(pts_t* in=0):culoc_t((loc_t*)in), dxsa(0), nxsa(0){
		if(in){
			dxsa=in->dx;
			nxsa=in->nxsa;
		}
	}
};
/**
   Compute the cubic spline influence function coefficients. See types.cu.

   @param[in] iac Inter-actuator coupling.
   @return The coefficients on GPU, or an empty array if iac is 0.
*/
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
	Real dratio;/**<M1 to DM demagnification ratio (>1 for shrinking beam)*/
	Real aoi;   /**<Angle of incidence in radian for ray tracing input */
	curmat cubic_cc; /*coefficients for cubic influence function. */
	//use default copy assignment operator and copy constructor

	/**
	   Copy the grid parameters (not the data) from a CPU map.

	   @param[in] in CPU map. Nothing is done if NULL.
	   @return Reference to this.
	*/
	cugrid_t& operator=(const map_t* in){
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
			dratio=in->dratio;
			aoi=in->aoi;
			cubic_cc=iac2cc(in->iac);
		}
		return *this;
	}
	/**
	   Create a grid from a CPU map.

	   @param[in] in CPU map. Default values are used if NULL.
	*/
	cugrid_t(const map_t* in=NULL):nx(0), ny(0), ox(0), oy(0), dx(0), dy(0), ht(0), vx(0), vy(0),dratio(1),aoi(0){
		if(in){
			(*this).operator=(in);
		}
	}
	/**
	   @param[in] sc Scaling factor.
	   @return A copy of the grid with origin and sampling scaled by sc.
	*/
	cugrid_t Scale(Real sc)const{
		cugrid_t tmp(*this);
		tmp.ox*=sc;
		tmp.oy*=sc;
		tmp.dx*=sc;
		tmp.dy*=sc;
		return tmp;
	}
	/**
	   @return Whether the grid is not empty.
	*/
	explicit operator bool(){
		return nx&&ny;
	}
};
class cumap_t:public cugrid_t{
public:
	curmat p;
	/**
	   Create an empty map.
	*/
	cumap_t(){}
	/**
	   Create a map with the grid and allocate the data.

	   @param[in] grid Grid parameters.
	*/
	cumap_t(cugrid_t& grid):cugrid_t(grid){
		if(nx>0&&ny>0){
			p=curmat(nx, ny);
		}
	}
	/**
	   Create a map with the grid referencing existing data.

	   @param[in] grid Grid parameters.
	   @param[in] pin Data to reference.
	*/
	cumap_t(cugrid_t &grid, curmat &pin):cugrid_t(grid), p(pin){
	}
	/**
	   Copy grid parameters from a CPU map and reallocate the data if the size
	   differs. The data is not copied.

	   @param[in] in CPU map. Nothing is done if NULL.
	   @return Reference to this.
	*/
	cumap_t& operator=(map_t* in){
		if(in){
			cugrid_t::operator=(in);
			if(p.Nx()!=nx||p.Ny()!=ny){
				p=curmat(nx, ny);
			}
		}
		return *this;
	}
	/**
	   @return Reference to the data.
	*/
	operator const curmat& ()const{
		return p;
	}
	/**
	   @return Reference to the data.
	*/
	operator curmat& (){
		return p;
	}
	/**
	   @return Device pointer to the data.
	*/
	Real* operator()(){
		return p();
	}
	/**
	   @return Device pointer to the data.
	*/
	const Real* operator()()const{
		return p();
	}
	/**
	   Zero the data.

	   @param[in] stream Stream to run in.
	*/
	void Zero(cudaStream_t stream=0){
		p.Zero(stream);
	}
	/**
	   Copy grid parameters and data from in.

	   @param[in] in Map to copy from.
	   @param[in] stream Stream to run in.
	*/
	void Copy(const cumap_t& in, cudaStream_t stream=0){
		cugrid_t::operator=(in);
		p.Copy(in.p, stream);
	}
};

typedef CellArray<cumap_t> cumapcell;
typedef CellArray<cugrid_t> cugridcell;
/**
   Copy cell array block by block. We do not use member method here because
   some type do not define a Copy method. out is zeroed if in is empty and
   reallocated if size mismatch.

   @param[in,out] out Destination cell array.
   @param[in] in Source cell array.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Copy(CellArray<T> &out, const CellArray<T>& in, cudaStream_t stream=0){
	if(&out!=&in){
		if(!in.N()){
			out.Zero(stream);
		}else{
			if(out.N()!=in.N()){
				out.init(in.Nx(), in.Ny());
			}
			for(long i=0; i<in.N(); i++){
				out[i].Copy(in[i], stream);
			}
		}
	}
}
/**
   Copy between same device. See NumArray::Copy().

   @param[in,out] out Destination array. Reallocated if size mismatch.
   @param[in] in Source array.
   @param[in] stream Stream to run in.
*/
template <typename T, typename T2, template<typename> class Dev>
void Copy(NumArray<T, Dev> &out, const NumArray<T2, Dev> &in, cudaStream_t stream=0){
	out.Copy(in, stream);
}

/**
   Copy between same device. See NumCell::Copy().

   @param[in,out] out Destination cell array. Reallocated if size mismatch.
   @param[in] in Source cell array.
   @param[in] stream Stream to run in.
*/
template <typename T, typename T2, template<typename> class Dev>
void Copy(NumCell<T, Dev> &out, const NumCell<T2, Dev> &in, cudaStream_t stream=0){
	out.Copy(in, stream);
}
//Copy between different device, same or different type
/*template <typename T, template<typename> class Dev, typename T2, template<typename> class Dev2>
void Copy(NumArray<T, Dev> &out, const NumArray<T2, Dev2> &in, cudaStream_t stream=0){
	cp2gpu(out(), in(), in.Nx(), in.Ny(), stream);
}*/
/**
   Copy from CPU array to GPU array with type conversion. See cp2gpu().

   @param[out] out Destination GPU array. Reallocated if size mismatch.
   @param[in] in Source CPU array.
   @param[in] stream Stream to copy in.
*/
template <typename T, typename T2>
void Copy(NumArray<T, Gpu> &out, const NumArray<T2> &in, cudaStream_t stream=0){
	if(out.N()!=in.N()){
		out.init(in.Nx(), in.Ny());
	}
	cp2gpu(out(), in(), in.Nx(), in.Ny(), stream);
}
/*static inline void Copy(cumap_t&out, const cumap_t &in, cudaStream_t stream=0){
	out.Copy(in, stream);
}*/
/**
   Scale a vector: out*=beta. Nothing is done if beta is 1.

   @param[in,out] out Array to scale. Nothing is done if empty.
   @param[in] beta Scaling factor.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Scale(NumArray<T, Gpu> &out, Real beta, cudaStream_t stream=0){
	if(!out) return;
	if(Z(fabs)(beta-(Real)1)>EPS){
		scale_do<<<DIM(out.N(), 256), 0, stream>>>
			(out(), out.N(), beta);
	}
}
/**
   Add of the same type: out=out*alpha+in*beta.

   @param[in,out] out Destination array. Set to in*beta if empty or alpha is 0.
   @param[in] alpha Scaling factor of out.
   @param[in] in Source array. Nothing is done if empty.
   @param[in] beta Scaling factor of in.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Add(NumArray<T, Gpu> &out, T alpha, const NumArray<T, Gpu> &in, T beta, cudaStream_t stream=0){
	if(!in) return;
	if(!out||alpha==0){
		Copy(out, in, stream);
		Scale(out, beta, stream);
	} else{
		assert(in.N()==out.N());
		add_do<<<DIM(in.N(), 256), 0, stream>>>
			(out(), (T*)NULL, alpha, in(), (T*)NULL, beta, in.N());
	}
}
/**
   Add a scalar to every element: A+=beta.

   @param[in,out] A Array.
   @param[in] beta Value to add.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Add(NumArray<T, Gpu> &A, T beta, cudaStream_t stream=0){
	add_do<<<DIM(A.N(), 256), 0, stream>>>(A(), beta, A.N());
}
/**
   out=out+in*(*alpha)*alpha2, where alpha lives in device memory.

   @param[in,out] out Destination array. Allocated (zeroed) if empty.
   @param[in] in Source array.
   @param[in] alpha Device pointer to scaling factor of in. NULL is treated as 1.
   @param[in] alpha2 Additional scaling factor of in.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Add(NumArray<T, Gpu> &out, const NumArray<T, Gpu> &in, T *alpha, T alpha2, cudaStream_t stream=0){
	if(!out){
		out=NumArray<T, Gpu>(in.Nx(), in.Ny());
	}
	add_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), in(), alpha, alpha2, in.N());
}
/**
   out=out*(*alpha1)+in, where alpha1 lives in device memory.

   @param[in,out] out Destination array. Allocated (zeroed) if empty.
   @param[in] alpha1 Device pointer to scaling factor of out. NULL is treated as 1.
   @param[in] in Source array.
   @param[in] stream Stream to run in.
*/
template <typename T>
void Add(NumArray<T, Gpu> &out, T *alpha1, const NumArray<T, Gpu> &in, cudaStream_t stream=0){
	if(!out){
		out=NumArray<T, Gpu>(in.Nx(), in.Ny());
	}
	add_do<<<DIM(in.N(), 256), 0, stream>>>
		(out(), alpha1, (T)1, in(), in.N());
}
/**
   Zero array if exist, otherwise allocate and zero.

   @param[in,out] A Array.
   @param[in] nx Number of rows if allocated.
   @param[in] ny Number of columns if allocated.
*/
template <typename T, template<typename> class Dev >
void initzero(NumArray<T, Dev>& A, long nx, long ny){
	if(A){
		A.Zero();
	} else{
		A.init(nx, ny);
	}
}

#define cuzero(A,B...) (A).Zero(B)
#define cuwrite(A,B...) A.Write(B)
#endif

