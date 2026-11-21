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
#include <errno.h>
#include "random.h"
#include "mathdef.h"
#include "defs.h"/*Defines T, X, etc */

/**
   The only function that allocates memory for the matrix object. It ensures that all
   fields are properly initialized. 
   - If p is NULL, mem should be NULL as well. memory is allocated. if mem is not null, it is ignored and warning is print. new mem_t is created.
   - If p is not NULL, mem is NULL, we act as a weak pointer. p is managed externally. 
   - If both p and mem are not NULL, we increment mem reference count.

   @param[in,out] pout The matrix. Allocated if NULL. If it already has data, only the dimension is checked.
   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @param[in] p Existing memory to wrap. NULL to allocate.
   @param[in] mem Reference counter of p. May be NULL.
   @return 0 on success, -1 on failure.
*/
static int X(init_do)(X(mat) **pout, long nx, long ny, T *p, mem_t *mem){
	if(!pout){
		warning("Null pointer for pout\n");
		return -1;
	}
	if(!*pout&&!(*pout=mycalloc(1, X(mat)))){
		warning("malloc for X(mat) failed\n");
		return -1;
	}
	X(mat)*out=*pout;
	if(nx<=0||ny<=0){
		nx=0;
		ny=0;
	}
	
	if(out->id && (out->id&M_T)!=M_T){
		error("Mismatch: pout is of type %u, expect %u.\n", out->id, M_T);
		return -1;
	}else if(out->id && out->p){//data exists. check dimension
		if(out->nx!=nx || out->ny!=ny){
			error("Mismatch: pout is %ldx%ld, expect %ldx%ld\n", out->nx, out->ny, nx, ny);
			return -1;
		}//else: just return
	}else{
		if(!p && nx && ny){//no pointer is supplied
			if(!(p=mycalloc((nx*ny), T))){
				error("malloc for %ldx%ld of size %ld failed.\n", nx, ny, sizeof(T));
			}
			if(mem){
				error("mem=%p should be NULL when p is NULL. It will be overriden.\n", mem);
			}
			mem=mem_new(p);
		}
		
		out->id=M_T;
		out->nx=nx;
		out->ny=ny;
		out->p=p;
		if(mem) out->mem=mem_ref(mem);
	}
	return 0;
}
/**
   check the size of matrix if exist. Otherwise create it. content is not zeroed if already exists.
   - empty matrix is initialized.
   - existing matrix is not resized.

   @param[in,out] A The matrix. Allocated if NULL; otherwise its dimension must match.
   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @return 0 on success, -1 on failure.
*/
int X(init)(X(mat) **A, long nx, long ny){
	return X(init_do)(A, nx, ny, NULL, 0);
}
/**
   Create a new X(mat) of nx*ny. If p is not NULL, wraps the existing memory p
   with reference counting by mem (may be NULL) instead of allocating.

   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @param[in] p Existing memory to wrap. NULL to allocate.
   @param[in] mem Reference counter of p. May be NULL for a weak reference.
   @return Newly allocated matrix; the caller frees it.
*/
X(mat) *X(new_do)(long nx, long ny, T *p, mem_t *mem){
	X(mat)* out=NULL;
	X(init_do)(&out, nx, ny, p, mem);
	return out;
}
/**
   Create a new X(mat). Initialized to zero.

   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @return Newly allocated matrix; the caller frees it.
*/
 X(mat)* X(new)(long nx, long ny){
	return X(new_do)(nx, ny, NULL, 0);
}
/**
   Calls X(new) with a filename to be saved to. The data is written to the file
   asynchronously, or when it is freed for fits files.

   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @param[in] keywords Keywords to save with the data. May be NULL.
   @param[in] format printf style format of the file name. A leading - disables saving.
   @return Newly allocated matrix; the caller frees it. NULL if nx or ny is 0.
*/
X(mat)* X(new_file)(long nx, long ny, const char* keywords, const char* format, ...){
	if(!nx||!ny) return NULL;
	format2fn;
	if(fn&&fn[0]=='-') fn=NULL;//leading - disables filename.
	X(mat)* out=X(new)(nx, ny);
	if(out&&keywords) out->keywords=strdup(keywords);
	if(disable_save&&!IS_SHM(fn))fn=NULL;
	if(out&&fn) {
		out->fp=zfopen(fn, "w");
		writedata(out->fp, out, -1);
		//out->async=async_init(out->fp, sizeof(T), M_T, out->keywords, P(out), out->nx, out->ny);
	}
	return out;
}
/**
   check and cast an object to matrix

   @param[in] A The object.
   @return A cast to matrix if it is a matrix of this type, otherwise NULL. Not a new reference.
*/
X(mat)* X(mat_cast)(const_anyarray A){
	return ismat(A.c)?(X(mat)*)A.c:NULL;
}
/**
   free content of a matrix object.

   @param[in,out] A The matrix. Its content is freed and the struct is kept.
*/
void X(free_content)(X(mat)* A){
	if(check_mat(A)){
		//only need to initiate write when fp is set.
		//if(A->fp) X(writedata)(A->fp, A, A->ny);//don't do this. 
		if(A->fp && !A->async && zfisfits(A->fp)){//async writing is not supported for fits. write all data now.
			X(writedata)(A->fp, A, 0);
		}
		if(A->async){//make sure entire data is written.
			async_write(A->async, PN(A)*sizeof(T), 1);
			async_free(A->async); A->async=NULL;
		}
		if(A->fp) {zfclose(A->fp); A->fp=NULL;}
		mem_unref(&A->mem);//takes care of freeing memory.
		if(A->deinit) A->deinit((cell*)A);//for derived types
#ifndef COMP_LONG
		if(A->fft) {X(fft_free_plan)(A->fft); A->fft=NULL;}
#endif
		if(A->keywords) {free(A->keywords); A->keywords=NULL;}
		A->nx=0;
		A->ny=0;
		A->p=0;
	}else{
		warning("check_mat returns false\n");
	}
}
/**
   free a matrix object.

   @param[in] A The matrix to free.
*/
void X(free_do)(X(mat) *A){
	X(free_content)(A);
	free(A);
}
/**
   create an new reference another with different shape.
   - total number of elements must be preserved
   - set ny to 1 or 0 to return a vector.

   @param[in] in The matrix to reference.
   @param[in] nx New number of rows.
   @param[in] ny New number of columns.
   @return Newly allocated matrix sharing data with in; the caller frees it. NULL if in is invalid.
*/
X(mat)* X(ref_reshape)(const X(mat)* in, long nx, long ny){
	X(mat) *out=NULL;
	if(check_mat(in)){
		out=X(ref)(in);
		if(reshape(out, nx, ny)){//failed
			error("Dimension mismatch: (%ldx%ld) vs (%ldx%ld)\n", NX(in), NY(in), nx, ny);
			X(free)(out);
		}
	}
	return out;
}

/**
   creat a reference to an existing matrix. header is duplicated if exists.

   @param[in] in The matrix to reference.
   @return Newly allocated matrix sharing data with in; the caller frees it. NULL if in is invalid.
*/
X(mat) *X(ref)(const X(mat) *in){
	if(!check_mat(in)) return NULL;
	X(mat) *out=X(new_do)(in->nx, in->ny, P(in), in->mem);
	if(out&&in->keywords) out->keywords=strdup(in->keywords);
	return out;
}

/**
   creat a new matrix referencing columns in existing
   matrix. reference counted. 

   @param[in] in The matrix to reference.
   @param[in] icol Index of the first column.
   @param[in] ncol Number of columns.
   @return Newly allocated matrix sharing data with in; the caller frees it. NULL if in is invalid.
*/
X(mat)* X(refcols)(const X(mat)* in, long icol, long ncol){
	if(!check_mat(in)) return NULL;
	if(icol<0 || ncol<0 || (icol+ncol)>NY(in)){
		error("Invalid parameters: icol=%ld, ncol=%ld, ny=%ld", icol, ncol, NY(in));
		return NULL;
	}
	return X(new_do)(in->nx, ncol, PCOL(in, icol), in->mem);
}
/**
   Override a stack matrix struct with pointers to columns of another matrix. Does not add reference count the original data.

   @param[in,out] out The matrix struct. Its existing content is freed if initialized.
   @param[in] in The matrix to reference.
   @param[in] icol Index of the first column.
   @param[in] ncol Number of columns.
*/
void X(cols)(X(mat) *out, const X(mat) *in, long icol, long ncol){
	if(!check_mat(in)||!out) return;
	if(icol<0||ncol<0||(icol+ncol)>NY(in)){
		error("Invalid parameters: icol=%ld, ncol=%ld, ny=%ld", icol, ncol, NY(in));
		return;
	}
	if(out->id) X(free_content)(out);	
	out->id=M_T;
	out->nx=in->nx;
	out->ny=ncol;
	out->p=PCOL(in, icol);
}
/**
   Create a new sub matrix of nx*ny starting from(sx,sy) and copy the data.

   @param[in] in The input matrix.
   @param[in] sx Starting row.
   @param[in] nx Number of rows. <=0 to use all remaining rows.
   @param[in] sy Starting column.
   @param[in] ny Number of columns. <=0 to use all remaining columns.
   @return Newly allocated matrix; the caller frees it. NULL if in is invalid or out of range.
*/
X(mat)* X(sub)(const X(mat)* in, long sx, long nx, long sy, long ny){
	if(!check_mat(in)) return NULL;
	if(nx<=0){
		nx=in->nx-sx;
	}
	if(ny<=0){
		ny=in->ny-sy;
	}
	if(sx+nx>in->nx||sy+ny>in->ny){
		error("Invalid parameter range: (%ld-%ld)x(%ld-%ld) is outside of %ldx%ld\n",
			sx, sx+nx, sy, sy+ny, in->nx, in->ny);
		return NULL;
	} else{
		X(mat)* out=X(new)(nx, ny);
		for(int iy=0; iy<ny; iy++){
			memcpy(PCOL(out, iy), &P(in, sx, iy+sy), sizeof(T)*nx);
		}
		return out;
	}
}
/**
   Resize a matrix by adding or removing columns or rows. Data is kept whenever
   possible.
   - -1 in dimension preserves original value.
   - 0 in dimension will free array.

   @param[in,out] A The matrix resized in place. Must not be referenced elsewhere.
   @param[in] nx New number of rows.
   @param[in] ny New number of columns.
   @return 0 on success, -1 on failure.
*/
int X(resize)(X(mat)* A, long nx, long ny){
	if(A){
		if(!nx || !ny){
			dbg("nx=%ld, ny=%ld will free array\n", nx, ny);
			X(free_content)(A);
		}
		if(!ismat(A)){
			error("Incorrect type: id=%d, cancelled\n", A->id);
			return -1;
		}
	}else{
		if(!nx || !ny){
			return 0;
		}else{
			error("Trying to resize NULL array to %ldx%ld, cancelled.\n", nx, ny);
			return -1;
		}
	}

	if(nx<0){
		nx=A->nx;
	}
	if(ny<0){
		ny=A->ny;
	}
	if(P(A) && mem_nref(A->mem)!=1){
		error("Resizing referenced matrix or weak reference is invalid. Cancelled.\n");
		return -1;
	}
	if(A->nx!=nx||A->ny!=ny){
		if(A->nx==nx||(A->ny==1&&(ny==1||nx>=A->nx))){//layout is preserved by linear realloc
			P(A)=myrealloc(P(A), nx*ny, T);
			if(nx*ny>PN(A)){
				memset(P(A)+PN(A), 0, (nx*ny-PN(A))*sizeof(T));
			}
		} else {//move memory content
			long nymin=MIN(ny, A->ny);
			if(A->nx<nx){//expand row count
				P(A)=myrealloc(P(A), nx*ny, T);//expand memory size
				for(long icol=nymin-1; icol>0; icol--){
					memmove(P(A)+icol*nx, P(A)+icol*A->nx, A->nx*sizeof(T));
				}
				for(int icol=0; icol<nymin; icol++){//zero the new rows
					memset(P(A)+icol*nx+A->nx, 0, (nx-A->nx)*sizeof(T));
				}
			}else{//reduce row count
				for(long icol=1; icol<nymin; icol++){
					memmove(P(A)+icol*nx, P(A)+icol*A->nx, nx*sizeof(T));
				}
				P(A)=myrealloc(P(A), nx*ny, T);//reduce memory size
			}
			if(ny>nymin){//zero the new columns
				memset(P(A)+nymin*nx, 0, (ny-nymin)*nx*sizeof(T));
			}
		} /*else{//
			T* p=mycalloc(nx*ny, T);
			if(!p){
				error("malloc for %ldx%ld of size %ld failed.\n", nx, ny, sizeof(T));
				return -1;
			}
			long minx=A->nx<nx?A->nx:nx;
			long miny=A->ny<ny?A->ny:ny;
			for(long iy=0; iy<miny; iy++){
				memcpy(p+iy*nx, PCOL(A, iy), sizeof(T)*minx);
			}
			free(P(A));
			P(A)=p;
		}*/
		if(A->mem){
			mem_replace(A->mem, P(A));
		} else{
			A->mem=mem_ref(mem_new(P(A)));
		}
		A->nx=nx;
		A->ny=ny;
	}
	return 0;
}

/**
   Concatenate two matrixes into 1 along "dim" (1 for x, 2 for y)

   @param[in] in1 The first matrix.
   @param[in] in2 The second matrix.
   @param[in] dim 1 to concatenate along x (rows), 2 along y (columns).
   @return Newly allocated matrix; the caller frees it. A copy of the other if one is empty.
*/
X(mat)* X(cat)(const X(mat)* in1, const X(mat)* in2, int dim){
	if(!check_mat(in2)){
		if(check_mat(in1)){
			return X(dup)(in1);
		} else{
			return NULL;
		}
	} else if(!check_mat(in1)){
		return X(dup)(in2);
	}
	X(mat)* out=NULL;

	if(dim==1){ /*along x. */
		if(in1->ny!=in2->ny){
			error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld), cancelled.\n",
				in1->nx, in1->ny, in2->nx, in2->ny);
		} else{
			out=X(new)(in1->nx+in2->nx, in1->ny);
			for(long iy=0; iy<in1->ny; iy++){
				memcpy(PCOL(out, iy), PCOL(in1, iy), in1->nx*sizeof(T));
				memcpy(PCOL(out, iy)+in1->nx, PCOL(in2, iy), in2->nx*sizeof(T));
			}
		}
	} else if(dim==2){	/*along y. */
		if(in1->nx!=in2->nx){
			error("Mismatch. in1 is (%ld, %ld), in2 is (%ld, %ld), cancelled.\n",
				in1->nx, in1->ny, in2->nx, in2->ny);
		} else{
			out=X(new)(in1->nx, in1->ny+in2->ny);
			memcpy(P(out), P(in1), in1->nx*in1->ny*sizeof(T));
			memcpy(P(out)+PN(in1), P(in2), in2->nx*in2->ny*sizeof(T));
		}
	} else{
		error("Invalid dim, excepts 1 or 2, cancelled.\n");
	}
	return out;
}
/**
   Set elements to zero.

   @param[in,out] A The matrix. Skipped if NULL.
*/
void X(zero)(X(mat)* A){
	if(A&&P(A)){
		memset(P(A), 0, sizeof(T)*PN(A));
	}
}
/**
   Set column elements to zero.

   @param[in,out] A The matrix.
   @param[in] icol Index of the column.
   @return 0 on success, -1 if A is invalid or icol is out of range.
*/
int X(zerocol)(X(mat)* A, int icol){
	int ans=0;
	if(check_mat(A)){
		if(icol>=0&&icol<A->ny){
			memset(PCOL(A, icol), 0, sizeof(T)*A->nx);
		} else{
			warning("Invalid range.\n");
			ans=-1;
		}
	}else{
		ans=-1;
	}
	return ans;
}
/**
   Compute the hashlittle of the matrix content using key as the seed. Returns
   the updated key. Returns key unchanged if A is invalid.

   @param[in] A The matrix.
   @param[in] key Initial hash key.
   @return Updated hash key.
*/
uint32_t X(hash)(const X(mat)* A, uint32_t key){
	if(check_mat(A)){
		key=hashlittle(P(A), PN(A)*sizeof(T), key);
	}
	return key;
}
/**
   copy the values from one matrix to another.

   @param[in,out] out0 Output matrix. Allocated if NULL; otherwise its dimension must match in. Zeroed if in is empty.
   @param[in] in Input matrix.
*/
void X(cp)(X(mat)** out0, const X(mat)* in){
	if(check_mat(in)){
		if(X(init)(out0, in->nx, in->ny)) return;
		X(mat)* out=*out0;
		if(in->keywords){
			free(out->keywords);
			out->keywords=strdup(in->keywords);
		}
		if(P(out)!=P(in)){
			memcpy(P(out), P(in), in->nx*in->ny*sizeof(T));
		}
	} else{
		X(zero)(*out0);
	}
}

/**
   duplicate a matrix including keywords. The caller owns the result.

   @param[in] in The matrix to copy.
   @return Newly allocated matrix; the caller frees it. NULL if in is invalid.
*/
X(mat)* X(dup)(const X(mat)* in){
	if(!check_mat(in)) return NULL;
	X(mat)* out=NULL;
	X(cp)(&out, in);
	return out;
}

/**
   transpose a matrix into a newly allocated one.

   @param[in] A The matrix.
   @return Newly allocated matrix; the caller frees it. NULL if A is invalid.
*/
X(mat)* X(trans)(const X(mat)* A){
	if(!check_mat(A)) return NULL;
	X(mat)* B=X(new)(A->ny, A->nx);
	if(A->nx==1||A->ny==1){
		memcpy(P(B), P(A), PN(A)*sizeof(T));
	} else{
		for(int iy=0; iy<A->ny; iy++){
			for(int ix=0; ix<A->nx; ix++){
				P(B, iy, ix)=P(A, ix, iy);
			}
		}
	}
	return B;
}
/**
   set values of each element in an matrix to val.

   @param[in,out] A The matrix.
   @param[in] val The value.
*/
void X(set)(X(mat)* A, const T val){
	if(!check_mat(A)) return;
	for(long i=0; i<PN(A); i++){
		P(A, i)=val;
	}
}

/**
   display the content of a matrix, labeled by the printf-style format, with
   info(). Output is truncated to 20 rows and 12 columns (100 elements for a vector).

   @param[in] A The matrix.
   @param[in] format printf style format of the label.
*/
void X(show)(const X(mat)* A, const char* format, ...){
	format2fn;
	if(!check_mat(A)){
		warning("Invalid: %s is not matrix.", fn);
		return;
	}
	if(PN(A)>20 || (NX(A)>1 && NY(A)>1)){
		info("Displaying content of %s (%ldx%ld):\n", fn, NX(A), NY(A));
	}else{
		info("%s:", fn);//display inline
	}
	int colmax=12;
	int rowmax=20;
	//int iset, i, j;
	int nset=1;//(A->ny+colmax-1)/colmax;
	for(int iset=0; iset<nset; iset++){
		int ncol=(iset+1)*colmax;
		if(A->ny==1){//Display horizontaly
			if(ncol>A->nx) ncol=(int)A->nx;
			if(ncol>100) ncol=100;
			for(int i=0; i<ncol; i++){
				PRINT(P(A,i));
				if ((i+1)%colmax==0) info("\n");
			}
			if(ncol%colmax!=0) info("\n");
		}else{
			if(ncol>A->ny) ncol=(int)A->ny;
			//info("Cols %d to %d\n", iset, ncol-1);
			int nxmax=(int)MIN(A->nx, rowmax);
			for(int j=0; j<nxmax; j++){
				for(int i=iset*colmax; i<ncol; i++){
					PRINT(P(A, j, i));
				}
				info("\n");
			}
		}
	}
}
/**
   display diagonal of a matrix. 

   @param[in] A The matrix.
   @param[in] format printf style format of the label.
*/
void X(show_diag)(const X(mat)* A, const char* format, ...){
	format2fn;
	if(!check_mat(A)){
		warning("Invalid: %s is not matrix.", fn);
		return;
	}
	info("%s diagonal:", fn);//display inline
	long n=MIN(NX(A), NY(A));
	long ns=MIN(n, 6);
	for(long i=0; i<ns; i++){
		PRINT(P(A,i,i));
	}
	if(n>ns*2){
		info("...");
		ns=n-ns;
	}
	for(long i=ns; i<n; i++){
		PRINT(P(A,i,i));
	}
	info("\n");
}
/**
   Permute the vector so that
   out(:)=in(p);

   @param[out] out Output vector with the same dimension as in.
   @param[in] in Input vector.
   @param[in] perm Permutation index. A negative (or zero) entry -j takes the conjugate of in(j).
*/
void X(vecperm)(X(mat)* out, const X(mat)* in, const long* perm){
	if(!check_mat(in)||!check_mat(out)||!check_match(in, out)) return;
	for(long i=0; i<in->nx*in->ny; i++){
		if(perm[i]>0){
			P(out, i)=P(in, perm[i]);
		} else{
			P(out, i)=CONJ(P(in, -perm[i]));
		}
	}
}
/**
   Reverse permute the vector so that
   out(p)=in(:);

   @param[out] out Output vector with the same dimension as in.
   @param[in] in Input vector.
   @param[in] perm Permutation index. A negative (or zero) entry -j stores the conjugate into out(j).
*/
void X(vecpermi)(X(mat)* out, const X(mat)* in, const long* perm){
	if(!check_mat(in)||!check_mat(out)||!check_match(in, out)) return;
	assert(in->nx*in->ny==out->nx*out->ny);
	for(long i=0; i<in->nx*in->ny; i++){
		if(perm[i]>0){
			P(out, perm[i])=P(in, i);
		} else{
			P(out, -perm[i])=CONJ(P(in, i));
		}
	}
}
/**
   Flip the matrix along the set axis. 0 to flip both, 1 along x, 2 along y.

   @param[in,out] A The matrix flipped in place.
   @param[in] axis 0 to flip both, 1 along x, 2 along y.
   @return 0 on success, -1 if axis is invalid.
*/
int X(flip)(X(mat)* A, int axis){
	if(!A || !PN(A)) return 0;
	if(axis<0 || axis>2) {
		dbg("flip: A is null or axis is invalid: %p %d\n", A, axis);
		return -1;
	}
	if(axis==0){//flipping both axes is equivalent to reversing the linear index
		const long nn=PN(A);
		for(long i=0; i<nn/2; i++){
			T tmp=P(A, nn-1-i);
			P(A, nn-1-i)=P(A, i);
			P(A, i)=tmp;
		}
		return 0;
	}
	const long nx=A->nx;
	const long ny=A->ny;
	const long nx2=axis==1?nx/2:nx;//only loop over half of the flipped axis
	const long ny2=axis==2?ny/2:ny;
	for(long iy=0; iy<ny2; iy++){
		const long jy=axis==2?(ny-1-iy):iy;
		for(long ix=0; ix<nx2; ix++){
			const long jx=axis==1?(nx-1-ix):ix;
			T tmp=P(A, jx, jy);
			P(A, jx, jy)=P(A, ix, iy);
			P(A, ix, iy)=tmp;
		}
	}
	return 0;
}
/**
   Calculate sum of all the elements in A.

   @param[in] A The matrix.
   @return The sum; 0 if A is invalid.
*/
T X(sum)(const X(mat)* A){
	T v=0;
	if(check_mat(A)){
		v=X(vecsum)(P(A), PN(A));
	}
	return v;
}
/**
   Calculate average of all the elements in A.

   @param[in] A The matrix.
   @return The average; 0 if A is NULL.
*/
T X(mean)(const X(mat) *A){
	if(!A) return 0;
	return X(sum)(A)/(R)PN(A);
}
/**
   compute the trace (sum of diagonal elements)

   @param[in] A The matrix.
   @return The trace; 0 if A is invalid.
*/
T X(trace)(const X(mat)* A){
	TD trace=0;
	if(check_mat(A)){
		long n=MIN(A->nx, A->ny);
		for(long i=0; i<n; i++){
			trace+=P(A, i, i);
		}
	}
	return (T)trace;
}

#ifndef COMP_COMPLEX
/**
   qsort comparison function for ascending order.

   @param[in] pA The first value.
   @param[in] pB The second value.
   @return 1, 0, or -1 if *A is greater than, equal to, or less than *B.
*/
static int sort_ascend(const void* pA, const void* pB){
	const T A=*(const T*)pA;
	const T B=*(const T*)pB;
	return (A>B)-(A<B);
}
/**
   qsort comparison function for descending order.

   @param[in] pA The first value.
   @param[in] pB The second value.
   @return -1, 0, or 1 if *A is greater than, equal to, or less than *B.
*/
static int sort_descend(const void* pA, const void* pB){
	const T A=*(const T*)pA;
	const T B=*(const T*)pB;
	return (A<B)-(A>B);
}
/**
  Sort all columns of A, in ascending order if ascend is non zero, otherwise in descending order.

  @param[in,out] A The matrix sorted in place.
  @param[in] ascend Non zero for ascending order, otherwise descending.
*/
void X(sort)(X(mat)* A, int ascend){
	if(!check_mat(A)) return;
	for(int i=0; i<A->ny; i++){
		qsort(PCOL(A, i), A->nx, sizeof(R), ascend?sort_ascend:sort_descend);
	}
}

/**
   Find unique values and sort them.

   @param[in] input Input array.
   @param[in] ascend Non zero for ascending order, otherwise descending.
   @return Newly allocated vector of unique values; the caller frees it.
*/
X(mat) *X(unique)(const X(mat) *input, int ascend){
	X(mat) *output=X(new)(PN(input), 1);//unique dtrats
	int nout=0;
	for(int ix=0; ix<PN(input); ix++){
		int found=0;
		for(int jx=0; jx<nout; jx++){
			if(P(output,jx)==P(input,ix)){
				found=1; break;
			}
		}
		if(!found){
			P(output,nout)=P(input,ix);
			nout++;
		}
	}
	X(resize)(output, nout, 1);
	X(sort)(output, ascend);//from fast to slow
	return output;
}
#endif
/**
   find the maximum and minimum value

   @param[in] A The matrix.
   @param[out] max The maximum value. May be NULL.
   @param[out] min The minimum value. May be NULL.
*/
void X(maxmin)(const X(mat)*A, R*max, R*min){
	if(!check_mat(A)) return ;
	X(vecmaxmin)(P(A), PN(A), max, min);
}
/**
   find the maximum value

   @param[in] A The matrix.
   @return The maximum value; 0 if A is invalid.
*/
R X(max)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	R max, min;
	X(vecmaxmin)(P(A), PN(A), &max, &min);
	return max;
}

/**
   find the minimum value

   @param[in] A The matrix.
   @return The minimum value; 0 if A is invalid.
*/
R X(min)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	R max, min;
	X(vecmaxmin)(P(A), PN(A), &max, &min);
	return min;
}
/**
   find the maximum of abs of all elements

   @param[in] A The matrix.
   @return The maximum absolute value; 0 if A is invalid.
*/
R X(maxabs)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	return X(vecmaxabs)(P(A), PN(A));
}
/**
   compute the sum of abs of all elements

   @param[in] A The matrix.
   @return The sum of absolute values; 0 if A is invalid.
*/
R X(sumabs)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	R out=0;
	for(long i=0; i<PN(A); i++){
		out+=ABS(P(A, i));
	}
	return out;
}
/**
   compute the sum of A.*A

   @param[in] A The matrix.
   @return The sum of squared absolute values; 0 if A is invalid.
*/
R X(sumsq)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	R out=0;
	for(long i=0; i<PN(A); i++){
		out+=ABS2(P(A, i));
	}
	return out;
}
/**
   compute the sum of (A-B)^2

   @param[in] A The first matrix.
   @param[in] B The second matrix with the same dimension.
   @return The sum of squared absolute differences; -1 if the inputs are invalid.
*/
R X(sumdiffsq)(const X(mat)* A, const X(mat)* B){
	if(!check_mat(A)||!check_mat(B)||!check_match(A, B)) return -1;
	RD out=0;
	for(long i=0; i<PN(A); i++){
		out+=ABS2(P(A, i)-P(B, i));
	}
	return (R)out;
}


/**
   shift frequency components by n/2

   @param[in,out] A The matrix shifted in place. Dimensions should be even.
*/
void X(fftshift)(X(mat)* A){
	if(!check_mat(A)) return;
	long i;
	const long nx=A->nx;
	const long ny=A->ny;
	if((nx&1)==1){
		warning_once("nx=%ld is not multiple of 2\n", nx);
	}
	const long nx2=nx/2;
	const long ny2=ny/2;
	const long nx2d=nx2*sizeof(T);
	T* tmp=(T*)malloc(nx2d);
	T* data=P(A);
	if(ny==1){
		memcpy(tmp, data, nx2d);
		memcpy(data, data+nx2, nx2d);
		memcpy(data+nx2, tmp, nx2d);
	} else{
		assert((ny&1)==0);
		for(i=0; i<ny2; i++){
			memcpy(tmp, data+i*nx, nx2d);
			memcpy(data+i*nx, data+(i+ny2)*nx+nx2, nx2d);
			memcpy(data+(i+ny2)*nx+nx2, tmp, nx2d);
			memcpy(tmp, data+i*nx+nx2, nx2d);
			memcpy(data+i*nx+nx2, data+(i+ny2)*nx, nx2d);
			memcpy(data+(i+ny2)*nx, tmp, nx2d);
		}
	}

	free(tmp);
}

/**
   reorder B and embed/crop into center of A .

   \verbatim
   4 * * 3
   * * * *
   * * * *
   2 * * 1
   \endverbatim
   to
   \verbatim
   1 2
   3 4
   \endverbatim

   @param[out] A Output matrix. Zeroed if larger than B.
   @param[in] B Input matrix. Must not share memory with A.
   @return 0 on success, -1 if the inputs are invalid or A and B share memory.
*/
int X(cpcorner2center)(X(mat)* A, const X(mat)* B){
	if(!check_mat(A)||!check_mat(B) || P(A)==P(B)) return -1;
	const long Anx=A->nx;
	const long Any=A->ny;
	const long Bnx=B->nx;
	const long Bny=B->ny;
	if(Anx>Bnx||Any>Bny){
		X(zero)(A);
	}
	if(!((Anx&1)==0&&(Any&1)==0&&(Bnx&1)==0&&(Bny&1)==0)){
		warning("Arrays with odd index may not work as expected\n");
	}
	const long ny2=MIN(Any, Bny)>>1;
	const long nx2=MIN(Anx, Bnx)>>1;
	const long xoff=Anx/2-nx2;
	const long yoff=Any/2-ny2;
	const size_t nx2d=nx2*sizeof(T);
	for(long i=0; i<ny2; i++){
		memcpy(&P(A, xoff, 		yoff+i), 	 &P(B, Bnx-nx2,	Bny-ny2+i), nx2d);
		memcpy(&P(A, xoff+nx2, 	yoff+i),	 &P(B, 0, 	  	Bny-ny2+i), nx2d);
		memcpy(&P(A, xoff, 		yoff+i+ny2), &P(B, Bnx-nx2,	i), 		nx2d);
		memcpy(&P(A, xoff+nx2, 	yoff+i+ny2), &P(B, 0, 		i), 		nx2d);
	}
	return 0;
}

/**
   cyclic shift A by nx and ny to B.
   \verbatim
   4   3     1   2

   2   1 to  3   4
   \endverbatim

   @param[in,out] B0 Output matrix. Allocated if NULL; otherwise its dimension must match A.
   @param[in] A Input matrix. Must not share memory with *B0.
   @param[in] sx Shift along x (rows).
   @param[in] sy Shift along y (columns).
   @return 0 on success, -1 if A is invalid or shares memory with *B0.
*/
int X(shift)(X(mat)** B0, const X(mat)* A, int sx, int sy){
	if(!check_mat(A)) return -1;
	if(X(init)(B0, A->nx, A->ny)) return -1;
	X(mat)* B=*B0;
	if(P(B)==P(A)) return -1;
	const long nx=A->nx;
	const long ny=A->ny;
	sx=(int)(sx%nx); if(sx<0) sx+=nx;
	sy=(int)(sy%ny); if(sy<0) sy+=ny;
	if(sx!=0||sy!=0){
		long dy=ny-sy;
		long dx=nx-sx;
		for(long iy=0; iy<sy; iy++){
			memcpy(&P(B, 0, iy), &P(A, dx, iy+dy), sx*sizeof(T));/*3 */
			memcpy(&P(B, sx, iy), &P(A, 0, iy+dy), dx*sizeof(T));/*4 */
		}
		for(long iy=sy; iy<ny; iy++){
			memcpy(&P(B, 0, iy), &P(A, dx, iy-sy), sx*sizeof(T));/*1 */
			memcpy(&P(B, sx, iy), &P(A, 0, iy-sy), dx*sizeof(T));/*2 */
		}
	} else{
		X(cp)(B0, A);
	}
	return 0;
}
/**
   cast a cell object to actual cell after checking.

   @param[in] A The cell object.
   @return A cast to X(cell) if all non-empty elements are matrices, otherwise NULL. Not a new reference.
*/
X(cell)* X(cell_cast)(const cell* A){
	if(!iscell(A)) return NULL;
	for(long iy=0; iy<NY(A); iy++){
		for(long ix=0; ix<NX(A); ix++){
			if(P(A, ix, iy)&&!ismat(P(A, ix, iy))){
				dbg("A(%ld,%ld) is not mat, return NULL.\n", ix, iy);
				return NULL;
			}
		}
	}
	return (X(cell)*)A;
}
/**
   create an new cell similar to A in shape.
   When a cell is empty, it is created with an emtry matrix and cannot be overriden.

   @param[in] A The cell to copy the shape from.
   @return Newly allocated zero cell with continuous memory; the caller frees it. NULL if A is NULL.
*/
X(cell)* X(cellnew2)(const X(cell)* A){
	if(!A) return NULL;
	X(cell)* out=X(cellnew)(A->nx, A->ny);
	long tot=0;
	for(long i=0; i<PN(A); i++){
		if(!isempty(P(A, i))){
			tot+=P(A, i)->nx*P(A, i)->ny;
		}
	}
	out->m=X(new)(tot, 1);
	tot=0;
	for(int i=0; i<PN(A); i++){
		if(!isempty(P(A, i))){
			P(out, i)=X(new_do)(P(A, i)->nx, P(A, i)->ny, P(out->m)+tot, out->m->mem);
			tot+=P(A, i)->nx*P(A, i)->ny;
		} else{
			P(out, i)=X(new)(0, 0);//place holder to avoid been overriden.
		}
	}
	return out;
}

/**
   Create an new cell with matrix specified. Each block is stored continuously in memory.

   @param[in] nx Number of rows of the cell.
   @param[in] ny Number of columns of the cell.
   @param[in] nnx Number of rows of each block, nx*ny entries. A negative pointer value -m means all blocks have m rows.
   @param[in] nny Number of columns of each block, nx*ny entries. A negative pointer value -m means all blocks have m columns. NULL means 1.
   @return Newly allocated zero cell; the caller frees it.
*/
X(cell)* X(cellnew3)(long nx, long ny, long* nnx, long* nny){
	long tot=0;
	for(long i=0; i<nx*ny; i++){
		long mx=(long)nnx>0?nnx[i]:(-(long)nnx);
		long my=nny?((long)nny>0?nny[i]:(-(long)nny)):1;
		tot+=mx*my;
	}
	X(cell)* out=X(cellnew)(nx, ny);
	
	out->m=X(new)(tot, 1);
	tot=0;
	for(long i=0; i<nx*ny; i++){
		long mx=(long)nnx>0?nnx[i]:(-(long)nnx);
		long my=nny?((long)nny>0?nny[i]:(-(long)nny)):1;
		P(out, i)=X(new_do)(mx, my, P(out->m)+tot, out->m->mem);
		tot+=mx*my;
	}
	
	return out;
}
/**
   Create an new cell with all blocks of the same size. Each block is stored continuously in memory.

   @param[in] nx Number of rows of the cell.
   @param[in] ny Number of columns of the cell.
   @param[in] mx Number of rows of each block.
   @param[in] my Number of columns of each block.
   @return Newly allocated zero cell; the caller frees it.
*/
X(cell)* X(cellnew_same)(long nx, long ny, long mx, long my){
	return X(cellnew3)(nx, ny, (long*)-mx, (long*)-my);
}
/**
   Calls cellnew3 with a filename to be saved to.

   @param[in] nx Number of rows of the cell.
   @param[in] ny Number of columns of the cell.
   @param[in] nnx Number of rows of each block. See X(cellnew3).
   @param[in] nny Number of columns of each block. See X(cellnew3).
   @param[in] keywords Keywords to save with the data. May be NULL.
   @param[in] format printf style format of the file name. A leading - disables saving.
   @return Newly allocated cell; the caller frees it. NULL if nx or ny is 0.
*/
X(cell)* X(cellnew_file)(long nx, long ny, long* nnx, long* nny, const char* keywords, const char* format, ...){
	if(!nx||!ny) return NULL;
	format2fn;
	if(fn&&fn[0]=='-') fn=NULL;//leading - disables filename.
	X(cell)* out=X(cellnew3)(nx, ny, nnx, nny);
	if(out&&keywords) out->keywords=strdup(keywords);
	if(disable_save&&!IS_SHM(fn))fn=NULL;
	if(out && fn) {
		out->fp=zfopen(fn, "w");
		writedata(out->fp, out, -1);
	}
	return out;
}
/**
   Calls cellnew_same with a filename to be saved to.

   @param[in] nx Number of rows of the cell.
   @param[in] ny Number of columns of the cell.
   @param[in] mx Number of rows of each block.
   @param[in] my Number of columns of each block.
   @param[in] keywords Keywords to save with the data. May be NULL.
   @param[in] format printf style format of the file name. A leading - disables saving.
   @return Newly allocated cell; the caller frees it. NULL if nx or ny is 0.
*/
X(cell)* X(cellnewsame_file)(long nx, long ny, long mx, long my, const char* keywords, const char* format, ...){
	if(!nx||!ny) return NULL;
	format2fn;
	if(fn&&fn[0]=='-') fn=NULL;//leading - disables filename.
	X(cell)* out=X(cellnew_same)(nx, ny, mx, my);
	if(out&&keywords) out->keywords=strdup(keywords);
	if(disable_save&&!IS_SHM(fn))fn=NULL;
	if(out && fn) {
		out->fp=zfopen(fn, "w");
		writedata(out->fp, out, -1);
	}
	return out;
}
/**
   creat a cell reference an existing cell by referencing the
   elements.

   @param[in] in The cell to reference.
   @return Newly allocated cell sharing data with in; the caller frees it. NULL if in is NULL.
*/
X(cell)* X(cellref)(const X(cell)* in){
	if(!in) return NULL;
	X(cell)* out=X(cellnew)(in->nx, in->ny);
	if(in->m){
		out->m=X(ref)(in->m);
	}
	for(int i=0; i<in->nx*in->ny; i++){
		P(out, i)=X(ref)(P(in, i));
	}
	if(in->keywords) out->keywords=strdup(in->keywords);
	return out;
}

/**
   setting all elements of a cell to val.

   @param[in,out] dc The cell.
   @param[in] val The value.
*/
void X(cellset)(X(cell)* dc, T val){
	if(dc){
		for(int ix=0; ix<dc->nx*dc->ny; ix++){
			if(P(dc, ix)) X(set)(P(dc, ix), val);
		}
	}
}

/**
   transpose a cell object, including each element.

   @param[in] A The cell.
   @return Newly allocated cell; the caller frees it. NULL if A is NULL.
*/
X(cell)* X(celltrans)(const X(cell)* A){
	if(!A) return NULL;
	X(cell)* B=X(cellnew)(A->ny, A->nx);
	for(int ix=0; ix<A->nx; ix++){
		for(int iy=0; iy<A->ny; iy++){
			P(B, iy, ix)=X(trans)(P(A, ix, iy));
		}
	}
	return B;
}

/**
   reduce nx*ny cell matrix to 1*ny if dim=1 and nx*1 if dim=2 by merging the cells.

   @param[in] A The cell.
   @param[in] dim 1 to merge along rows, 2 to merge along columns.
   @return Newly allocated cell; the caller frees it. NULL if A is empty.
*/
X(cell)* X(cellreduce)(const X(cell)* A, int dim){
	if(!A) return NULL;
	X(cell)* out=NULL;
	long nx, ny, * nxs, * nys;
	celldim(A, &nx, &ny, &nxs, &nys);
	if(nx==0||ny==0){
		free(nxs);
		free(nys);
		return NULL;
	}
	if(dim==1){
		out=X(cellnew)(1, A->ny);
		for(long iy=0; iy<A->ny; iy++){
			if(nys[iy]==0) continue;
			P(out, iy)=X(new)(nx, nys[iy]);
			for(long icol=0; icol<nys[iy]; icol++){
				long kr=0;
				for(long ix=0; ix<A->nx; ix++){
					if(!isempty(P(A, ix, iy))){
						memcpy(&P(P(out, iy), kr, icol), PCOL(P(A, ix, iy), icol),
							nxs[ix]*sizeof(T));
					}
					kr+=nxs[ix];
				}
			}
		}
	} else if(dim==2){
		out=X(cellnew)(A->nx, 1);
		for(long ix=0; ix<A->nx; ix++){
			if(nxs[ix]==0) continue;
			P(out, ix)=X(new)(nxs[ix], ny);
			long kr=0;
			for(long iy=0; iy<A->ny; iy++){
				if(!isempty(P(A, ix, iy))){
					memcpy(PCOL(P(out, ix), kr), P(P(A, ix, iy)), nxs[ix]*nys[iy]*sizeof(T));
				}
				kr+=nys[iy];
			}
		}
	} else{
		error("Invalid dim=%d\n", dim);
	}
	free(nxs);
	free(nys);
	return out;
}

/**
   drop empty rows or columns. size of *A0 is changed.

   @param[in,out] A0 The cell. Freed and set to NULL if all are dropped.
   @param[in] dim 1 to drop empty rows, 2 to drop empty columns.
*/
void X(celldropempty)(X(cell)** A0, int dim){
	if(!A0) return;
	X(cell)* A=*A0;
	if(!A) return;
	if(dim==1){
		/*drop rows */
		int keep[A->nx];
		int ndrop=0;
		for(int ix=0; ix<A->nx; ix++){
			keep[ix]=0;
			for(int iy=0; iy<A->ny; iy++){
				if(!isempty(P(A, ix, iy))){
					keep[ix]=1;
					break;
				}
			}
			if(keep[ix]==0){
				ndrop++;
			}
		}
		if(ndrop){
			if(ndrop==A->nx){
				X(cellfree)(A);
				*A0=NULL;
			} else{
				X(cell)* B=X(cellnew)(A->nx-ndrop, A->ny);
				int count=0;
				for(int ix=0; ix<A->nx; ix++){
					if(keep[ix]){
						for(int iy=0; iy<A->ny; iy++){
							P(B, count, iy)=P(A, ix, iy);
							P(A, ix, iy)=NULL;
						}
						count++;
					} else{
						warning("row %d dropped\n", ix);
					}
				}
				B->m=A->m; A->m=NULL;//transfer the contiguous storage
				X(cellfree)(A);//frees the dropped (empty) blocks and the container
				*A0=B;
			}
		}
	} else if(dim==2){
	   /*drop cols */
		int count=0;
		for(int iy=0; iy<A->ny; iy++){
			int keep=0;
			for(int ix=0; ix<A->nx; ix++){
				if(!isempty(P(A, ix, iy))){
					keep=1;
					break;
				}
			}
			if(keep){
				for(int ix=0; ix<A->nx; ix++){
					P(A, ix, count)=P(A, ix, iy);
				}
				count++;
			} else{
				for(int ix=0; ix<A->nx; ix++){//free the empty blocks
					X(free)(P(A, ix, iy));
				}
			}
		}
		A->ny=count;
		if(count==0){
			X(cellfree)(A);
			*A0=NULL;
		} else{
			P(A)=myrealloc(P(A), A->ny*A->nx, X(mat)*);
		}
	} else{
		error("Invalid dim: %d\n", dim);
	}

}


/**
   convert a vector to cell using dimensions specified in dims. Reference the vector

   @param[in] A The vector.
   @param[in] dims Number of elements of each block. Must sum to the length of A.
   @param[in] ndim Number of blocks.
   @return Newly allocated ndim*1 cell sharing data with A; the caller frees it.
*/
X(cell)* X(2cellref)(const X(mat)* A, long* dims, long ndim){
	if(!A||!dims) return NULL;
	long nx=0;
	for(int ix=0; ix<ndim; ix++){
		nx+=dims[ix];
	}
	if(A->ny!=1){
		error("Use d2cell2 instead for non vectors\n");
		return NULL;
	}
	if(nx!=A->nx){
		error("Shape doesn't agree. nx=%ld, nx=%ld\n", nx, A->nx);
		return NULL;
	}
	long kr=0;
	X(cell)* B=X(cellnew)(ndim, 1);
	B->m=X(ref)(A);
	for(long ix=0; ix<ndim; ix++){
		P(B, ix)=X(new_do)(dims[ix], 1, P(A)+kr, A->mem);
		kr+=dims[ix];
	}
	return B;
}

/**
   make A a cell array using shape information from ref if *B is NULL.

   Notice that empty cells may be replaced by zeros if it is within valid column or row.

   @param[in,out] B The output cell. Allocated using shape of ref if NULL; otherwise used as reference.
   @param[in] A The matrix with the same total dimension as the reference.
   @param[in] ref The reference cell for the shape. Used only if *B is NULL.
*/
void X(2cell)(X(cell)** B, const X(mat)* A, const X(cell)* ref){
	long nx, ny, * nxs, * nys;
	if(*B) ref=*B;/*use B as reference. */
	if(!ref) return;
	celldim(ref, &nx, &ny, &nxs, &nys);
	if(!nxs||!nys){//empty or invalid reference
		if(A&&PN(A)){
			error("Shape doesn't agree. Reference is empty or invalid but input is %ldx%ld\n", A->nx, A->ny);
		}else if(!*B && !PN(ref)){
			*B=X(cellnew)(ref->nx, ref->ny);
		}
		free(nxs); free(nys);
		return;
	}
	if(!A){
		if(nx!=0||ny!=0){
			error("Shape doesn't agree. A is 0x0. Reference is %ldx%ld\n",
				nx, ny);
			free(nxs); free(nys);
			return;
		}
	} else if(nx!=A->nx||ny!=A->ny){
		error("Shape doesn't agree. Reference is %ldx%ld but input is %ldx%ld\n",
			nx, ny, A->nx, A->ny);
		free(nxs); free(nys);
		return;
	}
	if(!*B){
		*B=X(cellnew)(ref->nx, ref->ny);
	}
	X(cell)* Bp=(*B);
	for(long iy=0; iy<ref->ny; iy++){
		for(long ix=0; ix<ref->nx; ix++){
			if(nxs[ix]&&nys[iy]&&!P(Bp, ix, iy)){
				P(Bp, ix, iy)=X(new)(nxs[ix], nys[iy]);
			}
		}
	}
	
	long jcol=0;
	for(long iy=0; iy<ref->ny; iy++){
		for(long icol=0; icol<nys[iy]; icol++){
			long kr=0;
			for(long ix=0; ix<ref->nx; ix++){
				if(nxs[ix]){
					memcpy(PCOL(P(Bp, ix, iy), icol), &P(A, kr, icol+jcol), nxs[ix]*sizeof(T));
					kr+=nxs[ix];
				}
			}
		}
		jcol+=nys[iy];
	}
	free(nxs);
	free(nys);
}



/**
   Create a new sub cell matrix of nx*ny starting from(sx,sy)

   @param[in] in The input cell.
   @param[in] sx Starting row.
   @param[in] nx Number of rows. <=0 to use all remaining rows.
   @param[in] sy Starting column.
   @param[in] ny Number of columns. <=0 to use all remaining columns.
   @return Newly allocated cell with elements referencing those of in; the caller frees it. NULL if in is NULL.
*/
X(cell)* X(cellsub)(const X(cell)* in, long sx, long nx, long sy, long ny){
	if(!in) return NULL;
	if(nx<=0){
		nx=in->nx-sx;
	}
	if(ny<=0){
		ny=in->ny-sy;
	}
	if(sx<0||sy<0||sx+nx>in->nx||sy+ny>in->ny){
		error("Invalid parameter range\n");
		return NULL;
	}
	X(cell)* out=X(cellnew)(nx, ny);
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(out, ix, iy)=X(ref)(P(in, ix+sx, iy+sy));
		}
	}
	return out;
}

/**
   input is nsa*ncol cell. each cell has npix=nx*ny elements. Extract icol of cell as npix*nsa array.

   @param[in] input The cell.
   @param[in] icol Index of the column.
   @return Newly allocated matrix; the caller frees it. If input->m is set, it is a weak reference into input->m and must not outlive input. NULL if input is NULL.
*/
X(mat)* X(cell_col)(X(cell)* input, long icol){
	if(!input) return NULL;
	if(icol<0||icol>=NY(input)){
		error("Column %ld exceeds number of columns %ld\n", icol, NY(input));
		return NULL;
	}
	long npix=PN(input, 0, 0);
	long nsa=input->nx;
	if(input->m){
		if(PN(input->m)!=npix*PN(input)){
			error("dimension mismatch.\n");
			return NULL;
		}
		return X(new_do)(npix, nsa, P(input->m)+(npix*nsa)*icol, NULL);
	} else{
		dbg("Consider using m to speed up\n");
		X(mat)* output=X(new)(npix, nsa);
		for(long isa=0; isa<nsa; isa++){
			if(PN(input, isa, icol)!=npix){
				error("dimension mismatch.\n");
				X(free)(output);
				return NULL;
			}
			memcpy(PCOL(output, isa), P(P(input, isa, icol)), sizeof(T)*npix);
		}
		return output;
	}
}
/**
   create sum of all the elements in A.

   @param[in] A The cell.
   @return The sum of all elements of all matrices; 0 if A is empty.
*/
T X(cellsum)(const X(cell)* A){
	if(isempty(A)) return 0;
	T v=0;
	for(long i=0; i<PN(A); i++){
		v+=X(sum)(P(A, i));
	}
	return v;
}
/**
   return sum of each cell in a X(mat)

   @param[in] A The cell.
   @return Newly allocated matrix of the same shape as A with the sum of each element; the caller frees it. NULL if A is empty.
*/
X(mat)* X(cellsum_each)(const X(cell)* A){
	if(isempty(A)) return 0;
	X(mat) *B=X(new)(NX(A), NY(A));
	for(long i=0; i<PN(A); i++){
		P(B,i)+=X(sum)(P(A, i));
	}
	return B;
}
