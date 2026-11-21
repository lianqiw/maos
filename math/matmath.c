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
#include "../sys/sys.h"
#include "random.h"
#include "type.h"
#include "mathdef.h"
#include "defs.h"

/*
	Math routines that are not applicable to integers and therefore not included in mat.c
*/
/**
   Scale each element of A by w.

   @param[in,out] A The matrix to scale in place.
   @param[in] w The scaling factor.
*/
void X(scale)(X(mat)* A, R w){
	if(!A || !check_mat(A)) return;
	if(w==(T)0){
		if(P(A)) memset(P(A), 0, sizeof(T)*A->nx*A->ny);
	} else{
		OMP_SIMD()
		for(long i=0; i<A->nx*A->ny; i++){
			P(A, i)*=w;
		}
	}
}
/**
   Normalize the sum of A to sum.

   @param[in,out] A The matrix to normalize in place.
   @param[in] sum The target sum.
*/
void X(normalize_sum)(X(mat) *A, T sum){
	if(A) X(scale)(A, REAL(sum/X(sum)(A)));
}
/**
   Normalize the sum of abs(A) to sum.

   @param[in,out] A The matrix to normalize in place.
   @param[in] sum The target sum of absolute values.
*/
void X(normalize_sumabs)(X(mat) *A, R sum){
	if(A) X(scale)(A, sum/X(sumabs)(A));
}
/**
   Check for NaN in elements.

   @param[in] A The matrix.
   @return 1 if any element (real part) is NaN, 0 otherwise or if A is invalid.
*/
int X(isnan)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	for(long i=0; i<A->nx*A->ny; i++){
		if(isnan(REAL(P(A, i)))){
			return 1;
		}
	}
	return 0;
}
/**
   Compute the 2-norm of A divided by sqrt of number of columns, i.e., sqrt(sum(abs(A).^2)/ny).

   @param[in] A The matrix.
   @return sqrt(sum(abs(A).^2)/ny), or 0 if A is invalid.
*/
R X(norm)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	return sqrt(X(sumsq)(A)/A->ny);
}
/**
   Compute the (sample) standard deviation of all elements of A.

   @param[in] A The matrix.
   @return The standard deviation, or 0 if A is invalid.
*/
R X(std)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	long N=A->nx*A->ny;
	T sum=X(sum)(A);
	R var=(X(sumsq)(A)-ABS2(sum)/N)/(N-1);
	return sqrt(var);
}
/**
   Fill A with random uniform numbers between [0, 1)*max.

   @param[out] A The matrix to fill.
   @param[in] max The upper bound.
   @param[in,out] rstat The random stream.
*/
void X(randu)(X(mat)* A, const T max, rand_t* rstat){
	if(!check_mat(A)){
		warning("Input is not valid");
		return;
	}
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=RANDU(rstat)*max;
	}
}
/**
   Fill A with random normal distribution numbers with
   standard deviation of sigma.

   @param[out] A The matrix to fill.
   @param[in] sigma The standard deviation.
   @param[in,out] rstat The random stream.
*/
void X(randn)(X(mat)* A, const T sigma, rand_t* rstat){
	if(!check_mat(A)){
		warning("Input is not valid");
		return;
	}
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=RANDN(rstat)*sigma;
	}
}

/**
   Compute the dot product (inner product) of A and B without conjugation.

   @param[in] A The first matrix.
   @param[in] B The second matrix. Must have the same size as A.
   @return sum(A.*B), or 0 if A or B is NULL or empty, or if they do not match (error is raised).
*/
T X(dot)(const X(mat)* A, const X(mat)* B){
	if(!A||!B||!PN(A)||!PN(B)||!check_match(A, B)) return 0;
	TD out=0;
	OMP_SIMD_R(reduction(+:out))
	for(long i=0; i<A->nx*A->ny; i++){
		out+=P(A, i)*P(B, i);
	}
	if(isnan(REAL(out))){
		X(show)(A, "inn_A");
		X(show)(B, "inn_B");
		error("X(dot): NaN found\n");
	}
	return (T)out;
}

/**
   Compute weighted dot product a'*(w*b).

   @param[in] a The left vector of length w->nx.
   @param[in] w The weighting matrix.
   @param[in] b The right vector of length w->ny.
   @return a'*w*b, or 0 if input is invalid.
*/
T X(wdot)(const T* a, const X(mat)* w, const T* b){
	if(!check_mat(w)||!a||!b) return 0;
	TD res=0;
	for(long j=0; j<w->ny; j++){
		OMP_SIMD_R(reduction(+:res))
		for(long i=0; i<w->nx; i++){
			res+=P(w, i, j)*a[i]*b[j];
		}
	}
	if(isnan(REAL(res))){
		writebin(w, "wdot_w");
		writearr("wdot_a", 1, sizeof(T), M_T, NULL, a, w->nx, 1);
		writearr("wdot_b", 1, sizeof(T), M_T, NULL, b, w->ny, 1);
		error("NaN found.\n");
	}
	return (T)res;
}

/**
   Compute component wise multiply A=A.*B

   @param[in,out] A The matrix to multiply in place.
   @param[in] B The multiplier. Must match A in size.
*/
void X(cwm)(X(mat)* restrict A, const X(mat)* restrict B){
	if(!check_mat(A, B)||!check_match(A, B)){
		error("Input is not valid\n");
		return;
	}
OMP_SIMD()
	for(long i=0; i<B->nx*B->ny; i++){
		P(A, i)*=P(B, i);
	}
}
/**
   Compute component wise multiply A=A.*(B1*wt1+B2*wt2)
*/
/*
void X(cwm2)(X(mat)* restrict A, const X(mat)* B1, R wt1, const X(mat)* B2, R wt2){
	if(!check_mat(A)){
		warning("A is not valid\n");
		return;
	}
	int has_b1=B1&&wt1&&check_match(A, B1);
	int has_b2=B2&&wt2&&check_match(A, B2);
	if(has_b1&&has_b2){
OMP_SIMD()
		for(long i=0; i<B1->nx*B1->ny; i++){
			P(A, i)*=(P(B1, i)*wt1+P(B2, i)*wt2);
		}
	} else if(has_b1){
OMP_SIMD()
		for(long i=0; i<B1->nx*B1->ny; i++){
			P(A, i)*=P(B1, i)*wt1;
		}
	} else if(has_b2){
OMP_SIMD()
		for(long i=0; i<B2->nx*B2->ny; i++){
			P(A, i)*=P(B2, i)*wt2;
		}
	}
}
*/

/**
   component-wise multiply of three matrices.
   A=A.*W.*(B1*wt1+B2*wt2)
*/
/*void X(cwm3)(X(mat)* restrict A, const X(mat)* restrict W,
	const X(mat)* restrict B1, R wt1, const X(mat)* restrict B2, R wt2){
	if(!check_mat(A)){
		warning("A is not valid\n");
		return;
	}
	if(!W){
		X(cwm2)(A, B1, wt1, B2, wt2);
	} else{
		int has_b1=B1&&wt1&&check_match(A, B1);
		int has_b2=B2&&wt2&&check_match(A, B2);
		if(has_b1&&has_b2){
			OMP_SIMD()
			for(long i=0; i<B1->nx*B1->ny; i++){
				P(A, i)*=P(W, i)*(P(B1, i)*wt1+P(B2, i)*wt2);
			}
		} else if(has_b1){
			OMP_SIMD()
			for(long i=0; i<B1->nx*B1->ny; i++){
				P(A, i)*=P(W, i)*P(B1, i)*wt1;
			}
		} else if(has_b2){
			OMP_SIMD()
			for(long i=0; i<B2->nx*B2->ny; i++){
				P(A, i)*=P(W, i)*P(B2, i)*wt2;
			}
		}
	}
}
*/
/**
   Component-wise multiply each column of A with B
   A(:,i)=A(:,i).*B;

   @param[in,out] A The matrix to multiply in place.
   @param[in] B The column vector. Must have A->nx rows. Nothing is done if it does not match.
*/
void X(cwmcol)(X(mat)* restrict A, const X(mat)* restrict B){
	if(check_mat(A, B)&&A->nx==B->nx&&B->ny==1){
		T* B1=P(B);
		for(long iy=0; iy<A->ny; iy++){
			OMP_SIMD()
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=B1[ix];
			}
		}
	}
}

/**
   component-wise multiply of columns of A with combination of B1 and B2:
   A(:,i)=A(:,i)*(B1*wt1+B2*wt2);
*/
/*void X(cwmcol2)(X(mat)* restrict A,
	const X(mat)* restrict B1, const R wt1,
	const X(mat)* restrict B2, const R wt2){
	if(!check_mat(A)){
		warning("A is not valid\n");
		return;
	}
	int has_b1=check_mat(B1)&&wt1&&(A->nx==B1->nx)&&(B1->ny==1);
	int has_b2=check_mat(B2)&&wt2&&(A->nx==B2->nx)&&(B2->ny==1);
	if(has_b1&&has_b2){
		for(long ix=0; ix<A->nx; ix++){
			T junk=P(B1, ix)*wt1+P(B2, ix)*wt2;
			//OMP_SIMD()
			for(long iy=0; iy<A->ny; iy++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b1){
		for(long ix=0; ix<A->nx; ix++){
			T junk=P(B1, ix)*wt1;
			//OMP_SIMD()
			for(long iy=0; iy<A->ny; iy++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b2){
		X(cwmcol2)(A, B2, wt2, B1, wt1);
	}
}
*/
/**
   component wise multiply of 2d complex matrix A,W and 1d vector B.
   A(:,i)=A(:,i).*W(:,i).*B;
*/
/*void X(cwm3col)(X(mat)* restrict A, const X(mat)* restrict W,
	const X(mat)* restrict B1, const R wt1,
	const X(mat)* restrict B2, const R wt2){

	if(!W){
		X(cwmcol2)(A, B1, wt1, B2, wt2);
	} else{
		int has_b1=check_mat(B1)&&wt1&&(A->nx==B1->nx)&&(B1->ny==1);
		int has_b2=check_mat(B2)&&wt2&&(A->nx==B2->nx)&&(B2->ny==1);
		if(has_b1&&has_b2){
			T* B1p=P(B1);
			T* B2p=P(B2);
			for(long iy=0; iy<A->ny; iy++){
				OMP_SIMD()
				for(long ix=0; ix<A->nx; ix++){
					P(A, ix, iy)=P(A, ix, iy)*P(W, ix, iy)*(B1p[ix]*wt1+B2p[ix]*wt2);
				}
			}
		} else if(has_b1){
			T* B1p=P(B1);
			for(long iy=0; iy<A->ny; iy++){
				OMP_SIMD()
				for(long ix=0; ix<A->nx; ix++){
					P(A, ix, iy)=P(A, ix, iy)*P(W, ix, iy)*B1p[ix]*wt1;
				}
			}
		} else if(has_b2){
			X(cwm3col)(A, W, B2, wt2, B1, wt1);
		}
	}
}*/
/**
   Component wise multiply each row of A with B.
   A(i,:)=A(i,:).*B'

   @param[in,out] A The matrix to multiply in place.
   @param[in] B The column vector with A->ny elements. Element i multiplies column i of A.
*/
void X(cwmrow)(X(mat)* restrict A, const X(mat)* restrict B){
	if(!A||!B||!check(A->ny==B->nx&&B->ny==1)){
		warning("Input is not valid\n");
		return;
	}
	T* B1=P(B);
	for(long iy=0; iy<A->ny; iy++){
		T junk=B1[iy];
		OMP_SIMD()
		for(long ix=0; ix<A->nx; ix++){
			P(A, ix, iy)*=junk;
		}
	}
}

/**
   component-wise multiply of rows of A with combination of B1 and B2:
   A(i,:)=A(i,:).*(B1*wt1+B2*wt2);
*/
/*void X(cwmrow2)(X(mat)* restrict A,
	const X(mat)* restrict B1, const R wt1,
	const X(mat)* restrict B2, const R wt2){
	if(!check_mat(A)){
		warning("A is not valid\n");
		return;
	}
	int has_b1=check_mat(B1)&&wt1&&check(A->nx==B1->ny&&B1->nx==1);
	int has_b2=check_mat(B2)&&wt2&&check(A->nx==B2->ny&&B2->nx==1);

	if(has_b1&&has_b2){
		for(long iy=0; iy<A->ny; iy++){
			T junk=P(B1, iy)*wt1+P(B2, iy)*wt2;
			OMP_SIMD()
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b1){
		for(long iy=0; iy<A->ny; iy++){
			T junk=P(B1, iy)*wt1;
			OMP_SIMD()
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b2){
		X(cwmrow2)(A, B2, wt2, B1, wt1);
	}
}
*/
/**
   Component wise division B=B./A. NaN results (e.g., 0/0) are replaced by 'value'.

   @param[in,out] B The numerator, replaced by the result.
   @param[in] A The denominator.
   @param[in] value The value to use when the result is NaN (e.g., 0/0).
*/
void X(cwdiv)(X(mat)* B, const X(mat)* A, T value){
	if(!check_match(A, B)){
		warning("Input does not match\n");
		return;
	}
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(B, i)/=P(A, i);
		if(isnan(REAL(P(B, i)))) P(B, i)=value;
	}
}
/**
   multiply a X(mat) matrix with a vector and accumulate to y:
   y+=A*x*alpha

   @param[in,out] y The output vector of length A->nx, accumulated into.
   @param[in] A The matrix.
   @param[in] x The input vector of length A->ny.
   @param[in] alpha The scaling factor.
*/
void X(mulvec)(T* restrict y, const X(mat)* restrict A, const T* restrict x, const T alpha){
	if(!check(y&&x&&A)){
		warning("input is not valid\n");
		return;
	}
	for(long iy=0; iy<A->ny; iy++){
		OMP_SIMD()
		for(long ix=0; ix<A->nx; ix++){
			y[ix]+=(x[iy]*alpha)*P(A, ix, iy);
		}
	}
}

/**
   Matrix vector multiply optimized for 3x3 matrix. y=A*x;

   @param[out] y The output vector of length 3.
   @param[in] A The 3x3 matrix.
   @param[in] x The input vector of length 3.
*/
void X(mulvec3)(T* y, const X(mat)* A, const T* x){
	if(!check_mat(A)||!check(A->nx==3&&A->ny==3)){
		warning("Input is not valid");
		return;
	}
	/*calculate y=A*x for 3. */
	y[0]=P(A, 0, 0)*x[0]+P(A, 0, 1)*x[1]+P(A, 0, 2)*x[2];
	y[1]=P(A, 1, 0)*x[0]+P(A, 1, 1)*x[1]+P(A, 1, 2)*x[2];
	y[2]=P(A, 2, 0)*x[0]+P(A, 2, 1)*x[1]+P(A, 2, 2)*x[2];
}

/**
   Compute (A'*W*A), where diag(W)=wt. No conjugation is applied.

   @param[in] A The matrix.
   @param[in] wt The weighting vector of length A->nx.
   @return The A->ny by A->ny result, or NULL if input is invalid. Newly allocated; the caller frees it.
*/
X(mat)* X(mcc)(const X(mat)* A, const X(mat)* wt){
	if(!check_mat(A)||!check(A->nx==wt->nx&&wt->ny==1)) return NULL;
	long nmod=A->ny;
	long nsa2=A->nx;
	X(mat)* ata=X(new)(nmod, nmod);;
	for(long imod=0; imod<nmod; imod++){
		for(long jmod=imod; jmod<nmod; jmod++){
			T tmp=0;
			OMP_SIMD_R(reduction(+:tmp))
			for(long ik=0; ik<nsa2; ik++){
				tmp+=P(A, ik, imod)*P(A, ik, jmod)*P(wt, ik);
			}
			P(ata, jmod, imod)=tmp;
			if(imod!=jmod){
				P(ata, imod, jmod)=P(ata, jmod, imod);
			}
		}
	}
	return ata;
}

/**
   Compute (A*W*A'), where diag(W)=wt. No conjugation is applied.

   @param[in] A The matrix.
   @param[in] wt The weighting vector of length A->ny.
   @return The A->nx by A->nx result, or NULL if input is invalid. Newly allocated; the caller frees it.
*/
X(mat)* X(tmcc)(const X(mat)* A, const X(mat)* wt){
	if(!check_mat(A)||!check(A->ny==wt->nx&&wt->ny==1)) return NULL;
	long nmod=A->nx;
	long nsa2=A->ny;
	X(mat)* ata=X(new)(nmod, nmod);;
	for(long imod=0; imod<nmod; imod++){
		for(long jmod=imod; jmod<nmod; jmod++){
			T tmp=0;
			OMP_SIMD_R(reduction(+:tmp))
			for(long k=0; k<nsa2; k++){
				tmp+=P(A, imod, k)*P(A, jmod, k)*P(wt, k);
			}
			P(ata, jmod, imod)=tmp;
			if(imod!=jmod)
				P(ata, imod, jmod)=P(ata, jmod, imod);
		}
	}
	return ata;
}

/**
   Compute the relative difference between two matrices as
   sqrt(sum(abs(A-B).^2)*2/(X(norm)(A)+X(norm)(B))). For debugging purpose.

   @param[in] A The first matrix.
   @param[in] B The second matrix.
   @return The relative difference, or -1 if A and B do not match or the result is NaN.
*/
T X(diff)(const X(mat)* A, const X(mat)* B){
	if(!check_match(A, B)) return -1;
	T d=sqrt(X(sumdiffsq)(A, B)*2/(X(norm)(A)+X(norm)(B)));
	return isnan(REAL(d))?-1:d;
}
/**
   Generate a new gray pixel map based on bilinear influence functions used in
   mkw.  It creates slightly larger map than an filled circle to account for the
   sampling effect.  The Center and Radius cx,cy,r are in unit of meter, The
   sampling dx,dy specify spacing of the points in meter along x or y dimension.
   Each full point received value of 'val'

   @param[in,out] A The map to accumulate into.
   @param[in] cx Center x coordinate (m).
   @param[in] cy Center y coordinate (m).
   @param[in] dx Sampling along x (m).
   @param[in] dy Sampling along y (m).
   @param[in] r Radius (m).
   @param[in] val Value of fully covered points.
   @return 0 if succeeded, -1 if input is invalid.
*/
int X(circle)(X(mat)* A, R cx, R cy, R dx, R dy, R r, T val){
	if(!check_mat(A)||r<0){
		warning("Input is not valid\n");
		return -1;
	}
	const int nres=10;
	const R res=(R)(2./nres);
	const R res2=(R)(res*res);
	R resm=(R)((nres-1)*0.5);
	R r2=r*r;
	const R dmax=MAX(fabs(dx), fabs(dy));//margins below are in unit of pixel
	R r2l=r>1.5*dmax?((r-1.5*dmax)*(r-1.5*dmax)):0;
	R r2u=(r+2.5*dmax)*(r+2.5*dmax);
	for(int iy=0; iy<A->ny; iy++){
		R r2y=(iy*dy-cy)*(iy*dy-cy);
		for(int ix=0; ix<A->nx; ix++){
			R r2r=(ix*dx-cx)*(ix*dx-cx)+r2y;
			T val2=0;
			if(r2r<r2l){
				val2=val;
			} else if(r2r<r2u){
				T tot=0.;
				for(int jy=0; jy<nres; jy++){
					R iiy=iy+(jy-resm)*res;
					R rr2y=(iiy*dy-cy)*(iiy*dy-cy);
					R wty=1.-fabs(iy-iiy);//takes into account bilinear influence function
OMP_SIMD()
					for(int jx=0; jx<nres; jx++){
						R iix=ix+(jx-resm)*res;
						R rr2r=(iix*dx-cx)*(iix*dx-cx)+rr2y;
						R wtx=1.-fabs(ix-iix);
						if(rr2r<r2){
							tot+=res2*wty*wtx;
						}
					}
				}
				val2=tot*val;
			}
			P(A, ix, iy)+=val2;
		}
	}
	return 0;
}
/**
   Generate a new gray pixel map based on bilinear influence functions used in
   mkw.  It creates slightly larger map than an filled rectangle to account for
   the sampling effect. The rectangle is rotated by theta. The center and half
   width cx,cy,rx,ry are in unit of pixel. Each full point received value of
   'val', which is added to A.

   @param[in,out] A The map to accumulate into.
   @param[in] cx Center x coordinate (pixel).
   @param[in] cy Center y coordinate (pixel).
   @param[in] rx Half width along x before rotation (pixel).
   @param[in] ry Half width along y before rotation (pixel).
   @param[in] theta Rotation angle (radian).
   @param[in] val Value of fully covered points.
   @return 0 if succeeded, -1 if input is invalid.
*/
int X(rectangle)(X(mat) *A, R cx, R cy, R rx, R ry, R theta, T val){
	if(!check_mat(A)||rx<0||ry<0){
		warning("Input is not valid\n");
		return -1;
	}
	const int nres=12;//super-resolution. Even and not a multiple of 4 so that no sub-sample falls on integer or half-integer offsets.
	const R res=(R)(2./nres);
	const R res2=(R)(res*res);
	const R cnt=cos(theta);
	const R snt=sin(theta);
	
	//The bilinear influence function has support [-1,1]x[-1,1] in pixel
	//coordinate, which extends by up to ext along the rotated axes.
	const R ext=fabs(cnt)+fabs(snt);
	R resm=(R)((nres-1)*0.5);
	for(int iy=0; iy<A->ny; iy++){
		for(int ix=0; ix<A->nx; ix++){
			R xi0=ix-cx;
			R yi0=iy-cy;
			R xi=xi0*cnt+yi0*snt;
			R yi=-xi0*snt+yi0*cnt;
			T val2=0;
			if(fabs(xi)<=rx-ext&&fabs(yi)<=ry-ext){//support completely inside
				val2=val;				
			} else if(fabs(xi)<rx+ext&&fabs(yi)<ry+ext){//partially inside
				T tot=0.;
				for(int jy=0; jy<nres; jy++){
					OMP_SIMD()
					for(int jx=0; jx<nres; jx++){
						R xj0=ix+(jx-resm)*res-cx;
						R yj0=iy+(jy-resm)*res-cy;
						R xj=xj0*cnt+yj0*snt;
						R yj=-xj0*snt+yj0*cnt;
						R wtx=1.-fabs(xi0-xj0);
						R wty=1.-fabs(yi0-yj0);

						if(fabs(xj)<rx && fabs(yj)<ry){
							tot+=res2*wty*wtx;
						}
					}
				}
				val2=tot*val;
			}
			P(A, ix, iy)+=val2;
		}
	}
	return 0;
}

/**
   rotate the column vectors CCW, equivalent as rotate coordinate theta CW.
   A is nx2 while A(:,0) is x, A(:,1) is y.

   @param[in,out] A The nx2 matrix of vectors, rotated in place.
   @param[in] theta The rotation angle (radian).
*/
void X(rotvec)(X(mat)* A, const R theta){
	if(!check_mat(A)||A->ny!=2){
		warning("A is not valid");
		return;
	}
	const R ctheta=cos(theta);
	const R stheta=sin(theta);
	OMP_SIMD()
	for(long i=0; i<A->nx; i++){
		T tmp=P(A, i, 0)*ctheta-P(A, i, 1)*stheta;
		P(A, i, 1)=P(A, i, 0)*stheta+P(A, i, 1)*ctheta;
		P(A, i, 0)=tmp;
	}
}

/**
   rotate the row vectors CCW. same as rotate coordinate theta CW.
   A is 2xn A(0, :) is x, A(1, :) is y.

   @param[in,out] A The 2xn matrix of vectors, rotated in place.
   @param[in] theta The rotation angle (radian).
*/
void X(rotvect)(X(mat)* A, const R theta){
	if(!check_mat(A)||A->nx!=2){
		warning("A has wrong dimension. Canceled.\n");
		return;
	}
	const R ctheta=cos(theta);
	const R stheta=sin(theta);
	OMP_SIMD()
	for(long i=0; i<A->ny; i++){
		T tmp=P(A, 0, i)*ctheta-P(A, 1, i)*stheta;
		P(A, 1, i)=P(A, 0, i)*stheta+P(A, 1, i)*ctheta;
		P(A, 0, i)=tmp;
	}
}

/**
   rotate a 2x2 covariance matrix A by theta CCW
   (coordinate rotate -theta CCW) or from ra to xy
   coordinate.  R*A*R';

   @param[in,out] B0 The 2x2 output matrix. Allocated if NULL. Overwritten.
   @param[in] A The 2x2 input matrix.
   @param[in] theta The rotation angle (radian).
*/
void X(rotvecnn)(X(mat)** B0, const X(mat)* A, R theta){
	if(!check_mat(A)||!check(A->nx==2&&A->ny==2)){
		warning("A has wrong dimension\n");
		return;
	}
	if(X(init)(B0, 2, 2)) return;
	X(mat)* B=*B0;
	const T ctheta=cos(theta);
	const T stheta=sin(theta);
	X(mat)* tmp=X(new)(2, 2);
	/*first apply left R */
	P(tmp, 0, 0)=ctheta*P(A, 0, 0)-stheta*P(A, 1, 0);
	P(tmp, 0, 1)=ctheta*P(A, 0, 1)-stheta*P(A, 1, 1);
	P(tmp, 1, 0)=stheta*P(A, 0, 0)+ctheta*P(A, 1, 0);
	P(tmp, 1, 1)=stheta*P(A, 0, 1)+ctheta*P(A, 1, 1);
	/*then apply right R' */

	P(B, 0, 0)=ctheta*P(tmp, 0, 0)-stheta*P(tmp, 0, 1);
	P(B, 0, 1)=stheta*P(tmp, 0, 0)+ctheta*P(tmp, 0, 1);
	P(B, 1, 0)=ctheta*P(tmp, 1, 0)-stheta*P(tmp, 1, 1);
	P(B, 1, 1)=stheta*P(tmp, 1, 0)+ctheta*P(tmp, 1, 1);
	X(free)(tmp);
}

/**
   Compute the 2D cross-correlation of A and B, out(offx,offy)=sum(A(x,y)*B(x-offx,y-offy)), with zero
   offset at the center of out.

   @param[in,out] pout The correlation output with odd dimensions. Allocated as (2nx-1)x(2ny-1) if NULL. Overwritten.
   @param[in] A The first matrix.
   @param[in] B The second matrix. Must match A in size.
*/
void X(corr)(X(mat)** pout, const X(mat)* A, const X(mat)* B){
	if(!check_match(A, B)){
		warning("Empty matrix or mismatch\n");
		return;
	}
	if(!*pout){
		*pout=X(new)(A->nx*2-1, A->ny*2-1);
	} else{
		if(((*pout)->nx&1)!=1||((*pout)->ny&1)!=1){
			error("output dimension shall be odd\n");
			return;
		}
	}
	X(mat)* out=*pout;
	const long offx2=(out->nx-1)>>1;
	const long offy2=(out->ny-1)>>1;

	for(long offy=-offy2; offy<=offy2; offy++){
		long sy1, nny;
#define SHIFT_PEX(sy1, nny, ny, offy) if(offy>0){sy1=offy; nny=ny;}else{sy1=0;nny=ny+offy;}
		SHIFT_PEX(sy1, nny, A->ny, offy);
		for(long offx=-offx2; offx<=offx2; offx++){
			long sx1, nnx;
			SHIFT_PEX(sx1, nnx, A->nx, offx);
			T tmp=0;
			for(long iy1=sy1; iy1<nny; iy1++){
				OMP_SIMD_R(reduction(+:tmp))
				for(long ix1=sx1; ix1<nnx; ix1++){
					tmp+=P(A, ix1, iy1)*P(B, ix1-offx, iy1-offy);
				}
			}
			P(out, offx+offx2, offy+offy2)=tmp;
		}
	}
}
/**
   Compute the parabolic fit center using 3x3 points around the peak.

   @param[out] grad The peak location (x, y) relative to the center of corr (pixel).
   @param[in] corr The correlation map.
*/
void X(para3)(R* grad, const X(mat)* corr){
	R valmax=0;
	int jy=0, jx=0;
	//Find Peak location (jx, jy)
	for(int iy=1; iy<corr->ny-1; iy++){
		for(int ix=1; ix<corr->nx-1; ix++){
			if(REAL(P(corr, ix, iy))>valmax){
				jy=iy; jx=ix;
				valmax=REAL(P(corr, ix, iy));
			}
		}
	}
	//Calculate 1d sum of 3 row/columns.
	R vx[3], vy[3];
	for(long iy=0; iy<3; iy++){
		vy[iy]=0; vx[iy]=0;
		for(long ix=0; ix<3; ix++){
			vy[iy]+=REAL(P(corr, ix+jx-1, iy+jy-1));
			vx[iy]+=REAL(P(corr, iy+jx-1, ix+jy-1));
		}
	}
	//Parabolic fit.
	R px[2], py[2];
	px[0]=(vx[0]+vx[2])*0.5-vx[1];
	py[0]=(vy[0]+vy[2])*0.5-vy[1];
	px[1]=(vx[2]-vx[0])*0.5;
	py[1]=(vy[2]-vy[0])*0.5;
	//Center
	grad[0]=px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
	grad[1]=py[0]==0?0:(-py[1]/(2*py[0])+jy-(corr->ny-1)*0.5);
}

/**
   Compute thresholded center of gravity. The threshold
   is absolute value. bkgrnd is removed from i0 when
   computing cog.  offset is the offset of the reference
   point (cog=0) from the physical center.
   all length are given in terms of pixel. grad is set to 0 if
   the sum is not above thres.

   @param[out] grad The computed CoG (x, y) (pixel).
   @param[in] im The image.
   @param[in] offsetx Offset of the reference point along x (pixel).
   @param[in] offsety Offset of the reference point along y (pixel).
   @param[in] thres The absolute threshold applied to the background removed image. Ignored if mask is set.
   @param[in] bkgrnd The background to remove from each pixel.
   @param[in] flux If nonzero, used instead of the sum of the image as the denominator to linearize CoG.
   @param[in] mask Optional weighting mask. Can be NULL.
*/
void X(cog)(R* grad, const X(mat)* im, R offsetx, R offsety, R thres, R bkgrnd, R flux, const X(mat) *mask){
	R sum=0, sumx=0, sumy=0;
	R iI;
	if(mask){
		if(thres){
			thres=0;
			warning_once("X(cog): threshold is ignored when mask is given\n");
		}
		for(long iy=0; iy<im->ny; iy++){
			OMP_SIMD_R(reduction(+:sum, sumx, sumy))
			for(long ix=0; ix<im->nx; ix++){
				iI=(REAL(P(im, ix, iy))-bkgrnd)*REAL(P(mask, ix, iy));
				sum+=iI;
				sumx+=iI*ix;
				sumy+=iI*iy;
			}
		}
	}else{
		for(long iy=0; iy<im->ny; iy++){
			OMP_SIMD_R(reduction(+:sum, sumx, sumy))
			for(long ix=0; ix<im->nx; ix++){
				iI=REAL(P(im, ix, iy))-bkgrnd;
				if(iI>thres){
					sum+=iI;
					sumx+=iI*ix;
					sumy+=iI*iy;
				}
			}
		}
	}
	if(flux){//linerize CoG by overriding sum
		sum=flux;
	}
	if(sum>thres){
		grad[0]=sumx/sum-((R)(im->nx-1)*0.5+offsetx);
		grad[1]=sumy/sum-((R)(im->ny-1)*0.5+offsety);
	} else{
		grad[0]=0;
		grad[1]=0;
	}
}


/**
   Shift the image in A to center on physical
   center+[offsetx,offsety] using cog and fft.

   @param[in,out] A The image, shifted in place.
   @param[in] offsetx Target offset along x (pixel).
   @param[in] offsety Target offset along y (pixel).
*/
void X(shift2center)(X(mat)* A, R offsetx, R offsety){
	R grad[2];
	R Amax=X(max)(A);
	X(cog)(grad, A, offsetx, offsety, Amax*0.1, Amax*0.2, 0, NULL);
	if(fabs(grad[0])>0.1||fabs(grad[1])>0.1){
	/*dbg("Before shift, residual grad is %g %g\n",grad[0],grad[1]); */
		XC(mat)* B=XC(new)(A->nx, A->ny);
		//XC(fft2plan)(B,-1);
		//XC(fft2plan)(B,1);
#ifdef COMP_COMPLEX
		XC(cp)(&B, A);
#else
		XC(cpd)(&B, A);
#endif
		R scale=1./(A->nx*A->ny);
		XC(fftshift)(B);
		XC(fft2)(B, -1);
		XC(tilt)(B, -grad[0], -grad[1], 0);
		XC(fft2)(B, 1);
		XC(fftshift)(B);
		XC(scale)(B, scale);
#ifdef COMP_COMPLEX
		XC(cp)(&A, B);
#else
		XC(real2d)(&A, 0, B, 1);
#endif
		X(cog)(grad, A, offsetx, offsety, Amax*0.1, Amax*0.2, 0, NULL);
		/*dbg("After shift, residual grad is %g %g\n",grad[0],grad[1]); */
		XC(free)(B);
	}
}


/**
   OrthNormalize column vector in Mod, with weighting from vector amp.
   <Mod|wt|Mod> is equal to sum(wt).
   2010-07-21: Bug found: The result was not orthonormal. cause: nonvalid was not initialized to 0.

   @param[in,out] Mod The column vectors, orthonormalized in place.
   @param[in] amp The weighting vector of length Mod->nx. Can be NULL for uniform weighting.
*/
void X(gramschmidt)(X(mat)*restrict Mod, R* restrict amp){
	const long nmod=Mod->ny;
	const long nx=Mod->nx;
	R wtsum=(R)nx;
	if(amp){
		wtsum=XR(vecsum)(amp, nx);
	}
	int nonvalid[nmod];
	memset(nonvalid, 0, sizeof(int)*nmod);
	for(long imod=0; imod<nmod; imod++){
		if(nmod>10){
			info("Gramschmidt: %ld of %ld\n", imod, nmod);
		}
		if(imod>0){/*orthogonalize */
			/*compute dot product. */
			for(long jmod=0; jmod<imod; jmod++){
				if(nonvalid[jmod]) continue;
				TD dot=0;/*<jmod|wt|imod>, conjugated for complex*/
				for(long ix=0; ix<nx; ix++){
					dot+=CONJ(P(Mod, ix, jmod))*P(Mod, ix, imod)*(amp?amp[ix]:1);
				}
				const T cross=(T)(-dot/(double)wtsum);
				OMP_SIMD()
				for(long ix=0; ix<nx; ix++){
					P(Mod, ix, imod)+=cross*P(Mod, ix, jmod);
				}
			}
		}
		/*normalize*/
		real norm2=0;/*<imod|wt|imod>*/
		for(long ix=0; ix<nx; ix++){
			const R mag=ABS(P(Mod, ix, imod));
			norm2+=mag*mag*(amp?amp[ix]:1);
		}
		R norm=(R)sqrt(norm2/wtsum);
		if(fabs(norm)>1.e-15){
			norm=1./norm;
			OMP_SIMD()
			for(long ix=0; ix<nx; ix++){
				P(Mod, ix, imod)*=norm;
			}
		} else{
			nonvalid[imod]=1;
			warning("Column %ld is not independent on other columns\n", imod);
		}
	}
}


/**
   Limit numbers in A to within [min, max]. used for DM clipping.

   @param[in,out] A The matrix to clip in place. Only the real part is compared.
   @param[in] min The lower limit.
   @param[in] max The upper limit.
   @return The number of elements clipped.
*/
int X(clip)(X(mat)* restrict A, R min, R max){
	if(!A) return 0;
	if(isinf(min)&&isinf(max)) return 0;
	if(max<=min){
		error("upper light should be larger than lower limit\n");
		return 0;
	}
	T* restrict Ap=P(A);
	int nclip=0;
	OMP_SIMD_R(reduction(+:nclip))
	for(long i=0; i<A->nx*A->ny; i++){
		R Ar=REAL(Ap[i]);
		if(Ar>max){
			Ap[i]=max;
			nclip++;
		} else if(Ar<min){
			Ap[i]=min;
			nclip++;
		}
	}
	return nclip;
}

/**
   A=A*B, where diag(B)=s.

   @param[in,out] A The matrix to multiply in place.
   @param[in] s The diagonal of B with A->ny elements.
*/
void X(muldiag)(X(mat) *restrict A, const XR(mat) *restrict s){
	assert(A->ny==s->nx&&s->ny==1);
	const R *restrict ps=P(s);
	for(long iy=0; iy<A->ny; iy++){
		OMP_SIMD()
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=ps[iy];
			}
	}
}
/**
   Raise all elements to power 'power'.

   @param[in,out] A The matrix, modified in place.
   @param[in] power The exponent.
*/
void X(cwpow)(X(mat)* A, R power){
	if(!A) return;
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=pow(P(A, i), power);
	}
}

/**
   Compute exponential of all elements after scaling by alpha.

   @param[in,out] A The matrix, modified in place.
   @param[in] alpha The scaling factor.
*/
void X(cwexp)(X(mat)* A, R alpha){
	if(!A) return;
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=exp(P(A, i)*alpha);
	}
}

/**
   Raise all elements above thres*maxabs(A) to pow power. Set others to zero.

   @param[in,out] A The matrix, modified in place.
   @param[in] power The exponent.
   @param[in] thres The threshold relative to the maximum absolute value.
*/
void X(cwpow_thres)(X(mat)* restrict A, R power, R thres){
	thres*=X(maxabs)(A);
	//OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		if(ABS(P(A, i))>thres){
			P(A, i)=pow(P(A, i), power);
		} else{
			P(A, i)=0;
		}
	}
}
/**
   Compute polynomial functional value inplace: A=p(0)*A^(np-1)+...+p(np-1).

   @param[in,out] A The input values, replaced by the polynomial value.
   @param[in,out] p The polynomial coefficients, highest order first. Reshaped in place to a column vector if it is a row vector.
*/
void X(polyval)(X(mat)* A, XR(mat)* p){
	if(p->nx==1&&p->ny>1){
		p->nx=p->ny;
		p->ny=1;
	}
	long np=p->nx;
	for(long i=0; i<A->nx*A->ny; i++){
		T tmp=0;
		OMP_SIMD_R(reduction(+:tmp))
		for(long ip=0; ip<np; ip++){
			tmp+=P(p, ip)*(T)pow(P(A, i), np-ip-1);
		}
		P(A, i)=tmp;
	}
}
/**
   Add val to diagonal values of A.

   @param[in,out] A The matrix, modified in place.
   @param[in] val The value to add.
*/
void X(addI)(X(mat)*restrict A, T val){
	if(!A) return;
	if(A->nx!=A->ny)
		warning("daddI: A is not square\n");
	long M=A->nx<A->ny?A->nx:A->ny;
	OMP_SIMD()
	for(long i=0; i<M; i++){
		P(A, i, i)+=val;
	}
}
/**
   compute B=bc*B+ac*A
   behavior changed on 2009-11-02. if A is NULL, don't do anything.

   @param[in,out] B0 The matrix B. Allocated if NULL.
   @param[in] bc The scaling factor of B.
   @param[in] A The matrix to add. Nothing is done if it is NULL or empty.
   @param[in] ac The scaling factor of A. If it is 0, B is only scaled by bc.
*/
void X(add)(X(mat)** B0, T bc, const X(mat)* A, const T ac){
	if(A&&A->nx){
		if(!*B0){
			bc=0;/*no bother to accumulate. */
			X(init)(B0, A->nx, A->ny);
		}else if(PN(A) != PN(*B0)){
			error("A is %ldx%ld, B is %ldx%ld. They should match in size.\n",
				NX(A), NY(A), NX(*B0), NY(*B0));
			return;
		}
		if(ac!=(T)0){
			X(mat) * restrict B=*B0;
			if(bc!=(T)0){
				OMP_SIMD()
				for(long i=0; i<PN(A); i++){
					P(B, i)=P(B, i)*bc+P(A, i)*ac;
				}
			} else{
				if(ac==(T)1){
					X(cp)(B0, A);
				} else{/*just assign */
					OMP_SIMD()
					for(long i=0; i<PN(A); i++){
						P(B, i)=P(A, i)*ac;
					}
				}
			}
		}else if(bc!=(T)1){/*B=bc*B */
			X(mat) * restrict B=*B0;
			OMP_SIMD()
			for(long i=0; i<PN(B); i++){
				P(B, i)*=bc;
			}
		}
	}
}
/**
   Compute B[:,icol]=bc*B[:,icol]+ac*A where A is vector.

   @param[in,out] B The matrix B. Must already exist.
   @param[in] icol The column index of B.
   @param[in] bc The scaling factor of B[:,icol].
   @param[in] A The column vector with B->nx elements.
   @param[in] ac The scaling factor of A. Nothing is done if it is zero.
*/
void X(addcol)(X(mat) *B, long icol, T bc, const X(mat) *A, const T ac){
	if(NX(A)&&ac!=(T)0){
		if(!B||NY(A)>1||icol<0){
			error("only support adding a vector to an existing matrix\n");
			return;
		}else if(NX(A)!=NX(B)||NY(B)<=icol){
			error("Mismatch: A is %ldx%ld, B is %ldx%ld. icol is %ld.\n",
				NX(A), NY(A), NX(B), NY(B), icol);
			return;
		}
		if(bc!=(T)0){
			OMP_SIMD()
				for(long i=0; i<NX(A); i++){
					P(B, i, icol)=P(B, i, icol)*bc+P(A, i)*ac;
				}
		} else{
			if(ac==(T)1){
				OMP_SIMD()
					for(long i=0; i<NX(A); i++){
						P(B, i, icol)=P(A, i);
					}
			} else{/*just assign */
				OMP_SIMD()
					for(long i=0; i<NX(A); i++){
						P(B, i, icol)=P(A, i)*ac;
					}
			}
		}
	}
}
/**
   Compute B=bc*B+ac*A. Use smallest dimension of both.
   behavior changed on 2009-11-02. if A is NULL or ac is 0, don't do anything.

   @param[in,out] B0 The matrix B. Allocated if NULL.
   @param[in] bc The scaling factor of B.
   @param[in] A The matrix to add.
   @param[in] ac The scaling factor of A.
*/
void X(add_relax)(X(mat)** B0, T bc, const X(mat)* A, const T ac){
	if(A&&A->nx&&ac!=(T)0){
		if(!*B0){
			bc=0;/*no bother to accumulate. */
			X(init)(B0, A->nx, A->ny);
		}
		X(mat)* restrict B=*B0;
		long nx=MIN(A->nx, B->nx);
		long ny=MIN(A->ny, B->ny);
		for(long iy=0; iy<ny; iy++){
			OMP_SIMD()
			for(long ix=0; ix<nx; ix++){
				P(B, ix, iy)=P(B, ix, iy)*bc+P(A, ix, iy)*ac;
			}
		}
	}
}
/**
   Add a scalar to each element of matrix.

   @param[in,out] A The matrix, modified in place.
   @param[in] ac The scalar to add.
*/
void X(adds)(X(mat)*restrict A, const T ac){
	if(!A||!A->nx||ac==(T)0) return;
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)+=ac;
	}
}

/**
   Create log spaced vector.

   @param[in] emin log10 of the first value.
   @param[in] emax log10 of the last value.
   @param[in] n Number of values.
   @return The nx1 vector 10.^linspace(emin, emax, n). Newly allocated; the caller frees it.
*/
X(mat)* X(logspace)(R emin, R emax, long n){
	X(mat)* out=X(new)(n, 1);
	R esep=(emax-emin)/(n-1);
	//OMP_SIMD()
	for(long i=0; i<n; i++){
		P(out, i)=pow(10, emin+esep*i);
	}
	return out;
}

/**
   Create linearly spaced vector.

   @param[in] min The first value.
   @param[in] dx The spacing.
   @param[in] n Number of values.
   @return The nx1 vector min+dx*(0:n-1). Newly allocated; the caller frees it.
*/
X(mat)* X(linspace)(R min, R dx, long n){
	X(mat)* out=X(new)(n, 1);
	OMP_SIMD()
	for(long i=0; i<n; i++){
		P(out, i)=min+dx*i;
	}
	return out;
}
#ifndef COMP_COMPLEX
/**
   Check whether xin is linearly spaced. Only the first two values and the last value are checked.

   @param[in] xin The vector of coordinates. May be increasing or decreasing.
   @return The spacing (negative for decreasing) if linearly spaced, 0 otherwise.
*/
static R X(islinear)(const X(mat)* xin){
	if(!xin||xin->nx<2) return 0;
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=(P(xin, 0));
	R xmaxl=(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	if(!(xsep!=0) || !(fabs(xsep+xminl-P(xin, 1))<=fabs(xsep)*1.e-3)){
		return 0;
	} else{
		return xsep;
	}
}
/**
   Check whether xin is logarithmically spaced. Only the first two values and the last value are checked.

   @param[in] xin The vector of coordinates. May be increasing or decreasing.
   @return The spacing in log10 (negative for decreasing) if logarithmically spaced, 0 otherwise.
*/
static R X(islog)(const X(mat)* xin){
	if(!xin||xin->nx<2) return 0;
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=log10(P(xin, 0));
	R x1=log10(P(xin, 1));
	R xmaxl=log10(P(xin, nmax1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	if(!(xsep!=0) || isinf(xsep) || !(fabs(xsep+xminl-x1)<=fabs(xsep)*1.e-3)){
		return 0;
	} else{
		return xsep;
	}
}
/**
   Select the function values for the interp1 routines.

   @param[in] xin The coordinate in the first column, optionally followed by function values.
   @param[in] yin_ The function values, with the same number of rows as xin. Can be NULL.
   @return yin_ if it is not NULL; otherwise a reference to columns 2 and later of xin, which the caller
   frees. NULL if yin_ is NULL and xin has one column, or if the number of rows does not match.
*/
static X(mat)* X(interp_check_y)(const X(mat)* xin, const X(mat)* yin_){
	X(mat)*yin=NULL;
	if(!yin_){
		if(NY(xin)>1){
			yin=X(refcols)(xin, 1, NY(xin)-1);
		}else{
			error("yin is NULL and xin has only one column.\n");
			return NULL;
		}
	}else if(yin_->nx!=xin->nx){
		error("xin and yin must have the same nx.\n");
		return NULL;
	}else{
		yin=(X(mat)*)yin_;
	}
	return yin;
}
/**
   Interpolate using linear interp. xin is the coordinate of yin. xnew is the
   coordinate of the output.

   @param[in] xin The linearly spaced (increasing or decreasing) coordinate in the first column. If yin_ is
   NULL, the other columns are the function values.
   @param[in] yin_ The function values with the same number of rows as xin. Each column is a function. Can be NULL.
   @param[in] xnew The new coordinates. Only the first column is used.
   @param[in] ydefault The value for points outside of xin. If NaN, the nearest end value is used.
   @return The interpolated values, xnew->nx by the number of functions. Newly allocated; the caller frees it. NULL if input is invalid.
*/
X(mat)* X(interp1linear)(const X(mat)* xin, const X(mat)* yin_, const X(mat)* xnew, T ydefault){
	if(!xin||!xnew||NX(xin)<2) return NULL;
	if(!X(islinear)(xin)){
		error("xin is not linearly spaced\n");
		return NULL;
	}
	X(mat)* yin=X(interp_check_y)(xin, yin_);
	if(!yin) return NULL;
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=(P(xin, 0));
	R xmaxl=(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	R xsep1=1./xsep;
	X(mat)* ynew=X(new)(xnew->nx, yin->ny);
	for(long iy=0; iy<ynew->ny; iy++){
		for(long ix=0; ix<ynew->nx; ix++){
			R xx=((P(xnew, ix))-xminl)*xsep1;
			long xxm=ifloor(xx);
			if(xxm==nmax1 && xx-nmax1<=1e-6*nmax1) xxm=nmax1-1;//the last grid point is within range
			if(xxm<0){
				P(ynew, ix, iy)=isnan(ydefault)?P(yin, 0, iy):ydefault;
			} else if(xxm>=nmax1){
				P(ynew, ix, iy)=isnan(ydefault)?P(yin, nmax1, iy):ydefault;
			} else{
				R xxw=xx-xxm;
				P(ynew, ix, iy)=xxw*P(yin, xxm+1, iy)+(1.-xxw)*P(yin, xxm, iy);
			}
		}
	}
	if(yin!=yin_){
		X(free)(yin);
	}
	return ynew;
}

/**
   Interpolate using log(xin) and log(xnew)
   xin is the coordinate of yin. xnew is the coordinate of the output.

   @param[in] xin The logarithmically spaced (increasing or decreasing) coordinate in the first column. If yin_
   is NULL, the other columns are the function values.
   @param[in] yin_ The function values with the same number of rows as xin. Each column is a function. Can be NULL.
   @param[in] xnew The new coordinates. Only the first column is used.
   @param[in] ydefault The value for points outside of xin. If NaN, the nearest end value is used.
   @return The interpolated values, xnew->nx by the number of functions. Newly allocated; the caller frees it. NULL if input is invalid.
*/
X(mat)* X(interp1log)(const X(mat)* xin, const X(mat)* yin_, const X(mat)* xnew, T ydefault){
	if(!xin||!xnew||NX(xin)<2) return NULL;
	if(!X(islog)(xin)){
		error("xin is not logrithmically spaced\n");
		return NULL;
	}
	X(mat)* yin=X(interp_check_y)(xin, yin_);
	if(!yin) return NULL;
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=log10(P(xin, 0));
	R xmaxl=log10(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	R xsep1=1./xsep;
	X(mat)* ynew=X(new)(xnew->nx, yin->ny);
	for(long iy=0; iy<ynew->ny; iy++){
		for(long ix=0; ix<ynew->nx; ix++){
			R xx=(log10(P(xnew, ix))-xminl)*xsep1;
			long xxm=ifloor(xx);
			if(xxm==nmax1 && xx-nmax1<=1e-6*nmax1) xxm=nmax1-1;//the last grid point is within range
			if(xxm<0){
				P(ynew, ix, iy)=isnan(ydefault)?P(yin, 0, iy):ydefault;
			} else if(xxm>=nmax1){
				P(ynew, ix, iy)=isnan(ydefault)?P(yin, nmax1, iy):ydefault;
			} else{
				R xxw=xx-xxm;
				P(ynew, ix, iy)=xxw*P(yin, xxm+1, iy)+(1.-xxw)*P(yin, xxm, iy);
			}
		}
	}
	if(yin!=yin_){
		X(free)(yin);
	}
	return ynew;
}
/**
   Interpolation of 1d array. Uses X(interp1linear) or X(interp1log) if xin is linearly or logarithmically
   spaced, otherwise searches the bin for each point.

   @param[in] xin The coordinate in the first column. If yin_ is NULL, the other columns are the function values.
   @param[in] yin_ The function values with the same number of rows as xin. Each column is a function. Can be NULL.
   @param[in] xnew The new coordinates. Only the first column is used.
   @param[in] ydefault The value for points outside of xin. If NaN, values are extrapolated (linear or log spacing uses the end values instead).
   @return The interpolated values, xnew->nx by the number of functions. Newly allocated; the caller frees it. NULL if input is invalid.
*/
X(mat)* X(interp1)(const X(mat)* xin, const X(mat)* yin_, const X(mat)* xnew, T ydefault){
	if(!xin||!xnew||NX(xin)<2) return NULL;
	X(mat)* restrict ynew=NULL;
	X(mat)* yin=X(interp_check_y)(xin, yin_);
	if(!yin) return NULL;
	if(X(islinear)(xin)){
		ynew=X(interp1linear)(xin, yin, xnew, ydefault);
	} else if(X(islog)(xin)){
		ynew=X(interp1log)(xin, yin, xnew, ydefault);
	} else{//arbitrary spacing
		ynew=X(new)(xnew->nx, yin->ny);
		const long nend=xin->nx-1;
		long curpos=0;
		for(long ix=0; ix<ynew->nx; ix++){
			const R xn=P(xnew, ix);
			int found=0;
			for(curpos=0; curpos<nend; curpos++){//works for increasing or decreasing xin
				if((xn-P(xin, curpos))*(xn-P(xin, curpos+1))<=0){
					found=1;
					break;
				}
			}
			if(!found){//extrapolate using the first or last interval, whichever end is nearer
				curpos=(fabs(xn-P(xin, 0))<fabs(xn-P(xin, nend)))?0:(nend-1);
			}
			if(found||isnan(ydefault)){
				R xx=((P(xnew, ix))-P(xin, curpos))/(P(xin, curpos+1)-P(xin, curpos));
				for(long iy=0; iy<ynew->ny; iy++){
					P(ynew, ix, iy)=xx*P(yin, curpos+1, iy)+(1.-xx)*P(yin, curpos, iy);
				}
			} else{
				for(long iy=0; iy<ynew->ny; iy++){
					P(ynew, ix, iy)=ydefault;
				}
			}
		}
	}
	if(yin!=yin_){
		X(free)(yin);
	}
	return ynew;
}

#endif
#ifndef COMP_COMPLEX
/**
   Blend B into center of A with width of overlap. The center
   (size is B->nx-overlap, B->ny-overlap) of A is replaced by
   center of B. The overlapping area is blended.

   @param[in,out] A The destination matrix.
   @param[in] B The matrix to blend in. Not modified.
   @param[in] overlap Width of the overlapping area (pixel).
*/
void X(blend)(X(mat)* restrict A, X(mat)* restrict B, int overlap){
	const long ninx=B->nx;
	const long niny=B->ny;
	const long noutx=A->nx;
	const long nouty=A->ny;
	const long skipx=(noutx-ninx)/2;
	const long skipy=(nouty-niny)/2;
	long ixstart=0, ixlen=ninx;
	long iystart=0, iylen=niny;

	if(skipx<0){
		ixstart=-skipx;
		ixlen+=2*skipx;
	}
	if(skipy<0){
		iystart=-skipy;
		iylen+=2*skipy;
	}
	X(mat)* pA=A;
	X(mat)* pB=B;
	R wty, wtx;
	for(long iy=0; iy<iylen; iy++){
		T* restrict outi=&P(pA, ixstart+skipx, iystart+skipy+iy);
		T* restrict ini=&P(pB, ixstart, iystart+iy);
		if(iy<overlap){
			wty=(R)iy/(R)(overlap-1);
		} else if(iylen-iy-1<overlap){
			wty=(R)(iylen-iy-1)/(R)(overlap-1);
		} else{
			wty=1;
		}
		OMP_SIMD()
		for(long ix=0; ix<ixlen; ix++){
			if(ix<overlap){
				wtx=(R)ix/(R)(overlap-1);
			} else if(ixlen-ix-1<overlap){
				wtx=(R)(ixlen-ix-1)/(R)(overlap-1);
			} else{
				wtx=1;
			}
			outi[ix]=outi[ix]*(1-wtx*wty)+ini[ix]*wtx*wty;
		}
	}
}
/**
   For each entry in A, call repeatly to collect its histogram, centered at
   center, spaced by spacing, for n bins in total. center is at bin n/2.

   @param[in,out] out The n by numel(A) histogram, accumulated into. Allocated if NULL.
   @param[in] A The values. Each element has its own histogram column.
   @param[in] center The value at the center bin n/2.
   @param[in] spacing The bin spacing.
   @param[in] n Number of bins. Values outside are put in the end bins.
*/
void X(histfill)(X(mat)** out, const X(mat)* A, R center, R spacing, int n){
	if(!A||!P(A)) return;
	long nn=A->nx*A->ny;
	if(X(init)(out, n, nn)) return;
	X(mat)* Op=*out;
	const T* restrict Ap=P(A);
	const R spacingi=1./spacing;
	const int noff=n/2;
	const int n1=n-1;
	for(long i=0; i<A->nx*A->ny; i++){
		int ind=(int)round(REAL(Ap[i]-center)*spacingi)+noff;
		if(ind<0) ind=0;
		if(ind>n1) ind=n1;
		P(Op, ind, i)++;
	}
}

/**
   1D Cubic spline interpolation preparation.
   if x has only 1 column: x is the coordinate. y is the function value.
   if x has two columns: first column is the coordinate, y is null.

   It is upto the user to make sure that the coordinate is increasingly ordered
   and evenly spaced.

   If the values of a function \f$f(x)\f$ and its derivative are know at x=0,
   and x=1 (normalized coordinate), then the function can be interpolated on the
   interval [0,1] using a third degree polynomial. This is called cubic
   interpolation. The formula of this polynomial can be easily derived.

   A third degree polynomial and its derivative:
   \f[
   f(x)=ax^3+bx^2+cx+d
   \f]
   \f[
   f^\prime(x)=3ax^2+2bx+c
   \f]
   The coefficients can be derived from the value and derivatives:
   \f{eqnarray*}{
   a&=&2f(0)-2f(1)+f^\prime (0)+f^\prime(1)\\
   b&=&-3f(0)+3f(1)-2f^\prime(0)-f^\prime(1)\\
   c&=&f^\prime(0)\\
   d&=&f(0)\\
   \f}
   the derivatives can be computed as
   \f{eqnarray*}{
   f^\prime(0)&=&\frac{f(1)-f(-1)}{2}\\
   f^\prime(1)&=&\frac{f(2)-f(0)}{2}\\
   \f}
   so we have the formula
   \f{eqnarray*}{
   a&=&-0.5 f(-1) + 1.5 f(0) - 1.5 f(1) + 0.5 f(2)\\
   b&=&     f(-1) - 2.5 f(0) + 2   f(1) - 0.5 f(2)\\
   c&=&-0.5 f(-1)            + 0.5 f(1)           \\
   d&=&                 f(0)                      \\
   \f}

   for the boundary pints, replace
   \f[f^\prime(0)=(f(1)-f(-1))/2\f] by
   \f[f^\prime(0)=(f(1)-f(0))\f]
   Otehr type of boundaries are handled in the same way.

   see www.paulinternet.nl/?page=bicubicx

   @param[in] x The coordinate vector, or two columns of coordinate and function value.
   @param[in] y The function value. Must be NULL if x has two columns.
   @return The 4xnx coefficients (a, b, c, d) of each interval. Newly allocated; the caller frees it.
*/
X(mat)* X(spline_prep)(X(mat)* x, X(mat)* y){
	T* px, * py;
	if(!x) return NULL;
	const long nx=x->nx;
	px=P(x);
	switch(x->ny){
	case 1:
		py=(y&&PN(y)==nx)?P(y):NULL;
		break;
	case 2:
		assert(y==NULL);
		py=P(x)+nx;
		break;
	default:
		py=NULL;
	}
	if(!py||nx<2){
		error("Invalid input\n");
		return NULL;
	}
	X(mat)* coeff=X(new)(4, nx);
	T xsep=(px[nx-1]-px[0])/(nx-1);
	R thres=fabs(xsep)*1.e-5;

	X(mat)* pc=coeff;
	T ypriv, ynext;
	for(long ix=0; ix<nx-1; ix++){
		if(fabs(px[ix+1]-px[ix]-xsep)>thres){
			error("The coordinate is not evenly spaced\n");
			X(free)(coeff);
			return NULL;
		}
		if(UNLIKELY(ix==0)){
			ypriv=2*py[ix]-py[ix+1];
		} else{
			ypriv=py[ix-1];
		}
		if(UNLIKELY(ix==nx-2)){
			ynext=2*py[ix+1]-py[ix];
		} else{
			ynext=py[ix+2];
		}
		P(pc, 0, ix)=-0.5*ypriv+1.5*py[ix]-1.5*py[ix+1]+0.5*ynext;/*a */
		P(pc, 1, ix)=ypriv-2.5*py[ix]+2.0*py[ix+1]-0.5*ynext;/*b */
		P(pc, 2, ix)=-0.5*ypriv+0.5*py[ix+1];/*c */
		P(pc, 3, ix)=py[ix];/*d */
		/*
		  For any point within this bin, with normalized coordinate t (0<t<1);
		  y(t)=a*pow(t,3)+b*pow(t,2)+c*t+d;
		*/
	}
	return coeff;
}
/**
   Evaluate the cubic spline represented by 4xnx matrix coeff, at location array xnew.

   @param[in] coeff The 4xnx coefficients from X(spline_prep).
   @param[in] x The original coordinate vector.
   @param[in] xnew The coordinate to evaluate at.
   @return The values at xnew with the shape of xnew. Newly allocated; the caller frees it.
*/
X(mat)* X(spline_eval)(X(mat)* coeff, X(mat)* x, X(mat)* xnew){
	assert(coeff->nx==4);
	const long nx=coeff->ny;
	X(mat)* pc=coeff;
	T xmin=P(x, 0);
	T xsep1=(T)(nx-1)/(P(x, nx-1)-xmin);
	X(mat)* restrict out=X(new)(xnew->nx, xnew->ny);
	OMP_SIMD()
	for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
		R xn=REAL((P(xnew, ix)-xmin)*xsep1);
		long xnf=floor(xn);
		if(xnf<0) xnf=0;
		if(xnf>nx-2) xnf=nx-2;
		xn=xn-xnf;
		T xn2=xn*xn;
		T xn3=xn2*xn;
		P(out, ix)=P(pc, 0, xnf)*xn3+P(pc, 1, xnf)*xn2+P(pc, 2, xnf)*xn+P(pc, 3, xnf);
	}
	return out;
}
/**
   Do 1D cubic spline all at once by calling X(spline_prep) and X(spline_eval).

   @param[in] x The coordinate vector, or two columns of coordinate and function value.
   @param[in] y The function value. Must be NULL if x has two columns.
   @param[in] xnew The coordinate to evaluate at.
   @return The values at xnew. Newly allocated; the caller frees it.
*/
X(mat)* X(spline)(X(mat)* x, X(mat)* y, X(mat)* xnew){
	X(mat)* coeff=X(spline_prep)(x, y);
	X(mat)* out=X(spline_eval)(coeff, x, xnew);
	X(free)(coeff);
	return out;
}
#endif
/**
   Do a component wise log10 on each element of A.

   @param[in,out] A The matrix, modified in place.
*/
void X(cwlog10)(X(mat)* A){
	R ratio=1./log(10);
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=log(P(A, i))*ratio;
	}
}
/**
   Do a component wise natural log on each element of A.

   @param[in,out] A The matrix, modified in place.
*/
void X(cwlog)(X(mat)* A){
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=log(P(A, i));
	}
}

#define EMBED_OP(out, in, theta, OPin)\
const long ninx=in->nx;\
const long niny=in->ny;\
const long noutx=out->nx;\
const long nouty=out->ny;\
const long ninx2=ninx/2;\
const long noutx2=noutx/2;\
const long niny2=niny/2;\
const long nouty2=nouty/2;\
X(zero)(out);\
if(theta==0){/*no rotation.*/ \
	long skipx, ixstart, ixend, skipy, iystart, iyend;\
	skipx=noutx2-ninx2;/*skip in output array*/\
	if(ninx>noutx){/*reduce*/\
		ixstart=-skipx;/*start in intput array*/\
		ixend=ixstart+noutx;/*ninx+skipx;*/\
	} else{\
		ixstart=0, ixend=ninx;\
	}\
	skipy=nouty2-niny2;\
	if(niny>nouty){\
		iystart=-skipy;\
		iyend=iystart+nouty;\
	} else{\
		iystart=0, iyend=niny;\
	}\
	for(long iy=iystart; iy<iyend; iy++){\
		for(long ix=ixstart; ix<ixend; ix++){\
			P(out, ix+skipx, iy+skipy)=OPin(P(in, ix, iy));\
		}\
	}\
} else{\
	const R ctheta=cos(theta);\
	const R stheta=sin(theta);\
	long ix2, iy2;\
	for(long iy=0; iy<nouty; iy++){\
		R y=(R)(iy-nouty2);\
		for(long ix=0; ix<noutx; ix++){\
			R x=(R)(ix-noutx2);\
			R x2=x*ctheta+y*stheta+ninx2;\
			R y2=-x*stheta+y*ctheta+niny2;\
			if(x2>=0&&x2<=ninx-1&&y2>=0&&y2<=niny-1){\
				ix2=ifloor(x2);\
				iy2=ifloor(y2);\
				x2=x2-ix2;\
				y2=y2-iy2;\
				P(out, ix, iy)=OPin(\
					+P(in, ix2, iy2)*((1.-x2)*(1.-y2))\
						+(x2>0?(P(in, ix2+1, iy2)*(x2*(1.-y2))):0)\
					+(y2>0?(P(in, ix2, iy2+1)*((1.-x2)*y2)\
						+(x2>0?P(in, ix2+1, iy2+1)*(x2*y2):0)):0));\
			}\
		}\
	}\
}
/**
   Embed a ninx*niny matrix in into out with optional rotation by theta CCW
   (output points are rotated by theta -CCW) around the fft center. Used to rotate the PSF
   from x-y to radial-azimuthal coordinate in radial format CCD. out may be bigger
   or smaller than in.

   @param[out] out The output matrix. Zeroed first.
   @param[in] in The input matrix.
   @param[in] theta The rotation angle (radian).
*/
void X(embed)(X(mat)* restrict out, const X(mat)* restrict in, const R theta){
	EMBED_OP(out, in, theta, );
}
/**
   Embed a real ninx*niny matrix in into out with optional rotation by -theta CCW
   (coordinate rotate theta CCW) around the fft center. Used to rotate the PSF
   from x-y to radial-azimuthal coordinate in radial format CCD. out may be bigger or smaller than in.

   @param[out] out The output matrix. Zeroed first.
   @param[in] in The real input matrix.
   @param[in] theta The rotation angle (radian).
*/
void X(embedd)(X(mat)* restrict out, const XR(mat)* restrict in, const R theta){
	EMBED_OP(out, in, theta, )
}
/**
   Calculate number of pixels having values larger than or equal to half of
   maximum and convert to diameter.

   For more accurate calculation of Gaussian beams, use the 2*sqrt(2*log(2))*X(gauss_sigma).

   @param[in] A The image.
   @return The FWHM diameter (pixel), or 0 if A is NULL.
*/
R X(fwhm)(X(mat)* A){
	if(!A) return 0;
	R hm=0.5*X(max)(A);
	long n=0;
	OMP_SIMD_R(reduction(+:n))
	for(long ix=0; ix<A->nx*A->ny; ix++){
		if(REAL(P(A, ix))>=hm){
			n++;
		}
	}
	return sqrt(n*4/M_PI);
}

/**
   Fit a 2D Gaussian function by computing the first and second moments of A. All outputs are set to 0
   if the sum of A is 0.

   @param[out] mr Equivalent sigma, (det(C))^(1/4) where C is the second moment matrix (pixel). Can be NULL.
   @param[out] ma Sigma along the major axis (pixel). Can be NULL.
   @param[out] mb Sigma along the minor axis (pixel). Can be NULL.
   @param[out] angle Angle of the major axis (radian). Can be NULL.
   @param[out] amp Peak amplitude of the fitted Gaussian. Can be NULL, but is required if img is set.
   @param[in,out] img Returns the fitted Gaussian image. Allocated if NULL. Can be NULL to skip.
   @param[in] A The irradiance (intensity).
   @param[in] thres The threshold relative to peak for computing the second moments.
*/
void X(gauss_fit)(R* mr, R* ma, R* mb, R* angle, R* amp, X(mat)**img, X(mat)*A, R thres){
	long ny=NY(A);
	long nx=NX(A);

	R sum=0, sumx=0, sumy=0, Amax=0;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			R Ai=REAL(P(A, ix, iy));
			sum+=Ai;
			sumx+=Ai*ix;
			sumy+=Ai*iy;
			if(Amax<Ai) Amax=Ai;
		}
	}
	if(!sum){
		if(mr) *mr=0;
		if(ma) *ma=0;
		if(mb) *mb=0;
		if(angle) *angle=0;
		if(amp) *amp=0;
		return;
	}
	thres*=Amax;
	sumx/=sum;
	sumy/=sum;
	R sumx2=0, sumy2=0, sumxy=0;
	for(long iy=0; iy<ny; iy++){
		R y=(R)iy-sumy;
		OMP_SIMD_R(reduction(+:sumx2, sumy2, sumxy))
		for(long ix=0; ix<nx; ix++){
			R x=(R)ix-sumx;
			R Ai=REAL(P(A, ix, iy));
			if (Ai>thres){
				sumx2+=Ai*x*x;
				sumy2+=Ai*y*y;
				sumxy+=Ai*x*y;
			}
		}
	}
	sumx2/=sum;
	sumy2/=sum;
	sumxy/=sum;
	R tmp=sqrt((sumx2-sumy2)*(sumx2-sumy2)+4*(sumxy*sumxy));
	//equivalent radius with same area is sqrt(ma*mb)
	const R det=sumx2*sumy2-sumxy*sumxy;
	if(mr) *mr=pow(det, 0.25);
	if(ma) *ma=sqrt(0.5*(sumx2+sumy2+tmp));
	if(mb) *mb=sqrt(0.5*(sumx2+sumy2-tmp));
	if(angle) *angle=0.5*atan2(2*sumxy,(sumx2-sumy2));
	const R amp0=sum/(2.*M_PI*sqrt(det));
	if(amp) *amp=amp0;
	if(img){
		if(X(init)(img, nx, ny)) return;
		// inverse covariance
        R inv00 =  sumy2 / det;
        R inv11 =  sumx2 / det;
        R inv01 = -sumxy / det;

        for(long iy=0; iy<ny; iy++){
            R y = (iy - sumy);
            for(long ix=0; ix<nx; ix++){
                R x = (ix - sumx);
                R Q = inv00*x*x + 2*inv01*x*y + inv11*y*y;
                P(*img, ix, iy) = amp0 * exp(-0.5 * Q);
            }
        }

	}
}
/**
   Determine the sigma_r of 2d gauss function using X(gauss_fit).

   @param[in] A The image.
   @param[in] thres The threshold relative to peak.
   @return The equivalent sigma (pixel).
*/
R X(gauss_sigma)(X(mat)*A, R thres){
	R mr;
	X(gauss_fit)(&mr, 0, 0, 0, 0, NULL, A, thres);
	return mr;
}
/**
   Use gauss fit to compute fwhm.

   @param[in] A The image.
   @return The FWHM (pixel) as 2*sqrt(2*log(2)) times the equivalent sigma.
*/
R X(fwhm_gauss)(X(mat)* A){
	return X(gauss_sigma)(A, 0.01)*2*1.17741;
}
#ifndef COMP_COMPLEX
struct enc_t{
	X(mat)* enc; /**<Output*/
	const X(mat)* dvec;/**<Radius wanted*/
	X(mat)* phat; /**<processed image.*/
	int type;
};

/**
   Thread worker for X(enc). Computes the result for a range of entries of dvec.

   @param[in,out] pdata The thread data with range of dvec indices and struct enc_t, whose enc is accumulated into.
   @return Always NULL.
*/
static void* X(enc_thread)(thread_t* pdata){
	struct enc_t* data=(struct enc_t*)pdata->data;
	const X(mat)* dvec=data->dvec;
	X(mat)* restrict enc=data->enc;
	const X(mat)* restrict ppsf=data->phat;
	int type=data->type;
	const R* restrict dr=P(dvec);
	const long ncomp2=data->phat->nx;
	const long ncomp=ncomp2/2;
	const R dk=1./ncomp2;
	const R pi2=2*M_PI;
	if(type==0 || type==2){
		X(mat)* ksinc=X(new)(PN(dvec), ncomp2);
		/*Cache the data. */
		for(long iy=0; iy<ncomp2; iy++){
			R ky=(iy<ncomp?iy:iy-ncomp2)*dk;
			OMP_SIMD()
			for(long ir=pdata->start; ir<pdata->end; ir++){
				P(ksinc, ir, iy)=sinc(ky*dr[ir])*dr[ir];
			}
		}
		for(long iy=0; iy<ncomp2; iy++){
			for(long ix=0; ix<ncomp2; ix++){
				OMP_SIMD()
				for(long ir=pdata->start; ir<pdata->end; ir++){
					R s=type==2?P(ksinc, ir, ix):P(ksinc, ir, iy)*P(ksinc, ir, ix);
					P(enc, ir)+=s*P(ppsf, ix, iy);
				}
			}
		}
		X(free)(ksinc);
	} else{
		for(long iy=0; iy<ncomp2; iy++){
			R ky=(iy<ncomp?iy:iy-ncomp2)*dk;
			for(long ix=0; ix<ncomp2; ix++){
				R kx=(ix<ncomp?ix:ix-ncomp2)*dk;
				switch(type){
				case -1: {/*azimuthal average. dr is radius */
					R k=sqrt(kx*kx+ky*ky);
					//OMP_SIMD()
					for(long ir=pdata->start; ir<pdata->end; ir++){
						R s=j0(k*pi2*dr[ir]);
						P(enc, ir)+=s*P(ppsf, ix, iy);
					}
				} break;
				case 0:/*Ensquared energy in a pixel*/
					break;
				case 1: {/*Encircled energy. dr is diameter */
					R k=sqrt(kx*kx+ky*ky);
					for(long ir=pdata->start; ir<pdata->end; ir++){
						const R r=dr[ir]*0.5;
						const R tmp=k*pi2*r;
						R s=j1(tmp)*r/k;
						if(!ix&&!iy) s=M_PI*r*r;/*special case when k=0. J1(x)/x=1/2 when x=0. */
						P(enc, ir)+=s*P(ppsf, ix, iy);
					}
				} break;
				case 2:/*Enslitted energe in a slit. */
					break;
				default:
					error("Not implemented\n");
					return NULL;
				}
			}
		}
	}
	return NULL;
}
/**
   Compute the enclosed energy or azimuthal average of a.

   @param[in] psf The input array.
   @param[in] dvec The diameter for enclosed energy, or radius for azimuthal average.
   @param[in] type The type. -1: azimuthal average, 0: within a square, 1: within a circle, 2: within a slit.
   @param[in] nthread Number of CPU threads to use. NTHREAD if not positive.
   @return The result for each entry of dvec. Newly allocated; the caller frees it.
*/
X(mat)* X(enc)(const X(mat)* psf, const X(mat)* dvec, int type, int nthread){
	if(type<-1||type>2){
		error("Usage: type= \n-1: azimuthal average, \n0: within a square, \n1: within a circle, \n2: within a slit\n");
		return NULL;
	}
	if(nthread<=0) nthread=NTHREAD;	
	const R rmax=ceil(X(max)(dvec))+1;
	long ncomp=nextfftsize(rmax*2);
	const long ncomp_max=MAX(psf->nx, psf->ny);
	X(mat)* psfc;
	if(ncomp_max>ncomp){//reduce the size of array
		psfc=X(new)(ncomp, ncomp);
		X(embed)(psfc, psf, 0);
		//info("1ncomp=%ld, %ld, nref=%d, %d, %p %p\n", ncomp, ncomp_max, mem_nref(psf->mem), mem_nref(psfc->mem), P(psf), P(psfc));
	} else{
		psfc=X(ref)(psf);
		//info("2ncomp=%ld, %ld, nref=%d\n", ncomp, ncomp_max, mem_nref(psf->mem));
		ncomp=ncomp_max;
	}
	long ncomp2=ncomp*2;//to avoid wrapping
	XC(mat)* psf2=XC(new)(ncomp2, ncomp2);
	XC(embedd)(psf2, psfc, 0);
	X(free)(psfc);
	XC(fftshift)(psf2);
	if(type==2){
		XC(fft2partial)(psf2, 0, -1);
	}else{
		XC(fft2)(psf2, -1);
	}
	X(mat)* phat=NULL;
	XC(real2d)(&phat, 0, psf2, (type==2)?(1./ncomp2):(1./ncomp2/ncomp2));
	XC(free)(psf2);
	X(mat)* enc=X(new)(PN(dvec), 1);
	struct enc_t data={enc, dvec, phat, type};
	thread_t* tdata=thread_prep(0, PN(dvec), nthread, X(enc_thread), &data);
	CALL_THREAD(tdata, 0);
	free(tdata);
	X(free)(phat);
	return enc;
}

#endif

/**
   Trapzoidal integration of each column of y. If x is NULL, unit spacing is assumed.

   @param[in] x The abscissa (one column or one per column of y). Can be NULL.
   @param[in] y The function values.
   @return The sum of integrations over all columns.
*/
T X(trapz)(const X(mat)* restrict x, const X(mat)* restrict y){
	if(!y) return 0;
	if(x&&x->nx!=y->nx){
		error("First dimension of x must match y\n");
		return 0;
	}
	T out=0;
	for(long icol=0; icol<y->ny; icol++){
		T* restrict py=PCOL(y, icol);
		T* restrict px=0;
		if(x){
			if(x->ny==y->ny){
				px=PCOL(x, icol);
			} else{
				px=P(x);
			}
		}
		T ans=0;
		if(px){
			OMP_SIMD_R(reduction(+:ans))
			for(long i=0; i<y->nx-1; i++){
			//notice use of abs here.
				ans+=ABS(px[i+1]-px[i])*(py[i+1]+py[i]);
			}
		} else{
			OMP_SIMD_R(reduction(+:ans))
			for(long i=0; i<y->nx; i++){
				ans+=py[i];
			}
			ans=(ans*2.-py[0]-py[y->nx-1]);
		}
		out+=ans*0.5;
	}
	return out;
}

/**
   duplicate a X(cell) object.
   Not implemented for lmat or imat as X(cellcp) is defined in spmm.c

   @param[in] in The cell to duplicate.
   @return The deep copy. Newly allocated; the caller frees it.
*/
X(cell)* X(celldup)(const X(cell)* in){
	X(cell)* out=NULL;
	X(cellcp)(&out, in);
	return out;
}
/**
   concatenate two cell matrices along dimenstion 'dim'.

   @param[in] A The first cell.
   @param[in] B The second cell.
   @param[in] dim 1: along x, 2: along y.
   @return The concatenated cell. Newly allocated; the caller frees it.
*/
X(cell)* X(cellcat)(const X(cell)* A, const X(cell)* B, int dim){
	if(!A){
		if(B) return X(celldup)(B); 
		else return NULL;
	}		
	X(cell)* out=X(celldup)(A);
	X(cellcat2)(&out, B, dim);
	return out;
}
/**
   Concatenate cell B into A along dimension 'dim'.

   @param[in,out] pA The cell A, extended in place. Set to a copy of B if NULL.
   @param[in] B The cell to append. Its blocks are duplicated.
   @param[in] dim 1: along x, 2: along y.
*/
void X(cellcat2)(X(cell)** pA, const X(cell)* B, int dim){
	if(!B) return;
	if(!pA){error("A must not be NULL\n"); return;}
	if(!*pA){*pA=X(celldup)(B);	return; }
	X(cell)*A=*pA;
	if(A->m){
		warning("cellcat2 require A->m to be unset. freed.");
		X(free)(A->m);
	}
	if(dim==1){/*along x. */
		if(A->ny!=B->ny){
			error("Mismatch: A is (%ld, %ld), B is (%ld, %ld)\n",
				A->nx, A->ny, B->nx, B->ny);
			return;
		}
		long xoff=A->nx;
		cellresize(A, NX(A)+NX(B),NY(A));
		for(long iy=0; iy<A->ny; iy++){
			for(long ix=0; ix<B->nx; ix++){
				P(A, ix+xoff, iy)=X(dup)(P(B, ix, iy));
			}
		}
	} else if(dim==2){/*along y. */
		if(A->nx!=B->nx){
			error("Mismatch. A is (%ld, %ld), B is (%ld, %ld)\n",
				A->nx, A->ny, B->nx, B->ny);
			return;
		}
		long yoff=A->ny;
		cellresize(A, NX(A), NY(A)+NY(B));
		for(long iy=0; iy<B->ny; iy++){
			for(long ix=0; ix<B->nx; ix++){
				P(A, ix, iy+yoff)=X(dup)(P(B, ix, iy));
			}
		}
	} else{
		error("Invalid dim\n");
	}
}

/**
   concatenate coresponding elements of each X(cell). They must
   have the same shape.

   @param[in] A The first cell.
   @param[in] B The second cell.
   @param[in] dim 1: along x, 2: along y.
   @return The cell of concatenated blocks. Newly allocated; the caller frees it.
*/
X(cell)* X(cellcat_each)(const X(cell)* A, const X(cell)* B, int dim){
	if(!A){
		if(!B){
			return NULL;
		} else{
			return X(celldup)(B);
		}
	} else if(!B){
		return X(celldup)(A);
	}
	if(A->nx!=B->nx||A->ny!=B->ny){
		error("Mismatch: (%ld %ld), (%ld %ld)\n", A->nx, A->ny, B->nx, B->ny);
		return NULL;
	}
	X(cell)* out=X(cellnew)(A->nx, A->ny);
	for(long ix=0; ix<A->nx*A->ny; ix++){
		P(out, ix)=X(cat)(P(A, ix), P(B, ix), dim);
	}
	return out;
}

/**
   Sum of X(norm) of each cell element.

   @param[in] A The cell.
   @return The sum of X(norm) of each block.
*/
R X(cellnorm)(const X(cell)* A){
	R out=0;
	for(int i=0; i<A->nx*A->ny; i++){
		out+=X(norm)(P(A, i));
	}
	return out;
}



/**
   Compute the dot product of two cells.

   @param[in] A The first cell, a column.
   @param[in] B The second cell, a column.
   @return The sum of X(dot) of each block, or 0 if A or B is NULL.
*/
T X(celldot)(const X(cell)* A, const X(cell)* B){
	if(!A||!B) return 0;
	if(A->nx!=B->nx||A->ny!=1||B->ny!=1){
		error("mismatch: A is %ldx%ld, B is %ldx%ld\n", A->nx, A->ny, B->nx, B->ny);
		return 0;
	}
	T out=0;
	for(int i=0; i<A->nx; i++){
		out+=X(dot)(P(A, i), P(B, i));
	}
	return out;
}

/**
   Component wise multiply of two cells.
   B=A.*B

   @param[in,out] B The cell B, modified in place.
   @param[in] A The cell A.
*/
void X(cellcwm)(X(cell)* B, const X(cell)* A){
	if(!A||!B) return;
	if(A->nx!=B->nx||A->ny!=B->ny){
		error("mismatch\n");
		return;
	}
	for(int i=0; i<A->nx*A->ny; i++){
		X(cwm)(P(B, i), P(A, i));
	}
}

/**
   Drop empty blocks (zero). Size of B is not modified.

   @param[in,out] B The cell. Blocks with all absolute values not above thres are freed.
   @param[in] thres The threshold.
*/
void X(celldropzero)(X(cell)* B, R thres){
	X(cell)* Bp=B;
	for(long iy=0; iy<B->ny; iy++){
		for(long ix=0; ix<B->nx; ix++){
			X(mat)* tmp=P(Bp, ix, iy);
			if(!tmp) continue;
			int hasnonzero=0;
			for(int ixy=0; ixy<tmp->nx*tmp->ny; ixy++){
				if(ABS(P(tmp, ixy))>thres){
					hasnonzero=1;
					break;
				}
			}
			if(!hasnonzero){
				X(free)(P(Bp, ix, iy));
				P(Bp, ix, iy)=NULL;
				/*warning("Dropped block (%ld, %ld)\n", ix, iy); */
			}
		}
	}
}

/**
   Compute the relative difference of two cells as
   sqrt(X(cellnorm)(A-B)*2/(X(cellnorm)(A)+X(cellnorm)(B))).

   @param[in] A The first cell.
   @param[in] B The second cell.
   @return The relative difference, or 0 if it is NaN.
*/
R X(celldiff)(const X(cell)* A, const X(cell)* B){
	X(cell)* C=NULL;
	X(cellcp)(&C, A);
	X(celladd)(&C, 1, B, -1);
	R d=sqrt(X(cellnorm)(C)*2/(X(cellnorm)(A)+X(cellnorm)(B)));
	X(cellfree)(C);
	return isnan(d)?0:d;
}

/**
   Clip a cell array to max at 'max', min at 'min'.

   @param[in,out] Ac The cell, clipped in place.
   @param[in] min The lower limit.
   @param[in] max The upper limit.
   @return The number of elements clipped.
*/
int X(cellclip)(X(cell)* Ac, R min, R max){
	if(!Ac||!P(Ac)) return 0;
	if(isinf(min)&&isinf(max)) return 0;
	int nclip=0;
	for(long i=0; i<Ac->nx*Ac->ny; i++){
		nclip+=X(clip)(P(Ac, i), min, max);
	}
	return nclip;
}

/**
   Raise each element of each block in the cell array to power of power.

   @param[in,out] A The cell, modified in place.
   @param[in] power The exponent.
*/
void X(cellcwpow)(X(cell)* A, R power){
	if(!A) return;
	for(long ib=0; ib<A->nx*A->ny; ib++){
		X(cwpow)(P(A, ib), power);
	}
}

/**
   2D cubic spline interpolation preparation. x is the x coordinate vector of
   the 2-d grid. y is the y coordinate vector of the 2-d grid. z is defined on the
   2-d grid.  It is upto the user to make sure that the coordinate is increasingly
   ordered and evenly spaced .

   The boundaries are handled in the same way is X(spline). i.e. replace
   \f[f^\prime(0)=(f(1)-f(-1))/2\f] by
   \f[f^\prime(0)=(f(1)-f(0))\f]
   Otehr type of boundaries are handled in the same way.

   @param[in] x The x coordinate vector.
   @param[in] y The y coordinate vector.
   @param[in] z The function values of size numel(x) by numel(y).
   @return The nx by ny cell of 4x4 coefficients. The last row and column are empty. Newly allocated; the caller frees it.
*/
X(cell)* X(bspline_prep)(X(mat)* x, X(mat)* y, X(mat)* z){
	const long nx=x->nx;
	const long ny=y->nx;
	assert(x->ny==1&&y->ny==1&&z->nx==nx&&z->ny==ny);
	X(cell)* coeff=X(cellnew)(nx, ny);
	X(cell)* pc=coeff;

	X(mat)* p=z;
	T p00, p01, p02, p03, p10, p11, p12, p13, p20, p21, p22, p23, p30, p31, p32, p33;
	for(long iy=0; iy<ny-1; iy++){
		for(long ix=0; ix<nx-1; ix++){
			if(iy==0){
				if(ix==0){
					p00=2.*(2.*P(p, ix, iy)-P(p, ix+1, iy))-(2.*P(p, ix, iy+1)-P(p, ix+1, iy+1));/*from a */
				} else{
					p00=2.*P(p, ix-1, iy)-P(p, ix-1, iy+1);/*from b */
				}
				p01=2.*P(p, ix, iy)-P(p, ix, iy+1);
				p02=2.*P(p, ix+1, iy)-P(p, ix+1, iy+1);
				if(ix==nx-2){
					p03=2.*(P(p, ix+1, iy)*2.-P(p, ix, iy))-(P(p, ix+1, iy+1)*2.-P(p, ix, iy+1));/*from n */
				} else{
					p03=2.*P(p, ix+2, iy)-P(p, ix+2, iy+1);/*from m */
				}
			} else{
				if(ix==0){
					p00=2.*P(p, ix, iy-1)-P(p, ix+1, iy-1);/*a from b */
				} else{
					p00=P(p, ix-1, iy-1);/*b */
				}
				p01=P(p, ix, iy-1);
				p02=P(p, ix+1, iy-1);
				if(ix==nx-2){
					p03=P(p, ix+1, iy-1)*2.-P(p, ix, iy-1);/*n from m */
				} else{
					p03=P(p, ix+2, iy-1);/*m */
				}
			}
			if(ix==0){
				p10=P(p, ix, iy)*2.-P(p, ix+1, iy);/*from c */
			} else{
				p10=P(p, ix-1, iy);/*c */
			}
			p11=P(p, ix, iy);
			p12=P(p, ix+1, iy);
			if(ix==nx-2){
				p13=P(p, ix+1, iy)*2.-P(p, ix, iy);/*from d */
			} else{
				p13=P(p, ix+2, iy);/*d */
			}
			if(ix==0){
				p20=P(p, ix, iy+1)*2.-P(p, ix+1, iy+1);/*from e */
			} else{
				p20=P(p, ix-1, iy+1);/*e */
			}
			p21=P(p, ix, iy+1);
			p22=P(p, ix+1, iy+1);
			if(ix==nx-2){
				p23=P(p, ix+1, iy+1)*2.-P(p, ix, iy+1);/*from f */
			} else{
				p23=P(p, ix+2, iy+1);/*f */
			}
			if(iy==ny-2){
				if(ix==0){
					p30=2.*(P(p, ix, iy+1)*2.-P(p, ix+1, iy+1))-(P(p, ix, iy)*2.-P(p, ix+1, iy));/*from h */
				} else{
					p30=2.*P(p, ix-1, iy+1)-P(p, ix-1, iy);/*from g */
				}
				p31=2.*P(p, ix, iy+1)-P(p, ix, iy);
				p32=2.*P(p, ix+1, iy+1)-P(p, ix+1, iy);
				if(ix==nx-2){
					p33=2.*(2.*P(p, ix+1, iy+1)-P(p, ix, iy+1))-(2.*P(p, ix+1, iy)-P(p, ix, iy));/*from j */
				} else{
					p33=2.*P(p, ix+2, iy+1)-P(p, ix+2, iy);/*from i */
				}
			} else{
				if(ix==0){
					p30=P(p, ix, iy+2)*2.-P(p, ix+1, iy+2);/*h from g */
				} else{
					p30=P(p, ix-1, iy+2);/*g */
				}
				p31=P(p, ix, iy+2);
				p32=P(p, ix+1, iy+2);
				if(ix==nx-2){
					p33=2.*P(p, ix+1, iy+2)-P(p, ix, iy+2);/*j from i */
				} else{
					p33=P(p, ix+2, iy+2);/*i */
				}
			}
			P(pc, ix, iy)=X(new)(4, 4);
			X(mat*) ppc=P(pc, ix, iy);
			P(ppc, 0, 0)=p11;
			P(ppc, 1, 0)=-.5*p10+.5*p12;
			P(ppc, 2, 0)=p10-2.5*p11+2.*p12-.5*p13;
			P(ppc, 3, 0)=-.5*p10+1.5*p11-1.5*p12+.5*p13;
			P(ppc, 0, 1)=-.5*p01+.5*p21;
			P(ppc, 1, 1)=.25*p00-.25*p02-.25*p20+.25*p22;
			P(ppc, 2, 1)=-.5*p00+1.25*p01-p02+.25*p03+.5*p20-1.25*p21+p22-.25*p23;
			P(ppc, 3, 1)=.25*p00-.75*p01+.75*p02-.25*p03-.25*p20+.75*p21-.75*p22+.25*p23;
			P(ppc, 0, 2)=p01-2.5*p11+2.*p21-.5*p31;
			P(ppc, 1, 2)=-.5*p00+.5*p02+1.25*p10-1.25*p12-p20+p22+.25*p30-.25*p32;
			P(ppc, 2, 2)=p00-2.5*p01+2.*p02-.5*p03-2.5*p10+6.25*p11-5.*p12+1.25*p13+2.*p20-5.*p21+4.*p22-p23-.5*p30+1.25*p31-p32+.25*p33;
			P(ppc, 3, 2)=-.5*p00+1.5*p01-1.5*p02+.5*p03+1.25*p10-3.75*p11+3.75*p12-1.25*p13-p20+3.*p21-3.*p22+p23+.25*p30-.75*p31+.75*p32-.25*p33;
			P(ppc, 0, 3)=-.5*p01+1.5*p11-1.5*p21+.5*p31;
			P(ppc, 1, 3)=.25*p00-.25*p02-.75*p10+.75*p12+.75*p20-.75*p22-.25*p30+.25*p32;
			P(ppc, 2, 3)=-.5*p00+1.25*p01-p02+.25*p03+1.5*p10-3.75*p11+3.*p12-.75*p13-1.5*p20+3.75*p21-3.*p22+.75*p23+.5*p30-1.25*p31+p32-.25*p33;
			P(ppc, 3, 3)=.25*p00-.75*p01+.75*p02-.25*p03-.75*p10+2.25*p11-2.25*p12+.75*p13+.75*p20-2.25*p21+2.25*p22-.75*p23-.25*p30+.75*p31-.75*p32+.25*p33;

		}
	}
	return coeff;
}

/**
   Evaluate 2D cubic spline at location defined 2-d arrays by xnew, ynew.

   @param[in] coeff The coefficients from X(bspline_prep).
   @param[in] x The x coordinate vector of the grid.
   @param[in] y The y coordinate vector of the grid.
   @param[in] xnew The x coordinate of the new points.
   @param[in] ynew The y coordinate of the new points. Must match xnew in size.
   @return The values with the shape of xnew. Newly allocated; the caller frees it.
*/
X(mat)* X(bspline_eval)(X(cell)* coeff, X(mat)* x, X(mat)* y, X(mat)* xnew, X(mat)* ynew){
	const long nx=x->nx;
	const long ny=y->nx;
	T xmin=P(x, 0);
	T ymin=P(y, 0);
	T xsep1=(R)(nx-1)/(P(x, nx-1)-xmin);
	T ysep1=(R)(ny-1)/(P(y, ny-1)-ymin);
	assert(xnew->nx==ynew->nx&&xnew->ny==ynew->ny);
	X(mat)* zz=X(new)(xnew->nx, xnew->ny);
	X(cell)* pc=coeff;
	for(long ix=0; ix<xnew->nx*xnew->ny; ix++){
		R xm=REAL((P(xnew, ix)-xmin)*xsep1);
		long xmf=floor(xm);
		if(xmf<0) xmf=0;
		if(xmf>nx-2) xmf=nx-2;
		xm=xm-xmf;

		R ym=REAL((P(ynew, ix)-ymin)*ysep1);
		long ymf=floor(ym);
		if(ymf<0) ymf=0;
		if(ymf>ny-2) ymf=ny-2;
		ym=ym-ymf;

		T xm2=xm*xm;
		T xm3=xm2*xm;
		T ym2=ym*ym;
		T ym3=ym2*ym;
		X(mat*)ppc=P(pc, xmf, ymf);
		P(zz, ix)=P(ppc, 0, 0)+P(ppc, 1, 0)*xm+P(ppc, 2, 0)*xm2+P(ppc, 3, 0)*xm3+
			P(ppc, 0, 1)*ym+P(ppc, 1, 1)*ym*xm+P(ppc, 2, 1)*ym*xm2+P(ppc, 3, 1)*ym*xm3+
			P(ppc, 0, 2)*ym2+P(ppc, 1, 2)*ym2*xm+P(ppc, 2, 2)*ym2*xm2+P(ppc, 3, 2)*ym2*xm3+
			P(ppc, 0, 3)*ym3+P(ppc, 1, 3)*ym3*xm+P(ppc, 2, 3)*ym3*xm2+P(ppc, 3, 3)*ym3*xm3;

	}
	return zz;
}
