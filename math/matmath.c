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

#include "../sys/sys.h"
#include "random.h"
#include "type.h"
#include "mathdef.h"
#include "defs.h"

/*
Math routines that are not included in mat.c
*/
/**
   scale each element of A by w.
*/
void X(scale)(X(mat)* A, R w){
	if(!A) return;
	if(!check_mat(A)){
		warning("Input is not valid");
		return;
	}
	if(w==(T)0){
		memset(P(A), 0, sizeof(T)*A->nx*A->ny);
	} else{
		OMP_SIMD()
		for(int i=0; i<A->nx*A->ny; i++){
			P(A, i)*=w;
		}
	}
}

/**
 * Check for NaN in elements
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
   compute the norm(2) of A
*/
R X(norm)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	return sqrt(X(sumsq)(A)/A->ny);
}
/**
   compute the standard deviation of A
*/
R X(std)(const X(mat)* A){
	if(!check_mat(A)) return 0;
	long N=A->nx*A->ny;
	T sum=X(sum)(A);
	R var=(X(sumsq)(A)-ABS2(sum)/N)/(N-1);
	return sqrt(var);
}
/**
   Fill A with random uniform numbers between [0, 1]*max
*/
void X(randu)(X(mat)* A, const T max, rand_t* rstat){
	if(!check_mat(A)){
		warning("Input is not valid");
		return;
	}
	for(int i=0; i<A->nx*A->ny; i++){
		P(A, i)=RANDU(rstat)*max;
	}
}
/**
   Fill A with random normal distribution numbers with
   standard deviation of sigma.
*/
void X(randn)(X(mat)* A, const T sigma, rand_t* rstat){
	if(!check_mat(A)){
		warning("Input is not valid");
		return;
	}
	for(int i=0; i<A->nx*A->ny; i++){
		P(A, i)=RANDN(rstat)*sigma;
	}
}

/**
   compute the inner product of A and B. (inner product)
*/

T X(inn)(const X(mat)* A, const X(mat)* B){
	if(!check_match(A, B)) return -1;
	TD out=0;
	OMP_SIMD(reduction(+:out))
	for(int i=0; i<A->nx*A->ny; i++){
		out+=P(A, i)*P(B, i);
	}
	if(isnan(REAL(out))){
		writebin(A, "inn_A");
		writebin(B, "inn_B");
		warning("NaN found\n");
	}
	return out;
}

/**
   compute weighted dot product a'*(w*b)
*/
T X(wdot)(const T* a, const X(mat)* w, const T* b){
	if(!check_mat(w)||!a||!b) return 0;
	TD res=0;
	OMP_SIMD(reduction(+:res) collapse(2))
	for(int j=0; j<w->ny; j++){
		for(int i=0; i<w->nx; i++){
			res+=P(w, i, j)*a[i]*b[j];
		}
	}
	if(isnan(REAL(res))){
		writebin(w, "wdot_w");
		writearr("wdot_a", 1, sizeof(real), M_T, NULL, a, w->nx, w->ny);
		writearr("wdot_b", 1, sizeof(real), M_T, NULL, b, w->nx, w->ny);
		error("NaN found.\n");
	}
	return res;
}

/**
   Compute component wise multiply A=A.*B
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


/**
   component-wise multiply of three matrices.
   A=A.*W.*(B1*wt1+B2*wt2)
*/
void X(cwm3)(X(mat)* restrict A, const X(mat)* restrict W,
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
			for(long i=0; i<B1->nx*B1->ny; i++){
				P(A, i)*=P(W, i)*(P(B1, i)*wt1+P(B2, i)*wt2);
			}
		} else if(has_b1){
			for(long i=0; i<B1->nx*B1->ny; i++){
				P(A, i)*=P(W, i)*P(B1, i)*wt1;
			}
		} else if(has_b2){
			for(long i=0; i<B2->nx*B2->ny; i++){
				P(A, i)*=P(W, i)*P(B2, i)*wt2;
			}
		}
	}
}

/**
   Component-wise multiply each column of A with B
   A(:,i)=A(:,i).*B;
*/
void X(cwmcol)(X(mat)* restrict A, const X(mat)* restrict B){
	if(check_mat(A, B)&&A->nx==B->nx&&B->ny==1){
		T* B1=P(B);
		for(long iy=0; iy<A->ny; iy++){
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
void X(cwmcol2)(X(mat)* restrict A,
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
			for(long iy=0; iy<A->ny; iy++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b1){
		for(long ix=0; ix<A->nx; ix++){
			T junk=P(B1, ix)*wt1;
			for(long iy=0; iy<A->ny; iy++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b2){
		X(cwmcol2)(A, B2, wt2, B1, wt1);
	}
}

/**
   component wise multiply of 2d complex matrix A,W and 1d vector B.
   A(:,i)=A(:,i).*W(:,i).*B;
*/
void X(cwm3col)(X(mat)* restrict A, const X(mat)* restrict W,
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
				for(long ix=0; ix<A->nx; ix++){
					P(A, ix, iy)=P(A, ix, iy)*P(W, ix, iy)*(B1p[ix]*wt1+B2p[ix]*wt2);
				}
			}
		} else if(has_b1){
			T* B1p=P(B1);

			for(long iy=0; iy<A->ny; iy++){
				for(long ix=0; ix<A->nx; ix++){
					P(A, ix, iy)=P(A, ix, iy)*P(W, ix, iy)*B1p[ix]*wt1;
				}
			}
		} else if(has_b2){
			X(cwm3col)(A, W, B2, wt2, B1, wt1);
		}
	}
}
/**
   Component wise multiply each row of A with B.
   A(i,:)=A(i,:)*B
*/
void X(cwmrow)(X(mat)* restrict A, const X(mat)* restrict B){
	if(!A||!B||!check(A->ny==B->nx&&B->ny==1)){
		warning("Input is not valid\n");
		return;
	}
	T* B1=P(B);
	for(long iy=0; iy<A->ny; iy++){
		T junk=B1[iy];
		for(long ix=0; ix<A->nx; ix++){
			P(A, ix, iy)*=junk;
		}
	}
}

/**
   component-wise multiply of rows of A with combination of B1 and B2:
   A(i,:)=A(i,:).*(B1*wt1+B2*wt2);
*/
void X(cwmrow2)(X(mat)* restrict A,
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
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b1){
		for(long iy=0; iy<A->ny; iy++){
			T junk=P(B1, iy)*wt1;
			for(long ix=0; ix<A->nx; ix++){
				P(A, ix, iy)*=junk;
			}
		}
	} else if(has_b2){
		X(cwmrow2)(A, B2, wt2, B1, wt1);
	}
}

/**
   Component wise division B=B./A. 0/0 is replace by 'value';
*/
void X(cwdiv)(X(mat)* B, const X(mat)* A, T value){
	if(!check_match(A, B)){
		warning("Input does not match\n");
		return;
	}
	for(int i=0; i<A->nx*A->ny; i++){
		P(B, i)/=P(A, i);
		if(isnan(REAL(P(B, i)))) P(B, i)=value;
	}
}
/**
   multiply a X(mat) matrix with a vector and accumulate to y:
   y+=A*x*alpha
*/
void X(mulvec)(T* restrict y, const X(mat)* restrict A,
	const T* restrict x, const T alpha){
	if(!check(y&&x&&A)){
		warning("input is not valid\n");
		return;
	}
	for(int iy=0; iy<A->ny; iy++){
		for(int ix=0; ix<A->nx; ix++){
			y[ix]+=(x[iy]*alpha)*P(A, ix, iy);
		}
	}
}

/**
   T matrix vector multiply optimized for just
   three values.  y=A*x;
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
   compute (A'*W*A); where diag(W)=wt
*/
X(mat)* X(mcc)(const X(mat)* A, const X(mat)* wt){
	if(!check_mat(A)||!check(A->nx==wt->nx&&wt->ny==1)) return NULL;
	int nmod=A->ny;
	int nsa2=A->nx;
	X(mat)* ata=X(new)(nmod, nmod);;
	for(int imod=0; imod<nmod; imod++){
		for(int jmod=imod; jmod<nmod; jmod++){
			T tmp=0;
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
   compute (A*W*A'); where diag(W)=wt
*/
X(mat)* X(tmcc)(const X(mat)* A, const X(mat)* wt){
	if(!check_mat(A)||!check(A->ny==wt->nx&&wt->ny==1)) return NULL;
	int nmod=A->nx;
	int nsa2=A->ny;
	X(mat)* ata=X(new)(nmod, nmod);;
	for(int imod=0; imod<nmod; imod++){
		for(int jmod=imod; jmod<nmod; jmod++){
			T tmp=0;
			for(int k=0; k<nsa2; k++){
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
   compute the relative difference betwee two vectors.
   sqrt(||A-B||/||A||) using norm2. for debugging purpose.
*/
T X(diff)(const X(mat)* A, const X(mat)* B){
	if(!check_match(A, B)) return -1;
	T d=sqrt(X(sumdiffsq)(A, B)*2/(X(norm)(A)+X(norm)(B)));
	return isnan(REAL(d))?-1:d;
}
/**
   Generate a new gray pixel map based on bilinear influence functions used in
   mkw.  It creates slightly larger map than an filled circle.  The Center and
   Radius cx,cy,r are in unit of meter, The sampling dx,dy specify spacing of
   the points in meter along x or y dimension. Each full point received value of
   'val'
*/
void X(circle)(X(mat)* A, R cx, R cy, R dx, R dy, R r, T val){
	if(!check_mat(A)||r<0){
		warning("Input is not valid");
		return;
	}
	const int nres=10;
	const R res=(R)(2./nres);
	const R res2=(R)(res*res);
	R resm=(R)((nres-1)*0.5);
	R r2=r*r;
	R r2l=r>1.5?((r-1.5)*(r-1.5)):0;
	R r2u=(r+2.5)*(r+2.5);
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
					R wty=1.-fabs(iy-iiy);
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
}

/**
   Mark valid grid points. If any direct neighbor of a point is within r, make
   the point valid. Parameters are the same as X(circle).
*/
void X(circle_symbolic)(X(mat)* A, R cx, R cy, R dx, R dy, R r){
	if(!check_mat(A)){
		warning("A is not valid");
		return;
	}
	R r2=r*r;
	R r2l=(r-1.5)*(r-1.5);//lower limit
	R r2u=(r+2.5)*(r+2.5);//upper limit
	for(int iy=0; iy<A->ny; iy++){
		R r2y=(iy*dy-cy)*(iy*dy-cy);
		for(int ix=0; ix<A->nx; ix++){
			R r2r=(ix*dx-cx)*(ix*dx-cx)+r2y;
			if(r2r<r2l){
				P(A, ix, iy)=1;
			} else if(r2r<r2u){
				for(R jy=-1; jy<=1; jy++){
					R iiy=iy+jy;
					R rr2y=(iiy*dy-cy)*(iiy*dy-cy);
					for(R jx=-1; jx<=1; jx++){
						R iix=ix+jx;
						R rr2r=(iix*dx-cx)*(iix*dx-cx)+rr2y;
						if(rr2r<=r2){
							P(A, ix, iy)=1;
							continue;
						}
					}
				}
			}
		}
	}
}

/**
   rotate the column vectors CCW, equivalent as rotate coordinate theta CW.
   A is nx2 while A(:,1) is x, A(:,2) is y.
*/
void X(rotvec)(X(mat)* A, const R theta){
	if(!check_mat(A)||A->ny!=2){
		warning("A is not valid");
		return;
	}
	const R ctheta=cos(theta);
	const R stheta=sin(theta);
	//OMP_SIMD()
	for(int i=0; i<A->nx; i++){
		T tmp=P(A, i, 0)*ctheta-P(A, i, 1)*stheta;
		P(A, i, 1)=P(A, i, 0)*stheta+P(A, i, 1)*ctheta;
		P(A, i, 0)=tmp;
	}
}

/**
   rotate the row vectors CCW. same as rotate coordinate theta CW.
   A is 2xn A(:,1) is x, A(:,2) is y.
*/
void X(rotvect)(X(mat)* A, const R theta){
	if(!check_mat(A)||A->nx!=2){
		warning("A has wrong dimension\n");
	}
	const R ctheta=cos(theta);
	const R stheta=sin(theta);
	OMP_SIMD()
	for(int i=0; i<A->ny; i++){
		T tmp=P(A, 0, i)*ctheta-P(A, 1, i)*stheta;
		P(A, 1, i)=P(A, 0, i)*stheta+P(A, 1, i)*ctheta;
		P(A, 0, i)=tmp;
	}
}

/**
   rotate a 2x2 covariance matrix A by theta CCW
   (coordinate rotate -theta CCW) or from ra to xy
   coordinate.  R*A*R';
*/
void X(rotvecnn)(X(mat)** B0, const X(mat)* A, R theta){
	if(!check_mat(A)||!check(A->nx==2&&A->ny==2)){
		warning("A has wrong dimension\n");
		return;
	}
	X(init)(B0, 2, 2);
	X(mat)* B=*B0;
	assert(B->nx==2&&B->ny==2);
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
   Compute the correlation matrix.
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
			OMP_SIMD(reduction(+:tmp) collapse(2))
			for(long iy1=sy1; iy1<nny; iy1++){
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
   all length are given in terms of pixel.
*/
void X(cog)(R* grad, const X(mat)* im, R offsetx, R offsety,
	R thres, R bkgrnd, R flux){
	R sum=0, sumx=0, sumy=0;
	R iI;
	OMP_SIMD(reduction(+:sum,sumx,sumy) collapse(2))
	for(int iy=0; iy<im->ny; iy++){
		for(int ix=0; ix<im->nx; ix++){
			iI=REAL(P(im, ix, iy))-bkgrnd;
			if(iI>thres){
				sum+=iI;
				sumx+=iI*ix;
				sumy+=iI*iy;
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
*/
void X(shift2center)(X(mat)* A, R offsetx, R offsety){
	R grad[2];
	R Amax=X(max)(A);
	X(cog)(grad, A, offsetx, offsety, Amax*0.1, Amax*0.2, 0);
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
		X(cog)(grad, A, offsetx, offsety, Amax*0.1, Amax*0.2, 0);
		/*dbg("After shift, residual grad is %g %g\n",grad[0],grad[1]); */
		XC(free)(B);
	}
}


/**
   OrthNormalize column vector in Mod, with weighting from vector amp.
   <Mod|wt|Mod> is equal to sum(wt).
   2010-07-21: Bug found: The result was not orthonormal. cause: nonvalid was not initialized to 0.
*/
void X(gramschmidt)(X(mat)*restrict Mod, R* restrict amp){
	const int nmod=Mod->ny;
	const long nx=Mod->nx;
	R wtsum=(R)nx;
	if(amp){
		wtsum=XR(vecsum)(amp, nx);
	}
	int nonvalid[nmod];
	memset(nonvalid, 0, sizeof(int)*nmod);
	for(int imod=0; imod<nmod; imod++){
		if(nmod>10){
			info("Gramschmidt: %d of %d\n", imod, nmod);
		}
		if(imod>0){/*orthogonalize */
			/*compute dot product. */
			for(int jmod=0; jmod<imod; jmod++){
				if(nonvalid[jmod]) continue;
				const T cross=-X(vecdot)(PCOL(Mod, imod), PCOL(Mod, jmod), amp, nx)/wtsum;
				OMP_SIMD()
				for(long ix=0; ix<nx; ix++){
					P(Mod, ix, imod)+=cross*P(Mod, ix, jmod);
				}
			}
		}
		/*normalize*/
		R norm=sqrt(REAL(X(vecdot)(PCOL(Mod, imod), PCOL(Mod, imod), amp, nx)/wtsum));
		if(fabs(norm)>1.e-15){
			norm=1./norm;
			OMP_SIMD()
			for(long ix=0; ix<nx; ix++){
				P(Mod, ix, imod)*=norm;
			}
		} else{
			nonvalid[imod]=1;
			warning("Column %d is not independent on other columns\n", imod);
		}
	}
}


/**
   Limit numbers in A to within [min, max]. used for DM clipping.
*/
int X(clip)(X(mat)* restrict A, R min, R max){
	if(!A) return 0;
	if(!isfinite(min)&&!isfinite(max)) return 0;
	if(max<=min){
		error("upper light should be larger than lower limit\n");
	}
	T* restrict Ap=P(A);
	int nclip=0;
	OMP_SIMD(reduction(+:nclip))
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
   A=A*B, where diag(B)=s
*/
void X(muldiag)(X(mat)* restrict A, const X(mat)* restrict s){
	assert(A->ny==s->nx&&s->ny==1);
	X(mat)* pA=A;
	const T* restrict ps=P(s);
	OMP_SIMD(collapse(2))
	for(long iy=0; iy<A->ny; iy++){
		for(long ix=0; ix<A->nx; ix++){
			P(pA, ix, iy)*=ps[iy];
		}
	}
}

/**
   Raise all elements to power 'power'
*/
void X(cwpow)(X(mat)* A, R power){
	if(!A) return;
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=pow(P(A, i), power);
	}
}

/**
   compute exponential of all elements after scaling by alpha
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
*/
void X(cwpow_thres)(X(mat)* restrict A, R power, R thres){
	thres*=X(maxabs)(A);
	for(long i=0; i<A->nx*A->ny; i++){
		if(fabs(P(A, i))>thres){
			P(A, i)=pow(P(A, i), power);
		} else{
			P(A, i)=0;
		}
	}
}
/**
   Compute polynomial functional value inplace
*/
void X(polyval)(X(mat)* A, XR(mat)* p){
	if(p->nx==1&&p->ny>1){
		p->nx=p->ny;
		p->ny=1;
	}
	int np=p->nx;
	for(long i=0; i<A->nx*A->ny; i++){
		T tmp=0;
		OMP_SIMD(reduction(+:tmp))
		for(long ip=0; ip<np; ip++){
			tmp+=P(p, ip)*(T)pow(P(A, i), np-ip-1);
		}
		P(A, i)=tmp;
	}
}
/**
   add val to diagonal values of A.
*/
void X(addI)(X(mat)*restrict A, T val){
	if(!A) return;
	if(A->nx!=A->ny)
		warning("daddI: A is not square\n");
	long M=A->nx<A->ny?A->nx:A->ny;
	for(long i=0; i<M; i++){
		P(A, i, i)+=val;
	}
}
/**
   compute B=bc*B+ac*A
   behavior changed on 2009-11-02. if A is NULL, don't do anything.
*/
void X(add)(X(mat)** B0, T bc, const X(mat)* A, const T ac){
	if(A&&A->nx){
		if(!*B0){
			bc=0;/*no bother to accumulate. */
			X(init)(B0, A->nx, A->ny);
		}else if(PN(A) != PN(*B0)){
			error("A is %ldx%ld, B is %ldx%ld. They should match in size.\n",
				NX(A), NY(A), NX(*B0), NY(*B0));
		}
		X(mat) * restrict B=*B0;
		if(bc){
			OMP_SIMD()
			for(int i=0; i<PN(A); i++){
				P(B, i)=P(B, i)*bc+P(A, i)*ac;
			}
		} else{
			if(ac==(T)1){
				X(cp)(B0, A);
			} else{/*just assign */
				OMP_SIMD()
				for(int i=0; i<PN(A); i++){
					P(B, i)=P(A, i)*ac;
				}
			}
		}
	}
}
/**
   compute B=bc*B+ac*A. Use smallest dimension of both.
   behavior changed on 2009-11-02. if A is NULL, don't do anything.
*/
void X(add_relax)(X(mat)** B0, T bc, const X(mat)* A, const T ac){
	if(A&&A->nx&&ac){
		if(!*B0){
			bc=0;/*no bother to accumulate. */
			X(init)(B0, A->nx, A->ny);
		}
		X(mat)* restrict B=*B0;
		long nx=MIN(A->nx, B->nx);
		long ny=MIN(A->ny, B->ny);
		OMP_SIMD(collapse(2))
		for(long iy=0; iy<ny; iy++){
			for(long ix=0; ix<nx; ix++){
				P(B, ix, iy)=P(B, ix, iy)*bc+P(A, ix, iy)*ac;
			}
		}
	}
}
/**
   Add a scalar to matrix
*/
void X(adds)(X(mat)*restrict A, const T ac){
	if(!A||!A->nx||ac==(T)0) return;
	OMP_SIMD()
	for(int i=0; i<A->nx*A->ny; i++){
		P(A, i)+=ac;
	}
}

/**
   Create log spaced vector.
*/
X(mat)* X(logspace)(R emin, R emax, long n){
	X(mat)* out=X(new)(n, 1);
	R esep=(emax-emin)/(n-1);
	for(long i=0; i<n; i++){
		P(out, i)=pow(10, emin+esep*i);
	}
	return out;
}

/**
   Create linearly spaced vector.
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
   Check whether xin is linearly spaced
*/
static int X(islinear)(const X(mat)* xin){
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=(P(xin, 0));
	R xmaxl=(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	if(fabs(xsep+xminl-P(xin, 1))>xsep*1.e-3){
		return 0;
	} else{
		return 1;
	}
}
/**
   Check whether xin is logrithmically spaced
*/
static int X(islog)(const X(mat)* xin){
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=log10(P(xin, 0));
	R x1=log10(P(xin, 1));
	R xmaxl=log10(P(xin, nmax1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	if(!isfinite(xsep)||fabs(xsep+xminl-x1)>xsep*1.e-3){
		return 0;
	} else{
		return 1;
	}
}
/**
   Interpolate using linear interp. xin is the coordinate of yin. xnew is the
   coordinate of the output.
*/
X(mat)* X(interp1linear)(const X(mat)* xin, const X(mat)* yin, const X(mat)* xnew, T ydefault){
	if(!X(islinear)(xin)){
		error("xin is not linearly spaced\n");
	}
	if(xin->ny!=1||xnew->ny!=1){
		error("Either xin or xnew is in wrong format\n");
	}
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=(P(xin, 0));
	R xmaxl=(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	R xsep1=1./xsep;
	X(mat)* ynew=X(new)(xnew->nx, xnew->ny);
	for(long iy=0; iy<ynew->ny; iy++){
		for(long ix=0; ix<ynew->nx; ix++){
			R xx=((P(xnew, ix))-xminl)*xsep1;
			long xxm=ifloor(xx);
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
	return ynew;
}

/**
   Interpolate using log(xin) and log(xnew)
   xin is the coordinate of yin. xnew is the coordinate of the output.
*/
X(mat)* X(interp1log)(const X(mat)* xin, const X(mat)* yin, const X(mat)* xnew, T ydefault){
	if(!X(islog)(xin)){
		error("xin is not logrithmically spaced\n");
	}
	if(xin->ny!=1||xnew->ny!=1){
		error("Either xin or xnew is in wrong format\n");
	}
	long nmax=xin->nx;
	long nmax1=nmax-1;
	R xminl=log10(P(xin, 0));
	R xmaxl=log10(P(xin, nmax-1));
	R xsep=(xmaxl-xminl)/(R)(nmax1);
	R xsep1=1./xsep;
	X(mat)* ynew=X(new)(xnew->nx, xnew->ny);
	for(long iy=0; iy<ynew->ny; iy++){
		for(long ix=0; ix<ynew->nx; ix++){
			R xx=(log10(P(xnew, ix))-xminl)*xsep1;
			long xxm=ifloor(xx);
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
	return ynew;
}
/**
   Interpolation of 1d array
 */
X(mat)* X(interp1)(const X(mat)* xin_, const X(mat)* yin_, const X(mat)* xnew, T ydefault){
	X(mat)* restrict ynew=NULL;
	X(mat) *restrict xin=X(refcols)(xin_, 0, 1);
	X(mat) *restrict yin=yin_?X(ref)(yin_):X(refcols)(xin_, 1, 1);

	if(X(islinear)(xin)){
		ynew=X(interp1linear)(xin, yin, xnew, ydefault);
	} else if(X(islog)(xin)){
		ynew=X(interp1log)(xin, yin, xnew, ydefault);
	} else{//arbitrary spacing
		if(xin->ny!=1||xnew->ny!=1){
			error("Either xin or xnew is in wrong format\n");
		}
		ynew=X(new)(xnew->nx, yin->ny);
		int curpos=0;
		for(long ix=0; ix<ynew->nx; ix++){
			int found=0;
			for(curpos=0; curpos<xin->nx-2; curpos++){
				if(P(xnew, ix)>P(xin, curpos)&&P(xnew, ix)<P(xin, curpos+1)){
					found=1;
					break;
				}
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
	X(free)(xin);
	X(free)(yin);
	return ynew;
}

#endif
#ifndef COMP_COMPLEX
/**blend B into center of A with width of overlap. The center
  (size is B->nx-overlap, B->ny-overlap) of A is replaced by
  center of B . The overlapping area is blended*/
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
   center, spaced by spacing, for n bins in total. center if at bin n/2.  */
void X(histfill)(X(mat)** out, const X(mat)* A,
	R center, R spacing, int n){
	if(!A||!P(A)) return;
	int nn=A->nx*A->ny;
	X(init)(out, n, nn);
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
   and evenly spaced .

   If the values of a function \f$f(x)\f$ and its derivative are know at x=0,
   and x=1 (normalized coordinate), then the function can be interpolated on the
   interval [0,1] using a third degree polynomial. This is called cubic
   interpolation. The formula of this polynomial can be easily derived.

   A third degree polynomial and its derivative:
   \f[
   f(x)=ax^3+bx^2+cx+d
   \f]
   \f[
   f(x)=3ax^3+2bx+c
   \f]
   The coefficients can be derived from the value and derivatives:
   \f{eqnarray*}{
   a&=&2f(0)-2f(1)+f^\prime (0)+f^\prime(1)\\
   b&=&-3f(0)+3f(1)-2f^\prime(0)-f^\prime(0)\\
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

   see www.paulinternet.nl/?page=bicubicx */
X(mat)* X(spline_prep)(X(mat)* x, X(mat)* y){
	T* px, * py;
	const long nx=x->nx;
	px=P(x);
	switch(x->ny){
	case 1:
		py=P(y);
		break;
	case 2:
		assert(y==NULL);
		py=P(x)+nx;
		break;
	default:
		py=NULL;
		error("Invalid input\n");
	}
	X(mat)* coeff=X(new)(4, nx);
	T xsep=(px[nx-1]-px[0])/(nx-1);
	R thres=fabs(xsep)*1.e-5;

	X(mat)* pc=coeff;
	T ypriv, ynext;
	OMP_SIMD()
	for(long ix=0; ix<nx-1; ix++){
		if(fabs(px[ix+1]-px[ix]-xsep)>thres){
			error("The coordinate is not evenly spaced\n");
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
   Evluate the cubic spline represented by nx5 matrix coeff, at location array xnew.
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
   Do 1D cubic spline all at once by calling X(spline_prep) and X(spline_evald)
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
*/
void X(cwlog10)(X(mat)* A){
	R ratio=1./log(10);
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=log(P(A, i))*ratio;
	}
}
/**
   Do a component wise log10 on each element of A.
*/
void X(cwlog)(X(mat)* A){
	OMP_SIMD()
	for(long i=0; i<A->nx*A->ny; i++){
		P(A, i)=log(P(A, i));
	}
}
/**
   embed a ninx*niny matrix in into A with optional rotation by -theta CCW
   (coordinate rotate theta CCW) around the fft center. Used to rotate the PSF
   from x-y to radial-azimuthal coordinate in radial format CCD. A may be bigger or smaller than B.
   \todo{
   merge this definition with cembed in cmat.c
   }
*/

void X(embed)(X(mat)* restrict A, const X(mat)* restrict B, const R theta){

	const long ninx=B->nx;
	const long niny=B->ny;
	const long noutx=A->nx;
	const long nouty=A->ny;
	memset(P(A), 0, sizeof(T)*noutx*nouty);
	if(fabs(theta)<1.e-10){/*no rotation. */
		const long skipx=(noutx-ninx-1)/2;//-1 to handle odd case
		const long skipy=(nouty-niny-1)/2;
		long ixstart=0, ixend=ninx;
		long iystart=0, iyend=niny;
		if(skipx<0){
			ixstart=-skipx;
			ixend=ninx+skipx;
		}
		if(skipy<0){
			iystart=-skipy;
			iyend=niny+skipy;
		}
		for(long iy=iystart; iy<iyend; iy++){
			T* outi=&P(A, skipx+ixstart, skipy+iy);
			T* ini=&P(B, ixstart, iy);
			memcpy(outi, ini, sizeof(T)*(ixend-ixstart));
		}
	} else{
		const R ctheta=cos(theta);
		const R stheta=sin(theta);
		R x2, y2;
		R x, y;
		long ninx2=ninx/2;
		long noutx2=noutx/2;
		long niny2=niny/2;
		long nouty2=nouty/2;
		long ix2, iy2;
		for(long iy=0; iy<nouty; iy++){
			y=(R)(iy-nouty2);
			OMP_SIMD()
			for(long ix=0; ix<noutx; ix++){
				x=(R)(ix-noutx2);
				x2=x*ctheta+y*stheta+ninx2;
				y2=-x*stheta+y*ctheta+niny2;
				if(x2>0&&x2<ninx-1&&y2>0&&y2<niny-1){
					ix2=ifloor(x2);
					iy2=ifloor(y2);
					x2=x2-ix2;
					y2=y2-iy2;
					P(A, ix, iy)=
						+P(B, ix2, iy2)*((1.-x2)*(1.-y2))
						+P(B, ix2+1, iy2)*(x2*(1.-y2))
						+P(B, ix2, iy2+1)*((1-x2)*y2)
						+P(B, ix2+1, iy2+1)*(x2*y2);
				}
			}
		}
	}
}

/**
   Calculate number of pixels having values larger than or equal to half of
   maximum and convert to diameter.

   For more accurate calculation of Gaussian beams, use the sqrt(2*log(2))*X(gauss_width).*/
R X(fwhm)(X(mat)* A){
	if(!A) return 0;
	R hm=0.5*X(max)(A);
	long n=0;
	OMP_SIMD(reduction(+:n))
	for(long ix=0; ix<A->nx*A->ny; ix++){
		if(fabs(P(A, ix))>=hm){
			n++;
		}
	}
	return sqrt(n*4/M_PI);
}

/**
 * Calculate D4s width of Gaussian 2D function
 * */
void X(gauss_fit)(
	R* mr, /**<Equivalent radius of the same area*/
	R* ma, /**<major axis*/
	R* mb, /**<minor axi*/
	R* angle, /**<angle*/
	X(mat)* A, /**<The irradiance (intensity)*/
	R thres  ///The threshold relative to peak.
	){
	long ny=NY(A);
	long nx=NX(A);

	R sum=0, sumx=0, sumy=0, Amax=0;
	//OMP_SIMD(reduction(+:sum,sumx,sumy) reduction(max:Amax) collapse(2))
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
		return;
	}
	thres*=Amax;
	sum=1./sum;
	sumx*=sum;
	sumy*=sum;
	R sumx2=0, sumy2=0, sumxy=0;
	sum=0;
	for(long iy=0; iy<ny; iy++){
		R y=(R)iy-sumy;
		OMP_SIMD(reduction(+:sumx2, sumy2, sumxy, sum))
		for(long ix=0; ix<nx; ix++){
			R x=(R)ix-sumx;
			R Ai=REAL(P(A, ix, iy));
			if (Ai>thres){
				sumx2+=Ai*x*x;
				sumy2+=Ai*y*y;
				sumxy+=Ai*x*y;
				sum+=Ai;
			}
		}
	}
	sum=1./sum;
	sumx2*=sum;
	sumy2*=sum;
	sumxy*=sum;
	R tmp=sqrt((sumx2-sumy2)*(sumx2-sumy2)+4*(sumxy*sumxy));
	//equivalent radius with same area is sqrt(ma*mb)
	if(mr) *mr=2*pow(sumx2*sumy2-sumxy*sumxy,1./4.);
	if(ma) *ma=sqrt(2*(sumx2+sumy2+tmp));
	if(mb) *mb=sqrt(2*(sumx2+sumy2-tmp));
	if(angle) *angle=0.5*atan2(2*sumxy,(sumx2-sumy2));
}
/**
 * A convenient wrapper
 * */
R X(gauss_width)(X(mat)*A, R thres){
	R mr;
	X(gauss_fit)(&mr, 0, 0, 0, A, thres);
	return mr;
}
/**
 * Use gauss fit to compute fwhm
 * */
R X(fwhm_gauss)(X(mat)* A){
	return X(gauss_width)(A, 0.01)*1.17741;
}
#ifndef COMP_COMPLEX
typedef struct{
	X(mat)* enc; /**<Output*/
	X(mat)* dvec;/**<Radius wanted*/
	X(mat)* phat; /**<processed image.*/
	int type;
}ENC_T;

static void X(enc_thread)(thread_t* pdata){
	ENC_T* data=(ENC_T*)pdata->data;
	const X(mat)* dvec=data->dvec;
	X(mat)* restrict enc=data->enc;
	const X(mat)* restrict ppsf=data->phat;
	int type=data->type;
	const R* restrict dr=P(dvec);
	const long ncomp2=data->phat->nx;
	const long ncomp=ncomp2/2;
	const R dk=1./ncomp2;
	const R pi2=2*M_PI;
	if(type==0){
		X(mat)* ksinc=X(new)(dvec->nx, ncomp2);
		X(mat)* pks=ksinc;
		/*Cache the data. */
		for(long iy=0; iy<ncomp2; iy++){
			R ky=(iy<ncomp?iy:iy-ncomp2)*dk;
			OMP_SIMD()
			for(long ir=pdata->start; ir<pdata->end; ir++){
				P(pks, ir, iy)=sinc(ky*dr[ir])*dr[ir];
			}
		}
		for(long iy=0; iy<ncomp2; iy++){
			for(long ix=0; ix<ncomp2; ix++){
				OMP_SIMD()
				for(long ir=pdata->start; ir<pdata->end; ir++){
					R s=P(pks, ir, iy)*P(pks, ir, ix);
					P(enc, ir)+=s*P(ppsf, ix, iy);
				}
			}
		}
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
				case 0:
					break;
				case 1: {/*Encircled energy. dr is diameter */
					R k=sqrt(kx*kx+ky*ky);
					//OMP_SIMD()
					for(long ir=pdata->start; ir<pdata->end; ir++){
						const R r=dr[ir]*0.5;
						const R tmp=k*pi2*r;
						R s=j1(tmp)*r/k;
						if(!ix&&!iy) s=pi2*r*r;/*special case. */
						P(enc, ir)+=s*P(ppsf, ix, iy);
					}
				} break;
				case 2:/*Enstripped energe in a slit. */
					error("To implement: Do FFT only along 1-d\n");
					break;
				default:
					error("Not implemented\n");
				}
			}
		}
	}
}
/**
   Compute the enclosed energy or azimuthal average of a.
*/
X(mat)* X(enc)(X(mat)* psf, /**<The input array*/
	X(mat)* dvec,/**<The diameter for enclosed energy, or radius for azimuthal average*/
	int type,  /**<The type. -1: azimuthal average, 0: within a square, 1: within a circle, 2: within a slit*/
	int nthread
	){
	if(type<-1||type>2){
		error("Usage: type= \n-1: azimuthal average, \n0: within a square, \n1: within a circle, \n2: within a slit\n");
	}
	R rmax=ceil(X(max)(dvec))+1;
	long ncomp;
	ncomp=nextfftsize(rmax*2);//avoid wrapping
	long ncomp_max=psf->nx>psf->ny?psf->nx:psf->ny;
	X(mat)* psfc;
	if(ncomp_max>ncomp){
		psfc=X(new)(ncomp, ncomp);
		X(embed)(psfc, psf, 0);
	} else{
		ncomp=ncomp_max;
		psfc=X(ref)(psf);
	}
	long ncomp2=ncomp*2;
	XC(mat)* psf2=XC(new)(ncomp2, ncomp2);
	//XC(fft2plan)(psf2, -1);
	XC(embedd)(psf2, psfc, 0);
	X(free)(psfc);
	XC(fftshift)(psf2);
	XC(fft2)(psf2, -1);
	X(mat)* phat=NULL;
	XC(real2d)(&phat, 0, psf2, 1);
	X(scale)(phat, pow((R)ncomp2, -2));
	XC(free)(psf2);
	X(mat)* enc=X(new)(dvec->nx, 1);
	ENC_T data={enc, dvec, phat, type};
	thread_t* tdata=thread_prep(0, dvec->nx, nthread, X(enc_thread), &data);
	CALL_THREAD(tdata, 0);
	free(tdata);
	X(free)(phat);
	return enc;
}

#endif

/**
   Trapzoidal integration
*/
T X(trapz)(const X(mat)* restrict x, const X(mat)* restrict y){
	if(!y) return 0;
	if(x&&x->nx!=y->nx){
		error("First dimension of x must match y\n");
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
			OMP_SIMD(reduction(+:ans))
			for(long i=0; i<y->nx-1; i++){
			//notice use of abs here.
				ans+=fabs(px[i+1]-px[i])*(py[i+1]+py[i]);
			}
		} else{
			OMP_SIMD(reduction(+:ans))
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
*/
X(cell)* X(celldup)(const X(cell)* in){
	X(cell)* out=NULL;
	X(cellcp)(&out, in);
	return out;
}

/**
   concatenate two cell matrices along dimenstion 'dim'.
*/
X(cell)* X(cellcat)(const X(cell)* A, const X(cell)* B, int dim){
	if(!A){
		if(!B){
			return NULL;
		} else{
			return X(celldup)(B);
		}
	} else if(!B){
		return X(celldup)(A);
	}

	X(cell)* out=NULL;

	if(dim==1){
	/*along x. */
		if(A->ny!=B->ny){
			error("Mismatch: A is (%ld, %ld), B is (%ld, %ld)\n",
				A->nx, A->ny, B->nx, B->ny);
		}
		out=X(cellnew)(A->nx+B->nx, A->ny);
		for(long iy=0; iy<A->ny; iy++){
			for(long ix=0; ix<A->nx; ix++){
				P(out, ix, iy)=X(dup)(P(A, ix, iy));
			}
			for(long ix=0; ix<B->nx; ix++){
				P(out, ix+A->nx, iy)=X(dup)(P(B, ix, iy));
			}
		}
	} else if(dim==2){
	/*along y. */
		if(A->nx!=B->nx){
			error("Mismatch. A is (%ld, %ld), B is (%ld, %ld)\n",
				A->nx, A->ny, B->nx, B->ny);
		}
		out=X(cellnew)(A->nx, A->ny+B->ny);
		for(long iy=0; iy<A->ny; iy++){
			for(long ix=0; ix<A->nx; ix++){
				P(out, ix, iy)=X(dup)(P(A, ix, iy));
			}
		}
		for(long iy=0; iy<B->ny; iy++){
			for(long ix=0; ix<B->nx; ix++){
				P(out, ix, iy+A->ny)=X(dup)(P(B, ix, iy));
			}
		}
	} else{
		error("Invalid dim\n");
	}
	return out;
}

/**
   concatenate two cell matrices along dimenstion 'dim'.
*/
void X(cellcat2)(X(cell)** A, const X(cell)* B, int dim){
	if(!B) return;
	X(cell)* Anew=X(cellcat)(*A, B, dim);
	X(cellfree)(*A);
	*A=Anew;
}
/**
   concatenate coresponding elements of each X(cell). They must
   have the same shape.
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
	}
	X(cell)* out=X(cellnew)(A->nx, A->ny);
	for(long ix=0; ix<A->nx*A->ny; ix++){
		P(out, ix)=X(cat)(P(A, ix), P(B, ix), dim);
	}
	return out;
}

/**
   compute norm2.
*/
R X(cellnorm)(const X(cell)* A){
	R out=0;
	for(int i=0; i<A->nx*A->ny; i++){
		out+=X(norm)(P(A, i));
	}
	return out;
}



/**
   Compute the inner produce of two dcell.
*/
T X(cellinn)(const X(cell)* A, const X(cell)* B){
	if(!A||!B) return 0;
	if(A->nx!=B->nx||A->ny!=1||B->ny!=1) error("mismatch\n");
	T out=0;
	for(int i=0; i<A->nx; i++){
		out+=X(inn)(P(A, i), P(B, i));
	}
	return out;
}

/**
   Component wise multiply of two dcell
   B=A.*B*alpha
*/
void X(cellcwm)(X(cell)* B, const X(cell)* A){
	if(A->nx!=B->nx||A->ny!=B->ny) error("mismatch\n");
	for(int i=0; i<A->nx*A->ny; i++){
		X(cwm)(P(B, i), P(A, i));
	}
}

/**
   drop empty blocks (zero). Size of B is not modified.
*/
void X(celldropzero)(X(cell)* B, R thres){
	X(cell)* Bp=B;
	for(long iy=0; iy<B->ny; iy++){
		for(long ix=0; ix<B->nx; ix++){
			X(mat)* tmp=P(Bp, ix, iy);
			if(!tmp) continue;
			int hasnonzero=0;
			for(int ixy=0; ixy<tmp->nx*tmp->ny; ixy++){
				if(fabs(P(tmp, ixy))>thres){
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
   compute ||A-B||/||A||
   use mean.
*/
R X(celldiff)(const X(cell)* A, const X(cell)* B){
	X(cell)* C=NULL;
	X(cellcp)(&C, A);
	X(celladd)(&C, 1, B, -1);
	R d=sqrt(X(cellnorm)(C)*2/(X(cellnorm)(C)+X(cellnorm)(B)));
	return isnan(d)?0:d;
}

/**
   clip a X(cell) array to max at 'max', min at 'min'
*/
int X(cellclip)(X(cell)* Ac, R min, R max){
	if(!Ac||!P(Ac)) return 0;
	if(!isfinite(min)&&!isfinite(max)) return 0;
	int nclip=0;
	for(long i=0; i<Ac->nx*Ac->ny; i++){
		nclip+=X(clip)(P(Ac, i), min, max);
	}
	return nclip;
}

/**
   raise each cell in the cell array to power of power.
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
   Evaluate 2D cubic spline at location defined 2-d arrays by xnew, ynew
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
