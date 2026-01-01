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
#include "mathdef.h"
#include "defs.h"/*Defines T, X, etc */

/*
   The following are functions that are only useful for complex mat only

   notice that fma and modf are slower than native functions.
 */

/**
   component-wise multiply A with conjugate of B: A=A.*conj(B);
 */
void X(cwmc)(X(mat)* restrict A, const X(mat)* restrict B){
	if(!check_mat(A, B) || !B) return;
	const long ntot=A->nx*A->ny;
	for(long i=0; i<ntot; i++){
		P(A,i)*=conj(P(B,i));
	}
}
/**
   component multiply X(mat) A with XR(mat) B. A=A.*B;
 */
void X(cwmd)(X(mat)* restrict A, const XR(mat)* restrict B){
	if(!check_mat(A) || !B) return;
	const long ntot=A->nx*A->ny;
	for(long i=0; i<ntot; i++){
		P(A,i)*=P(B,i);
	}
}

/**
   reshape amp.*exp(-2*pi/lambda*opd) into square and embed into center of A,
   with optional rotation of theta CW.

   benchmark: embed 32x32 into 64x64,
   without rotation takes 0.161 ms.
   with rotation takes 0.330 ms.
*/
void X(embed_wvf)(X(mat)* restrict A, const R* opd, const R* amp,
	const int nopdx, const int nopdy,
	const R wvl, const R theta){
	if(!check_mat(A)) return;
	T* psf=P(A);
	const int npsfx=A->nx;
	const int npsfy=A->ny;

	R wvk=TWOPI/wvl;
	memset(psf, 0, sizeof(T)*npsfx*npsfy);
	if(fabs(theta)<1.e-10){/*no rotation. */
		const int skipx=(npsfx-nopdx)/2;
		const int skipy=(npsfy-nopdy)/2;
		assert(skipx>=0&&skipy>=0);
		T* psf0=psf+skipy*npsfx+skipx;
		for(long iy=0; iy<nopdy; iy++){
			T* psfi=psf0+iy*npsfx;
			const R* opdi=opd+iy*nopdx;
			const R* ampi=amp+iy*nopdx;
			for(int ix=0; ix<nopdx; ix++){
			/*sqrt(1-cos^2)!=sin due to sign problem.  */
				psfi[ix]=ampi[ix]*EXPI(wvk*opdi[ix]);
				/*most time is spend in EXP. */
			}
		}
	} else{
	/*rotated for LGS. */
	/*rotate image CCW to -theta. coordinate rotate in reverse way. */
	/*xnew=x*cos(theta)-y*sin(theta) */
	/*ynew=x*sin(theta)+y*cos(theta); */
	/*rotation is around the fft center: (nx/2, ny/2); */
	/*
	  The original method of adding in to out is not right.
	*/
		const R ctheta=cos(theta);
		const R stheta=sin(theta);
		const int nopdx2=nopdx/2;
		const int npsfx2=npsfx/2;
		const int nopdy2=nopdy/2;
		const int npsfy2=npsfy/2;
		const int maxr=iceil(sqrt(nopdx*nopdx+nopdy*nopdy));
		const int xskip=npsfx>maxr?(npsfx-maxr)/2:0;
		const int yskip=npsfy>maxr?(npsfy-maxr)/2:0;
		for(long iy=yskip; iy<npsfy-yskip; iy++){
			R y=iy-npsfy2;
			for(int ix=xskip; ix<npsfx-xskip; ix++){
				R x=ix-npsfx2;
				R x2=x*ctheta+y*stheta+nopdx2;
				R y2=-x*stheta+y*ctheta+nopdy2;
				if(x2>0&&x2<nopdx-1&&y2>0&&y2<nopdy-1){
					int ix2=ifloor(x2);
					int iy2=ifloor(y2);
					x2=x2-ix2;
					y2=y2-iy2;
					R iopd=opd[iy2*nopdx+ix2]*(1.-x2)*(1.-y2)
						+opd[iy2*nopdx+ix2+1]*(x2*(1.-y2))
						+opd[(iy2+1)*nopdx+ix2]*((1-x2)*y2)
						+opd[(iy2+1)*nopdx+ix2+1]*(x2*y2);
					R iamp=amp[iy2*nopdx+ix2]*(1.-x2)*(1.-y2)
						+amp[iy2*nopdx+ix2+1]*(x2*(1.-y2))
						+amp[(iy2+1)*nopdx+ix2]*((1-x2)*y2)
						+amp[(iy2+1)*nopdx+ix2+1]*(x2*y2);
					psf[iy*npsfx+ix]=iamp*EXPI(wvk*iopd);
				}
			}
		}
	}
}
#define cmpcpy(A,B,S) memcpy(A,B,S*sizeof(T)); /**<macro to do copying*/
/**
   do a modules square and copy to real part of output
*/
static inline void sq2cpy(T* out, const T* in, const long length){
	for(long i=0; i<length; i++){
		out[i]=ABS2(in[i]);
	}
}
/**
   Copy the real part.
*/
static inline void realcpy(T* out, const T* in, const long length){
	for(long i=0; i<length; i++){
		out[i]=REAL(in[i]);
	}
}
/**
   Copy the absolute value to output.
 */
static inline void abscpy(T* out, const T* in, const long length){
	for(long i=0; i<length; i++){
		out[i]=ABS(in[i]);
	}
}
static inline T RA2XY(T A){
	//convert r/a to x/y
	return COMPLEX(REAL(A)*cos(IMAG(A)), REAL(A)*sin(IMAG(A)));
}
static inline T XY2RA(T A){
	//convert x/y to r/a*/
	return COMPLEX(ABS(A), atan2(IMAG(A), REAL(A)));
}

/**
   copy and embed/crop psfin into psfout
   \verbatim
   4 * * 3
   * * * *
   * * * *
   2 * * 1
   \endverbatim
   to
   \verbatim
   4 * * * 3
   * * * * *
   * * * * *
   * * * * *
   2 * * * 1
   \endverbatim
   }
*/
void X(cpcorner)(X(mat)* A, const X(mat)* restrict B, CEMBED flag){
	if(!check_mat(A, B)) return;
	const long nx=A->nx;
	const long ny=A->ny;
	const long ninx=B->nx;
	const long niny=B->ny;
	T* psfout=P(A);
	const T* restrict psfin=P(B);
	assert((nx&1)==0&&(ny&1)==0&&(ninx&1)==0&&(niny&1)==0);
	memset(psfout, 0, sizeof(T)*nx*ny);
	const int ny2=(ny<niny)?ny/2:niny/2;
	const int nx2=(nx<ninx)?nx/2:ninx/2;

#define DO_COPY2(CCMD)			\
	CCMD(psfout+(ny-ny2+i)*nx+(nx-nx2),		\
	     psfin+(niny-ny2+i)*ninx+(ninx-nx2),nx2);	\
	CCMD(psfout+(ny-ny2+i)*nx,			\
	     psfin+(niny-ny2+i)*ninx, nx2);		\
	CCMD(psfout+i*nx+(nx-nx2),			\
	     psfin+i*ninx+(ninx-nx2), nx2);		\
	CCMD(psfout+i*nx,				\
	     psfin+i*ninx, nx2); 
	for(int i=0; i<ny2; i++){
		switch(flag){
		case C_FULL:
			DO_COPY2(cmpcpy);
			break;
		case C_ABS2:
			DO_COPY2(sq2cpy);
			break;
		case C_REAL:
			DO_COPY2(realcpy);
			break;
		case C_ABS:
			DO_COPY2(abscpy);
			break;
		default:
			error("Invalid flag\n");
		}
	}
}
/**
   Put abs square of each element into its realpart.
   A=abs(A).^2;
*/
void X(abs2toreal)(X(mat)* A, R scale){
	if(!check_mat(A)) return;
	for(int i=0; i<A->nx*A->ny; i++){
		P(A,i)=ABS2(P(A,i))*scale;
	}
}
/**
   Put abs of each elemnt into its realpart. A=abs(A);
*/
void X(abstoreal)(X(mat)* A){
	if(!check_mat(A)) return;
	/*put abs to real */
	for(int i=0; i<A->nx*A->ny; i++){
		P(A,i)=ABS(P(A,i));
	}
}
/**
   Copy a XR(mat) into real part of X(mat).
 */
void X(cpd)(X(mat)** restrict A0, const XR(mat)* restrict B){
	if(!B) return;
	X(mat)* restrict A=*A0;

	if(!A){
		*A0=A=X(new)(B->nx, B->ny);
	} else{
		assert((A->nx==B->nx&&A->ny==B->ny));
	}
	for(long i=0; i<B->nx*B->ny; i++){
		P(A,i)=P(B,i);
	}
}
/**
   Copy real part of a X(mat) to XR(mat) with optional scaling:
   A0=A0.*alpha+real(B)*beta
*/

#define CMAT_TO_DMAT(name, opt)\
void X(name)(XR(mat)** restrict A0, R alpha, const X(mat)* restrict B, R beta){\
	if(!check_mat(B)) return;\
	XR(init)(A0, B->nx,B->ny);\
	XR(mat)* restrict A=*A0;\
	if(fabs(alpha)<EPS){\
		for(long i=0; i<B->nx*B->ny; i++){\
			P(A,i)=opt(P(B,i))*beta;\
		}\
	} else{\
		for(long i=0; i<B->nx*B->ny; i++){\
			P(A,i)=P(A,i)*alpha+opt(P(B,i))*beta;\
		}\
	}\
}
CMAT_TO_DMAT(real2d, REAL)
CMAT_TO_DMAT(imag2d, IMAG)
CMAT_TO_DMAT(abs22d, ABS2)
CMAT_TO_DMAT(arg2d, CARG)

/**
   Tilt the otf to make the image shift.

   apply exp(-2*pi*sx/nx)*exp(-2*pi*sy/ny) to otf to move the image
   sx, sy are shifts in terms of pixel.
   sx/nx is equivalent to sx*dtheta*du.
   peak_corner=1: peak is in center
   peak_corner=0: peak is in corner
*/
void X(tilt2)(X(mat)* otf, X(mat)* otfin, R sx, R sy, int peak_corner){
	if(!check_mat(otf, otfin)) return;
	int nx=otf->nx;
	int ny=otf->ny;
	R dux=1./(R)nx;
	R duy=1./(R)ny;
	T ux[nx];
	T uy[ny];
	T cx=EXPI(-TWOPI*dux*sx);
	T cy=EXPI(-TWOPI*duy*sy);
	//warning_once("Consider caching ux, uy\n");
	if(peak_corner==1){/*peak in center */
		ux[0]=EXPI(-TWOPI*dux*sx*(-nx/2));
		for(int i=1; i<nx; i++){
			ux[i]=ux[i-1]*cx;
		}
		uy[0]=EXPI(-TWOPI*duy*sy*(-ny/2));
		for(int i=1; i<ny; i++){
			uy[i]=uy[i-1]*cy;
		}
	} else{
		ux[0]=1;
		for(int i=1; i<nx/2; i++){
			ux[i]=ux[i-1]*cx;
		}
		ux[nx/2]=EXPI(-TWOPI*dux*sx*(-nx/2));
		for(int i=nx/2+1; i<nx; i++){
			ux[i]=ux[i-1]*cx;
		}
		uy[0]=1;
		for(int i=1; i<ny/2; i++){
			uy[i]=uy[i-1]*cy;
		}
		uy[ny/2]=EXPI(-TWOPI*duy*sy*(-ny/2));
		for(int i=ny/2+1; i<ny; i++){
			uy[i]=uy[i-1]*cy;
		}
	}
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(otf, ix, iy)=P(otfin, ix, iy)*ux[ix]*uy[iy];
		}
	}
}
/**
   Inplace tilt the otf to make the image shift.
*/
void X(tilt)(X(mat)* otf, R sx, R sy, int peak_corner){
	X(tilt2)(otf, otf, sx, sy, peak_corner);
}
