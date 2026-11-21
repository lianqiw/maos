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
#include "kernel.h"

/* A few kernels. */

/**
   Add piston and tip/tilt to opd defined on loc: opd+=pis+x*tx+y*ty. Grid-stride loop.

   @param[in,out] opd OPD array of length n.
   @param[in] loc Coordinates (x, y) of the n points.
   @param[in] n Number of points.
   @param[in] pis Piston.
   @param[in] tx Tip (x tilt) coefficient.
   @param[in] ty Tilt (y tilt) coefficient.
*/
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2],
	int n, Real pis, Real tx, Real ty){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		opd[i]+=pis+loc[i][0]*tx+loc[i][1]*ty;
	}
}
/**
   Same as above, but pis, tx, ty are further scaled by ptt[0], ptt[1], ptt[2],
   which live in device memory.

   @param[in,out] opd OPD array of length n.
   @param[in] loc Coordinates (x, y) of the n points.
   @param[in] n Number of points.
   @param[in] ptt Device array of 3 scaling factors for pis, tx, ty.
   @param[in] pis Piston.
   @param[in] tx Tip (x tilt) coefficient.
   @param[in] ty Tilt (y tilt) coefficient.
*/
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2],
	int n, Real* ptt, Real pis, Real tx, Real ty){
	pis*=ptt[0];
	tx*=ptt[1];
	ty*=ptt[2];
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		opd[i]+=pis+loc[i][0]*tx+loc[i][1]*ty;
	}
}
/**
   Add focus mode to opd: opd+=focus*((x-ox)^2+(y-oy)^2), centered at (ox, oy).

   @param[in,out] opd OPD array of length n.
   @param[in] loc Coordinates (x, y) of the n points.
   @param[in] n Number of points.
   @param[in] focus Focus coefficient.
   @param[in] ox Center x coordinate.
   @param[in] oy Center y coordinate.
*/
__global__ void add_focus_do(Real* restrict opd, Real(*restrict loc)[2],
	int n, Real focus, Real ox, Real oy){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		opd[i]+=((loc[i][0]-ox)*(loc[i][0]-ox)+(loc[i][1]-oy)*(loc[i][1]-oy))*focus;
	}
}
/**
   Add NGS modes (tip/tilt, plate scale ps1-ps3, astigmatism, focus) to opd:
   opd+=alpha*sum(mode*coeff). The plate scale modes depend on the direction
   (thetax, thetay), conjugation height ht and cone scaling factor scale.

   @param[in,out] opd OPD array of length n.
   @param[in] loc Coordinates (x, y) of the n points.
   @param[in] n Number of points.
   @param[in] ttx Tip coefficient (mode x).
   @param[in] tty Tilt coefficient (mode y).
   @param[in] ps1 Plate scale mode 1 coefficient.
   @param[in] ps2 Plate scale mode 2 coefficient.
   @param[in] ps3 Plate scale mode 3 coefficient.
   @param[in] astig1 Astigmatism coefficient (mode x^2-y^2).
   @param[in] astig2 Astigmatism coefficient (mode x*y).
   @param[in] focus Focus coefficient (mode x^2+y^2).
   @param[in] thetax Direction x angle (rad).
   @param[in] thetay Direction y angle (rad).
   @param[in] scale Cone effect scaling factor.
   @param[in] ht Conjugation height (m).
   @param[in] alpha Overall scaling factor.
*/
__global__ void add_ngsmod_do(Real* restrict opd, Real(*restrict loc)[2], int n,
	Real ttx, Real tty,
	Real ps1, Real ps2, Real ps3,
	Real astig1, Real astig2, Real focus,
	Real thetax, Real thetay, Real scale, Real ht, Real alpha
){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		Real x=loc[i][0];
		Real y=loc[i][1];
		Real xy=x*y;
		Real x2=x*x;
		Real y2=y*y;
		opd[i]+=alpha*(+x*ttx
			+y*tty
			+focus*(x2+y2)
			+astig1*(x2-y2)
			+astig2*(xy)
			+ps1*(-2*scale*ht*(thetax*x+thetay*y))
			+ps2*(-2*scale*ht*(thetax*x-thetay*y))
			+ps3*(  -scale*ht*(thetay*x+thetax*y))
		);
	}
}

/**
   a=a*alpha+abs(b)^2*beta, where b is complex.

   @param[in,out] a Real array of length n.
   @param[in] alpha Scaling factor of a.
   @param[in] b Complex array of length n.
   @param[in] beta Scaling factor of abs(b)^2.
   @param[in] n Number of elements.
*/
__global__ void addcabs2_do(Real* restrict a, Real alpha,
	const Comp* restrict b, Real beta, int n){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]=a[i]*alpha+CABS2(b[i])*beta;
	}
}
/**
   a+=abs(b)^2*beta, where b is complex.

   @param[in,out] a Real array of length n.
   @param[in] b Complex array of length n.
   @param[in] beta Scaling factor of abs(b)^2.
   @param[in] n Number of elements.
*/
__global__ void addcabs2_do(Real *restrict a, 
	const Comp *restrict b, Real beta, int n){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]+=+CABS2(b[i])*beta;
	}
}
/*reduction routines*/

/**
   Reduction: *res+=sum(a), accumulated in double precision. *res must be
   initialized by the caller. Needs blockDim.x*sizeof(double) dynamic shared
   memory with blockDim.x a power of 2; see sum_wrap().

   @param[in,out] res Result in device memory, accumulated into.
   @param[in] a Input array of length n.
   @param[in] n Number of elements.
*/
__global__ void sum_do(Real* restrict res, const Real* a, const int n){
	extern __shared__ double sbd[];
	sbd[threadIdx.x]=0;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		sbd[threadIdx.x]+=a[i];
	}
	/*
	  Was using step for above loop. Some threads may proceed faster and
	  replaces step with new value, thus corrfsming the computation.
	 */
	for(int step=(blockDim.x>>1);step>0;step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
			sbd[threadIdx.x]+=sbd[threadIdx.x+step];
		}
	}
	if(threadIdx.x==0){
		atomicAdd(res, sbd[0]);
	}
}
/**
  In each block, we first do the reduction in each warp. This avoid syncthreads
  and if test. Then we copy results from each wrap to the first wrap and do the
  reduction again. *res+=sum(a); *res must be initialized by the caller.

  @param[in,out] res Result in device memory, accumulated into.
  @param[in] a Input array of length n.
  @param[in] n Number of elements.
*/
__global__ void sum2_do(Real* restrict res, const Real* a, const int n){
	typedef double sumtype;//use double to reduce truncation errors
	__shared__ sumtype ssb[REDUCE_WRAP*REDUCE_STRIDE];
	const int idx=threadIdx.x;
	const int wrap=idx/WRAP_SIZE; //which wrap
	const int jdx=(WRAP_SIZE-1)&idx;//index within this wrap
	volatile sumtype* s=ssb+REDUCE_STRIDE*wrap+jdx+WRAP_SIZE/2;
	s[-16]=0;
	//Read in vector from global mem
	sumtype sum=0;
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+idx; i<n; i+=step){
		sum+=a[i];
	}
	s[0]=sum;
	//Handle each wrap without sync
#pragma unroll
	for(int i=0; i<5; i++){
		int offset=1<<i;
		sum+=s[-offset];//every thread retrives current value
		s[0]=sum;//every thread write new value.
	}
	__syncthreads();//synchronize different wraps
	if(idx<REDUCE_WRAP){//use a few threads for reduce
		sumtype sum2=ssb[REDUCE_STRIDE*idx+WRAP_SIZE/2+WRAP_SIZE-1];
		//reuse sb for size of REDUCE_WRAP+REDUCE_WRAP/2;
		ssb[idx]=0;
		volatile sumtype* s2=ssb+REDUCE_WRAP/2+idx;
		s2[0]=sum2;
#pragma unroll	
		for(int i=0; i<REDUCE_WRAP_LOG2; i++){
			int offset=1<<i;
			sum2+=s2[-offset];
			s2[0]=sum2;
		}
		if(idx+1==REDUCE_WRAP){
			atomicAdd(res, sum2);
		}
	}
}

/**
   Reduction: *res=max(*res, max(a)). Values below 0 are ignored since the
   per-thread maximum starts at 0. *res must be initialized by the caller.
   Needs blockDim.x*sizeof(Real) dynamic shared memory with blockDim.x a power
   of 2; see max_wrap().

   @param[in,out] res Result in device memory, updated.
   @param[in] a Input array of length n.
   @param[in] n Number of elements.
*/
__global__ void max_do(Real* restrict res, const Real* a, const int n){
	extern __shared__ Real sb[];
	sb[threadIdx.x]=0;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		if(sb[threadIdx.x]<a[i]) sb[threadIdx.x]=a[i];
	}
	for(int step=(blockDim.x>>1);step>0;step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
			if(sb[threadIdx.x]<sb[threadIdx.x+step]){
				sb[threadIdx.x]=sb[threadIdx.x+step];
			}
		}
	}
	if(threadIdx.x==0){
		atomicMax(res, sb[0]);
	}
}

/**
   Reduction: *res=max(*res, max(abs(a))). *res must be initialized by the
   caller. Needs blockDim.x*sizeof(Real) dynamic shared memory with blockDim.x a
   power of 2; see maxabs_wrap().

   @param[in,out] res Result in device memory, updated.
   @param[in] a Input array of length n.
   @param[in] n Number of elements.
*/
__global__ void maxabs_do(Real* restrict res, const Real* a, const int n){
	extern __shared__ Real sb[];
	sb[threadIdx.x]=0;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		if(sb[threadIdx.x]<Z(fabs)(a[i])) sb[threadIdx.x]=Z(fabs)(a[i]);
	}
	for(int step=(blockDim.x>>1);step>0;step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
			if(sb[threadIdx.x]<sb[threadIdx.x+step]){
				sb[threadIdx.x]=sb[threadIdx.x+step];
			}
		}
	}
	if(threadIdx.x==0){
		atomicMax(res, sb[0]);
	}
}
/**
   Reduction: *res_add+=sum(a.*b), accumulated in double precision. Needs
   blockDim.x*sizeof(double) dynamic shared memory with blockDim.x a power of
   2; see inn_wrap().

   @param[in,out] res_add Result in device memory, accumulated into.
   @param[in] a First array of length n.
   @param[in] b Second array of length n.
   @param[in] n Number of elements.
*/
__global__ void inn_do(Real* res_add, const Real* a, const Real* b, const int n){
	extern __shared__ double sbd[];
	sbd[threadIdx.x]=0;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		sbd[threadIdx.x]+=a[i]*b[i];
	}
	for(int step=(blockDim.x>>1);step>0;step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
			sbd[threadIdx.x]+=sbd[threadIdx.x+step];
		}
	}
	if(threadIdx.x==0){
		atomicAdd(res_add, sbd[0]);
	}
}
/**
   Embed real to complex data with zero imaginary part.

   @param[out] out Complex array of length nx.
   @param[in] in Real array of length nx.
   @param[in] nx Number of elements.
*/
__global__ void embed_do(Comp* out, Real* in, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[ix].x=in[ix];
		out[ix].y=0;
	}
}
/**
   Extract real part from complex data.

   @param[out] out Real array of length nx.
   @param[in] in Complex array of length nx.
   @param[in] nx Number of elements.
*/
__global__ void extract_do(Real* out, Comp* in, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[ix]=Z(cuCreal)(in[ix]);
	}
}
/**
   Forward permutation of complex data: out[i]=in[perm[i]].

   @param[out] out Output array of length nx.
   @param[in] in Input array.
   @param[in] perm Permutation vector of length nx.
   @param[in] nx Number of elements.
*/
__global__ void perm_f_do(Comp* restrict out, const Comp* restrict in, int* restrict perm, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[ix]=in[perm[ix]];
	}
}
/**
   Inverse permutation of complex data: out[perm[i]]=in[i].

   @param[out] out Output array.
   @param[in] in Input array of length nx.
   @param[in] perm Permutation vector of length nx.
   @param[in] nx Number of elements.
*/
__global__ void perm_i_do(Comp* restrict out, const Comp* restrict in, int* restrict perm, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[perm[ix]]=in[ix];
	}
}
/**
   Forward permutation of real data: out[i]=in[perm[i]].

   @param[out] out Output array of length nx.
   @param[in] in Input array.
   @param[in] perm Permutation vector of length nx.
   @param[in] nx Number of elements.
*/
__global__ void perm_f_do(Real* restrict out, const Real* restrict in, int* restrict perm, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[ix]=in[perm[ix]];
	}
}
/**
   Inverse permutation of real data: out[perm[i]]=in[i].

   @param[out] out Output array.
   @param[in] in Input array of length nx.
   @param[in] perm Permutation vector of length nx.
   @param[in] nx Number of elements.
*/
__global__ void perm_i_do(Real* restrict out, const Real* restrict in, int* restrict perm, int nx){
	const int step=blockDim.x*gridDim.x;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=step){
		out[perm[ix]]=in[ix];
	}
}

/**
   Convert opd and amplitude defined on nloc points to complex wavefront
   amp*exp(i*2*pi/wvl*opd) and embed it into wvf at indices embed[]. Only the
   embedded points of wvf are written; the caller zeros it beforehand.

   @param[in,out] wvf Complex wavefront array; points at embed[] are overwritten.
   @param[in] opd OPD of the nloc points (m).
   @param[in] amp Amplitude of the nloc points.
   @param[in] embed Index of each point in wvf.
   @param[in] nloc Number of points.
   @param[in] wvl Wavelength (m).
*/
__global__ void embed_wvf_do(Comp* restrict wvf,
	const Real* restrict opd, const Real* restrict amp,
	const int* embed, const int nloc, const Real wvl){
	const Real pi2l=2.f*M_PI/wvl;
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nloc; ix+=blockDim.x*gridDim.x){
		Real s, c;
		Z(sincos)(pi2l*opd[ix], &s, &c);
		wvf[embed[ix]].x=amp[ix]*c;
		wvf[embed[ix]].y=amp[ix]*s;
	}
}

/**
   Embed or crop an array to another array, moving the corners of in (FFT
   ordering) to the center of out. Only the overlapping region is written.

   @param[out] out Output array of noutx*nouty.
   @param[in] noutx Size along x (contiguous dimension) of out.
   @param[in] nouty Size along y of out.
   @param[in] in Input array of ninx*niny.
   @param[in] ninx Size along x (contiguous dimension) of in.
   @param[in] niny Size along y of in.
*/
__global__ void corner2center_do(Comp* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny){
	int nx, ny;
	ny=MIN(niny, nouty)>>1;
	nx=MIN(ninx, noutx)>>1;
	int noutx2=noutx>>1;
	int nouty2=nouty>>1;
	for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
			out[(iy+nouty2)*noutx+(ix+noutx2)]=in[iy*ninx+ix];
			out[(iy+nouty2)*noutx+(noutx2-1-ix)]=in[iy*ninx+(ninx-1-ix)];
			out[(nouty2-1-iy)*noutx+(noutx2-1-ix)]=in[(niny-1-iy)*ninx+(ninx-1-ix)];
			out[(nouty2-1-iy)*noutx+(ix+noutx2)]=in[(niny-1-iy)*ninx+(ix)];
		}
	}
}

/**
   Same as corner2center_do(), but accumulates abs(in)^2 into out.

   @param[in,out] out Real output array of noutx*nouty, accumulated into.
   @param[in] noutx Size along x (contiguous dimension) of out.
   @param[in] nouty Size along y of out.
   @param[in] in Complex input array of ninx*niny.
   @param[in] ninx Size along x (contiguous dimension) of in.
   @param[in] niny Size along y of in.
*/
__global__ void corner2center_abs2_do(Real* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny){
	int nx, ny;
	ny=MIN(niny, nouty)>>1;
	nx=MIN(ninx, noutx)>>1;
	int noutx2=noutx>>1;
	int nouty2=nouty>>1;
	for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
			out[(iy+nouty2)*noutx+(ix+noutx2)]+=CABS2(in[iy*ninx+ix]);
			out[(iy+nouty2)*noutx+(noutx2-1-ix)]+=CABS2(in[iy*ninx+(ninx-1-ix)]);
			out[(nouty2-1-iy)*noutx+(noutx2-1-ix)]+=CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]);
			out[(nouty2-1-iy)*noutx+(ix+noutx2)]+=CABS2(in[(niny-1-iy)*ninx+(ix)]);
		}
	}
}

/**
   Same as corner2center_abs2_do(), but accumulates with atomicAdd so that
   multiple kernels can add to out concurrently.

   @param[in,out] out Real output array of noutx*nouty, accumulated into.
   @param[in] noutx Size along x (contiguous dimension) of out.
   @param[in] nouty Size along y of out.
   @param[in] in Complex input array of ninx*niny.
   @param[in] ninx Size along x (contiguous dimension) of in.
   @param[in] niny Size along y of in.
*/
__global__ void corner2center_abs2_atomic_do(Real* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny){
	int nx, ny;
	ny=MIN(niny, nouty)>>1;
	nx=MIN(ninx, noutx)>>1;
	int noutx2=noutx>>1;
	int nouty2=nouty>>1;
	for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
			atomicAdd(&out[(iy+nouty2)*noutx+(ix+noutx2)], CABS2(in[iy*ninx+ix]));
			atomicAdd(&out[(iy+nouty2)*noutx+(noutx2-1-ix)], CABS2(in[iy*ninx+(ninx-1-ix)]));
			atomicAdd(&out[(nouty2-1-iy)*noutx+(noutx2-1-ix)], CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]));
			atomicAdd(&out[(nouty2-1-iy)*noutx+(ix+noutx2)], CABS2(in[(niny-1-iy)*ninx+(ix)]));
		}
	}
}
/**
   FFT shift of an nx*ny complex array in place by swapping quadrants. nx and
   ny must be even. Uses a 2-d grid-stride loop over the first quadrant.

   @param[in,out] wvf Complex array of nx*ny.
   @param[in] nx Size along x (contiguous dimension).
   @param[in] ny Size along y.
*/
__global__ void fftshift_do(Comp* wvf, const int nx, const int ny){
	int nx2=nx>>1;
	int ny2=ny>>1;
	for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny2; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx2; ix+=blockDim.x*gridDim.x){
			Comp tmp;
			tmp=wvf[ix+iy*nx];
			wvf[ix+iy*nx]=wvf[(ix+nx2)+(iy+ny2)*nx];
			wvf[(ix+nx2)+(iy+ny2)*nx]=tmp;
			tmp=wvf[ix+(iy+ny2)*nx];
			wvf[ix+(iy+ny2)*nx]=wvf[(ix+nx2)+iy*nx];
			wvf[(ix+nx2)+iy*nx]=tmp;
		}
	}
}
/**
   Add tip/tilt to the array. OPD=OPD+x*ttx+y*tty, where x=ix*dx+ox, y=iy*dx+oy.
   Uses only a single block.

   @param[in,out] opd OPD array of nx*ny.
   @param[in] nx Size along x (contiguous dimension).
   @param[in] ny Size along y.
   @param[in] ox x coordinate of the first point.
   @param[in] oy y coordinate of the first point.
   @param[in] dx Grid spacing in both x and y.
   @param[in] ttx Tip (x tilt) coefficient.
   @param[in] tty Tilt (y tilt) coefficient.
*/
__global__ void add_tilt_do(Real* opd, int nx, int ny, Real ox, Real oy, Real dx, Real ttx, Real tty){
	for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
		Real vty=(oy+iy*dx)*tty;
		for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
			opd[ix+iy*nx]+=vty+(ox+ix*dx)*ttx;
		}
	}
}

/**
   Component wise multiply of complex array by real array: dest=dest.*from.
   Specialization of cwm_do() in kernel.h.

   @param[in,out] dest Complex array of length n.
   @param[in] from Real array of length n.
   @param[in] n Number of elements.
*/
template <>
__global__ void cwm_do(Comp* dest, const Real* from, long n){
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<n; i+=blockDim.x*gridDim.x){
		dest[i].x*=from[i];
		dest[i].y*=from[i];
	}
}

/**
   Component wise multiply of 2-d arrays: dest=dest.*from.

   @param[in,out] dest Complex array with leading dimension lda.
   @param[in] from Complex array with leading dimension ldb.
   @param[in] lda Leading dimension of dest.
   @param[in] ldb Leading dimension of from.
   @param[in] nx Size along x to process.
   @param[in] ny Size along y to process.
*/
__global__ void cwm_do(Comp* dest, Comp* from, int lda, int ldb, int nx, int ny){
	for(int iy=threadIdx.y+blockIdx.y*blockDim.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockIdx.x*blockDim.x; ix<nx; ix+=blockDim.x*gridDim.x){
			dest[ix+iy*lda]*=from[ix+iy*ldb];
		}
	}
}


/**
   Component wise multiply of 2-d arrays: dest=from1.*from2. from1 has same
   shape as dest.

   @param[out] dest Complex array with leading dimension lda.
   @param[in] from1 Complex array with leading dimension lda.
   @param[in] from2 Complex array with leading dimension ldb.
   @param[in] lda Leading dimension of dest and from1.
   @param[in] ldb Leading dimension of from2.
   @param[in] nx Size along x to process.
   @param[in] ny Size along y to process.
*/
__global__ void cwm_do(Comp* dest, Comp* from1, Comp* from2, int lda, int ldb, int nx, int ny){
	for(int iy=threadIdx.y+blockIdx.y*blockDim.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockIdx.x*blockDim.x; ix<nx; ix+=blockDim.x*gridDim.x){
			dest[ix+iy*lda]=from1[ix+iy*lda]*from2[ix+iy*ldb];
		}
	}
}
/**
   Unwrap the phase of wvf back to opd, assuming the result is within wvl/2 of
   the original opd. Each point is handled by exactly one thread, with any 1-d
   launch configuration.

   @param[in] wvf Complex wavefront; point i is at wvf[embed[i]].
   @param[in,out] opd OPD of the n points (m). Input is the reference; output
   is the unwrapped OPD.
   @param[in] embed Index of each point in wvf.
   @param[in] n Number of points.
   @param[in] wvl Wavelength (m).
*/
__global__ void unwrap_phase_do(Comp* wvf, Real* opd, int* embed, int n, Real wvl){
	Real kki=wvl/(2*M_PI);
	Real wvlh=wvl*0.5;
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<n; i+=blockDim.x*gridDim.x){
		Real val=atan2(wvf[embed[i]].y, wvf[embed[i]].x)*kki;
		Real diff=Z(fmod)(val-opd[i]+wvlh, wvl);
		if(diff<0) diff+=wvl;
		opd[i]+=diff-wvlh;
	}
}
/**
   Matrix vector multiply a+=mvm*g, one thread per actuator. The number of
   total threads should not be less than, ideally equal to, nact. The caller
   should setup shared memory to contain blockDim.x*sizeof(Real).

   @param[in] mvm Matrix of size nact*ng.
   @param[in,out] a Vector of size nact*1, accumulated into.
   @param[in] g Vector of size ng*1.
   @param[in] nact Number of rows (actuators).
   @param[in] ng Number of columns (gradients).
*/
__global__ void
mvm_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng){
	extern __shared__ Real acc[];
	int iact=threadIdx.x+blockIdx.x*blockDim.x;
	if(iact<nact){
		acc[threadIdx.x]=0;
		for(int ig=0; ig<ng; ig++){
			Real mvmi=mvm[nact*ig+iact];
			acc[threadIdx.x]+=mvmi*g[ig];
		}
		a[iact]+=acc[threadIdx.x];
	}
}
/**
   Matrix vector multiply a+=mvm*g with the gradients split among multiple
   threads per actuator, combined with atomicAdd. The number of total threads
   should be multiple times of nact, to maximize occupancy. The caller should
   setup shared memory to contain blockDim.x*sizeof(Real).

   @param[in] mvm Matrix of size nact*ng.
   @param[in,out] a Vector of size nact*1, accumulated into.
   @param[in] g Vector of size ng*1.
   @param[in] nact Number of rows (actuators).
   @param[in] ng Number of columns (gradients).
*/
__global__ void
multimv_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng){
	extern __shared__ Real acc[];
	const int ind=threadIdx.x+blockIdx.x*blockDim.x;//max at blockDim.x*gridDim.x
	const int nset=(blockDim.x*gridDim.x)/nact;//total number of over runs
	const int iset=ind/nact;
	if(iset>=nset) return;
	const int iact=ind-nact*iset;//actual actuator
	acc[threadIdx.x]=0;
	/*This block handle grads from igi to ngi.*/
	const int igi=(iset*ng)/nset;
	const int ngi=((iset+1)*ng)/nset;
	for(int ig=igi; ig<ngi; ig++){
		acc[threadIdx.x]+=mvm[nact*ig+iact]*g[ig];
	}
	atomicAdd(&a[iact], acc[threadIdx.x]);
}
