/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#include <fftw3.h>
#include "mathdef.h"
#include "defs.h"
//This file is only used by cmath and zmath, which defines COMP_COMPLEX
#if defined(COMP_COMPLEX) && (!defined(COMP_SINGLE) || HAS_FFTWF==1)

/**
   load FFT wisdom from file.
*/
static char fnwisdom[64];
static void load_wisdom(){
	FILE* fpwisdom;
	if((fpwisdom=fopen(fnwisdom, "r"))){
		FFTW(import_wisdom_from_file)(fpwisdom);
		fclose(fpwisdom);
	}
}
/**
   save FFT wisdom to file.
*/

//Put the following within COMP_COMPLEX to avoid run multiple copies of it.
static void save_wisdom(){
	FILE* fpwisdom;
	if((fpwisdom=fopen(fnwisdom, "w"))){
		FFTW(export_wisdom_to_file)(fpwisdom);
		fclose(fpwisdom);
	}
}
/**
   executed after main() exits.
*/
static __attribute__((destructor))void deinit(){
	save_wisdom();
}
static void (*init_threads())(int){
	void (*p_fftw_plan_with_nthreads)(int)=NULL;
#if HAS_FFTW_THREADS
	const char *suffix="_threads";
#ifdef COMP_SINGLE
	fftwf_init_threads();
	p_fftw_plan_with_nthreads=(void(*)(int))fftwf_plan_with_nthreads;
#else
	fftw_init_threads();
	p_fftw_plan_with_nthreads=(void(*)(int))fftw_plan_with_nthreads;
#endif
#else
	dbg("FFTW thread library is not evailable\n");
	const char *suffix="_serial";
#endif	
#ifdef COMP_SINGLE
	const char *libname="fftwf";
#else
	const char *libname="fftw";
#endif
	sprintf(fnwisdom, "%s/%s_wisdom%s", DIRCACHE, libname, suffix);
	load_wisdom();
	return p_fftw_plan_with_nthreads;
}



#if HAS_FFTW_CALLBACK
//sometimes the callback function is not included in the header.
typedef void callback_fun(void *(*work)(char *), char *jobdata, size_t elsize, int njobs, void *data);
extern void fftw_threads_set_callback(callback_fun callback, void *callback_data);
static void FFTW(parallel_callback)(void *(*work)(char *), char *jobdata, size_t elsize, int njobs, void *data){
	(void)data;
	//dbg("fft_task_callback: %d threads\n", njobs);
#if _OPENMP
	//Adaptively handling in or outside of parallel region.
	if(NTHREAD>1&&!omp_in_parallel()){
		#pragma omp parallel for default(shared) 
		for(int i=0; i<njobs; ++i){
			work(jobdata+elsize*i);
		}
	}else{
		#pragma omp taskloop default(shared) num_tasks(njobs) priority(1)
		for(int i=0; i<njobs; ++i){
			work(jobdata+elsize*i);
		}
	}
#else
	tp_counter_t group={0};
	for(int i=0; i<njobs; ++i){
		QUEUE(&group, work, jobdata+elsize*i, 1, 1);
	}
	WAIT(&group, 1);
#endif
}
#endif
static void FFTW(fft_threads)(long nx, long ny){
	static int fft_has_threads=-1;
	static int last_nthread=-1;
	static void (*p_fftw_plan_with_nthreads)(int)=NULL;
	if(fft_has_threads==-1){//first initialization
		if((p_fftw_plan_with_nthreads=init_threads())){
#if HAS_FFTW_CALLBACK
			fftw_threads_set_callback(FFTW(parallel_callback), NULL); //since version 3.3.9
			dbg("fftw_thread_set_callback is set.\n");
#else			
			dbg("fftw_thread_set_callback is not available.\n");
#endif
			fft_has_threads=1;
		} else{
			fft_has_threads=0;
		}
	}
	if(fft_has_threads==1){
		int nth=(nx*ny>256*256)?NTHREAD:1;
		if(nth!=last_nthread){
			dbg3("FFTW %ldx%ld using %d threads \n", nx, ny, nth);
			p_fftw_plan_with_nthreads(nth);
			last_nthread=nth;
		}
	}
}

#ifdef __cplusplus
#define COMP(A) reinterpret_cast<fftw_complex*>(A)
#else
#define COMP(A) A
#endif
/**
   An arrays of 1-d plans that are used to do 2-d FFTs only over specified region.
*/
typedef struct PLAN1D_T{
	int ncomp;         /**< For a NxN array, only convert center ncomp*ncomp to Fourier space. */
	FFTW(plan) plan[3]; /**< Array of plans for 1-d FFT */
}PLAN1D_T;

struct fft_t{
	FFTW(plan) plan[3];
	PLAN1D_T *plan1d[3];
};
/**
   Free FFTW plans.
*/
void X(fft_free_plan)(fft_t *fft){
	if(!fft) return;
	for(int idir=-1; idir<2; idir+=2){
		LOCK_FFT;
		if(fft->plan1d[idir+1]){
			FFTW(destroy_plan)(fft->plan1d[idir+1]->plan[0]);
			FFTW(destroy_plan)(fft->plan1d[idir+1]->plan[1]);
			FFTW(destroy_plan)(fft->plan1d[idir+1]->plan[2]);
			dbg3("Plan %p destroyed\n", fft->plan1d[idir+1]);
			free(fft->plan1d[idir+1]);
		}
		if(fft->plan[idir+1]){
			FFTW(destroy_plan)(fft->plan[idir+1]);
			dbg3("Plan %p destroyed\n", fft->plan[idir+1]);

		}
		UNLOCK_FFT;
	}
	free(fft);
}
/**
   Create FFTW plans for 2d FFT transforms. This operation destroyes the data in
   the array if FFTW_FLAGS is not FFTW_ESTIMATE. So do it before filling in data.
*/
static void X(fft2plan)(X(mat)* A, int dir){
	assert(abs(dir)==1&&A&&P(A));
	if(!A->fft){
		A->fft=mycalloc(1, fft_t);
	} else if(A->fft->plan[dir+1]) return;
	int FFTW_FLAGS;
	FFTW_FLAGS=FFTW_ESTIMATE;//Always use ESTIMATE To avoid override data.
	LOCK_FFT;
	FFTW(fft_threads)(A->nx, A->ny);
	/*!!fft uses row major mode. so need to reverse order */
	if(A->nx==1||A->ny==1){
		A->fft->plan[dir+1]=FFTW(plan_dft_1d)(A->ny*A->nx, COMP(P(A)), COMP(P(A)), dir, FFTW_FLAGS);
	} else{
		A->fft->plan[dir+1]=FFTW(plan_dft_2d)(A->ny, A->nx, COMP(P(A)), COMP(P(A)), dir, FFTW_FLAGS);
	}
	if(!A->fft->plan[dir+1]){
		error("Plan is empty\n");
	}
	UNLOCK_FFT;
	dbg3("Plan %p created\n", A->fft->plan[dir+1]);
}

/**
   make plans for cfft2partial
*/
static void X(fft2partialplan)(X(mat)* A, int ncomp, int dir){
	assert(abs(dir)==1);
	if(!A->fft){
		A->fft=mycalloc(1, fft_t);
	} else if(A->fft->plan1d[dir+1]) return;
	const int nx=A->nx;
	const int ny=A->ny;
	int FFTW_FLAGS;
	FFTW_FLAGS=FFTW_ESTIMATE;
	PLAN1D_T* plan1d=A->fft->plan1d[dir+1]=mycalloc(1, PLAN1D_T);
	LOCK_FFT;
	FFTW(fft_threads)(A->nx, A->ny);
	/*along columns for all columns. */
	plan1d->plan[0]=FFTW(plan_many_dft)(1, &nx, ny,
		COMP(P(A)), NULL, 1, nx,
		COMP(P(A)), NULL, 1, nx,
		dir, FFTW_FLAGS);
	/*selected along rows, beginning */
	plan1d->plan[1]=FFTW(plan_many_dft)(1, &ny, ncomp/2,
		COMP(P(A)), NULL, nx, 1,
		COMP(P(A)), NULL, nx, 1,
		dir, FFTW_FLAGS);
	/*selected along rows, end */
	plan1d->plan[2]=FFTW(plan_many_dft)(1, &ny, ncomp/2,
		COMP(P(A))+nx-ncomp/2, NULL, nx, 1,
		COMP(P(A))+nx-ncomp/2, NULL, nx, 1,
		dir, FFTW_FLAGS);
	if(!plan1d->plan[0]||!plan1d->plan[1]||!plan1d->plan[2]){
		error("Plan is empty\n");
	}
	UNLOCK_FFT;
	plan1d->ncomp=ncomp;
	dbg3("Plan %p created\n", A->fft->plan1d[dir+1]);
}

/**
   Do 2d FFT transforms.
*/
void X(fft2)(X(mat)* A, int dir){
	assert(abs(dir)==1); assert(A&&P(A));
	/*do 2d FFT on A. */
	if(!A->fft||!A->fft->plan[dir+1]){
		X(fft2plan)(A, dir);//Uses FFTW_ESTIMATE to avoid override data.
	}
	FFTW(execute)(A->fft->plan[dir+1]);
}

/**
   Do 2d inverse FFT (scaling factor of 1/(nx*ny) is applied)
*/
void X(fft2i)(X(mat)* A, int dir){
	/*Apply scaling factor*/
	X(fft2)(A, dir);
	X(scale)(A, 1./(R)(A->nx*A->ny));
}

/**
   Do 2d FFT transforms and scale the output by 1/sqrt(nx*ny)
*/
void X(fft2s)(X(mat)* A, int dir){/*symmetrical cfft2. */
	X(fft2)(A, dir);
	X(scale)(A, 1./sqrt((R)(A->nx*A->ny)));
}

/**
   Apply 2d FFT partially over ncomp \f$\times\f$ncomp region use two 1d plans that
   takes 1d fft through its column and selectly for its rows to produce smaller
   psf.
*/
void X(fft2partial)(X(mat)* A, int ncomp, int dir){
	assert(abs(dir)==1);
	assert(A&&P(A));
	if(!A->fft||!A->fft->plan1d[dir+1]){
		X(fft2partialplan)(A, ncomp, dir);
	}
	PLAN1D_T* plan1d=A->fft->plan1d[dir+1];
	if(ncomp!=plan1d->ncomp) error("Plan and fft mismatch\n");
	for(int i=0; i<3; i++){
		FFTW(execute)(plan1d->plan[i]);
	}
}

/**
   returns IFFT(fftshift(FFT(A)))
*/
X(mat)* X(ffttreat)(X(mat)* A){
	if(!A) return NULL;
	X(mat)* B=X(new)(A->nx, A->ny);
	X(fft2plan)(B, 1);
	X(fft2plan)(B, -1);
	X(cp)(&B, A);
	X(fft2)(B, -1);
	X(fftshift)(B);
	X(fft2)(B, 1);
	X(scale)(B, 1./(R)(A->nx*A->ny));
	return B;
}
void XR(fft_free_plan)(fft_t *fft){
	if(fft) X(fft_free_plan)(fft);
}
/**
   Create a fftw plan based on a 2 element cell array that contains real/imaginary parts respectively.
 */
static void XR(cell_fft2plan)(XR(cell)* dc, int dir){
	assert(abs(dir)==1);
	if(dc->nx*dc->ny!=2){
		error("XR(cell) of two elements is required\n");
	}
	int nx=P(dc,0)->nx;
	int ny=P(dc,0)->ny;
	if(P(dc,1)->nx!=nx||P(dc,1)->ny!=ny){
		error("The two elements in XR(cell) must be of the same size\n");
	}
	fftw_iodim dims[2]={{nx,1,1},{ny,nx,nx}};
	fftw_iodim howmany_dims={1,1,1};
	R* restrict p1=P(P(dc,0));
	R* restrict p2=P(P(dc,1));
	/*Use FFTW_ESTIMATE since the size may be large, and measuring takes too long. */
	fft_t* fft=mycalloc(1, fft_t);
	if(!fft->plan[dir+1]){
		LOCK_FFT;
		FFTW(fft_threads)(nx, ny);
		fft->plan[dir+1]=FFTW(plan_guru_split_dft)
			(2, dims, 1, &howmany_dims, p1, p2, p1, p2, FFTW_ESTIMATE);
		if(!fft->plan[dir+1]){
			error("Plan is empty\n");
		}
		UNLOCK_FFT;
	}
	dc->fft=fft;
}
/**
   Do FFT based on a 2 element cell array that contains real/imaginary parts respectively.
 */
void XR(cell_fft2)(XR(cell)* dc, int dir){
	assert(abs(dir)==1);
	if(!dc->fft||!dc->fft->plan[dir+1]){
		XR(cell_fft2plan)(dc, dir);
	}
	FFTW(execute)(dc->fft->plan[dir+1]);
}
/**
   Create a fftw plan for 1d real to real FFT.
 */
void XR(fft1plan_r2hc)(XR(mat)* A, int dir){
	if(A->nx!=1&&A->ny!=1){
		error("not supported\n");
	}
	assert(abs(dir)==1&&A&&P(A));
	if(!A->fft) A->fft=mycalloc(1, fft_t);
	int FFTW_FLAGS;
	FFTW_FLAGS=FFTW_ESTIMATE;
	LOCK_FFT;
	if(!A->fft->plan[dir+1]){
		if(dir==-1){
			A->fft->plan[dir+1]=FFTW(plan_r2r_1d)(A->nx*A->ny, P(A), P(A), FFTW_R2HC, FFTW_FLAGS);
		} else{
			A->fft->plan[dir+1]=FFTW(plan_r2r_1d)(A->nx*A->ny, P(A), P(A), FFTW_HC2R, FFTW_FLAGS);
		}
		if(!A->fft->plan[dir+1]){
			error("fftw_plan_r2r_1d: Plan is empty. Please check FFT library.\n");
		}
	}
	UNLOCK_FFT;
}
/**
   Do 1d real to real FFT.
 */
void XR(fft1)(XR(mat)* A, int dir){
	assert(A->fft&&abs(dir)==1);
	FFTW(execute)(A->fft->plan[dir+1]);
}

#endif //if defined(COMP_COMPLEX) && (!defined(COMP_SINGLE) || HAS_FFTWF==1) 
