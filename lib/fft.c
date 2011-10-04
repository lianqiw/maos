/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#define USE_COMPLEX /**<for complex data*/
#include <sys/file.h>
#include "common.h"
#include <fftw3.h>
#include "thread.h"
#include "process.h"
#include "misc.h"
#include "fft.h"
#include "cmat.h"
PNEW(mutex_fftw);
#define LOCK_FFT LOCK(mutex_fftw)
#define UNLOCK_FFT UNLOCK(mutex_fftw)
/**
   \file fft.c
   Routines to do FFT on cmat.
*/
/**
   An arrays of 1-d plans that are used to do 2-d FFTs only over specified region.
 */
typedef struct PLAN1D_T{
    int ncomp;         /**< For a NxN array, only convert center ncomp*ncomp to Fourier space. */
    fftw_plan plan[3]; /**< Array of plans for 1-d FFT */
}PLAN1D_T;

struct fft_t{
    fftw_plan plan[3];
    PLAN1D_T *plan1d[3];
};
static char fnwisdom[64];
/**
   load FFT wisdom from file.
 */
static void load_wisdom(){
    FILE *fpwisdom;
    if((fpwisdom=fopen(fnwisdom,"r"))){
	int fd=fileno(fpwisdom);
	if(!flock(fd,LOCK_SH)){
	    fftw_import_wisdom_from_file(fpwisdom);
	    flock(fd,LOCK_UN);
	}else{
	    perror("flock");
	}
	fclose(fpwisdom);
    }
}
/**
   save FFT wisdom to file.
 */
static void save_wisdom(){
    FILE *fpwisdom;
    if((fpwisdom=fopen(fnwisdom,"w"))){
	int fd=fileno(fpwisdom);
	if(!flock(fd,LOCK_EX)){
	    fftw_export_wisdom_to_file(fpwisdom);
	    flock(fd,LOCK_UN);
	}else{
	    perror("flock");
	}
	fclose(fpwisdom);
    }
}
#if USE_FFTW_THREADS == 1
#define FFTW_THREADS(n) fftw_plan_with_nthreads(n)
#else
#define FFTW_THREADS(n)
#endif

/**
   executed before main(). FFTW website says "We should also mention one other
   restriction: if you save wisdom from a program using the multi-threaded FFTW,
   that wisdom cannot be used by a program using only the single-threaded FFTW
   (i.e. not calling fftw_init_threads)." So here we use different names.
 */
static __attribute__((constructor))void init(){
#if USE_FFTW_THREADS == 1
    fftw_init_threads();
    sprintf(fnwisdom, "%s/.aos/fftw_wisdom_thread",HOME);
#else
    sprintf(fnwisdom, "%s/.aos/fftw_wisdom_serial",HOME);
#endif
    load_wisdom();
}
/**
   executed after main() exits.
 */
static __attribute__((destructor))void deinit(){
    save_wisdom();
#if USE_FFTW_THREADS == 1
    fftw_cleanup_threads();
#endif
}
/**
   Create FFTW plans for 2d FFT transforms. This operation destroyes the data in
   the array. So do it before filling in data.
 */
void cfft2plan(cmat *A, int dir){
    assert(abs(dir)==1 && A && A->p);
    if(!A->fft) A->fft=calloc(1, sizeof(fft_t));
    int FFTW_FLAGS;
    if(A->nx<1024){
	FFTW_FLAGS=FFTW_MEASURE;
    }else{
	FFTW_FLAGS=FFTW_ESTIMATE;
    }
    LOCK_FFT;
    /*!!fft uses row major mode. so need to reverse order */
    if(A->nx==1 || A->ny==1){
	A->fft->plan[dir+1]=fftw_plan_dft_1d(A->ny*A->nx, A->p, A->p, dir, FFTW_FLAGS);
    }else{
	A->fft->plan[dir+1]=fftw_plan_dft_2d(A->ny, A->nx, A->p, A->p, dir, FFTW_FLAGS);
    }
    UNLOCK_FFT;  
    /*info("Plan %p created\n", A->fft->plan[dir+1]); */
}

/**
   make plans for cfft2partial
*/
void cfft2partialplan(cmat *A, int ncomp, int dir){
    assert(abs(dir)==1);
    if(!A->fft)  A->fft=calloc(1, sizeof(fft_t));
    const int nx=A->nx;
    const int ny=A->ny;
    int FFTW_FLAGS;
    if(A->nx<1024){
	FFTW_FLAGS=FFTW_MEASURE;
    }else{
	FFTW_FLAGS=FFTW_ESTIMATE;
    }
    PLAN1D_T *plan1d=A->fft->plan1d[dir+1]=calloc(1, sizeof(PLAN1D_T));
    LOCK_FFT;
    /*along columns for all columns. */
    plan1d->plan[0]=fftw_plan_many_dft(1, &nx, ny,
				       A->p,NULL,1,nx,
				       A->p,NULL,1,nx,
				       dir,FFTW_FLAGS);
    /*selected along rows, beginning */
    plan1d->plan[1]=fftw_plan_many_dft(1, &ny, ncomp/2,
				       A->p,NULL,nx,1,
				       A->p,NULL,nx,1,
				       dir,FFTW_FLAGS);
    /*selected along rows, end */
    plan1d->plan[2]=fftw_plan_many_dft(1,&ny,ncomp/2, 
				 A->p+nx-ncomp/2,NULL,nx,1,
				 A->p+nx-ncomp/2,NULL,nx,1,
				 dir,FFTW_FLAGS);
    UNLOCK_FFT; 
    plan1d->ncomp=ncomp;
    /*info("Plan %p created\n", A->fft->plan1d[dir+1]); */
}
/**
   Free FFTW plans.
*/
void fft_free_plan(fft_t *fft){
    if(!fft) return;
    for(int idir=-1; idir<2; idir+=2){
	LOCK_FFT;
	if(fft->plan1d[idir+1]){
	    fftw_destroy_plan(fft->plan1d[idir+1]->plan[0]);
	    fftw_destroy_plan(fft->plan1d[idir+1]->plan[1]);
	    fftw_destroy_plan(fft->plan1d[idir+1]->plan[2]);
	    /*info("Plan %p destroyed\n", fft->plan1d[idir+1]); */
	    free(fft->plan1d[idir+1]);
	}
	if(fft->plan[idir+1]){
	    fftw_destroy_plan(fft->plan[idir+1]);
	    /*info("Plan %p destroyed\n", fft->plan[idir+1]); */

	}
	UNLOCK_FFT;
    }
    free(fft);
}
/**
   Do 2d FFT transforms.
 */
void cfft2(cmat *A, int dir){
    assert(abs(dir)==1); assert(A && A->p);
    /*do 2d FFT on A. */
    dir++;
    /*can not do planning here because planning will override the data. */
    if(!A->fft || !A->fft->plan[dir]) error("Please run cfft2plan first\n");
    fftw_execute(A->fft->plan[dir]);
}

/**
   Do 2d inverse FFT (scaling factor of 1/(nx*ny) is applied)
*/
void cifft2(cmat *A, int dir){
    /*Apply scaling factor*/
    cfft2(A,dir);
    cscale(A,1./(double)(A->nx*A->ny));
}

/**
   Do 2d FFT transforms and scale the output by 1/sqrt(nx*ny)
 */
void cfft2s(cmat *A, int dir){/*symmetrical cfft2. */
    cfft2(A,dir);
    cscale(A,1./sqrt((double)(A->nx*A->ny)));
}

/**
   Apply 2d FFT partially over ncomp \f$\times\f$ncomp region use two 1d plans that
   takes 1d fft through its column and selectly for its rows to produce smaller
   psf.  
*/
void cfft2partial(cmat *A, int ncomp, int dir){
    assert(abs(dir)==1);
    assert(A && A->p);
    if(!A->fft) error("Please run cfft2partialplan first\n");
    PLAN1D_T *plan1d=A->fft->plan1d[dir+1];
    if(!plan1d) error("Please run cfft2partialplan first\n");
    if(ncomp!=plan1d->ncomp) error("Plan and fft mismatch\n");
    for(int i=0; i<3; i++){
	fftw_execute(plan1d->plan[i]);
    }
}

/**
   returns IFFT(fftshift(FFT(A)))
 */
cmat *cffttreat(cmat *A){
    if (!A) return NULL;
    cmat *B=cnew(A->nx, A->ny);
    cfft2plan(B,1);
    cfft2plan(B,-1);
    ccp(&B, A);
    cfft2(B,-1);
    cfftshift(B);
    cfft2(B,1);
    cscale(B,1./(double)(A->nx*A->ny));
    return B;
}

/**
 * Create a fftw plan based on a 2 element dmat cell array that contains
 * real/imaginary parts respectively
 */
fft_t* dcell_fft2plan(dcell *dc, int dir, int nthreads){
    assert(abs(dir)==1);
    if(dc->nx*dc->ny!=2){
	error("dcell of two elements is required\n");
    }
    long nx=dc->p[0]->nx;
    long ny=dc->p[0]->ny;
    if(dc->p[1]->nx!=nx || dc->p[1]->ny!=ny){
	error("The two elements in dcell must be of the same size\n");
    }
    fftw_iodim dims[2]={{nx,1,1},{ny,nx,nx}};
    fftw_iodim howmany_dims={1,1,1};
    double *restrict p1=dc->p[0]->p;
    double *restrict p2=dc->p[1]->p;
    /*Use FFTW_ESTIMATE since the size may be large, and measuring takes too long. */
    fft_t *fft=calloc(1, sizeof(fft_t));
    if(!fft->plan[dir+1]){
	TIC;tic;
	LOCK_FFT;
	if(nthreads<1) nthreads=1;
	info2("Creating fft plan with %d threads ...", nthreads);
	FFTW_THREADS(nthreads);
	fft->plan[dir+1]=fftw_plan_guru_split_dft
	    (2, dims, 1, &howmany_dims, p1, p2, p1, p2, FFTW_ESTIMATE);
	FFTW_THREADS(1);
	UNLOCK_FFT;
	toc2("done");
    }
    return fft;
}

/**
 * Free the plan in dcell
 */
void fft2(fft_t *fft, int dir){
    assert(fft && abs(dir)==1);
    fftw_execute(fft->plan[dir+1]);
}
