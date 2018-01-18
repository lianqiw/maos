/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include <sys/file.h>
#include "../sys/sys.h"
#include <fftw3.h>
#include "mathdef.h"
#include "defs.h"
#if !defined(USE_SINGLE) || HAS_FFTWF==1
extern int has_threads;
extern void (*p_fftw_plan_with_nthreads)(int n);
static int FFTW_VERBOSE=0;
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
	    /*info("Plan %p destroyed\n", fft->plan1d[idir+1]); */
	    free(fft->plan1d[idir+1]);
	}
	if(fft->plan[idir+1]){
	    FFTW(destroy_plan)(fft->plan[idir+1]);
	    /*info("Plan %p destroyed\n", fft->plan[idir+1]); */

	}
	UNLOCK_FFT;
    }
    free(fft);
}
/**
   load FFT wisdom from file.
*/
#if defined(USE_COMPLEX) && !defined (USE_SINGLE)
int has_threads=-1;
void (*p_fftw_plan_with_nthreads)(int n)=NULL;
//Put the following within USE_COMPLEX to avoid run multiple copies of it.
static char fnwisdom[64];
static void load_wisdom(){
    FILE *fpwisdom;
    if((fpwisdom=fopen(fnwisdom,"r"))){
	int fd=fileno(fpwisdom);
	if(!flock(fd,LOCK_SH)){
	    FFTW(import_wisdom_from_file)(fpwisdom);
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
	    FFTW(export_wisdom_to_file)(fpwisdom);
	    flock(fd,LOCK_UN);
	}else{
	    perror("flock");
	}
	fclose(fpwisdom);
    }
}
/**
   executed after main() exits.
*/
static __attribute__((destructor))void deinit(){
    save_wisdom();
}
/**
   executed before main()
*/
static __attribute__((constructor))void init(){
    READ_ENV_INT(FFTW_VERBOSE, 0, 1);
}
#include <dlfcn.h>
static void init_threads(){
    void *libfftw_threads=NULL;
    const char *fn=0;
    //quitfun is not null when run in matlab.
#if BUILTIN_FFTW_THREADS == 0
#if _OPENMP>200805
#if defined(USE_SINGLE)
    fn="libfftw3f_omp." LDSUFFIX;
#else
    fn="libfftw3_omp." LDSUFFIX;
#endif
#else //else 
#if defined(USE_SINGLE)
    fn="libfftw3f_threads." LDSUFFIX;
#else
    fn="libfftw3_threads." LDSUFFIX;
#endif
#endif //if OPENMP
#endif //if FFTW_THREADS
    
    if(!fn || (libfftw_threads=dlopen(fn, RTLD_LAZY))){
	if(!fn){
	    info2("FFTW thread library is built in\n");
	}else{
	    info2("Open FFTW thread library %s: success\n", fn);
	}
	int (*p_fftw_init_threads)(void)=NULL;
	has_threads=1;
#ifdef USE_SINGLE
	sprintf(fnwisdom, "%s/.aos/fftwf_wisdom_thread",HOME);
#if BUILTIN_FFTW_THREADS
	p_fftw_init_threads=(int(*)())fftwf_init_threads;
	p_fftw_plan_with_nthreads=(void(*)(int))fftwf_plan_with_nthreads;
#else
	p_fftw_init_threads=(int(*)())dlsym(libfftw_threads, "fftwf_init_threads");
	p_fftw_plan_with_nthreads=(void(*)(int))dlsym(libfftw_threads, "fftwf_plan_with_nthreads");
#endif
#else
	sprintf(fnwisdom, "%s/.aos/fftw_wisdom_thread",HOME);
#if BUILTIN_FFTW_THREADS
	p_fftw_init_threads=(int(*)())fftw_init_threads;
	p_fftw_plan_with_nthreads=(void(*)(int))fftw_plan_with_nthreads;
#else
	p_fftw_init_threads=(int(*)())dlsym(libfftw_threads, "fftw_init_threads");
	p_fftw_plan_with_nthreads=(void(*)(int))dlsym(libfftw_threads, "fftw_plan_with_nthreads");
#endif
#endif
	p_fftw_init_threads();
    }else{
	if(!quitfun){
	    info2("Open FFTW thread library %s: failed\n", fn);
	}
	has_threads=0;
#ifdef USE_SINGLE
	sprintf(fnwisdom, "%s/.aos/fftwf_wisdom_serial",HOME);
#else
	sprintf(fnwisdom, "%s/.aos/fftw_wisdom_serial",HOME);
#endif
    }
    load_wisdom();
}
#else //if USE_COMPLEX
static void init_threads(){}//avoid multiple init.
#endif
static void fft_execute(FFTW(plan) plan){
/*
  #if _OPENMP>200805
  if(has_threads && !omp_in_parallel()){//testing purpose.
  OMPTASK_SINGLE{
  FFTW(execute)(plan);
  }
  }else
  #endif
*/
    FFTW(execute)(plan);
}
static void fft_threads(long nx, long ny){
    if(has_threads==-1){
	init_threads();
	if(FFTW_VERBOSE){
	    info2("FFTW: has_threads=%d\n", has_threads);
	}
    }
    if(has_threads==1){
	int nth=(nx*ny>256*256)?NTHREAD:1;
	if(FFTW_VERBOSE){
	    info2("FFTW %ldx%ld using %d threads \n", nx, ny, nth);
	}
	p_fftw_plan_with_nthreads(nth);
    }
}
#ifdef __cplusplus
#define COMP(A) reinterpret_cast<fftw_complex*>(A)
#else
#define COMP(A) A
#endif

#ifdef USE_COMPLEX
/**
   Create FFTW plans for 2d FFT transforms. This operation destroyes the data in
   the array. So do it before filling in data.
*/
static void X(fft2plan)(X(mat) *A, int dir){
    assert(abs(dir)==1 && A && A->p);
    if(!A->fft) {
	A->fft=mycalloc(1,fft_t);
    }else if(A->fft->plan[dir+1]) return;
    int FFTW_FLAGS;
    FFTW_FLAGS=FFTW_ESTIMATE;//Always use ESTIMATE To avoid override data.
    LOCK_FFT;
    fft_threads(A->nx, A->ny);
    /*!!fft uses row major mode. so need to reverse order */
    if(A->nx==1 || A->ny==1){
	A->fft->plan[dir+1]=FFTW(plan_dft_1d)(A->ny*A->nx, COMP(A->p), COMP(A->p), dir, FFTW_FLAGS);
    }else{
	A->fft->plan[dir+1]=FFTW(plan_dft_2d)(A->ny, A->nx,COMP(A->p), COMP(A->p), dir, FFTW_FLAGS);
    }
    if(!A->fft->plan[dir+1]){
	error("Plan is empty\n");
    }
    UNLOCK_FFT;  
    /*info("Plan %p created\n", A->fft->plan[dir+1]); */
}

/**
   make plans for cfft2partial
*/
static void X(fft2partialplan)(X(mat) *A, int ncomp, int dir){
    assert(abs(dir)==1);
    if(!A->fft){
	A->fft=mycalloc(1,fft_t);
    } else if(A->fft->plan1d[dir+1]) return;
    const int nx=A->nx;
    const int ny=A->ny;
    int FFTW_FLAGS;
    FFTW_FLAGS=FFTW_ESTIMATE;
    PLAN1D_T *plan1d=A->fft->plan1d[dir+1]=mycalloc(1,PLAN1D_T);
    LOCK_FFT;
    fft_threads(A->nx, A->ny);
    /*along columns for all columns. */
    plan1d->plan[0]=FFTW(plan_many_dft)(1, &nx, ny,
					COMP(A->p),NULL,1,nx,
					COMP(A->p),NULL,1,nx,
					dir,FFTW_FLAGS);
    /*selected along rows, beginning */
    plan1d->plan[1]=FFTW(plan_many_dft)(1, &ny, ncomp/2,
					COMP(A->p),NULL,nx,1,
					COMP(A->p),NULL,nx,1,
					dir,FFTW_FLAGS);
    /*selected along rows, end */
    plan1d->plan[2]=FFTW(plan_many_dft)(1,&ny,ncomp/2, 
					COMP(A->p)+nx-ncomp/2,NULL,nx,1,
					COMP(A->p)+nx-ncomp/2,NULL,nx,1,
					dir,FFTW_FLAGS);
    UNLOCK_FFT; 
    plan1d->ncomp=ncomp;
    /*info("Plan %p created\n", A->fft->plan1d[dir+1]); */
}

/**
   Do 2d FFT transforms.
*/
void X(fft2)(X(mat) *A, int dir){
    assert(abs(dir)==1); assert(A && A->p);
    /*do 2d FFT on A. */
    if(!A->fft || !A->fft->plan[dir+1]) {
	X(fft2plan)(A, dir);//Uses FFTW_ESTIMATE to avoid override data.
    }
    fft_execute(A->fft->plan[dir+1]);
}

/**
   Do 2d inverse FFT (scaling factor of 1/(nx*ny) is applied)
*/
void X(fft2i)(X(mat) *A, int dir){
    /*Apply scaling factor*/
    X(fft2)(A,dir);
    X(scale)(A,1./(R)(A->nx*A->ny));
}

/**
   Do 2d FFT transforms and scale the output by 1/sqrt(nx*ny)
*/
void X(fft2s)(X(mat) *A, int dir){/*symmetrical cfft2. */
    X(fft2)(A,dir);
    X(scale)(A,1./sqrt((R)(A->nx*A->ny)));
}

/**
   Apply 2d FFT partially over ncomp \f$\times\f$ncomp region use two 1d plans that
   takes 1d fft through its column and selectly for its rows to produce smaller
   psf.  
*/
void X(fft2partial)(X(mat) *A, int ncomp, int dir){
    assert(abs(dir)==1);
    assert(A && A->p);
    if(!A->fft){
	X(fft2partialplan)(A, ncomp, dir);
    }
    PLAN1D_T *plan1d=A->fft->plan1d[dir+1];
    if(!plan1d) error("Please run //cfft2partialplan first\n");
    if(ncomp!=plan1d->ncomp) error("Plan and fft mismatch\n");
    for(int i=0; i<3; i++){
	fft_execute(plan1d->plan[i]);
    }
}

/**
   returns IFFT(fftshift(FFT(A)))
*/
X(mat) *X(ffttreat)(X(mat) *A){
    if (!A) return NULL;
    X(mat) *B=X(new)(A->nx, A->ny);
    X(fft2plan)(B,1);
    X(fft2plan)(B,-1);
    X(cp)(&B, A);
    X(fft2)(B,-1);
    X(fftshift)(B);
    X(fft2)(B,1);
    X(scale)(B,1./(R)(A->nx*A->ny));
    return B;
}
#else
/**
 * Create a fftw plan based on a 2 element X(mat) cell array that contains
 * real/imaginary parts respectively
 */
static void X(cell_fft2plan)(X(cell) *dc, int dir){
    assert(abs(dir)==1);
    if(dc->nx*dc->ny!=2){
	error("X(cell) of two elements is required\n");
    }
    int nx=dc->p[0]->nx;
    int ny=dc->p[0]->ny;
    if(dc->p[1]->nx!=nx || dc->p[1]->ny!=ny){
	error("The two elements in X(cell) must be of the same size\n");
    }
    fftw_iodim dims[2]={{nx,1,1},{ny,nx,nx}};
    fftw_iodim howmany_dims={1,1,1};
    T *restrict p1=dc->p[0]->p;
    T *restrict p2=dc->p[1]->p;
    /*Use FFTW_ESTIMATE since the size may be large, and measuring takes too long. */
    fft_t *fft=mycalloc(1,fft_t);
    if(!fft->plan[dir+1]){
	TIC;tic;
	LOCK_FFT;
	fft_threads(nx, ny);
	fft->plan[dir+1]=FFTW(plan_guru_split_dft)
	    (2, dims, 1, &howmany_dims, p1, p2, p1, p2, FFTW_ESTIMATE);
	if(!fft->plan[dir+1]){
	    error("Plan is empty\n");
	}
	UNLOCK_FFT;
	toc2("done");
    }
    dc->fft=fft;
}
void X(cell_fft2)(X(cell) *dc, int dir){
    assert(abs(dir)==1);
    if(!dc->fft || !dc->fft->plan[dir+1]){
	X(cell_fft2plan)(dc, dir);
    }
    fft_execute(dc->fft->plan[dir+1]);
}

void X(fft1plan_r2hc)(X(mat) *A, int dir){
    if(A->nx!=1 && A->ny!=1){
	error("not supported\n");
    }
    assert(abs(dir)==1 && A && A->p);
    if(!A->fft) A->fft=mycalloc(1,fft_t);
    int FFTW_FLAGS;
    FFTW_FLAGS=FFTW_ESTIMATE;
    LOCK_FFT;
    if(!A->fft->plan[dir+1]){
	if(dir==-1){
	    A->fft->plan[dir+1]=FFTW(plan_r2r_1d)(A->nx*A->ny, A->p, A->p, FFTW_R2HC, FFTW_FLAGS);
	}else{
	    A->fft->plan[dir+1]=FFTW(plan_r2r_1d)(A->nx*A->ny, A->p, A->p, FFTW_HC2R, FFTW_FLAGS);
	}
	if(!A->fft->plan[dir+1]){
	    error("fftw_plan_r2r_1d: Plan is empty. Please check FFT library.\n");
	}
    }
    UNLOCK_FFT;
}

void X(fft1)(X(mat) *A, int dir){
    assert(A->fft && abs(dir)==1);
    fft_execute(A->fft->plan[dir+1]);
}

#endif //#ifdef USE_COMPLEX
#else
void X(fft_free_plan)(fft_t *fft){
    if(fft){
	error("libfftw3f is not available\n");
    }
}
void X(fft2)(X(mat) *A, int dir){
    (void)A; (void)dir;
    error("libfftw3f is not available\n");
}
#endif //!defined(USE_SINGLE) || HAS_FFTWF==1
