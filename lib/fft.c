/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <complex.h>
//need to include complex.h before fftw3.h
#include <fftw3.h>
#include "common.h"
#include "misc.h"
#include "fft.h"
#include "cmat.h"
PNEW(mutex_fftw);
#define LOCK_FFT LOCK(mutex_fftw)
#define UNLOCK_FFT UNLOCK(mutex_fftw)
#define FFTW_FLAGS FFTW_MEASURE
/**
   \file fft.c
   Routines to do FFT on cmat.
*/

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
/**
   executed before main().
 */
static __attribute__((constructor))void init(){
    sprintf(fnwisdom, "%s/.aos/fftw_wisdom",HOME);
    load_wisdom();
}
/**
   executed after main() exits.
 */
static __attribute__((destructor))void deinit(){
    save_wisdom();
}
/**
   Create FFTW plans for 2d FFT transforms. This operation destroyes the data in
   the array. So do it before filling in data.
 */
void cfft2plan(cmat *A, int dir){
    assert(abs(dir)==1);
    assert(A && A->p && !A->plan[dir+1]);
    LOCK_FFT;
    //!!fft uses row major mode. so need to reverse order
    if(A->nx==1 || A->ny==1){
	A->plan[dir+1]=fftw_plan_dft_1d(A->ny*A->nx, A->p, A->p, dir, FFTW_FLAGS);
    }else{
	A->plan[dir+1]=fftw_plan_dft_2d(A->ny, A->nx, A->p, A->p, dir, FFTW_FLAGS);
    }
    UNLOCK_FFT;  
}

/**
   make plans for cfft2partial
*/
void cfft2partialplan(cmat *A, int ncomp, int dir){
    assert(abs(dir)==1);

    const int nx=A->nx;
    const int ny=A->ny;
    PLAN1D_T *plan1d=A->plan1d[dir+1]=calloc(1, sizeof(PLAN1D_T));
    LOCK_FFT;
    //along columns for all columns.
    plan1d->plan[0]=fftw_plan_many_dft(1, &nx, ny,
				       A->p,NULL,1,nx,
				       A->p,NULL,1,nx,
				       dir,FFTW_FLAGS);
    //selected along rows, beginning
    plan1d->plan[1]=fftw_plan_many_dft(1, &ny, ncomp/2,
				       A->p,NULL,nx,1,
				       A->p,NULL,nx,1,
				       dir,FFTW_FLAGS);
    //selected along rows, end
    plan1d->plan[2]=fftw_plan_many_dft(1,&ny,ncomp/2, 
				 A->p+nx-ncomp/2,NULL,nx,1,
				 A->p+nx-ncomp/2,NULL,nx,1,
				 dir,FFTW_FLAGS);
    UNLOCK_FFT; 
    plan1d->ncomp=ncomp;
}
/**
   Free FFTW plans.
*/
void cfree_plan(cmat *A){
    for(int idir=-1; idir<2; idir+=2){
	if(A->plan1d[idir+1]){
	    LOCK_FFT;
	    fftw_destroy_plan(A->plan1d[idir+1]->plan[0]);
	    fftw_destroy_plan(A->plan1d[idir+1]->plan[1]);
	    fftw_destroy_plan(A->plan1d[idir+1]->plan[2]);
	    free(A->plan1d[idir+1]);
	    UNLOCK_FFT;
	}
	if(A->plan[idir+1]){
	    LOCK_FFT;
	    fftw_destroy_plan(A->plan[idir+1]);
	    UNLOCK_FFT;
	}
    }
}
/**
   Do 2d FFT transforms.
 */
void cfft2(cmat *A, int dir){
    assert(abs(dir)==1); assert(A && A->p);
    //do 2d FFT on A.

    dir++;
    //can not do planning here because planning will override the data.
    if(!A->plan[dir]) error("Please run cfft2plan first\n");
    fftw_execute(A->plan[dir]);
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
void cfft2s(cmat *A, int dir){//symmetrical cfft2.
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

    PLAN1D_T *plan1d=A->plan1d[dir+1];
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

