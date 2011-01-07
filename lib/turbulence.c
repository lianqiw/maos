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

#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "common.h"
#include "thread.h"
#include "shm.h"
#include "process.h"
#include "turbulence.h"
#include "path.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
#include "fft.h"
#include "hashlittle.h"
#include "fractal.h"
#include "mathmisc.h"
#include "nr/nr.h"

int disable_atm_shm=0;
/**
 *  \file turbulence.c
 *  Contains routines to generate atmospheric turbulence screens
 */


/**
 *  map the shm to map_t array.
 */
void map_shm(map_t **screen, long totmem, int nlayer, int fd, int rw){
    
    /*
      set access, modification time to current. This information is used to
      detect unused screens
    */

    futimes(fd, NULL);
    int prot = rw ? PROT_READ|PROT_WRITE : PROT_READ;
    int op = rw ? LOCK_EX : LOCK_SH;
    
    /*
      In rw mode, apply an exclusive lock so no other process can use it.
      In ro mode, apply an shared lock.
    */
    if(flock(fd, op)){
	error("Failed to lock file\n");
    }
    /*
      Allocate memory by calling ftruncate in rw mode.
    */
    if(rw && ftruncate(fd, totmem)){
	error("Failed to allocate memory\n");
    }
    if((screen[0]->p=mmap(NULL, totmem, prot, MAP_SHARED ,fd, 0))<0){
	error("Unable to mmap\n");
    }
    /*
      since fd might be 0. we add fd by 1 and assign to shm. This will later be
      used to close the fd and release the lock.
    */
    screen[0]->shm=fd +1;
    int m=screen[0]->nx;
    int n=screen[0]->ny;
    for(int ilayer=1; ilayer<nlayer; ilayer++){
	screen[ilayer]->shm = -1;
	screen[ilayer]->p = screen[0]->p+m*n*ilayer;
    }
}

/**
 *  Allocate for memory for atmosphere. If shm is enabled and has enough shared
 *  memory, will allocate memory in shm, otherwise, allocate memory in
 *  heap. inshm is returned value, and could take the following values:
 *
 * 0: not in shm. 
 * 1: existing screen is resued.
 * 2: in shm and need to generate screen
 *
 * dx, r0, L0, wt, are used to make the key.
 */
map_t **atmnew_shm(int *fd0, int *inshm, rand_t *rstat, long nx, long ny, double dx, 
		   double r0, double L0, double *wt, int nlayer, int method){
    map_t**screen=calloc(nlayer,sizeof(map_t*));
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	screen[ilayer]=mapnew(nx, ny, dx, (void*)1);//give it 1 temporarily.
    }
#if USE_POSIX_SHM == 1
    int use_shm=0;
    char fnshm[NAME_MAX];
    uint32_t key;
    long totmem;
    if(disable_atm_shm){
	use_shm=0;
	shm_free_unused("",0);
    }else{
	key=hashlittle(rstat, sizeof(rand_t), 0);//contains seed
	key=hashlittle(&nx, sizeof(long), key);
	key=hashlittle(&ny, sizeof(long), key);
	key=hashlittle(&dx, sizeof(double), key);
	key=hashlittle(&r0, sizeof(double), key);
	key=hashlittle(&L0, sizeof(double), key);
	key=hashlittle(wt, sizeof(double)*nlayer, key);
	key=hashlittle(&method, sizeof(int), key);
	key=hashlittle(&nlayer, sizeof(int), key);
	snprintf(fnshm,NAME_MAX,"/maos_atm_%ud_%d_%ldx%ld_%g",key,nlayer,nx,ny,dx);
	long already_exist = shm_free_unused(fnshm, 0);
	totmem=(nx*ny*nlayer+1)*sizeof(double);//+1 for sanity check.
	long shm_avail=shm_get_avail();
	if(!already_exist && shm_avail<totmem*1.1){
	    warning2("Need %ld MiB, only %ld MiB available", totmem/1024/1024, shm_avail/1024/1024);
	    warning2("Fall back to non-shared atm");
	    use_shm=0;
	}else{
	    use_shm=1;
	}
    }
    if(use_shm){
	int fd;
	*inshm=1;
    retry:
	//sleep so that the creation daemonize.has enough time to make an exclusive lock.
	fd=shm_open(fnshm, O_RDONLY, 00777); usleep(100);//umask'ed

	if(fd<0){//unable to open. We then try to create
	    fd=shm_open(fnshm, O_RDWR|O_CREAT|O_EXCL, 00777);
	    if(fd<0){//some other process created it in the mean time
		fd=shm_open(fnshm, O_RDONLY, 00777); usleep(100);
		if(fd<0){
		    error("Unable to open shared segment for read.\n");
		}else{
		    if(flock(fd, LOCK_SH)){//wait to make a shared lock
			error("Failed to apply a lock on shared segment.\n");
		    }	
		}
	    }else{
		*inshm=2;
		info2("\nCreating %s ...",fnshm);
		map_shm(screen, totmem, nlayer, fd, 1);
		screen[0]->p[nx*ny*nlayer]=0;//say the data is not valid.
		*fd0=fd;
		return screen;
	    }
	}else{
	    info2("\nReusing %s ...",fnshm);
	}
	*fd0=fd;
	map_shm(screen, totmem, nlayer, fd, 0);
	if(fabs(screen[0]->p[nx*ny*nlayer])<1.e-15){//test validity of data.
	    warning2("Data in shm %s is invalid, free and redo it.\n", fnshm);
	    flock(fd, LOCK_UN);//release lock.
	    close(fd);//close descriptor.
	    if(munmap(screen[0]->p, totmem)){//unmap.
		error("munmap failed\n");
	    }
	    shm_unlink(fnshm);//destroy the map.
	    goto retry;
	}
    }else{
#endif
	for(int ilayer=0; ilayer<nlayer; ilayer++){
	    screen[ilayer]->p=malloc(nx*ny*sizeof(double));
	    screen[ilayer]->shm=0;
	}
	*inshm=0;
#if USE_POSIX_SHM == 1
    }
#endif
    return screen;
}

typedef struct GENSCREEN_T{
    rand_t *rstat;
    dmat *spect;
    map_t **screen;
    double *wt;
    double r0;
    double L0;
    long nlayer;
    long ilayer;
    long ninit;
    pthread_mutex_t mutex_ilayer;
}GENSCREEN_T;
/**
 *   Generate turbulence screens.
 */
static void *genscreen_do(GENSCREEN_T *data){
    const dmat *spect=data->spect;
    rand_t *rstat=data->rstat;
    
    const long m=spect->nx;
    const long n=spect->ny;
    const long nlayer=data->nlayer;
    map_t** screen=data->screen;
    const double *wt=data->wt;
    cmat* cspect=cnew(m,n);
    cfft2plan(cspect,-1);
 repeat:
    LOCK(data->mutex_ilayer);
    int ilayer=data->ilayer;
    data->ilayer+=2;//generate 2 layers at a time.
    if(ilayer>=nlayer){
	UNLOCK(data->mutex_ilayer);
	cfree(cspect);
	return NULL;
    }
    /*We generate random numbers inside mutex lock to make
      sure the random sequence is repeatable.*/
    
    for(long i=0; i<m*n; i++){
	cspect->p[i]=(randn(rstat)+I*randn(rstat))*spect->p[i];
    }
    UNLOCK(data->mutex_ilayer);
    cfft2(cspect,-1); 	
    double wt1=sqrt(wt[ilayer]);
    double *p=screen[ilayer]->p;

    for(long i=0; i<m*n; i++){
	p[i]=creal(cspect->p[i])*wt1;
    }
    
    if(ilayer+1<nlayer){
	double *q=screen[ilayer+1]->p;
	double wt2=sqrt(wt[ilayer+1]);
	for(long i=0; i<m*n; i++){
	    q[i]=cimag(cspect->p[i])*wt2;
	}
    }
    goto repeat;
}

/**
 * Generates multiple screens from spectrum.
 */
map_t** genscreen_from_spect(rand_t *rstat, dmat *spect, double dx,double r0, double L0,
			     double* wt, int nlayer, int nthread){
    
    GENSCREEN_T screendata;
    memset(&screendata, 0, sizeof(screendata));
    screendata.rstat=rstat;
    screendata.spect=spect;
    screendata.wt=wt;
    screendata.nlayer=nlayer;
    screendata.ilayer=0;
    PINIT(screendata.mutex_ilayer);
    spect->p[0]=0;//zero out piston.
    long nx=spect->nx;
    long ny=spect->ny;
    int inshm=-1;
    long totmem=(nx*ny*nlayer+1)*sizeof(double);//+1 for sanity check.
    int fd;
    map_t **screen=screendata.screen=atmnew_shm(&fd, &inshm, rstat, nx, ny, 
						dx, r0, L0, wt, nlayer, 1);
    if(inshm != 1){
	CALL(genscreen_do, &screendata, nthread);
	if(inshm == 2){//need to remap to r/o mode.
	    screen[0]->p[nx*ny*nlayer]=1;//say the data is valid.
	    if(munmap(screen[0]->p, totmem)){//unmap read/write. remap as read only later
		error("munmap failed\n");
	    }
	    //the exclusive lock will be converted to shared lock automatically.
	    map_shm(screen, totmem, nlayer, fd, 0);
	}
    }
    return screen;
}
/**
 *   Generate vonkarman screens from turbulence statistics.
 */
map_t** vonkarman_screen(ATM_ARGS){
    (void)ninit;
    char fnspect[PATH_MAX];
    dmat *spect;
    mymkdir("%s/.aos/spect/",HOME);
    snprintf(fnspect,PATH_MAX,"%s/.aos/spect/spect_%ldx%ld_dx1_%g_r0%g_L0%g.bin",
	     HOME,nx,ny,1./dx,r0,L0);
    if(exist(fnspect)){
	spect=({char fn[PATH_MAX];strcpy(fn,fnspect);dread("%s",fn);});
    }else{
	info2("\nGenerating spect...");
	TIC;
	tic;
	spect=turbpsd(nx,ny,dx,r0,L0,0.5);
	toc2("done");
	dwrite(spect,"%s",fnspect);
    }
    map_t **screen=genscreen_from_spect(rstat,spect,dx,r0,L0,wt,nlayer,nthread);
    dfree(spect);
    return(screen);
}

/**
 *  Generate screens from PSD with power of 12/3 instead of 11/3.
 */
map_t** biharmonic_screen(ATM_ARGS){
    (void)ninit;
    /**
       Generate vonkarman screen.
    */
    dmat *spect;
  
    info2("\nGenerating spect...");
    TIC;
    tic;
    spect=turbpsd_full(nx,ny,dx,r0,L0,-2,0.5);
    toc2("done");

    map_t **screen=genscreen_from_spect(rstat,spect,dx,r0,L0,wt,nlayer,nthread);
    dfree(spect);
    return(screen);
}
static void fractal_screen_do(GENSCREEN_T *data){
    rand_t *rstat=data->rstat;
    map_t** screen=data->screen;
    const double *wt=data->wt;
    long nx=screen[0]->nx;
    long ny=screen[0]->ny;
 repeat:
    LOCK(data->mutex_ilayer);
    int ilayer=data->ilayer;
    data->ilayer++;
    if(ilayer>=data->nlayer){
	UNLOCK(data->mutex_ilayer);
	return;
    }
    for(long i=0; i<nx*ny; i++){
	screen[ilayer]->p[i]=randn(rstat);
    }
    UNLOCK(data->mutex_ilayer);
    double r0i=data->r0*pow(wt[ilayer], -3./5.);
    //info("r0i=%g\n", r0i);
    fractal(screen[ilayer]->p, nx, ny, screen[0]->dx, r0i, data->L0, data->ninit);
    remove_piston(screen[ilayer]->p, nx*ny);
    goto repeat;
}

/**
 * Generate Fractal screens. Not good statistics.
 */

map_t **fractal_screen(ATM_ARGS){
    GENSCREEN_T screendata;
    memset(&screendata, 0, sizeof(screendata));
    screendata.rstat=rstat;
    screendata.wt=wt;
    screendata.nlayer=nlayer;
    screendata.ilayer=0;
    screendata.r0=r0;
    screendata.L0=L0;
    screendata.ninit=ninit;
    PINIT(screendata.mutex_ilayer);
    int inshm=-1;
    long totmem=(nx*ny*nlayer+1)*sizeof(double);
    int fd;
    map_t **screen=screendata.screen=atmnew_shm(&fd, &inshm, rstat, nx, ny, 
						dx, r0, L0, wt, nlayer, 1);
    if(inshm != 1){
	CALL(fractal_screen_do,&screendata, nthread);
	if(inshm == 2){
	    screen[0]->p[nx*ny*nlayer]=1;//say the data is valid.
	    if(munmap(screen[0]->p, totmem)){//unmap read/write. remap as read only later
		error("munmap failed\n");
	    }
	    //the exclusive lock will be converted to shared lock automatically.
	    map_shm(screen, totmem, nlayer, fd, 0);
	}
    }
    return screen;
}

/**
 *  Compute the covariance for separation of r, and put the values in cov. In
 *  kolmogorov spectrum, the variance are defined as half of the structure
 *  function between two points separated by rmax.
 */

dmat* turbcov(dmat *r, double rmax, double r0, double L0){
    double tg1=tgamma(11./6) * pow(24./5 * tgamma(6./5.), 5./6.)
	* pow(2 * M_PI/0.5e-6, -2) / pow(M_PI, 8./3.);
    double vkcoeff  = tg1 / pow(2, 5./6.);    
    double vkcoeff0 = tg1 * tgamma(5./6.) / 2 ;//for variance
    dmat *cov=dnew(r->nx, r->ny);
    long n=r->nx*r->ny;
    if(isinf(L0)){//kolmogorov.
	const double power=5./3.;
	double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
	double sigma2=0.5*coeff*pow(rmax, power);
	for(long i=0; i<n; i++){
	    cov->p[i]=sigma2-0.5*coeff*pow(r->p[i], power);
	}
    }else{//von karman.
	const double f0=1./L0;
	const double r0f0p=pow(r0*f0, -5./3.);
	double ri, rk, rip, rkp;
	double r2pif0;	
	for(long i=0; i<n; i++){
	    if(fabs(r->p[i])<EPS){
		cov->p[i]=vkcoeff0*r0f0p;
	    }else{
		r2pif0=r->p[i]*2*M_PI*f0;
		bessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
		cov->p[i]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	    }
	}
    }
    return cov;
}

/**
 * Compute the turbulence spectrum at size nx*ny, with spacing dx. Notice that
 * the zero frequency component is in the corner psd->p[0].
 */

dmat *turbpsd_full(long nx,      /**<The size*/
		   long ny,      /**<The size*/
		   double dx,    /**<The sampling of spatial coordinate.*/
		   double r0,    /**<The Fried parameter*/
		   double L0,    /**<The outer scale*/
		   double slope, /**<should be -11/6 for von karman or kolmogorov
				    screens, or -2 for biharmonic screen (just
				    testing only).*/
		   double power  /**< optionally do a power of psd.*/
		   ){
    if(nx & 1 || ny & 1){
	warning("Screen is odd size.");
    }
    slope*=power;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double dfx2=dfx*dfx;
    const double L02=pow(L0,-2);
    const double scrnstr=pow(0.0229*pow(r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2)*(dfx*dfy),power);
    const int nx2=nx/2;
    const int ny2=ny/2;
    dmat *psd=dnew(nx,ny);
    for(int i=0;i<ny;i++){
	double r2y=pow((i<ny2?i:i-ny)*dfy,2);// to avoid fft shifting.
	double *psd1=psd->p+i*nx;
	for(int j=0;j<nx2;j++){
	    double r2x=j*j*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
	for(int j=nx2;j<nx;j++){
	    double r2x=(j-nx)*(j-nx)*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
    }
    return psd;
}

