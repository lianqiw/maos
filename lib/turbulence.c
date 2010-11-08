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
#include "common.h"
#include "turbulence.h"
#include "path.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
#include "fft.h"
#include "hashlittle.h"
int disable_atm_shm=0;
/**
   \file turbulence.c
   Contains routines to generate atmospheric turbulence screens
*/
/**
   Create a vonkarman spectrum. This is the square root of PSD. Good for
   rectangular screen. piston is in (0,0)
*/
dmat *vonkarman_spect(int nx, int ny, double dx, double r0, double L0){
    if(nx & 1 || ny & 1){
	error("Screen is odd size.");
    }

    const double power=-11./12.;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double dfx2=dfx*dfx;
    const double L02=pow(L0,-2);
    const double scrnstr=sqrt(0.0229)*pow(r0,-5./6.)*(0.5e-6)/(2.*M_PI)*sqrt(dfx*dfy);
    const int nx2=nx/2;
    const int ny2=ny/2;
    dmat *spect=dnew(nx,ny);
    for(int i=0;i<ny;i++){
	double r2y=pow((i<ny2?i:i-ny)*dfy,2);// to avoid fft shifting.
	double *spect1=spect->p+i*nx;
	for(int j=0;j<nx2;j++){
	    double r2x=j*j*dfx2;
	    spect1[j] = pow(r2x+r2y+L02,power)*scrnstr;
	}
	for(int j=nx2;j<nx;j++){
	    double r2x=(j-nx)*(j-nx)*dfx2;
	    spect1[j] = pow(r2x+r2y+L02,power)*scrnstr;
	}
    }
    //spect->p[0]=0;//piston is in corner
    return spect;
}

/**
   inverse of PSD. zero frequency in lower left corner.
*/
dmat *vonkarman_invpsd(int nx, int ny, double dx, double r0, double L0){
    dmat *psd=vonkarman_spect(nx, ny, dx, r0, L0);
    double scale=pow((double)(nx*ny),-2);
    for(long i=0; i<psd->nx*psd->ny; i++){
	psd->p[i]=pow(psd->p[i],-2)*scale;
    }
    return psd;
}
typedef struct GENSCREEN_T{
    rand_t *rstat;
    dmat *spect;
    map_t **screen;
    double *wt;
    int nlayer;
    int ilayer;
    pthread_mutex_t mutex_ilayer;
}GENSCREEN_T;
/**
   Generate turbulence screens.
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
   Generates multiple screens from spectrum.
*/
map_t** genscreen_from_spect(rand_t *rstat, dmat *spect, double dx,
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
    long m=spect->nx;
    long n=spect->ny;
    map_t** screen;

    screen=calloc(nlayer,sizeof(map_t*));
    screendata.screen=screen; 
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	screen[ilayer]=calloc(1, sizeof(map_t));
	screen[ilayer]->nx=m;
	screen[ilayer]->ny=n;
	screen[ilayer]->ox=-m*dx/2;
	screen[ilayer]->oy=-n*dx/2;
	screen[ilayer]->dx=dx;
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
	key=hashlittle(spect->p, m*n*sizeof(double), key);
	key=hashlittle(wt, nlayer*sizeof(double), key);
	snprintf(fnshm,NAME_MAX,"/maos_atm_%ud_%d_%ld_%g",key,nlayer,m,dx);
	long already_exist = shm_free_unused(fnshm, 0);
	totmem=(m*n*nlayer+1)*sizeof(double);//+1 for sanity check.
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
    retry:
	//sleep so that the creation process has enough time to make an exclusive lock.
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
		info2("\nCreating %s ...",fnshm);
		/*apply an exclusive lock so no other process can make a shared
		  lock thus preventing them from reading the data.*/
		if(flock(fd, LOCK_EX)){
		    error("Failed to apply an exclusive lock");
		}
		if(ftruncate(fd, totmem)){//allocate the size.
		    perror("ftruncate");
		    error("Failed to allocate memory. check kernel.shmmax in /etc/sysctl.conf");
		}
		fchmod(fd, 00777);//make it globally writable. (not working as what I wanted)
		if((screen[0]->p=mmap(NULL, totmem, PROT_READ|PROT_WRITE, MAP_SHARED ,fd, 0))<0){
		    error("Unable to mmap for read/write\n");
		}
		screen[0]->shm=fd;//use the file descriptor so we can release the lock.
		for(int ilayer=1; ilayer<nlayer; ilayer++){
		    screen[ilayer]->shm = -1;
		    screen[ilayer]->p   = screen[0]->p+m*n*ilayer;
		}
		screen[0]->p[m*n*nlayer]=0;//say the data is not valid.
		CALL(genscreen_do, &screendata, nthread);
		screen[0]->p[m*n*nlayer]=1;//say the data is valid.
		if(munmap(screen[0]->p, totmem)){//unmap read/write. remap as read only later
		    error("munmap failed\n");
		}
		flock(fd, LOCK_UN);//release exclusive lock.
		if(close(fd)){//close file descriptor
		    warning("Error closing %d\n",fd);
		}
		//Reopen as read only and do shared lock.
		if((fd=shm_open(fnshm, O_RDONLY, 00777))<0){
		    error("Failed to open shared segment.\n");
		}
	    }
	}else{
	    info2("\nReusing %s ...",fnshm);
	}
	//fd might be zero.
	futimes(fd, NULL);//set access, modification time to current.
	if(flock(fd, LOCK_SH)){
	    error("Failed to apply a shared lock on shared segment.\n");
	}
	//we map it in read only mode.
	screen[0]->p=mmap(NULL, totmem, PROT_READ, MAP_SHARED, fd, 0);
	if(fabs(screen[0]->p[m*n*nlayer])<1.e-15){
	    warning2("Data in shm %s is invalid, free and redo it.\n", fnshm);
	    flock(fd, LOCK_UN);//release lock.
	    close(fd);//close descriptor.
	    if(munmap(screen[0]->p, totmem)){//unmap.
		error("munmap failed\n");
	    }
	    shm_unlink(fnshm);//destroy the map.
	    goto retry;
	}
	if(screen[0]->p<0){
	    error("Unable to mmap for read/write\n");
	}								
	//since fd might be 0. we add fd by 1 and assign to shm.
	screen[0]->shm=fd+1;//save file descriptor of shared mem for easy munmap later and close fd
	for(int ilayer=1; ilayer<nlayer; ilayer++){
	    screen[ilayer]->shm = -1;//say this is a part of the shared mem
	    screen[ilayer]->p = screen[0]->p+m*n*ilayer;
	}  
	/*In the end, a shared lock is placed on the shm
	  memory. if the screen is freed or program exited, the
	  shared lock will be released. THis is the method I
	  invented to tell how many process is using this shm.*/
    }else{
#endif
	for(int ilayer=0; ilayer<nlayer; ilayer++){
	    screen[ilayer]->p=malloc(m*n*sizeof(double));
	    screen[ilayer]->shm=0;
	}
	CALL(genscreen_do, &screendata, nthread);
#if USE_POSIX_SHM == 1
    }
#endif
    return screen;
}
/**
   Generate vonkarman screens from turbulence statistics.
*/
map_t** vonkarman_screen(rand_t *rstat, int m, int n, double dx, 
			   double r0, double L0, double* wt, int nlayer, int nthread){
   
    char fnspect[PATH_MAX];
    dmat *spect;
    mymkdir("%s/.aos/spect/",HOME);
    snprintf(fnspect,PATH_MAX,"%s/.aos/spect/spect_%dx%d_dx1_%g_r0%g_L0%g.bin",
	     HOME,m,n,1./dx,r0,L0);
    if(exist(fnspect)){
	spect=({char fn[PATH_MAX];strcpy(fn,fnspect);dread("%s",fn);});
    }else{
	info2("\nGenerating spect...");
	TIC;
	tic;
	spect=vonkarman_spect(m,n,dx,r0,L0);
	toc2("done");
	dwrite(spect,"%s",fnspect);
    }
    map_t **screen=genscreen_from_spect(rstat,spect,dx,wt,nlayer,nthread);
    dfree(spect);
    return(screen);
}

/**
   Create a vonkarman spectrum. This is the square root of PSD. Good for
   rectangular screen. piston is in (0,0)
*/
static dmat *biharmonic_spect(int nx, int ny, double dx, double r0, double L0){
  
    if(nx & 1 || ny & 1){
	error("Screen is odd size.");
    }

    const double power=-1.;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double dfx2=dfx*dfx;
    const double L02=pow(L0,-2);
    const double scrnstr=sqrt(0.0229)*pow(r0,-5./6.)*(0.5e-6)/(2.*M_PI)*sqrt(dfx*dfy);
    const int nx2=nx/2;
    const int ny2=ny/2;
    dmat *spect=dnew(nx,ny);
    for(int i=0;i<ny;i++){
	double r2y=pow((i<ny2?i:i-ny)*dfy,2);// to avoid fft shifting.
	double *spect1=spect->p+i*nx;
	for(int j=0;j<nx2;j++){
	    double r2x=j*j*dfx2;
	    spect1[j] = pow(r2x+r2y+L02,power)*scrnstr;
	}
	for(int j=nx2;j<nx;j++){
	    double r2x=(j-nx)*(j-nx)*dfx2;
	    spect1[j] = pow(r2x+r2y+L02,power)*scrnstr;
	}
    }
    //spect->p[0]=0;//piston is in corner
    return spect;
}
/**
   Generate screens from PSD with power of 12/3 instead of 11/3.
*/
map_t** biharmonic_screen(rand_t *rstat, int m, int n, double dx, 
			    double r0, double L0, double* wt, int nlayer, int nthread){
    /**
       Generate vonkarman screen.
    */

    dmat *spect;
   
    info2("\nGenerating spect...");
    TIC;
    tic;
    spect=biharmonic_spect(m,n,dx,r0,L0);
    toc2("done");

    map_t **screen=genscreen_from_spect(rstat,spect,dx,wt,nlayer, nthread);
    dfree(spect);
    return(screen);
}
