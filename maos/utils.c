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

#include <math.h>
#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <signal.h>
#include <fcntl.h>           /* For O_* constants */
#include <errno.h>
#include <getopt.h>
#include "maos.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

char *dirsetup=NULL;
char *dirskysim=NULL;
/**
   \file maos/utils.c
   A few utility routines
 */
/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
out. pcalib2 part of bkgrnd2 is calibrated out.  */
void addnoise(dmat *A,              /**<The pixel intensity array*/
	      rand_t* rstat,   /**<The random stream*/
	      const double bkgrnd,  /**<background in PDEs per pixel per frame*/
	      const double pcalib,  /**<Fraction of bkgrnd that is calibrated out*/
	      const double *bkgrnd2, /**<background in PDEs of each pixel per frame.*/
	      const double *bkgrnd2c,/**<calibration bkgrnd2*/
	      const double rne      /**<Read out noise per pixel per read*/
	      ){
 
    if(bkgrnd2){
	for(int ix=0; ix<A->nx*A->ny; ix++){
	    A->p[ix]=randp(rstat,A->p[ix]+bkgrnd+bkgrnd2[ix]) +rne*randn(rstat)
		-bkgrnd*pcalib ;
	}
    }else{
	for(int ix=0; ix<A->nx*A->ny; ix++){
	    A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)
			   -bkgrnd*pcalib+rne*randn(rstat);
	}
    }
    if(bkgrnd2c){
	for(int ix=0; ix<A->nx*A->ny; ix++){
	    A->p[ix]-=bkgrnd2c[ix];
	}
    }
}
/**
   Calls create_metapupil is simplified interface by returning a map_t object.
 */
map_t * create_metapupil_wrap(const PARMS_T *parms, double ht,double dx,
			      double offset,double guard, long ninx,long niny,
			      int pad,int square){
    long nx, ny;
    double ox, oy, *p;
    create_metapupil(parms,ht,dx,offset,&(nx),&(ny),
		     &(ox),&(oy),&(p),guard, ninx, niny, pad,square);
    map_t *amp=mapnew(nx, ny, dx, p);
    amp->ox=ox;
    amp->oy=oy;
    amp->h=ht;
    return amp;
}
/**
   create a metapupil map, with size nx*ny, origin at (ox,oy), sampling of dx,
   height of ht, that can cover all the WFS and science beams.  

   offset: distance in pixel from the point closest to the origin to origin (right
   side).  
   0: there is a point on the origin.  
   1/2: the closest point to origin is 1/2 pixel.  
   
   pad: 1: round nx, ny to power of 2.  */

void create_metapupil(const PARMS_T *parms, double ht,double dx,
		      double offset, long* nxout, long* nyout,
		      double *oxout, double *oyout, double**map, 
		      double guard, long ninx, long niny, int pad,int square){
    double R=parms->aper.d/2;
    double maxx=0,maxy=0;
    double sx,sy;/*temporary variables */
    int i;
    int use_wfs_hi=1;
    int use_wfs_lo=1;
    int use_evl=1;
    int use_fit=1;
   
    /*find minimum map size to cover all the beams */
    for(i=0; i<parms->nwfs; i++){
	int ipowfs=parms->wfs[i].powfs;
	if((parms->powfs[ipowfs].lo && !use_wfs_lo) 
	   ||(!parms->powfs[ipowfs].lo && !use_wfs_hi)){
	    continue;
	}
	sx=fabs(parms->wfs[i].thetax*ht)
	    +(1.-ht/parms->powfs[ipowfs].hs)*R+guard;
	sy=fabs(parms->wfs[i].thetay*ht)
	    +(1.-ht/parms->powfs[ipowfs].hs)*R+guard;
	if(sx>maxx) maxx=sx;
	if(sy>maxy) maxy=sy;
    }
    if(use_evl){
	for(i=0; i<parms->evl.nevl; i++){
	    sx=fabs(parms->evl.thetax[i]*ht)+R+guard;
	    sy=fabs(parms->evl.thetay[i]*ht)+R+guard;
	    if(sx>maxx) maxx=sx;
	    if(sy>maxy) maxy=sy;
	}
    }
    if(use_fit){
	for(i=0; i<parms->fit.nfit; i++){
	    sx=fabs(parms->fit.thetax[i]*ht)+R+guard;
	    sy=fabs(parms->fit.thetay[i]*ht)+R+guard;
	    if(sx>maxx) maxx=sx;
	    if(sy>maxy) maxy=sy;
	}
    }
    /*normalized grid size +3 is extra guard band */
    maxx=ceil(maxx/dx)+2;
    maxy=ceil(maxy/dx)+2;
    long nx,ny;
    nx=iceil(maxx+offset+1)*2;
    ny=iceil(maxy+offset+1)*2;
    if(pad){/*pad to power of 2 */
	nx=1<<iceil(log2((double)nx));
	ny=1<<iceil(log2((double)ny));
    }
    /*Make it square */
    nx=(nx<ny)?ny:nx;
    ny=nx;
    if(ninx>1){
	if(ninx<nx) warning("ninx=%ld is too small\n",ninx);
	nx=ninx;
    }
    if(niny>1){
	if(niny<ny)  warning("niny=%ld is too small\n",niny);
	ny=niny;
    }
    double ox,oy;
    ox=((nx/2)-offset);
    oy=((ny/2)-offset);

    if(nxout)
	*nxout=nx;
    if(nyout)
	*nyout=ny;
    if(oxout)
	*oxout=-ox*dx;
    if(oyout)
	*oyout=-oy*dx;

    if(map && square){
	dmat *dmap=dnew(nx,ny);
	*map=dmap->p;
	dset(dmap,1);
	dfree_keepdata(dmap);
    }else if(map){
	/*guard+=0.707;//necessary with dcircle_deprecated. not with current dcircle */
	dmat *dmap=dnew(nx,ny);
	*map=dmap->p;
	double Rn=R/dx;
	double Rg=guard/dx;
	for(i=0; i<parms->nwfs; i++){
	    int ipowfs=parms->wfs[i].powfs;
	    if((parms->powfs[ipowfs].lo && !use_wfs_lo) 
	       ||(!parms->powfs[ipowfs].lo && !use_wfs_hi)){
		continue;
	    }
	    sx=ox+(parms->wfs[i].thetax*ht)/dx;
	    sy=oy+(parms->wfs[i].thetay*ht)/dx;
	    double RR=Rn*(1.-ht/parms->powfs[ipowfs].hs)+Rg;
	    dcircle_symbolic(dmap,sx,sy,RR);
	}

	if(use_evl){
	    for(i=0; i<parms->evl.nevl; i++){
		sx=ox+(parms->evl.thetax[i]*ht)/dx;
		sy=oy+(parms->evl.thetay[i]*ht)/dx;
		double RR=Rn+Rg;
		dcircle_symbolic(dmap,sx,sy,RR);
	    }
	}
	if(use_fit){
	    for(i=0; i<parms->fit.nfit; i++){
		sx=ox+(parms->fit.thetax[i]*ht)/dx;
		sy=oy+(parms->fit.thetay[i]*ht)/dx;
		double RR=Rn+Rg;
		dcircle_symbolic(dmap,sx,sy,RR);
	    }
	}

	for(i=0; i<nx*ny; i++){
	    dmap->p[i]=(dmap->p[i])>1.e-15?1:0;
	}
	dfree_keepdata(dmap);
    }
}
/**
   Plot the projection of all the different FoVs on a certain grid. used in setup_recon
 */
void plotloc(char *fig, const PARMS_T *parms, 
	     loc_t *loc, double ht, char *format,...){
    format2fn;
    int ncir=parms->evl.nevl + parms->nwfs;
    double (*cir)[4];
    cir=(double(*)[4])calloc(ncir*4,sizeof(double));
    int count=0;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	cir[count][0]=ht*parms->evl.thetax[ievl];
	cir[count][1]=ht*parms->evl.thetay[ievl];
	cir[count][2]=parms->aper.d*0.5;
	cir[count][3]=0xFF0000;/*rgb color */
	count++;
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	double hs=parms->powfs[parms->wfs[iwfs].powfs].hs;
	cir[count][0]=parms->wfs[iwfs].thetax*ht;
	cir[count][1]=parms->wfs[iwfs].thetay*ht;
	cir[count][2]=parms->aper.d*0.5*(1.-ht/hs);
	if(isinf(hs)){
	    cir[count][3]=0x44FF00;/*rgb color */
	}else{
	    cir[count][3]=0xFFFF33;/*rgb color */
	}
	count++;
    }
    plot_points(fig, 1, &loc, NULL ,NULL, NULL,ncir, cir, NULL,
	       "Coordinate","x (m)","y (m)", "%s",fn);
    free(cir);
}
/**
   ploted all the different beam directions as points. used in setup_parms */
void plotdir(char *fig, const PARMS_T *parms, double totfov, char *format,...){
    format2fn;
    int ncir=1;
    double (*cir)[4];
    cir=(double(*)[4])calloc(ncir*4,sizeof(double));
    cir[0][0]=0;
    cir[0][1]=0;
    cir[0][2]=totfov/2;
    cir[0][3]=0x000000;/*rgb color */
    int ngroup=2+parms->npowfs;
    loc_t **locs=calloc(ngroup, sizeof(loc_t*));
    int32_t *style=calloc(ngroup, sizeof(int32_t));

    style[0]=(0xFF0000<<8)+(4<<4)+3;
    locs[0]=locnew(parms->evl.nevl, 0);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	locs[0]->locx[ievl]=parms->evl.thetax[ievl]*206265;
	locs[0]->locy[ievl]=parms->evl.thetay[ievl]*206265;
    }

    style[1]=(0xFF22DD<<8)+(4<<4)+3;
    locs[1]=locnew(parms->fit.nfit, 0);
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	locs[1]->locx[ifit]=parms->fit.thetax[ifit]*206265;
	locs[1]->locy[ifit]=parms->fit.thetay[ifit]*206265;
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	locs[ipowfs+2]=locnew(parms->powfs[ipowfs].nwfs, 0);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs[jwfs];
	    locs[ipowfs+2]->locx[jwfs]=parms->wfs[iwfs].thetax*206265;
	    locs[ipowfs+2]->locy[jwfs]=parms->wfs[iwfs].thetay*206265;
	}
	if(!isinf(parms->powfs[ipowfs].hs)){
	    style[ipowfs+2]=(0xFF8800<<8)+(4<<4)+2;
	}else if(!parms->powfs[ipowfs].lo){
	    style[ipowfs+2]=(0xFFFF00<<8)+(4<<4)+1;
	}else if(parms->powfs[ipowfs].order>1){
	    style[ipowfs+2]=(0x0000FF<<8)+(4<<4)+4;
	}else{
	    style[ipowfs+2]=(0x0000FF<<8)+(4<<4)+1;
	}
    }
    double limit[4];
    limit[0]=limit[2]=-totfov/2;
    limit[1]=limit[3]=totfov/2;
    plot_points(fig, ngroup, locs, NULL, style,limit,ncir,cir, NULL,
		"Asterism","x (arcsec)", "y (arcsec)", "%s",fn);
    free(cir);
    locarrfree(locs, ngroup);
    free(style);
}
/**
   Rename the log files when simulation exits.
 */
void rename_file(int sig){
    char fnnew[256];
    char fnold[256];
    char suffix[16];
    int pid=getpid();
    switch(sig){
    case 0:
	sprintf(suffix,"done");
	break;
    case SIGBUS:
    case SIGILL:
    case SIGSEGV:
    case SIGABRT:
	sprintf(suffix,"err");
	break;
    case SIGKILL:
    case SIGINT: /*Ctrl-C */
    case SIGTERM:
    case SIGQUIT: /*Ctrl-'\' */
	sprintf(suffix,"killed");
	break;
    default:
	sprintf(suffix,"unknown");
    }
    snprintf(fnnew,256,"kill_%d",pid);
    if(exist(fnnew)) remove(fnnew);
    
    snprintf(fnold,256,"run_%d.log",pid);
    snprintf(fnnew,256,"run_%d.%s", pid,suffix);
    rename(fnold,fnnew);
    mysymlink(fnnew, "run_recent.log");

    if(curparms && sig!=0){
	char fn[80];
	const PARMS_T *parms=curparms;
	for(int iseed=curiseed; iseed<parms->sim.nseed; iseed++){
	    if(parms->fdlock[iseed]>=0){
		close(parms->fdlock[iseed]);
		int seed=parms->sim.seeds[iseed];
		snprintf(fn, 80, "Res_%d.lock",seed);
		(void) remove(fn);
	    }
	}
    }
}
/**
   Handles signals. We don't want to exit the simulation when SIGPIPE happens
   (when we are writing to closed sockets)
 */
void maos_signal_handler(int sig){
    if(sig==SIGPIPE){
	warning3("Program received signal SIGPIPE, broken pipe.\n");
	return;
    }
    disable_signal_handler;
    rename_file(sig);/*handles signal */
    if(sig!=0){
	char *info="Unknown";
	switch(sig){
	case SIGBUS:
	case SIGILL:
	case SIGSEGV:
	case SIGABRT:
	    info="Segmentation falt";
	    break;
	case SIGKILL:
	case SIGINT: /*Ctrl-C */
	case SIGTERM:
	case SIGQUIT: /*Ctrl-'\' */
	    info="Killed";
	case SIGUSR1:/*user defined */
	    info="Quit";
	    break;
	}
	warning2("Signal %d: %s\n", sig, info);
	fflush(stderr);
	fflush(stdout);
	if(sig !=0){
	    print_backtrace(0);
	}
	scheduler_finish(1);
	_Exit(sig);/*don't call clean up functions, but does flush files. */
    }
}
/**
   Print out usage information.
 */
static void print_usage(void){
    info2(
"Usage: maos [OPTION...] [FILE]...\n"
"maos is a simulation tool developed to simulate scao/mcao systems\n\n"
"Examples:\n"
"maos   # Run the default configuration of NFIRAOS: nfiaros.conf as the baseline.\n"
"maos -c scao_ngs.conf -s 2 -n 2 -d -o scao_ngs override.conf chol.conf\n"
"       # Run a single conjugate natural guide star case, with seed 2, 2 threads\n"
"       # detach from the terminal and output results to folder scao_ngs\n"
"       # and read in overriding parameters stored in override.conf and chol.conf\n"
"\n"
"Options: \n"
"-h, --help        to print out this message\n"
"-d, --detach      to detach from terminal and run in background\n"
"-f, --force       force starting simulation without scheduler\n"
"-n, --nthread=N   Use N threads, default is 1\n"
"-s, --seed=N      Use seed N instead of the numbers specified in .conf files\n"
"-o, --output=DIR  output the results to DIR.\n"
"-c, --conf=FILE.conf\n"
"                  Use FILE.conf as the baseline config instead of nfiraos.conf\n"
"-p, --path=dir    Add dir to the internal PATH\n"
"-P, --pause       paulse simulation in the end of every time step\n"
"-g, --gpu=N       Use the N'th gpu. 0 for the first. -1 to disable. default: automatic"
	  );
    exit(0);
}
  
#if USE_STATIC
static __attribute__((destructor)) void deinit(){
    char tmp[PATH_MAX];
    snprintf(tmp, PATH_MAX, "rm -rf %s/config-%ld/", TEMP, (long)getpid());
    if(system(tmp)){}
}
#endif

/**
   Parse command line arguments argc, argv
 */
ARG_T * parse_args(int argc, char **argv){
    ARG_T *arg=calloc(1, sizeof(ARG_T));
    int *seeds=NULL; int nseed=0;
    ARGOPT_T options[]={
	{"help",   'h',T_INT, 2, print_usage, NULL},
	{"detach", 'd',T_INT, 0, &arg->detach, NULL},
	{"force",  'f',T_INT, 0, &arg->force, NULL},
	{"output", 'o',T_STR, 1, &arg->dirout, NULL},
	{"nthread",'n',T_INT, 1, &arg->nthread,NULL},
	{"gpu",    'g',T_INTARR, 1, &arg->gpus, &arg->ngpu},
	{"conf",   'c',T_STR, 1, &arg->conf, NULL},
	{"seed",   's',T_INTARR, 1, &seeds, &nseed},
	{"path",   'p',T_STR, 3, addpath, NULL},
	{"pause",  'P',T_INT, 0, &arg->pause, NULL},
	{NULL, 0,0,0, NULL, NULL}
    };
    char *cmds=parse_argopt(argc, argv, options);
    if(arg->nthread>NCPU2 || arg->nthread<=0){
	arg->nthread=NCPU2;
    }
#if USE_PTHREAD == 0
    arg->nthread=1;
#endif
    char fntmp[PATH_MAX];
    snprintf(fntmp,PATH_MAX,"%s/maos_%ld.conf",TEMP,(long)getpid());
    FILE *fptmp=fopen(fntmp,"w");
    fputs(cmds, fptmp);
    free(cmds); cmds=NULL;
    if(nseed){
	fprintf(fptmp, "\nsim.seeds=[");
	for(int iseed=0; iseed<nseed; iseed++){
	    fprintf(fptmp, "%d ", seeds[iseed]);
	}
	fprintf(fptmp, "]\n");
	free(seeds); seeds=NULL;
    }
    fclose(fptmp);
    arg->confcmd=strdup(fntmp);
    
    if(arg->pause && arg->detach){
	warning("-p and -d are both specified, disable -d\n");
	arg->detach=0;
    }

    if(!arg->conf){ /*If -c is not specifid in path, will use default.conf*/
	arg->conf=strdup("default.conf");
    }
    /*Setup PATH and result directory so that the config_path is in the back of path */
    char *config_path=find_config("maos");
#if USE_STATIC && 0
    if(!exist(config_path)){
	/*We are binding the config files in the executable, extract and use it. */
	char *start=&_binary____config_tar_gz_start;
	char *end=&_binary____config_tar_gz_end;
	long len=end-start;
	assert(len>0);
	char fn[PATH_MAX];
	mymkdir("%s/config-%d", TEMP, (int)getpid());
	snprintf(fn, PATH_MAX, "%s/config-%d/config.tar.gz", TEMP, (int)getpid());
	FILE *fp=fopen(fn, "w");
	if(fp){
	    fwrite(start, 1, end-start, fp);
	    fclose(fp);
	    snprintf(fn, PATH_MAX, "cd %s/config-%d/ && tar xf config.tar.gz", TEMP, (int)getpid());
	    if(system(fn)<0){
		perror("system");
		warning("Unable to execute %s\n", fn);
	    }else{
		snprintf(fn, PATH_MAX, "%s/config-%d/config", TEMP, (int)getpid());
		config_path=stradd(fn,"/maos",NULL);
	    }
	}else{
	    warning("Unable to open %s\n", fn);
	}
    }
#endif
    if(!config_path || !exist(config_path)){
	error("Unable to find usable configuration file\n");
    }
    /*info2("Using config files found in %s\n", config_path); */
    char *bin_path=stradd(config_path, "/bin", NULL);
    addpath(config_path);
    addpath(bin_path);
    free(bin_path);
    free(config_path);

    if(!arg->dirout){
	arg->dirout=strtime();
    }
    addpath(".");
    mymkdir("%s",arg->dirout);
    if(chdir(arg->dirout)){
	error("Unable to chdir to %s\n", arg->dirout);
    }

    return arg;
}
/**
   Computes strehl from OPD without doing FFT. The strehl is simply 

   \f$s=\sum(A*exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]) \f$

   where A is the amplitude map.
 */
cmat *strehlcomp(const dmat *iopdevl, const double *amp, const double wvl){
    dcomplex i2pi=I*2*M_PI/wvl;
    dcomplex strehl=0;
    for(int iloc=0; iloc<iopdevl->nx; iloc++){
	strehl+=amp[iloc]*cexp(i2pi*iopdevl->p[iloc]);
    }
    cmat *psf2=cnew(1,1);
    psf2->p[0]=strehl;
    return psf2;
}
typedef struct PSFCOMP_T{
    ccell *psf2s;/*output psf. */
    const dmat *iopdevl;
    const double *amp;
    long **embeds;
    const long *nembeds;
    const int *psfsize;
    const int nwvl;
    const double *wvl;
}PSFCOMP_T;
/**
   Call psfcomp_iwvl in parallel to compute PSF in each wvl.
*/
ccell *psfcomp(const dmat *iopdevl, const double *amp,
	       long **embeds, const long *nembeds, const int *psfsize,
	       const int nwvl, const double *wvl){
    ccell *psf2s=ccellnew(nwvl, 1);
    PSFCOMP_T data={psf2s, iopdevl, amp, embeds, nembeds, psfsize, nwvl, wvl};
    thread_t dopsf[nwvl];
    thread_prep(dopsf, 0, nwvl, nwvl, (thread_wrapfun)psfcomp_iwvl, &data);
    CALL_THREAD(dopsf, nwvl, 0);
    return psf2s;
}
/**
   Computes PSF from OPD by doing FFT. The PSF is computed as

   \f$\textrm{PSF}=\mathcal{F}[A\times exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]]\f$

   The peak value (center) in the computed PSF is normalized by the peak value
   in the differaction limited PSF. In other words, the peak value in the
   computed PSF is the Strehl. Keep this in mind when you compute enclosed
   energy.
*/
void psfcomp_iwvl(thread_t *thdata){
    PSFCOMP_T *data=thdata->data;
    ccell *psf2s=data->psf2s;
    const dmat *iopdevl=data->iopdevl;
    const double *amp=data->amp;
    long **embeds=data->embeds;
    const long *nembeds=data->nembeds;
    const int *psfsize=data->psfsize;
    const double *wvl=data->wvl;
    
    for(int iwvl=thdata->start; iwvl<thdata->end; iwvl++){
	if(psfsize[iwvl]==1){
	    psf2s->p[iwvl]=strehlcomp(iopdevl, amp, wvl[iwvl]);
	}else{
	    long nembed=nembeds[iwvl];
	    long *embed=embeds[iwvl];
	    cmat *psf2=cnew(nembed,nembed);
	    int use1d;
	    int use1d_enable=0;
	    if(psfsize[iwvl]<nembed && use1d_enable){/*Want smaller PSF. */
		use1d=1;
		cfft2partialplan(psf2, psfsize[iwvl], -1);
	    }else{
		use1d=0;
		cfft2plan(psf2, -1);
	    }
	 
	    /*
	      //The following makes sum(psf)=1; 
	      double psfnorm=1./(sqrt(aper->sumamp2)*aper->nembed);
	    */
	    
	    /*
	      since sum(aper->amp)=1, this makes psf normalized by diffraction
	      limited PSF. max(psf) is strehl.
	    */
	    double psfnorm=1; 
	    dcomplex i2pi=I*2*M_PI/wvl[iwvl];
	    for(int iloc=0; iloc<iopdevl->nx; iloc++){
		psf2->p[embed[iloc]]=amp[iloc]*cexp(i2pi*iopdevl->p[iloc]);
	    }
	    if(use1d==1){
		cfft2partial(psf2, psfsize[iwvl], -1);
	    }else{
		cfft2(psf2,-1);
	    }
	    if(psfsize[iwvl]==nembed){/*just reference */
		cfftshift(psf2);
		psf2s->p[iwvl]=cref(psf2);
	    }else{/*create a new array, smaller. */
		psf2s->p[iwvl]=cnew(psfsize[iwvl], psfsize[iwvl]);
		ccpcorner2center(psf2s->p[iwvl], psf2);
	    }
	    if(fabs(psfnorm-1)>1.e-15) {
		cscale(psf2s->p[iwvl], psfnorm);
	    }
	    cfree(psf2);
	}
    }
}
/**
   Simple embed and accumulation. 
 */
void embed_in(double *out, const double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	out[embed[i]]+=in[i];
    }
}
/**
   Simple embed and accumulation
 */
void embed_out(const double *out, double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	in[i]+=out[embed[i]];
    }
}
/**
   Simple embed and accumulation.  for complex output
 */
void embedc_in(dcomplex *out, const double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	out[embed[i]]+=in[i];
    }
}
/**
   Simple embed and accumulation. for real output.
 */
void embedc_out(const dcomplex *out, double *in, long nin, long *embed){
    for(long i=0; i<nin; i++){
	in[i]+=creal(out[embed[i]]);
    }
}

/**
   Estimate anisoplanatic angle theta0 from Fried parameter r0, layer height and
   weights.  */
double calc_aniso(double r0, int nht, double *ht, double *wt){
    double wh=0;
    for(int iht=0; iht<nht; iht++){
	wh+=pow(ht[iht],5./3.)*wt[iht];
    }
    return 0.3144*r0*pow(wh,-3./5.);
}
/**
   Estimate generalized aniso angle theta2 from Fried parameter r0, and layer
   height and weights, and deformable mirror conjugation heights hc1 hc2 of the
   ground and altitude DMs. */
double calc_aniso2(double r0, int nht, double *ht, double *wt, double hc1, double hc2){
    double wh=0;
    double hh=pow(hc2-hc1,5./3.);
    for(int iht=0; iht<nht; iht++){
	double t1=0.5*pow(fabs(ht[iht]-hc1),5./3.)+0.5*pow(fabs(ht[iht]-hc2),5./3.);
	double t2=-0.25*hh-0.25/hh*pow(pow(fabs(ht[iht]-hc1),5./3.)-pow(fabs(ht[iht]-hc2),5./3.),2);
	wh+=wt[iht]*(t1+t2);
    }
    return 0.3144*r0*pow(wh,-3./5.);
}

/**
   prepare the integrator by shifting commands. similar to laos.
   inte->p[0]=inte->p[0]*ap[0]+inte->p[1]*ap[1]+...
*/
void shift_inte(int nap, double *ap, dcell **inte){
    dcell *tmp=NULL;
    dcell *keepjunk=inte[nap-1];
    for(int iap=nap-1; iap>=0; iap--){
	dcelladd(&tmp,1,inte[iap],ap[iap]);
	if(iap>0){
	    inte[iap]=inte[iap-1];/*shifting */
	}else{
	    inte[iap]=tmp;/*new command. */
	}
    }
    dcellfree(keepjunk);
}
char *evl_header(const PARMS_T *parms, const APER_T *aper, int ievl, int iwvl){
    char header[320];
    int nembed=aper->nembed[iwvl];
    double wvl=parms->evl.wvl[iwvl];
    double sumamp2=aper->sumamp2;
    int npos=parms->evl.psfmean;
    if(npos==1) {
	npos=parms->sim.end-parms->evl.psfisim;
    }
    snprintf(header, 320, 
	     "Science PSF at (%.15g, %.15g) arcsec\n"
	     "Turbulence: r0=%g, l0=%g\n"
	     "Wavelength: %.15gm\n"
	     "OPD Sampling: %.15gm\n"
	     "FFT Grid: %dx%d\n"
	     "PSF Sampling: %.15g arcsec\n"
	     "PSF Sum to %.15g\n"
	     "Exposure: %gs", 
	     ievl<0?0:parms->evl.thetax[ievl]*206265, ievl<0?0:parms->evl.thetay[ievl]*206265,
	     parms->atm.r0, parms->atm.l0,
	     wvl, parms->aper.dx, nembed, nembed, wvl/(nembed*parms->aper.dx)*206265,
	     sumamp2*nembed*nembed, parms->sim.dt*npos);
    return strdup(header);
}
