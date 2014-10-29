/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <signal.h>
#include <fcntl.h>           /* For O_* constants */
#include <errno.h>
#include <getopt.h>
#include "common.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif

/**
   \file maos/utils.c
   A few utility routines
 */

/**
   Plot the loc, together with all beams
 */
void plotloc(char *fig, const PARMS_T *parms, 
	     loc_t *loc, double ht, char *format,...){
    format2fn;
    int ncir=parms->evl.nevl + parms->fit.nfit + parms->nwfs;
    if(parms->sim.ncpa_calib){
	ncir+=parms->sim.ncpa_ndir;
    }
    double (*cir)[4];
    cir=(double(*)[4])calloc(ncir*4,sizeof(double));
    int count=0;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double hs=parms->evl.hs->p[ievl];
	cir[count][0]=ht*parms->evl.thetax->p[ievl];
	cir[count][1]=ht*parms->evl.thetay->p[ievl];
	cir[count][2]=parms->aper.d*0.5*(1-ht/hs);
	cir[count][3]=0xFF0000;/*rgb color */
	count++;
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	double hs=parms->fit.hs->p[ifit];
	cir[count][0]=ht*parms->fit.thetax->p[ifit];
	cir[count][1]=ht*parms->fit.thetay->p[ifit];
	cir[count][2]=parms->aper.d*0.5*(1-ht/hs);
	cir[count][3]=0xFF22DD;/*rgb color */
	count++;
    }
    for(int idir=0; idir<parms->sim.ncpa_ndir; idir++){
	double hs=parms->sim.ncpa_hs->p[idir];
	cir[count][0]=ht*parms->sim.ncpa_thetax->p[idir];
	cir[count][1]=ht*parms->sim.ncpa_thetay->p[idir];
	cir[count][2]=parms->aper.d*0.5*(1-ht/hs);
	cir[count][3]=0x22FF00;/*rgb color */
	count++;
    }

    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	double hs=parms->wfs[iwfs].hs;
	int ipowfs=parms->wfs[iwfs].powfs;
	cir[count][0]=parms->wfs[iwfs].thetax*ht;
	cir[count][1]=parms->wfs[iwfs].thetay*ht;
	cir[count][2]=parms->aper.d*0.5*(1.-ht/hs);
	if(isfinite(hs)){//LGS
	    cir[count][3]=0xFF8800;
	}else if(!parms->powfs[ipowfs].lo){//Hi NGS
	    cir[count][3]=0xFFFF00;
	}else if(parms->powfs[ipowfs].order>1){//TTF
	    cir[count][3]=0x0000FF;//TTF
	}else{
	    cir[count][3]=0x0000FF;//TT
	}
	count++;
    }
    plot_points(fig, 1, &loc, NULL ,NULL, NULL,NULL,ncir, cir, NULL,
	       "Coordinate","x (m)","y (m)", "%s",fn);
    free(cir);
}
/**
   ploted all the different beam directions as points. */
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
    ngroup+=1;
    loccell *locs=cellnew(ngroup, 1);
    int32_t *style=calloc(ngroup, sizeof(int32_t));
    int count=0;
    style[count]=(0xFF0000<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->evl.nevl, 0, 0);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	locs->p[count]->locx[ievl]=parms->evl.thetax->p[ievl]*206265;
	locs->p[count]->locy[ievl]=parms->evl.thetay->p[ievl]*206265;
    }
    count++;

    style[count]=(0xFF22DD<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->fit.nfit, 0, 0);
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	locs->p[count]->locx[ifit]=parms->fit.thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->fit.thetay->p[ifit]*206265;
    }
    count++;
    style[count]=(0x22FF00<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->sim.ncpa_ndir, 0, 0);
    for(int ifit=0; ifit<parms->sim.ncpa_ndir; ifit++){
	locs->p[count]->locx[ifit]=parms->sim.ncpa_thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->sim.ncpa_thetay->p[ifit]*206265;
    }
    count++;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	locs->p[count]=locnew(parms->powfs[ipowfs].nwfs, 0, 0);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    locs->p[count]->locx[jwfs]=parms->wfs[iwfs].thetax*206265;
	    locs->p[count]->locy[jwfs]=parms->wfs[iwfs].thetay*206265;
	}
	if(isfinite(parms->powfs[ipowfs].hs)){
	    style[count]=(0xFF8800<<8)+(4<<4)+2;
	}else if(!parms->powfs[ipowfs].lo){
	    style[count]=(0xFFFF00<<8)+(4<<4)+1;
	}else if(parms->powfs[ipowfs].order>1){
	    style[count]=(0x0000FF<<8)+(4<<4)+4;
	}else{
	    style[count]=(0x0000FF<<8)+(4<<4)+1;
	}
	count++;
    }
    if(count!=ngroup){
	error("count=%d, ngroup=%d. they should equal.\n", count, ngroup);
    }
    double limit[4];
    limit[0]=limit[2]=-totfov/2;
    limit[1]=limit[3]=totfov/2;
    plot_points(fig, ngroup, locs->p, NULL, style,limit,NULL,ncir,cir, NULL,
		"Asterism","x (arcsec)", "y (arcsec)", "%s",fn);
    free(cir);
    cellfree(locs);
    free(style);
}
/**
   Rename the log files when simulation exits.
 */
void rename_file(int sig){
    draw_final(0);
    if(disable_save) return;
    if(sig==0){
	rename("run_recent.log", "run_done.log");
	mysymlink("run_done.log", "run_recent.log");
	rename("maos_recent.conf", "maos_done.conf");
	mysymlink("maos_done.conf", "maos_recent.conf");
    }
    if(global && global->parms && global->parms->fdlock && sig!=0){
	char fn[80];
	const PARMS_T *parms=global->parms;
	for(int iseed=global->iseed; iseed<parms->sim.nseed; iseed++){
	    if(parms->fdlock->p[iseed]>0){
		close(parms->fdlock->p[iseed]);
		int seed=parms->sim.seeds->p[iseed];
		snprintf(fn, 80, "Res_%d.lock",seed);
		if(exist(fn)){
		    (void) remove(fn);
		}
	    }
	}
    }
}
/**
   Handles signals.
 */
int maos_signal_handler(int sig){
    psignal(sig, "maos");
    rename_file(sig);/*handles signal */
    if(global->parms->sim.mvmport){
	mvm_client_close();
    }
    scheduler_finish(sig);
    return 0;
}
/**
   Print out usage information.
 */
static void print_usage(void){
    info2(
"Usage: maos [OPTION...] [FILE]...\n"
"maos is a simulation tool developed to adaptive optics systems\n\n"
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
"-o, --output=DIR  output the results to DIR.\n"
"-c, --conf=FILE.conf\n"
"                  Use FILE.conf as the baseline config instead of nfiraos.conf\n"
"-p, --path=dir    Add dir to the internal PATH\n"
"-g, --gpu=i       Use the i'th gpu. 0 for the first. -1 to disable. default: automatic\n"
"-G, --ngpu=N'     Use a total of N gpus.\n"
"-r, --run=host    Run the job in another host.\n"
	  );
    exit(0);
}

/**
   Parse command line arguments argc, argv
 */
ARG_T * parse_args(int argc, const char *argv[]){
    ARG_T *arg=calloc(1, sizeof(ARG_T));
    char *host=NULL;
    int nthread=0;
    ARGOPT_T options[]={
	{"help",   'h',T_INT, 2, print_usage, NULL},
	{"detach", 'd',T_INT, 0, &arg->detach, NULL},
	{"force",  'f',T_INT, 0, &arg->force, NULL},
	{"override",'O',T_INT,0, &arg->override, NULL},
	{"output", 'o',T_STR, 1, &arg->dirout, NULL},
	{"nthread",'n',T_INT, 1, &nthread,NULL},
	{"gpu",    'g',T_INTARR, 1, &arg->gpus, &arg->ngpu},
	{"ngpu",   'G',T_INT, 1, &arg->ngpu2, NULL},
	{"conf",   'c',T_STR, 1, &arg->conf, NULL},
	{"path",   'p',T_STR, 3, addpath, NULL},
	{"run",    'r',T_STR, 1, &host, NULL},
	{"server", 'S',T_INT, 0, &arg->server,NULL},
	{NULL, 0,0,0, NULL, NULL}
    };
    char *cmds=strnadd(argc-1, argv+1, " ");
    parse_argopt(cmds, options);
    if((host || arg->detach) && !arg->dirout){
	error("detach mode requires specifying -o\n");
    }
    if(!host && !arg->detach){//forground running.
	arg->force=1;
    }else if(getenv("MAOS_DIRECT_LAUNCH")){
	/*lanched through scheduler to run locally. We are already detached, so
	  don't daemonize again.*/
	arg->detach=0;
	arg->force=0;
	detached=1;
    }else{
#ifndef MAOS_DISABLE_SCHEDULER
	/*Detached version. Always launch through scheduler if available.*/
	int locally=0;
	if(!host){
	    locally=1;
	    host=strdup("localhost");
	}
	if(scheduler_launch_exe(host, argc, argv)){
	    warning("Launch locally\n");
	    if(!locally){
		error2("Unable to launch maos at server %s\n", host);
	    }
	}else{
	    exit(EXIT_SUCCESS);
	}
#endif
    }
    free(host);
    if(nthread<NTHREAD && nthread>0){
        NTHREAD=nthread;
    }
    if(arg->ngpu2>0){
	if(!arg->gpus || arg->ngpu==0){
	    arg->ngpu=arg->ngpu2;
	}else{
	    error("-g and -G cannot both be specified\n");
	}
    }
    arg->confcmd=cmds;
    
    addpath(".");
    if(arg->dirout){
	mymkdir("%s",arg->dirout);
	if(chdir(arg->dirout)){
	    error("Unable to chdir to %s\n", arg->dirout);
	}
    }else{
	warning2("Disable saving when no -o is supplied.\n");
	disable_save=1;
    }
    return arg;
}

char *evl_header(const PARMS_T *parms, const APER_T *aper, int ievl, int iwvl){
    char header[320];
    int nembed=aper->embed->nembed->p[iwvl];
    double wvl=parms->evl.wvl->p[iwvl];
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
	     "PSF Sum to: %.15g\n"
	     "Exposure: %gs\n", 
	     ievl<0?0:parms->evl.thetax->p[ievl]*206265, ievl<0?0:parms->evl.thetay->p[ievl]*206265,
	     parms->atm.r0, parms->atm.l0,
	     wvl, parms->evl.dx, nembed, nembed, wvl/(nembed*parms->evl.dx)*206265,
	     sumamp2*nembed*nembed, parms->sim.dt*npos);
    return strdup(header);
}
void apply_fieldstop(dmat *opd, dmat *amp, lmat *embed, long nembed, dmat *fieldstop, double wvl){
    cmat *wvf=cnew(nembed, nembed);
    //cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
    double kk=2*M_PI/wvl;
    double kki=1./kk;
    double wvlh=wvl*0.5;
    dcomplex i2pi=I*kk;
    for(int iloc=0; iloc<opd->nx; iloc++){
	wvf->p[embed->p[iloc]]=amp->p[iloc]*cexp(i2pi*opd->p[iloc]);
    }
    cfft2(wvf, -1);
    ccwmd(wvf, fieldstop, 1);
    cfft2(wvf, 1);
    for(int iloc=0; iloc<opd->nx; iloc++){
	double val=carg(wvf->p[embed->p[iloc]])*kki;
	double diff=fmod(val-opd->p[iloc]+wvlh, wvl);
	if(diff<0) diff+=wvl;
	opd->p[iloc]+=diff-wvlh;
    }
    cfree(wvf);
}
void plot_setup(const PARMS_T *parms, const POWFS_T *powfs,
		const APER_T *aper, const RECON_T *recon){
    if(!parms->plot.setup) return;
    plotdir("FoV",parms,parms->sim.fov*206265,"fov");/*plot wfs/evaluation direction */
    plotloc("FoV",parms,recon->ploc,0, "ploc");
    plotloc("FoV",parms,recon->floc,0, "floc");
    for(int idm=0; idm<parms->ndm; idm++){
	double ht=parms->dm[idm].ht;
	plotloc("FoV", parms, recon->aloc->p[idm], ht, "aloc%d", idm);
    }
    for(int ips=0; ips<recon->npsr; ips++){
	const double ht=recon->ht->p[ips];
	plotloc("FoV",parms,recon->xloc->p[ips],ht, "xloc%d",ips);
    }
    drawopd("amp",aper->locs,aper->amp1->p,NULL,"Aperture Amplitude Map",
	    "x (m)","y (m)","aper");
    
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	drawopd("amp", powfs[ipowfs].loc, powfs[ipowfs].amp->p,NULL,
		"WFS Amplitude Map","x (m)","y (m)","powfs %d", ipowfs);
	if(powfs[ipowfs].amp_tel){
	    for(int wfsind=0; wfsind<powfs[ipowfs].nwfs; wfsind++){
		drawopd("amp", powfs[ipowfs].loc, powfs[ipowfs].amp_tel->p[wfsind]->p,NULL,
			"WFS Amplitude Map","x (m)","y (m)","powfs %d tel2wfs", ipowfs);
	    }
	}
    }
}


/**
   Create WFS amplitude map from coordinate, masked with annular defined by (D,Din). 
*/
dmat *mkamp(loc_t *loc, map_t *ampground, double misregx, double misregy, double D, double Din){
    dmat *amp=dnew(loc->nloc, 1);
    if(ampground){
	prop_grid(ampground, loc, NULL, amp->p, 1, misregx, misregy,1,0,0,0);
    }else{
	locannular(amp->p, loc, -misregx, -misregy, D*0.5, Din*0.5, 1);
    }
    return amp;
}
