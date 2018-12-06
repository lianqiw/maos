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


#include <search.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include <fcntl.h>           /* For O_* constants */
#include <errno.h>
#include <getopt.h>
#include "common.h"
#include "mvm_client.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#include "setup_powfs.h"
#include "genseotf.h"
/**
   \file maos/utils.c
   A few utility routines
*/

/**
   Plot the loc, together with all beams
*/
void plotloc(const char *fig, const PARMS_T *parms, 
	     loc_t *loc, double ht, const char *format,...){
    format2fn;
    int ncir=parms->evl.nevl + parms->fit.nfit + parms->nwfs;
    if(parms->sim.ncpa_calib){
	ncir+=parms->sim.ncpa_ndir;
    }
    dmat *cir=dnew(4, ncir);
    int count=0;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double hs=parms->evl.hs->p[ievl];
	IND(cir,0,count)=ht*parms->evl.thetax->p[ievl];
	IND(cir,1,count)=ht*parms->evl.thetay->p[ievl];
	IND(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	IND(cir,3,count)=0xFF0000;/*rgb color */
	count++;
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	double hs=parms->fit.hs->p[ifit];
	IND(cir,0,count)=ht*parms->fit.thetax->p[ifit];
	IND(cir,1,count)=ht*parms->fit.thetay->p[ifit];
	IND(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	IND(cir,3,count)=0xFF22DD;/*rgb color */
	count++;
    }
    for(int idir=0; idir<parms->sim.ncpa_ndir; idir++){
	double hs=parms->sim.ncpa_hs->p[idir];
	IND(cir,0,count)=ht*parms->sim.ncpa_thetax->p[idir];
	IND(cir,1,count)=ht*parms->sim.ncpa_thetay->p[idir];
	IND(cir,2,count)=parms->aper.d*0.5*(1-ht/hs);
	IND(cir,3,count)=0x22FF00;/*rgb color */
	count++;
    }

    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	double hs=parms->wfs[iwfs].hs;
	int ipowfs=parms->wfs[iwfs].powfs;
	IND(cir,0,count)=parms->wfs[iwfs].thetax*ht;
	IND(cir,1,count)=parms->wfs[iwfs].thetay*ht;
	IND(cir,2,count)=parms->aper.d*0.5*(1.-ht/hs);
	if(isfinite(hs)){//LGS
	    IND(cir,3,count)=0xFF8800;
	}else if(!parms->powfs[ipowfs].lo){//Hi NGS
	    IND(cir,3,count)=0xFFFF00;
	}else if(parms->powfs[ipowfs].order>1){//TTF
	    IND(cir,3,count)=0x0000FF;//TTF
	}else{
	    IND(cir,3,count)=0x0000FF;//TT
	}
	count++;
    }
    plot_points(fig, 1, &loc, NULL ,NULL, NULL,NULL, cir, NULL,
		"Coordinate","x (m)","y (m)", "%s",fn);
    dfree(cir);
}
/**
   ploted all the different beam directions as points. */
void plotdir(const char *fig, const PARMS_T *parms, double totfov, const char *format,...){
    format2fn;
    int ncir=1;
    dmat *cir=dnew(4, ncir);
    IND(cir,0,0)=0;
    IND(cir,1,0)=0;
    IND(cir,2,0)=totfov/2;
    IND(cir,3,0)=0x000000;/*rgb color */
    int ngroup=2+parms->npowfs;
    ngroup+=1;
    const char *legend[ngroup];
    loccell *locs=(loccell*)cellnew(ngroup, 1);
    int32_t *style=mycalloc(ngroup,int32_t);
    int count=0;
    legend[count]="Evaluation";
    style[count]=(0xFF0000<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->evl.nevl, 0, 0);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	locs->p[count]->locx[ievl]=parms->evl.thetax->p[ievl]*206265;
	locs->p[count]->locy[ievl]=parms->evl.thetay->p[ievl]*206265;
    }
    count++;
    legend[count]="DM Fitting";
    style[count]=(0xFF22DD<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->fit.nfit, 0, 0);
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	locs->p[count]->locx[ifit]=parms->fit.thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->fit.thetay->p[ifit]*206265;
    }
    count++;
    legend[count]="NCPA";
    style[count]=(0x22FF00<<8)+(4<<4)+3;
    locs->p[count]=locnew(parms->sim.ncpa_ndir, 0, 0);
    for(int ifit=0; ifit<parms->sim.ncpa_ndir; ifit++){
	locs->p[count]->locx[ifit]=parms->sim.ncpa_thetax->p[ifit]*206265;
	locs->p[count]->locy[ifit]=parms->sim.ncpa_thetay->p[ifit]*206265;
    }
    count++;
    char* const legwfs[]={
	"LGS WFS",
	"NGS WFS",
	"PWFS",
	"TTF WFS",
	"TT WFS",
	"Other WFS",
    };
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	int ilegwfs=6;
	if(parms->powfs[ipowfs].lo){
	    if(parms->powfs[ipowfs].order==1){
		ilegwfs=4;
	    }else{
		ilegwfs=3;
	    }
	}else{
	    if(parms->powfs[ipowfs].trs){
		ilegwfs=0;
	    }else{
		if(parms->powfs[ipowfs].type==1){
		    ilegwfs=2;
		}else{
		    ilegwfs=1;
		}
	    }
	}
	legend[count]=legwfs[ilegwfs];
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
    plot_points(fig, ngroup, locs->p, NULL, style,limit,NULL, cir, legend,
		"Asterism","x (arcsec)", "y (arcsec)", "%s",fn);
    dfree(cir);
    cellfree(locs);
    free(style);
}
/**
   Rename the log files when simulation exits.
*/
void rename_file(int sig){
    draw_final(1);
    if(disable_save) return;
    if(sig==0){
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
	remove("run_done.log");
	mysymlink(fn, "run_done.log");
	snprintf(fn, PATH_MAX, "maos_%s_%ld.conf", HOST, (long)getpid());
	remove("maos_done.conf");
	mysymlink(fn, "maos_done.conf");
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
    info("maos: %s", sys_siglist[sig]);
    rename_file(sig);/*handles signal */
    if(global && global->parms && global->parms->sim.mvmport){
	mvm_client_close();
    }
    scheduler_finish(sig);
    return 0;
}
/**
   Print out usage information.
*/
static void print_usage(void){
    info(
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
    ARG_T *arg=mycalloc(1,ARG_T);
    char *host=NULL;
    int nthread=0;
    ARGOPT_T options[]={
	{"help",   'h',M_INT, 0, 1, (void*)print_usage, NULL},
	{"detach", 'd',M_INT, 0, 0, &arg->detach, NULL},
	{"force",  'f',M_INT, 0, 0, &arg->force, NULL},
	{"override",'O',M_INT,0, 0, &arg->override, NULL},
	{"output", 'o',M_STR, 1, 0, &arg->dirout, NULL},
	{"nthread",'n',M_INT, 1, 0, &nthread,NULL},
	{"gpu",    'g',M_INT, 2, 0, &arg->gpus, &arg->ngpu},
	{"ngpu",   'G',M_INT, 1, 0, &arg->ngpu2, NULL},
	{"conf",   'c',M_STR, 1, 0, &arg->conf, NULL},
	{"path",   'p',M_STR, 1, 1, (void*)addpath, NULL},
	{"run",    'r',M_STR, 1, 0, &host, NULL},
	{"server", 'S',M_INT, 0, 0, &arg->server,NULL},
	{NULL,     0,  0,     0, 0, NULL, NULL}
    };
    char *cmds=strnadd(argc-1, argv+1, " ");
    parse_argopt(cmds, options);
    if((host || arg->detach) && !arg->dirout){
	error("detach mode requires specifying -o\n");
    }
    if(!host && !arg->detach){//forground running.
	arg->force=1;
    }else if(getenv("MAOS_DIRECT_LAUNCH")){
	/*being lanuched by scheduler. We are already detached, so don't daemonize again.*/
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
	    if(!locally){
		error("Unable to launch maos at server %s.\n", host);
	    }else{
		warning("Launch locally without scheduler.\n");
	    }
	}else{
	    exit(EXIT_SUCCESS);
	}
#else
	arg->force=1;//launch directly when scheduler is disabled.
#endif
    }
    free(host);
    if(nthread<MAXTHREAD && nthread>0){
        NTHREAD=nthread;
    }
    if(arg->ngpu2>0){
	if(!arg->gpus || arg->ngpu==0){
	    arg->ngpu=arg->ngpu2;
	}else{
	    error("-g and -G cannot be specified simultaneously\n");
	}
    }else if(arg->ngpu){//check for -g-1
	for(int ig=arg->ngpu-1; ig>=0; ig--){
	    if(arg->gpus[ig]<0){
		if(ig+1==arg->ngpu){//-g-1 appear last
		    arg->gpus[0]=-1;
		    arg->ngpu=1;
		}else{
		    //-g-1 is not last. It invalides previous -g's
		    arg->ngpu=arg->ngpu-(1+ig);
		    memcpy(arg->gpus, arg->gpus+ig+1, arg->ngpu*sizeof(int));
		}
		break;
	    }
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
	warning("Disable saving when no -o is supplied.\n");
	disable_save=1;
    }
    return arg;
}

char *evl_header(const PARMS_T *parms, const APER_T *aper, int ievl, int iwvl, int isim){
    char header[320];
    int nembed=aper->embed->nembed->p[iwvl];
    double wvl=parms->evl.wvl->p[iwvl];
    double sumamp2=aper->sumamp2;
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
	     parms->atm.r0, parms->atm.L0->p[0],
	     wvl, parms->evl.dx, nembed, nembed, wvl/(nembed*parms->evl.dx)*206265,
	     sumamp2*nembed*nembed, parms->sim.dt*(isim-parms->evl.psfisim+1));
    return strdup(header);
}
void apply_fieldstop(dmat *opd, dmat *amp, lmat *embed, long nembed, dmat *fieldstop, double wvl){
    cmat *wvf=cnew(nembed, nembed);
    //cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
    double kk=2*M_PI/wvl;
    double kki=1./kk;
    double wvlh=wvl*0.5;
    dcomplex i2pi=COMPLEX(0, kk);
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
    plotdir("Aperture",parms,parms->sim.fov*206265,"fov");/*plot wfs/evaluation direction */
    plotloc("Aperture",parms,recon->ploc,0, "ploc");
    plotloc("Aperture",parms,recon->floc,0, "floc");
    for(int idm=0; idm<parms->ndm; idm++){
	double ht=parms->dm[idm].ht;
	plotloc("Aperture", parms, recon->aloc->p[idm], ht, "aloc%d", idm);
    }
    for(int ips=0; ips<recon->npsr; ips++){
	const double ht=recon->ht->p[ips];
	plotloc("Aperture",parms,recon->xloc->p[ips],ht, "xloc%d",ips);
    }
    drawopd("Aperture",aper->locs,aper->amp1->p,NULL,"Aperture Amplitude Map",
	    "x (m)","y (m)","aper");
    
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	drawopd("Aperture", powfs[ipowfs].loc, powfs[ipowfs].amp->p,NULL,
		"WFS Amplitude Map","x (m)","y (m)","powfs %d", ipowfs);
	if(powfs[ipowfs].amp_tel){
	    for(int wfsind=0; wfsind<parms->powfs[ipowfs].nwfs; wfsind++){
		drawopd("Aperture", powfs[ipowfs].loc, powfs[ipowfs].amp_tel->p[wfsind]->p,NULL,
			"WFS Amplitude Map","x (m)","y (m)","powfs %d tel2wfs", ipowfs);
	    }
	}
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    if(powfs[ipowfs].gradncpa){
		drawopd("Goffx",powfs[ipowfs].saloc, powfs[ipowfs].gradncpa->p[jwfs]->p,NULL,
			"WFS Offset (x)","x (m)", "y (m)", "x %d",  iwfs);
		drawopd("Goffy",powfs[ipowfs].saloc, powfs[ipowfs].gradncpa->p[jwfs]->p+nsa, NULL,
			"WFS Offset (y)","x (m)", "y (m)", "y %d",  iwfs);
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
	prop_grid(ampground, loc, amp->p, 1, misregx, misregy,1,0,0,0);
    }else{
	locannular(amp->p, loc, -misregx, -misregy, D*0.5, Din*0.5, 1);
    }
    return amp;
}

void wfslinearity(const PARMS_T *parms, POWFS_T *powfs, const int iwfs){
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int nsa=powfs[ipowfs].saloc->nloc;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    ccell *potf=intstat->potf->p[intstat->nsepsf>1?wfsind:0];
    cmat *potf2=0;
    ccell *otf=ccellnew(nwvl,1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	otf->p[iwvl]=cnew(potf->p[0]->nx, potf->p[0]->ny);
    }
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    dmat **mtche=NULL;
    if(parms->powfs[ipowfs].phytype_sim==1){
	if(powfs[ipowfs].intstat->mtche->ny==1){
	    mtche=powfs[ipowfs].intstat->mtche->p;
	}else{
	    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
	}
    }
    int nllt=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->n:0;
    double *srot=NULL;
    cmat ***petf=NULL;
    void (*pccwm)(cmat*,const cmat*)=NULL;
    if(nllt){
	srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p;
	petf=mycalloc(nwvl,cmat**);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(powfs[ipowfs].etfsim[iwvl].p1){
		pccwm=ccwmcol;
		if(powfs[ipowfs].etfsim[iwvl].p1->ny==1)
		    petf[iwvl]=powfs[ipowfs].etfsim[iwvl].p1->p;
		else
		    petf[iwvl]=powfs[ipowfs].etfsim[iwvl].p1->p+wfsind*nsa;
	    }else{
		pccwm=ccwm;
		if(powfs[ipowfs].etfsim[iwvl].p2->ny==1)
		    petf[iwvl]=powfs[ipowfs].etfsim[iwvl].p2->p;
		else
		    petf[iwvl]=powfs[ipowfs].etfsim[iwvl].p2->p+wfsind*nsa;
	    }
	}
    }

    const int nsep=41;
    const double dg=0.1;
    double gx=0, gy=0, dgx=0, dgy=0;
    dmat *ints=dnew(powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
    double theta=0, cx=1, sx=0;
    double isep0=-(nsep-1)*0.5;
    dmat *xg=dlinspace(isep0*dg, dg, nsep);
    writebin(xg, "wfslinearity_wfs%d_sep", iwfs);
    dfree(xg);
    dmat *gnfra=0;
    if(srot){
	gnfra=dnew(nsep,nsa*2);
    }
    dmat *gnfxy=dnew(nsep,nsa*2);
    dmat* pgnfxy=gnfxy/*PDMAT*/;
    dmat* pgnfra=gnfra/*PDMAT*/;
    const int ndir=4;
    const char *dirs[]={"x", "y", "r", "a"};
    const char *types[]={"","MF", "CoG", "MAP"};
    if(parms->powfs[ipowfs].mtchcr){
	types[1]="MFC";
    }
    int type=parms->powfs[ipowfs].phytype_sim;
    for(int dir=0; dir<ndir; dir++){
	dzero(gnfra);
	dzero(gnfxy);
	for(int isa=0; isa<nsa; isa++){
	    info2("isa=%4d\b\b\b\b\b\b\b\b", nsa);
	    if(srot){
		theta=srot[isa];
		cx=cos(theta);
		sx=sin(theta);
	    }
	    switch(dir){
	    case 0://x
		dgx=dg*pixthetax;
		dgy=0;
		break;
	    case 1://y
		dgx=0;
		dgy=dg*pixthetay;
		break;
	    case 2://r
		dgx=cx*dg*pixthetax;
		dgy=sx*dg*pixthetax;
		break;
	    case 3://a
		dgx=-sx*dg*pixthetay;
		dgy=cx*dg*pixthetay;
		break;
	    default:
		error("Invalid dir\n");
	    }
	    for(int isep=0; isep<nsep; isep++){
		gx=dgx*(isep+isep0);
		gy=dgy*(isep+isep0);
		dzero(ints);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    double wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
			*parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
		    int idtf=powfs[ipowfs].dtf[iwvl].si->ny>1?wfsind:0;
		    int idtfsa=powfs[ipowfs].dtf[iwvl].si->nx>1?isa:0;
		    dspcell*  psi=powfs[ipowfs].dtf[iwvl].si/*PDSPCELL*/;
		    dsp *sis=IND(psi,idtfsa,idtf);
		    double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
		    double dtheta1=powfs[ipowfs].pts->nx*powfs[ipowfs].pts->dx*parms->powfs[ipowfs].embfac/wvl;
		    if(petf){
			ccp(&potf2, potf->p[isa+nsa*iwvl]);
			(*pccwm)(potf2,petf[iwvl][isa]);
		    }else{
			potf2=potf->p[isa+nsa*iwvl];
		    }
		    ctilt2(otf->p[iwvl], potf2, gx*dtheta1, gy*dtheta1, 0);
		    cfft2(otf->p[iwvl], 1);
		    dspmulcreal(ints->p, sis, otf->p[iwvl]->p, wvlsig);
		}
		//ddraw("ints", ints, NULL, NULL, "ints", "x", "y", "ints"); PAUSE;
		double g[3]={0,0,0};
		//notice that all the following gives gradients in x/y coord only.
		switch(type){
		case 0://no-op
		    break;
		case 1:{/*(constraint) Matched filter give gradients along x/y*/
		    dmulvec(g, mtche[isa], ints->p,1.);
		}
		    break;
		case 2:{/*tCoG gives gradients along r/a*/
		    dcog(g,ints,0.,0.,
			 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2],
			 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1], 0);
		    g[0]*=pixthetax;
		    g[1]*=pixthetay;
		}
		    break;
		case 3:{/*MAP gives gradient along r/a(?)*/
		    g[0]=gx+pixthetax*0.1;
		    g[1]=gy+pixthetay*0.1;
		    g[2]=1;
		    maxapriori(g, ints, parms, powfs, iwfs, isa, 1, 0, 1);
		}
		    break;
		default:
		    error("Invalid");
		}
		if(type==1 || !srot){
		    IND(pgnfxy,isep,isa)=g[0]/pixthetax;
		    IND(pgnfxy,isep,isa+nsa)=g[1]/pixthetay;
		    
		    if(srot){/*obtain gradients in r/a coord*/
			IND(pgnfra,isep,isa)=(g[0]*cx+g[1]*sx)/pixthetax;
			IND(pgnfra,isep,isa+nsa)=(-g[0]*sx+g[1]*cx)/pixthetay;
		    }
		}else{
		    IND(pgnfra,isep,isa)=g[0]/pixthetax;
		    IND(pgnfra,isep,isa+nsa)=g[1]/pixthetay;
		    IND(pgnfxy,isep,isa)=(g[0]*cx-g[1]*sx)/pixthetax;
		    IND(pgnfxy,isep,isa+nsa)=(g[0]*sx+g[1]*cx)/pixthetay;
		}
	    }
	}/*for isa*/
	writebin(gnfxy, "wfslinearity_xy_wfs%d_%s_%s", iwfs, types[type],dirs[dir]);
	if(srot){
	    writebin(gnfra, "wfslinearity_ra_wfs%d_%s_%s", iwfs, types[type],dirs[dir]);
	}
    }
    dfree(gnfxy);
    dfree(gnfra);
    dfree(ints);
    ccellfree(otf);
    if(petf){
	free(petf);
	cfree(potf2);
    }
}
/**
   Compute spherical aberration
*/
void lgs_wfs_sph_psd(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon, const int iwfs){
    int ipowfs=parms->wfs[iwfs].powfs;
    /*First save matched filter and i0*/
    dcell *mtche=powfs[ipowfs].intstat->mtche;
    dcell *i0=powfs[ipowfs].intstat->i0;
    dmat *i0sum=powfs[ipowfs].intstat->i0sum;
    powfs[ipowfs].intstat->mtche=0;
    powfs[ipowfs].intstat->i0=0;
    powfs[ipowfs].intstat->i0sum=0;
    //writebin(i0, "i0_initial");
    //writebin(mtche, "mtche_initial");
    dmat *opd=zernike(recon->ploc, parms->aper.d, 2, 15, 1);
    //writebin(opd, "ropd");
    dmat *GR=0;
    dspmm(&GR, recon->GP->p[iwfs], opd, "nn", 1);
    dfree(opd);
    dmat *RR=dpinv(GR, recon->saneai->p[iwfs+iwfs*parms->nwfsr]);
    //writebin(GR, "GR");
    //writebin(RR, "RR");
    int nsa=powfs[ipowfs].saloc->nloc;
    dmat *gradmf=dnew(nsa, 2);
    dmat *gradcg=dnew(nsa, 2);
    int ncol=1000;
    int dtrat=parms->powfs[ipowfs].llt->coldtrat;
    dmat *rmodmf=dnew(RR->nx, ncol/dtrat);
    dmat *rmodcg=dnew(RR->nx, ncol/dtrat);
    double scale=1;
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    double *srot=(parms->powfs[ipowfs].radpix)?
	powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p:NULL;
    for(int icol=0; icol<1000; icol+=dtrat){
	setup_powfs_etf(powfs, parms, ipowfs, 0, icol);
	gensei(parms, powfs, ipowfs);
	dcell *i0_new=powfs[ipowfs].intstat->i0;
	//writebin(i0_new, "i0_%d", icol);
	for(int isa=0; isa<nsa; isa++){
	    double geach[3]={0,0,1};
	    dmulvec(geach, mtche->p[isa], i0_new->p[isa]->p, 1);
	    if(parms->powfs[ipowfs].sigmatch){
		scale=i0sum->p[isa]/dsum(i0_new->p[isa]);
	    }
	    gradmf->p[isa]=geach[0]*scale;
	    gradmf->p[isa+nsa]=geach[1]*scale;
	    {
		dcog(geach, i0_new->p[isa], 0, 0, 
		     powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2],
		     powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1], 0);
		geach[0]*=pixthetax;
		geach[1]*=pixthetay;
		if(srot){
		    double theta=srot[isa];
		    double cx=cos(theta);
		    double sx=sin(theta);
		    double tmp=geach[0]*cx-geach[1]*sx;
		    geach[1]=geach[0]*sx+geach[1]*cx;
		    geach[0]=tmp;
		}
		gradcg->p[isa]=geach[0];
		gradcg->p[isa+nsa]=geach[1];
	    }
	    
	}
	dmulvec(rmodmf->p+icol/dtrat*rmodmf->nx, RR, gradmf->p, 1);
	dmulvec(rmodcg->p+icol/dtrat*rmodcg->nx, RR, gradcg->p, 1);
    }
    dfree(GR);
    dfree(RR);
    writebin(rmodmf, "rmod_mf");
    writebin(rmodcg, "rmod_cg");
    dfree(rmodmf);
    dfree(rmodcg);
    dfree(gradmf);
    dfree(gradcg);
    dcellfree(mtche);
    dcellfree(i0);
    dfree(i0sum);
}
typedef struct {
    const PARMS_T *parms;
    const POWFS_T *powfs;
    const dmat *ints;
    ccell *fotf;
    ccell *otf;//temporary.
    double bkgrnd;
    double rne;
    int noisy;
    int iwfs;
    int isa;
}mapdata_t;
/**
   The function to evaluate the result at x.
*/
static double mapfun(double *x, mapdata_t *info){
    const dmat *ints=info->ints;
    ccell *fotf=info->fotf;
    ccell *otf=info->otf;
    const PARMS_T *parms=info->parms;
    const POWFS_T *powfs=info->powfs;
    int iwfs=info->iwfs;
    int ipowfs=parms->wfs[iwfs].powfs;
    int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    int isa=info->isa;
    int nsa=fotf->nx;
    int nwvl=fotf->ny;
    if(!otf){
	info->otf=ccellnew(nwvl,1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    info->otf->p[iwvl]=cnew(info->fotf->p[0]->nx, info->fotf->p[0]->ny);
	    //cfft2plan(info->otf->p[iwvl], 1);
	    //cfft2plan(info->otf->p[iwvl], -1);
	}
	otf=info->otf;
    }
    dmat *ints2=dnew(ints->nx, ints->ny);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
	    *parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	dspcell*  psi=powfs[ipowfs].dtf[iwvl].si/*PDSPCELL*/;
	int idtf=powfs[ipowfs].dtf[iwvl].si->ny>1?wfsind:0;
	int idtfsa=powfs[ipowfs].dtf[iwvl].si->nx>1?isa:0;
	dsp *sis=IND(psi,idtfsa,idtf);
	double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	double dtheta1=powfs[ipowfs].pts->nx*powfs[ipowfs].pts->dx*parms->powfs[ipowfs].embfac/wvl;
	ctilt2(info->otf->p[iwvl], info->fotf->p[isa+nsa*iwvl], x[0]*dtheta1, x[1]*dtheta1, 0);
	cfft2(info->otf->p[iwvl], 1);
	dspmulcreal(ints2->p, sis, info->otf->p[iwvl]->p, wvlsig*x[2]);
    }
 
    double sigma=0;
    if(info->noisy){
	double noise=info->rne*info->rne+info->bkgrnd;
	for(int i=0; i<ints->nx*ints->ny; i++){
	    sigma+=pow(ints->p[i]-ints2->p[i],2)/(ints2->p[i]+noise);
	}
    }else{
	for(int i=0; i<ints->nx*ints->ny; i++){
	    sigma+=pow(ints->p[i]-ints2->p[i],2);
	}
    }
    /*dbg("Map fun called with [%g %g] %g, sigma=%g. noisy=%d\n", x[0], x[1], x[2], sigma, info->noisy);*/
    dfree(ints2);
    return sigma;
}
/**
   Implements MAP tracking algorithm. The polar coordinate is implicitly taken care of in mapfun if parms->powfs.radrot=0;
*/
void maxapriori(double *g, const dmat *ints, const PARMS_T *parms, 
		const POWFS_T *powfs, int iwfs, int isa, int noisy,
		double bkgrnd, double rne){
    int ipowfs=parms->wfs[iwfs].powfs;
    int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    ccell *fotf=intstat->fotf->p[intstat->nsepsf>1?wfsind:0];
    mapdata_t data={parms, powfs, ints, fotf, NULL, bkgrnd, rne, noisy, iwfs, isa};
    //info("isa %d: %.4e %.4e %.2f", isa, g[0], g[1], g[2]);
    int ncall=dminsearch(g, 3, MIN(pixthetax, pixthetay)*1e-2, 5000, (dminsearch_fun)mapfun, &data);
    ccellfree(data.otf);
    /* convert to native format along x/y or r/a to check for overflow*/
    if(parms->powfs[ipowfs].radpix && !parms->powfs[ipowfs].radrot){
	double theta=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p[isa];
	double cx=cos(theta);
	double sx=sin(theta);
	double tmp=g[0]*cx+g[1]*sx;
	g[1]=-g[0]*sx+g[1]*cx;
	g[0]=tmp;
    }
    double gx=g[0]/pixthetax*2./ints->nx;
    double gy=g[1]/pixthetay*2./ints->ny;
    if(fabs(gx)>0.55||fabs(gy)>0.55){
	warning("sa %4d iter %3d: wrapped: gx=%6.3f, gy=%6.3f ==> ", isa, ncall, gx, gy);
	gx=gx-floor(gx+0.5);
	gy=gy-floor(gy+0.5);
	warning("gx=%6.3f, gy=%6.3f\n", gx, gy);
	g[0]=pixthetax*ints->nx/2*gx;
	g[1]=pixthetay*ints->ny/2*gy;
    }
    //info("==> %.4e %.4e %.2f after %d iter\n", g[0], g[1], g[2], ncall);
    if(parms->powfs[ipowfs].radpix){
	double theta=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p[isa];
	double cx=cos(theta);
	double sx=sin(theta);
	double tmp=g[0]*cx-g[1]*sx;
	g[1]=g[0]*sx+g[1]*cx;
	g[0]=tmp;
    }
}

/**
   Compute the focus adjustment need to apply to OPD of wfs. Used in both CPU and GPU code.
*/
double wfsfocusadj(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    const POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const int isim=simu->isim;
    double focus=0;
    if(parms->powfs[ipowfs].llt){
	if(powfs[ipowfs].focus){
	    focus+=INDR(powfs[ipowfs].focus, isim, wfsind);
	}
	if(simu->zoomreal && parms->powfs[ipowfs].llt){
	    if(simu->zoompos && simu->zoompos->p[iwfs]){
		simu->zoompos->p[iwfs]->p[isim]=simu->zoomreal->p[iwfs];
	    }
	    focus-=simu->zoomreal->p[iwfs];
	}
    }
    if(simu->telfocusreal){
	focus-=simu->telfocusreal->p[0]->p[0];
    }
    return focus;
}
/**
   Expected averaged position of dithering signal during WFS integration. Called when (isim+1)%dtrat=0
*/
void dither_position(double *cs, double *ss, const PARMS_T *parms, int ipowfs, int isim, double deltam){
    //adjust for delay due to propagation, and computation delay. no effect when al+1=dtrat, which makes 2 wfs frame delay.
    const int adjust=parms->sim.alfsm+1-parms->powfs[ipowfs].dtrat;
    //adjust to get delay at beginning of integration
    const int adjust2=parms->powfs[ipowfs].dtrat-1;
    const double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
    const double angle=((isim-adjust-adjust2)/parms->powfs[ipowfs].dtrat)*anglei+deltam;
    const double angle2=((isim-adjust)/parms->powfs[ipowfs].dtrat)*anglei+deltam;
    const double delay=(double)adjust/parms->powfs[ipowfs].dtrat;
    const double beta=1+delay+floor(-delay);
    const double scale=1./(beta*beta+(1-beta)*(1-beta));
    //use average of two places during accumulation and scale
    *cs=(beta*cos(angle)+(1-beta)*cos(angle2))*scale;
    *ss=(beta*sin(angle)+(1-beta)*sin(angle2))*scale;
}
/**
   Find peak, then using parabolic fit on 3x3 window around it.
*/
/*
  void parabolic_peak_fit(double *grad, dmat *corr){
  double valmax=0;
  int jy=0, jx=0;
  //Find Peak location (jx, jy)
  for(int iy=1; iy<corr->ny-1; iy++){
  for(int ix=1; ix<corr->nx-1; ix++){
  if(IND(corr,ix,iy)>valmax){
  jy=iy; jx=ix;
  valmax=IND(corr,ix,iy);
  }
  }
  }
  //Calculate 1d sum of 3 row/columns.
  double vx[3], vy[3];
  for(long iy=0; iy<3; iy++){
  vy[iy]=0; vx[iy]=0;
  for(long ix=0; ix<3; ix++){
  vy[iy]+=IND(corr, ix+jx-1, iy+jy-1);
  vx[iy]+=IND(corr, iy+jx-1, ix+jy-1);
  }
  }
  //Parabolic fit.
  double px[2], py[2];
  px[0]=(vx[0]+vx[2])*0.5-vx[1];
  py[0]=(vy[0]+vy[2])*0.5-vy[1];
  px[1]=(vx[2]-vx[0])*0.5;
  py[1]=(vy[2]-vy[0])*0.5;
  //Center
  grad[0]=px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
  grad[1]=py[0]==0?0:(-py[1]/(2*py[0])+jy-(corr->ny-1)*0.5);
  }
*/
/**
   Fit 3 points around the peak of a 1-d profile.
*/
double parabolic_peak_1d(dmat *corr){
    double valmax=0;
    int jx=0;
    for(int ix=1; ix<corr->nx-1; ix++){
	if(IND(corr,ix)>valmax){
	    jx=ix;
	    valmax=IND(corr,ix);
	}
    }
    //Parabolic fit.
    double px[2];
    px[0]=(IND(corr, jx+1)+IND(corr, jx-1))*0.5-IND(corr, jx);
    px[1]=(IND(corr, jx+1)-IND(corr, jx-1))*0.5;
    return px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
}
/**
   First sum along 1 dimension, then fit 3 points around the peak. More robust than the old method.
*/
void parabolic_peak_sum(double *grad, dmat *corr, int nbox){
    const long nx=corr->nx;
    const long ny=corr->ny;
    if(nbox<=0 || nbox>nx) nbox=nx;
    if(nbox<=0 || nbox>ny) nbox=ny;
    dmat *corrx=dnew(nbox,1);
    dmat *corry=dnew(nbox,1);
    const long offx=(nx-nbox)/2;
    const long offy=(ny-nbox)/2;
    for(int iy=0; iy<nbox; iy++){
	for(int ix=0; ix<nbox; ix++){
	    IND(corrx, ix)+=IND(corr, ix+offx, iy+offy);
	    IND(corry, iy)+=IND(corr, ix+offx, iy+offy);
	}
    }
    grad[0]=parabolic_peak_1d(corrx);
    grad[1]=parabolic_peak_1d(corry);
    if (0){
	info("grad=%g, %g\n", grad[0], grad[1]);
	writebin(corr, "corr");
	exit(0);
    }
}
/**
   Calculate gradients using current specified algorithm
*/
void calc_phygrads(dmat **pgrad, dmat *ints[], const PARMS_T *parms, const POWFS_T *powfs, const int iwfs, const int phytype){
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    dmat **mtche=NULL;
    double *i0sum=NULL;
    double i0sumg=0;
    if(phytype==1){
	mtche=PINDR(powfs[ipowfs].intstat->mtche, 0, wfsind);
    }
    if(powfs[ipowfs].intstat->i0sum){
	i0sum=PINDR(powfs[ipowfs].intstat->i0sum, 0, wfsind);
	i0sumg=INDR(powfs[ipowfs].intstat->i0sumsum, wfsind, 0);
    }
    /*if(phytype==1){
	if(powfs[ipowfs].intstat->mtche->ny==1){
	    mtche=powfs[ipowfs].intstat->mtche->p;
	    i0sum=powfs[ipowfs].intstat->i0sum->p;
	    i0sumg=powfs[ipowfs].intstat->i0sumsum->p[0];
	}else{
	    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
	    i0sum=powfs[ipowfs].intstat->i0sum->p+nsa*wfsind;
	    i0sumg=powfs[ipowfs].intstat->i0sumsum->p[wfsind];
	    }
	    }*/
    const double *srot=(parms->powfs[ipowfs].radpix)?INDR(powfs[ipowfs].srot, wfsind, 0)->p:NULL;
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    /*output directly to simu->gradcl. replace */
    if(!*pgrad){
	*pgrad=dnew(nsa*2, 1);
    }
    double *pgradx=(*pgrad)->p;
    double *pgrady=pgradx+nsa;
    double i1sum=0;
    dmat *corr=0;
    if(parms->powfs[ipowfs].sigmatch==2){
	for(int isa=0; isa<nsa; isa++){
	    i1sum+=dsum(ints[isa]);//Disable the following to match GPU code.
	    /*const double thres=powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2];
	      const double offset=powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1];
	      for(int ip=0; ip<ints[isa]->nx*ints[isa]->ny; ip++){
	      double sig=ints[isa]->p[ip]-offset;
	      if(sig>thres){
	      i1sum+=sig;
	      }
	      }*/
	}
    }
    double sigtot=parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
    for(int isa=0; isa<nsa; isa++){
	double geach[3]={0,0,1};
	switch(phytype){
	case 1:{//matched filter
	    double scale=1.;
	    switch(parms->powfs[ipowfs].sigmatch){
	    case 0://no normalization
		break;
	    case 1://match per subaperture
		scale=i0sum[isa]/dsum(ints[isa]);
		break;
	    case 2://match globally.
		scale=i0sumg/i1sum;
		break;
	    }
	    dmulvec(geach, mtche[isa],ints[isa]->p,scale);
	}
	    break;
	case 2:{//CoG
	    double sumi=0;
	    switch(parms->powfs[ipowfs].sigmatch){
	    case 0://normalization use model intensity (linear model)
		if(i0sum){
		    sumi=i0sum[isa];
		}else{
		    sumi=sigtot*powfs[ipowfs].saa->p[isa];
		}
		break;
	    case 1://normalization use current intensity (non-linear)
		break;
	    case 2://normalized use scaled current intensity (non-linear)
		sumi=powfs[ipowfs].saa->p[isa]*(i1sum/powfs[ipowfs].saasum);
		warning_once("Please verify implementation is correct\n");
		break;
	    }
	    dcog(geach,ints[isa],0.,0.,
		 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2],
		 powfs[ipowfs].cogcoeff->p[wfsind]->p[isa*2+1], sumi);
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	  
	}
	    break;
	case 3:{//MAP: (to be removed. This algorithm is not very useful.)
	    geach[0]=pgradx[isa];//warm restart
	    geach[1]=pgrady[isa];
	    maxapriori(geach, ints[isa], parms, powfs, iwfs, isa, 1, bkgrnd, rne);
	}
	    break;
	case 4:{//Correlation. Peak first.
	    dcorr(&corr, ints[isa], powfs[ipowfs].intstat->i0->p[isa]);
	    dpara3(geach, corr);
	    if((fabs(geach[0])+fabs(geach[1]))>powfs[ipowfs].pixpsax){
		warning_once("wfs %d, subaperture %d has incorrect gradient measurement\n", iwfs, isa);
		geach[0]=0;
		geach[1]=0;
	    }
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	}
	    break;
	case 5:{//Correlation. Sum first (to be removed. Not working well)
	    dcorr(&corr, ints[isa], powfs[ipowfs].intstat->i0->p[isa]);
	    parabolic_peak_sum(geach, corr, 5);
	    geach[0]*=pixthetax;
	    geach[1]*=pixthetay;
	}
	    break;
	default:
	    error("Invalid");
	}
	if(phytype>1 && srot){
	    double theta=srot[isa];
	    double cx=cos(theta);
	    double sx=sin(theta);
	    double tmp=geach[0]*cx-geach[1]*sx;
	    geach[1]=geach[0]*sx+geach[1]*cx;
	    geach[0]=tmp;
	}
	pgradx[isa]=geach[0];
	pgrady[isa]=geach[1];
    }//for isa
    dfree(corr);
}
/**
   Read cell array from file specified by whole name or prefix.
*/
dcell *dcellread_prefix(const char *file, const PARMS_T *parms, int ipowfs){
    dcell *nea=0;
    if(!file){
	nea=0;
    }else if(zfexist(file)){
	//info("using %s\n", file);
	nea=dcellread(file);
    }else if(zfexist("%s_powfs%d.bin", file, ipowfs)){
	//info("using %s_powfs%d.bin\n", file, ipowfs);
	nea=dcellread("%s_powfs%d.bin", file, ipowfs);
    }else{
	nea=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<nea->nx; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    //info("using %s_wfs%d.bin\n", file, iwfs);
	    nea->p[jwfs]=dread("%s_wfs%d.bin", file, iwfs);
	}
    }
    return nea;
}
