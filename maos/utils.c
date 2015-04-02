/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    draw_final(1);
    if(disable_save) return;
    if(sig==0){
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "run_%s_%ld.log", myhostname(), (long)getpid());
	remove("run_done.log");
	rename(fn, "run_done.log");
	mysymlink("run_done.log", fn);
	remove("maos_done.conf");
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
		error("Unable to launch maos at server %s\n", host);
	    }
	}else{
	    exit(EXIT_SUCCESS);
	}
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
	     parms->atm.r0, parms->atm.L0,
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
	    for(int wfsind=0; wfsind<parms->powfs[ipowfs].nwfs; wfsind++){
		drawopd("amp", powfs[ipowfs].loc, powfs[ipowfs].amp_tel->p[wfsind]->p,NULL,
			"WFS Amplitude Map","x (m)","y (m)","powfs %d tel2wfs", ipowfs);
	    }
	}
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    if(powfs[ipowfs].gradoff){
		drawopd("Goffx",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p,NULL,
			"WFS Offset (x)","x (m)", "y (m)", "x %d",  iwfs);
		drawopd("Goffy",(loc_t*)powfs[ipowfs].pts, powfs[ipowfs].gradoff->p[jwfs]->p+nsa, NULL,
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
	prop_grid(ampground, loc, NULL, amp->p, 1, misregx, misregy,1,0,0,0);
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
    ccell *otf=cellnew(nwvl,1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	otf->p[iwvl]=cnew(potf->p[0]->nx, potf->p[0]->ny);
    }
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    dmat **mtche=NULL;
    if(parms->powfs[ipowfs].phytypesim==1){
	if(powfs[ipowfs].intstat->mtche->ny==1){
	    mtche=powfs[ipowfs].intstat->mtche->p;
	}else{
	    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
	}
    }
    int nllt=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->n:0;
    double *srot=NULL;
    const cmat ***petf=NULL;
    void (*pccwm)(cmat*,const cmat*)=NULL;
    if(nllt){
	srot=powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p;
	petf=calloc(nwvl, sizeof(void*));
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    if(powfs[ipowfs].etfsim[iwvl].p1){
		pccwm=ccwmcol;
		if(powfs[ipowfs].etfsim[iwvl].p1->ny==1)
		    petf[iwvl]=(void*)powfs[ipowfs].etfsim[iwvl].p1->p;
		else
		    petf[iwvl]=(void*)powfs[ipowfs].etfsim[iwvl].p1->p+wfsind*nsa;
	    }else{
		pccwm=ccwm;
		if(powfs[ipowfs].etfsim[iwvl].p2->ny==1)
		    petf[iwvl]=(void*)powfs[ipowfs].etfsim[iwvl].p2->p;
		else
		    petf[iwvl]=(void*)powfs[ipowfs].etfsim[iwvl].p2->p+wfsind*nsa;
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
    PDMAT(gnfxy,pgnfxy);
    PDMAT(gnfra,pgnfra);
    const int ndir=4;
    char *dirs[]={"x", "y", "r", "a"};
    char *types[]={"","MF", "CoG", "MAP"};
    if(parms->powfs[ipowfs].mtchcr){
	types[1]="MFC";
    }
    int type=parms->powfs[ipowfs].phytypesim;
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
		    PDSPCELL(powfs[ipowfs].dtf[iwvl].si, psi);
		    dsp *sis=psi[idtf][idtfsa];
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
		case 1:{/*(constraint) Matched filter give gradients along x/y*/
		    dmulvec(g, mtche[isa], ints->p,1.);
		}
		    break;
		case 2:{/*tCoG gives gradients along r/a*/
		    dcog(g,ints,0.,0.,
			 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2],
			 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2+1]);
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
		    pgnfxy[isa][isep]=g[0]/pixthetax;
		    pgnfxy[isa+nsa][isep]=g[1]/pixthetay;
		    
		    if(srot){/*obtain gradients in r/a coord*/
			pgnfra[isa][isep]=(g[0]*cx+g[1]*sx)/pixthetax;
			pgnfra[isa+nsa][isep]=(-g[0]*sx+g[1]*cx)/pixthetay;
		    }
		}else{
		    pgnfra[isa][isep]=g[0]/pixthetax;
		    pgnfra[isa+nsa][isep]=g[1]/pixthetay;
		    pgnfxy[isa][isep]=(g[0]*cx-g[1]*sx)/pixthetax;
		    pgnfxy[isa+nsa][isep]=(g[0]*sx+g[1]*cx)/pixthetay;
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
    dspmm(&GR, recon->GP->p[ipowfs], opd, 'n', 1);
    dfree(opd);
    dmat *RR=dpinv(GR, recon->saneai->p[iwfs+iwfs*parms->nwfsr]);
    //writebin(GR, "GR");
    //writebin(RR, "RR");
    int nsa=powfs[ipowfs].saloc->nloc;
    dmat *gradmf=dnew(nsa, 2);
    dmat *gradcg=dnew(nsa, 2);
    int ncol=1000;
    int dtrat=parms->powfs[ipowfs].llt->colsimdtrat;
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
	    if(parms->powfs[ipowfs].mtchscl){
		scale=i0sum->p[isa]/dsum(i0_new->p[isa]);
	    }
	    gradmf->p[isa]=geach[0]*scale;
	    gradmf->p[isa+nsa]=geach[1]*scale;
	    {
		dcog(geach, i0_new->p[isa], 0, 0, 
		     powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2],
		     powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2+1]);
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
	info->otf=cellnew(nwvl,1);
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
	PDSPCELL(powfs[ipowfs].dtf[iwvl].si, psi);
	int idtf=powfs[ipowfs].dtf[iwvl].si->ny>1?wfsind:0;
	int idtfsa=powfs[ipowfs].dtf[iwvl].si->nx>1?isa:0;
	dsp *sis=psi[idtf][idtfsa];
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
    /*info("Map fun called with [%g %g] %g, sigma=%g. noisy=%d\n", x[0], x[1], x[2], sigma, info->noisy);*/
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
    //info2("isa %d: %.4e %.4e %.2f", isa, g[0], g[1], g[2]);
    int ncall=dminsearch(g, 3, MIN(pixthetax, pixthetay)*1e-2, (dminsearch_fun)mapfun, &data);
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
	warning2("sa %4d iter %3d: wrapped: gx=%6.3f, gy=%6.3f ==> ", isa, ncall, gx, gy);
	gx=gx-floor(gx+0.5);
	gy=gy-floor(gy+0.5);
	warning2("gx=%6.3f, gy=%6.3f\n", gx, gy);
	g[0]=pixthetax*ints->nx/2*gx;
	g[1]=pixthetay*ints->ny/2*gy;
    }
    //info2("==> %.4e %.4e %.2f after %d iter\n", g[0], g[1], g[2], ncall);
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
    if(powfs[ipowfs].focus){
	const long nx=powfs[ipowfs].focus->nx;
	focus+=powfs[ipowfs].focus->p[(isim%nx)+nx*(powfs[ipowfs].focus->ny==parms->powfs[ipowfs].nwfs?wfsind:0)];
    }
    if(simu->zoomreal && parms->powfs[ipowfs].llt){
	if(simu->zoompos && simu->zoompos->p[iwfs]){
	    simu->zoompos->p[iwfs]->p[isim]=simu->zoomreal->p[iwfs];
	}
	focus-=simu->zoomreal->p[iwfs];
    }
    return focus;
}
/**
   Expected averaged position of dithering signal during WFS integration
*/
void dither_position(double *cs, double *ss, const PARMS_T *parms, int ipowfs, int isim, double deltam){
    //adjust for delay due to propagation
    const int adjust=(parms->powfs[ipowfs].llt?parms->sim.alfsm:0)+1-parms->powfs[ipowfs].dtrat;
    //adjust to get delay at beginning of integration
    const int adjust2=parms->powfs[ipowfs].llt?(parms->powfs[ipowfs].dtrat-1):0;
    const double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
    const double angle=anglei*((isim-adjust-adjust2)/parms->powfs[ipowfs].dtrat)+deltam;
    const double angle2=anglei*((isim-adjust)/parms->powfs[ipowfs].dtrat)+deltam;
    const double delay=(double)adjust/parms->powfs[ipowfs].dtrat;
    const double beta=1+delay+floor(-delay);
    const double scale=1./(beta*beta+(1-beta)*(1-beta));
    //use average of two places during accumulation and scale
    *cs=(beta*cos(angle)+(1-beta)*cos(angle2))*scale;
    *ss=(beta*sin(angle)+(1-beta)*sin(angle2))*scale;
}

/**
   Calculate gradients using current specified algorithm
*/
void calc_phygrads(dmat **pgrad, dmat *ints[], const PARMS_T *parms, const POWFS_T *powfs, int iwfs, int phytype){
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    dmat **mtche=NULL;
    double *i0sum=NULL;
    if(phytype==1){
	if(powfs[ipowfs].intstat->mtche->ny==1){
	    mtche=powfs[ipowfs].intstat->mtche->p;
	    i0sum=powfs[ipowfs].intstat->i0sum->p;
	}else{
	    mtche=powfs[ipowfs].intstat->mtche->p+nsa*wfsind;
	    i0sum=powfs[ipowfs].intstat->i0sum->p+nsa*wfsind;
	}
    }
    const double *srot=(parms->powfs[ipowfs].radpix)?powfs[ipowfs].srot->p[powfs[ipowfs].srot->ny>1?wfsind:0]->p:NULL;
    double pixthetax=parms->powfs[ipowfs].radpixtheta;
    double pixthetay=parms->powfs[ipowfs].pixtheta;
    /*output directly to simu->gradcl. replace */
    if(!*pgrad){
	*pgrad=dnew(nsa*2, 1);
    }
    double *pgradx=(*pgrad)->p;
    double *pgrady=pgradx+nsa;
    for(int isa=0; isa<nsa; isa++){
	double geach[3]={0,0,1};
	switch(phytype){
	case 1:{
	    dmulvec(geach, mtche[isa],ints[isa]->p,1.);
	    if(parms->powfs[ipowfs].mtchscl){
		double scale=i0sum[isa]/dsum(ints[isa]);
		geach[0]*=scale;
		geach[1]*=scale;
	    }
	}
	    break;
	case 2:{
	    dcog(geach,ints[isa],0.,0.,
		 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2],
		 powfs[ipowfs].intstat->cogcoeff->p[wfsind]->p[isa*2+1]);
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
	}
	    break;
	case 3:{
	    geach[0]=pgradx[isa];//warm restart
	    geach[1]=pgrady[isa];
	    maxapriori(geach, ints[isa], parms, powfs, iwfs, isa, 1, bkgrnd, rne);
	}
	    break;
	default:
	    error("Invalid");
	}
	
	pgradx[isa]=geach[0];
	pgrady[isa]=geach[1];
    }//for isa
}
