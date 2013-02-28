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
#include <getopt.h>
#include "skyc.h"
#include "parms.h"
#include "utils.h"
char *dirsetup;
/**
   Print out usage information.
 */
static void print_usage(void){
    info2(
	  "Usage: skyc [OPTION...] [FILE]...\n"
	  "skyc is a simulation tool developed to do sky coveragepostprocessing\n\n"
	  "Options: \n"
	  "-h, --help        to print out this message\n"
	  "-d, --detach      to detach from terminal and run in background\n"
	  "-f, --force       force starting simulation without scheduler\n"
	  "-n N, --nthread=N Use N threads, default is 1\n"
	  "-o DIR, --output=DIR\n"
	  "                  output the results to DIR.\n"
	  "-c FILE.conf, --conf=FILE.conf\n"
	  "                  Use FILE.conf as the baseline config instead of maos.conf\n"
	  );
    exit(0);
}
/**
   Parse command line arguments argc, argv
 */
ARG_S *parse_args(int argc, const char *argv[]){
    ARG_S *arg=calloc(1, sizeof(ARG_S));
    char *host=NULL; int local=0;
    ARGOPT_T options[]={
	{"help", 'h', T_INT, 2, print_usage, NULL},
	{"detach", 'd',T_INT, 0, &arg->detach, NULL},
	{"override",'O',T_INT,0, &arg->override, NULL},
	{"force",  'f',T_INT, 0, &arg->force, NULL},
	{"output", 'o',T_STR, 1, &arg->dirout, NULL},
	{"nthread",'n',T_INT, 1, &arg->nthread,NULL},
	{"conf",   'c',T_STR, 1, &arg->conf, NULL},
	{"path",   'P',T_STR, 3, addpath, NULL},
	{"run",    'r',T_STR, 1, &host, NULL},
	{"local",  'l',T_INT, 0, &local, NULL},
	{NULL, 0,0,0, NULL, NULL}
    };
    char *cmds=parse_argopt(argc, argv, options);
    if(!local && host){ //run through scheduler
	if(scheduler_launch_exe(host, argc, argv)){
	    error2("Unable to launch skyc at server %s\n", host);
	}
	exit(EXIT_SUCCESS);
    }
    /*do not use NTHREAD here. causes too much cpu load*/
    if(arg->nthread>NCPU || arg->nthread<=0){
	arg->nthread=NCPU;
    }else{
        NTHREAD=arg->nthread;
    }
    char fntmp[PATH_MAX];
    snprintf(fntmp,PATH_MAX,"%s/skyc_%ld.conf",TEMP,(long)getpid());
    FILE *fptmp=fopen(fntmp,"w");
    fputs(cmds, fptmp);
    fclose(fptmp);
    free(cmds); cmds=NULL;
    arg->confcmd=strdup(fntmp);
    if(!arg->detach){/*foreground task will start immediately. */
	arg->force=1;
    }
    if(!arg->dirout){
	arg->dirout=strtime();
    }
    if(!arg->conf){ /*If -c is not specifid in path, will use maos.conf*/
	arg->conf=strdup("maos.conf");
    }
    /*Setup PATH and result directory */
    char *config_path=find_config("skyc");
    addpath(config_path);
    free(config_path);

    mymkdir("%s",arg->dirout);
    addpath(".");
    if(chdir(arg->dirout)){
	error("Unable to chdir to %s\n", arg->dirout);
    }
    return arg;
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
    }
    snprintf(fnnew,256,"kill_%d",pid);
    if(exist(fnnew)) remove(fnnew);
    
    snprintf(fnold,256,"run_%d.log",pid);
    snprintf(fnnew,256,"run_%d.%s",pid,suffix);
    rename(fnold,fnnew);
    mysymlink(fnnew, "run_recent.log");
}
/**
   Handles signals. 
 */
void skyc_signal_handler(int sig){
    disable_signal_handler;
    rename_file(sig);/*handles signal */
    if(sig!=0){
	info2("Caught signal %d\n",sig);
	if(sig == SIGSEGV){
	    print_backtrace();
	}
	scheduler_finish(1);
	raise(sig);
	exit(sig);
    }
}

/**
   Add two PSDs that doesn't have the same frequency. the first column of each
   dmat is the frequency nu, and the second column is PSD*/
static dmat *add_psd_nomatch(dmat *psd1,dmat *psd2){
    dmat *nu1=dsub(psd1,0,psd1->nx,0,1);
    dmat *psd2x=dnew_ref(psd2->nx, 1, psd2->p);
    dmat *psd2y=dnew_ref(psd2->nx,1,psd2->p+psd2->nx);
    dmat *psd2new=dinterp1(psd2x,psd2y,nu1);
    dfree(psd2x); dfree(psd2y);
    dmat *psd=dnew(nu1->nx,2);
    double *ppsd=psd->p+psd->nx;
    for(long i=0; i<psd->nx; i++){
	psd->p[i]=nu1->p[i];
	ppsd[i]=psd2new->p[i];
    }
    dfree(nu1);
    dfree(psd2new);
    return psd;
}
/**
   Add two PSDs. the first column of each dmat is the frequency nu, and the
   second column is PSD*/
dmat *add_psd(dmat *psd1, dmat *psd2){
    if(psd1->nx!=psd2->nx){
	//warning("The two PSDs have different length\n");
	return add_psd_nomatch(psd1, psd2);
    }
    dmat *psd=dnew(psd1->nx,2);
    double *restrict pp=psd->p+psd->nx;
    const double *p1=psd1->p+psd1->nx;
    const double *p2=psd2->p+psd2->nx;
    for(long i=0; i<psd->nx; i++){
	if(fabs(psd1->p[i]-psd2->p[i])>1.e-2){
	    warning("The two PSDs doesn't have the same freq.");
	    dfree(psd);
	    return add_psd_nomatch(psd1,psd2);
	    /*todo: apply interp1 to interpolate the second PSD. */
	}
	psd->p[i]=psd1->p[i];
	pp[i]=p1[i]+p2[i];
    }
    return psd;
}
/*
  Add a PSD to another.
*/
void add_psd2(dmat **out, dmat *in){
    if(!*out){
	*out=ddup(in);
    }else{
	dmat *tmp=*out;
	*out=add_psd(tmp, in);
	dfree(tmp);
    }
}
