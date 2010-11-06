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
	  "                  Use FILE.conf as the baseline config instead of nfiraos.conf\n"
	  );
}
/**
   Parse command line arguments argc, argv
 */
ARG_S *parse_args(int argc, char **argv){
    ARG_S *arg=calloc(1, sizeof(ARG_S));
    static struct option long_options[]={
	{"help",0,0,'h'},
	{"detach",0,0,'d'},
	{"force",0,0,'f'},
	{"output",1,0,'o'},
	{"nthread",1,0,'n'},
	{"conf",1,0,'c'},
	{"path",1,0,'p'},
	{NULL,0,0,0}
    };
    arg->nthread=1;
    while(1){
	int option_index = 0;
	int c = getopt_long(argc, argv, "hdfo:n:c:p:",
			    long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'h':
	    print_usage();
	    exit(0);
	    break;
	case 'd':
	    arg->detach=1;
	    break;
	case 'f':
	    arg->force=1; 
	    break;
	case 'o':
	    if(arg->dirout){
		error("Duplicate argument for dirout\n");
	    }
	    arg->dirout=strdup(optarg);
	    break;
	case 'n':
	    arg->nthread=strtol(optarg,NULL,10);
	    if(arg->nthread<=0){
		warning("illigal nthread. set to 0.\n");
		arg->nthread=1;
	    }else if(arg->nthread>NCPU){
		warning("nthread is larger than number of cpus, reset to %d\n",NCPU);
		arg->nthread=NCPU;
	    }
	    break;
	case 'c':
	    arg->conf=strdup(optarg);
	    info2("Main config file changed to %s\n",arg->conf);
	    break;
	case 'p':{
	    addpath(optarg);
	}
	case '?':
	    warning("Unregonized option. exit.\n");exit(1);
	    break;
	default:
	    printf("?? getopt returned 0%o ??\n", c);exit(1);
	    break;
	}
    }
    if(!arg->detach){//foreground task will start immediately.
	arg->force=1;
    }
    arg->iconf=optind;
    arg->argc=argc;
    arg->argv=argv;
    if(!arg->dirout){
	arg->dirout=strtime();
    }
    if(!arg->conf){
	/*If -c is not specifid in path, will use nfiraos.conf*/
	arg->conf=strdup("maos.conf");
    }
    //Setup PATH and result directory
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
static void rename_file(int sig){
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
    case SIGINT: //Ctrl-C
    case SIGTERM:
    case SIGQUIT: //Ctrl-'\'
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
   Handles signals. We don't want to exit the simulation when SIGPIPE happens
   (when we are writing to closed sockets)
 */
void skyc_signal_handler(int sig){
    if(sig==SIGPIPE){
	warning3("Program received signal SIGPIPE, broken pipe.\n");
	return;
    }
    disable_signal_handler;
    rename_file(sig);//handles signal
    if(sig!=0){
	info2("Caught signal %d\n",sig);
	if(sig == SIGSEGV){
	    print_backtrace(0);
	}
	scheduler_finish(1);
	exit(sig);
    }
}

/**
   Add two PSDs that doesn't have the same frequency. the first column of each
   dmat is the frequency nu, and the second column is PSD*/
static dmat *add_psd_nomatch(dmat *psd1,dmat *psd2){
    dmat *nu1=dsub(psd1,0,psd1->nx,0,1);
    dmat *psd2x=dnew_ref(psd2->p, psd2->nx, 1);
    dmat *psd2y=dnew_ref(psd2->p+psd2->nx,psd2->nx,1);
    dmat *psd2new=dinterp1log(psd2x,psd2y,nu1);
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
	warning("The two PSDs have different length\n");
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
	    //todo: apply interp1 to interpolate the second PSD.
	}
	psd->p[i]=psd1->p[i];
	pp[i]=p1[i]+p2[i];
    }
    return psd;
}
