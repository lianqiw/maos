/**
   Wrap of maos in mex so we can run full simulations in matlab, interrupt
   simulation, and examine data.

   maos() is the single entrace. Different functionalities are provided by the
   first parameter as a string.
*/

#include "interface.h"
#include "../maos/maos.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static PARMS_T *parms=0;
    static POWFS_T *powfs=0;
    static RECON_T *recon=0;
    static APER_T  *aper=0;
    static int iseed=0;
    static int iseed_last=-1;
    static int isim=-1;
    static SIM_T *simu=0;    
    char *cmd=" ";
    if(nrhs>0){
	cmd=mxArrayToString(prhs[0]);
    }
    if(!strcmp(cmd, "reset") && parms){
	maos_reset();
	parms=0;
	aper=0;
	powfs=0;
	recon=0;
	iseed=0;
	isim=-1;
	simu=0;
	return;
    }
    if(!strcmp(cmd, "setup") || !parms){
	int override=0;
	char *conf=0;
	char *dirout=0;
	char *mainconf=0;
	int nthread=0;
	int *gpus=0;
	int ngpu=0;
	int ngpu2=0;
	if(nrhs>1){
	    conf=mxArrayToString(prhs[1]);
	    ARGOPT_T options[]={
		{"override",'O',T_INT,0, &override, NULL},
		{"output", 'o',T_STR, 1, &dirout, NULL},
		{"nthread",'n',T_INT, 1, &nthread,NULL},
		{"conf",   'c',T_STR, 1, &mainconf, NULL},
		{"path",   'P',T_STR, 3, addpath, NULL},
		{"gpu",    'g',T_INTARR, 1, &gpus, &ngpu},
		{"ngpu",   'G',T_INT, 1, &ngpu2, NULL},
		{NULL, 0,0,0, NULL, NULL}
	    };
	    parse_argopt(conf, options);
	    if(dirout){
		mymkdir("%s",dirout);
		if(chdir(dirout)){
		    error("Unable to chdir to %s\n", dirout);
		}
	    }
	    if(nthread<NTHREAD && nthread>0){
		NTHREAD=nthread;
	    }
	    if(ngpu2>0){
		if(!gpus || ngpu==0){
		    ngpu=ngpu2;
		}else{
		    error("-g and -G cannot both be specified\n");
		}
	    }
	}
	if(parms){
	    fprintf(stderr, "Already setup. Please call 'reset' first\n");
	    return;
	}
	addpath(".");
	parms=setup_parms(mainconf, conf, override);
	setup_parms_gpu(parms, gpus, ngpu);
	maos_setup(parms);
	powfs=global->powfs;
	recon=global->recon;
	aper=global->aper;
    }
    int nstep=1;
    if(!strcmp(cmd, "sim")){
	if(nrhs>1){
	    nstep=(int)mxGetScalar(prhs[1]);
	    if(nstep<=0){
		nstep=parms->sim.end-parms->sim.start;
	    }
	}
    }
    if(isim==-1){
	isim=parms->sim.start;
    }
    if(iseed<parms->sim.nseed){
	if(!simu){
	    while(!(simu=maos_iseed(iseed))){
		iseed++;
		if(iseed==parms->sim.nseed){
		    info2("All seeds are finished\n");
		    return;
		}
	    }
	}
	while(nstep--){
	    if(isim<parms->sim.end){
		maos_isim(isim);
		isim++;
	    }else{//one seed finished
		free_simu(simu);simu=0;global->simu=0;
		iseed++;
		isim=parms->sim.start;
		break;
	    }
	}
    }else{
	info2("Simulation finished\n");
    }
}
