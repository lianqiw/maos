/**
   Wrap of maos in mex so we can run full simulations in matlab, interrupt
   simulation, and examine data.

   maos() is the single entrace. Different functionalities are provided by the
   first parameter as a string.
*/

#include "interface.h"
#include "../maos/maos.h"
#include "maos2mex.h"
static __attribute__((destructor)) void deinit_maos(){
    maos_reset();
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    exception_env=malloc(sizeof(jmp_buf));
    if(setjmp(*exception_env)){
	//We use longjump because calling mexErrMsgTxt causing matlab to crash (bug?)
	info2("Exception happened\n");
	return;
    }else{
	const PARMS_T *parms=0;
	static int isim=0;
	static int iseed=0;
	int nstep=0;
	char *cmd=0;//default action is sim
	int free_cmd=0;
	if(nrhs>0){
	    cmd=mxArrayToString(prhs[0]);
	    free_cmd=1;
	}else{
	    cmd="sim";
	}
	if(!strcmp(cmd, "reset")){
	    if(global) maos_reset();
	    iseed=0;
	    isim=0;
	    goto end;
	}
	if(!strcmp(cmd, "setup") || !global){
	    if(global){
		fprintf(stderr, "Already setup. Please call 'reset' first\n");
		goto end;
	    }
	    int override=0;
	    char *conf=0;
	    char *dirout=0;
	    char *mainconf=0;
	    int nthread=0;
	    int *gpus=0;
	    int ngpu=0;
	    int ngpu2=0;
	    if(!strcmp(cmd, "setup") && nrhs>1){
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
		if(nrhs>2){
		    if(!mxIsDouble(prhs[2])){
			error("The second parameter should be an integer\n");
		    }
		    nstep=(int)mxGetScalar(prhs[2]);
		}
	    }
	    addpath(".");
	    if(dirout){
		info2("dirout=%s\n", dirout);
		mymkdir("%s",dirout);
		if(chdir(dirout)){
		    error("Unable to chdir to %s\n", dirout);
		}
	    }else{
		warning2("Disable saving when no -o is supplied.\n");
		disable_save=1;
	    }
	    parms=setup_parms(mainconf, conf, override);
	    info2("setup_parms done\n");
	    setup_parms_gpu((PARMS_T*)parms, gpus, ngpu);
	    maos_setup(parms);//sets global
	}
	parms=global->parms;
	if(!strcmp(cmd, "sim")){
	    if(nrhs>1){
		if(!mxIsDouble(prhs[1])){
		    error("The second parameter should be an integer\n");
		}
		nstep=(int)mxGetScalar(prhs[1]);
	    }
	}
	if(nstep<0){
	    nstep=parms->sim.end-parms->sim.start;
	}
	if(nstep>0){
	    SIM_T *simu=global->simu;
	    if(iseed<parms->sim.nseed){
		if(!simu){
		    while(!(simu=maos_iseed(iseed))){
			iseed++;
			if(iseed==parms->sim.nseed){
			    info2("All seeds are finished\n");
			    goto end;
			}
		    }
		    isim=parms->sim.start;
		}
		while(nstep--){
		    if(isim<parms->sim.end){
			maos_isim(isim);
			isim++;
		    }else{//one seed finished
			free_simu(simu);simu=0;
			iseed++;
			break;
		    }
		    extern int utIsInterruptPending();
		    if(utIsInterruptPending()){
			info2("Simulation interrupted\n");
			goto end;
		    }
		}
	    }else{
		info2("Simulation finished\n");
	    }
	}
	int free_valname=0;
	char *valname=NULL;
	if(!strcmp(cmd, "get")){
	    if(nrhs>1){
		free_valname=1;
		valname=mxArrayToString(prhs[1]);
	    }else{
		valname="";
	    }
	}else if(nlhs>0){
	    valname="simu";
	}
	if(valname){
	    if(global){
		plhs[0]=get_data(global->simu, valname);
	    }else{
		plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
	    }
	    if(free_valname) mxFree(valname);
	}
      end:;
	if(free_cmd) mxFree(cmd);
    }
}
    
