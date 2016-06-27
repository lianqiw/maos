/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifdef __cplusplus
extern "C" {
#endif
    int utIsInterruptPending();
#ifdef __cplusplus
}
#endif
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    exception_env=(int(*)[37])malloc(sizeof(jmp_buf));
    if(setjmp(*exception_env)){
	//We use longjump because calling mexErrMsgTxt causing matlab to crash (bug?)
	info2("Exception happened\n");
	return;
    }else{
	const PARMS_T *parms=0;
	static int isim=0;
	static int iseed=0;
	int nstep=0;
	const char *cmd=0;//default action is sim
	int free_cmd=0;
	if(nrhs>0 && mxIsChar(prhs[0])){
	    cmd=mxArrayToString(prhs[0]);
	    free_cmd=1;
	}else{
	    cmd="sim";
	    if(nrhs>0 && mxIsDouble(prhs[0])){
		nstep=(int)mxGetScalar(prhs[0]);
	    }
	}
	if(!strcmp(cmd, "-h")){
	    printf("Usage: maos('setup', '-o dirout -n N -c scao_ngs.conf -g0')\n"
		 "  simu=maos('sim', nstep)\n"
		 "       maos('reset')\n"
		 "  simu=maos('get','simu')\n"
		 "dmreal=maos('get','dmreal')\n"
		);
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
		    {"override",'O',M_INT,0, 0, &override, NULL},
		    {"output", 'o',M_STR, 1, 0, &dirout, NULL},
		    {"nthread",'n',M_INT, 1, 0, &nthread,NULL},
		    {"gpu",    'g',M_INT, 2, 0, &gpus, &ngpu},
		    {"ngpu",   'G',M_INT, 1, 0, &ngpu2, NULL},
		    {"conf",   'c',M_STR, 1, 0, &mainconf, NULL},
		    {"path",   'P',M_STR, 1, 1, (void*)addpath, NULL},
		    {NULL,     0,  0,     0, 0, NULL, NULL}
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
		if(nstep<0){
		    nstep=parms->sim.end-parms->sim.start;
		}
	    }else{
		nstep=1;
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
			if(utIsInterruptPending()){
			    info2("Simulation interrupted\n");
			    goto end;
			}
		    }
		}else{
		    info2("Simulation finished\n");
		}
	    }
	}
	if(nlhs==1){
	    int free_valname=0;
	    const char *valname=NULL;
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
		if(free_valname) mxFree((char*)valname);
	    }
	}
      end:;
	if(free_cmd) mxFree((char*)cmd);
    }
}
    
