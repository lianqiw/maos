/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <unistd.h>
#include "interface.h"
#include "../maos/maos.h"
#include "maos2mex.h"
static __attribute__((constructor)) void init_maos(){
	register_deinit(maos_reset, NULL);
}
#ifdef __cplusplus
extern "C" {
#endif
	int utIsInterruptPending();
#ifdef __cplusplus
}
#endif
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	//exception_env=(jmp_buf*)malloc(sizeof(jmp_buf));
	/*if(setjmp(*exception_env)){
		//We use longjump because calling mexErrMsgTxt causing matlab to crash (bug?)
		printf("Exception happened\n");
		return;
	} else{*/
		static int isim=0;
		static int iseed=0;
		int nstep=0;
		const char* cmd="help";//default action is sim
		int free_cmd=0;
		if(nrhs>0&&mxIsChar(prhs[0])){
			cmd=mxArrayToString(prhs[0]);
			free_cmd=1;
		} else{
			cmd="sim";
			if(nrhs>0&&mxIsDouble(prhs[0])){
				nstep=(int)mxGetScalar(prhs[0]);
			}
		}
		if(nrhs==0||(strcmp(cmd, "setup")&&strcmp(cmd, "reset")&&strcmp(cmd, "sim")&&strcmp(cmd, "get"))){
			printf("Usage: \n"
				"       maos('setup', '-o dirout -n N -c scao_ngs.conf -g0')\n"
				"       maos('reset')\n"
				"  simu=maos('sim', nstep)\n"
				"  simu=maos('get','sim')\n"
				" parms=maos('get','parms')\n"
			);
			return;
		}
		if(!strcmp(cmd, "reset")){
			if(global) maos_reset();
			iseed=0;
			isim=0;
			goto end;
		}
		if(!strcmp(cmd, "setup")||!global){
			if(global){
				fprintf(stderr, "Already setup. Please call 'reset' first\n");
				goto end;
			}
			int over_ride=0;
			char* conf=0;
			char* dirout=0;
			char* mainconf=0;
			int nthread=0;
			int* gpus=0;
			int ngpu=0;
			int ngpu2=0;

			if(!strcmp(cmd, "setup")&&nrhs>1){
				conf=mxArrayToString(prhs[1]);
				argopt_t options[]={
					{"override",'O',M_INT,0, 0, &over_ride, NULL},
					{"output", 'o',M_STR, 1, 0, &dirout, NULL},
					{"nthread",'n',M_INT, 1, 0, &nthread,NULL},
					{"gpu",    'g',M_INT, 2, 0, &gpus, &ngpu},
					{"ngpu",   'G',M_INT, 1, 0, &ngpu2, NULL},
					{"conf",   'c',M_STR, 1, 0, &mainconf, NULL},
					{"path",   'P',M_STR, 1, 1, (void*)addpath, NULL},
					{NULL,     0,  0,     0, 0, NULL, NULL}
				};
				parse_argopt(conf, options);
				if(nthread<MAXTHREAD && nthread>0){
					NTHREAD=nthread;
				}
				if(ngpu2>0){
					if(!gpus||ngpu==0){
						ngpu=ngpu2;
					} else{
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
				printf("dirout=%s\n", dirout);
				mymkdir("%s", dirout);
				if(chdir(dirout)){
					error("Unable to chdir to %s\n", dirout);
				}
			} else{
				warning("Disable saving when no -o is supplied.\n");
				disable_save=1;
			}
			const parms_t *parms=setup_parms(mainconf, conf, over_ride);
			printf("setup_parms done\n");
			setup_parms_gpu((parms_t*)parms, gpus, ngpu);
			maos_setup(parms);//sets global
		}
		if(!strcmp(cmd, "sim")){
			const parms_t *parms=global->parms;
			if(nrhs>1){
				if(!mxIsDouble(prhs[1])){
					error("The second parameter should be an integer\n");
				}
				nstep=(int)mxGetScalar(prhs[1]);
				if(nstep<0){
					nstep=parms->sim.end-parms->sim.start;
				}
			} else if(!nstep){
				nstep=1;
			}

			if(nstep>0){
				sim_t* simu=global->simu;
				if(iseed<parms->sim.nseed){
					if(!simu){
						while(!(simu=maos_iseed(iseed))){
							iseed++;
							if(iseed==parms->sim.nseed){
								printf("All seeds are finished\n");
								goto end;
							}
						}
						isim=parms->sim.start;
					}
					while(nstep--){
						if(isim<parms->sim.end){
							maos_isim(isim);
							isim++;
						} else{//one seed finished
							free_simu(simu);simu=0;
							iseed++;
							break;
						}
					}
				} else{
					printf("Simulation finished\n");
				}
			}
		}
		if(nlhs==1){
			int free_valname=0;
			const char* valname=NULL;
			if(!strcmp(cmd, "get")){
				if(nrhs>1){
					free_valname=1;
					valname=mxArrayToString(prhs[1]);
				} else{
					valname="";
				}
			} else if(nlhs>0){
				valname="sim";
			}
			if(valname){
				if(global){
					plhs[0]=get_data(global->simu, valname);
				} else{
					plhs[0]=mxCreateDoubleMatrix(0, 0, mxREAL);
				}
				if(free_valname) mxFree((char*)valname);
			}
		}
end:;
		if(free_cmd) mxFree((char*)cmd);
	//}
}

