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
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE //DT_DIR #glibc > 2.19
#endif
#ifndef _BSD_SOURCE
#define _BSD_SOURCE //DT_DIR glibc < 2.19
#endif
#include <dirent.h>
#include <getopt.h>
#include <errno.h>
#include <unistd.h>
#include "../lib/aos.h"
static void usage(){
	fprintf(stderr, "Usage:\n"
		"drawres [-s 1] folder1 folder2 ...\n");
}
/**
   \file drawres.c
   Plot MAOS Results
*/
typedef struct arg_t{
	int iarg;
	int nseed;
	long* seeds;
}arg_t;
static arg_t* parse_args(int argc, char** argv){
	arg_t* arg=mycalloc(1, arg_t);
	static struct option long_options[]={
	{"help",0,0,'h'},
	{"seed",1,0,'s'},
	{NULL,0,0,0}
	};
	while(1){
		int option_index=0;
		int c=getopt_long(argc, argv, "hdfo:n:c:s:p:",
			long_options, &option_index);
		if(c==-1) break;
		switch(c){
		case 'h':
			usage();
			exit(0);
			break;
		case 's':{
			arg->nseed++;
			arg->seeds=myrealloc(arg->seeds, arg->nseed, long);
			arg->seeds[arg->nseed-1]=strtol(optarg, NULL, 10);
		}
				break;
		}
	}
	arg->iarg=optind;
	return arg;
}
void fixnan(dmat* res){
	for(long i=0; i<res->nx*res->ny; i++){
		if(isnan(P(res, i))){
			P(res, i)=0;
		}
	}
}
/**
   The main.
*/
int main(int argc, char* argv[]){
	int DRAWRES_HI=0;
	READ_ENV_INT(DRAWRES_HI, 0, 1);
	int drawres_tot=1;
	int drawres_lo=1;
	int drawres_hi=1;
	int drawres_ol=1;
	if(check_suffix(argv[0], "drawreshi")||DRAWRES_HI==1){
		drawres_lo=0;
		drawres_tot=0;
		drawres_ol=0;
		dbg("MAOS_DRAWRES_HI=1: Only plotting high order closed loop results.\n");
	} else{
		dbg("set MAOS_DRAWRES_HI=1 to only plot high order closed loop results.\n");
	}
	TIC;tic;
	arg_t* arg=parse_args(argc, argv);
	/*use the parent pid so same bash session has the same drawdaemon. */
	draw_id=DRAW_ID_RES;
	draw_direct=1;/*launch drawdaemon directly, without going through server. */
	draw_single=-1;//disable draw_single support.
	int has_ahst=0;//where there are ahst results
	char** path;
	int npath;
	if(arg->iarg<argc){
		npath=argc-arg->iarg;
		path=mycalloc(npath, char*);
		for(int ipath=0; ipath<npath; ipath++){
			path[ipath]=strdup(argv[ipath+arg->iarg]);
		}
	} else{
		npath=1;
		path=mycalloc(npath, char*);
		path[0]=mygetcwd();
	}
	int mpath=npath;
	long nseed=0;
	long* seed=NULL;
	long* seed2=NULL;
	/*Find all valid path and seeds. */
	int jpath=0;
	int restype=-1;//1: maos, 2: skyc
	for(int ipath=0; ipath<npath; ipath++){
		DIR* dir=opendir(path[ipath]);
		if(!dir){
			if(errno!=ENOTDIR){
				warning("Unable to read directory %s\n", path[ipath]);
			}
			continue;
		}
		struct dirent* dp;
		int nseedfound=0;
		while((dp=readdir(dir))){
			if(dp->d_type==DT_DIR){
				if(dp->d_name[0]!='.'&&strcmp(dp->d_name, "skysim")){
					if(npath==mpath){
						mpath*=2;
						path=myrealloc(path, mpath, char*);
					}
					path[npath]=stradd(path[ipath], "/", dp->d_name, NULL);
					npath++;
				}
			} else if(!strncmp(dp->d_name, "Res", 3)&&check_suffix(dp->d_name, ".bin")){
				long seedi, seedi2;
				int found=0;
				if(sscanf(dp->d_name, "Res_%ld.bin", &(seedi))==1){
					//MAOS results
					if(restype!=-1&&restype!=1){
						error("Only maos or skyc results can display together\n");
					}
					restype=1;
					if(arg->nseed){//filter seeds.
						int wanted=0;
						for(int is=0; is<arg->nseed; is++){
							if(arg->seeds[is]==seedi){
								wanted=1;
								break;
							}
						}
						if(!wanted){
							warning("Skip seed %ld\n", seedi);
							continue;
						}
					}
					for(int is=0; is<nseed; is++){
						if(seed[is]==seedi){
							found=1;
						}
					}
					if(!found){
						nseed++;
						seed=myrealloc(seed, nseed, long);
						seed[nseed-1]=seedi;
					}
					nseedfound++;
				} else if(sscanf(dp->d_name, "Res%ld_%ld.bin", &seedi, &seedi2)==2){
					//Sky coverage results
					if(restype!=-1&&restype!=2){
						warning("Only maos or skyc results can display together\n");
						continue;
					}
					restype=2;

					for(int is=0; is<nseed; is++){
						if(seed[is]==seedi&&seed2[is]==seedi2){
							found=1;
						}
					}
					if(!found){
						dbg("Found seed %ld %ld.\n", seedi, seedi2);
						nseed++;
						seed=myrealloc(seed, nseed, long);
						seed[nseed-1]=seedi;
						seed2=myrealloc(seed2, nseed, long);
						seed2[nseed-1]=seedi2;
					}
					nseedfound++;
				}
			}
		}
		closedir(dir);
		if(nseedfound){ /*record the path as valid */
			path[jpath]=path[ipath];
			info("Found path: %s\n", path[jpath]);
			jpath++;
		}else{
			free(path[ipath]);path[ipath]=NULL;
		}
	}
	npath=jpath;
	if(npath==0){
		info("Nothing to display\n");
		return 1;
	}

	info("Found seed:");
	for(int i=0; i<nseed; i++){
		info("%ld ", seed[i]);
	}
	info("\n");

	enum{
		P_TOT,
		P_HI,
		P_LO,
		P_TT,
		P_PS,
		P_F,
		P_OLTOT,
		P_OLHI,
		P_OLLO,
		N_ALL,
	};
	//name of the top tab.
	const char* toptab[]={
	"CL",
	"CL",
	"CL",
	"CL",
	"CL",
	"CL",
	"OL",//Openloop
	"OL",
	"OL"
	};
	//name of the side tab
	const char* sidetab[]={
	"Total",
	"High",
	"Low",
	"TT",
	"PS",
	"Focus",
	"Total",//OpenLoop
	"High",
	"Low"
	};
	//title of each plot
	const char* title[]={
	"Total Wavefront Error",
	"High Order Wavefront Error",
	"Low Order Wavefront Error",
	"Tip/Tilt Wavefront Error",
	"Plate Scale Wavefront Error",
	"Focus Mode Wavefront Error",
	"Total Wavefront Error",
	"High Order Wavefront Error",
	"Low Order Wavefront Error",
	};


	dccell* res=dccellnew(N_ALL, 1);
	dccell* resm=dccellnew(N_ALL, 1);
	for(int i=0; i<res->nx; i++){
		P(res, i)=dcellnew(npath, nseed);
		P(resm, i)=dcellnew(npath, 1);//seed averaged
	}

	const char* xlabel, * ylabel;
	dmat* ysky=NULL;
	const char* xylog;
	if(restype==1){
		xylog="nn";
		xlabel="Steps";
		ylabel="Wavefront Error (nm)";
	} else{
		xylog="yn";
		ylabel="Sky Coverage";
		xlabel="Wavefront Error (nm)";
	}

	/*dcell *upterr=dcellnew(npath, nseed);
	  dcell *uptcmd=dcellnew(npath, nseed);*/
	for(int ipath=0; ipath<npath; ipath++){
		for(int iseed=0; iseed<nseed; iseed++){
			char fn[PATH_MAX];
			switch(restype){
			case 1:
				snprintf(fn, PATH_MAX, "%s/Res_%ld.bin", path[ipath], seed[iseed]); break;
			case 2:
				snprintf(fn, PATH_MAX, "%s/Res%ld_%ld.bin", path[ipath], seed[iseed], seed2[iseed]); break;
			default:
				error("Invalid restype=%d\n", restype);
			}
			if(!zfexist("%s", fn)) continue;
			if(restype==1){//MAOS results.
				dcell* ires=dcellread("%s", fn);
				if(!ires||ires->nx<3||!ires->p){
					continue;
				}
				int indcl=0;
				int indlo=-1;
				int indhi=-1;
				int indtt=-1;
				int indfocus=-1;
				if(P(ires, 3)&&P(ires, 3)->nx>0){/*split tomography. */
					indcl=3;//cell index for CL results
					indhi=0;//high
					indtt=1;//tt
					indlo=2;//low order total
					if(P(ires, 3)->nx>3){
						indfocus=3;
					}
				} else{
					indcl=2;//cell index for CL results
					//0 is pr.
					indlo=1;//low order
					indtt=1;//tt
					indhi=2;//pttr 
				}
				long nstep=NY(P(ires, 0));//valid step;
				for(long i=0; i<nstep; i++){
					real val=P(P(ires, 0), 0, i);
					if(isnan(val)||val==0){
						//warning("Number of steps reduced from %ld to %ld\n", nstep, i);
						nstep=i;
						break;
					}
				}
				dmat* tmp;
				if(indhi!=-1 && (drawres_hi||drawres_tot)){
					tmp=dsub(P(ires, indcl), indhi, 1, 0, nstep);
					P(P(res, P_HI), ipath, iseed)=dtrans(tmp);
					dfree(tmp);
				}
				if(indlo!=-1 && (drawres_lo||drawres_tot)){
					tmp=dsub(P(ires, indcl), indlo, 1, 0, nstep);
					fixnan(tmp);
					P(P(res, P_LO), ipath, iseed)=dtrans(tmp);
					dfree(tmp);
				}
				if(drawres_tot){
					dadd(&P(P(res, P_TOT), ipath, iseed), 1, P(P(res, P_LO), ipath, iseed), 1);
					dadd(&P(P(res, P_TOT), ipath, iseed), 1, P(P(res, P_HI), ipath, iseed), 1);
					if(!drawres_hi){
						dfree(P(P(res, P_HI), ipath, iseed));
					}
					if(!drawres_lo){
						dfree(P(P(res, P_LO), ipath, iseed));
					}
				}
				if(drawres_lo){
					if(indfocus>-1){
						tmp=dsub(P(ires, indcl), indfocus, 1, 0, nstep);
						fixnan(tmp);
						P(P(res, P_F), ipath, iseed)=dtrans(tmp);
						dfree(tmp);
					}
					if(indtt>-1){
						tmp=dsub(P(ires, indcl), indtt, 1, 0, nstep);
						fixnan(tmp);
						P(P(res, P_TT), ipath, iseed)=dtrans(tmp);
						dfree(tmp);
					}
					if(indfocus>-1 && indtt>-1){
						dadd(&P(P(res, P_PS), ipath, iseed), 1, P(P(res, P_LO), ipath, iseed), 1);
						dadd(&P(P(res, P_PS), ipath, iseed), 1, P(P(res, P_TT), ipath, iseed), -1);
						dadd(&P(P(res, P_PS), ipath, iseed), 1, P(P(res, P_F), ipath, iseed), -1);
						if(dmax(P(res, P_PS, ipath, iseed))<1e-24){//<1e-3 nm
							dfree(P(res, P_PS, ipath, iseed));
						}
					}
				}
				if(indlo!=indtt){
					has_ahst++;
				}
				if(drawres_ol){
					tmp=dsub(P(ires, 0), 2, 1, 0, nstep);
					P(P(res, P_OLHI), ipath, iseed)=dtrans(tmp);
					dfree(tmp);
					tmp=dsub(P(ires, 0), 1, 1, 0, nstep);
					P(P(res, P_OLLO), ipath, iseed)=dtrans(tmp);
					dfree(tmp);
					dadd(&P(P(res, P_OLTOT), ipath, iseed), 1, P(P(res, P_OLLO), ipath, iseed), 1);
					dadd(&P(P(res, P_OLTOT), ipath, iseed), 1, P(P(res, P_OLHI), ipath, iseed), 1);
				}
				dcellfree(ires);
			} else if(restype==2){//Skycoverage results
				snprintf(fn, PATH_MAX, "%s/Res%ld_%ld.bin", path[ipath], seed[iseed], seed2[iseed]);
				dmat* res0=dread("%s", fn);
				if(NY(res0)==0){
					dbg("Empty result: %s, skip.\n", fn);
					continue;
				}
				if(ysky && NX(ysky)!=NY(res0)){
					dfree(ysky);
				}
				if(!ysky){
					int nsky=NY(res0);
					ysky=dnew(nsky, 1);
					for(int i=0; i<nsky; i++){
						P(ysky, i)=pow((double)i/(double)(nsky-1)*1e-9, 2);
					}
				}
				dmat *tmp=dsub(res0, 0, 1, 0, 0);//keep only first row
				tmp->nx=NY(tmp); tmp->ny=1; //make it column vector.
				dsort(tmp, 1);
				P(P(res, P_TOT), ipath, iseed)=dcat(tmp, ysky, 2);
				dfree(tmp);

				/*tmp=dsub(ires, 2, 1, 0, 0);
				tmp->nx=NY(tmp); tmp->ny=1; //make it column vector.
				dsort(tmp, 1);
				P(P(res,P_LO),ipath, iseed)=dcat(tmp, ysky, 2);
				dfree(tmp);*/
			} else{
				error("Invalid restype=%d\n", restype);
			}
		}//for seed
	}//for path
	for(int imod=0; imod<NX(res); imod++){//mode
		for(int ipath=0; ipath<npath; ipath++){
			int seedcount=0;
			int nsim=0;
			for(int iseed=0; iseed<nseed; iseed++){
				dmat* resi=P(P(res, imod), ipath, iseed);
				if(resi && NX(resi)){
					dadd_relax(&P(P(resm, imod), ipath), 1, resi, 1);
					nsim=MIN(NX(P(P(resm, imod), ipath)), NX(resi));
					seedcount++;
				}
			}
			if(seedcount>0){
				dscale(P(P(resm, imod), ipath), 1./seedcount);
				dresize(P(P(resm, imod), ipath), nsim, NY(P(resm, imod), ipath));
			}
		}
	}
	if(restype==1 && has_ahst==0){//no ahst results, TT and LOW are the same. 
		dcellfree(P(res, P_LO));
	}
	int nvalid=0;//number of valid modes in CL type
	for(int imod=0; imod<NX(res); imod++){
		if(dcellsum(P(res, imod))==0){
			dcellfree(P(res, imod));
			dcellfree(P(resm, imod));
		} else{
			dcellcwpow(P(res, imod), 0.5);
			dcellscale(P(res, imod), 1e9);
			dcellcwpow(P(resm, imod), 0.5);
			dcellscale(P(resm, imod), 1e9);
			if(toptab[imod][0]=='C'){
				nvalid++;
			}
		}
	}
	char* pathtag0[npath];
	char prefix[5]="A: ";
	for(int ipath=0; ipath<npath; ipath++){
		if(ipath<26){
			prefix[0]='A'+ipath;
		}else{
			prefix[0]='A'+ipath/26;
			prefix[1]='A'+ipath%26;
			prefix[2]=':';
		}
		pathtag0[ipath]=stradd(prefix, path[ipath], NULL);
	}
	if(nseed>1){//more than one seed is available. Every over the seeds
		if(npath==1 && nvalid>1){//plot different seeds when only 1 path is found.
			char* legs0[nseed+1];
			for(int iseed=0; iseed<nseed; iseed++){
				legs0[iseed]=mymalloc(50, char);
				snprintf(legs0[iseed], 50, "Seed %ld", seed[iseed]);
			}
			legs0[nseed]=mystrdup("Seed RMS");

			for(int ic=0; ic<res->nx; ic++){
				if(P(res, ic)){
					if(nseed>1){//seed average
						cellresize(P(res, ic), 1, nseed+1);
						P(P(res, ic), nseed)=dref(P(P(resm, ic), 0));
					}
					draw(toptab[ic], (plot_opts){.dc=P(res, ic), .xylog=xylog, .legend=(const char* const*)legs0},
						title[ic], xlabel, ylabel, "%s", sidetab[ic]);
					if(nseed>1){//seed average
						dfree(P(P(res, ic), nseed));
						cellresize(P(res, ic), 1, nseed);
					}
				}
			}

			for(int iseed=0; iseed<=nseed; iseed++){
				free(legs0[iseed]);
			}
			{//plot all modes together for seed average
				dcell *restmp=dcellnew(N_ALL, 1);
				for(int ic=0; ic<restmp->nx; ic++){
					if(P(resm, ic)){
						P(restmp, ic, 0)=dref(P(P(resm, ic),0));
					}
				}
				dcell *restmpcl=dcellsub(restmp, 0, 6, 0, 1);
				draw("CL", (plot_opts){ .ngroup=(int)restmpcl->nx, .dc=restmpcl, .xylog=xylog, .legend=(const char *const *)sidetab },
					"Closed loop Wavefront Error", xlabel, ylabel, "All");
				dcellfree(restmpcl);
				if(restype==1){
					dcell *restmpol=dcellsub(restmp, 6, 3, 0, 1);
					draw("OL", (plot_opts){ .ngroup=(int)restmpol->nx, .dc=restmpol, .xylog=xylog, .legend=(const char *const *)sidetab },
						"Open loop Wavefront Error", xlabel, ylabel, "All");
						dcellfree(restmpol);
				}
				
				dcellfree(restmp);
			}
		} else{//plot seed average for each mode.
			for(int ic=0; ic<resm->nx; ic++){
				if(P(resm, ic)){
					draw(toptab[ic], (plot_opts){.ngroup=(int)npath,.dc=P(resm, ic), .xylog=xylog, .legend=(const char* const*)pathtag0},
						title[ic], xlabel, ylabel, "%s", sidetab[ic]);
				}
			}
		}
	}
	for(int iseed=0; iseed<nseed; iseed++){//each seed
		if(npath==1 && nvalid>1){//plot all modes together for each seed
			dcell *restmp=dcellnew(N_ALL,1);
			for(int ic=0; ic<restmp->nx; ic++){
				if(P(res, ic)){
					P(restmp, ic, 0)=dref(P(P(res, ic),0,iseed));
				}
			}
			dcell* restmpcl=dcellsub(restmp, 0, 6, 0, 1);
			draw("CL", (plot_opts){ .ngroup=(int)restmpcl->nx, .dc=restmpcl, .xylog=xylog, .legend=(const char *const *)sidetab },
				"Closed loop Wavefront Error", xlabel, ylabel, "%s:%-4ld", "All", seed[iseed]);
			dcellfree(restmpcl);
			if(restype==1){
				dcell *restmpol=dcellsub(restmp, 6, 3, 0, 1);
				draw("OL", (plot_opts){ .ngroup=(int)restmpol->nx, .dc=restmpol, .xylog=xylog, .legend=(const char *const *)sidetab },
				"Open loop Wavefront Error", xlabel, ylabel, "%s:%-4ld", "All", seed[iseed]);
				dcellfree(restmpol);
			}
			dcellfree(restmp);
		}else{//plot all path together for each seed
			for(int ic=0; ic<res->nx; ic++){
				if(P(res, ic)){
					dcell* tmp=dcellsub(P(res, ic), 0, 0, iseed, 1);
					draw(toptab[ic], (plot_opts){.ngroup=(int)npath, .dc=tmp, .xylog=xylog, .legend=(const char* const*)pathtag0},
						title[ic], xlabel, ylabel, "%s:%-4ld", sidetab[ic], seed[iseed]);
					dcellfree(tmp);
				}
			}
		}
	}

	draw_final(1);
	cellfree(ysky);
	cellfree(res);
	cellfree(resm);
	for(int ipath=0; ipath<npath; ipath++){
		free(path[ipath]);
		free(pathtag0[ipath]);
	}
	free(path);
	free(seed);
	free(arg->seeds);
	free(arg);
	toc("drawres");
}
