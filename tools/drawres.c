/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <dirent.h>
#include <getopt.h>

#include "../lib/aos.h"
static void usage(){
    fprintf(stderr,"Usage:\n"
	    "drawres [-s 1] folder1 folder2 ...\n");
}
/**
   \file drawres.c
   Plot MAOS Results
*/
typedef struct ARG_T{
    int iarg;
    int nseed;
    long *seeds;
}ARG_T;
static ARG_T * parse_args(int argc, char **argv){
    ARG_T *arg=mycalloc(1,ARG_T);
    static struct option long_options[]={
	{"help",0,0,'h'},
	{"seed",1,0,'s'},
	{NULL,0,0,0}
    };
    while(1){
	int option_index = 0;
	int c = getopt_long(argc, argv, "hdfo:n:c:s:p:",
			    long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'h':
	    usage();
	    exit(0);
	    break;
	case 's':{
	    arg->nseed++;
	    arg->seeds=myrealloc(arg->seeds,arg->nseed,long);
	    arg->seeds[arg->nseed-1]=strtol(optarg,NULL,10);
	}
	    break;
	}
    }
    arg->iarg=optind;
    return arg;
}

/**
   The main.
*/
int main(int argc, char *argv[]){
    ARG_T *arg=parse_args(argc, argv);
    /*use the parent pid so same bash session has the same drawdaemon. */
    DRAW_ID=getsid(0)+2e6;/*variables in draw.c */
    DRAW_DIRECT=1;/*launch drawdaemon directly, without going through server. */
    char **path;
    int npath;

    if(arg->iarg<argc){
	npath=argc-arg->iarg;
	path=mycalloc(npath,char*);
	for(int ipath=0; ipath<npath; ipath++){
	    path[ipath]=argv[ipath+arg->iarg];
	}
    }else{
	npath=1;
	path=mycalloc(npath,char*);
	path[0]=mygetcwd();
    }
    long nseed=0;
    long *seed=NULL;
    long *seed2=NULL;
    /*Find all valid path and seeds. */
    int jpath=0;
    int restype=-1;
    for(int ipath=0; ipath<npath; ipath++){
	if(!isdir(path[ipath])){
	    continue;
	}
	/*dbg("Try path %s\n", path[ipath]); */
	DIR *dir=opendir(path[ipath]);
	if(!dir){
	    warning("Unable to read directory %s\n", path[0]);
	    continue;
	}
	struct dirent *dp;
	int nseedfound=0;
	while((dp=readdir(dir))){
	    if(!strncmp(dp->d_name, "Res", 3) && check_suffix(dp->d_name, ".bin")){
		long seedi, seedi2;
		int found=0;
		if(sscanf(dp->d_name, "Res_%ld.bin",&(seedi))==1){
		    //MAOS results
		    if(restype!=-1 && restype!=1){
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
			if(!wanted) {
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
			seed=myrealloc(seed,nseed,long);
			seed[nseed-1]=seedi;
		    }
		    nseedfound++;
		}else if(sscanf(dp->d_name, "Res%ld_%ld.bin",&seedi, &seedi2)==2){
		    //Sky coverage results
		    if(restype!=-1 && restype!=2){
			error("Only maos or skyc results can display together\n");
		    }
		    restype=2;
		    
		    for(int is=0; is<nseed; is++){
			if(seed[is]==seedi && seed2[is]==seedi2){
			    found=1;
			}
		    }
		    if(!found){
			dbg("Found seed %ld %ld.\n", seedi, seedi2);
			nseed++;
			seed=myrealloc(seed,nseed,long);
			seed[nseed-1]=seedi;
			seed2=myrealloc(seed2,nseed,long);
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
	P_OLTOT,
	P_OLHI,
	P_OLLO,
	P_TOT,
	P_HI,
	P_LO,
	P_TT,
	P_PS,
	P_F,
	N_ALL,
    };
    const char *toptab[]={
	"OL",
	"OL hi",
	"OL lo",
	"CL",
	"CL hi",
	"CL lo",
	"CL lo",
	"CL lo",
	"CL lo"
    };
    const char *sidetab[]={
	"Total",
	"High",
	"Low",
	"Total",
	"High",
	"Low",
	"Low_TT",
	"Low_PS",
	"Low_Focus"
    };
    const char *title[]={
	"Total Wavefront Error",
	"High Order Wavefront Error",
	"Low Order Wavefront Error",
	"Total Wavefront Error",
	"High Order Wavefront Error",
	"Low Order Wavefront Error",
	"Tip/Tilt Wavefront Error",
	"Plate Scale Wavefront Error",
	"Focus Mode Wavefront Error",
    };
    
	
    dccell *res=dccellnew(N_ALL,1);
    dccell *resm=dccellnew(N_ALL,1);
    for(int i=0; i<res->nx; i++){
	res->p[i]=dcellnew(npath,nseed);
	resm->p[i]=dcellnew(npath,1);
    }
  
    const char *xlabel, *ylabel;
    dmat *ysky=NULL;
    const char* xylog;
    if(restype==1){
	xylog="nn";
	xlabel="Steps";
	ylabel="Wavefront Error (nm)";
    }else{
	xylog="yn";
	ylabel="Sky Coverage";
	xlabel="Wavefront Error (nm)";
    }

    /*dcell *upterr=dcellnew(npath, nseed);
      dcell *uptcmd=dcellnew(npath, nseed);*/
    for(int ipath=0; ipath<npath; ipath++){
	int seedcount=0;
	for(int iseed=0; iseed<nseed; iseed++){
	    char fn[PATH_MAX];
	    switch(restype){
	    case 1:
		snprintf(fn,PATH_MAX,"%s/Res_%ld.bin", path[ipath], seed[iseed]); break;
	    case 2:
		snprintf(fn,PATH_MAX,"%s/Res%ld_%ld.bin", path[ipath], seed[iseed], seed2[iseed]); break;
	    default:
		error("Invalid restype=%d\n", restype);
	    }
	    if(!zfexist(fn)) continue;
	    int ii=ipath+npath*iseed;
	    if(restype==1){//MAOS results.
		dcell *ires;
		ires=dcellread("%s",fn);
		if(ires->nx<3 || !ires->p){
		    continue;
		}
		int ind=0;
		int indlo=0;
		int indhi=0;
		int indtt=-1;
		int indfocus=-1;
		if(ires->p[3] && ires->p[3]->nx>0){/*split tomography. */
		    ind=3;
		    indlo=2;//total ngs
		    indhi=0;//high
		    indtt=1;//tt
		    if(ires->p[3]->nx>3){
			indfocus=3;
		    }
		}else{
		    ind=2;
		    indlo=1;/*tt */
		    indhi=2;/*pttr */
		}
		dmat *tmp;
		tmp=dsub(ires->p[ind], indhi, 1, 0, 0);
		res->p[P_HI]->p[ii]=dtrans(tmp);
		dfree(tmp);
		tmp=dsub(ires->p[ind], indlo, 1, 0, 0);
		res->p[P_LO]->p[ii]=dtrans(tmp);
		dfree(tmp);
		dadd(&res->p[P_TOT]->p[ii], 1, res->p[P_LO]->p[ii], 1);
		dadd(&res->p[P_TOT]->p[ii], 1, res->p[P_HI]->p[ii], 1);
		    
		if(indfocus>-1){
		    tmp=dsub(ires->p[ind], indfocus, 1, 0, 0);
		    res->p[P_F]->p[ii]=dtrans(tmp);
		    dfree(tmp);
		}
		if(indtt>-1){
		    tmp=dsub(ires->p[ind], indtt, 1, 0, 0);
		    res->p[P_TT]->p[ii]=dtrans(tmp);
		    dfree(tmp);
		    dadd(&res->p[P_PS]->p[ii], 1, res->p[P_LO]->p[ii], 1);
		    dadd(&res->p[P_PS]->p[ii], 1, res->p[P_TT]->p[ii], -1);
		    dadd(&res->p[P_PS]->p[ii], 1, res->p[P_F]->p[ii], -1);		    
		}
		
		tmp=dsub(ires->p[0], 2, 1, 0, 0);
		res->p[P_OLHI]->p[ii]=dtrans(tmp);
		dfree(tmp);
		tmp=dsub(ires->p[0], 1, 1, 0, 0);
		res->p[P_OLLO]->p[ii]=dtrans(tmp);
		dfree(tmp);
		dcellfree(ires);
		dadd(&res->p[P_OLTOT]->p[ii], 1, res->p[P_OLLO]->p[ii], 1);
		dadd(&res->p[P_OLTOT]->p[ii], 1, res->p[P_OLHI]->p[ii], 1);
	    }else if(restype==2){//Skycoverage results
		snprintf(fn,PATH_MAX,"%s/Res%ld_%ld.bin", path[ipath], seed[iseed], seed2[iseed]);
		dmat *res0=dread("%s", fn);
		dmat *ires=dtrans(res0); dfree(res0);
		dmat *tmp;
		if(!ysky){
		    int nsky=ires->nx;
		    ysky=dnew(nsky,1);
		    for(int i=0; i<nsky; i++){
			ysky->p[i]=pow((double)i/(double)(nsky-1)*1e-9,2);
		    }
		}
		
		if(ires->nx!=ysky->nx){
		    warning("Mismatch: %ld vs %ld\n", ires->nx, ysky->nx);
		    dfree(ires);
		    continue;
		}
		tmp=dsub(ires, 0, 0, 0, 1);
		dsort(tmp, 1);
		res->p[P_TOT]->p[ii]=dcat(tmp, ysky, 2); 
		dfree(tmp);

		/*tmp=dsub(ires, 0, 0, 2, 1);
		dsort(tmp, 1);
		res->p[P_LO]->p[ii]=dcat(tmp, ysky, 2); 
		dfree(tmp);*/
		
		dfree(ires);
	    }else{
		error("Invalid restype=%d\n", restype);
	    }
	    for(int i=0; i<res->nx; i++){
		dadd_relax(&resm->p[i]->p[ipath], 1, res->p[i]->p[ii], 1);
	    }
	  
	    seedcount++;
	}
	if(seedcount>0){
	    for(int i=0; i<res->nx; i++){
		dscale(resm->p[i]->p[ipath], 1./seedcount);
	    }
	  
	}
    }
  
    for(int i=0; i<res->nx; i++){
	dcellcwpow(res->p[i], 0.5);
	dcellscale(res->p[i], 1e9);
	dcellcwpow(resm->p[i], 0.5);
	dcellscale(resm->p[i], 1e9);
    }
   
 
    if(npath==1){
	char *legs0[nseed];
	const char *legs[nseed];//work around C warning.
	for(int iseed=0; iseed<nseed; iseed++){
	    legs0[iseed]=mymalloc(50,char);
	    snprintf(legs0[iseed], 50, "Seed %ld", seed[iseed]);
	    legs[iseed]=legs0[iseed];
	}
	for(int ic=0; ic<res->nx; ic++){
	    if(res->p[ic]){
		plot_points(toptab[ic], nseed, NULL, res->p[ic], NULL, NULL, xylog, NULL, legs,
			    title[ic], xlabel,ylabel, "%s",sidetab[ic]);
	    }
	}
	for(int iseed=0; iseed<nseed; iseed++){
	    free(legs0[iseed]);
	}
    }else{
	char *pathtag0[npath];
	const char *pathtag[npath];
	char prefix[4]="A: ";
	for(int ipath=0; ipath<npath; ipath++){
	    prefix[0]='A'+ipath;
	    pathtag0[ipath]=stradd(prefix, path[ipath], NULL);
	    pathtag[ipath]=pathtag0[ipath];
	}
	for(int ic=0; ic<res->nx; ic++){
	    if(res->p[ic]){
		plot_points(toptab[ic], npath, NULL, resm->p[ic], NULL, NULL, xylog, NULL, pathtag,
			    title[ic], xlabel,ylabel, "%s",sidetab[ic]);
	    }
	}

	for(int iseed=0; iseed<nseed && nseed>1; iseed++){
	    for(int ic=0; ic<res->nx; ic++){
		if(res->p[ic]){
		    dcell *tmp=dcellsub(res->p[ic], 0, 0, iseed, 1);
		    plot_points(toptab[ic], npath, NULL, tmp , NULL, NULL, xylog, NULL, pathtag,
				title[ic], xlabel,ylabel, "%s_%ld", sidetab[ic], seed[iseed]);
		    dcellfree(tmp);
		}
	    }
	}
    }
    draw_final(1);
    cellfree(res);
    cellfree(resm);
    /*
      writebin(upterr, "upterr");
      if(upterr && upterr->p[0]){
      for(int iseed=0; iseed<nseed; iseed++){

      plot_points("upterr", nseed, NULL, upterr, NULL, NULL, xylog, NULL, NULL, 
      "Uplink error", xlabel, "Error (rad)", "%d", iseed);
      }
      }*/
}
