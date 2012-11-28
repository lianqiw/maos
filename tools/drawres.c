/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <unistd.h>
#include "../lib/aos.h"
static void usage(){
    fprintf(stderr,"Usage:\n"
	    "drawres [-s 1] folder1 folder2 ...\n");
}
/**
   Plot MAOS Results
*/
typedef struct ARG_T{
    int iarg;
    int nseed;
    long *seeds;
}ARG_T;
static ARG_T * parse_args(int argc, char **argv){
    ARG_T *arg=calloc(1, sizeof(ARG_T));
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
	    arg->seeds=realloc(arg->seeds, sizeof(long)*arg->nseed);
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
    DRAW_ID=getppid();/*variables in draw.c */
    DRAW_DIRECT=1;/*launch drawdaemon directly, without going through server. */
    char **path;
    int npath;

    if(arg->iarg<argc){
	npath=argc-arg->iarg;
	path=calloc(npath, sizeof(char*));
	for(int ipath=0; ipath<npath; ipath++){
	    path[ipath]=argv[ipath+arg->iarg];
	}
    }else{
	npath=1;
	path=calloc(npath, sizeof(char*));
	path[0]=mygetcwd();
    }
    long nseed=0;
    long *seed=NULL;
    long *seed2=NULL;
    /*Find all valid path and seeds. */
    int jpath=0, mpath=0;
    int restype=-1;
    for(int ipath=0; ipath<npath; ipath++){
	if(!isdir(path[ipath])){
	    continue;
	}
	/*info("Try path %s\n", path[ipath]); */
	DIR *dir=opendir(path[ipath]);
	if(!dir){
	    warning("Unable to read directory %s\n", path[0]);
	    continue;
	}
	struct dirent *dp;
	int nseedfound=0;
	while((dp=readdir(dir))){
	    if(strncmp(dp->d_name, "Res", 3) || !check_suffix(dp->d_name, ".bin")){
		continue;
	    }else{
		long seedi, seedi2;
		int found=0;
		if(sscanf(dp->d_name, "Res_%ld.bin",&(seedi))==1){
		    nseedfound++;
		    if(restype!=-1 && restype!=1){
			error("Only maos or skyc results can display together\n");
		    }
		    restype=1;
		    if(arg->nseed){
			int wanted=0;
			for(int is=0; is<arg->nseed; is++){
			    if(arg->seeds[is]==seedi){
				wanted=1;
				break;
			    }
			}
			if(!wanted) {
			    warning2("Skip seed %ld\n", seedi);
			    continue;
			}
		    }
		    for(int is=0; is<nseed; is++){
			if(seed[is]==seedi){
			    found=1;
			}
		    }
		    if(!found){
			info2("Found unique seed %ld\n", seedi);
			nseed++;
			seed=realloc(seed, sizeof(long)*nseed);
			seed[nseed-1]=seedi;
		    }
		}else if(sscanf(dp->d_name, "Res%ld_%ld.bin",&seedi, &seedi2)==2){
		    nseedfound++;/*sky coverage*/
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
			info("Found unique seed %ld %ld\n", seedi, seedi2);
			nseed++;
			seed=realloc(seed, sizeof(long)*nseed);
			seed[nseed-1]=seedi;
			seed2=realloc(seed2, sizeof(long)*nseed);
			seed2[nseed-1]=seedi2;
		    }
		}
	    }
	}
	closedir(dir);
	if(nseedfound){ /*record the path as valid */
	    path[jpath]=path[ipath];
	    info2("Found path: %s\n", path[jpath]);
	    jpath++; mpath++;
	}
    }
    npath=mpath;
    if(npath==0){
	info2("Nothing to display\n");
	return 1;
    }
  
    info2("Found seed:");
    for(int i=0; i<nseed; i++){
	info2("%ld ", seed[i]);
    }
    info2("\n");

    dcell *resolhi=dcellnew(npath,nseed);
    dcell *resollo=dcellnew(npath,nseed);
    dcell *reshi=dcellnew(npath,nseed);
    dcell *reslo=dcellnew(npath,nseed);
    dcell *reshim=dcellnew(npath,1);
    dcell *reslom=dcellnew(npath,1);
    dcell *resolhim=dcellnew(npath,1);
    dcell *resollom=dcellnew(npath,1);
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
	    seedcount++;
	    int ii=ipath+npath*iseed;
	    if(restype==1){
		dcell *res;
		res=dcellread("%s",fn);
		int ind=0;
		int indlo=0;
		int indhi=0;
		if(res->p[3] && res->p[3]->nx>0){/*split tomography. */
		    ind=3;
		    indlo=2;
		    indhi=0;
		}else{
		    ind=2;
		    indlo=1;/*tt */
		    indhi=2;/*pttr */
		}
		dmat *tmp;
		tmp=dsub(res->p[ind], indhi, 1, 0, 0);
		reshi->p[ii]=dtrans(tmp);
		dfree(tmp);
		tmp=dsub(res->p[ind], indlo, 1, 0, 0);
		reslo->p[ii]=dtrans(tmp);
		dfree(tmp);

		tmp=dsub(res->p[0], 2, 1, 0, 0);
		resolhi->p[ii]=dtrans(tmp);
		dfree(tmp);
		tmp=dsub(res->p[0], 1, 1, 0, 0);
		resollo->p[ii]=dtrans(tmp);
		dfree(tmp);
		dcellfree(res);
	    }else if(restype==2){
		snprintf(fn,PATH_MAX,"%s/Res%ld_%ld.bin", path[ipath], seed[iseed], seed2[iseed]);
		dmat *res0=dread("%s", fn);
		dmat *res=dtrans(res0); dfree(res0);
		dmat *tmp;
		if(!ysky){
		    int nsky=res->nx;
		    ysky=dnew(nsky,1);
		    for(int i=0; i<nsky; i++){
			ysky->p[i]=pow((double)i/(double)(nsky-1)*1e-9,2);
		    }
		}
		
		tmp=dsub(res, 0, 0, 0, 1);
		dsort(tmp, 1);
		reshi->p[ii]=dcat(tmp, ysky, 2); 
		dfree(tmp);

		tmp=dsub(res, 0, 0, 2, 1);
		dsort(tmp, 1);
		reslo->p[ii]=dcat(tmp, ysky, 2); 
		dfree(tmp);
	    }else{
		error("Invalid restype=%d\n", restype);
	    }

	    dadd(&reshim->p[ipath], 1, reshi->p[ii], 1);
	    dadd(&reslom->p[ipath], 1, reslo->p[ii], 1);
	    if(restype==1){
		dadd(&resolhim->p[ipath], 1, resolhi->p[ii], 1);
		dadd(&resollom->p[ipath], 1, resollo->p[ii], 1);
	    }
	    /*
	      snprintf(fn, PATH_MAX, "%s/Resuptcmd_%ld.bin", path[ipath], seed[iseed]);
	      dcell *upt=dcellread("%s",fn);
	    tmp=dcell2m(upt);
	    uptcmd->p[ii]=dtrans(tmp);
	    dfree(tmp);
	    dcellfree(upt);
	    snprintf(fn, PATH_MAX, "%s/Resupterr_%ld.bin", path[ipath], seed[iseed]);
	    upt=dcellread("%s",fn);
	    tmp=dcell2m(upt);
	    upterr->p[ii]=dtrans(tmp);
	    dfree(tmp);
	    dcellfree(upt);*/
	}
	if(seedcount>0){
	    dscale(reshim->p[ipath], 1./seedcount);
	    dscale(reslom->p[ipath], 1./seedcount);
	    if(restype==1){
		dscale(resolhim->p[ipath], 1./seedcount);
		dscale(resollom->p[ipath], 1./seedcount);
	    }
	}
    }
    dcellcwpow(reshi, 0.5); dcellscale(reshi, 1e9);
    dcellcwpow(reshim, 0.5); dcellscale(reshim, 1e9);
    dcellcwpow(reslo, 0.5); dcellscale(reslo, 1e9);
    dcellcwpow(reslom, 0.5); dcellscale(reslom, 1e9);
    if(restype==1){
	dcellcwpow(resolhi, 0.5); dcellscale(resolhi, 1e9);
	dcellcwpow(resollo, 0.5); dcellscale(resollo, 1e9);

	dcellcwpow(resolhim, 0.5); dcellscale(resolhim, 1e9);
	dcellcwpow(resollom, 0.5); dcellscale(resollom, 1e9);
    }
 
    if(npath==1){
	char *legs[nseed];
	for(int iseed=0; iseed<nseed; iseed++){
	    legs[iseed]=malloc(50*sizeof(char));
	    snprintf(legs[iseed], 50, "Seed %ld", seed[iseed]);
	}
	if(restype==1){
	plot_points("Reshi", nseed, NULL, reshi, NULL, NULL, xylog, 0, NULL, legs,
		    "High order Wavefront Error", xlabel,ylabel, "High");
	plot_points("Reslo", nseed, NULL, reslo, NULL, NULL, xylog, 0, NULL,legs,
		    "Low order Wavefront Error", xlabel,ylabel, "Low");
	plot_points("ResOLhi", nseed, NULL, resolhi, NULL, NULL, xylog, 0, NULL, legs,
		    "High order Openloop Wavefront Error", xlabel,ylabel, "High");
	plot_points("ResOLlo", nseed, NULL, resollo, NULL, NULL, xylog, 0, NULL, legs,
		    "Low order Openloop Wavefront Error", xlabel,ylabel, "Low");
	}else{
	    plot_points("Tot", nseed, NULL, reshi, NULL, NULL, xylog, 0, NULL, legs,
			"Total OIWFS Mode Wavefront Error", xlabel,ylabel, "All");
	    plot_points("Atm TT", nseed, NULL, reslo, NULL, NULL, xylog, 0, NULL,legs,
			"ATM T/T Wavefront Error", xlabel,ylabel, "TT");
	}
	for(int iseed=0; iseed<nseed; iseed++){
	    free(legs[iseed]);
	}
    }else{
	if(restype==1){
	    plot_points("Reshi", npath, NULL, reshim, NULL, NULL, xylog, 0, NULL, path,
			"High order Wavefront Error", xlabel,ylabel, "High");
	    plot_points("Reslo", npath, NULL, reslom, NULL, NULL, xylog, 0, NULL, path,
			"Low order Wavefront Error", xlabel,ylabel, "Low");
	    plot_points("ResOLhi", npath, NULL, resolhim, NULL, NULL, xylog, 0, NULL, path,
			"High order Openloop Wavefront Error", xlabel,ylabel, "High");
	    plot_points("ResOLlo", npath, NULL, resollom, NULL, NULL, xylog, 0, NULL, path,
			"Low order Openloop Wavefront Error", xlabel,ylabel, "Low");
	}else{
	    plot_points("Tot", npath, NULL, reshim, NULL, NULL, xylog, 0, NULL, path,
			"Total OIWFS Mode Wavefront Error", xlabel,ylabel, "All");
	    plot_points("Atm TT", npath, NULL, reslom, NULL, NULL, xylog, 0, NULL,path,
			"ATM T/T Wavefront Error", xlabel,ylabel, "TT");
	}
	for(int iseed=0; iseed<nseed; iseed++){
	    dcell *reshi_i=dcellsub(reshi, 0,0,iseed, 1);
	    dcell *reslo_i=dcellsub(reslo, 0,0,iseed, 1);
	    dcell *resolhi_i=dcellsub(resolhi, 0,0,iseed, 1);
	    dcell *resollo_i=dcellsub(resollo, 0,0,iseed, 1);
	    if(restype==1){
	    plot_points("Reshi", npath, NULL, reshi_i, NULL, NULL, xylog, 0, NULL, path,
			"High order Wavefront Error", xlabel,ylabel, "High_%ld",seed[iseed]);
	    plot_points("Reslo", npath, NULL, reslo_i, NULL, NULL, xylog, 0, NULL, path,
			"Low order Wavefront Error", xlabel,ylabel, "Low_%ld",seed[iseed]);
	    plot_points("ResOLhi", npath, NULL, resolhi_i, NULL, NULL, xylog, 0, NULL, path,
			"High order Openloop Wavefront Error", xlabel,ylabel, "High_%ld",seed[iseed]);
	    plot_points("ResOLlo", npath, NULL, resollo_i, NULL, NULL, xylog, 0, NULL, path,
			"Low order Openloop Wavefront Error", xlabel,ylabel, "Low_%ld",seed[iseed]);
	    }else{
		plot_points("Tot", npath, NULL, reshi_i, NULL, NULL, xylog, 0, NULL, path,
			    "Total OIWFS Mode Wavefront Error", xlabel,ylabel, "All_%ld",seed[iseed]);
		plot_points("Atm TT", npath, NULL, reslo_i, NULL, NULL, xylog, 0, NULL,path,
			    "ATM T/T Wavefront Error", xlabel,ylabel, "TT_%ld",seed[iseed]);

	    }
	    dcellfree(reshi_i);dcellfree(reslo_i);dcellfree(resolhi_i);dcellfree(resollo_i);
	}
    }
    /*
    dcellwrite(upterr, "upterr");
    if(upterr && upterr->p[0]){
	for(int iseed=0; iseed<nseed; iseed++){

	    plot_points("upterr", nseed, NULL, upterr, NULL, NULL, xylog, 0, NULL, NULL, 
			"Uplink error", xlabel, "Error (rad)", "%d", iseed);
	}
	}*/
}
