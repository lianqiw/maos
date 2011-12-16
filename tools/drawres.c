/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "aos.h"
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
    long nseed=arg->nseed;
    long *seed=arg->seeds;
    /*Find all valid path and seeds. */
    int jpath=0, mpath=0;
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
	if(!nseed){/*find available seeds */
	    long nmax=20;
	    seed=calloc(nmax, sizeof(long));
	    struct dirent *dp;
	    while((dp=readdir(dir))){
		if(strncmp(dp->d_name, "Res_", 4) || !check_suffix(dp->d_name, ".bin")){
		    continue;
		}else{
		    if(sscanf(dp->d_name, "Res_%ld.bin",&(seed[nseed]))==1){
			/*info2("found seed %ld\n", seed[nseed]); */
			nseed++;
			if(nseed>nmax){
			    nmax*=2;
			    seed=realloc(seed, nmax*sizeof(long));
			}
		    }
		}
	    }
	    closedir(dir);
	    if(nseed>0){
		seed=realloc(seed, nseed*sizeof(long));
	    }else{
		free(seed);/*try the next folder. */
	    }
	}
	/*record the path as valid */
	path[jpath]=path[ipath];
	info2("Found path: %s\n", path[jpath]);
	jpath++; mpath++;
    }
    npath=mpath;
    if(npath==0){
	info2("Nothing to display\n");
	return 1;
    }
    int jseed=0, mseed=0;
    /*Find seeds that are available in all folders. */
    for(int iseed=0; iseed<nseed; iseed++){
	int ipath;
	for(ipath=0; ipath<npath; ipath++){
	    char fn[PATH_MAX];
	    if(!path[ipath]){
		warning("path is NULL\n");
		break;
	    }
	    snprintf(fn,PATH_MAX,"%s/Res_%ld.bin", path[ipath], seed[iseed]);
	    if(!exist(fn)){
		break;
	    }
	}
	if(ipath==npath){
	    seed[jseed]=seed[iseed];
	    jseed++;
	    mseed++;
	}
    }
    nseed=mseed;
    if(nseed==0){
	info2("Nothing to display\n");
	return 1;
    }else{
	info2("Found seed:");
	for(int i=0; i<nseed; i++){
	    info2("%ld ", seed[i]);
	}
	info2("\n");
    }
    dcell *resolhi=dcellnew(npath,nseed);
    dcell *resollo=dcellnew(npath,nseed);
    dcell *reshi=dcellnew(npath,nseed);
    dcell *reslo=dcellnew(npath,nseed);
    dcell *reshim=dcellnew(npath,1);
    dcell *reslom=dcellnew(npath,1);
    dcell *resolhim=dcellnew(npath,1);
    dcell *resollom=dcellnew(npath,1);

    dcell *upterr=dcellnew(npath, nseed);
    dcell *uptcmd=dcellnew(npath, nseed);
    for(int iseed=0; iseed<nseed; iseed++){
	for(int ipath=0; ipath<npath; ipath++){
	    char fn[PATH_MAX];
	    snprintf(fn,PATH_MAX,"%s/Res_%ld.bin", path[ipath], seed[iseed]);
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
	    int ii=ipath+npath*iseed;
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
	    dadd(&reshim->p[ipath], 1, reshi->p[ii], 1);
	    dadd(&reslom->p[ipath], 1, reslo->p[ii], 1);
	    dadd(&resolhim->p[ipath], 1, resolhi->p[ii], 1);
	    dadd(&resollom->p[ipath], 1, resollo->p[ii], 1);
	    
	    snprintf(fn, PATH_MAX, "%s/Resuptcmd_%ld.bin", path[ipath], seed[iseed]);
	    dcell *upt=dcellread(fn);
	    tmp=dcell2m(upt);
	    uptcmd->p[ii]=dtrans(tmp);
	    dfree(tmp);
	    dcellfree(upt);
	    snprintf(fn, PATH_MAX, "%s/Resupterr_%ld.bin", path[ipath], seed[iseed]);
	    upt=dcellread(fn);
	    tmp=dcell2m(upt);
	    upterr->p[ii]=dtrans(tmp);
	    dfree(tmp);
	    dcellfree(upt);
	    info("upterr->p[%d]=%p. %ldx%ld\n", 
		 ii, upterr->p[ii], upterr->p[ii]->nx, upterr->p[ii]->ny);
	}
    }
    dcellcwpow(reshi, 0.5); dcellscale(reshi, 1e9);
    dcellcwpow(reslo, 0.5); dcellscale(reslo, 1e9);
    dcellcwpow(resolhi, 0.5); dcellscale(resolhi, 1e9);
    dcellcwpow(resollo, 0.5); dcellscale(resollo, 1e9);
    dcellscale(reshim, 1./nseed);
    dcellscale(reslom, 1./nseed);
    dcellscale(resolhim, 1./nseed);
    dcellscale(resollom, 1./nseed);
    dcellcwpow(reshim, 0.5); dcellscale(reshim, 1e9);
    dcellcwpow(reslom, 0.5); dcellscale(reslom, 1e9);
    dcellcwpow(resolhim, 0.5); dcellscale(resolhim, 1e9);
    dcellcwpow(resollom, 0.5); dcellscale(resollom, 1e9);


    if(npath==1){
	char *legs[nseed];
	for(int iseed=0; iseed<nseed; iseed++){
	    legs[iseed]=malloc(50*sizeof(char));
	    snprintf(legs[iseed], 50, "Seed %ld", seed[iseed]);
	}
	plot_points("Reshi", nseed, NULL, reshi, NULL, NULL, 0, NULL, legs,
		    "High order wavefront Error", "Steps","Error (nm)", "High");
	plot_points("Reslo", nseed, NULL, reslo, NULL, NULL, 0, NULL,legs,
		    "Low order wavefront Error", "Steps","Error (nm)", "Low");
	plot_points("ResOLhi", nseed, NULL, resolhi, NULL, NULL, 0, NULL, legs,
		    "High order open loop wavefront Error", "Steps","Error (nm)", "High");
	plot_points("ResOLlo", nseed, NULL, resollo, NULL, NULL, 0, NULL, legs,
		    "Low order open loop wavefront Error", "Steps","Error (nm)", "Low");
	for(int iseed=0; iseed<nseed; iseed++){
	    free(legs[iseed]);
	}
    }else{
	plot_points("Reshi", npath, NULL, reshim, NULL, NULL, 0, NULL, path,
		    "High order wavefront Error", "Steps","Error (nm)", "High");
	plot_points("Reslo", npath, NULL, reslom, NULL, NULL, 0, NULL, path,
		    "Low order wavefront Error", "Steps","Error (nm)", "Low");
	plot_points("ResOLhi", npath, NULL, resolhim, NULL, NULL, 0, NULL, path,
		    "High order open loop wavefront Error", "Steps","Error (nm)", "High");
	plot_points("ResOLlo", npath, NULL, resollom, NULL, NULL, 0, NULL, path,
		    "Low order open loop wavefront Error", "Steps","Error (nm)", "Low");

	for(int iseed=0; iseed<nseed; iseed++){
	    dcell *reshi_i=dcellsub(reshi, 0,0,iseed, 1);
	    dcell *reslo_i=dcellsub(reslo, 0,0,iseed, 1);
	    dcell *resolhi_i=dcellsub(resolhi, 0,0,iseed, 1);
	    dcell *resollo_i=dcellsub(resollo, 0,0,iseed, 1);
	
	    plot_points("Reshi", npath, NULL, reshi_i, NULL, NULL, 0, NULL, path,
			"High order wavefront Error", "Steps","Error (nm)", "High_%ld",seed[iseed]);
	    plot_points("Reslo", npath, NULL, reslo_i, NULL, NULL, 0, NULL, path,
			"Low order wavefront Error", "Steps","Error (nm)", "Low_%ld",seed[iseed]);
	    plot_points("ResOLhi", npath, NULL, resolhi_i, NULL, NULL, 0, NULL, path,
			"High order open loop wavefront Error", "Steps","Error (nm)", "High_%ld",seed[iseed]);
	    plot_points("ResOLlo", npath, NULL, resollo_i, NULL, NULL, 0, NULL, path,
			"Low order open loop wavefront Error", "Steps","Error (nm)", "Low_%ld",seed[iseed]);
	    dcellfree(reshi_i);dcellfree(reslo_i);dcellfree(resolhi_i);dcellfree(resollo_i);
	}
    }
    /*
    dcellwrite(upterr, "upterr");
    if(upterr && upterr->p[0]){
	for(int iseed=0; iseed<nseed; iseed++){

	    plot_points("upterr", nseed, NULL, upterr, NULL, NULL, 0, NULL, NULL, 
			"Uplink error", "Steps", "Error (rad)", "%d", iseed);
	}
	}*/
}
