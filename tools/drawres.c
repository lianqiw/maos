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
    return arg;
}
/**
   The main.
*/
int main(int argc, char *argv[]){
    if(mystrcmp(argv[0], "scheduler")==0){//launch the scheduler.
	scheduler();
	exit(0);
    }
    ARG_T *arg=parse_args(argc, argv);
    //use the parent pid so same bash session has the same drawdaemon.
    DRAW_ID=getppid();
    char **path;
    int npath;
    if(argc>=2){
	npath=argc-1;
	path=calloc(npath, sizeof(char*));
	for(int ipath=0; ipath<npath; ipath++){
	    path[ipath]=argv[ipath+1];
	}
    }else{
	npath=1;
	path=calloc(npath, sizeof(char*));
	path[0]=mygetcwd();
    }
    long nseed=arg->nseed;
    long *seed=arg->seeds;
    if(arg->nseed==0){//do all seeds
	DIR *dir=opendir(path[0]);
	if(!dir){
	    error("Unable to directory %s\n", path[0]);
	}
	long nmax=20;
	seed=calloc(nmax, sizeof(long));
	struct dirent *dp;
	while((dp=readdir(dir))){
	    if(strncmp(dp->d_name, "Res_", 4) || !check_suffix(dp->d_name, ".bin")){
		continue;
	    }else{
		if(sscanf(dp->d_name, "Res_%ld.bin",&(seed[nseed]))==1){
		    info2("Processing seed %ld, %s\n", seed[nseed], dp->d_name);
		    nseed++;
		    if(nseed>nmax){
			nmax*=2;
			seed=realloc(seed, nmax*sizeof(long));
		    }
		}
	    }
	}
    	closedir(dir);
	seed=realloc(seed, nseed*sizeof(long));
    }
    info2("nseed=%ld\n", nseed);
    if(nseed==0){
	info2("Nothing to display\n");
	return 1;
    }
    dcell *resolhi=dcellnew(npath,nseed);
    dcell *resollo=dcellnew(npath,nseed);
    dcell *reshi=dcellnew(npath,nseed);
    dcell *reslo=dcellnew(npath,nseed);
    int skip[nseed];
    memset(skip, 0, sizeof(int)*nseed);
    for(int iseed=0; iseed<nseed; iseed++){
	for(int ipath=0; ipath<npath; ipath++){
	    char fn[PATH_MAX];
	    snprintf(fn,PATH_MAX,"%s/Res_%ld.bin", path[ipath], seed[iseed]);
	    dcell *res;
	    if(exist(fn)){
		res=dcellread("%s",fn);
	    }else{
		skip[iseed]=1;
		break;
	    }
	    int ind=0;
	    int indlo=0;
	    if(res->p[3]){//split tomography.
		ind=3;
		indlo=2;
	    }else{
		ind=2;
		indlo=1;
	    }
	    int index=ipath+npath*iseed;
	    dmat *tmp;
	    tmp=dsub(res->p[ind], 0, 1, 0, 0);
	    reshi->p[index]=dtrans(tmp);
	    dfree(tmp);
	    tmp=dsub(res->p[ind], indlo, 1, 0, 0);
	    reslo->p[index]=dtrans(tmp);
	    dfree(tmp);

	    tmp=dsub(res->p[0], 0, 1, 0, 0);
	    resolhi->p[index]=dtrans(tmp);
	    dfree(tmp);
	    tmp=dsub(res->p[0], 1, 1, 0, 0);
	    resollo->p[index]=dtrans(tmp);
	    dfree(tmp);
	}
    }
    dcellcwpow(reshi, 0.5); dcellscale(reshi, 1e9);
    dcellcwpow(reslo, 0.5); dcellscale(reslo, 1e9);
    dcellcwpow(resolhi, 0.5); dcellscale(resolhi, 1e9);
    dcellcwpow(resollo, 0.5); dcellscale(resollo, 1e9);
    if(npath==1){
	plot_points("Res", nseed, NULL, reshi, NULL, NULL, 0, NULL, 
		    "High order wavefront Error", "Steps","Error (nm)", "High");
	plot_points("Res", nseed, NULL, reslo, NULL, NULL, 0, NULL,
		    "Low order wavefront Error", "Steps","Error (nm)", "Low");
	plot_points("ResOL", nseed, NULL, resolhi, NULL, NULL, 0, NULL, 
		    "High order open loop wavefront Error", "Steps","Error (nm)", "High");
	plot_points("ResOL", nseed, NULL, resollo, NULL, NULL, 0, NULL, 
		    "Low order open loop wavefront Error", "Steps","Error (nm)", "Low");
    }else{
	for(int iseed=0; iseed<nseed; iseed++){
	    if(skip[iseed]) continue;
	    dcell *reshi_i=dcellsub(reshi, 0,0,iseed, 1);
	    dcell *reslo_i=dcellsub(reslo, 0,0,iseed, 1);
	    dcell *resolhi_i=dcellsub(resolhi, 0,0,iseed, 1);
	    dcell *resollo_i=dcellsub(resollo, 0,0,iseed, 1);
	    plot_points("Res", npath, NULL, reshi_i, NULL, NULL, 0, NULL, 
			"High order wavefront Error", "Steps","Error (nm)", "High_%ld",seed[iseed]);
	    plot_points("Res", npath, NULL, reslo_i, NULL, NULL, 0, NULL,
			"Low order wavefront Error", "Steps","Error (nm)", "Low_%ld",seed[iseed]);
	    plot_points("ResOL", npath, NULL, resolhi_i, NULL, NULL, 0, NULL, 
			"High order open loop wavefront Error", "Steps","Error (nm)", "High_%ld",seed[iseed]);
	    plot_points("ResOL", npath, NULL, resollo_i, NULL, NULL, 0, NULL, 
			"Low order open loop wavefront Error", "Steps","Error (nm)", "Low_%ld",seed[iseed]);
	    dcellfree(reshi_i);dcellfree(reslo_i);dcellfree(resolhi_i);dcellfree(resollo_i);
	}
    }
}
