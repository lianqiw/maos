#include <dirent.h>
#include "aos.h"
/**
   Plot MAOS Results
*/
/**
   The main.
*/
int main(int argc, char *argv[]){
    if(mystrcmp(argv[0], "scheduler")==0){//launch the scheduler.
	scheduler();
	exit(0);
    }
    //use the parent pid so same bash session has the same drawdaemon.
    DRAW_ID=getppid();
    long nseed=0;
    long *seed;
    if(argc>1){
	nseed=1;
	seed=calloc(1, sizeof(long));
	seed[0]=strtol(argv[1], NULL, 10);
    }
    char *path;
    if(argc>2){
	path=argv[2];
    }else{
	path=mygetcwd();
    }
    if(nseed==0){//do all seeds
	DIR *dir=opendir(path);
	if(!dir){
	    error("Unable to directory %s\n", path);
	}
	long nmax=20;
	seed=calloc(nmax, sizeof(long));
	struct dirent *dp;
	while((dp=readdir(dir))){
	    if(strncmp(dp->d_name, "Res_", 4) || !check_suffix(dp->d_name, ".bin")){
		continue;
	    }else{
		if(sscanf(dp->d_name, "Res_%ld.bin",&seed[nseed])==1){
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
    dcell *resolhi=dcellnew(nseed,1);
    dcell *resollo=dcellnew(nseed,1);
    dcell *reshi=dcellnew(nseed,1);
    dcell *reslo=dcellnew(nseed,1);
    for(int iseed=0; iseed<nseed; iseed++){
	dcell *res=dcellread("%s/Res_%ld.bin", path, seed[iseed]);
	int ind=0;
	int indlo=0;
	if(res->p[3]){//split tomography.
	    ind=3;
	    indlo=2;
	}else{
	    ind=2;
	    indlo=1;
	}
	dmat *tmp;
	tmp=dsub(res->p[ind], 0, 1, 0, 0);
	reshi->p[iseed]=dtrans(tmp);
	dfree(tmp);
	tmp=dsub(res->p[ind], indlo, 1, 0, 0);
	reslo->p[iseed]=dtrans(tmp);
	dfree(tmp);

	tmp=dsub(res->p[0], 0, 1, 0, 0);
	resolhi->p[iseed]=dtrans(tmp);
	dfree(tmp);
	tmp=dsub(res->p[0], 1, 1, 0, 0);
	resollo->p[iseed]=dtrans(tmp);
	dfree(tmp);
    }

    dcellcwpow(reshi, 0.5);
    dcellscale(reshi, 1e9);
    dcellcwpow(reslo, 0.5);
    dcellscale(reslo, 1e9);
    dcellcwpow(resolhi, 0.5);
    dcellscale(resolhi, 1e9);
    dcellcwpow(resollo, 0.5);
    dcellscale(resollo, 1e9);
    plot_points("Res", nseed, NULL, reshi, NULL, NULL, 0, NULL, 
		"High order wavefront Error", "Steps","Error (nm)", "High");
    plot_points("Res", nseed, NULL, reslo, NULL, NULL, 0, NULL,
		"Low order wavefront Error", "Steps","Error (nm)", "Low");
    plot_points("ResOL", nseed, NULL, resolhi, NULL, NULL, 0, NULL, 
		"High order open loop wavefront Error", "Steps","Error (nm)", "High");
    plot_points("ResOL", nseed, NULL, resollo, NULL, NULL, 0, NULL, 
		"Low order open loop wavefront Error", "Steps","Error (nm)", "Low");
    free(seed);
}
