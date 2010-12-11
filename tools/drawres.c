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

    DIR *dir=opendir(path);
    if(!dir){
	error("Unable to directory %s\n", path);
    }
    
    if(nseed==0){//do all seeds
	long nmax=20;
	seed=calloc(nmax, sizeof(long));
	struct dirent *dp;
	while((dp=readdir(dir))){
	    if(strncmp(dp->d_name, "Res_", 4)){
		continue;
	    }else{
		sscanf(dp->d_name, "Res_%ld.bin",&seed[nseed]);
		nseed++;
		if(nseed>nmax){
		    nmax*=2;
		    seed=realloc(seed, nmax*sizeof(long));
		}
	    }
	}
    	closedir(dir);
	seed=realloc(seed, nseed*sizeof(long));
    }
    for(int iseed=0; iseed<nseed; iseed++){
	info2("Processing seed %ld\n", seed[iseed]);
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "%s/Res_%ld.bin", path, seed[iseed]);
	dmat *res=dcellread(fn);
	if(res->p[3]){//split tomography.
	    dmat *tmp=dnew(
	}
    }
    free(seed);
}
