#include "aos.h"

#if USE_DAEMON==1

/**
   A recursive opd drawing routine. The two files must contains matched information of loc grid and opd.
*/
static void draw_opd(file_t *fp1, file_t *fp2, int id){
    uint32_t magic1, magic2;
    zfread(&magic1, sizeof(uint32_t), 1, fp1);
    zfread(&magic2, sizeof(uint32_t), 1, fp2);
    uint64_t nx1,ny1,nx2,ny2;
    zfreadlarr(fp1, 2, &nx1, &ny1);
    zfreadlarr(fp2, 2, &nx2, &ny2);
    char *name=strdup(fp2->fn);
    char *dot=index(name,'.');
    if(dot) *dot='\0';
    if(magic1==MC_DBL && magic2==MC_DBL){
	if(nx1*ny1!=nx2*ny2){
	    error("cell arrays does have the same length.\n");
	}
	for(int ii=0; ii<nx1*ny1; ii++){
	    draw_opd(fp1, fp2, ii);
	}
    }else if(magic1==M_DBL && magic2==M_DBL){
	if(nx1!=nx2 || ny1!=2 || ny2!=1){
	    info("nx1=%lu, ny1=%lu, nx2=%lu, ny2=%lu\n",nx1,ny1,nx2,ny2);
	    error("we expect matching loc_t and a double vector.\n");
	}
	loc_t *loc=locreaddata2(fp1, nx1, ny1);
	double *opd=malloc(sizeof(double)*nx2);
	zfread(opd, sizeof(double), nx2, fp2);
	drawopd("opd", loc, opd, name,"x","y","%s:%d",name,id);
    }else{
	error("Invalid input\n");
    }
}
/**
   A recursive loc drawing routine
*/
static void draw_loc(file_t *fp, int id){
    uint32_t magic;
    zfread(&magic, sizeof(uint32_t), 1, fp);
    uint64_t nx,ny;
    zfreadlarr(fp, 2, &nx, &ny);
    char *name=strdup(fp->fn);
    char *dot=index(name,'.');
    if(dot) *dot='\0';
    switch(magic){
    case MC_DBL:
	for(int ii=0; ii<nx*ny; ii++){
	    draw_loc(fp, ii);
	}
    	break;
    case M_DBL:
	{
	    if(ny==2 && nx>2){//plot of coordinates
		loc_t *loc=locreaddata2(fp, nx, ny);
		info2("drawloc");
		if(loc->nloc>100000){//if too many points, we draw it.
		    drawloc(fp->fn,loc,fp->fn,"x","y","%s",fp->fn);
		}else{//we plot individual points.
		    plot_coord("loc",loc->nloc,loc->locx,loc->locy,
			       NULL,NULL,0,NULL,name,"x","y","%s:%d",name,id);
		}
		locfree(loc);
	    }else{
		warning("skip %d\n",id);
	    }
	}
	break;
    default:
	error("Not implemented yet\n");
    }
}
static void usage(){
	info2("Usage:\n"
"drawbin loc ploc.bin\n"
"drawbin opd powfs0_loc.bin powfs0_amp.bin\n"
);
}
/**
   The main.
*/
int main(int argc, char *argv[]){
    if(mystrcmp(argv[0], "scheduler")==0){//launch the scheduler.
	scheduler();
	exit(0);
    }
    if(argc<2){
	usage();
	return 0;
    }
    //use the parent pid so same bash session has the same drawdaemon.
    DRAW_ID=getppid();
    scheduler_launch();
    //launch scheduler if it is not already running.
    if(!strcmp(argv[1],"loc")){
	if(argc!=3){
	    error("Invalid number of input\n");
	}
	file_t *fp=zfopen(argv[2],"r");
	draw_loc(fp, 0);
	
    }else if(!strcmp(argv[1],"opd")){
	if(argc!=4){
	    error("Invalid number of input\n");
	}
	file_t *fp1=zfopen(argv[2],"r");
	file_t *fp2=zfopen(argv[3],"r");
	draw_opd(fp1, fp2, 0);
    }else{
	error("Invalid arguments\n");
    }
}

#endif
