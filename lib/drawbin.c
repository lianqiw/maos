#include "aos.h"

#if USE_DAEMON==1
/**
   draw square map.
 */
static void draw_map(file_t *fp, int id){
    char *header=NULL;
    uint32_t magic=read_magic(fp, &header);
    uint64_t nx,ny;
    char *name=mybasename(fp->fn);
    char *suf=strstr(name, ".bin");
    suf[0]='\0';
    switch(magic){
    case MC_DBL:
    case MCC_ANY:
	zfreadlarr(fp, 2, &nx, &ny);
	for(int ii=0; ii<nx*ny; ii++){
	    draw_map(fp, ii);
	}
    	break;
    case M_DBL:
	{
	    map_t *data=sqmapreaddata(fp, magic, header);
	    drawmap("map", data, name, "x", "y", "%s:%d", name, id);
	    sqmapfree(data);
	}
	break;
    default:
	error("magic=%x, Not implemented yet\n", magic);
    }
    free(header);
    free(name);
}
/** 
   A recursive opd drawing routine. The two files must contains matched information of loc grid and opd.
*/
static void draw_opd(file_t *fp1, file_t *fp2, int id){
    uint32_t magic1, magic2;
    uint64_t nx1,ny1,nx2,ny2;
    char *name=mybasename(fp2->fn);
    char *suf=strstr(name, ".bin");
    suf[0]='\0';
    char *header1=NULL, *header2=NULL;
    magic1=read_magic(fp1, &header1);
    magic2=read_magic(fp2, &header2);
    if((magic1==MC_DBL || magic1==MCC_ANY) && (magic2==MC_DBL || magic2==MCC_ANY)){
	//cells
	zfreadlarr(fp1, 2, &nx1, &ny1);
	zfreadlarr(fp2, 2, &nx2, &ny2);
	if(nx1*ny1!=nx2*ny2){
	    error("cell arrays does have the same length.\n");
	}
	for(int ii=0; ii<nx1*ny1; ii++){
	    draw_opd(fp1, fp2, ii);
	}
    }else if(magic1==M_DBL && magic2==M_DBL){//single
	loc_t *loc=locreaddata(fp1, magic1, header1);
	dmat *opd=dreaddata(fp2, magic2);
	if(loc->nloc!=opd->nx){
	    error("we expect matching loc_t and a double vector.\n");
	}
	drawopd("opd", loc, opd->p, name,"x","y","%s:%d",name,id);
	dfree(opd);
	locfree(loc);
    }else{
	error("Invalid input\n");
    }
    free(header1);
    free(header2);
    free(name);
}
/**
   A recursive loc drawing routine
*/
static void draw_loc(file_t *fp, int id){
    char *header=NULL;
    uint32_t magic=read_magic(fp, &header);
    uint64_t nx,ny;
    char *name=mybasename(fp->fn);
    char *suf=strstr(name, ".bin");
    suf[0]='\0';
    switch(magic){
    case MC_DBL:
    case MCC_ANY:
	zfreadlarr(fp, 2, &nx, &ny);
   	for(int ii=0; ii<nx*ny; ii++){
	    draw_loc(fp, ii);
	}
    	break;
    case M_DBL:
	{
	    loc_t *loc=locreaddata(fp, magic, header);
	    if(loc->nloc>100000){//if too many points, we draw it.
		drawloc("loc",loc,fp->fn,"x","y","%s",fp->fn);
	    }else{//we plot individual points.
		plot_coord("loc",loc->nloc,loc->locx,loc->locy,
			   NULL,NULL,0,NULL,name,"x","y","%s:%d",name,id);
	    }
	    locfree(loc);
	}
	break;
    default:
	error("magic=%x, nx=%lu, ny=%lu. Not implemented yet\n", magic, nx, ny);
    }
    free(name);
    free(header);
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
    //launch scheduler if it is not already running.
    scheduler_launch();
    if(!strcmp(argv[1],"loc")){//draw coordinate grid
	if(argc!=3){
	    error("Invalid number of input\n");
	}
	file_t *fp=zfopen(argv[2],"r");
	draw_loc(fp, 0);
	zfclose(fp);
    }else if(!strcmp(argv[1],"opd")){//draw OPD with coordinate
	if(argc==3){
	    file_t *fp1=zfopen(argv[2],"r");
	    draw_map(fp1, 0);
	    zfclose(fp1);
	}else if(argc==4){
	    file_t *fp1=zfopen(argv[2],"r");
	    file_t *fp2=zfopen(argv[3],"r");
	    draw_opd(fp1, fp2, 0);
	    zfclose(fp1);
	    zfclose(fp2);
	}else{
	    error("Invalid number of input\n");
	}
    }else{
	error("Invalid arguments\n");
    }
    exit_success=1;
}

#endif
