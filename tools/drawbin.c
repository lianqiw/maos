/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "aos.h"
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
    if(iscell(magic)){
	zfreadlarr(fp, 2, &nx, &ny);
	for(int ii=0; ii<nx*ny; ii++){
	    draw_map(fp, ii);
	}
    }else{
	map_t *data=sqmapreaddata(fp, magic, header);
	drawmap("map", data, name, "x", "y", "%s:%d", name, id);
	sqmapfree(data);
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
    if(iscell(magic1) && iscell(magic2)){ //cells
	zfreadlarr(fp1, 2, &nx1, &ny1);
	zfreadlarr(fp2, 2, &nx2, &ny2);
	if(nx1*ny1!=nx2*ny2){
	    error("cell arrays does have the same length.\n");
	}
	for(int ii=0; ii<nx1*ny1; ii++){
	    draw_opd(fp1, fp2, ii);
	}
    }else{
	loc_t *loc=locreaddata(fp1, magic1, header1);
	dmat *opd=dreaddata(fp2, magic2);
	if(loc->nloc!=opd->nx){
	    error("we expect matching loc_t and a double vector.\n");
	}
	drawopd("opd", loc, opd->p, name,"x","y","%s:%d",name,id);
	dfree(opd);
	locfree(loc);
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
    if(iscell(magic)){
	zfreadlarr(fp, 2, &nx, &ny);
   	for(int ii=0; ii<nx*ny; ii++){
	    draw_loc(fp, ii);
	}
    }else{
	loc_t *loc=locreaddata(fp, magic, header);
	if(loc->nloc>100000){//if too many points, we draw it.
	    drawloc("loc",loc,fp->fn,"x","y","%s",fp->fn);
	}else{//we plot individual points.
	    plot_points("loc",1, &loc, NULL, NULL,NULL,0,NULL,name,"x","y","%s:%d",name,id);
	}
	locfree(loc);
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

