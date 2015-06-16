/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "../lib/aos.h"
/**
   draw square map.
 */
static void draw_map(file_t *fp, int id){
    header_t header={0};
    read_header(&header, fp);
    char *name=mybasename(zfname(fp));
    char *suf=strstr(name, ".bin");
    if(!suf) suf=strstr(name, ".fits");
    if(suf) suf[0]='\0';
    if(iscell(&header.magic)){
	for(int ii=0; ii<header.nx*header.ny; ii++){
	    draw_map(fp, ii);
	}
	free(header.str);
    }else{
	dmat *in=dreaddata(fp, &header);
	map_t *data=d2map(in);
	drawmap("map", data, NULL, name,"x", "y", "%s:%d", name, id);
	mapfree(data);
	dfree(in);
    }
    free(name);
}
/** 
   A recursive opd drawing routine. The two files must contains matched information of loc grid and opd.
*/
static void draw_opd(file_t *fp1, file_t *fp2, int id){
    char *name=mybasename(zfname(fp2));
    char *suf=strstr(name, ".bin");
    if(!suf) suf=strstr(name, ".fits");
    if(suf) suf[0]='\0';
    header_t header1={0}, header2={0};
    read_header(&header1, fp1);
    read_header(&header2, fp2);
    if(iscell(&header1.magic) && iscell(&header2.magic)){ /*cells */
	if(header1.nx*header1.ny!=header2.nx*header2.ny){
	    error("cell arrays does have the same length.\n");
	}
	for(int ii=0; ii<header1.nx*header1.ny; ii++){
	    draw_opd(fp1, fp2, ii);
	}
	free(header1.str);
	free(header2.str);
    }else{
	loc_t *loc=locreaddata(fp1, &header1);
	dmat *opd=dreaddata(fp2, &header2);
	if(loc->nloc!=opd->nx){
	    error("we expect matching loc_t and a double vector.\n");
	}
	drawopd("opd", loc, opd->p, NULL, name,"x","y","%s:%d",name,id);
	dfree(opd);
	locfree(loc);
    }
    free(name);
}
/**
   A recursive loc drawing routine
*/
static void draw_loc(file_t *fp, int id){
    header_t header={0};
    read_header(&header, fp);
    char *name=mybasename(zfname(fp));
    char *suf=strstr(name, ".bin");
    if(!suf) suf=strstr(name, ".fits");
    if(suf) suf[0]='\0';
    if(iscell(&header.magic)){
   	for(int ii=0; ii<header.nx*header.ny; ii++){
	    draw_loc(fp, ii);
	}
	free(header.str);
    }else{
	loc_t *loc=locreaddata(fp, &header);
	if(loc->nloc>100000){/*if too many points, we draw it. */
	    drawloc("loc",loc,NULL,zfname(fp),"x","y","%s",zfname(fp));
	}else{/*we plot individual points. */
	    plot_points("loc",1, &loc, NULL, NULL,NULL,NULL,0,NULL,NULL,name,"x","y","%s:%d",name,id);
	}
	locfree(loc);
    }
    free(name);
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
    if(argc<2){
	usage();
	return 0;
    }
    /*use the parent pid so same bash session has the same drawdaemon. */
    DRAW_ID=getppid();
    DRAW_DIRECT=1;
    /*launch scheduler if it is not already running. */
    if(!strcmp(argv[1],"loc")){/*draw coordinate grid */
	if(argc!=3){
	    error("Invalid number of input\n");
	}
	file_t *fp=zfopen(argv[2],"r");
	draw_loc(fp, 0);
	zfclose(fp);
    }else if(!strcmp(argv[1],"opd")){/*draw OPD with coordinate */
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
    draw_final(1);
}

