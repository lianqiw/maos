/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#include <unistd.h>
#include "../lib/aos.h"

static void usage(){
	info("Usage:\n"
		"drawbin loc.bin\n"
		"drawbin loc.bin opd.bin\n"
		"drawbin map.bin\n"
	);
}
/**
   The main.
*/
int main(int argc, char* argv[]){
	if(argc<2){
		usage();
		return 0;
	}
	/*use the parent pid so same bash session has the same drawdaemon. */
	draw_id=getsid(0)+1e6;
	draw_direct=1;
	draw_single=-1;//disable draw_single support;
	/*launch scheduler if it is not already running. */
	dcell *arg1=0, *arg2=0;
	if(argc>1){
		arg1=dcellread("%s", argv[1]);
		if(!arg1){
			error("Unable to read %s\n", argv[1]);
		}else{
			info("%s is %ldx%ld\n", argv[1], arg1->nx, arg1->ny);
		}
	}
	if(argc>2){
		arg2=dcellread("%s",argv[2]);
		if(!arg2){
			error("Unable to read %s\n", argv[2]);
		}else{
			info("%s is %ldx%ld\n", argv[2], arg2->nx, arg2->ny);
		}
	}
	int mx=0;
	int my=0;
	if(argc>3){
		mx=strtol(argv[3], NULL, 10);
	}
	if(argc>4){
		my=strtol(argv[4], NULL, 10);
	}
	char *title1=argv[argc>2?2:1];
	char *title2=NULL;
	if((title2=strrchr(title1, '/'))){
		title1=title2+1;
	}
	if((title2=strrchr(title1, '.'))){
		title2[0]='\0';
	}
	
	long nx=MAX(NX(arg1), arg2?NX(arg2):0); nx=MIN(nx, mx+20);
	long ny=MAX(NY(arg1), arg2?NY(arg2):0); ny=MIN(ny, my+20);
	int id=0;
	loc_t *loc_save=0;
	dmat* p1_save=0;
	for(long iy=my; iy<ny; iy++){
		for(long ix=mx; ix<nx; ix++){
			dmat *p1=(dmat*)PR(arg1, ix, iy);
			dmat *p2=arg2?(dmat*)PR(arg2, ix, iy):NULL;
			if(!p1 || (!p2 && (ix>=NX(arg1) || iy>=NY(arg1)))) continue;
			loc_t* loc=0;
			if(NY(p1)==2&&NX(p1)>2){//first parameter is loc
				if(p1==p1_save){
					loc=loc_save;
				}
				if(!loc){
					loc=d2loc(p1);
					if(!loc_save){
						loc_save=loc;
						p1_save=p1;
					}
				}
			}
			if(!p2){//single parameter
				if(loc){//is loc
					if(loc->nloc>100000){/*if too many points, we draw it. */
						drawloc("loc", loc, NULL, title1, "x", "y", "%s[%02d]", title1, id++);
					} else{/*we plot individual points. */
						plot_points("loc", (plot_opts){.ngroup=1,.loc=&loc},  title1, "x", "y", "%s[%02d]", title1, id++);
					}
					if(loc!=loc_save){
						locfree(loc);
					}
				}else if(p1->nx>1 && p1->ny>1){//map
					map_t* data=d2map(p1);
					drawmap("map", data, NULL, title1, "x", "y", "%s[%02d]", title1, id++);
					mapfree(data);
				}else{
					plot_points("points", (plot_opts){.ngroup=1, .dc=arg1} , title1, "x", "y", "%s[%02d]", title1, id++);
				}
			}else{//two parameter
				if(loc&&p2&&p2->nx&&p2->ny){
					drawopd("opd", loc, p2, NULL, title1, "x", "y", "%s[%02d]", title1, id++);
				}
			}
		}
	}
	if(loc_save){
		locfree(loc_save);
	}
	cellfree(arg2);
	cellfree(arg1);
	draw_final(1);
}

