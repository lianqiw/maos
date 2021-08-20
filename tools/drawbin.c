/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	/*launch scheduler if it is not already running. */
	dcell *arg1=0, *arg2=0;
	const char *title1=NULL, *title2=NULL;
	if(argc>1){
		arg1=dcellread("%s", argv[1]);
		info("arg1 is %ldx%ld\n", arg1->nx, arg1->ny);
		if((title1=strrchr(argv[1],'/'))){
			title1++;
		}else{
			title1=argv[1];
		}
	}
	if(argc>2){
		arg2=dcellread("%s",argv[2]);
		info("arg2 is %ldx%ld\n", arg2->nx, arg2->ny);
		if((title2=strrchr(argv[2], '/'))){
			title2++;
		} else{
			title2=argv[2];
		}
	}

	long nx=MAX(NX(arg1), arg2?NX(arg2):0);
	long ny=MAX(NY(arg1), arg2?NY(arg2):0);
	int id=0;
	loc_t *loc_save=0;
	dmat* p1_save=0;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			dmat *p1=(dmat*)PR(arg1, ix, iy);
			dmat *p2=arg2?(dmat*)PR(arg2, ix, iy):NULL;
			if(!p1) continue;
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
						drawloc("loc", loc, NULL, title1, "x", "y", "%s[%d]", title1, id++);
					} else{/*we plot individual points. */
						plot_points("loc", 1, &loc, NULL, NULL, NULL, NULL, NULL, NULL, title1, "x", "y", "%s[%d]", title1, id++);
					}
					if(loc!=loc_save){
						locfree(loc);
					}
				}else if(p1->nx>1 && p1->ny>1){//map
					map_t* data=d2map(p1);
					drawmap("map", data, NULL, title1, "x", "y", "%s[%d]", title1, id++);
					mapfree(data);
				}else{
					plot_points("points", 1, NULL, arg1, NULL, NULL, NULL, NULL, NULL, title1, "x", "y", "%s[%d]", title1, id++);
				}
			}else{//two parameter
				if(loc&&p2&&p2->nx&&p2->ny){
					drawopd("opd", loc, p2, NULL, title2, "x", "y", "%s[%d]", title2, id++);
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

