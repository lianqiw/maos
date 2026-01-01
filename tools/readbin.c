/*
  Copyright 2009-2026 Lianqi Wang
  
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
		"readbin file.bin [file2.bin ...] \n"
	);
}
CHECK_ARG(2)
void cell_show(cell *var, const char* format, ...){
	format2fn;
	if(var){
		info("%s is %ldx%ld of type %u\n", fn, NX(var), NY(var), var->id);
		if(iscell(var)){
			for(long iy=0; iy<NY(var); iy++){
				for(long ix=0; ix<NX(var); ix++){
					cell_show(P(var,ix,iy), "%s[%ld,%ld]", fn, ix, iy);
				}
			}
		}
	}
}
/**
   wrap readbin() to an executable for debugging.
*/
int main(int argc, char *argv[]){
	if(argc<2){
		usage();
		return 0;
	}
	for(int iarg=1; iarg<argc; iarg++){
		cell* var=readbin("%s", argv[iarg]);
		cell_show(var, "%s", argv[iarg]);
		cellfree(var);
	}
}
