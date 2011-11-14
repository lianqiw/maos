/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

/**
   Convert bin files to fits files. Supports matrix and cellarray of matrix.
*/

#include "../lib/aos.h"

int main(int argc, char *argv[]){
    if(argc==1){
	info("Usage: %s file1.bin file2.bin \n", argv[0]);
	exit(0);
    }
    for(int iarg=1; iarg<argc; iarg++){
	const char *fn=argv[iarg];
	char fn2[strlen(fn)+1]; strcpy(fn2, fn);
	char *tmp=strstr(fn2, ".bin");
	if(!tmp){
	    warning("Entry %s has invalid format, skip it\n", fn);
	    continue;
	}
	strcpy(tmp, ".fits");
	info("Copying from %s to %s\n", fn, fn2);
	dcell *temp=dcellread(fn);
	dcellwrite(temp, fn2);
	dcellfree(temp);
    }
}
