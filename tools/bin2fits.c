/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file bin2fits.c
   Convert bin files to fits files. Supports matrix and zfarray of matrix.
*/

#include "../lib/aos.h"

int main(int argc, char *argv[]){
    if(argc==1){
	dbg("Usage: %s file1.bin file2.bin \n", argv[0]);
	exit(0);
    }
    int jarg=1;
    char **header=NULL;
    if(!strcmp(argv[jarg], "gc")){
	jarg++;
	dbg("Creating headers for GC\n");
	header=mymalloc(25,char*);
	char tmp[320];
	double wvl[5]={0.9, 0.975, 1, 1.025, 1.1};
	for(int iwvl=0; iwvl<5; iwvl++){
	    wvl[iwvl]*=2.12e-6;
	}
	int psfgrid=1024;
	double psfsum=6.34772385106148;
	double dx=1./16.;
	for(int istep=0; istep<5; istep++){
	    for(int iwvl=0; iwvl<5; iwvl++){
		snprintf(tmp, 320, 
			 "Turbulence: r0=0.1987m, l0=30m\n"
			 "Wavelength:   %gm\n"
			 "OPD Sampling: %gm\n"
			 "FFT Grid: %dx%d\n"
			 "PSF Sampling: %g arcsec\n"
			 "PSF Sum to %g\n"
			 "Exposure: %gs\n",
			 wvl[iwvl], dx, psfgrid, psfgrid,  wvl[iwvl]/(dx*psfgrid)*206265,
			 psfsum,  (double)((istep+1)*4));
		header[iwvl+istep*5]=strdup(tmp);
	    }
	}
    }
    for(int iarg=jarg; iarg<argc; iarg++){
	const char *fn=argv[iarg];
	char fn2[strlen(fn)+1]; strcpy(fn2, fn);
	char *tmp=strstr(fn2, ".bin");
	if(!tmp){
	    warning("Entry %s has invalid format, skip it\n", fn);
	    continue;
	}
	strcpy(tmp, ".fits");
	info("Copying from %s to %s\n", fn, fn2);
	dcell *temp=dcellread("%s", fn);
	if(header){
	    for(int i=0; i<temp->nx*temp->ny; i++){
		temp->p[i]->header=strdup(header[i]);
	    }
	}
	writebin(temp, "%s", fn2);
	dcellfree(temp);
    }
}
