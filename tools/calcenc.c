/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
 * Standalone to postprocess MAOS PSF output and compute enclosed energe
 */
#include <unistd.h>
#include <getopt.h>

#include "../lib/aos.h"

static void calcenc(const char *fn, dmat *dvec, int type, int nthread){
    char fnout[PATH_MAX];
    char *suffix=strstr(fn, ".bin");
    if(!suffix) suffix=strstr(fn, ".fits");
    if(!suffix) error("%s has wrong suffix\n", suffix);
    memcpy(fnout, fn, suffix-fn+1);
    fnout[suffix-fn]='\0';
    switch(type){
    case -1:
	strcat(fnout, "_azmean.bin");
	break;
    case 0:
	strcat(fnout, "_ensquare.bin");
	break;
    case 1:
	strcat(fnout, "_encircle.bin");
	break;
    case 2:
	strcat(fnout, "_enstrip.bin");
	break;
    default:
	error("Not implemented\n");
    }
    if(exist(fnout)){
	if(fmtime(fn)<fmtime(fnout)){
	    info2("%s: skip\n", fn);
	    return;
	}
    }
    file_t *fp=zfopen(fn, "r");
    file_t *fpout=zfopen(fnout, "w");
    info2("%s --> %s\n", fn, fnout);
    header_t header;
    read_header(&header, fp);
    int free_dvec=0;
    if(iscell(header.magic)){
	write_header(&header, fpout);
	for(long i=0; i<header.nx*header.ny; i++){
	    info2("%ld of %ld\n", i, (long)(header.nx*header.ny));
	    dmat *psf=dreaddata(fp, &header);
	    if(!dvec){
		free_dvec=1;
		dvec=dlinspace(0, 1, psf->nx);
	    }
	    dmat *enc=denc(psf, dvec, type, nthread);
	    dwritedata(fpout, enc);
	    dfree(psf); dfree(enc);
	}
    }else{
	int nenc=0;
	dmat **encs=NULL;
	do{
	    dmat *psf=dreaddata(fp, &header);
	    if(!dvec){
		free_dvec=1;
		dvec=dlinspace(0, 1, psf->nx);
	    }
	    dmat *enc=denc(psf, dvec, type, nthread);
	    nenc++;
	    encs=realloc(encs, nenc*sizeof(dmat*));
	    encs[nenc-1]=enc;
	}while(!read_header2(&header, fp));
	dcell *encs2=dcellnew(nenc, 1);
	memcpy(encs2->p, encs, sizeof(dmat*)*nenc);
	dcellwritedata(fpout, encs2);
    }
    if(free_dvec) dfree(dvec);
    zfclose(fp);
    zfclose(fpout);
}
static void usage(){
    info2("Usage: calcenc [options] PSF1.bin PSF2.bin \n"
	  "Options are:\n"
	  "-n num --nthread: Maximum number of threads\n"
	  "-d num --diam   : Maximum diameter to compute enclosed energy. radius if type==-1\n"
	  "-s num --step   : Step in diamter\n"
	  "-t num --type   : type of computation:\n"
	  "\t\t\t0: ensquare energy inside a square\n"
	  "\t\t\t1: encircle energy inside a circle\n"
	  "\t\t\t2: enstriped energy inside a strip\n"
	  "\t\t\t-1: azimuthal average\n"
	  "-h         : Print this help"
	  );
}
int main(int argc, char *argv[]){
    if(argc==1){
	usage();
    }
    double rmax=0;
    double dstep=1;
    int nthread=NCPU;
    int ipos=0;
    int type=0;/*default is square */
    /*
    ARGOPT_T options[]={
	{"help",    'h', T_INT, 2, usage, NULL},
	{"nthread", 'n', T_INT, 1, &nthread, NULL},
	{"diam",    'd', T_DBL, 1, &rmax, NULL},
	{"step",    's', T_INT, 1, &dstep, NULL},
	{"type",    't', T_INT, 1, &type, NULL},
	{NULL, 0,0,0,NULL,NULL}};
    char *cmds=parse_argopt(argc, argv, options);
    */
    static struct option long_options[]={
	{"nthread", 1, 0, 'n'},
	{"diam",    1, 0, 'd'},
	{"step",    1, 0, 's'},
	{"type",    1, 0, 't'},
	{"help",    0, 0, 'h'},
	{NULL,0,0,0}
    };
    int option_index = 0;
    while(1){
	int c=getopt_long(argc, argv, "n:d:s:t:h",
			  long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'n':
	    nthread=strtol(optarg, NULL, 10);
	    break;
	case 'd':
	    rmax=strtod(optarg, NULL);
	    break;
	case 's':
	    dstep=strtod(optarg, NULL);
	    break;
	case 't':
	    type=strtol(optarg, NULL, 10);
	    break;
	case 'h':
	    usage();
	    break;
	default:
	    error("Invalid option.\n");
	}
    }
    ipos=optind;
    dmat *dvec=NULL;
    if(rmax>0){
	dvec=dlinspace(0, dstep, (long)ceil(rmax/dstep));
    }
    THREAD_POOL_INIT(nthread); ;
    for(; ipos<argc; ipos++){
	calcenc(argv[ipos], dvec, type, nthread);
    }
    dfree(dvec);
}
