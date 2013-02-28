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
#if defined(__linux__)
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <dirent.h>
#include <signal.h>
#include "../sys/sys.h"
/**
   \file record_cpu.c
   
   Records the cpu consumption of a process.

   Usage: record_cpu pid
*/
int main(int argc, char**argv){
    if(argc<2){
	error("Usage: record_cpu PID second\n");
    }
    long pid=strtol(argv[1],NULL,10);
    double sec;
    if(argc>2){
	sec=strtod(argv[2], NULL);
    }else{
	sec=0.1;
    }
    unsigned int usec=(int)(sec*1e6);
    char fn[PATH_MAX];
    snprintf(fn,PATH_MAX,"/proc/%ld/stat", pid);
    FILE *fp=fopen(fn,"r");
    if(!fp){
	return 1;
    }
    long stime,utime;
    double tck=sysconf(_SC_CLK_TCK)*sec;
    snprintf(fn,PATH_MAX,"%ld.cpu",pid);
    info("Recording PID %ld every %g second, and save to %s\n", pid, sec, fn);
    FILE *fpout=fopen(fn,"w");
    if(fscanf(fp,"%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %ld %ld",
	      &stime, &utime)!=2){
	error("Unable to read\n");
    }
    long last=utime+stime;
    rewind(fp);
    usleep(usec);
    setbuf(fp,NULL);
    for(int count=0; count<1000;count++){
	if(fscanf(fp,"%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %ld %ld",
		  &stime, &utime)!=2){
	    error("Unable to read\n");
	}
	fprintf(fpout,"%g\n",(double)(stime+utime-last)/tck);
	fflush(fpout);
	rewind(fp);
	last=stime+utime;
	usleep(usec);
    }
    fclose(fpout);
}
#endif
