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


#ifndef AOS_PROC_H
#define AOS_PROC_H
#include "common.h"
/**
   \file process.h
   A few routines handles process resource.
*/
#if defined(__CYGWIN__)
void GetTempPath(long, char*);
#endif
void init_process(void);
void free_process(void);
int    get_usage_running(void);
double get_usage_load(void);
double get_usage_mem(void);
double get_usage_cpu(void);

int get_job_progname(char *res, int nres, int pid);
size_t get_job_mem(void);
double get_job_launchtime(int pid);

int get_cpu_avail(void);
int read_cpu_counter(long *user, long *tot);
double read_self_cpu(void);
extern int NCPU;  /**<True number of cores*/
extern int MAXTHREAD;/**<Number of hyper threads. may be larger than NCPU*/
extern int NTHREAD;  /**<Number of threads for simulation*/
extern int TCK;
extern int PID;
extern long NMEM;
extern const char *HOME;/*the user home */
extern const char *USER;/*the user name */
extern char* HOST;/*the host name */
extern char* TEMP;/*the temporary folder */
extern char* DIRCACHE; /*the caching folder */
extern char* DIRLOCK; /*the lock file folder*/
extern char* DIREXE;/*the exe path*/
extern char* DIRSTART; /*the startup directory*/
extern char* DIROUT; /*the result directory (set by maos)*/
extern char *EXENAME;/*the name of the executable.*/
int get_ncpu(void);
void set_dirout(const char* dir);
#define print_mem(A,...) info(A ":\t%.1f MiB.\n", ##__VA_ARGS__, get_job_mem()/1024.);
#endif
