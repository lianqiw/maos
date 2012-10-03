/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#if defined(__CYGWIN__)
void GetTempPath(long, char*);
#endif
void init_process(void);
int    get_usage_running(void);
double get_usage_load(void);
double get_usage_mem(void);
double get_usage_cpu(void);

const char *get_job_progname(void);
int get_job_mem(void);
double get_job_launchtime(int pid);

int get_cpu_avail(void);
int read_usage_cpu(long *user, long *tot);
void wait_cpu(int nthread);
double read_self_cpu(void);
extern int NCPU;  /**<True number of cores*/
extern int NTHREAD; /**<Number of hyper threads. may be larger than NCPU*/
extern int TCK;
extern long NMEM;
extern const char *HOME;/*the user home */
extern const char *TEMP;/*the temporary folder */
extern const char *USER;/*the user name */
int get_ncpu(void);
#endif
