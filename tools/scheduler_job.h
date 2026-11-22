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
#ifndef AOS_SCHEDULER_JOB_H
#define AOS_SCHEDULER_JOB_H
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include "../sys/sys.h"
/**
   Struct to hold information of running jobs.
*/
typedef struct RUN_T{
	struct RUN_T* next;
	char* exe; /*Path to the executable.*/
	char* path;/*Job path and Job arguments.*/
	char* path0;/*same as path, with fields separated by \n instead of space*/
	status_t status;
	time_t last_time; //time of last report
	int pid;
	int pidnew;//the new pid
	int sock;//socket for primary maos connection. -1 if not connected.
	int sock_cmd;//socket for secondary maos connection for maos_command(). -1 if not connected.
	int nthread;//number of threads
	int ngpu;//number of gpus requested.
	int counted;//nthread and ngpu are counted in nused_cpu and nused_gpu.
}RUN_T;
extern int NGPU;
void nrun_add(RUN_T *irun);
void nrun_sub(RUN_T *irun);
int  nrun_get(int type);
RUN_T* running_add(int pid, int sock);
void running_remove(int pid, int status);
RUN_T* running_get_by_pid(int pid);
RUN_T* running_get_by_sock(int sock);
RUN_T* running_head(void);
void running_kill(int pid);
void runned_remove(int pid);
void runned_restart(int pid);
RUN_T* runned_head(void);
int  check_job_exists(RUN_T *irun);
void check_jobs(void);
void process_queue(void);
char* remove_endl(const char* path0);
void queue_new_job(const char* exename, const char* execmd);
/*Implemented by the user of this file (scheduler.c)*/
void monitor_send(RUN_T* irun, const char* path);
#endif
