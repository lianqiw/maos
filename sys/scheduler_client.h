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
#ifndef AOS_SCHEDULER_CLIENT_H
#define AOS_SCHEDULER_CLIENT_H
#include "scheduler.h"
/**
   \file scheduler_client.h
   Contains routines that will be used to talk to the scheduler.
*/
extern int nhost;
extern char** hosts;
extern int hid;
extern uint16_t PORT;
/*called by maos */
int scheduler_listen(void(*fun)(int));
int scheduler_start(char *path, int nthread, int waiting);
int scheduler_wait(void);
int scheduler_finish(int status);
int scheduler_report(STATUS_T *status);
/*called by monitor */
int scheduler_launch_drawdaemon(char *fifo);
char* scheduler_get_drawdaemon(int pid, int direct);
int scheduler_launch_exe(const char *host, int argc, const char *argv[]);
/*Handling backtrace*/
char* call_addr2line(const char *buf);
void print_backtrace_symbol(void *const *buffer, int size);
#endif
