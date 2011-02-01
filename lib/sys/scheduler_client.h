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
#ifndef AOS_SCHEDULER_CLIENT_H
#define AOS_SCHEDULER_CLIENT_H
#define MAXMSG  512
#include "scheduler_server.h"
void start_scheduler(int argc, char **argv);
int init_sockaddr (struct sockaddr_in *name, const char *hostname, uint16_t port);
int scheduler_connect(int ihost, int block,int mode);
int scheduler_connect_self(int block,int mode);
#ifdef __CYGWIN__
#define print_backtrace(A)
#endif
//called by maos
int scheduler_start(char *path, int nthread, int waiting);
int scheduler_wait(void);
void scheduler_finish(int status);
void scheduler_report(STATUS_T *status);
//called by monitor
int scheduler_kill_job(int ihost,int pid);
int scheduler_remove_job(int ihost, int pid);
int scheduler_launch_drawdaemon(char *fifo);
char* scheduler_get_drawdaemon(int pid, int direct);
void  scheduler_shutdown(int *sock, int mode);
#endif
